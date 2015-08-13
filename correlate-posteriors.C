/****************************************************************
 correlate-posteriors.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Regex.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/VectorSorter.H"
#include "AVES.H"
#include "SparseMatrix3D.H"
#include "BranchHMM.H"
#include "ResidueOrthologyGraph.H"
#include "CollapsedOrthologyMatrix.H"
using namespace std;
using namespace BOOM;

typedef SparseMatrix3D::EntryList EntryList;

class Application : public AVES {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  BranchHMM *hmm;
  STATE matchState;
  String outfile;
  Array2D<SparseMatrix3D*> matrices;
  Regex filenameRegex, pathRegex;
  bool shouldDumpMatrices;
  float contrastRatio;
  float zero;
  void loadMatrices(const String &dirname);
  virtual void initBranches();
  void computeCorrelations(const String &filename);
  void dumpMatrices(const String &filename);
  void dumpMatrices2(const String &filename);
};



int main(int argc,char *argv[])
{
  try
    {
      Application app;
      return app.main(argc,argv);
    }
  catch(const char *p)
    {
      cerr << p << endl;
    }
  catch(const string &msg)
    {
      cerr << msg.c_str() << endl;
    }
  catch(const exception &e)
    {
      cerr << "STL exception caught in main:\n" << e.what() << endl;
    }
  catch(...)
    {
      cerr << "Unknown exception caught in main" << endl;
    }
  return -1;
}



Application::Application()
  : pathRegex("(.*\/)([^\/]+).post2"),
    filenameRegex("([^\/-]+)-([^\/\\.-]+)\\.post2")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"dc:z:");
  if(cmd.numArgs()!=3)
    throw String("\n\
correlate-posteriors [s] <*.aln> <library> <*.lambda>\n\
   where:\n\
     -d = dump posterior matrices\n\
     -c X = use contrast ratio X (0.0 - 1.0), default=1 (unenhanced)\n\
     -z X = use X as sparsity threshold when dumping matrices\n\
");
  String alignmentFile=cmd.arg(0);
  String libraryDir=cmd.arg(1);
  String lambdaFile=cmd.arg(2);
  shouldSeed=false;
  baselineMode=false;
  usePosteriors=false;
  shouldHillClimb=false;
  wantFuncDollo=false;
  explicitHistories=true;
  shouldDumpModel=false;
  useHirschberg=true;
  useContentSensor=false;
  usePrescan=true;
  noPrescanOnDownpass=true;
  shouldDumpMatrices=cmd.option('d');
  contrastRatio=cmd.option('c') ? cmd.optParm('c').asFloat() : 1.0;
  zero=cmd.option('z') ? cmd.optParm('z').asFloat() : 0.0;

  modelCompiler->setGCthreshold(LARGEST_INTEGER);
  LambdaAPI &lambda=modelCompiler->getLambda();
  lambda.getGC().setSilence(true);
  modelCompiler->parse(lambdaFile);
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  transducerTemplate=modelCompiler->instantiate(ctor);
  hmm=templateInstantiator->instantiate(transducerTemplate,1.0,NULL);
  int numStates=hmm->getNumStates();
  for(int i=0 ; i<numStates ; ++i)
    if(hmm->getStateType(i)==PHMM_MATCH) matchState=i;
  tree=modelCompiler->getGuideTree(transducerTemplate);
  tree->gatherNodes(phylogenyNodes);
  tree->constructBranches();
  tree->computePairwiseDistances();

  // Link taxa to phylogeny nodes
  numTaxa=phylogenyNodes.size();
  taxa.resize(numTaxa);
  for(int i=0 ; i<numTaxa ; ++i) {
    PhylogenyNode *node=phylogenyNodes[i];
    int ID=node->getID();
    nameToTaxon[node->getName()]=ID;
    Taxon &taxon=taxa[ID];
    taxon.setNode(node);
    node->getDecoration()=&taxon;
  }
  alignmentBuilder=
    new AlignmentBuilder(tree->getRoot(),alphabet,gapSymbol,numTaxa);
  gatherCladeMembers(); // ### ?

  // Instantiate BranchHMM's and substitution matrices on tree
  initBranches();

  // Load posterior matrices
  loadMatrices(libraryDir);

  computeCorrelations(alignmentFile);
  cerr<<"dumping matrices"<<endl;
  if(shouldDumpMatrices) dumpMatrices(alignmentFile);
  cerr<<"done."<<endl;

  return 0;
}



/****************************************************************
                     Application::loadMatrices()
 ****************************************************************/
void Application::loadMatrices(const String &inDir)
{  
  matrices.resize(numTaxa,numTaxa);
  Array1D<int> lengths(numTaxa);
  Vector<String> files;
  if(!File::getFileList(inDir,files)) 
    throw String("Can't get dir listing for ")+inDir;
  Vector<String>::iterator cur=files.begin(), end=files.end();
  for(; cur!=end ; ++cur) {
    const String &filename=*cur;
    String baseName, path;
    if(pathRegex.match(filename)) {
      path=pathRegex[1];
      baseName=pathRegex[2];
    }
    else baseName=filename;
    if(!filenameRegex.match(baseName)) continue;
    String firstSpecies=filenameRegex[1], secondSpecies=filenameRegex[2];
    if(!nameToTaxon.isDefined(firstSpecies)) continue;
    if(!nameToTaxon.isDefined(secondSpecies)) continue;
    int tax1=nameToTaxon[firstSpecies], tax2=nameToTaxon[secondSpecies];
    cout<<"loading "<<filename<<endl;
    SparseMatrix3D *M=SparseMatrix3D::loadBinary(inDir+"/"+filename);
    matrices[tax1][tax2]=M;
    matrices[tax2][tax1]=M->transpose();
    lengths[tax1]=M->getFirstDim();
    lengths[tax2]=M->getSecondDim();
  }
  for(int tax=0 ; tax<numTaxa ; ++tax) {
    if(!taxa[tax].isLeaf()) continue;
    int L=lengths[tax];//taxa[tax].getSeqLen();
    taxa[tax].setSeqLen(L);
    SparseMatrix3D &M=*new SparseMatrix3D(L+1,L+1,3);
    matrices[tax][tax]=&M;
    for(int i=1 ; i<=L ; ++i) 
      for(int j=1 ; j<=L ; ++j)
	M.addEntry(i,j,PHMM_MATCH,LOG_1);
  }
}



/****************************************************************
                         Application::initBranches()
 ****************************************************************/
void Application::initBranches()
{
  //tree->gatherBranches(branches);
  Vector<PhylogenyBranch*> tmp;
  tree->collectBranches(tmp);
  int n=tmp.size();
  branches.resize(n);
  for(int i=0 ; i<n ; ++i) branches[i]=*tmp[i];
  branchAttributes.resize(n);
  for(int i=0 ; i<n ; ++i) {
    PhylogenyBranch *branch=&branches[i];
    PhylogenyNode *parentNode=branch->getParent();
    PhylogenyNode *childNode=branch->getChild();
    branchAttributes[i]=BranchAttributes(NULL,branch);
    //cout<<"assigning "<<branch<<" -> "<<&branchAttributes[i]<<endl;
    branch->getDecoration()=&branchAttributes[i];
    tmp[i]->getDecoration()=&branchAttributes[i];
    Taxon *parent=&taxa[parentNode->getID()];
    parentNode->getDecoration()=parent;
    switch(parentNode->getNodeType()) 
      {
      case ROOT_NODE: parent->getIthBranch(0)=&branchAttributes[i]; break;
      case LEAF_NODE: break;
      case INTERNAL_NODE: {
	InternalNode *in=static_cast<InternalNode*>(parentNode);
	WhichChild c=in->whichChildIsThis(childNode);
	parent->getIthBranch(c)=&branchAttributes[i];
	}
	break;
      }
  }
}



/****************************************************************
                 Application::computeCorrelations()
 ****************************************************************/
void Application::computeCorrelations(const String &alignFile)
{
  ResidueOrthologyGraph rog(tree,taxa);
  Array2D< Array2D<float> > counts(numTaxa,numTaxa);
  for(int i=0 ; i<numTaxa ; ++i)
    for(int j=i+1 ; j<numTaxa ; ++j) {
      Array2D<float> &cnt=counts[i][j];
      cnt.resize(taxa[i].getSeqLen(),taxa[j].getSeqLen());
      cnt.setAllTo(0);
    }
  File f(alignFile);
  int numAlignments;
  numAlignments=f.readInt();


  //numAlignments=10000;


  cerr<<numAlignments<<" alignments"<<endl;
  for(int i=0 ; i<numAlignments ; ++i) {
    //cout<<"i="<<i<<endl;
    rog.loadConnections(f);
    if(f.eof()) break;//cerr<<"EOF"<<endl;
    for(int j=0 ; j<numTaxa ; ++j) {
      //cout<<"j="<<j<<endl;
      if(!taxa[j].isLeaf()) continue;
      for(int k=j+1 ; k<numTaxa ; ++k) {
	//cout<<"k="<<k<<endl;
	if(!taxa[k].isLeaf()) continue;
	IndexMap forward, backward;
	rog.getAlignment(j,k,forward,backward);
	//cout<<"a"<<endl;
	Array2D<float> &cnt=counts[j][k];
	//cout<<"b"<<endl;
	int L=forward.size();
	//cout<<"c "<<L<<endl;
	for(int x=0 ; x<L ; ++x) {
	  //cout<<"x="<<x<<endl;
	  int y=forward[x];
	  //cout<<"y="<<y<<endl;
	  if(y!=IndexMap::UNDEFINED) ++cnt[x][y];
	}
	//cout<<"A"<<endl;
      }
      //cout<<"B"<<endl;
    }
    //cout<<"C"<<endl;
  }
  f.close();
  for(int i=0 ; i<numTaxa ; ++i) {
    if(!taxa[i].isLeaf()) continue;
    for(int j=i+1 ; j<numTaxa ; ++j) {
      if(!taxa[j].isLeaf()) continue;
      Array2D<float> &cnt=counts[i][j];
      SparseMatrix3D *M=matrices[i][j];
      if(!M) {cout<<"i="<<i<<" "<<taxa[i].getName()<<" j="<<j<<" "<<taxa[j].getName()<<endl; INTERNAL_ERROR;}
      int L1=cnt.getFirstDim()-1, L2=cnt.getSecondDim()-1;
      for(int x=0 ; x<=L1 ; ++x) {
	Array2D<float>::RowIn2DArray<float> row=cnt[x];
	for(int y=0 ; y<=L2 ; ++y) {
	  row[y]/=numAlignments;
	  //if(row[y]>0) osObs<<x<<" "<<y<<" 0 "<<powf(row[y],0.2)<<endl;
	}
      }
      Array2D<float> expected(L1,L2);
      expected.setAllTo(0.0);
      for(int x=0 ; x<L1 ; ++x) {
	Array2D<float>::RowIn2DArray<float> row=cnt[x];
	SparseMatrix3D::EntryList el=(*M)(x+1,PHMM_MATCH);
	SparseMatrix3D::EntryList::iterator cur=el.begin(), end=el.end();
	for(; cur!=end ; ++cur) {
	  SparseMatrix3D::Entry e=*cur;
	  float obs=row[e.y-1];
	  float exponent=1.0;//0.5;
	  float value=exp(e.value*exponent);
	  if(value>1) value=1;
	  //if(value>0) osExp<<x<<" "<<e.y-1<<" 0 "<<powf(value,0.2)<<endl;
	  expected[x][e.y-1]=value;//exp(e.value);//value;
	}
      }
      for(int x=0 ; x<L1 ; ++x)
	for(int y=0 ; y<L2 ; ++y)
	  cout<<cnt[x][y]<<"\t"<<expected[x][y]<<endl;
    }
  }
}



void Application::dumpMatrices(const String &alignFile)
{
  ResidueOrthologyGraph rog(tree,taxa);
  for(PHMM_StateType t=0 ; t<3 ; ++t) {
    Array2D< Array2D<float> > counts(numTaxa,numTaxa);
    for(int i=0 ; i<numTaxa ; ++i)
      for(int j=i+1 ; j<numTaxa ; ++j) {
	Array2D<float> &cnt=counts[i][j];
	cnt.resize(taxa[i].getSeqLen(),taxa[j].getSeqLen());
	cnt.setAllTo(0);
      }
    File f(alignFile);
    int numAlignments;
    numAlignments=f.readInt();
    for(int i=0 ; i<numAlignments ; ++i) {
      rog.loadConnections(f);
      if(f.eof()) break;
      for(int j=0 ; j<numTaxa ; ++j) {
	if(!taxa[j].isLeaf()) continue;
	for(int k=j+1 ; k<numTaxa ; ++k) {
	  if(!taxa[k].isLeaf()) continue;
	  IndexMap forward, backward;
	  rog.getAlignment(j,k,forward,backward);
	  Array2D<float> &cnt=counts[j][k];
	  switch(t)
	    {
	    case PHMM_MATCH:
	      int L=forward.size();
	      for(int x=0 ; x<L ; ++x) {
		int y=forward[x];
		if(y!=IndexMap::UNDEFINED) ++cnt[x][y];
	      } break;
	    case PHMM_DELETE: {
	      int L=forward.size();
	      int prevY=0;
	      for(int x=0 ; x<L ; ++x) {
		int y=forward[x];
		if(y==IndexMap::UNDEFINED) ++cnt[x][prevY];
		else prevY=y;
	      }
	    } break;
	    case PHMM_INSERT: {
	      int L=backward.size();
	      int prevX=0;
	      for(int y=0 ; y<L ; ++y) {
		int x=backward[y];
		if(x==IndexMap::UNDEFINED) ++cnt[prevX][y];
		else prevX=x;
	      }
	    } break;
	    }
	}
      }
    }
    f.close();
    for(int i=0 ; i<numTaxa ; ++i) {
      if(!taxa[i].isLeaf()) continue;
      String name1=taxa[i].getName();
      for(int j=i+1 ; j<numTaxa ; ++j) {
	if(!taxa[j].isLeaf()) continue;
	String name2=taxa[j].getName();
	String obsName=String("obs-")+name1+"-"+name2+"-"+toString(t)+".txt";
	String expName=String("exp-")+name1+"-"+name2+"-"+toString(t)+".txt";
	ofstream osObs(obsName.c_str());
	ofstream osExp(expName.c_str());
	Array2D<float> &cnt=counts[i][j];
	SparseMatrix3D *M=matrices[i][j];
	int L1=cnt.getFirstDim()-1, L2=cnt.getSecondDim()-1;
	for(int x=0 ; x<=L1 ; ++x) {
	  Array2D<float>::RowIn2DArray<float> row=cnt[x];
	  for(int y=0 ; y<=L2 ; ++y) {
	    row[y]/=numAlignments;
	    if(row[y]>zero) osObs<<x<<" "<<y<<" 0 "
				 <<powf(row[y],contrastRatio)<<endl;
	  }
	}
	for(int x=0 ; x<L1 ; ++x) {
	  SparseMatrix3D::EntryList el=(*M)(x+1,t);
	  SparseMatrix3D::EntryList::iterator cur=el.begin(), end=el.end();
	  for(; cur!=end ; ++cur) {
	    SparseMatrix3D::Entry e=*cur;
	    float value=exp(e.value);
	    if(value>zero) osExp<<x<<" "<<e.y-1<<" 0 "
			     <<powf(value,contrastRatio)<<endl;
	  }
	}
      }
    }
  }
}



void Application::dumpMatrices2(const String &alignFile)
{
  ResidueOrthologyGraph rog(tree,taxa);
  Array2D< Array2D<float> > counts(numTaxa,numTaxa);
  for(int i=0 ; i<numTaxa ; ++i)
    for(int j=i+1 ; j<numTaxa ; ++j) {
      Array2D<float> &cnt=counts[i][j];
      cnt.resize(taxa[i].getSeqLen(),taxa[j].getSeqLen());
      cnt.setAllTo(0);
    }
  int numAlignments;
  for(PHMM_StateType t=0 ; t<3 ; ++t) {
    File f(alignFile);
    numAlignments=f.readInt();
    for(int i=0 ; i<numAlignments ; ++i) {
      rog.loadConnections(f);
      for(int j=0 ; j<numTaxa ; ++j) {
	if(!taxa[j].isLeaf()) continue;
	for(int k=j+1 ; k<numTaxa ; ++k) {
	  if(!taxa[k].isLeaf()) continue;
	  IndexMap forward, backward;
	  rog.getAlignment(j,k,forward,backward);
	  Array2D<float> &cnt=counts[j][k];
	  switch(t)
	    {
	    case PHMM_MATCH:
	      int L=forward.size();
	      for(int x=0 ; x<L ; ++x) {
		int y=forward[x];
		if(y!=IndexMap::UNDEFINED) ++cnt[x][y];
	      } break;
	    case PHMM_DELETE: {
	      int L=forward.size();
	      int prevY=0;
	      for(int x=0 ; x<L ; ++x) {
		int y=forward[x];
		if(y==IndexMap::UNDEFINED) ++cnt[x][prevY];
		else prevY=y;
	      }
	    } break;
	    case PHMM_INSERT: {
	      int L=backward.size();
	      int prevX=0;
	      for(int y=0 ; y<L ; ++y) {
		int x=backward[y];
		if(x==IndexMap::UNDEFINED) ++cnt[prevX][y];
		else prevX=x;
	      }
	    } break;
	    }
	}
      }
    }
    f.close();
  }
  for(int i=0 ; i<numTaxa ; ++i) {
    if(!taxa[i].isLeaf()) continue;
    String name1=taxa[i].getName();
    for(int j=i+1 ; j<numTaxa ; ++j) {
      if(!taxa[j].isLeaf()) continue;
      String name2=taxa[j].getName();
      String obsName=String("obs-")+name1+"-"+name2+"-"+toString(t)+".txt";
      String expName=String("exp-")+name1+"-"+name2+"-"+toString(t)+".txt";
      ofstream osObs(obsName.c_str());
      ofstream osExp(expName.c_str());
      Array2D<float> &cnt=counts[i][j];
      SparseMatrix3D *M=matrices[i][j];
      int L1=cnt.getFirstDim()-1, L2=cnt.getSecondDim()-1;
      for(int x=0 ; x<=L1 ; ++x) {
	Array2D<float>::RowIn2DArray<float> row=cnt[x];
	for(int y=0 ; y<=L2 ; ++y) {
	  row[y]/=numAlignments;
	  if(row[y]>zero) osObs<<x<<" "<<y<<" 0 "
			       <<powf(row[y],contrastRatio)<<endl;
	}
      }
      for(int x=0 ; x<L1 ; ++x) {
	SparseMatrix3D::EntryList el=(*M)(x+1,t);
	SparseMatrix3D::EntryList::iterator cur=el.begin(), end=el.end();
	for(; cur!=end ; ++cur) {
	  SparseMatrix3D::Entry e=*cur;
	  float value=exp(e.value);
	  if(value>zero)
	    osExp<<x<<" "<<e.y-1<<" 0 "<<powf(value,contrastRatio)<<endl;
	}
      }
    }
  }
}



