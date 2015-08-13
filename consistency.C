/****************************************************************
 consistency.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Regex.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Array3D.H"
#include "BOOM/Time.H"
#include "SparseMatrix3D.H"
#include "PhyLib/Phylogeny.H"
#include "PairHMM/State.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  Application();
  int main(int argc,char *argv[]);
  typedef SparseMatrix3D::EntryList EntryList;
  typedef SparseMatrix3D::Entry Entry;
protected:
  Phylogeny *tree;
  String inDir, outDir;
  int numIterations;
  int numLeaves;
  Map<String,int> taxonToID; // not applicable to the phylogeny's ID's!
  Vector<String> taxonName;
  Vector<PhylogenyNode*> nodes;
  Regex filenameRegex, pathRegex;
  float threshold;
  Array2D<SparseMatrix3D*> matrices;
  Array1D<int> lengths;
  int getTaxonID(const String &taxon);
  void loadMatrices();
  void transform();
  SparseMatrix3D *transform(int taxonId1,int taxonId2);
  void normalize(SparseMatrix3D &);
  void save();
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
    filenameRegex("([^\/]+)-([^\/\\.]+)\\.post2")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line`
  CommandLine cmd(argc,argv,"i:s:");
  if(cmd.numArgs()!=3)
    throw String(
"\n\
consistency [options] <phylogeny.phy> <in-dir> <out-dir>\n\
   options are:\n\
      -i k = use k iterations of consistency transformation (default=2)\n\
      -s k = use k as sparseness threshold for matrices (default=0.001)\n\
\n\
");
  String phyFile=cmd.arg(0);
  inDir=cmd.arg(1);
  outDir=cmd.arg(2);
  numIterations=cmd.option('i') ? cmd.optParm('i').asInt() : 2;
  threshold=cmd.option('s') ? cmd.optParm('s').asFloat() : 0.001;
  //threshold=log(threshold);

  // Load phylogeny
  tree=new Phylogeny(phyFile);
  numLeaves=tree->getNumLeaves();
  matrices.resize(numLeaves,numLeaves);
  matrices.setAllTo(NULL);
  lengths.resize(numLeaves);
  lengths.setAllTo(0);

  Vector<PhylogenyNode*> allNodes;
  tree->gatherNodes(allNodes);
  int numNodes=allNodes.size();
  for(int i=0 ; i<numNodes ; ++i) {
    PhylogenyNode *node=allNodes[i];
    if(node->getNodeType()!=LEAF_NODE) continue;
    getTaxonID(node->getName());
  }

  // Load matrices
  loadMatrices();

  // Do It.
  transform();

  // Save
  save();

  return 0;
}



void Application::save()
{
  for(int i=0 ; i<numLeaves ; ++i) {
    for(int j=i ; j<numLeaves ; ++j) {
      String name1=taxonName[i], name2=taxonName[j];
      String outfile=outDir+"/"+name1+"-"+name2+".post2";
      matrices[i][j]->saveBinary(outfile);
    }
  }
}




int Application::getTaxonID(const String &taxon)
{
  if(!taxonToID.isDefined(taxon)) {
    taxonToID[taxon]=taxonName.size();
    taxonName.push_back(taxon);
    nodes.push_back(tree->findNode(taxon));
    if(tree->findNode(taxon)==NULL) {
      cout<<"can't find node "<<taxon<<endl;
      INTERNAL_ERROR;
    }
  }
  return taxonToID[taxon];
}



void Application::loadMatrices()
{  
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
    if(!taxonToID.isDefined(firstSpecies) ||
       !taxonToID.isDefined(secondSpecies)) continue;
    int tax1=getTaxonID(firstSpecies), tax2=getTaxonID(secondSpecies);
    SparseMatrix3D *M=SparseMatrix3D::loadBinary(inDir+"/"+filename);
    matrices[tax1][tax2]=M;
    matrices[tax2][tax1]=M->transpose();
    lengths[tax1]=M->getFirstDim()-1;
    lengths[tax2]=M->getSecondDim()-1;
  }
  for(int tax=0 ; tax<numLeaves ; ++tax) {
    int L=lengths[tax];
    SparseMatrix3D &M=*new SparseMatrix3D(L+1,L+1,3);
    matrices[tax][tax]=&M;
    for(int i=1 ; i<=L ; ++i) 
      //for(int j=1 ; j<=L ; ++j)
      M.addEntry(i,i,PHMM_MATCH,LOG_1);
  }
}



void Application::transform()
{
  for(int iter=0 ; iter<numIterations ; ++iter) {
    Array2D<SparseMatrix3D*> newM(numLeaves,numLeaves);
    newM.setAllTo(NULL);
    for(int tax1=0 ; tax1<numLeaves ; ++tax1)
      for(int tax2=0 ; tax2<numLeaves ; ++tax2) {
	//if(tax2==tax1) continue;
	newM[tax1][tax2]=transform(tax1,tax2);
      }
    for(int tax1=0 ; tax1<numLeaves ; ++tax1)
      for(int tax2=0 ; tax2<numLeaves ; ++tax2) {
	delete matrices[tax1][tax2];
	matrices[tax1][tax2]=newM[tax1][tax2];
      }	
    //matrices=newM;
  }
}



SparseMatrix3D *Application::transform(int xID,int yID)
{
  bool useDelta=false;//true;
  cout<<"transforming "<<taxonName[xID]<<":"<<taxonName[yID]<<endl;
  Time timer;
  timer.startCounting();
  PhylogenyNode *x=nodes[xID], *y=nodes[yID];
  if(x->getNodeType()!=LEAF_NODE || y->getNodeType()!=LEAF_NODE) return NULL;
  int Lx=lengths[xID], Ly=lengths[yID];
  Array3D<float> M(Lx+1,Ly+1,3);
  M.setAllTo(LOG_0);

  // Match states:
  for(int i=1 ; i<=Lx ; ++i) {
    Array1D< Vector<float> > sums(Ly+1);
    for(int zID=0 ; zID<numLeaves ; ++zID) {
      PhylogenyNode *z=nodes[zID];
      int Lz=lengths[zID];
      float delta=useDelta ? 
	tree->getSpannedDistance(x->getID(),y->getID(),z->getID()) : 0;
      SparseMatrix3D &Mxz=*matrices[xID][zID], &Mzy=*matrices[zID][yID];
      EntryList &row=Mxz(i,PHMM_MATCH);
      EntryList::iterator cur=row.begin(), end=row.end();
      for(; cur!=end ; ++cur) {
	Entry &e1=*cur;
	int k=e1.y;
	float v1=e1.value;
	EntryList &col=Mzy(k,PHMM_MATCH);
	EntryList::iterator cur=col.begin(), end=col.end();
	for(; cur!=end ; ++cur) {
	  Entry &e2=*cur;
	  int j=e2.y;
	  float v2=e2.value;
	  float product=v1+v2-delta;
	  sums[j].push_back(product);
	}
      }
    }
    for(int j=0 ; j<=Ly ; ++j)
      M[i][j][PHMM_MATCH]=sumLogProbs(sums[j])-numLeaves;
  }
  //cout<<M<<endl;

  // Insertion states:
  for(int i=0 ; i<=Lx ; ++i) {
    Array1D< Vector<float> > sums(Ly+1);
    for(int zID=0 ; zID<numLeaves ; ++zID) {
      PhylogenyNode *z=nodes[zID];
      int Lz=lengths[zID];
      float delta=useDelta ? 
	tree->getSpannedDistance(x->getID(),y->getID(),z->getID()) : 0;
      SparseMatrix3D &Mxz=*matrices[xID][zID], &Mzy=*matrices[zID][yID];
      EntryList &row=Mxz(i,PHMM_INSERT);
      EntryList::iterator cur=row.begin(), end=row.end();
      for(; cur!=end ; ++cur) {
	Entry &e1=*cur;
	int k=e1.y;
	float v1=e1.value;
	EntryList &col=Mzy(k,PHMM_MATCH);
	EntryList::iterator cur=col.begin(), end=col.end();
	for(; cur!=end ; ++cur) {
	  Entry &e2=*cur;
	  int j=e2.y;
	  float v2=e2.value;
	  float product=v1+v2-delta;
	  sums[j].push_back(product);
	}
      }
    }
    for(int j=0 ; j<=Ly ; ++j) 
      M[i][j][PHMM_INSERT]=sumLogProbs(sums[j])-numLeaves;
  }

  // Deletion states:
  for(int i=0 ; i<=Lx ; ++i) {
    Array1D< Vector<float> > sums(Ly+1);
    for(int zID=0 ; zID<numLeaves ; ++zID) {
      PhylogenyNode *z=nodes[zID];
      int Lz=lengths[zID];
      float delta=useDelta ? 
	tree->getSpannedDistance(x->getID(),y->getID(),z->getID()) : 0;
      SparseMatrix3D &Mxz=*matrices[xID][zID], &Mzy=*matrices[zID][yID];
      EntryList &row=Mxz(i,PHMM_MATCH);
      EntryList::iterator cur=row.begin(), end=row.end();
      for(; cur!=end ; ++cur) {
	Entry &e1=*cur;
	int k=e1.y;
	float v1=e1.value;
	EntryList &col=Mzy(k,PHMM_DELETE);
	EntryList::iterator cur=col.begin(), end=col.end();
	for(; cur!=end ; ++cur) {
	  Entry &e2=*cur;
	  int j=e2.y;
	  float v2=e2.value;
	  float product=v1+v2-delta;
	  sums[j].push_back(product);
	}
      }
    }
    for(int j=0 ; j<=Ly ; ++j) 
      M[i][j][PHMM_DELETE]=sumLogProbs(sums[j])-numLeaves;
  }
  //cout<<M<<endl;
  timer.stopCounting();
  cout<<timer.elapsedTime()<<endl;
  if(useDelta) {
    SparseMatrix3D *newM=new SparseMatrix3D(M,0.0);
    normalize(*newM);
    SparseMatrix3D *newerM=new SparseMatrix3D(*newM,threshold);
    delete newM;
    return newerM;
  }	
  SparseMatrix3D *newM=new SparseMatrix3D(M,threshold);

  {//###
    normalize(*newM);
    SparseMatrix3D *newerM=new SparseMatrix3D(*newM,threshold);
    delete newM;
    return newerM;
  }//###

  return newM;
}



void Application::normalize(SparseMatrix3D &M)
{
  int L1=M.getFirstDim()-1, L2=M.getSecondDim()-1;

  // Match and insertion states:
  Array1D< Vector<float> > V(L2+1);
  for(int i=0 ; i<=L1 ; ++i) {
    EntryList &matches=M(i,PHMM_MATCH);
    EntryList &insertions=M(i,PHMM_INSERT);
    EntryList::iterator curM=matches.begin(), endM=matches.end();
    EntryList::iterator curI=matches.begin(), endI=matches.end();
    for(; curM!=endM ; ++curM) {
      Entry &e=*curM;
      V[e.y].push_back(e.value);
    }
    for(; curI!=endI ; ++curI) {
      Entry &e=*curI;
      V[e.y].push_back(e.value);
    }
  }
  Array1D<float> sums(L2+1);
  for(int j=0 ; j<=L2 ; ++j) sums[j]=sumLogProbs(V[j]);
  for(int i=0 ; i<=L1 ; ++i) {
    EntryList &matches=M(i,PHMM_MATCH);
    EntryList &insertions=M(i,PHMM_INSERT);
    EntryList::iterator curM=matches.begin(), endM=matches.end();
    EntryList::iterator curI=matches.begin(), endI=matches.end();
    for(; curM!=endM ; ++curM) {
      Entry &e=*curM;
      e.value-=sums[e.y];
    }
    for(; curI!=endI ; ++curI) {
      Entry &e=*curI;
      e.value-=sums[e.y];
    }
  }

  // Match and deletion states:
  for(int i=0 ; i<=L1 ; ++i) {
    Vector<float> V;
    EntryList &matches=M(i,PHMM_MATCH);
    EntryList &deletions=M(i,PHMM_DELETE);
    EntryList::iterator curM=matches.begin(), endM=matches.end();
    EntryList::iterator curD=matches.begin(), endD=matches.end();
    for(; curM!=endM ; ++curM) {
      Entry &e=*curM;
      V.push_back(e.value);
    }
    for(; curD!=endD ; ++curD) {
      Entry &e=*curD;
      V.push_back(e.value);
    }
    float sum=sumLogProbs(V);
    for(; curM!=endM ; ++curM) {
      Entry &e=*curM;
      e.value-=sum;
    }
    for(; curD!=endD ; ++curD) {
      Entry &e=*curD;
      e.value-=sum;
    }
  }

  // ### THIS MAY NEED TO BE DONE  ITERATIVELY...SHOULD ADD SOME CODE HERE
  //     TO DUMP THE SUMS AND SEE IF THE SECOND NORMALIZATION SCREWED UP THE
  //     RESULTS OF THE FIRST ONE...
}



