


OBSOLETE




/****************************************************************
 marginalize-class.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Regex.H"
#include "BOOM/SumLogProbs.H"
#include "AVES.H"
#include "SparseMatrix3D.H"
using namespace std;
using namespace BOOM;

enum CombiningStrategy {
  CS_MEAN,
  CS_POWER_MEAN,
  CS_WEIGHTED_POWER_MEAN,
  CS_PHYLO_SPAN,
  CS_PHYLO_LOGISTIC,
  CS_PHYLO_CRF
};


class Application : public AVES {
public:
  Application();
  int main(int argc,char *argv[]);
  typedef SparseMatrix3D::EntryList EntryList;
  typedef SparseMatrix3D::Entry Entry;
protected:
  CombiningStrategy combiningStrategy;
  String targetSpecies;
  String outDir;
  String wantClass;
  String wantFactor;
  bool wantForegroundDump;
  String foregroundDumpFile;
  String restrictInformant;
  int numAlpha;
  Regex filenameRegex, pathRegex;
  BranchHMM *hmm;
  bool skipBackground;
  float power; // used in Power Mean
  void topLevel(const String &infile,const String &outfile);
  void processDir();
  float combine_Ave(Vector<float> &);
  float combine_powerMean(Vector<float> &);
  float combine_phyloLogistic(Vector<float> &);
  float combine_phyloCRF(Vector<float> &);
  float combine_phyloSpan(Vector<float> &,Vector<int> &);
  float combine_weightedPower(Vector<float> &,Vector<float> &weights);

  //debugging:
  float combine_phyloSpan_old(Vector<float> &,Vector<int> &);
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
  : pathRegex("(.*\/)([^\/]+).post"),
    filenameRegex("([^\/]+)-([^\/]+)-([^\/\\.]+)\\.post")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line`
  CommandLine cmd(argc,argv,"dGBn:f:s:p:F:I:");
  if(cmd.numArgs()!=4)
    throw String(
"\n\
marginalize-class [options] <model.lambda> <target-species> <infile> <outfile>\n\
   options are:\n\
      -d = filename is actually a directory\n\
      -G = no garbage collection\n\
      -B = skip background states\n\
      -f file = dump graph of foreground classes (pooled) into file\n\
      -F factor = ignore all factors but this one\n\
      -n name = dump graph for functional class with given name\n\
                (e.g., 'B_0')\n\
      -s # = stategy for combining evidence: (default=1)\n\
             1 : arithmetic mean          2 : power mean\n\
             3 : weighted power mean      4 : phylo span\n\
             5 : phylo logistic           6 : phylo CRF\n\
      -p k = use k as power in power mean (default=2)\n\
      -I X = use only one informant, X\n\
\n\
");
  String lambdaFile=cmd.arg(0);
  targetSpecies=cmd.arg(1);
  String infile=cmd.arg(2);
  String outfile=cmd.arg(3);
  if(cmd.option('d')) outDir=infile;
  disableGC=cmd.option('G');
  if(disableGC) modelCompiler->setGCthreshold(LARGEST_INTEGER);
  numAlpha=alphabet.size();
  skipBackground=cmd.option('B');
  wantClass=cmd.option('n') ? cmd.optParm('n') : "";
  wantForegroundDump=cmd.option('f');
  if(wantForegroundDump) foregroundDumpFile=cmd.optParm('f');
  combiningStrategy=cmd.option('s') ? cmd.optParm('s').asInt()-1 : 0;
  power=cmd.option('p') ? cmd.optParm('p').asFloat() : 2;
  if(cmd.option('F')) wantFactor=cmd.optParm('F');
  if(cmd.option('I')) restrictInformant=cmd.optParm('I');
  
  // Load HMM
  //cout<<"executing "<<lambdaFile<<"..."<<endl;
  LambdaAPI &lambda=modelCompiler->getLambda();
  lambda.getGC().setSilence(true);
  modelCompiler->parse(lambdaFile);
  modelCompiler->deleteForeignObjects();
  lambda.checkGarbageLevel();
  
  //cout<<"instantiating model"<<endl;
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  transducerTemplate=modelCompiler->instantiate(ctor);
  //cout<<transducerTemplate->getNumStates()<<" states"<<endl;

  // Get phylogeny
  tree=modelCompiler->getGuideTree(transducerTemplate);
  tree->gatherNodes(phylogenyNodes);
  tree->constructBranches();

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
  alignmentBuilder=NULL;
  gatherCladeMembers();
  initBranches();

  // Instantiate a dummy HMM so we can lookup state classes
  hmm=templateInstantiator->instantiate(transducerTemplate,1.0,NULL);

  // Perform the marginalization
  if(cmd.option('d')) processDir();
  else topLevel(infile,outfile);

  return 0;
}



void Application::topLevel(const String &infile,const String &outfile)
{
  //cout<<"Loading sparse matrix "<<infile<<endl;
  SparseMatrix3D *M=SparseMatrix3D::loadBinary(infile);

  String baseName;
  if(pathRegex.match(infile)) baseName=pathRegex[2];
  else baseName=infile;
  if(!filenameRegex.match(baseName)) 
    throw String("Can't parse species from filename: ")+baseName;
  BranchEnd branchEnd;
  String factor=filenameRegex[1];
  String firstSpecies=filenameRegex[2], secondSpecies=filenameRegex[3];
  if(targetSpecies==secondSpecies) {
    SparseMatrix3D *Mt=M->transpose();
    delete M;
    M=Mt;
    branchEnd=CHILD;
  }
  else if(targetSpecies==firstSpecies) {
    branchEnd=PARENT;
  }
  else throw String("Target species \"")+targetSpecies+" does not match "
    +firstSpecies+" or "+secondSpecies;

  File out(outfile,"w");
  short L1=M->getFirstDim(), L2=M->getSecondDim(), numStates=M->getThirdDim();
  int numClasses=FunctionalClass::numClasses();

  int numForegroundClasses=0;
  for(short fc=0 ; fc<numClasses ; ++fc)  //###
    if(FunctionalClass(fc).fg_or_bg()==FOREGROUND) ++numForegroundClasses;

  ofstream *os=NULL;
  if(wantForegroundDump || !wantClass.isEmpty())
    os=new ofstream(foregroundDumpFile.c_str());
  Array1D<float> sums(numClasses);
  for(short x=0 ; x<L1 ; ++x) {
    sums.setAllTo(0.0);
    for(short q=0 ; q<numStates ; ++q) {
      FunctionalClass fc=hmm->getFunctionalClass(q,branchEnd);
      if(skipBackground && fc.isBackground()) continue;
      EntryList &entries=(*M)(x,q);
      EntryList::iterator cur=entries.begin(), end=entries.end();
      for(; cur!=end ; ++cur) {
	Entry &entry=*cur;
	//sums[fc]+=entry.value;
	sums[fc]=sumLogProbs(sums[fc],entry.value);
      }
    }
    out.write(x);
    for(short fc=0 ; fc<numClasses ; ++fc) {
      out.write(fc);
      out.write(sums[fc]);
    }	
    if(!wantClass.isEmpty()) {
      short fc=FunctionalClass::getClassByName(wantClass);
      (*os)<<x<<"\t"<<sums[fc]<<endl;
    }
    if(wantForegroundDump) {
      float sum=0.0;
      for(short fc=0 ; fc<numClasses ; ++fc)
	if(FunctionalClass(fc).fg_or_bg()==FOREGROUND) 
	  sum=sumLogProbs(sum,sums[fc]);
      (*os)<<x<<"\t"<<sum/numForegroundClasses<<endl;
    }
  }
  delete os;
  delete M;
}



void Application::processDir()
{
  FunctionalClass fcWant;
  if(!wantClass.isEmpty())
    fcWant=FunctionalClass::getClassByName(wantClass);
  int numClasses=FunctionalClass::numClasses();
  ofstream *os=NULL;
  if(wantForegroundDump || !wantClass.isEmpty())
    os=new ofstream(foregroundDumpFile.c_str());
  Vector<String> files;
  if(!File::getFileList(outDir,files)) 
    throw String("Can't get dir listing for ")+outDir;
  int targetLength=-1;
  Array1D< Vector<float> > uberSums;
  Vector<int> taxonIDs;
  Vector<float> weights;
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
    BranchEnd branchEnd;
    String factor=filenameRegex[1];
    /*
    if(wantFactor.length()>0 && !(factor==wantFactor || 
      factor==wantFactor+"fwd" || factor==wantFactor+"rev")) continue;
    */
    if(wantFactor.length()>0 && !(factor==wantFactor)) continue;
    String firstSpecies=filenameRegex[2], secondSpecies=filenameRegex[3];
    if(!(targetSpecies==firstSpecies || targetSpecies==secondSpecies)) 
      continue;
    if(!(targetSpecies==firstSpecies)) continue;//###
    int tax1id=nameToTaxon[firstSpecies], tax2id=nameToTaxon[secondSpecies];
    int otherTaxonID;
    SparseMatrix3D *M=SparseMatrix3D::loadBinary(outDir+"/"+filename);
    if(targetSpecies==secondSpecies) {
      if(!restrictInformant.empty() && firstSpecies!=restrictInformant) {
	delete M;
	continue;
      }
      SparseMatrix3D *Mt=M->transpose();
      delete M;
      M=Mt;
      branchEnd=CHILD;
      otherTaxonID=tax1id;
    }
    else if(targetSpecies==firstSpecies) {
      if(!restrictInformant.empty() && secondSpecies!=restrictInformant) {
	delete M;
	continue;
      }
      branchEnd=PARENT;
      otherTaxonID=tax2id;
    }
    PhylogenyNode *node1=taxa[tax1id].getNode(), *node2=taxa[tax2id].getNode();
    float branchLen=tree->distanceBetween(node1,node2);
    float weight=1/branchLen;
    if(weight<0) INTERNAL_ERROR;
    weights.push_back(weight);
    taxonIDs.push_back(otherTaxonID);
    String outfile=outDir+"/"+factor+"-"+targetSpecies+".marg";
    File out(outfile,"w");
    short L1=M->getFirstDim(), L2=M->getSecondDim();
    short numStates=M->getThirdDim();
    if(uberSums.size()<L1) {
      uberSums.resize(L1);
      targetLength=L1;
    }
    for(short x=0 ; x<L1 ; ++x) {
      Vector<float> logValues;
      for(short q=0 ; q<numStates ; ++q) {
	FunctionalClass fc=hmm->getFunctionalClass(q,branchEnd);
	if(fc.fg_or_bg()!=FOREGROUND) continue;
	if(!wantClass.isEmpty() && fc!=fcWant) continue;
	EntryList &entries=(*M)(x,q);
	EntryList::iterator cur=entries.begin(), end=entries.end();
	for(; cur!=end ; ++cur) {
	  Entry &entry=*cur;
	  logValues.push_back(entry.value);
	}
      }
      uberSums[x].push_back(sumLogProbs(logValues));
    }
    delete M;
  }
  for(short x=0 ; x<targetLength ; ++x) {
    float combined;
    switch(combiningStrategy)
      {
      case CS_MEAN:
	combined=combine_Ave(uberSums[x]);
	break;
      case CS_POWER_MEAN:
	combined=combine_powerMean(uberSums[x]);
	break;
      case CS_WEIGHTED_POWER_MEAN:
	combined=combine_weightedPower(uberSums[x],weights);
	break;
      case CS_PHYLO_LOGISTIC:
	combined=combine_phyloLogistic(uberSums[x]);
	break;
      case CS_PHYLO_CRF:
	combined=combine_phyloCRF(uberSums[x]);
	break;
      case CS_PHYLO_SPAN:
	combined=combine_phyloSpan(uberSums[x],taxonIDs);
	break;
      default: INTERNAL_ERROR;
      }
    //if(os) (*os)<<x<<"\t"<<exp(combined)<<endl;
    if(os) (*os)<<x<<"\t"<<exp(combined)<<endl;
  }
  delete os;
}



float Application::combine_Ave(Vector<float> &x)
{
  return sumLogProbs(x)-log(x.size());
}



float vectMax(Vector<float> &v) {
  float m=NEGATIVE_INFINITY;
  Vector<float>::iterator cur=v.begin(), end=v.end();
  for(; cur!=end ; ++cur)
    if(*cur>m) m=*cur;
  return m;
}


float Application::combine_powerMean(Vector<float> &v)
{
  const float k=power;
  if(k>=100) return vectMax(v);

  int nn=v.size();
  Vector<float> x;
  for(int i=0 ; i<nn ; ++i) {
    float f=v[i];
    if(isFinite(f)) x.push_back(f);
  }
  int n=x.size();
  if(n==0) return NEGATIVE_INFINITY;
  if(n==1) return x[0];
  double sum=1.0;
  if(k==0) { // geometric mean
    sum=0;
    for(int i=0 ; i<n ; ++i) sum+=x[i];
    return sum/n;
  }
  const float x0=x[0];
  for(int i=1 ; i<n ; ++i)
    sum+=exp(k*(x[i]-x0));
  sum=log(sum);
  sum/=k;
  sum+=x0;
  sum-=log(nn)/k;
  //if(sum>0) throw "Error in power mean: exponent may be too high?";
  return sum;
}



float Application::combine_weightedPower(Vector<float> &v,
					 Vector<float> &weights)
{
  const float k=power;
  if(k>=100) return vectMax(v);

  int nn=v.size();
  Vector<float> x, w;
  for(int i=0 ; i<nn ; ++i) {
    float f=v[i];
    if(isFinite(f)) { x.push_back(f); w.push_back(weights[i]); }
  }
  int n=x.size();
  if(n==0) return NEGATIVE_INFINITY;
  if(n==1) return x[0];
  double sum=0.0;
  if(k==0) { // geometric mean
    for(int i=0 ; i<n ; ++i)
      sum+=x[i]+log(w[i]); // ?
    return sum/n;
  }
  const float x0=x[0];
  float sumW=w[0];
  for(int i=1 ; i<n ; ++i) {
    sum+=exp(k*(x[i]-x0+log(w[i])));
    sumW+=w[i];
  }
  sum/=exp(k*log(w[0]));
  sum+=1;
  sum=log(sum);
  sum/=k;
  sum+=x0;
  sum+=log(w[0])/k;
  sum-=log(sumW)/k;
  //if(sum>0) throw "Error in power mean: exponent may be too high?";
  return sum;  
}



float Application::combine_phyloSpan_old(Vector<float> &X,
					 Vector<int> &taxonIDs)
{
  // First, get the minimum values for the ancestral nodes, to serve as
  // estimates for those unobservable values

  int targetID=nameToTaxon[targetSpecies];
  PhylogenyNode *targetNode=taxa[targetID].getNode();
  Array1D<float> inferred(numTaxa); // for internal nodes only (non-root)
  inferred.setAllTo(NEGATIVE_INFINITY);
  int n=X.size();
  float combined=0.0;
  for(int i=0 ; i<n ; ++i) {
    float x=exp(X[i]);
    if(!isFinite(x)) continue;
    combined+=x;
    //cout<<"adding "<<x<<endl;
    Vector<PhylogenyNode*> path;
    int informantID=taxonIDs[i];
    PhylogenyNode *informantNode=taxa[informantID].getNode();
    tree->getPath(targetNode,informantNode,path);
    int pathLen=path.size();
    float fullDist=tree->distanceBetween(targetNode,informantNode);
    for(int j=1 ; j<pathLen-1 ; ++j) {
      PhylogenyNode *node=path[j];
      float partialDist=tree->distanceBetween(targetNode,node);
      float prorated=x*partialDist/fullDist;
      int ID=static_cast<Taxon*>(node->getDecoration())->getID();
      float &inf=inferred[ID];
      if(!isFinite(inf) || prorated<inf) inf=prorated;
    }
  }

  // Now form the sum of the observables minus the unobservables
  // (to negate the "double counting" that results from the dependency
  // structure when summing the observables)

  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    if(taxon.getNodeType()==INTERNAL_NODE &&
       taxon.getNode()->getParent()!=NULL) {
      //cout<<"subtracting "<<inferred[i]<<endl;
      combined-=inferred[i];
    }	
  }

  return log(combined/tree->getNumLeaves());
}



float Application::combine_phyloLogistic(Vector<float> &x)
{
  INTERNAL_ERROR;
}



float Application::combine_phyloCRF(Vector<float> &x)
{
  INTERNAL_ERROR;
}



float Application::combine_phyloSpan(Vector<float> &X,
				     Vector<int> &taxonIDs)
{
  int targetID=nameToTaxon[targetSpecies];
  PhylogenyNode *targetNode=taxa[targetID].getNode();
  int n=X.size();
  Array1D< Vector<PhylogenyBranch*> > paths(n);
  Map<PhylogenyBranch*,int> counts;
  float sum=0;
  for(int i=0 ; i<n ; ++i) {
    float x=exp(X[i]);
    if(!isFinite(x)) continue;
    int informantID=taxonIDs[i];
    PhylogenyNode *informantNode=taxa[informantID].getNode();
    Vector<PhylogenyBranch*> &path=paths[i];
    tree->getPath(targetNode,informantNode,path);
    int pathLen=path.size();
    for(int j=0 ; j<pathLen ; ++j) {
      PhylogenyBranch *branch=path[j];
      if(!counts.isDefined(branch)) counts[branch]=0;
      ++counts[branch];
      sum+=branch->getLength();
    }
  }
  float combined=0;
  for(int i=0 ; i<n ; ++i) {
    float x=exp(X[i]);
    if(!isFinite(x)) continue;
    Vector<PhylogenyBranch*> &path=paths[i];
    int pathLen=path.size();
    float normLength=0;
    for(int j=0 ; j<pathLen ; ++j) {
      PhylogenyBranch *branch=path[j];
      if(!counts.isDefined(branch)) INTERNAL_ERROR;
      int count=counts[branch];
      normLength+=branch->getLength()/count;
    }
    cout<<"sum="<<sum<<" normLength="<<normLength<<endl;
    normLength/=sum;
    cout<<"    now norm="<<normLength<<endl;
    combined+=x*normLength;
    cout<<"    combined="<<combined<<endl;
  }
  return combined;
}

