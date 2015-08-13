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
#include "AVES.H"
#include "SparseMatrix3D.H"
#include "SparseMatrix2D.H"
using namespace std;
using namespace BOOM;


class Application : public AVES {
public:
  Application();
  int main(int argc,char *argv[]);
  typedef SparseMatrix3D::EntryList EntryList;
  typedef SparseMatrix3D::Entry Entry;
protected:
  String inDir, outDir;
  String wantFactor;
  int numAlpha, numClasses, numStates;
  Regex filenameRegex, pathRegex;
  BranchHMM *hmm;
  float threshold;
  void sumOutStates();
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
}



int Application::main(int argc,char *argv[])
{
  // Process command line`
  CommandLine cmd(argc,argv,"Gs:");
  if(cmd.numArgs()!=4)
    throw String(
"\n\
marginalize-states [options] <factor-model.lambda> <factor> <in-dir> <out-dir>\n\
   options are:\n\
      -G = no garbage collection\n\
      -s k = use k as sparseness threshold for matrices (default=0.001)\n\
\n\
");
  String lambdaFile=cmd.arg(0);
  wantFactor=cmd.arg(1);
  inDir=cmd.arg(2);
  outDir=cmd.arg(3);
  disableGC=cmd.option('G');
  if(disableGC) modelCompiler->setGCthreshold(LARGEST_INTEGER);
  numAlpha=alphabet.size();
  threshold=cmd.option('s') ? cmd.optParm('s').asFloat() : 0.001;
  
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
  numStates=transducerTemplate->getNumStates();

  // Get phylogeny
  tree=modelCompiler->getGuideTree(transducerTemplate);
  tree->gatherNodes(phylogenyNodes);

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
  numClasses=FunctionalClass::numClasses();

  // Do It.
  sumOutStates();

  return 0;
}



void Application::sumOutStates()
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
    String factor=filenameRegex[1];
    if(!(factor==wantFactor)) continue;
    String firstSpecies=filenameRegex[2], secondSpecies=filenameRegex[3];
    int tax1id=nameToTaxon[firstSpecies], tax2id=nameToTaxon[secondSpecies];
    String outfile=
      outDir+"/"+factor+"-"+firstSpecies+"-"+secondSpecies+".post2";
    SparseMatrix3D &M3=*SparseMatrix3D::loadBinary(inDir+"/"+filename);
    int L1=M3.getFirstDim(), L2=M3.getSecondDim();
    Array3D<float> M2(L1,L2,3);
    M2.setAllTo(LOG_0);
    for(int i=0 ; i<L1 ; ++i) {
      for(int q=0 ; q<numStates ; ++q) {
	PairHMMStateType stateType=hmm->getStateType(q);
	if(stateType>2) continue; // ###
	EntryList &row=M3(i,q);
	EntryList::iterator cur=row.begin(), end=row.end();
	for(; cur!=end ; ++cur) {
	  Entry &e3=*cur;
	  float &e2=M2[i][e3.y][stateType]; // FIXED: 12/9/2010
	  e2=sumLogProbs(e2,e3.value); 
	}
      }
    }
    delete &M3;
    SparseMatrix3D Ms(M2,threshold);
    //cout<<M2<<endl;
    Ms.saveBinary(outfile);
  }
}





