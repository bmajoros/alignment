/****************************************************************
 library-info.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Regex.H"
#include "BOOM/SummaryStats.H"
#include "SparseMatrix3D.H"
#include "AVES.H"
using namespace std;
using namespace BOOM;

typedef SparseMatrix3D::EntryList EntryList;

class Application : public AVES {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Regex filenameRegex, pathRegex;
  bool shouldDumpMatrices;
  void loadMatrices(const String &dirname);
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
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("\n\
library-info <library-dir> <model.lambda>\n\
");
  String libraryDir=cmd.arg(0);
  String lambdaFile=cmd.arg(1);

  LambdaAPI &lambda=modelCompiler->getLambda();
  lambda.getGC().setSilence(true);
  modelCompiler->parse(lambdaFile);
  lambda.checkGarbageLevel();
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  transducerTemplate=modelCompiler->instantiate(ctor);

  // Get phylogeny
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

  loadMatrices(libraryDir);
  return 0;
}



/****************************************************************
                     Application::loadMatrices()
 ****************************************************************/
void Application::loadMatrices(const String &inDir)
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
    int tax1=nameToTaxon[firstSpecies], tax2=nameToTaxon[secondSpecies];
    float D=tree->distanceBetweenLeaves(tax1,tax2);
    SparseMatrix3D *M=SparseMatrix3D::loadBinary(inDir+"/"+filename);
    int L1=M->getFirstDim();
    int L2=M->getSecondDim();
    int numStates=M->getThirdDim();
    //cout<<"tax1="<<tax1<<"="<<taxa[tax1].getName()<<"="<<firstSpecies<<", tax2="<<tax2<<"="<<taxa[tax2].getName()<<"="<<secondSpecies<<endl;
    cout<<baseName<<": "<<L1<<"x"<<L2<<", #states="<<numStates
	<<" phyletic distance="<<D<<endl;
    Vector<int> sizes;
    for(int x=0 ; x<L1 ; ++x)
      for(int q=0 ; q<numStates ; ++q) {
	SparseMatrix3D::EntryList &row=(*M)(x,q);
	sizes.push_back(row.size());
      }
    SummaryStats stats(sizes);
    cout<<"  ave row size="<<stats.getMean()<<" ("<<stats.getMin()
	<<"-"<<stats.getMax()<<", sd="<<stats.getStdDev()<<" N="
	<<stats.getN()<<" total="<<stats.getSum()<<")"<<endl;
    delete M;
  }
}




