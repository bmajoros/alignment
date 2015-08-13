/****************************************************************
 combine-factors.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Regex.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/PowerMean.H"
#include "AVES.H"
#include "SparseMatrix3D.H"
#include "SparseMatrix2D.H"
using namespace std;
using namespace BOOM;

enum CombiningStrategy {
  CS_POWER_MEAN
};

class Application : public AVES {
public:
  Application();
  int main(int argc,char *argv[]);
  typedef SparseMatrix3D::EntryList EntryList;
  typedef SparseMatrix3D::Entry Entry;
protected:
  CombiningStrategy strategy;
  String inDir, outDir;
  int numAlpha;
  Regex filenameRegex, pathRegex;
  float power; // used in Power Mean
  float threshold;
  void combineFactors();
  void combineFactors(Taxon &x,Taxon &y);
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
    filenameRegex("([^\/]+)-([^\/]+)-([^\/\\.]+)\\.post2")
{
}



int Application::main(int argc,char *argv[])
{
  // Process command line`
  CommandLine cmd(argc,argv,"Gp:s:t:");
  if(cmd.numArgs()!=2)
    throw String(
"\n\
combine-factors [options] <in-dir> <out-dir>\n\
   options are:\n\
      -G = no garbage collection\n\
      -p k = use k as power in power mean (default=1.5)\n\
      -s k = use combining strategy k (default=1):\n\
             1 = power mean\n\
      -t k = use threshold k for sparseness (default=0.001)\n\
\n\
");
  inDir=cmd.arg(0);
  outDir=cmd.arg(1);
  disableGC=cmd.option('G');
  if(disableGC) modelCompiler->setGCthreshold(LARGEST_INTEGER);
  numAlpha=alphabet.size();
  power=cmd.option('p') ? cmd.optParm('p').asFloat() : 1.5;
  strategy=(CombiningStrategy) cmd.option('s') ? cmd.optParm('s').asInt() : 1;
  threshold=cmd.option('t') ? cmd.optParm('t').asFloat() : 0.001;
  
  // Do It.
  combineFactors();

  return 0;
}



void Application::combineFactors()
{
  Vector<String> files;
  if(!File::getFileList(inDir,files)) 
    throw String("Can't get dir listing for ")+inDir;
  Vector<String>::iterator curF=files.begin(), endF=files.end();
  //Map<String,Vector<SparseMatrix3D*> > matrices;
  Set<String> keys;
  Map<String,Vector<String> > filenames;
  for(; curF!=endF ; ++curF) {
    const String &filename=*curF;
    String baseName, path;
    if(pathRegex.match(filename)) {
      path=pathRegex[1];
      baseName=pathRegex[2];
    }
    else baseName=filename;
    if(!filenameRegex.match(baseName)) continue;
    String factor=filenameRegex[1];
    String firstSpecies=filenameRegex[2], secondSpecies=filenameRegex[3];
    //SparseMatrix3D *factorModel=SparseMatrix3D::loadBinary(inDir+"/"+filename);
    String key=firstSpecies+"-"+secondSpecies;
    keys.insert(key);
    filenames[key].push_back(filename);
    //matrices[key].push_back(factorModel);
  }
  //Set<String> *keys=matrices.getKeys();
  Set<String>::iterator cur=keys.begin(), end=keys.end();
  for(; cur!=end ; ++cur) {
    int L1, L2;
    Array3D<Vector<float> > VM;
    Array3D<float> M;
    const String &key=*cur;
    String outfile=outDir+"/"+key+".post2";
    //Vector<SparseMatrix3D*> &V=matrices[key];
    //Vector<SparseMatrix3D*> matrices;
    Vector<String> &fnames=filenames[key];
    int numFiles=fnames.size();
    for(int i=0 ; i<numFiles ; ++i) {
      String filename=fnames[i];
      String path=inDir+"/"+filename;
      SparseMatrix3D &F=*SparseMatrix3D::loadBinary(path);
      //Vector<SparseMatrix3D*>::iterator cur=V.begin(), end=V.end();
      //for(; cur!=end ; ++cur) {
      //SparseMatrix3D &F=**cur;
      L1=F.getFirstDim(), L2=F.getSecondDim();
      if(M.getFirstDim()==0) {
	M.resize(L1,L2,3); 
	M.setAllTo(LOG_0);
	VM.resize(L1,L2,3);
      }
      for(PHMM_StateType q=0 ; q<3 ; ++q)
	for(int x=0 ; x<L1 ; ++x) {
	  EntryList &row=F(x,q);
	  EntryList::iterator cur=row.begin(), end=row.end();
	  for(; cur!=end ; ++cur) {
	    Entry &entry=*cur;
	    int y=entry.y;
	    float value=entry.value;
	    VM[x][y][q].push_back(value);
	  }
      }
      delete &F;
    }
    for(int x=0 ; x<L1 ; ++x)
      for(int y=0 ; y<L2 ; ++y) 
	for(PHMM_StateType q=0 ; q<3 ; ++q)
	  M[x][y][q]=PowerMean::compute_log(VM[x][y][q],power);
    SparseMatrix3D S(M,threshold);
    S.saveBinary(outfile);
  }
}










