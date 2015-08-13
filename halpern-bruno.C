/****************************************************************
 halpern-bruno.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/File.H"
#include "BOOM/PureDnaAlphabet.H"
#include "PhyLib/HalpernBruno.H"
using namespace std;
using namespace BOOM;

typedef Vector< Array1D<double> > WeightMatrix;

class Application
{
  PureDnaAlphabet alphabet;
  bool separateFiles, forceOverwrite;
  int numColumns;
  void loadWMM(String filename,WeightMatrix &);
  void revCompWMM(const WeightMatrix &from,WeightMatrix &to);
  void output(WeightMatrix &wmm,RateMatrix *Q,String filestem);
public:
  Application();
  int main(int argc,char *argv[]);
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
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"fsr:");
  int numParms=cmd.numArgs();
  if(numParms!=3)
    throw String("\n\
halpern-bruno [options] <bg-rate-matrix> <PSSM> <outfile>\n\
     where: bg-rate-matrix = a substitution rate matrix\n\
            PSSM = position-specific scoring matrix (\"weight matrix\")\n\
                   (NOTE: do not include an N column)\n\
            -s = produce separate output files for each position\n\
            -f = force overwrite of output files\n\
            -r <file2> = output reverse-complement matrix into file2\n\
");
  String rateMatrixFile=cmd.arg(0);
  String wmmFile=cmd.arg(1);
  String filestem=cmd.arg(2);
  forceOverwrite=cmd.option('f');
  separateFiles=cmd.option('s');

  // Load weight matrix
  WeightMatrix wmm; // the weight matrix
  loadWMM(wmmFile,wmm);

  // Load background matrix
  RateMatrix *Q=RateMatrix::load(rateMatrixFile);

  // Generate output
  output(wmm,Q,filestem);

  if(cmd.option('r')) {
    WeightMatrix revWMM;
    revCompWMM(wmm,revWMM);
    RateMatrix *rcQ=Q->complement();
    output(revWMM,rcQ,cmd.optParm('r'));
  }
  
  return 0;
}



void Application::output(WeightMatrix &wmm,RateMatrix *Q,String filestem)
{
  // Generate output filenames
  Array1D<String> outfiles(numColumns);
  if(separateFiles)
    for(int i=0 ; i<numColumns ; ++i) {
      outfiles[i]=filestem+"_pos"+i;
      if(File::exists(outfiles[i]) && !forceOverwrite) 
	throw outfiles[i]+" already exists -- delete it or use -f option";
    }

  // Process each site separately
  Array1D<HalpernBruno*> matrices(numColumns);
  for(int col=0 ; col<numColumns ; ++col)
    matrices[col]=new HalpernBruno(*Q,wmm[col]);

  // Save the models
  if(separateFiles)
    for(int col=0 ; col<numColumns ; ++col)
      static_cast<RateMatrix*>(matrices[col])->save(outfiles[col]);
  else {
    ofstream os(filestem.c_str());
    os<<numColumns<<endl;
    for(int col=0 ; col<numColumns ; ++col)
      matrices[col]->save(os);
  }
}



void Application::loadWMM(String filename,WeightMatrix &wmm)
{
  // Check the type of file first
  bool isModelFile=false;
  File test(filename);
  String line=test.getline();
  test.close();
  isModelFile=line=="WMM";
  line.trimWhitespace();

  // Now process the file
  File file(filename);
  if(isModelFile) for(int i=0 ; i<4 ; ++i) file.getline(); // skip four lines
  int col=0;
  Array1D<double> probs(4);
  while(!file.eof()) {
    String line=file.getline();
    line.trimWhitespace();
    Vector<String> &fields=*line.getFields();
    int numFields=fields.size();
    if(numFields<=1) continue;
    if(numFields!=4) throw String("Bad line in weight matrix file: \n")+line;
    for(int i=0 ; i<numFields ; ++i)
      probs[i]=fields[i].asDouble();
    if(isModelFile) {
      double sum=0.0;
      for(int i=0 ; i<4 ; ++i) {
	probs[i]=exp(probs[i]);
	if(probs[i]==0.0) probs[i]=0.001;
	if(probs[i]>1.05) throw "matrix entries should be in log space";
	sum+=probs[i];
      }	
      for(int i=0 ; i<4 ; ++i) probs[i]/=sum;
    }
    wmm.push_back(probs);
    delete &fields;
    ++col;
  }
  numColumns=wmm.size();
}



void Application::revCompWMM(const WeightMatrix &from,WeightMatrix &to)
{
  int L=from.size();
  to.resize(L);
  for(int i=0 ; i<L ; ++i) {
    const Array1D<double> &fromCol=from[i];
    Array1D<double> &toCol=to[L-i-1];
    toCol.resize(4);
    for(Symbol a=0 ; a<4 ; ++a) {
      Symbol aC=alphabet.complement(a);
      toCol[a]=fromCol[aC];
    }
  }
}




