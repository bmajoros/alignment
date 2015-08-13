/****************************************************************
 sample-wmm.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GSL/DirichletDistribution.H"
#include "BOOM/GSL/Random.H"
#include "BOOM/File.H"
using namespace std;
using namespace BOOM;
using namespace GSL;

class Application
{
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
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("sample-wmm <counts-file.txt> <random-seed>");
  String infile=cmd.arg(0);
  int seed=cmd.arg(1).asInt();
  GSL::Random::seed(seed);

  // Read counts file and sample from Dirichlet
  File f(infile);
  Vector< Array1D<double> > matrix;
  while(!f.eof()) {
    String line=f.getline();
    BOOM::Vector<BOOM::String> &fields=*line.getFields();
    if(fields.size()<4) continue;
    BOOM::Array1D<double> counts(4), probs(4);
    for(int i=0 ; i<4 ; ++i) counts[i]=fields[i].asDouble();
    delete &fields;
    DirichletDistribution d(counts);
    d.generate(probs);
    matrix.push_back(probs);
  }
  int N=matrix.size();

  // Write out sampled matrix
  cout<<"WMM\nPROMOTER\n0 "<<N<<" 4\n"<<N<<" 0 0 +\n";
  for(int i=0 ; i<N ; ++i) {
    Array1D<double> row=matrix[i];
    for(int j=0 ; j<4 ; ++j) cout<<log(row[j])<<(j<3?"\t":"\n");
  }

  return 0;
}

