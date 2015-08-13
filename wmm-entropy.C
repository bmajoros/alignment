/****************************************************************
 wmm-entropy.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Random.H"
#include "BOOM/Entropy.H"
#include "BOOM/PureDnaAlphabet.H"
#include "EGGS/WMM.H"
using namespace std;
using namespace BOOM;

const int NUM_ITERATIONS=5000;
const float MAX_INF_DIFF=0.4;
BOOM::Alphabet alphabet;

class Application
{
  float infContent(int pos,BOOM::Array2D<float> &);
  void compute(WMM &);
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
    randomize();
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=1)
    throw String("wmm-entropy <infile>");
  String infile=cmd.arg(0);
  GarbageIgnorer GC;

  // Load matrix
  WMM *wmm=dynamic_cast<WMM*>(SignalSensor::load(infile,GC));
  if(!wmm) throw "Error loading matrix file";
  
  // Compute entropy
  compute(*wmm);

  return 0;
}



float Application::infContent(int pos,BOOM::Array2D<float> &M)
{
  Vector<float> row;
  M.getRow(pos,row);
  int n=row.size();
  Vector<double> dRow(n);
  for(int i=0 ; i<n ; ++i) dRow[i]=exp(row[i]);
  return Entropy::informationContent(dRow);
}



void Application::compute(WMM &wmm)
{
  BOOM::Array2D<float> &M=wmm.getMatrix(); // indexed as M[pos][symbol]
  int numPositions=M.getFirstDim();
  float totalInf=0.0;
  for(int pos=0 ; pos<numPositions ; ++pos) {
    float inf=infContent(pos,M);
    totalInf+=inf;
  }
  float infPerPos=totalInf/numPositions;
  cout<<totalInf<<" bits total, "<<infPerPos<<" bits per position"<<endl;
}


