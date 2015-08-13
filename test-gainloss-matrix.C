/****************************************************************
 test-gainloss-matrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BirthDeathMatrix.H"
using namespace std;
using namespace BOOM;


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


void dump(Vector<float> V,float inc)
{
  int n=V.size();
  float t=0;
  for(int i=0 ; i<n ; ++i) {
    cout<<t<<"\t"<<V[i]<<endl;
    t+=inc;
  }
}

int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw String("test-gainloss-matrix <lambda> <mu> <max-time>");
  float lambda=cmd.arg(0).asFloat();
  float mu=cmd.arg(1).asFloat();
  float maxT=cmd.arg(2).asFloat();
  const float DELTA=0.01;

  BirthDeathMatrix Q(lambda,mu);
  float birthProb, deathProb;
  Vector<float> births, deaths;
  for(float t=0 ; t<maxT ; t+=DELTA) {
    SubstitutionMatrix &Pt=*Q.instantiate(t);
    births.push_back(Pt(0,1));
    deaths.push_back(Pt(1,0));
    delete &Pt;
  }
  dump(births,DELTA);
  cout<<endl;
  dump(deaths,DELTA);
  return 0;
}

