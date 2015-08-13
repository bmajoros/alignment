/****************************************************************
 matrix-entropy.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/Vector.H"
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



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=1)
      throw String("matrix-entropy <*.wmm>");
    String infile=cmd.arg(0);
    
    // Load matrix
    Vector< Array1D<float> > M;
    File file(infile);
    for(int i=0 ; i<4 ; ++i) String line=file.readLine();
    int nucIndex=0, L=0;
    bool areLogs=false;
    while(!file.eof()) {
      String line=file.readLine();
      Vector<String> &fields=*line.getFields();
      int n=fields.size();
      if(n<4) break;
      Array1D<float> dummy(4);
      M.push_back(dummy);
      for(int i=0 ; i<4 ; ++i) {
	float x=fields[i].asFloat();
	if(x<0.0) areLogs=true;
	M[nucIndex][i]=x;
      }
      ++nucIndex;
      delete &fields;
    }
    L=M.size();
    Array2D<float> matrix(4,L);
    for(int i=0 ; i<4 ; ++i)
      for(int j=0 ; j<L ; ++j) 
	matrix[i][j]=areLogs ? exp(M[j][i]) : M[j][i];
    
    // Produce output
    double totalBits=0.0, maxH=log2(4), maxBits=0.0;
    for(int j=0 ; j<L ; ++j) {
      double H=0.0;
      for(int i=0 ; i<4 ; ++i) {
	double P=matrix[i][j];
	if(P>0)	H-=P*log2(P);
      }
      double IC=maxH-H;
      totalBits+=IC;
      maxBits+=maxH;
    }
    cout<<totalBits<<" bits (max="<<maxBits<<" bits over "<<L
	<<" positions)"<<endl;

    return 0;
  }

