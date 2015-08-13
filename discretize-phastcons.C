/****************************************************************
 discretize-phastcons.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Regex.H"
using namespace std;
using namespace BOOM;


class Application
{
  Regex start;
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
  : start("start=(\\\d+)")
  {
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=2)
      throw String("discretize-phastcons <phastcons.pp> <outfile>");
    String infile=cmd.arg(0), outfile=cmd.arg(1);
    
    ofstream os(outfile.c_str());
    File in(infile);
    int pos=0;
    while(!in.eof()) {
      String line=in.getline();
      if(start.search(line)) {
	int newPos=start[1].asInt()-1;
	for(int i=pos ; i<newPos ; ++i) os<<'0';
	pos=newPos;
      }
      else {
	float f=line.asFloat();
	int i=int(f*10);
	if(i>9) i=9;
	++pos;
	os<<i;
      }
    }

    return 0;
  }

