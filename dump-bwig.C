/****************************************************************
 dump-bwig.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
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
  CommandLine cmd(argc,argv,"b:e:");
  if(cmd.numArgs()!=1)
    throw String("\n\
dump-bwig [options] <infile.bwig>\n\
   options:\n\
      -b begin\n\
      -e end\n\
");
  String infile=cmd.arg(0);
  File inFile(infile);
  int begin=cmd.option('b') ? cmd.optParm('b').asInt() : 0;
  int end=cmd.option('e') ? cmd.optParm('e').asInt() : 
    inFile.getSize()/sizeof(float);
  if(begin>0) inFile.seek(begin*sizeof(float));
  for(int pos=begin ; pos<=end ; ++pos) {
    //while(!inFile.eof())
    cout<<inFile.readFloat()<<endl;
  }
  return 0;
}

