/****************************************************************
 compile-wig.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Regex fixedStep, number;
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
  : fixedStep("fixedStep\\s+chrom=(\\S+)\\s+start=(\\d+)\\s+step=(\\d+)"),
    number("\\s*((\\d|\\.)+)\\s*")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("compile-wig <infile.wig> <outfile.bwig>");
  String infile=cmd.arg(0), outfile=cmd.arg(1);

  File inFile(infile), outFile(outfile,"w");
  String chrom;
  int step=1, pos=0;
  float value=0.0;
  while(!inFile.eof()) {
    String line=inFile.getline();
    if(fixedStep.search(line)) {
      chrom=fixedStep[1];
      int start=fixedStep[2].asInt();
      step=fixedStep[3].asInt();
      value=0.0;
      for(; pos<start ; ++pos) {
	outFile.write(value);
      }
      //cout<<chrom<<" "<<start<<" "<<step<<endl;
    }
    else if(number.match(line)) {
      float value=number[1].asFloat();
      for(int i=0 ; i<step ; ++i) {
	outFile.write(value);
	++pos;
      }
      //cout<<"float="<<value<<endl;
    }
  }

  return 0;
}

