/****************************************************************
 merge-slave-alignments.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/File.H"
#include "BOOM/Exceptions.H"
using namespace std;
using namespace BOOM;

int master(int argc,char *argv[]);

int main(int argc,char *argv[])
  {
    try
      {
	return master(argc,argv);
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
    catch(const RootException &e)
      {
	cerr << "Exception: "<<e.getMessage()<<endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



int master(int argc,char *argv[])
{
  int numSlaves=199; // ###
  String outfile="sampled.aln"; // ###
  int numSamples=20000;
  int skipSamples=50;
  int samplesPerSlave=numSamples/numSlaves;
  File out(outfile,"w");
  int n=numSamples-numSlaves*skipSamples;
  out.write(n);
  for(int i=0 ; i<numSlaves ; ++i) {
    int slaveID=i+1;
    String slaveFile=outfile+"-"+slaveID;
    File in(slaveFile);
    long size=in.getSize();
    for(int j=0 ; j<size ; ++j) {
      unsigned char c=in.readByte();
      out.write(c);
    }
  }
  return 0;
}




