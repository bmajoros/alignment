/****************************************************************
 find-pileups.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Map.H"
#include "BOOM/Constants.H"
#include "BOOM/Histogram.H"
#include "BOOM/Regex.H"
using namespace std;
using namespace BOOM;

const int NUM_BINS=100;
const float PSEUDOCOUNT=0;
const float DEFAULT_TAIL_AREA=0.01;

class Application {
  int minDiff;
  float tailArea;
  Regex factorRegex;
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
  : factorRegex("-([^-]+).bed")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"f:d:t:");
  if(cmd.numArgs()!=1)
    throw String("\n\
find-pileups [options] <infile.bed>\n\
   options:\n\
      -f # : flag any position having this many reads or more\n\
      -d # : min difference between neighboring counts to be flagged\n\
      -t P : flag using specified P-value\n\
");
  String filename=cmd.arg(0);
  int flag=cmd.option('f') ? cmd.optParm('f').asInt() : 0;
  int minDiff=cmd.option('d') ? cmd.optParm('d').asInt() : 10;
  tailArea=cmd.option('t') ? cmd.optParm('t').asFloat() : DEFAULT_TAIL_AREA;

  if(!factorRegex.search(filename)) throw filename;
  String factor=factorRegex[1];

  Map<String,Map<int,int> > chrPosCount;
  File f(filename);
  long double lengthSum=0, sumN=0;
  while(!f.eof()) {
    String line=f.getline();
    if(f.eof()) break;
    Vector<String> &fields=*line.getFields();
    int n=fields.size();
    if(n>=6) {
      const String &chr=fields[0];
      int begin=fields[1].asInt(), end=fields[2].asInt();
      lengthSum+=abs(end-begin);
      ++sumN;
      Map<int,int> &count=chrPosCount[chr];
      if(!count.isDefined(begin)) {count[begin]=1;}
      else count[begin]++;
    }
    delete &fields;
  }
  int aveLen=lengthSum/sumN;
  
  Set<String> chroms;
  chrPosCount.getKeys(chroms);
  Set<String>::iterator cur=chroms.begin(), end=chroms.end();
  /*
  int min=LARGEST_INTEGER, max=SMALLEST_INTEGER;
  for(; cur!=end ; ++cur) {
    const String &key=*cur;
    Map<int,int> &counts=chrPosCount[key];
    Set<int> keys;
    counts.getKeys(keys);
    Set<int>::iterator cur=keys.begin(), end=keys.end();
    for(; cur!=end ; ++cur) {
      int pos=*cur;
      int count=counts[pos];
      if(count<min) min=count;
      else if(count>max) max=count;
    }
  }
  cout<<"min="<<min<<" max="<<max<<endl;
  Histogram<float> H(min,max,NUM_BINS,PSEUDOCOUNT);
  cur=chroms.begin();
  for(; cur!=end ; ++cur) {
    const String &key=*cur;
    Map<int,int> &counts=chrPosCount[key];
    Set<int> keys;
    counts.getKeys(keys);
    Set<int>::iterator cur=keys.begin(), end=keys.end();
    for(; cur!=end ; ++cur) {
      int pos=*cur;
      int count=counts[pos];
      H.addPoint(count);
    }
  }
  H.normalize();
  cout<<H<<endl;
  exit(0);
  int threshold=(int) H.getRightTailThreshold(tailArea);
  cout<<"threshold="<<threshold<<endl;
  exit(0);
  */
  for(; cur!=end ; ++cur) {
    const String &key=*cur;
    Map<int,int> &counts=chrPosCount[key];
    Set<int> keys;
    counts.getKeys(keys);
    Set<int>::iterator cur=keys.begin(), end=keys.end();
    for(; cur!=end ; ++cur) {
      int pos=*cur;
      int count=counts[pos];
      if(flag>0) {
	if(count>flag) {
	  int diffL=counts.isDefined(pos-1) ? count-counts[pos-1] : count;
	  int diffR=counts.isDefined(pos+1) ? count-counts[pos+1] : count;
	  int diffL2=counts.isDefined(pos-2) ? count-counts[pos-2] : count;
	  int diffR2=counts.isDefined(pos+2) ? count-counts[pos+2] : count;
	  //if(diffL>=minDiff && diffR>=minDiff) 
	  float thresh=count/2.0;
	  if(diffL>=thresh && diffR>=thresh && diffL2>=thresh 
	     && diffR2>=thresh)
	    cout<<key<<"\tpileup\t"<<factor<<"\t"<<pos<<"\t"<<pos+aveLen
		<<"\t0\t+\t0\t"<<"/count="<<count<<";/dL="<<diffL<<";/dR="<<diffR<<endl;
	}
      }
      else cout<<pos<<"\t"<<count<<endl;
    }
  }

  return 0;
}

