/****************************************************************
 scan-wig.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/List.H"
#include "BOOM/PriorityTree.H"
#include "WigBinary.H"
using namespace std;
using namespace BOOM;


class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  DirectComparator<float> cmp;
  PriorityTree<float> wTree;
  List<float> wList;
  void init(WigBinary &,int windowSize);
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
  : wTree(cmp)
{
  wTree.enableDuplicates();
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4)
    throw String("scan-wig <*.bwig> <min> <max> <window-size>");
  String infile=cmd.arg(0);
  float minValue=cmd.arg(1).asFloat();
  float maxValue=cmd.arg(2).asFloat();
  int windowSize=cmd.arg(3).asInt();

  WigBinary wig(infile);
  int L=wig.getLength();
  int lastPos=L-windowSize-1;
  init(wig,windowSize);
  for(int pos=0 ; pos<lastPos ; ++pos) {
    float min=wTree.peekMin(), max=wTree.peekMax();
    if(min<=minValue && max>=maxValue) {
      cout<<pos<<"\t"<<min<<"\t"<<max<<endl;
      //pos+=windowSize;
    }
    float x=wig.read(pos+windowSize);
    wTree.insert(x);
    wList.append(x);
    float f=wList.firstElem();
    wList.removeFirstElem();
    wTree.remove(f);
  }
  
  return 0;
}



void Application::init(WigBinary &wig,int windowSize)
{
  for(int i=0 ; i<windowSize ; ++i) {
    float x=wig.read(i);
    wTree.insert(x);
    wList.append(x);
  }
}


