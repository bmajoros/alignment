/****************************************************************
 alternate-maf-elements.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/Vector.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultiAlignment.H"
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
    if(cmd.numArgs()!=6)
      throw String("alternate-maf-elements <in-pos.maf> <in-neg.maf> <out.maf> <out.gff> <max-elements> <skip>");
    String infile1=cmd.arg(0);
    String infile2=cmd.arg(1);
    String outfileMaf=cmd.arg(2);
    String outfileGff=cmd.arg(3);
    int maxElements=cmd.arg(4).asInt();
    int skip=cmd.arg(5).asInt();

    // First, read alignment elements and combine into composite MAF file
    ifstream mafFile1(infile1.c_str()), mafFile2(infile2.c_str());
    if(!mafFile1.good()) throw String("Error opening ")+infile1;
    if(!mafFile2.good()) throw String("Error opening ")+infile2;
    int elementsRead=0;
    Vector<MultiAlignment*> elements;
    for(int i=0 ; i<skip ; ++i) {
      MultiAlignment *a1=MultiAlignment::nextAlignmentFromMAF(mafFile1);
      MultiAlignment *a2=MultiAlignment::nextAlignmentFromMAF(mafFile2);
      delete a1; delete a2;
    }
    while(elementsRead<maxElements) {
      MultiAlignment *a1=MultiAlignment::nextAlignmentFromMAF(mafFile1);
      MultiAlignment *a2=MultiAlignment::nextAlignmentFromMAF(mafFile2);
      if(!a1 || !a2) break;
      a1->ensureAllSameLength(); a2->ensureAllSameLength();
      elements.push_back(a1); elements.push_back(a2);
      elementsRead+=2;
    }
    MultiAlignment *combined=MultiAlignment::combine(elements,false);
    combined->toupper();
    mafFile1.close(); mafFile2.close();
    ofstream outMaf(outfileMaf.c_str());
    outMaf<<*combined<<endl;

    // Now generate the GFF file
    int n=elements.size();
    ofstream outGff(outfileGff.c_str());
    int begin=1;
    for(int i=0 ; i<n ; ++i) {
      MultiAlignment *element=elements[i];
      int len=element->getLength();
      int end=begin+len;
      String type=(i%2 ? "neg" : "pos");
      outGff<<"1\tknown\t"<<type<<"\t"<<begin<<"\t"<<end<<"\t.\t+\t.\n";
      begin=end;
    }

    return 0;
  }

