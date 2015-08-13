/****************************************************************
 generate-background-seq.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaWriter.H"
#include "BOOM/Sequence.H"
#include "BOOM/RouletteWheel.H"
#include "PhyLib/RateMatrix.H"
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
    if(cmd.numArgs()!=2)
      throw String("generate-background-seq <rate-matrix> <seq-length>");
    String filename=cmd.arg(0);
    int length=cmd.arg(1).asInt();

    // Load rate matrix
    RateMatrix *M=RateMatrix::load(filename);

    // Get EQ freqs
    SubstitutionMatrix *Pt=M->instantiate(1.0);
    Array1D<double> eqFreqs;
    Pt->getEqFreqs(eqFreqs);

    // Generate sequence
    Sequence S;
    RouletteWheel W;
    for(int i=0 ; i<4 ; ++i) W.addSector(eqFreqs[i]);
    W.doneAddingSectors;
    for(int i=0 ; i<length ; ++i) {
      S+=Symbol(W.spin());
    }

    // Emit
    FastaWriter writer;
    String *s=S.toString(Pt->getAlphabet());
    writer.addToFasta(String(">background"),*s,cout);

    return 0;
  }

