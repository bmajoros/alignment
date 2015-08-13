/****************************************************************
 debug-lengths.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaDashAlphabet.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/Map.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Alphabet.H"
#include "PairHMM/PairHMM.H"
#include "PairHMM/Transducer.H"
#include "PairHMM/Sampler.H"
#include "PhyLib/Phylogeny.H"
#include "ModelCompiler.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;

class App
{
  ModelCompiler *modelCompiler; // parses lambda program -> TransducerTemplate
  TransducerTemplate *transducerTemplate; // transducer with free variable t
  TemplateInstantiator *templateInstantiator; // makes transducer for fixed t
  Phylogeny *tree;
  Vector<PhylogenyBranch> branches;
  DnaDashAlphabet alphabet;
  DropGapMapping alphabetMap;
  PureDnaAlphabet pureDnaAlphabet;
  BOOM::Symbol gapSymbol;
public:
  App();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	App app;
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
    catch(const BOOM::RootException &e)
      {
	cerr<<e.getMessage()<<endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



App::App()
  : alphabetMap(DnaDashAlphabet::global(),PureDnaAlphabet::global())
{
  gapSymbol=alphabet.lookup('-');
  modelCompiler=new ModelCompiler;
  templateInstantiator=modelCompiler->getInstantiator();
}



int App::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2) throw String(
"\n\
debug-lengths [options] <in.lambda> <branch-length>\n\
\n\
");
  String lambdaFile=cmd.arg(0);
  double branchLength=cmd.arg(1).asDouble();
  modelCompiler->parse(lambdaFile);
  if(modelCompiler->getNumModels()<1)
    throw "use (register-transducer ...) to register a transducer";
  transducerTemplate=modelCompiler->getIthModel(0);
  BranchHMM *hmm=
    templateInstantiator->instantiate(transducerTemplate,branchLength,&cout);
  GenerativeSampler sampler(*hmm);
  for(int i=0 ; i<1000 ; ++i) {
    StatePath *path=sampler.samplePath();
    int parentLen, childLen;
    path->getSeqLengths(parentLen,childLen);
    cout<<path->length()<<"\t"<<parentLen<<"\t"<<childLen<<endl;
    delete path;
  }
  return 0;
}



