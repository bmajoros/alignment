/****************************************************************
 pair-align.C
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaDashAlphabet.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Array2D.H"
#include "BOOM/Random.H"
#include "Viterbi.H"
#include "BackwardAlgorithm.H"
#include "Sampler.H"
using namespace std;
using namespace BOOM;


class Application
{
  DnaDashAlphabet alphabet;
  Symbol gapSymbol;
  bool shouldSample;
  int numSamples;
  void emit(StatePath &path,const PairHMM &hmm,const String &S1,
	    const String &S2,ostream &);
  Array2D<char> *pathToAlignment(StatePath &path,const PairHMM &hmm,
				 const String &S1,const String &S2);
  void doSampling(const Sequence &,const Sequence &,const String &seq1,
		  const String &seq2,const PairHMM &);
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

    gapSymbol=alphabet.lookup('-');
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"s:");
  if(cmd.numArgs()!=3)
    throw String(
"\npair-align [options] <*.phmm> <*.fasta> <*.fasta>\n\
    where -s = sample N alignments\n\
");
  String hmmFilename=cmd.arg(0);
  String fastaFile1=cmd.arg(1);
  String fastaFile2=cmd.arg(2);
  shouldSample=cmd.option('s');
  if(shouldSample) numSamples=cmd.optParm('s');
  randomize();
  
  // Load the HMM
  PairHMM hmm(alphabet,gapSymbol,hmmFilename);
  hmm.convertToLogs();
  
  // Iterate over pairs of sequences from the two files
  FastaReader reader1(fastaFile1,alphabet), reader2(fastaFile2,alphabet);
  String seq1, seq2, def1, def2;
  while(reader1.nextSequence(def1,seq1) &&
	reader2.nextSequence(def2,seq2)) {
    Sequence S1(seq1,alphabet), S2(seq2,alphabet);
    S1.replaceAll(INVALID_SYMBOL,gapSymbol);
    S2.replaceAll(INVALID_SYMBOL,gapSymbol);
    if(shouldSample) { doSampling(S1,S2,seq1,seq2,hmm); continue; }
    Viterbi viterbi(hmm,S1,S2);
    StatePath *path=viterbi.getPath();
    emit(*path,hmm,seq1,seq2,cout);
    delete path;
  }
  return 0;
}



void Application::emit(StatePath &path,const PairHMM &hmm,
		       const String &S1,const String &S2,
		       ostream &os) 
{
  Array2D<char> &alignment=*pathToAlignment(path,hmm,S1,S2);
  int L=alignment.getFirstDim();
  int numRows=L/60;
  if(L%60>0) ++numRows;
  for(int r=0 ; r<numRows ; ++r) {
    int begin=r*60;
    int end=begin+60;
    if(end>=L) end=L-1;
    for(int i=0 ; i<3 ; ++i) {
      for(int j=begin ; j<=end ; ++j)
	os<<alignment[j][i];
      os<<endl;
    }
    os<<endl;
  }
  delete &alignment;
}


Array2D<char> *Application::pathToAlignment(StatePath &path,
					    const PairHMM &hmm,
					    const String &S1,
					    const String &S2)
{
  int L=path.size(), i=0, j=0;
  Array2D<char> &alignment=*new Array2D<char>(L,3);
  for(int k=0 ; k<L ; ++k) {
    STATE q=path[k];
    switch(hmm.getStateType(q))
      {
      case PHMM_MATCH:
	alignment[k][0]=S1[i++];
	alignment[k][1]='|';
	alignment[k][2]=S2[j++];
	break;
      case PHMM_INSERT:
	alignment[k][0]='-';
	alignment[k][1]='v';
	alignment[k][2]=S2[j++];
	break;
      case PHMM_DELETE:
	alignment[k][0]=S1[i++];
	alignment[k][1]='^';
	alignment[k][2]='-';
	break;
      }
  }
  return &alignment;
}



void Application::doSampling(const Sequence &S1,const Sequence &S2,
			     const String &seq1,const String &seq2,
			     const PairHMM &hmm)
{
  for(int i=0 ; i<numSamples ; ++i) {
    BackwardAlgorithm B(hmm,S1,S2);
    Sampler sampler(hmm,B,S1,S2);
    StatePath *path=sampler.samplePath();
    emit(*path,hmm,seq1,seq2,cout);
    delete path;
  }
}


