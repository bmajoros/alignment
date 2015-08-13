/****************************************************************
 train-3state.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/FastaReader.H"
#include "BOOM/DnaDashAlphabet.H"
#include "BOOM/Sequence.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array3D.H"
#include "PairHMM.H"
using namespace std;
using namespace BOOM;


enum States {
  MATCH_STATE=1,
  INSERTION_STATE=2,
  DELETION_STATE=3
};


class Application
{
  DnaDashAlphabet alphabet;
  Symbol gapSymbol;
  int nAlpha;
  Array2D<int> transCounts; // fromState x toState
  Array3D<int> emitCounts;  // fromState x symbol1 x symbol2
  void updateCounts(Sequence &,Sequence &);
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
  nAlpha=alphabet.size();
}



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=2)
      throw String("train-3state <alignment-file> <outfile>");
    String infile=cmd.arg(0);
    String outfile=cmd.arg(1);

    // Initialization
    transCounts.resize(4,4);
    emitCounts.resize(4,nAlpha,nAlpha);
    transCounts.setAllTo(0);
    emitCounts.setAllTo(0);

    // Process all alignments (pairs of sequences in fasta file)
    FastaReader reader(infile,alphabet);
    String def1, seq1, def2, seq2;
    while(reader.nextSequence(def1,seq1) &&
	  reader.nextSequence(def2,seq2)) {
      Sequence S1(seq1,alphabet), S2(seq2,alphabet);
      updateCounts(S1,S2);
    }

    // Compute HMM parameters
    PairHMM hmm(alphabet,gapSymbol,4);
    hmm.setStateType(INSERTION_STATE,PHMM_INSERT);
    hmm.setStateType(DELETION_STATE,PHMM_DELETE);
    hmm.setStateType(MATCH_STATE,PHMM_MATCH);
    hmm.setStateType(0,PHMM_OTHER);
    for(STATE i=0 ; i<4 ; ++i) {
      double sum=0;
      for(STATE j=0 ; j<4 ; ++j)
	sum+=transCounts[i][j];
      for(STATE j=0 ; j<4 ; ++j)
	hmm.setTransP(i,j,transCounts[i][j]/sum);
    }
    for(STATE k=0 ; k<4 ; ++k) {
      double sum=0;
      for(Symbol i=0 ; i<nAlpha ; ++i)
	for(Symbol j=0 ; j<nAlpha ; ++j)
	  sum+=emitCounts(k,i,j);
      if(sum==0) sum=1;
      for(Symbol i=0 ; i<nAlpha ; ++i)
	for(Symbol j=0 ; j<nAlpha ; ++j)
	  hmm.setEmitP(k,i,j,emitCounts(k,i,j)/sum);
    }      

    // Save HMM
    hmm.save(outfile);

    return 0;
  }



void Application::updateCounts(Sequence &S1,Sequence &S2)
{
  // First, infer state sequence
  int L=S1.getLength();
  if(S2.getLength()!=L) throw "unequal sequence lengths";
  Array1D<STATE> path(L);
  for(int i=0 ; i<L ; ++i)
    if(S1[i]==gapSymbol) path[i]=INSERTION_STATE;
    else if(S2[i]==gapSymbol) path[i]=DELETION_STATE;
    else path[i]=MATCH_STATE;
  
  // Now count state transitions
  for(int i=0 ; i<L-1 ; ++i)
    ++transCounts[path[i]][path[i+1]];
  ++transCounts[0][path[0]];
  ++transCounts[path[L-1]][0];

  // Make sure there are no ambiguity codes
  S1.replaceAll(INVALID_SYMBOL,gapSymbol);
  S2.replaceAll(INVALID_SYMBOL,gapSymbol);

  // Now count state emissions
  for(int i=0 ; i<L ; ++i)
    ++emitCounts(path[i],S1[i],S2[i]);
}


