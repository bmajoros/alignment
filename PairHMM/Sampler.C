/****************************************************************
 Sampler.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "Sampler.H"
#include "BOOM/RouletteWheel.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


Sampler::Sampler(const PairHMM &hmm,const BackwardAlgorithm &B,
		 const Sequence &S1,const Sequence &S2)
  : hmm(hmm), B(B), S1(S1), S2(S2)
{
  // ctor
}



StatePath *Sampler::samplePath()
{
  Symbol gapSymbol=hmm.getGapSymbol();
  StatePath *path=new StatePath;
  int L1=S1.getLength(), L2=S2.getLength(), numStates=hmm.getNumStates();
  int i=0, j=0;
  STATE prevState=0;
  do {
    RouletteWheel wheel;
    wheel.addSector(0.0);
    for(STATE k=1 ; k<numStates ; ++k) {
      double P=hmm.getTransP(prevState,k)-B(i,j,prevState);
      switch(hmm.getStateType(k))
	{
	case PHMM_MATCH:
	  if(i<L1 && j<L2)
	    P+=B(i+1,j+1,k)+hmm.getEmitP(k,S1[i],S2[j]);
	  else P=NEGATIVE_INFINITY;
	  break;
	case PHMM_INSERT: 
	  if(j<L2)
	    P+=B(i,j+1,k)+hmm.getEmitP(k,gapSymbol,S2[j]);
	  else P=NEGATIVE_INFINITY;
	  break;
	case PHMM_DELETE: 
	  if(i<L1)
	    P+=B(i+1,j,k)+hmm.getEmitP(k,S1[i],gapSymbol);
	  else P=NEGATIVE_INFINITY;
	  break;
	default: throw "error in Sampler.C";
	}
      wheel.addSector(exp(P));
    }
    wheel.doneAddingSectors();
    STATE k=wheel.spin(); // could be made faster
    path->push_back(k);
    switch(hmm.getStateType(k)) 
      {
      case PHMM_MATCH:  ++i; ++j; break;
      case PHMM_INSERT:      ++j; break;
      case PHMM_DELETE: ++i;      break;
      default: throw "error in Sampler.C";
      }
    prevState=k;
    if(i==L1 && j==L2) break;
  }
  while(i<L1 || j<L2);
  return path;
}




/****************************************************************
                    GenerativeSampler methods
 ****************************************************************/

GenerativeSampler::GenerativeSampler(const PairHMM &hmm)
  : hmm(hmm)
{
  if(!hmm.isInLogSpace()) 
    throw String("GenerativeSampler requires HMM in log space");
}



StatePath *GenerativeSampler::samplePath()
{
  double score=0.0;
  STATE oldState=0;
  StatePath *path=new StatePath(&hmm);
  while(true) {
    STATE newState=hmm.sampleNextState(oldState);
    score+=hmm.getTransP(oldState,newState);
    if(newState) {
      path->push_back(newState);
      oldState=newState;
    }
    else break;
  } 
  path->setScore(score);
  return path;
}


