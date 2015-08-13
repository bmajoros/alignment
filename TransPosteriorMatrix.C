/****************************************************************
 TransPosteriorMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "TransPosteriorMatrix.H"
using namespace std;
using namespace BOOM;


TransPosteriorMatrix::TransPosteriorMatrix(BandedForward &F,BandedBackward &B,
					   PairHMM &hmm,const Sequence &S1,
					   const Sequence &S2,
					   PosteriorMatrix &statePosteriors)
{
  hmm.enumerateTransitions(transitions);
  compute(F,B,hmm,S1,S2,statePosteriors);
}



TransPosteriorMatrix *TransPosteriorMatrix::loadBinary(const String &filename)
{
  TransPosteriorMatrix *PM=new TransPosteriorMatrix;
  PM->loadFromBinary(filename);
  return PM;
}



const Vector<PHMM_Trans> &TransPosteriorMatrix::getTransitions() const 
{
  return transitions;
}



void TransPosteriorMatrix::compute(BandedForward &F,BandedBackward &B,
				   PairHMM &hmm,const Sequence &S1,
				   const Sequence &S2,
				   PosteriorMatrix &statePosteriors)
{
  int L1=F.getFirstDim(), L2=F.getSecondDim(), numStates=F.getThirdDim();
  int nTrans=transitions.size();
  double LL=F.getLikelihood();
  for(int i=0 ; i<L1 ; ++i) {
    Symbol s1=S1[i];
    for(int j=0 ; j<L2 ; ++j) {
      Symbol s2=S2[j];
      for(int it=0 ; it<nTrans ; ++it) {
	const PHMM_Trans &trans=transitions[it];
	const STATE &h=trans.from, k=trans.to;
	const float transP=trans.prob;
	float statePosterior=statePosteriors(i,j,h);
	if(!isFinite(statePosterior)) continue; // sparseness
	int newI=i, newJ=j;
	hmm.updateColumnsFwd(k,newI,newJ);
	M(i,j,it)=
	  F(i,j,h)                   // Forward
	  +transP                    // transition
	  +hmm.hmm.getEmitP(k,s1,s2) // emission
	  +B(newI,newJ,k)            // Backward
	  -statePosterior;           // denominator
      }
    }
  }
}


