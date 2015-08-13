/****************************************************************
 BackwardAlgorithm.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BackwardAlgorithm.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


const double log0=NEGATIVE_INFINITY;
const double log1=log(1.0);


BackwardAlgorithm::BackwardAlgorithm(const PairHMM &hmm,const Sequence &S1,
				     const Sequence &S2)
  : hmm(hmm), S1(S1), S2(S2)
{
  // ctor

  fillMatrix();
}



double BackwardAlgorithm::getLikelihood() const 
{
  return B(0,0,0);
}



double BackwardAlgorithm::operator()(int i,int j,STATE k) const
{
  return B(i,j,k);
}



int BackwardAlgorithm::getFirstDim() const
{
  return B.getFirstDim();
}



int BackwardAlgorithm::getSecondDim() const
{
  return B.getSecondDim();
}



int BackwardAlgorithm::getThirdDim() const
{
  return B.getThirdDim();
}



void BackwardAlgorithm::fillMatrix()
{
  if(!hmm.isInLogSpace()) throw "HMM is not in log space";
  
  // Initialization:
  
  Symbol gapSymbol=hmm.getGapSymbol();
  Vector<STATE> I, D, M;
  hmm.getStatesOfType(PHMM_INSERT,I);
  hmm.getStatesOfType(PHMM_DELETE,D);
  hmm.getStatesOfType(PHMM_MATCH,M);
  int nI=I.size(), nD=D.size(), nM=M.size();
  int m=S1.getLength(), n=S2.getLength(), numStates=hmm.getNumStates();
  B.resize(m+1,n+1,numStates);
  B.setAllTo(NEGATIVE_INFINITY);
  for(STATE k=1 ; k<numStates ; ++k)
    B(m,n,k)=hmm.getTransP(k,0);
  Array1D<double> logProbs(nI);
  for(int j=n-1 ; j>=0 ; --j) {
    for(STATE k=1 ; k<numStates ; ++k) {
      for(int ih=0 ; ih<nI ; ++ih) {
	STATE h=I[ih];
	logProbs[ih]=
	  hmm.getTransP(k,h)+hmm.getEmitP(h,gapSymbol,S2[j])+B(m,j+1,h);
      }
      B(m,j,k)=sumLogProbs<double>(logProbs);
    }
  }
  logProbs.resize(nD);
  for(int i=m-1 ; i>=0 ; --i) {
    for(STATE k=1 ; k<numStates ; ++k) {
      for(int ih=0 ; ih<nD ; ++ih) {
	STATE h=D[ih];
	logProbs[ih]=
	  hmm.getTransP(k,h)+hmm.getEmitP(h,S1[i],gapSymbol)+B(i+1,n,h);
      }
      B(i,n,k)=sumLogProbs<double>(logProbs);
    }
  }

  // Recurrence:

  logProbs.resize(numStates);
  for(int i=m-1 ; i>=0 ; --i) {
    for(int j=n-1 ; j>=0 ; --j) {
      for(STATE k=1 ; k<numStates ; ++k) {
	logProbs.setAllTo(log0);
	for(STATE h=1 ; h<numStates ; ++h) {
	  double trans=hmm.getTransP(k,h);
	  switch(hmm.getStateType(h)) 
	    {
	    case PHMM_MATCH:
	      logProbs[h]=
		trans+hmm.getEmitP(h,S1[i],S2[j])+B(i+1,j+1,h);
	      break;
	    case PHMM_INSERT:
	      logProbs[h]=
		trans+hmm.getEmitP(h,gapSymbol,S2[j])+B(i,j+1,h);
	      break;
	    case PHMM_DELETE:
	      logProbs[h]=
		trans+hmm.getEmitP(h,S1[i],gapSymbol)+B(i+1,j,h);
	      break;
	    default: throw "error in BackwardAlgorithm.C";
	    }
	}
	B(i,j,k)=sumLogProbs<double>(logProbs);
      }
    }
  }

  // Termination:

  for(STATE h=1 ; h<numStates ; ++h) {
    double trans=hmm.getTransP(0,h);
    switch(hmm.getStateType(h)) 
      {
      case PHMM_MATCH:
	logProbs[h]=
	  trans+hmm.getEmitP(h,S1[0],S2[0])+B(1,1,h);
	break;
      case PHMM_INSERT:
	logProbs[h]=
	  trans+hmm.getEmitP(h,gapSymbol,S2[0])+B(0,1,h);
	break;
      case PHMM_DELETE:
	logProbs[h]=
	  trans+hmm.getEmitP(h,S1[0],gapSymbol)+B(1,0,h);
	break;
      default: throw "error in BackwardAlgorithm.C";
      }
  }
  B(0,0,0)=sumLogProbs<double>(logProbs);
}


