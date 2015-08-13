/****************************************************************
 BandedBackward.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BandedBackward.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Constants.H"
#include "BOOM/PureDnaAlphabet.H"
using namespace std;
using namespace BOOM;


const double log0=NEGATIVE_INFINITY;
const double log1=log(1.0);



/****************************************************************
                      BandedBackward methods
 ****************************************************************/

BandedBackward::BandedBackward(const BranchHMM &hmm,Sequence &left,
			     Sequence &right,const AlphabetMap &alphabetMap,
			     BandingPattern &bandingPattern)
  : BandedFB_Base(PureDnaAlphabet::global(),alphabetMap,hmm.getGapSymbol(),
		  hmm,bandingPattern,left,right),
    B(M)
{
  fillMatrix();
}



double BandedBackward::getLikelihood() const
{
  return B(0,0,0);
}



void BandedBackward::fillMatrix()
{
  if(!hmm.isInLogSpace()) throw "HMM is not in log space";
  
  // Initialization:
  
  int m=left.getLength(), n=right.getLength();
  int minY, maxY;
  bandingPattern.getBounds(m,minY,maxY);
  if(maxY>n-1) maxY=n-1; // ### necessary?
  B.resize(m+1,n+1,numStates);
  B.setAllTo(NEGATIVE_INFINITY);
  for(STATE k=1 ; k<numStates ; ++k) B(m,n,k)=hmm.getTransP(k,0);
  Array1D<double> logProbs(nI);
  //for(int j=n-1 ; j>=0 ; --j) {
  for(int j=maxY ; j>=minY ; --j) {
    for(STATE k=1 ; k<numStates ; ++k) {
      for(int ih=0 ; ih<nI ; ++ih) {
	STATE h=QI[ih];
	logProbs[ih]=safeAdd(hmm.getTransP(k,h),getEmitP(h,m,j),B(m,j+1,h));
      }
      B(m,j,k)=sumLogProbs<double>(logProbs);
    }
  }
  logProbs.resize(nD);
  for(int i=m-1 ; i>=0 ; --i) {
    bandingPattern.getBounds(i,minY,maxY);
    if(maxY<n) break;
    for(STATE k=1 ; k<numStates ; ++k) {
      for(int ih=0 ; ih<nD ; ++ih) {
	STATE h=QD[ih];
	logProbs[ih]=safeAdd(hmm.getTransP(k,h),B(i+1,n,h),
			     getEmitP(h,i,n));
      }
      B(i,n,k)=sumLogProbs<double>(logProbs);
    }
  }

  // Recurrence:

  logProbs.resize(numStates);
  logProbs[0]=NEGATIVE_INFINITY;
  for(int i=m-1 ; i>=0 ; --i) {
    bandingPattern.getBounds(i,minY,maxY);
    if(maxY>n-1) maxY=n-1;
    //for(int j=n-1 ; j>=0 ; --j) {
    for(int j=maxY ; j>=minY ; --j) {
      for(STATE k=1 ; k<numStates ; ++k) {
	Vector<double> logProbs;
	const Array1D<STATE> &succ=hmm.getSuccessors(k);
	STATE *p=&succ[0];
	int numSucc=succ.size();
	for(int s=0 ; s<numSucc ; ++s) {
	  STATE h=*p;
	  double trans=hmm.getTransP(k,h);
	  switch(hmm.getStateType(h)) 
	    {
	    case PHMM_MATCH:
	      logProbs.push_back(safeAdd(trans,getEmitP(h,i,j),B(i+1,j+1,h)));
	      break;
	    case PHMM_INSERT:
	      logProbs.push_back(safeAdd(trans,getEmitP(h,m,j),B(i,j+1,h)));
	      break;
	    case PHMM_DELETE:
	      logProbs.push_back(safeAdd(trans,getEmitP(h,i,n),B(i+1,j,h)));
	      break;
	    }
	  ++p;
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
	logProbs[h]=safeAdd(trans,getEmitP(h,0,0),B(1,1,h));
	break;
      case PHMM_INSERT:
	logProbs[h]=safeAdd(trans,getEmitP(h,m,0),B(0,1,h));
	break;
      case PHMM_DELETE:
	logProbs[h]=safeAdd(trans,getEmitP(h,0,n),B(1,0,h));
	break;
      }
  }
  B(0,0,0)=sumLogProbs<double>(logProbs);
}




