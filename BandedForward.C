/****************************************************************
 BandedForward.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BandedForward.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Constants.H"
#include "BOOM/PureDnaAlphabet.H"
using namespace std;
using namespace BOOM;




/****************************************************************
                      BandedForward methods
 ****************************************************************/

BandedForward::BandedForward(const BranchHMM &hmm,Sequence &left,
			     Sequence &right,const AlphabetMap &alphabetMap,
			     BandingPattern &bandingPattern)
  : BandedFB_Base(PureDnaAlphabet::global(),alphabetMap,hmm.getGapSymbol(),
		  hmm,bandingPattern,left,right),
    F(M)
{
  fillMatrix();
}



double BandedForward::getLikelihood() const
{
  return likelihood;
}



void BandedForward::fillMatrix()
{
  if(!hmm.isInLogSpace()) throw "HMM is not in log space";
  
  // Initialization:

  int minY, maxY;
  bandingPattern.getBounds(0,minY,maxY);
  if(minY<1) minY=1;
  int m=left.getLength(), n=right.getLength();
  F.resize(m+1,n+1,numStates);
  F.setAllTo(LOG_0); // ### redundant
  F(0,0,0)=LOG_1;
  for(STATE k=1 ; k<numStates ; ++k) F(0,0,k)=LOG_0;
  Array1D<double> logProbs(numStates);
  //for(int j=1 ; j<=n ; ++j) {
    for(int j=minY ; j<=maxY ; ++j) {
    //cout<<"OK1"<<endl;
    for(int ik=0 ; ik<nI ; ++ik) { // insertion states
      STATE k=QI[ik];
      double emitP=getEmitP(k,m,j-1);
      for(STATE h=0 ; h<numStates ; ++h) {
	logProbs[h]=safeAdd(hmm.getTransP(h,k),emitP,F(0,j-1,h));
      }
      F(0,j,k)=sumLogProbs<double>(logProbs);
    }
  }
  for(int i=1 ; i<=m ; ++i) {
    bandingPattern.getBounds(i,minY,maxY);
    if(minY>0) break;
    //cout<<"OK2"<<endl;
    for(int ik=0 ; ik<nD ; ++ik) { // deletion states
      STATE k=QD[ik];
      double emitP=getEmitP(k,i-1,n);
      for(STATE h=0 ; h<numStates ; ++h) {
	logProbs[h]=safeAdd(hmm.getTransP(h,k),emitP,F(i-1,0,h));
      }
      F(i,0,k)=sumLogProbs<double>(logProbs);
    }
  }

  // Recurrence:

  for(int i=1 ; i<=m ; ++i) { 
    bandingPattern.getBounds(0,minY,maxY);
    if(minY<1)  minY=1;
    //for(int j=1 ; j<=n ; ++j) { 
    //cout<<"minY="<<minY<<" maxY="<<maxY<<" i="<<i<<endl;
    for(int j=minY ; j<=maxY ; ++j) { 
      //cout<<"OK3"<<endl;
      F(i,j,0)=LOG_0;
      for(STATE k=1 ; k<numStates ; ++k) {
	Vector<double> logProbs;
	const Array1D<STATE> &pred=hmm.getPredecessors(k);
	STATE *p=&pred[0];
	int numPred=pred.size();
	for(int s=0 ; s<numPred ; ++s) {
	  STATE h=*p;
	  double trans=hmm.getTransP(h,k);
	  switch(hmm.getStateType(k)) 
	    {
	    case PHMM_MATCH:
	      logProbs.push_back(safeAdd(trans,getEmitP(k,i-1,j-1),
					 F(i-1,j-1,h)));
	      break;
	    case PHMM_INSERT:
	      logProbs.push_back(safeAdd(trans,getEmitP(k,m,j-1),
					 F(i,j-1,h)));
	      break;
	    case PHMM_DELETE:
	      logProbs.push_back(safeAdd(trans,getEmitP(k,i-1,n),
					 F(i-1,j,h)));
	      break;
	    }
	  ++p;
	}
	F(i,j,k)=sumLogProbs<double>(logProbs);
      }
    }
  }

  // Termination:
  logProbs[0]=LOG_0;
  for(STATE k=1 ; k<numStates ; ++k)
    logProbs[k]=F(m,n,k)+hmm.getTransP(k,0);
  likelihood=sumLogProbs<double>(logProbs);
}


