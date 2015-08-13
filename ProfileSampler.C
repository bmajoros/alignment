/****************************************************************
 ProfileSampler.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "ProfileSampler.H"
#include "BOOM/RouletteWheel.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


ProfileSampler::ProfileSampler(const PairHMM &hmm,
			       const ProfileBackward &B,
			       const AlignmentView &S1,
			       const AlignmentView &S2)
  : hmm(hmm), B(B), S1(S1), S2(S2)
{
  // ctor
}



StatePath *ProfileSampler::samplePath(double &pathScore,
					  double &oldPathScore)
{
  pathScore=0;
  Symbol gapSymbol=hmm.getGapSymbol();
  StatePath *path=new StatePath(&hmm);
  int L1=S1.getLength(), L2=S2.getLength(), numStates=hmm.getNumStates();
  int i=0, j=0;
  STATE prevState=0;
  do {
    RouletteWheel wheel;
    Vector<double> wheelScores;
    wheel.addSector(0.0); 
    wheelScores.push_back(NEGATIVE_INFINITY);
    for(STATE k=1 ; k<numStates ; ++k) {
      double transP=hmm.getTransP(prevState,k);
      double P=transP-B(i,j,prevState);
      switch(hmm.getStateType(k))
	{
	case PHMM_MATCH:
	  if(i<L1 && j<L2) {
	    double emitP=B.getCachedEmitP(k,i,j);
	    P+=B(i+1,j+1,k)+emitP;
	  }
	  else P=NEGATIVE_INFINITY; 
	  break;
	case PHMM_INSERT: 
	  if(j<L2) {
	    double emitP=B.getCachedEmitP(k,L1,j);
	    P+=B(i,j+1,k)+emitP;
	  }
	  else P=NEGATIVE_INFINITY;
	  break;
	case PHMM_DELETE: 
	  if(i<L1) {
	    double emitP=B.getCachedEmitP(k,i,L2);
	    P+=B(i+1,j,k)+emitP;
	  }
	  else P=NEGATIVE_INFINITY;
	  break;
	default: throw "error in ProfileSampler.C";
	}
      wheel.addSector(exp(P));
      wheelScores.push_back(P);
    }
    wheel.doneAddingSectors();
    STATE k=wheel.spin(); // could be made faster

    /* ### A FASTER WAY: generate random value (r) first, then while
           assembling the roulette wheel, stop when it sums to >r
     */

    pathScore+=wheelScores[k];
    if(!isFinite(pathScore)) throw "bad";
    path->push_back(k);
    switch(hmm.getStateType(k)) 
      {
      case PHMM_MATCH:  ++i; ++j; break;
      case PHMM_INSERT:      ++j; break;
      case PHMM_DELETE: ++i;      break;
      default: 
	cout<<"k="<<k<<" wheel="<<wheel<<endl;
	throw "error in ProfileSampler.C";
      }
    prevState=k;
    if(i==L1 && j==L2) break;
  }
  while(i<L1 || j<L2);
  pathScore+=hmm.getTransP((*path)[path->size()-1],0);
  if(!isFinite(pathScore)) {cout<<"hmm.getTransP()="<<hmm.getTransP((*path)[path->size()-1],0)<<endl;throw "bad";}
  oldPathScore=computeOldPathScore();
  return path;
}



double ProfileSampler::computeOldPathScore()
{
  Vector<STATE> matchStates, insertionStates, deletionStates;
  hmm.getStatesOfType(PHMM_MATCH,matchStates);
  hmm.getStatesOfType(PHMM_INSERT,insertionStates);
  hmm.getStatesOfType(PHMM_DELETE,deletionStates);
  STATE qMatch=matchStates[0], qInsert=insertionStates[0], 
    qDelete=deletionStates[0];
  double logP=-B.getLikelihood();
  int m=S1.getLength(), n=S2.getLength();
  int i=0, j=0, prevState=0, nextState;
  while(i<m && j<n) {
    int I=(i<m ? S1.mapColumn(i) : LARGEST_INTEGER);
    int J=(j<n ? S2.mapColumn(j) : LARGEST_INTEGER);
    if(I==J && I<LARGEST_INTEGER) { // MATCH
      nextState=qMatch;
      logP+=B.getCachedEmitP(qMatch,i,j);
      ++i; ++j;
    }
    else if(I<J) { // DELETION
      nextState=qDelete;
      logP+=B.getCachedEmitP(qDelete,i,n);
      ++i;
    }
    else { // INSERTION
      nextState=qInsert;
      logP+=B.getCachedEmitP(qInsert,m,j);
      ++j;
    }
    logP+=hmm.getTransP(prevState,nextState);
    prevState=nextState;
  }
  logP+=hmm.getTransP(prevState,0);
  return logP;
}


