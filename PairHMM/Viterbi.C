/****************************************************************
 Viterbi.C
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "Viterbi.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;

const log0=NEGATIVE_INFINITY;
const log1=0.0;


Viterbi::Viterbi(const PairHMM &hmm,const Sequence &S1,const Sequence &S2)
  : hmm(hmm), S1(S1), S2(S2),
    V(S1.getLength()+1,S2.getLength()+1,hmm.getNumStates()),
    T(S1.getLength()+1,S2.getLength()+1,hmm.getNumStates()),
    L1(S1.getLength()), L2(S2.getLength()),
    gap(hmm.getGapSymbol()), numStates(hmm.getNumStates())
{
  if(!hmm.isInLogSpace()) throw "HMM not in log space";
  fillMatrix();
}



double Viterbi::operator()(int i,int j,STATE k) const
{
  return V(i,j,k);
}



int Viterbi::getFirstDim() const
{
  return V.getFirstDim();
}



int Viterbi::getSecondDim() const
{
  return V.getSecondDim();
}



int Viterbi::getThirdDim() const
{
  return V.getThirdDim();
}



void Viterbi::fillMatrix()
{
  // Misc initialization:

  Vector<STATE> I, D, M;
  hmm.getStatesOfType(PHMM_INSERT,I);
  hmm.getStatesOfType(PHMM_DELETE,D);
  hmm.getStatesOfType(PHMM_MATCH,M);
  int nI=I.size(), nD=D.size(), nM=M.size();

  // Initialization of DP matrix:

  V.setAllTo(log0);
  T.setAllTo(-1);
  V(0,0,0)=log1;
  for(int i=1 ; i<=L1 ; ++i)
    for(int ik=0 ; ik<nD ; ++ik) { // deletion states
      STATE k=D[ik], bestH=0;
      double Pe=hmm.getEmitP(k,S1[i-1],gap);
      double bestScore=V(i-1,0,0)+hmm.getTransP(0,k);
      for(int ih=0 ; ih<nD ; ++ih) {
	STATE h=D[ih];
	double score=V(i-1,0,h)+hmm.getTransP(h,k);
	if(score>bestScore) { bestH=h; bestScore=score; }
      }
      V(i,0,k)=bestScore+Pe;
      T(i,0,k)=bestH;
    }
  for(int j=1 ; j<=L2 ; ++j) 
    for(int ik=0 ; ik<nI ; ++ik) { // insertion states
      STATE k=I[ik], bestH=0;
      double Pe=hmm.getEmitP(k,gap,S2[j-1]);
      double bestScore=V(0,j-1,0)+hmm.getTransP(0,k);
      for(int ih=0 ; ih<nI ; ++ih) {
	STATE h=I[ih];
	double score=V(0,j-1,h)+hmm.getTransP(h,k);
	if(score>bestScore) { bestH=h; bestScore=score; }
      }
      V(0,j,k)=bestScore+Pe;
      T(0,j,k)=bestH;
    }
  
  // Now for the recursion:

  for(int i=1 ; i<=L1 ; ++i)
    for(int j=1 ; j<=L2 ; ++j) {
      for(int ik=0 ; ik<nM ; ++ik) { // match states
	STATE k=M[ik], bestH=-2;
	double Pe=hmm.getEmitP(k,S1[i-1],S2[j-1]);
	double bestScore=NEGATIVE_INFINITY;
	for(STATE h=0 ; h<numStates ; ++h) {
	  double score=V(i-1,j-1,h)+hmm.getTransP(h,k);
	  if(score>bestScore) { bestH=h; bestScore=score; }
	}
	V(i,j,k)=bestScore+Pe;
	T(i,j,k)=bestH;
      }
      for(int ik=0 ; ik<nD ; ++ik) { // deletion states
	STATE k=D[ik], bestH=-2;
	double Pe=hmm.getEmitP(k,S1[i-1],gap);
	double bestScore=NEGATIVE_INFINITY;
	for(STATE h=0 ; h<numStates ; ++h) {
	  double score=V(i-1,j,h)+hmm.getTransP(h,k);
	  if(score>bestScore) { bestH=h; bestScore=score; }
	}
	V(i,j,k)=bestScore+Pe;
	T(i,j,k)=bestH;
      }
      for(int ik=0 ; ik<nI ; ++ik) { // insertion states
	STATE k=I[ik], bestH=-2;
	double Pe=hmm.getEmitP(k,gap,S2[j-1]);
	double bestScore=NEGATIVE_INFINITY;
	for(STATE h=0 ; h<numStates ; ++h) {
	  double score=V(i,j-1,h)+hmm.getTransP(h,k);
	  if(score>bestScore) { bestH=h; bestScore=score; }
	}
	V(i,j,k)=bestScore+Pe;
	T(i,j,k)=bestH;
      }
    }
}



StatePath *Viterbi::getPath()
{
  // Initialization:

  StatePath *path=new StatePath;
  int bestH=-1;
  double bestScore=NEGATIVE_INFINITY;
  for(STATE h=1 ; h<numStates ; ++h) {
    double score=V(L1,L2,h)+hmm.getTransP(h,0);
    if(score>bestScore) { bestH=h; bestScore=score; }
  }
  cout<<"bestScore at end="<<bestScore<<", L1="<<L1<<", L2="<<L2<<endl;

  // Recursion:

  int i=L1, j=L2;
  STATE t=bestH;
  while(t>0) {
    path->push_back(t);
    STATE newT=T(i,j,t);
    switch(hmm.getStateType(t)) 
      {
      case PHMM_MATCH:  --i; --j; break;
      case PHMM_INSERT:      --j; break;
      case PHMM_DELETE: --i;      break;
      default:  throw "error";
      }
    t=newT;
  }

  // Termination:
  StatePath *rev=path->getReverse();
  delete path;
  return rev;
}
