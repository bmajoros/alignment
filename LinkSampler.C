/****************************************************************
 LinkSampler.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "LinkSampler.H"
#include "BOOM/RouletteWheel.H"
#include "BOOM/Constants.H"
#include "BOOM/Random.H"
#include "BOOM/Exceptions.H"
#include "BOOM/SumLogProbs.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                       LinkSampler methods
 ****************************************************************/

LinkSampler::LinkSampler(const BranchHMM &hmm,const LinkBackward &B,
			 Taxon &parent, Taxon &child,BranchAttributes *branch)
  : hmm(hmm), B(B), parent(parent), child(child), branch(branch)
{
  // ctor
}



StatePath *LinkSampler::samplePath(double &pathScore,
				   double &oldPathScore)
{
  pathScore=0;
  Symbol gapSymbol=hmm.getGapSymbol();
  StatePath *path=new StatePath(&hmm);
  int L1=parent.getSeqLen(), L2=child.getSeqLen(), 
    numStates=hmm.getNumStates();
  int i=0, j=0;
  STATE prevState=0, k;
  do {
    double r=Random0to1(), cumulative=0.0, P;
    while(r==1) r=Random0to1();
    for(k=1 ; k<numStates ; ++k) {
      double transP=hmm.getTransP(prevState,k);
      P=transP-B(i,j,prevState);
      switch(hmm.getStateType(k)) {
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
      default: INTERNAL_ERROR;
      }
      if(finite(P)) cumulative+=exp(P);
      if(cumulative>r) break; // ###
    }
    if(k==numStates) --k;
    pathScore+=P;
    if(!isFinite(pathScore)) INTERNAL_ERROR;
    path->push_back(k);
    updateCols(k,i,j);
    prevState=k;
  }
  while(i<L1 || j<L2);
  if(!isFinite(pathScore)) 
    {cout<<"hmm.getTransP()="<<hmm.getTransP((*path)[path->size()-1],0)<<endl;
      INTERNAL_ERROR;}
  oldPathScore=computePathScore(*branch->getStatePath());
  return path;
}



double LinkSampler::computePathScore(StatePath &path)
{
  int L1=parent.getSeqLen(), L2=child.getSeqLen();
  int len=path.length();
  STATE prev=0;
  double logP=0;
  int colI=0, colJ=0;
  for(int i=0 ; i<len ; ++i) {
    STATE q=path[i];
    logP+=hmm.getTransP(prev,q);
    switch(hmm.getStateType(q)) 
      {
      case PHMM_MATCH:  logP+=B.getCachedEmitP(q,colI,colJ); break;
      case PHMM_INSERT: logP+=B.getCachedEmitP(q,L1,colJ); break;
      case PHMM_DELETE: logP+=B.getCachedEmitP(q,colI,L2); break;
      }
    updateCols(q,colI,colJ);
    prev=q;
  }
  logP+=hmm.getTransP(prev,0);
  logP-=B.getLikelihood(); // =B(0,0,0)
  return logP;
}



void LinkSampler::updateCols(STATE q,int &i,int &j)
{
  switch(hmm.getStateType(q)) 
    {
    case PHMM_MATCH:  ++i; ++j; break;
    case PHMM_INSERT:      ++j; break;
    case PHMM_DELETE: ++i;      break;
    default: 
      cout<<"state="<<q<<", stateType="<<hmm.getStateType(q)<<endl;
      INTERNAL_ERROR;
    }
}






/****************************************************************
                       RootSampler methods
 ****************************************************************/

RootSampler::RootSampler(const BranchHMM &hmm,const RootBackward &B,
			 Taxon &parent, Taxon &child,BranchAttributes *branch)
  : hmm(hmm), B(B), parent(parent), child(child), branch(branch),
    gapSymbol(hmm.getGapSymbol()), numStates(hmm.getNumStates())
{
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  nI=QI.size(); nD=QD.size(); nM=QM.size();
}



StatePath *RootSampler::samplePath(int L,double &pathScore)
{
  pathScore=0;
  StatePath *path=new StatePath(&hmm);
  int i=0, j=0;
  STATE prevState=0, k;
  do {
    //cout<<"pos "<<i<<endl;
    double b=B(i,prevState);
    double r=Random0to1(), cumulative=NEGATIVE_INFINITY, P=NEGATIVE_INFINITY;
    while(r==1) r=Random0to1();
    const Array1D<STATE> &succ=hmm.getSuccessors(prevState);
    int numSucc=succ.size();
    r=log(r);
    Vector<double> temp;
    for(int ik=0 ; ik<numSucc ; ++ik) {
      k=succ[ik];
      if(k==0) continue;
      double transP=hmm.getTransP(prevState,k);
      P=transP-b;
      switch(hmm.getStateType(k)) {
      case PHMM_MATCH:  // fall through...
      case PHMM_DELETE: P+=B(i+1,k); break;
      case PHMM_INSERT: P+=B(i,k); break;
      default: INTERNAL_ERROR;
      }
      //cout<<"k="<<k<<" transP="<<transP<<" P="<<P<<" B="<<b<<endl;
      if(finite(P)) cumulative=sumLogProbs(cumulative,P);//cumulative+=exp(P);
      if(finite(P)) {temp.push_back(P);cumulative=sumLogProbs(temp);}
      //cout<<"cumulative="<<cumulative<<" r="<<r<<endl;
      if(cumulative>r) break; // DO NOT change to >=
    }
    //if(k==numStates) --k;
    //cout<<"end of loop: k="<<k<<" prev="<<prevState<<" numSucc="<<numSucc<<" P="<<P<<" pathScore="<<pathScore<<" r="<<r<<" cum="<<cumulative<<endl;
    if(ik==numSucc) {
      Vector<double> all, some;
      for(STATE i=0 ; i<numStates ; ++i) {
	double P=hmm.getTransP(prevState,i);
	if(finite(P)) all.push_back(P);
      }
      for(STATE i=0 ; i<numSucc ; ++i) {
	double P=hmm.getTransP(prevState,i);
	if(finite(P)) some.push_back(P);
      }
      double sum=sumLogProbs(all), sumSucc=sumLogProbs(some);
      //cout<<"trans sum="<<sum<<" expSum="<<exp(sum)<<" sumSucc="<<sumSucc<<" expSumSucc="<<exp(sumSucc)<<endl;
      throw "didn't select a state";
    }
    pathScore+=P;
    if(!isFinite(pathScore)) INTERNAL_ERROR; // ### happened on March 18, 2009
    path->push_back(k);
    updateCol(k,i);
    prevState=k;
  }
  while(i<L);

  // Extra step: at end of parent sequence, we have to sample 0 or more
  // insertion states:
  GSL::Vector &TstarV=B.getTstarV();
  while(1) {
    double r=Random0to1(), P;
    while(r==1) r=Random0to1();
    double cumulative=exp(hmm.getTransP(prevState,0)-B(L,prevState));
    if(r<=cumulative) break;
    for(int i=0 ; i<nI ; ++i) {
      k=QI[i];
      P=hmm.getTransP(prevState,k)+TstarV[i]-B(L,prevState);
      if(finite(P)) cumulative+=exp(P);
      if(cumulative>=r) break;
    }
    pathScore+=P; // ###
    if(!isFinite(pathScore)) INTERNAL_ERROR;
    path->push_back(k);
    prevState=k;    
  }	
  if(!isFinite(pathScore)) 
    {cout<<"hmm.getTransP()="<<hmm.getTransP((*path)[path->size()-1],0)<<endl;
      INTERNAL_ERROR;}
  return path;
}



void RootSampler::updateCol(STATE q,int &i)
{
  switch(hmm.getStateType(q)) 
    {
    case PHMM_MATCH:  ++i; break;
    case PHMM_INSERT:      break;
    case PHMM_DELETE: ++i; break;
    default: 
      cout<<"state="<<q<<", stateType="<<hmm.getStateType(q)<<endl;
      INTERNAL_ERROR;
    }
}





/****************************************************************
                       ChildSampler methods
 ****************************************************************/

ChildSampler::ChildSampler(const BranchHMM &hmm,const ChildBackward &B,
			 Taxon &parent, Taxon &child,BranchAttributes *branch)
  : hmm(hmm), B(B), parent(parent), child(child), branch(branch),
    gapSymbol(hmm.getGapSymbol()), numStates(hmm.getNumStates())
{
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  nI=QI.size(); nD=QD.size(); nM=QM.size();
}



StatePath *ChildSampler::samplePath(double &pathScore)
{
  pathScore=0;
  StatePath *path=new StatePath(&hmm);
  int i=0, j=0, L=parent.getSeq().getLength();
  STATE prevState=0, k;
  do {
    double r=Random0to1(), cumulative=0.0, P;
    while(r==1) r=Random0to1();
    for(k=1 ; k<numStates ; ++k) {
      double transP=hmm.getTransP(prevState,k);
      P=transP-B(i,prevState);
      switch(hmm.getStateType(k)) {
      case PHMM_MATCH:  // fall through...
      case PHMM_DELETE: 
	P+=safeAdd(B(i+1,k),B.getCachedEmitP(k,i)); 
	break;
      case PHMM_INSERT: 
	P+=B(i,k); 
	break;
      default: INTERNAL_ERROR;
      }
      if(finite(P)) cumulative+=exp(P);
      if(cumulative>r) break; // DO NOT change to >=
    }
    if(k==numStates) {
      if(cumulative<0.99) {
	cout<<"i="<<i<<" r="<<r<<" fc="<<parent.getFunctionalParse()[i]<<" trans="<<hmm.getTransP(prevState,k)<<" B("<<i<<","<<prevState<<")="<<B(i,prevState)<<" B("<<i+1<<","<<k<<")="<<B(i+1,k)<<" emit="<<B.getCachedEmitP(k,i)<<"L="<<L<<" i="<<i<<" k="<<k<<" P="<<P<<" exp(P)="<<exp(P)<<" cum="<<cumulative<<endl;
	INTERNAL_ERROR;
      }
      --k;
    }
    pathScore+=P;
    if(!isFinite(pathScore)) INTERNAL_ERROR; //throw "bad 6";
    path->push_back(k);
    updateCol(k,i);
    prevState=k;
  }
  while(i<L);

  // Extra step: at end of parent sequence, we have to sample 0 or more
  // insertion states:
  GSL::Vector &TstarV=B.getTstarV();
  while(1) {
    double r=Random0to1(), P;
    while(r==1) r=Random0to1();
    double cumulative=exp(hmm.getTransP(prevState,0)-B(L,prevState));
    if(r<=cumulative) break;
    for(int i=0 ; i<nI ; ++i) {
      k=QI[i];
      P=hmm.getTransP(prevState,k)+TstarV[i]-B(L,prevState);
      if(finite(P)) cumulative+=exp(P);
      if(cumulative>=r) break;
      //cout<<"P="<<P<<" cumulative="<<cumulative<<" r="<<r<<endl;
    }
    pathScore+=P;
    if(!isFinite(pathScore)) INTERNAL_ERROR;
    path->push_back(k);
    prevState=k;    
  }	
  if(!isFinite(pathScore)) 
    {cout<<"hmm.getTransP()="<<hmm.getTransP((*path)[path->size()-1],0)<<endl;
      INTERNAL_ERROR;}
  return path;
}



void ChildSampler::updateCol(STATE q,int &i)
{
  switch(hmm.getStateType(q)) 
    {
    case PHMM_MATCH:  ++i; break;
    case PHMM_INSERT:      break;
    case PHMM_DELETE: ++i; break;
    default: 
      cout<<"state="<<q<<", stateType="<<hmm.getStateType(q)<<endl;
      INTERNAL_ERROR;
    }
}



