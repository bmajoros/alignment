/****************************************************************
 LinkForward.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "LinkForward.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Constants.H"
#include "BOOM/PureDnaAlphabet.H"
#include "LinkFelsenstein.H"
using namespace std;
using namespace BOOM;




/****************************************************************
                      LinkForward methods
 ****************************************************************/

LinkForward::LinkForward(const BranchHMM &hmm,Taxon &parent,Taxon &child,
			 const AlphabetMap &alphabetMap,int numTaxa,
			 LinkFelsenstein &fels)
  : hmm(hmm), parent(parent), child(child), alphabetMap(alphabetMap),
    alphabet(PureDnaAlphabet::global()), gapSymbol(hmm.getGapSymbol()),
    numTaxa(numTaxa), felsenstein(fels), numStates(hmm.getNumStates())
{
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  fillMatrix();
}



double LinkForward::getLikelihood() const
{
  return likelihood;
}



double LinkForward::operator()(int i,int j,STATE k) const
{
  return F(i,j,k);
}



int LinkForward::getFirstDim() const
{
  return F.getFirstDim();
}



int LinkForward::getSecondDim() const
{
  return F.getSecondDim();
}



int LinkForward::getThirdDim() const
{
  return F.getThirdDim();
}



double LinkForward::getCachedEmitP(STATE k,int col1,int col2)
{
  return cachedEmitP(k,col1,col2);
}



double LinkForward::getEmitP(STATE k,int parentResidueIndex,
			      int childResidueIndex)
{
  throw "needs to be updated";

  FunctionalClass fc=hmm.getFunctionalClass(k,PARENT);
  double &cachedValue=cachedEmitP(k,parentResidueIndex,childResidueIndex);
  if(finite(cachedValue)) return cachedValue;
  switch(hmm.getStateType(k)) 
    {
    case PHMM_MATCH:
      cachedValue=
	felsenstein.ancestralLikelihood(parentResidueIndex,parent,fc);
      break;
    case PHMM_INSERT:
      cachedValue=felsenstein.logLikelihood(childResidueIndex,child,fc);
      break;
    case PHMM_DELETE:
      cachedValue=
	felsenstein.outsideLikelihood(parentResidueIndex,parent,child,fc);
      break;
    default: throw "error in LinkForward::getEmitP";
    }
  return cachedValue;
}



void LinkForward::fillMatrix()
{
  if(!hmm.isInLogSpace()) throw "HMM is not in log space";
  
  // Initialization:
  
  int nI=QI.size(), nD=QD.size(), nM=QM.size();
  int m=parent.getSeqLen(), n=child.getSeqLen();
  cachedEmitP.resize(numStates,m+1,n+1);
  cachedEmitP.setAllTo(NEGATIVE_INFINITY);
  F.resize(m+1,n+1,numStates);
  F.setAllTo(LOG_0); // ### redundant
  F(0,0,0)=LOG_1;
  for(STATE k=1 ; k<numStates ; ++k) F(0,0,k)=LOG_0;
  Array1D<double> logProbs(numStates);
  for(int j=1 ; j<=n ; ++j) {
    for(int ik=0 ; ik<nI ; ++ik) { // insertion states
      STATE k=QI[ik];
      for(STATE h=0 ; h<numStates ; ++h) {
	logProbs[h]=safeAdd(hmm.getTransP(h,k),getEmitP(k,m,j-1),
			    F(0,j-1,h));
      }
      F(0,j,k)=sumLogProbs<double>(logProbs);
    }
  }
  for(int i=1 ; i<=m ; +i) {
    for(int ik=0 ; ik<nD ; ++ik) { // deletion states
      STATE k=QD[ik];
      for(STATE h=0 ; h<numStates ; ++h) {
	logProbs[h]=safeAdd(hmm.getTransP(h,k),getEmitP(k,i-1,n),
			    F(i-1,0,h));
      }
      F(i,0,k)=sumLogProbs<double>(logProbs);
    }
  }

  // Recurrence:

  for(int i=m-1 ; i>=0 ; --i) {
    for(int j=n-1 ; j>=0 ; --j) {
      F(i,j,0)=LOG_0;
      for(STATE k=1 ; k<numStates ; ++k) {
	Vector<double> logProbs;
	const Array1D<STATE> &pred=hmm.getPredecessors(k);
	STATE *p=&pred[0];
	int numPred=pred.size();
	for(int s=0 ; s<numPred ; ++s) {
	  STATE h=*p;
	  double trans=hmm.getTransP(h,k);
	  switch(hmm.getStateType(h)) 
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
	    default: throw "error in LinkForward::fillMatrix (recurrence)";
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


