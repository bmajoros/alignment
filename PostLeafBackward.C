/****************************************************************
 PostLeafBackward.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "PostLeafBackward.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Constants.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BranchAttributes.H"
using namespace std;
using namespace BOOM;


const double log0=NEGATIVE_INFINITY;
const double log1=0.0;
typedef SparseMatrix3D::EntryList EntryList;



/****************************************************************
                      PostLeafBackward methods
 ****************************************************************/

PostLeafBackward::PostLeafBackward(const BranchHMM &hmm,Taxon &parent,
				   Taxon &child,SparseMatrix3D &M,
				   double logPseudocount,float indelCoef)
  : hmm(hmm), parent(parent), child(child), numStates(hmm.getNumStates()), 
    M(M), logPseudocount(logPseudocount), indelCoef(indelCoef)
{
  parentL=parent.getSeqLen(), childL=child.getSeqLen();
  cachedEmitP.resize(numStates,parentL+1,childL+1);
  cachedEmitP.setAllTo(log0);
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  nI=QI.size(); nD=QD.size(); nM=QM.size();
  fillMatrix();
}



double PostLeafBackward::getLikelihood() const
{
  return B(0,0,0);
}



double PostLeafBackward::operator()(int i,int j,STATE k) const
{
  return B(i,j,k);
}



int PostLeafBackward::getFirstDim() const
{
  return B.getFirstDim();
}



int PostLeafBackward::getSecondDim() const
{
  return B.getSecondDim();
}



int PostLeafBackward::getThirdDim() const
{
  return B.getThirdDim();
}



float PostLeafBackward::getCachedEmitP(STATE k,int col1,int col2)
{
  return cachedEmitP(k,col1,col2);
}



double PostLeafBackward::getEmitP(STATE q,int newI,int newJ)
{
  float &cachedValue=cachedEmitP(q,newI,newJ);
  if(finite(cachedValue)) return cachedValue;
  bool justPseudo;
  PHMM_StateType stateType=hmm.getStateType(q);
  Vector<float> sumV;
  EntryList &row=M(newI,stateType);
  EntryList::iterator cur=row.begin(), end=row.end();
  for(; cur!=end ; ++cur) {
    SparseMatrix3D::Entry e=*cur;
    if(e.y<newI) continue;
    if(e.y>newI) break;
    float value=e.value;
    if(stateType!=PHMM_MATCH) value+=log(indelCoef); // ###
    return cachedValue=value;
  }
  return logPseudocount;
}




void PostLeafBackward::fillMatrix()
{
  // Initialization:
  
  int m=parent.getSeqLen(), n=child.getSeqLen();
  B.resize(m+1,n+1,numStates);
  B.setAllTo(NEGATIVE_INFINITY);
  for(STATE k=1 ; k<numStates ; ++k) B(m,n,k)=LOG_1;
  Array1D<double> logProbs(nI);
  for(int j=n-1 ; j>=0 ; --j) {
    for(STATE k=1 ; k<numStates ; ++k) {
      for(int ih=0 ; ih<nI ; ++ih) {
	STATE h=QI[ih];
	logProbs[ih]=safeAdd(getEmitP(h,m,j),B(m,j+1,h));
      }
      B(m,j,k)=sumLogProbs<double>(logProbs);
    }
  }
  logProbs.resize(nD);
  for(int i=m-1 ; i>=0 ; --i) {
    for(STATE k=1 ; k<numStates ; ++k) {
      for(int ih=0 ; ih<nD ; ++ih) {
	STATE h=QD[ih];
	logProbs[ih]=safeAdd(B(i+1,n,h),getEmitP(h,i,n));
      }
      B(i,n,k)=sumLogProbs<double>(logProbs);
    }
  }

  // Recurrence:

  logProbs.resize(numStates);
  logProbs[0]=NEGATIVE_INFINITY;
  for(int i=m-1 ; i>=0 ; --i) {
    for(int j=n-1 ; j>=0 ; --j) {
      for(STATE k=1 ; k<numStates ; ++k) {
	Vector<double> logProbs;
	const Array1D<STATE> &succ=hmm.getSuccessors(k);
	STATE *p=&succ[0];
	int numSucc=succ.size();
	for(int s=0 ; s<numSucc ; ++s) {
	  STATE h=*p;
	  switch(hmm.getStateType(h)) 
	    {
	    case PHMM_MATCH:
	      logProbs.push_back(safeAdd(getEmitP(h,i,j),B(i+1,j+1,h)));
	      break;
	    case PHMM_INSERT:
	      logProbs.push_back(safeAdd(getEmitP(h,m,j),B(i,j+1,h)));
	      break;
	    case PHMM_DELETE:
	      logProbs.push_back(safeAdd(getEmitP(h,i,n),B(i+1,j,h)));
	      break;
	    default: 
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
    switch(hmm.getStateType(h)) 
      {
      case PHMM_MATCH:
	logProbs[h]=safeAdd(getEmitP(h,0,0),B(1,1,h));
	break;
      case PHMM_INSERT:
	logProbs[h]=safeAdd(getEmitP(h,m,0),B(0,1,h));
	break;
      case PHMM_DELETE:
	logProbs[h]=safeAdd(getEmitP(h,0,n),B(1,0,h));
	break;
      default: throw "error in PostLeafBackward::fillMatrix (termination)";
      }
  }
  B(0,0,0)=sumLogProbs<double>(logProbs);
}

