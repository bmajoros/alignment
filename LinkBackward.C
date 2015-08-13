/****************************************************************
 LinkBackward.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "LinkBackward.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Constants.H"
#include "BOOM/PureDnaAlphabet.H"
#include "LinkFelsenstein.H"
#include "BranchAttributes.H"
using namespace std;
using namespace BOOM;


const double log0=NEGATIVE_INFINITY;
const double log1=log(1.0);



/****************************************************************
                      LinkBackward methods
 ****************************************************************/

LinkBackward::LinkBackward(const BranchHMM &hmm,Taxon &parent,Taxon &child,
			   const AlphabetMap &alphabetMap,int numTaxa,
			   LinkFelsenstein &F)
  : hmm(hmm), parent(parent), child(child), alphabetMap(alphabetMap),
    alphabet(PureDnaAlphabet::global()), gapSymbol(hmm.getGapSymbol()),
    numTaxa(numTaxa), F(F), numStates(hmm.getNumStates())
{
  branch=child.getBranchToParent();
  upMap=&branch->getUpMap();
  downMap=&branch->getDownMap();
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  nI=QI.size(); nD=QD.size(); nM=QM.size();
  fillMatrix();
}



double LinkBackward::getLikelihood() const
{
  return B(0,0,0);
}



double LinkBackward::operator()(int i,int j,STATE k) const
{
  return B(i,j,k);
}



int LinkBackward::getFirstDim() const
{
  return B.getFirstDim();
}



int LinkBackward::getSecondDim() const
{
  return B.getSecondDim();
}



int LinkBackward::getThirdDim() const
{
  return B.getThirdDim();
}



double LinkBackward::getCachedEmitP(STATE k,int col1,int col2)
{
  return cachedEmitP(k,col1,col2);
}



double LinkBackward::getEmitP(STATE k,int parentResidueIndex,
			      int childResidueIndex)
{
  FunctionalClass fc=hmm.getFunctionalClass(k,PARENT);
  double &cachedValue=cachedEmitP(k,parentResidueIndex,childResidueIndex);
  if(finite(cachedValue)) return cachedValue;
  int tempU, tempD;
  switch(hmm.getStateType(k)) 
    {
    case PHMM_MATCH:
      tempD=(*downMap)[parentResidueIndex];
      tempU=(*upMap)[childResidueIndex];
      (*downMap)[parentResidueIndex]=childResidueIndex;
      (*upMap)[childResidueIndex]=parentResidueIndex;
      cachedValue=F.ancestralLikelihood(parentResidueIndex,parent,fc);
      (*downMap)[parentResidueIndex]=tempD;
      (*upMap)[childResidueIndex]=tempU;
      break;
    case PHMM_INSERT:
      cachedValue=F.logLikelihood(childResidueIndex,child,fc);
      break;
    case PHMM_DELETE:
      cachedValue=F.outsideLikelihood(parentResidueIndex,parent,child,fc);
      break;
    default: throw "error in LinkBackward::getEmitP";
    }
  return cachedValue;
}



void LinkBackward::fillMatrix()
{
  if(!hmm.isInLogSpace()) throw "HMM is not in log space";
  
  // Initialization:
  
  int m=parent.getSeqLen(), n=child.getSeqLen();
  cachedEmitP.resize(numStates,m+1,n+1);
  cachedEmitP.setAllTo(NEGATIVE_INFINITY);
  B.resize(m+1,n+1,numStates);
  B.setAllTo(NEGATIVE_INFINITY);
  for(STATE k=1 ; k<numStates ; ++k) B(m,n,k)=hmm.getTransP(k,0);
  Array1D<double> logProbs(nI);
  for(int j=n-1 ; j>=0 ; --j) {
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
    for(int j=n-1 ; j>=0 ; --j) {
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
	    default: throw "error in LinkBackward::fillMatrix (recurrence)";
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
      default: throw "error in LinkBackward::fillMatrix (termination)";
      }
  }
  B(0,0,0)=sumLogProbs<double>(logProbs);
}




/****************************************************************
                      RootBackward methods
 ****************************************************************/

RootBackward::RootBackward(const BranchHMM &hmm,Taxon &parent,Taxon &child,
			   const AlphabetMap &alphabetMap,int numTaxa,
			   int parentLength)
  : hmm(hmm), parent(parent), child(child), alphabetMap(alphabetMap),
    alphabet(PureDnaAlphabet::global()), gapSymbol(hmm.getGapSymbol()),
    numTaxa(numTaxa), numStates(hmm.getNumStates()),
    parentLength(parentLength)
{
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  QMD=QM;
  QMD.append(QD);
  nI=QI.size(); nD=QD.size(); nM=QM.size(); nMD=QMD.size();
  h_to_ih.resize(numStates);
  for(int i=0 ; i<nI ; ++i) h_to_ih[QI[i]]=i;
  sumW.resize(nMD);
  computeTstar();
  fillMatrix();
}



double RootBackward::getLikelihood() const
{
  return B(0,0);
}



double RootBackward::operator()(int i,STATE k) const
{
  return B(i,k);
}



int RootBackward::getFirstDim() const
{
  return B.getFirstDim();
}



int RootBackward::getSecondDim() const
{
  return B.getSecondDim();
}



void RootBackward::computeTstar()
{
  // First, compute matrix T
  T.resize(nI,nI);
  for(int i=0 ; i<nI ; ++i) {
    STATE qi=QI[i];
    for(int j=0 ; j<nI ; ++j) {
      STATE qj=QI[j];
      T(i,j)=exp(hmm.getTransP(qi,qj));
    }
  }
  //cout<<"T="<<T<<endl;

  // Now compute Tstar
  GSL::Matrix I(nI,nI), I_T(nI,nI);
  I.becomeIdentity();
  I.subtract(T,I_T);
  if(!I_T.invert(Tstar)) 
    throw "Failed to invert I-T in RootBackward::computeTstar()";
  //cout<<"T*="<<Tstar<<endl;

  // Also compute V (trans probs from insertion states to q0)
  V.resize(nI);
  for(int i=0 ; i<nI ; ++i)
    V[i]=exp(hmm.getTransP(QI[i],0));
  //cout<<"V="<<V<<endl;

  // Now compute TstarV = log(Tstar * V)
  Tstar.multiply(V,TstarV);
  //cout<<"T*V="<<TstarV<<endl;
  TstarV.convertToLogs();
  Tstar.convertToLogs();
}



void RootBackward::fillMatrix()
{
  if(!hmm.isInLogSpace()) throw "HMM is not in log space";
  
  // Initialization:

  int m=parentLength;
  B.resize(m+1,numStates);
  B.setAllTo(NEGATIVE_INFINITY); // ### could be made faster
  //Array1D<double> logProbs(nI+1);
  Vector<double> logProbs;
  for(STATE k=1 ; k<numStates ; ++k) {
    logProbs.clear();
    logProbs.push_back(hmm.getTransP(k,0));
    //cout<<"trans("<<k<<",0)="<<hmm.getTransP(k,0)<<endl;
    for(int ih=0 ; ih<nI ; ++ih) {
      logProbs.push_back(hmm.getTransP(k,QI[ih])+TstarV[ih]);
      //cout<<"trans("<<k<<","<<QI[ih]<<")="<<hmm.getTransP(k,QI[ih])<<endl;
      //cout<<"T*V="<<TstarV[ih]<<endl;
    }
    B(m,k)=sumLogProbs<double>(logProbs);
    //cout<<"B("<<m<<","<<k<<") = "<<B(m,k)<<endl;
  }

  // Recurrence:

  GSL::Vector W(nI), TstarW(nI);
  for(int i=m-1 ; i>=0 ; --i) { // iterate over positions
    computeW(W,TstarW,i);
    for(STATE k=1 ; k<numStates ; ++k) {
      logProbs.clear();
      const Array1D<STATE> &succ=hmm.getSuccessors(k);
      const STATE *pH=&succ[0];
      int numSucc=succ.size();
      for(int j=0 ; j<numSucc ; ++j, ++pH) {
	STATE h=*pH;
	double transP=hmm.getTransP(k,h);
	logProbs.push_back(hmm.getStateType(h)==PHMM_INSERT ?
			   safeAdd(transP,TstarW[h_to_ih[h]]) :
			   safeAdd(transP,B(i+1,h)));
      }
      /*
      for(int ih=0 ; ih<nMD ; ++ih) {
	STATE h=QMD[ih];
	logProbs[h]=safeAdd(hmm.getTransP(k,h),B(i+1,h));
      }
      for(int ih=0 ; ih<nI ; ++ih) {
	STATE h=QI[ih];
	logProbs[h]=safeAdd(hmm.getTransP(k,h),TstarW[ih]);
      }
      */
      B(i,k)=sumLogProbs<double>(logProbs);
      //cout<<"B("<<i<<","<<k<<")="<<B(i,k)<<endl;
      //cout<<"array="<<logProbs<<endl;
    } // STATE k
  } // positions

  // Termination:

  computeW(W,TstarW,0);
  logProbs.clear();
  for(int ih=0 ; ih<nMD ; ++ih) {
    STATE h=QMD[ih];
    logProbs.push_back(safeAdd(hmm.getTransP(0,h),B(1,h)));
  }
  for(int ih=0 ; ih<nI ; ++ih) {
    STATE h=QI[ih];
    logProbs.push_back(safeAdd(hmm.getTransP(0,h),TstarW[ih]));
  }
  B(0,0)=sumLogProbs<double>(logProbs);
  //cout<<"0,0 = "<<B(0,0)<<endl;
}



void RootBackward::computeW(GSL::Vector &W,GSL::Vector &TstarW,int i)
{
  for(int j=0 ; j<nI ; ++j) { // iterate over insertion states
    STATE qj=QI[j];
    for(int k=0 ; k<nMD ; ++k) { // iterate over match & delete states
      STATE qk=QMD[k];
      sumW[k]=safeAdd(hmm.getTransP(qj,qk),B(i+1,qk));
    }
    W[j]=sumLogProbs<double>(sumW);
  }
  Tstar.multiplyInLogSpace(W,TstarW);
}



GSL::Vector &RootBackward::getTstarV()
{
  return TstarV;
}



/****************************************************************
                      ChildBackward methods
 ****************************************************************/

ChildBackward::ChildBackward(const BranchHMM &hmm,Taxon &parent,
			       Taxon &child,const AlphabetMap &alphabetMap,
			       int numTaxa)
  : hmm(hmm), parent(parent), child(child), alphabetMap(alphabetMap),
    alphabet(PureDnaAlphabet::global()), gapSymbol(hmm.getGapSymbol()),
    numTaxa(numTaxa), numStates(hmm.getNumStates()), 
    parentFuncParse(parent.getFunctionalParse())
{
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  QMD=QM;
  QMD.append(QD);
  nI=QI.size(); nD=QD.size(); nM=QM.size(); nMD=QMD.size();
  h_to_ih.resize(numStates);
  for(int i=0 ; i<nI ; ++i) h_to_ih[QI[i]]=i;
  sumW.resize(nMD);
  computeTstar();
  fillMatrix();
}



double ChildBackward::getLikelihood() const
{
  return B(0,0);
}



double ChildBackward::operator()(int i,STATE k) const
{
  return B(i,k);
}



int ChildBackward::getFirstDim() const
{
  return B.getFirstDim();
}



int ChildBackward::getSecondDim() const
{
  return B.getSecondDim();
}



void ChildBackward::computeTstar()
{
  // First, compute matrix T
  T.resize(nI,nI);
  for(int i=0 ; i<nI ; ++i) {
    STATE qi=QI[i];
    for(int j=0 ; j<nI ; ++j) {
      STATE qj=QI[j];
      T(i,j)=exp(hmm.getTransP(qi,qj));
    }
  }
  //cout<<"T="<<T<<endl;

  // Now compute Tstar
  GSL::Matrix I(nI,nI), I_T(nI,nI);
  I.becomeIdentity();
  I.subtract(T,I_T);
  if(!I_T.invert(Tstar)) 
    throw "Failed to invert I-T in ChildBackward::computeTstar()";

  // Also compute V (trans probs from insertion states to q0)
  V.resize(nI);
  for(int i=0 ; i<nI ; ++i)
    V[i]=exp(hmm.getTransP(QI[i],0));

  // Now compute TstarV = log(Tstar * V)
  Tstar.multiply(V,TstarV);
  TstarV.convertToLogs();
  Tstar.convertToLogs();
}



void ChildBackward::fillMatrix()
{
  if(!hmm.isInLogSpace()) throw "HMM is not in log space";
  
  // Initialization:

  int m=parent.getSeq().getLength();
  cachedEmitP.resize(numStates,m+1);
  cachedEmitP.setAllTo(NEGATIVE_INFINITY);
  B.resize(m+1,numStates);
  B.setAllTo(NEGATIVE_INFINITY); // ### could be faster
  Array1D<double> logProbs(nI+1);
  for(STATE k=1 ; k<numStates ; ++k) {
    logProbs[nI]=hmm.getTransP(k,0);
    for(int ih=0 ; ih<nI ; ++ih)
      logProbs[ih]=hmm.getTransP(k,QI[ih])+TstarV[ih];
    B(m,k)=sumLogProbs<double>(logProbs);
  }

  // Recurrence:

  GSL::Vector W(nI), TstarW(nI);
  Vector<double> vLogProbs;
  for(int i=m-1 ; i>=0 ; --i) { // iterate over positions
    computeW(W,TstarW,i);
    for(STATE k=1 ; k<numStates ; ++k) {
      vLogProbs.clear();
      const Array1D<STATE> &succ=hmm.getSuccessors(k);
      const STATE *pH=&succ[0];
      int numSucc=succ.size();
      for(int j=0 ; j<numSucc ; ++j, ++pH) {
	STATE h=*pH;
	double transP=hmm.getTransP(k,h);
	if(hmm.getStateType(h)==PHMM_INSERT)
	  vLogProbs.push_back(safeAdd(transP,TstarW[h_to_ih[h]]));
	else
	  vLogProbs.push_back(hmm.getFunctionalClass(h,PARENT)==
			     parentFuncParse[i] ? 
			     safeAdd(transP,getEmitP(h,i),B(i+1,h))
			     : NEGATIVE_INFINITY);
      }
      B(i,k)=sumLogProbs<double>(vLogProbs);
    } // STATE k
  } // positions

  // Termination:

  logProbs.resize(numStates);
  logProbs[0]=NEGATIVE_INFINITY;
  computeW(W,TstarW,0);
  for(int ih=0 ; ih<nMD ; ++ih) {
    STATE h=QMD[ih];
    logProbs[h]=
      hmm.getFunctionalClass(h,PARENT)==parentFuncParse[0] ? 
      safeAdd(hmm.getTransP(0,h),getEmitP(h,0),B(1,h)) : NEGATIVE_INFINITY;
  }
  for(int ih=0 ; ih<nI ; ++ih) {
    STATE h=QI[ih];
    logProbs[h]=safeAdd(hmm.getTransP(0,h),TstarW[ih]);
  }
  B(0,0)=sumLogProbs<double>(logProbs);
}



void ChildBackward::computeW(GSL::Vector &W,GSL::Vector &TstarW,int i)
{
  for(int j=0 ; j<nI ; ++j) { // iterate over insertion states
    STATE qj=QI[j];
    for(int k=0 ; k<nMD ; ++k) { // iterate over match & delete states
      STATE qk=QMD[k];
      sumW[k]=
	hmm.getFunctionalClass(qk,PARENT)==parentFuncParse[i] ?
	safeAdd(hmm.getTransP(qj,qk),B(i+1,qk),getEmitP(qk,i)) :   //###
	NEGATIVE_INFINITY;
      //cout<<"\tk="<<k<<" fc(qk)="<<hmm.getFunctionalClass(qk)<<" par_fc="<<parentFuncParse[i]<<" trans="<<hmm.getTransP(qj,qk)<<" B="<<B(i+1,qk)<<" emit="<<getEmitP(qk,i)<<endl;
    }
    //cout<<"i="<<i<<" fc="<<parentFuncParse[i]<<" sumW="<<sumW<<endl;
    W[j]=sumLogProbs<double>(sumW);
  }
  Tstar.multiplyInLogSpace(W,TstarW);
  //cout<<"Tstar="<<Tstar<<"\nW="<<W<<"\nTstarW="<<TstarW<<endl;
}



GSL::Vector &ChildBackward::getTstarV()
{
  return TstarV;
}



double ChildBackward::getEmitP(STATE k,int parentResidueIndex)
{
  double &cachedValue=cachedEmitP(k,parentResidueIndex);
  if(finite(cachedValue)) return cachedValue;
  if(hmm.getFunctionalClass(k,PARENT)!=parentFuncParse[parentResidueIndex]) {
    cachedValue=NEGATIVE_INFINITY;
    //if(hmm.crossFunctionalType(k)==GLT_RETENTION) throw "aha!";
  }
  else switch(hmm.getStateType(k)) 
    {
    case PHMM_MATCH:  // fall through...
    case PHMM_DELETE:
      {
	Symbol s=parent.getSeq()[parentResidueIndex];
	cachedValue=log(hmm.getEqFreqs(k)[s]);
	// ### NOTE: the line above may cause unexpected behavior for
	//           cross-functional states, since the substitution mixture
	//           model's EQ freqs correspond to the descendent...
      }
      break;
    case PHMM_INSERT: 
      cachedValue=0; // =log(1)
      break;
    default: throw "error in LinkBackward::getEmitP";
    }
  return cachedValue;
}



double ChildBackward::getCachedEmitP(STATE k,int pos)
{
  return cachedEmitP(k,pos);
}


