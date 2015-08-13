/****************************************************************
 ConstrainedBackward.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "ConstrainedBackward.H"
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
                      ConstrainedBackward methods
 ****************************************************************/

ConstrainedBackward::ConstrainedBackward(const BranchHMM &hmm,Taxon &parent,
					 Taxon &child,
					 const AlphabetMap &alphabetMap,
					 int numTaxa,LinkFelsenstein &F,
					 const FunctionalParse &parentParse)
  : hmm(hmm), parent(parent), child(child), alphabetMap(alphabetMap),
    alphabet(PureDnaAlphabet::global()), gapSymbol(hmm.getGapSymbol()),
    numTaxa(numTaxa), F(F), numStates(hmm.getNumStates()),
    parentParse(parentParse)
{
  branch=child.getBranchToParent();
  upMap=&branch->getUpMap();
  downMap=&branch->getDownMap();
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  nI=QI.size(); nD=QD.size(); nM=QM.size();
  path=branch->inferStatePath();
  //cout<<"PATH "<<parent.getName()<<" "<<child.getName()<<" "<<*path<<endl;
  pathLen=path->size();
  constrainParent=
    !(parent.getBranchToParent()==NULL && child.whichChild()==LEFT);
  fillMatrix();
}



Vector<PHMM_StateType> *ConstrainedBackward::getPath()
{
  return path;
}



double ConstrainedBackward::getLikelihood() const
{
  return B(0,0);
}



double ConstrainedBackward::operator()(int i,STATE k) const
{
  return B(i,k);
}



int ConstrainedBackward::getFirstDim() const
{
  return B.getFirstDim();
}



int ConstrainedBackward::getSecondDim() const
{
  return B.getSecondDim();
}



double ConstrainedBackward::getCachedEmitP(STATE k,int col)
{
  return cachedEmitP(k,col);
}



double ConstrainedBackward::getEmitP(STATE k,int  parentResidueIndex,
				     int childResidueIndex,int col)
{
  FunctionalClass fc=hmm.getFunctionalClass(k,PARENT);
  double &cachedValue=cachedEmitP(k,col);
  if(finite(cachedValue)) return cachedValue;
  int tempU, tempD;
  switch(hmm.getStateType(k)) 
    {
    case PHMM_MATCH:
      tempD=(*downMap)[parentResidueIndex];
      tempU=(*upMap)[childResidueIndex];
      (*downMap)[parentResidueIndex]=childResidueIndex;
      (*upMap)[childResidueIndex]=parentResidueIndex;
      cachedValue=F.ancestralLikelihood(parentResidueIndex,parent,fc,false);
      (*downMap)[parentResidueIndex]=tempD;
      (*upMap)[childResidueIndex]=tempU;
      break;
    case PHMM_INSERT:
      cachedValue=F.logLikelihood(childResidueIndex,child,fc,false);
      break;
    case PHMM_DELETE:
      cachedValue=F.outsideLikelihood(parentResidueIndex,parent,child,fc,false);
      break;
    default: throw "error in ConstrainedBackward::getEmitP";
    }
  return cachedValue;
}



bool ConstrainedBackward::ok()
{
  return !aborted;
}



ostream &operator<<(ostream &os,const ConstrainedBackward &b)
{
  b.printOn(os);
  return os;
}



void ConstrainedBackward::printOn(ostream &os) const
{
  os<<B<<endl;
}



void ConstrainedBackward::fillMatrix()
{
  if(!hmm.isInLogSpace()) throw "HMM is not in log space";
  Vector<PHMM_StateType> &path=*this->path;
  
  // Initialization:
  
  int m=parent.getSeqLen(), n=child.getSeqLen();
  cachedEmitP.resize(numStates,pathLen+1);
  cachedEmitP.setAllTo(NEGATIVE_INFINITY);
  B.resize(pathLen+1,numStates);
  B.setAllTo(NEGATIVE_INFINITY);
  for(STATE k=1 ; k<numStates ; ++k) {
    if(constrainParent && hmm.getFunctionalClass(k,PARENT)!=parentParse[m-1])
      continue;
    B(pathLen,k)=hmm.getTransP(k,0);
  }

  // Recurrence:

  int i=m-1, j=n-1;
  for(int c=pathLen-1 ; c>0 ; --c) {
    PHMM_StateType t=path[c-1], tNext=path[c];
    Vector<STATE> &Q=(t==PHMM_MATCH ? QM : (t==PHMM_INSERT ? QI : QD));
    int nQ=Q.size();
    bool found=false;
    for(int iq=0 ; iq<nQ ; ++iq) {
      STATE k=Q[iq];
      Vector<double> logProbs;
      const Array1D<STATE> &succ=hmm.getSuccessors(k);
      STATE *p=&succ[0];
      int numSucc=succ.size();
      for(int s=0 ; s<numSucc ; ++s, ++p) {
	STATE h=*p;
	PHMM_StateType hType=hmm.getStateType(h);
	if(hType!=tNext) continue;
	if(constrainParent && emitsParent(hType) &&
	   hmm.getFunctionalClass(h,PARENT)!=parentParse[i]) continue;
	double trans=hmm.getTransP(k,h);
	switch(hType) 
	  {
	  case PHMM_MATCH:
	    logProbs.push_back(safeAdd(trans,getEmitP(h,i,j,c),B(c+1,h)));
	    break;
	  case PHMM_INSERT:
	    logProbs.push_back(safeAdd(trans,getEmitP(h,m,j,c),B(c+1,h)));
	    break;
	  case PHMM_DELETE:
	    logProbs.push_back(safeAdd(trans,getEmitP(h,i,n,c),B(c+1,h)));
	    break;
	  default: break;
	  }
      }
      B(c,k)=sumLogProbs<double>(logProbs);
      if(isFinite(B(c,k))) found=true;
    }
    if(!found) {
      cout<<"ABORT "<<parent.getName()<<":"<<child.getName()<<" "<<path<<endl;
      aborted=true;
      //return;

      cout<<"c="<<c<<" tNext="<<tNext<<" t="<<t<<" m="<<m<<" n="<<n<<" pathLen="<<pathLen<<" constrain="<<constrainParent<<" class="<<parentParse[i]<<endl;
      cout<<"parentParse="<<parent.getFunctionalParse()<<Endl;
      for(int iq=0 ; iq<nQ ; ++iq) {
	STATE k=Q[iq];
	const Array1D<STATE> &succ=hmm.getSuccessors(k);
	STATE *p=&succ[0];
	int numSucc=succ.size();
	cout<<"state "<<k<<" has "<<numSucc<<" successors; ktype="<<hmm.getStateType(k)<<" kclass="<<hmm.getFunctionalClass(k,PARENT)<<"/"<<hmm.getFunctionalClass(k,CHILD)<<endl;
	for(int s=0 ; s<numSucc ; ++s, ++p) {
	  STATE h=*p;
	  PHMM_StateType hType=hmm.getStateType(h);
	  //if(hType!=tNext) continue;
	  //if(constrainParent && emitsParent(hType) && 
	  //  hmm.getFunctionalClass(h,PARENT)!=parentParse[i]) continue;
	  double trans=hmm.getTransP(k,h);
	  cout<<"h="<<h<<" trans="<<trans<<" hType="<<hType<<" class="<<hmm.getFunctionalClass(h,PARENT)<<endl;
	  switch(hType) 
	    {
	    case PHMM_MATCH:
	      cout<<"emit("<<h<<","<<i<<","<<j<<","<<c<<")="<<getEmitP(h,i,j,c)<<" B="<<B(c+1,h)<<endl;
	      break;
	    case PHMM_INSERT:
	      cout<<"emit("<<h<<","<<m<<","<<j<<","<<c<<")="<<getEmitP(h,m,j,c)<<" B="<<B(c+1,h)<<endl;
	      break;
	    case PHMM_DELETE:
	      cout<<"emit("<<h<<","<<i<<","<<n<<","<<c<<")="<<getEmitP(h,i,n,c)<<" B="<<B(c+1,h)<<endl;
	      break;
	    default: 
	      break;
	    }
	}
      }
      INTERNAL_ERROR;
    }
    updateCoordsRev(tNext,i,j);
  }

  // Termination:

  Array1D<double> logProbs(numStates);
  logProbs[0]=NEGATIVE_INFINITY;
  for(STATE h=1 ; h<numStates ; ++h) {
    double trans=hmm.getTransP(0,h);
    switch(hmm.getStateType(h)) 
      {
      case PHMM_MATCH:
	logProbs[h]=safeAdd(trans,getEmitP(h,0,0,0),B(1,h));
	break;
      case PHMM_INSERT:
	logProbs[h]=safeAdd(trans,getEmitP(h,m,0,0),B(1,h));
	break;
      case PHMM_DELETE:
	logProbs[h]=safeAdd(trans,getEmitP(h,0,n,0),B(1,h));
	break;
      default: throw "error in ConstrainedBackward::fillMatrix (termination)";
      }
  }
  B(0,0)=sumLogProbs<double>(logProbs);
  if(!isFinite(B(0,0))) {
    cout<<"XXX"<<endl;
    cout<<parent.getFunctionalParse()<<endl;
    for(STATE h=1 ; h<numStates ; ++h) {
      double trans=hmm.getTransP(0,h);
      cout<<h<<" "<<trans<<" "<<B(1,h)<<" ";
      switch(hmm.getStateType(h)) {
      case PHMM_MATCH: cout<<getEmitP(h,0,0,0)<<endl; break;
      case PHMM_INSERT: cout<<getEmitP(h,m,0,0)<<endl; break;
      case PHMM_DELETE: cout<<getEmitP(h,0,n,0)<<endl; break;
      }
    }
    cout<<path<<endl;
    cout<<B<<endl;
    INTERNAL_ERROR;
  }
  //cout<<"B(0,0)="<<B(0,0)<<endl;
  aborted=false;
}


