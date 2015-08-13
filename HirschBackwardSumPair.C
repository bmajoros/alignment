/****************************************************************
 HirschBackwardSumPair.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "HirschBackwardSumPair.H"
#include "BOOM/Constants.H"
#include "BOOM/Exceptions.H"
#include "BOOM/SumLogProbs.H"
using namespace std;
using namespace BOOM;

HirschBackwardSumPair::HirschBackwardSumPair(int minX,int maxX,int minY,int maxY,
			       const BandingPattern &bp,BranchHMM &hmm,
			       const Vector<STATE> &QI,int numI,
			       const Vector<STATE> &QD,int numD,
			       const Vector<STATE> &QM,int numM,
			       ViterbiConstraint vc,
			       Taxon *parent,Taxon *child,
			       FunctionalParse *parentFP, 
			       FunctionalParse *childFP,
			       IndexMap *upMap,IndexMap *downMap,
			       Hirschberg *hirschberg,
				     ContentSensor *contentSensor)
  : HirschPass(minX,maxX,minY,maxY,bp,hmm,QI,numI,QD,numD,QM,numM,
	       BACKWARD,vc,NULL,parentFP,childFP,upMap,downMap,hirschberg,
	       contentSensor), 
    parent(parent), child(child), 
    parentSeq(parent->getSeq()), childSeq(child->getSeq())
{
  frame.initWindow(maxX);
}



void HirschBackwardSumPair::run()
{
  WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();

  // FIRST INITIALIZATION STEP -- fill in remainder of thisCol:
  int lowerY=max(minY,thisCol.minY), upperY=min(maxY,thisCol.maxY);
  for(int y=upperY-1 ; y>=lowerY ; --y) {
    Array2D<float>::RowIn2DArray<float> pillar=thisCol[y];
    Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[y+1];
    pillar[0]=LOG_0;
    for(STATE k=1 ; k<numStates ; ++k) {
      float bestP=LOG_0;
      const Array1D<STATE> &succQ=hmm.getSuccessors(k);
      int numSuccQ=succQ.size();
      for(int ih=0 ; ih<numSuccQ ; ++ih) {
	STATE h=succQ[ih];
	if(hmm.getStateType(h)!=PHMM_DELETE) continue;
	float emitP=emission(h,y,IndexMap::UNDEFINED,PHMM_DELETE);
	if(emitP==LOG_0) continue;
	float inductiveP=prevPillar[h];
	if(inductiveP==LOG_0) continue;
	float transP=getTransP(k,h);
	float logP=transP+emitP+inductiveP;
	//if(transP>0 || emitP>0 || inductiveP>0) INTERNAL_ERROR; // ###
	if(logP>bestP) bestP=logP; 
      }
      pillar[k]=bestP;
    }
  }

  for(int j=maxX-1 ; j>=minX ; --j) 
    {
      // SECOND INITIALIZATION STEP -- fill in top pillar of nextCol:
      
      WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();
      int lowerY=max(minY,nextCol.minY), upperY=min(maxY,nextCol.maxY);
      Array2D<float>::RowIn2DArray<float> pillar=nextCol[upperY];
      //if(upperY<thisCol.minY || upperY>thisCol.maxY) INTERNAL_ERROR;//###DEBUG
      Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[upperY];
      bool matchValid=upperY<maxY && upperY<thisCol.maxY;
      pillar.setAllTo(LOG_0);//frame.zeroOut(pillar);//
      for(STATE k=1 ; k<numStates ; ++k) {
	float bestP=LOG_0;
	const Array1D<STATE> &succQ=hmm.getSuccessors(k);
	int numSuccQ=succQ.size();
	for(int ih=0 ; ih<numSuccQ ; ++ih) {
	  STATE h=succQ[ih];
	  float logP, emitP, transP, inductiveP;
	  switch(hmm.getStateType(h)) {
	  case PHMM_INSERT:
	    emitP=emission(h,IndexMap::UNDEFINED,j,PHMM_INSERT);
	    if(emitP==LOG_0) continue;
	    inductiveP=prevPillar[h];
	    if(inductiveP==LOG_0) continue;
	    transP=getTransP(k,h);
	    logP=transP+emitP+inductiveP;
	    break;
	  case PHMM_MATCH:
	    if(matchValid) {
	      emitP=emission(h,upperY,j,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      inductiveP=thisCol[upperY+1][h];
	      if(inductiveP==LOG_0) continue;
	      transP=getTransP(k,h);
	      logP=transP+emitP+inductiveP;
	    }
	    else logP=LOG_0;
	    break;
	  default: logP=LOG_0; break;
	  }
	  if(logP>bestP) bestP=logP; 
	}
	pillar[k]=bestP;
      }
      
    // RECURRENCE -- fill in the rest of nextCol:

      for(int y=upperY-1 ; y>=lowerY ; --y) {
	Array2D<float>::RowIn2DArray<float> pillar=nextCol[y];
	Array2D<float>::RowIn2DArray<float> rightPillar=thisCol[y];
	Array2D<float>::RowIn2DArray<float> upPillar=nextCol[y+1];
	Array2D<float>::RowIn2DArray<float> diagPillar=thisCol[y+1];
	bool rightValid=y>=thisCol.minY, diagValid=y+1>=thisCol.minY;
	//if(rightValid && y>thisCol.maxY) INTERNAL_ERROR; // ### DEBUGGING
	//if(diagValid && y+1>thisCol.maxY) INTERNAL_ERROR; // ### DEBUGGING
	pillar[0]=LOG_0; // ### ?
	for(STATE k=1 ; k<numStates ; ++k) {
	  float bestP=LOG_0;
	  const Array1D<STATE> &succQ=hmm.getSuccessors(k);
	  int numSuccQ=succQ.size();
	  for(int ih=0 ; ih<numSuccQ ; ++ih) {
	    STATE h=succQ[ih];
	    float emitP, inductiveP;
	    switch(hmm.getStateType(h)) {
	    case PHMM_MATCH:
	      emitP=emission(h,y,j,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      if(!diagValid) continue;
	      inductiveP=diagPillar[h];
	      break;
	    case PHMM_INSERT:
	      emitP=emission(h,IndexMap::UNDEFINED,j,PHMM_INSERT);
	      if(emitP==LOG_0) continue;
	      if(!rightValid) continue;
	      inductiveP=rightPillar[h];
	      break;
	    case PHMM_DELETE:
	      emitP=emission(h,y,IndexMap::UNDEFINED,PHMM_DELETE);
	      if(emitP==LOG_0) continue;
	      inductiveP=upPillar[h];
	      break; 
	    default:
	      emitP=inductiveP=LOG_0;
	      break;
	    }
	    if(inductiveP==LOG_0) continue;
	    float transP=getTransP(k,h);
	    float logP=transP+emitP+inductiveP;
	    if(logP>bestP) bestP=logP; 
	  }
	  pillar[k]=bestP;
	}
      }
      if(j>minX) frame.advanceWindow();
    }
}



/****************************************************************
                         runPrescan()
 ****************************************************************/
void HirschBackwardSumPair::runPrescan()
{
  WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();

  // FIRST INITIALIZATION STEP -- fill in remainder of thisCol:
  //cout<<"first initialization"<<endl;
  Vector<float> logProbs;
  Vector<STATE> states;
  int lowerY=max(minY,thisCol.minY), upperY=min(maxY,thisCol.maxY);
  for(int y=upperY-1 ; y>=lowerY ; --y) {
    Array2D<float>::RowIn2DArray<float> pillar=thisCol[y];
    Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[y+1];
    pillar.setAllTo(LOG_0);
    states.clear();
    getPrescanStates(y-1,maxX-1,states);
    Vector<STATE>::iterator cur=states.begin(), end=states.end();
    for(; cur!=end ; ++cur) {
      logProbs.clear();
      STATE k=*cur;
      const Array1D<STATE> &succQ=hmm.getSuccessors(k);
      int numSuccQ=succQ.size();
      for(int ih=0 ; ih<numSuccQ ; ++ih) {
	STATE h=succQ[ih];
	if(hmm.getStateType(h)!=PHMM_DELETE) continue;
	float inductiveP=prevPillar[h];
	if(inductiveP==LOG_0) continue;
	float emitP=emission(h,y,IndexMap::UNDEFINED,PHMM_DELETE);
	if(emitP==LOG_0) continue;
	float transP=hmm.getTransP(k,h);
	float logP=transP+emitP+inductiveP;
	logProbs.push_back(logP);
      }
      pillar[k]=sumLogProbs<float>(logProbs);
    }
  }

  for(int j=maxX-1 ; j>=minX ; --j) 
    {
      // SECOND INITIALIZATION STEP -- fill in top pillar of nextCol:
      //cout<<"second initialization, minX="<<minX<<" j="<<j<<endl;

      WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();
      int lowerY=max(minY,nextCol.minY), upperY=min(maxY,nextCol.maxY);
      Array2D<float>::RowIn2DArray<float> pillar=nextCol[upperY];
      Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[upperY];
      bool matchValid=upperY<maxY && upperY<thisCol.maxY;
      pillar.setAllTo(LOG_0);
      states.clear();
      getPrescanStates(upperY-1,j-1,states);
      Vector<STATE>::iterator cur=states.begin(), end=states.end();
      for(; cur!=end ; ++cur) {
	logProbs.clear();
	STATE k=*cur;
	const Array1D<STATE> &succQ=hmm.getSuccessors(k);
	int numSuccQ=succQ.size();
	for(int ih=0 ; ih<numSuccQ ; ++ih) {
	  STATE h=succQ[ih];
	  float logP, emitP, transP, inductiveP;
	  switch(hmm.getStateType(h)) {
	  case PHMM_INSERT:
	    inductiveP=prevPillar[h];
	    if(inductiveP==LOG_0) continue;
	    emitP=emission(h,IndexMap::UNDEFINED,j,PHMM_INSERT);
	    if(emitP==LOG_0) continue;
	    transP=hmm.getTransP(k,h);
	    logP=transP+emitP+inductiveP;
	    break;
	  case PHMM_MATCH:
	    if(matchValid) {
	      inductiveP=thisCol[upperY+1][h];
	      if(inductiveP==LOG_0) continue;
	      emitP=emission(h,upperY,j,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      transP=hmm.getTransP(k,h);
	      logP=transP+emitP+inductiveP;
	    }
	    else continue;//logP=LOG_0;
	    break;
	  default: continue; //logP=LOG_0; break;
	  }
	  logProbs.push_back(logP);
	}
	pillar[k]=sumLogProbs<float>(logProbs);
      }
      
    // RECURRENCE -- fill in the rest of nextCol:
      //cout<<"recurrence"<<endl;

      for(int y=upperY-1 ; y>=lowerY ; --y) {
	Array2D<float>::RowIn2DArray<float> pillar=nextCol[y];
	Array2D<float>::RowIn2DArray<float> rightPillar=thisCol[y];
	Array2D<float>::RowIn2DArray<float> upPillar=nextCol[y+1];
	Array2D<float>::RowIn2DArray<float> diagPillar=thisCol[y+1];
	bool rightValid=y>=thisCol.minY, diagValid=y+1>=thisCol.minY;
	pillar.setAllTo(LOG_0);
	//pillar.printOn(cout);cout<<endl;throw "OK";
	states.clear();
	getPrescanStates(y-1,j-1,states);
	Vector<STATE>::iterator cur=states.begin(), end=states.end();
	for(; cur!=end ; ++cur) {
	  logProbs.clear();
	  STATE k=*cur;
	  const Array1D<STATE> &succQ=hmm.getSuccessors(k);
	  int numSuccQ=succQ.size();
	  //if(k==0) {cout<<numSuccQ<<" SUCCESSORS OF STATE 0"<<endl;}
	  for(int ih=0 ; ih<numSuccQ ; ++ih) {
	    STATE h=succQ[ih];
	    float emitP, inductiveP;
	    switch(hmm.getStateType(h)) {
	    case PHMM_MATCH:
	      inductiveP=diagValid ? diagPillar[h] : LOG_0;
	      if(inductiveP==LOG_0) continue;
	      emitP=emission(h,y,j,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      break;
	    case PHMM_INSERT:
	      inductiveP=rightValid ? rightPillar[h] : LOG_0;
	      if(inductiveP==LOG_0) continue;
	      emitP=emission(h,IndexMap::UNDEFINED,j,PHMM_INSERT);
	      if(emitP==LOG_0) continue;
	      break;
	    case PHMM_DELETE:
	      inductiveP=upPillar[h];
	      if(inductiveP==LOG_0) continue;
	      emitP=emission(h,y,IndexMap::UNDEFINED,PHMM_DELETE);
	      if(emitP==LOG_0) continue;
	      break; 
	    default:
	      continue;
	    }
	    float transP=hmm.getTransP(k,h);
	    float logP=transP+emitP+inductiveP;
	    logProbs.push_back(logP);
	  }
	  pillar[k]=sumLogProbs<float>(logProbs);
	}
      }

      /* // DEBUGGING
      if(columnIsEmpty(j,lowerY,upperY)) {
	for(int y=upperY-1 ; y>=lowerY ; --y) {
	  states.clear();
	  getPrescanStates(y-1,j-1,states);
	  cout<<"j="<<j<<"y="<<y<<" : ";
	  Vector<STATE>::iterator cur=states.begin(), end=states.end();
	  for(; cur!=end ; ++cur) {
	    STATE q=*cur;
	    cout<<q<<":"<<hmm.getFunctionalClass(q,PARENT)<<" ";
	  }
	  cout<<endl;
	  cout<<"parent FC="<<(*parentFP)[y]<<endl;
	}
	cout<<firstTaxon->getName()<<" vs. "<<secondTaxon->getName()<<endl;
	INTERNAL_ERROR;
      }
      */

      //cout<<"elephant j="<<j<<endl<<frame.getThisCol()<<endl;
      if(j>minX) frame.advanceWindow();
    }
}




double HirschBackwardSumPair::getEmitP(STATE k,int p,int c,PairHMMStateType t)
{
  switch(vc)
    {
    case VC_UNCONSTRAINED: 
      break;
    case VC_PARENT_PARSE:
      if(emitsParent(t) && hmm.getFunctionalClass(k,PARENT)!=(*parentFP)[p])
	return LOG_0 ;
      break;
    case VC_CHILD_PARSE:
      if(emitsChild(t) && hmm.getFunctionalClass(k,CHILD)!=(*childFP)[c])
	return LOG_0;
      break;
    case VC_PARENT_AND_CHILD_PARSE:
      if(emitsParent(t) && hmm.getFunctionalClass(k,PARENT)!=(*parentFP)[p] ||
	 emitsChild(t) && hmm.getFunctionalClass(k,CHILD)!=(*childFP)[c])
	return LOG_0;
      break;
    case VC_INDEL_HISTORY: // ### ?
      if(emitsParent(t) && (*downMap)[p]!=c ||
	 emitsChild(t) && (*upMap)[c]!=p)
	return LOG_0;
      break;
    case VC_INDEL_AND_PARENT_PARSE:
      if(emitsParent(t))
	if((*downMap)[p]!=c || 
	   hmm.getFunctionalClass(k,PARENT)!=(*parentFP)[p]) return LOG_0;
      if(emitsChild(t) && (*upMap)[c]!=p) return LOG_0 ;
      break;
    case VC_INDEL_AND_CHILD_PARSE:
      if(emitsParent(t) && (*downMap)[p]!=c) return LOG_0;
      if(emitsChild(t))
	if((*upMap)[c]!=p ||
	   hmm.getFunctionalClass(k,CHILD)!=(*childFP)[c]) return LOG_0;
      break;
    }
  return emission(k,p,c,t);
}



double HirschBackwardSumPair::emission(STATE k,int parentResidueIndex,
			     int childResidueIndex,
			     PairHMMStateType stateType)
{
  switch(stateType)
    {
    case PHMM_MATCH:{
      //if(parentResidueIndex<0 || childResidueIndex<0) INTERNAL_ERROR; // ###
      FunctionalClass fcParent=hmm.getFunctionalClass(k,PARENT);
      FunctionalClass fcChild=hmm.getFunctionalClass(k,CHILD);
      const Array1D<double> &eqFreqs=hmm.getEqFreqs(k);
      SubstitutionMatrix &Pt=*hmm.getSubstMatrix(k);
      Symbol s=parentSeq[parentResidueIndex];
      double eq=log(eqFreqs[s]);
      return Pt(s,childSeq[childResidueIndex])+eq;
      }
    case PHMM_INSERT:{
      //if(childResidueIndex<0) INTERNAL_ERROR; // ###
      FunctionalClass fc=hmm.getFunctionalClass(k,CHILD);
      Array1D<double> &eqFreqs=fc.getEqFreqs(); // ###
      return log(eqFreqs[childSeq[childResidueIndex]]);
      }
    case PHMM_DELETE:{
      //if(parentResidueIndex<0) INTERNAL_ERROR; // ###
      FunctionalClass fc=hmm.getFunctionalClass(k,PARENT);
      Array1D<double> &eqFreqs=fc.getEqFreqs();
      return log(eqFreqs[parentSeq[parentResidueIndex]]);
      }
    default:
      INTERNAL_ERROR;
    }
}



