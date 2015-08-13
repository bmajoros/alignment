/****************************************************************
 HirschBackwardSum.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "HirschBackwardSum.H"
#include "BOOM/Constants.H"
#include "BOOM/Exceptions.H"
#include "BOOM/SumLogProbs.H"
using namespace std;
using namespace BOOM;

HirschBackwardSum::HirschBackwardSum(int minX,int maxX,int minY,int maxY,
			       const BandingPattern &bp,BranchHMM &hmm,
			       const Vector<STATE> &QI,int numI,
			       const Vector<STATE> &QD,int numD,
			       const Vector<STATE> &QM,int numM,
			       ViterbiConstraint vc,
			       PrecomputedEmissions *pe,
			       FunctionalParse *parentFP, 
			       FunctionalParse *childFP,
			       IndexMap *upMap,IndexMap *downMap,
			       Hirschberg *hirschberg,
				     ContentSensor *contentSensor)
  : HirschPass(minX,maxX,minY,maxY,bp,hmm,QI,numI,QD,numD,QM,numM,
	       BACKWARD,vc,pe,parentFP,childFP,upMap,downMap,hirschberg,
	       contentSensor)
{
  frame.initWindow(maxX);
}



void HirschBackwardSum::run()
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
	float emitP=getEmitP(h,y,IndexMap::UNDEFINED,PHMM_DELETE);
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
	    emitP=getEmitP(h,IndexMap::UNDEFINED,j,PHMM_INSERT);
	    if(emitP==LOG_0) continue;
	    inductiveP=prevPillar[h];
	    if(inductiveP==LOG_0) continue;
	    transP=getTransP(k,h);
	    logP=transP+emitP+inductiveP;
	    break;
	  case PHMM_MATCH:
	    if(matchValid) {
	      emitP=getEmitP(h,upperY,j,PHMM_MATCH);
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
	      emitP=getEmitP(h,y,j,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      if(!diagValid) continue;
	      inductiveP=diagPillar[h];
	      break;
	    case PHMM_INSERT:
	      emitP=getEmitP(h,IndexMap::UNDEFINED,j,PHMM_INSERT);
	      if(emitP==LOG_0) continue;
	      if(!rightValid) continue;
	      inductiveP=rightPillar[h];
	      break;
	    case PHMM_DELETE:
	      emitP=getEmitP(h,y,IndexMap::UNDEFINED,PHMM_DELETE);
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
void HirschBackwardSum::runPrescan()
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
      float bestP=LOG_0;
      const Array1D<STATE> &succQ=hmm.getSuccessors(k);
      int numSuccQ=succQ.size();
      for(int ih=0 ; ih<numSuccQ ; ++ih) {
	STATE h=succQ[ih];
	if(hmm.getStateType(h)!=PHMM_DELETE) continue;
	float inductiveP=prevPillar[h];
	if(inductiveP==LOG_0) continue;
	float emitP=getEmitP(h,y,IndexMap::UNDEFINED,PHMM_DELETE);
	if(emitP==LOG_0) continue;
	float transP=getTransP(k,h);
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
	float bestP=LOG_0;
	const Array1D<STATE> &succQ=hmm.getSuccessors(k);
	int numSuccQ=succQ.size();
	for(int ih=0 ; ih<numSuccQ ; ++ih) {
	  STATE h=succQ[ih];
	  float logP, emitP, transP, inductiveP;
	  switch(hmm.getStateType(h)) {
	  case PHMM_INSERT:
	    inductiveP=prevPillar[h];
	    if(inductiveP==LOG_0) continue;
	    emitP=getEmitP(h,IndexMap::UNDEFINED,j,PHMM_INSERT);
	    if(emitP==LOG_0) continue;
	    transP=getTransP(k,h);
	    logP=transP+emitP+inductiveP;
	    break;
	  case PHMM_MATCH:
	    if(matchValid) {
	      inductiveP=thisCol[upperY+1][h];
	      if(inductiveP==LOG_0) continue;
	      emitP=getEmitP(h,upperY,j,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      transP=getTransP(k,h);
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
	  float bestP=LOG_0;
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
	      emitP=getEmitP(h,y,j,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      break;
	    case PHMM_INSERT:
	      inductiveP=rightValid ? rightPillar[h] : LOG_0;
	      if(inductiveP==LOG_0) continue;
	      emitP=getEmitP(h,IndexMap::UNDEFINED,j,PHMM_INSERT);
	      if(emitP==LOG_0) continue;
	      break;
	    case PHMM_DELETE:
	      inductiveP=upPillar[h];
	      if(inductiveP==LOG_0) continue;
	      emitP=getEmitP(h,y,IndexMap::UNDEFINED,PHMM_DELETE);
	      if(emitP==LOG_0) continue;
	      break; 
	    default:
	      continue;
	    }
	    float transP=getTransP(k,h);
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

