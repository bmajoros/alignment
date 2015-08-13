/****************************************************************
 HirschBackwardMax.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "HirschBackwardMax.H"
#include "BOOM/Constants.H"
#include "BOOM/Exceptions.H"
using namespace std;
using namespace BOOM;

HirschBackwardMax::HirschBackwardMax(int minX,int maxX,int minY,int maxY,
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



void HirschBackwardMax::run()
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
void HirschBackwardMax::runPrescan()
{
  WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();

//###
/*
FunctionalParse &fp=hirschberg->getFirstTaxon()->getFunctionalParse();
int fpLen=fp.size();
cout<<"YYY ";
for(int i=0 ; i<fpLen ; ++i) cout<<i<<"="<<fp[i]<<" ";
cout<<endl;
*/
//###

  // FIRST INITIALIZATION STEP -- fill in remainder of thisCol:
  //cout<<"first initialization"<<endl;
  Vector<STATE> states;
  int lowerY=max(minY,thisCol.minY), upperY=min(maxY,thisCol.maxY);
  for(int y=upperY-1 ; y>=lowerY ; --y) {
    Array2D<float>::RowIn2DArray<float> pillar=thisCol[y];
    Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[y+1];
    pillar.setAllTo(LOG_0);//frame.zeroOut(pillar);
    states.clear();
    getPrescanStates(y-1,maxX-1,states);
    //if(states.size()==0) INTERNAL_ERROR; // ###
    Vector<STATE>::iterator cur=states.begin(), end=states.end();
    for(; cur!=end ; ++cur) {
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
	if(logP>bestP) bestP=logP; 
      }
      pillar[k]=bestP;
    }
  }

  //cout<<"second initialization"<<endl;
  for(int j=maxX-1 ; j>=minX ; --j) 
    {
      // SECOND INITIALIZATION STEP -- fill in top pillar of nextCol:
      
      WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();
      int lowerY=max(minY,nextCol.minY), upperY=min(maxY,nextCol.maxY);
      Array2D<float>::RowIn2DArray<float> pillar=nextCol[upperY];
      Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[upperY];
      bool matchValid=upperY<maxY && upperY<thisCol.maxY;
      pillar.setAllTo(LOG_0);//frame.zeroOut(pillar);
      states.clear();
      getPrescanStates(upperY-1,j-1,states);
      //if(states.size()==0) INTERNAL_ERROR; // ###
      Vector<STATE>::iterator cur=states.begin(), end=states.end();
      for(; cur!=end ; ++cur) {
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
	      emitP=getEmitP(h,upperY/*-1*/,j,PHMM_MATCH);//####
	      if(emitP==LOG_0) continue;
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

      //cout<<"recurrence"<<endl;
      for(int y=upperY-1 ; y>=lowerY ; --y) {
	Array2D<float>::RowIn2DArray<float> pillar=nextCol[y];
	Array2D<float>::RowIn2DArray<float> rightPillar=thisCol[y];
	Array2D<float>::RowIn2DArray<float> upPillar=nextCol[y+1];
	Array2D<float>::RowIn2DArray<float> diagPillar=thisCol[y+1];
	bool rightValid=y>=thisCol.minY, diagValid=y+1>=thisCol.minY;
	pillar.setAllTo(LOG_0);//frame.zeroOut(pillar);
	PHMM_StateType debug=(PHMM_StateType)-1;
	states.clear();
	getPrescanStates(y-1,j-1,states);
	//cout<<states<<endl;
	
	/*
	if(j==1494 && y>=1235 && y<=1238 && hirschberg->getSecondTaxon()->getName()=="vir3") {
	  cout<<"bbb ["<<y<<"] ";
	  Vector<STATE>::iterator cur=states.begin(), end=states.end();
	  for(; cur!=end ; ++cur) {
	    STATE k=*cur;
	    cout<<k<<",";
	  }
	  cout<<endl;
	}
	*/
	
	//if(states.size()==0) INTERNAL_ERROR; // ###
	Vector<STATE>::iterator cur=states.begin(), end=states.end();
	for(; cur!=end ; ++cur) {
	  STATE k=*cur;
	  float bestP=LOG_0;
	  const Array1D<STATE> &succQ=hmm.getSuccessors(k);
	  int numSuccQ=succQ.size();

	  //if(k==721 && j==1494) {for(int ih=0 ; ih<numSuccQ ; ++ih) {STATE h=succQ[ih];cout<<h<<",";}cout<<endl;exit(0);} //### DEBUGGING

	  for(int ih=0 ; ih<numSuccQ ; ++ih) {
	    STATE h=succQ[ih];
	    float emitP, inductiveP;
	    switch(hmm.getStateType(h)) {
	    case PHMM_MATCH:
	      inductiveP=diagValid ? diagPillar[h] : LOG_0;
		
		//if(k==726 && h==724 && j==1495 && y==1235) cout<<"ccc trans="<<getTransP(k,h)<<" ind="<<inductiveP<<" valid="<<diagValid<<" emit="<<getEmitP(h,y,j,PHMM_MATCH)<<endl;
		//if(k>=710 && k<735 && h>=712 && h<735 && j>=1490 && j<=1501 && y>=1231 && y<=1240 && isFinite(inductiveP)) cout<<"ddd trans="<<getTransP(k,h)<<" ind="<<inductiveP<<" valid="<<diagValid<<" emit="<<getEmitP(h,y,j,PHMM_MATCH)<<" k="<<k<<" h="<<h<<" j="<<j<<" y="<<y<<endl;
	      if(inductiveP==LOG_0) continue;
	      emitP=getEmitP(h,y,j,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      break;
	    case PHMM_INSERT:
	      inductiveP=rightValid ? rightPillar[h] : LOG_0;
		//if(k==718 && h==723 && j==1494 && y==1235) cout<<"aaa trans="<<getTransP(k,h)<<" ind="<<inductiveP<<" valid="<<rightValid<<" emit="<<getEmitP(h,IndexMap::UNDEFINED,j,PHMM_INSERT)<<endl;

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
	      emitP=inductiveP=LOG_0;
	      break;
	    }
	    float transP=getTransP(k,h);
	    //if(k==721 && h==726 && j==1494) cout<<"AAA "<<transP<<endl;
	    float logP=transP+emitP+inductiveP;
	    if(logP>bestP) {
	      bestP=logP; 
	      debug=hmm.getStateType(h);
	      //cout<<"ZZZ j="<<j<<" y="<<y<<" "<<debug<<" k="<<k<<endl;
	    }	
	  }
	  pillar[k]=bestP;
	}
	//if(debug>=0) cout<<"ZZZ j="<<j<<" y="<<y<<" "<<debug<<endl;
      }

      // /*
      //if(0)
      //if(frame.getThisColPos()!=j) {
      //cout<<"j="<<j<<" thiscolpos="<<frame.getThisColPos()<<endl;
      //INTERNAL_ERROR;
      //}
      if(columnIsEmpty(frame.getThisColPos()-1,lowerY,upperY)) {
	IndexMap &up=*upMap, &down=*downMap;
	cout<<"\nUPMAP: "<<up<<"\n\nDOWNMAP: "<<down<<"\n"<<endl;
	int parentPos=(*upMap)[j];
	//cout<<"child="<<j<<" parent="<<parentPos<<" upperY="<<upperY<<" lowerY="<<lowerY<<endl;
	//for(int y=prevParentPos ; y<nextParentPos ; ++y) {
	y=parentPos;//+1;
	  states.clear();
	  getPrescanStates(y-1,j-1,states);
	  cout<<"j="<<j<<" y="<<y<<" : ";
	  Vector<STATE>::iterator cur=states.begin(), end=states.end();
	  for(; cur!=end ; ++cur) {
	    STATE q=*cur;
	    if(hmm.getStateType(q)==PHMM_MATCH) {
	      int qq;
	      for(int z=5 ; z<872 ; ++z)
		if(isFinite(getTransP(q,z))) 
		  // && isFinite(frame.getNextCol()[y][z]))
		  cout<<q<<"->"<<z<<" trans="<<getTransP(q,z)<<" emit="<<getEmitP(z,y,j,PHMM_MATCH)<<" induc="<<frame.getThisCol()[y][z]<<" "<<hmm.getStateType(z)<<" "<<hmm.getFunctionalClass(z,PARENT)<<endl;
	     //cout<<q<<":"<<hmm.getStateType(q)<<":"<<hmm.getFunctionalClass(q,PARENT)<<",emit="<<getEmitP(q,IndexMap::UNDEFINED,j,PHMM_INSERT)<<",trans="<<getTransP(97,q)<<",induc="<<frame.getNextCol()[y][j]<<endl;
	    }
	  }
	  cout<<endl;
	  cout<<"parent FC="<<(*parentFP)[y/*-1*/]<<endl;
	  //}
	cout<<firstTaxon->getName()<<" vs. "<<secondTaxon->getName()<<endl;

	WindowColumn &col=frame.getThisCol();//frame.getCol(j+1);
	int minY=col.minY, maxY=col.maxY;
	STATE Q;
	for(int y=minY ; y<maxY ; ++y) {
	  Array2D<float>::RowIn2DArray<float> pillar=col[y];
	  bool found=false;
	  for(int q=5 ; q<numStates ; ++q) {
	    if(isFinite(pillar[q])) {
	      cout<<q<<"="<<hmm.getStateType(q)<<"="<<hmm.getFunctionalClass(q,PARENT)<<"="<<pillar[q]<<" y="<<y<<endl;
	      found=true;
	      Q=q;
	      cout<<"NONZERO "<<Q<<" "<<y<<endl;
	    }
	  }
	  if(found) cout<<"nonzero: "<<y<<endl;
	}
	//cout<<"transP="<<getTransP(2,Q)<<" emitP="<<getEmitP(2,IndexMap::UNDEFINED,j,PHMM_INSERT)<<" inductive="<<col[1235][Q]<<endl;
	cout<<"j="<<j<<endl;
	INTERNAL_ERROR;
      }
      // */

      if(j>minX) frame.advanceWindow();
    }
  //cout<<"/runPrescan"<<endl;
}

