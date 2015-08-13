/****************************************************************
 HirschForwardSum.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "HirschForwardSum.H"
#include "BOOM/Constants.H"
#include "Hirschberg.H"
#include "BOOM/Exceptions.H"
#include "BOOM/SumLogProbs.H"
using namespace std;
using namespace BOOM;



HirschForwardSum::HirschForwardSum(int minX,int maxX,int minY,int maxY,
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
	       FORWARD,vc,pe,parentFP,childFP,upMap,downMap,hirschberg,
	       contentSensor)
{
  //cout<<"minX="<<minX<<" maxX="<<maxX<<" minY="<<minY<<" maxY="<<maxY<<endl;
  frame.initWindow(minX);
}



void HirschForwardSum::run()
{
  WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();

  //thisCol.pred.setAllTo(-1); nextCol.pred.setAllTo(-1); // ### DEBUGGING

  // FIRST INITIALIZATION STEP -- fill in remainder of thisCol:

  int lowerY=max(minY,thisCol.minY), upperY=min(maxY,thisCol.maxY);
  for(int y=lowerY+1 ; y<=upperY ; ++y) {
    Array2D<float>::RowIn2DArray<float> pillar=thisCol[y];
    Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[y-1];
    Array2D<STATE>::RowIn2DArray<STATE> preds=thisCol.indexPred(y);
    pillar.setAllTo(LOG_0);
    for(int ik=0 ; ik<numD ; ++ik) {
      STATE k=QD[ik];
      float emitP=getEmitP(k,y-1,IndexMap::UNDEFINED,PHMM_DELETE);
      if(emitP==LOG_0) continue;
      float bestP=LOG_0;
      STATE bestH;
      const Array1D<STATE> &predQ=hmm.getPredecessors(k);
      int numPredQ=predQ.size();
      for(int ih=0 ; ih<numPredQ ; ++ih) {
	STATE h=predQ[ih];
	float inductiveP=prevPillar[h];
	if(inductiveP==LOG_0) continue;
	float transP=getTransP(h,k);
	float logP=transP+emitP+inductiveP;
	//if(transP>0 || emitP>0 || inductiveP>0) INTERNAL_ERROR; // ###
	if(logP>bestP) { bestP=logP; bestH=h; }
      }
      pillar[k]=bestP;
      preds[k]=bestH;
    }
  }

  for(int j=minX+1 ; j<=maxX ; ++j) 
    {
      // SECOND INITIALIZATION STEP -- fill in bottom pillar of nextCol:
      
      WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();
      int lowerY=max(minY,nextCol.minY), upperY=min(maxY,nextCol.maxY);
      Array2D<float>::RowIn2DArray<float> pillar=nextCol[lowerY];
      Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[lowerY];
      //if(lowerY>thisCol.maxY || lowerY<thisCol.minY) INTERNAL_ERROR;//###DEBUG
      Array2D<STATE>::RowIn2DArray<STATE> preds=nextCol.indexPred(lowerY);
      pillar.setAllTo(LOG_0);
      for(int ik=0 ; ik<numI ; ++ik) {
	STATE k=QI[ik];
	float emitP=getEmitP(k,IndexMap::UNDEFINED,j-1,PHMM_INSERT);
	if(emitP==LOG_0) continue;
	float bestP=LOG_0;
	STATE bestH;
	const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	int numPredQ=predQ.size();
	for(int ih=0 ; ih<numPredQ ; ++ih) {
	  STATE h=predQ[ih];
	  float inductiveP=prevPillar[h];
	  if(inductiveP==LOG_0) continue;
	  float transP=getTransP(h,k);
	  float logP=transP+emitP+inductiveP;
	  if(logP>bestP) { bestP=logP; bestH=h; }
	}
	pillar[k]=bestP;
	preds[k]=bestH;	
      }
      if(lowerY>minY && lowerY>thisCol.minY) {
	//if(lowerY-1<thisCol.minY || lowerY-1>thisCol.maxY) INTERNAL_ERROR;//###DEBUGGING
	Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[lowerY-1];
	for(int ik=0 ; ik<numM ; ++ik) {
	  STATE k=QM[ik];
	  float emitP=getEmitP(k,lowerY-1,j-1,PHMM_MATCH);
	  if(emitP==LOG_0) continue;
	  float bestP=LOG_0;
	  STATE bestH;
	  const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	  int numPredQ=predQ.size();
	  for(int ih=0 ; ih<numPredQ ; ++ih) {
	    STATE h=predQ[ih];
	    float inductiveP=prevPillar[h];
	    if(inductiveP==LOG_0) continue;
	    float transP=getTransP(h,k);
	    float logP=transP+emitP+inductiveP;
	    if(logP>bestP) { bestP=logP; bestH=h; }
	  }
	  pillar[k]=bestP;
	  preds[k]=bestH;	
	}
      }
      
    // RECURRENCE -- fill in the rest of nextCol:

      for(int y=lowerY+1 ; y<=upperY ; ++y) {
	Array2D<float>::RowIn2DArray<float> pillar=nextCol[y];
	Array2D<float>::RowIn2DArray<float> leftPillar=thisCol[y];
	Array2D<float>::RowIn2DArray<float> downPillar=nextCol[y-1];
	Array2D<float>::RowIn2DArray<float> diagPillar=thisCol[y-1];
	Array2D<STATE>::RowIn2DArray<STATE> preds=nextCol.indexPred(y);
	bool leftValid=y<=thisCol.maxY, diagValid=y-1<=thisCol.maxY;
	pillar[0]=LOG_0;
	for(int ik=0 ; ik<numM ; ++ik) { // MATCH STATES
	  STATE k=QM[ik];
	  float emitP=getEmitP(k,y-1,j-1,PHMM_MATCH);
	  if(emitP==LOG_0) continue;
	  float bestP=LOG_0;
	  STATE bestH;
	  const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	  int numPredQ=predQ.size();
	  for(int ih=0 ; ih<numPredQ ; ++ih) {
	    STATE h=predQ[ih];
	    float inductiveP=diagValid ? diagPillar[h] : LOG_0;
	    if(inductiveP==LOG_0) continue;
	    float transP=getTransP(h,k);
	    float logP=transP+emitP+inductiveP;
	    if(logP>bestP) { bestP=logP; bestH=h; }
	  }
	  pillar[k]=bestP;
	  preds[k]=bestH;	
	}
	for(int ik=0 ; ik<numI ; ++ik) { // INSERT STATES
	  STATE k=QI[ik];
	  float emitP=getEmitP(k,IndexMap::UNDEFINED,j-1,PHMM_INSERT);
	  if(emitP==LOG_0) continue;
	  float bestP=LOG_0;
	  STATE bestH;
	  const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	  int numPredQ=predQ.size();
	  for(int ih=0 ; ih<numPredQ ; ++ih) {
	    STATE h=predQ[ih];
	    float inductiveP=leftValid ? leftPillar[h] : LOG_0;
	    if(inductiveP==LOG_0) continue;
	    float transP=getTransP(h,k);
	    float logP=transP+emitP+inductiveP;
	    if(logP>bestP) { bestP=logP; bestH=h; }
	  }
	  pillar[k]=bestP;
	  preds[k]=bestH;	
	}
	for(int ik=0 ; ik<numD ; ++ik) { // DELETE STATES
	  STATE k=QD[ik];
	  float emitP=getEmitP(k,y-1,IndexMap::UNDEFINED,PHMM_DELETE);
	  if(emitP==LOG_0) continue;
	  float bestP=LOG_0;
	  STATE bestH=0; //-1; //### DEBUGGING
	  const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	  int numPredQ=predQ.size();
	  for(int ih=0 ; ih<numPredQ ; ++ih) {
	    STATE h=predQ[ih];
	    float inductiveP=downPillar[h];
	    if(inductiveP==LOG_0) continue;
	    float transP=getTransP(h,k);
	    float logP=transP+emitP+inductiveP;
	    if(logP>bestP) { bestP=logP; bestH=h; }
	  }
	  pillar[k]=bestP;
	  preds[k]=bestH;	
	}
      }
      if(j<maxX) frame.advanceWindow();
    }
}



/****************************************************************
                         runPrescan()
 ****************************************************************/
void HirschForwardSum::runPrescan()
{
  WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();

  // FIRST INITIALIZATION STEP -- fill in remainder of thisCol:

  Vector<float> logProbs;
  Vector<STATE> states;
  int lowerY=max(minY,thisCol.minY), upperY=min(maxY,thisCol.maxY);
  for(int y=lowerY+1 ; y<=upperY ; ++y) {
    Array2D<float>::RowIn2DArray<float> pillar=thisCol[y];
    Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[y-1];
    pillar.setAllTo(LOG_0);
    states.clear();
    getPrescanStates(y-1,minX-1,states);
    Vector<STATE>::iterator cur=states.begin(), end=states.end();
    for(; cur!=end ; ++cur) {
      logProbs.clear();
      STATE k=*cur;
      if(hmm.getStateType(k)!=PHMM_DELETE) continue;
      float emitP=getEmitP(k,y-1,IndexMap::UNDEFINED,PHMM_DELETE);
      if(emitP==LOG_0) continue;
      const Array1D<STATE> &predQ=hmm.getPredecessors(k);
      int numPredQ=predQ.size();
      for(int ih=0 ; ih<numPredQ ; ++ih) {
	STATE h=predQ[ih];
	float inductiveP=prevPillar[h];
	if(inductiveP==LOG_0) continue;
	float transP=getTransP(h,k);
	float logP=transP+emitP+inductiveP;
	logProbs.push_back(logP);
      }
      pillar[k]=sumLogProbs<float>(logProbs);
    }
  }

  for(int j=minX+1 ; j<=maxX ; ++j) 
    {
      // SECOND INITIALIZATION STEP -- fill in bottom pillar of nextCol:
      
      WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();
      int lowerY=max(minY,nextCol.minY), upperY=min(maxY,nextCol.maxY);
      Array2D<float>::RowIn2DArray<float> pillar=nextCol[lowerY];
      Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[lowerY];
      Array2D<STATE>::RowIn2DArray<STATE> preds=nextCol.indexPred(lowerY);
      pillar.setAllTo(LOG_0);
      states.clear();
      getPrescanStates(lowerY-1,j-1,states);
      Vector<STATE>::iterator cur=states.begin(), end=states.end();
      for(; cur!=end ; ++cur) {
	logProbs.clear();
	STATE k=*cur;
	if(hmm.getStateType(k)!=PHMM_INSERT) continue;
	float emitP=getEmitP(k,IndexMap::UNDEFINED,j-1,PHMM_INSERT);
	if(emitP==LOG_0) continue;
	//float bestP=LOG_0;
	//STATE bestH;
	const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	int numPredQ=predQ.size();
	for(int ih=0 ; ih<numPredQ ; ++ih) {
	  STATE h=predQ[ih];
	  float inductiveP=prevPillar[h];
	  if(inductiveP==LOG_0) continue;
	  float transP=getTransP(h,k);
	  float logP=transP+emitP+inductiveP;
	  logProbs.push_back(logP);
	}
	pillar[k]=sumLogProbs<float>(logProbs);
      }
      if(lowerY>minY && lowerY>thisCol.minY) {
	Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[lowerY-1];
	cur=states.begin();
	for(; cur!=end ; ++cur) {
	  logProbs.clear();
	  STATE k=*cur;
	  if(hmm.getStateType(k)!=PHMM_MATCH) continue;
	  float emitP=getEmitP(k,lowerY-1,j-1,PHMM_MATCH);
	  if(emitP==LOG_0) continue;
	  //float bestP=LOG_0;
	  //STATE bestH;
	  const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	  int numPredQ=predQ.size();
	  for(int ih=0 ; ih<numPredQ ; ++ih) {
	    STATE h=predQ[ih];
	    float inductiveP=prevPillar[h];
	    if(inductiveP==LOG_0) continue;
	    float transP=getTransP(h,k);
	    float logP=transP+emitP+inductiveP;
	    logProbs.push_back(logP);
	  }
	  pillar[k]=sumLogProbs<float>(logProbs);
	}
      }
      
      // RECURRENCE -- fill in the rest of nextCol:

      for(int y=lowerY+1 ; y<=upperY ; ++y) {
	Array2D<float>::RowIn2DArray<float> pillar=nextCol[y];
	Array2D<float>::RowIn2DArray<float> leftPillar=thisCol[y];
	Array2D<float>::RowIn2DArray<float> downPillar=nextCol[y-1];
	Array2D<float>::RowIn2DArray<float> diagPillar=thisCol[y-1];
	Array2D<STATE>::RowIn2DArray<STATE> preds=nextCol.indexPred(y);
	bool leftValid=y<=thisCol.maxY, diagValid=y-1<=thisCol.maxY;
	pillar.setAllTo(LOG_0);//frame.zeroOut(pillar);
	states.clear();
	getPrescanStates(y-1,j-1,states);
	//if(states.size()==0) INTERNAL_ERROR; // ###
	Vector<STATE>::iterator cur=states.begin(), end=states.end();
	for(; cur!=end ; ++cur) {
	  logProbs.clear();
	  STATE k=*cur;
	  switch(hmm.getStateType(k)) 
	    {
	    case PHMM_MATCH: {
	      float emitP=getEmitP(k,y-1,j-1,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      //float bestP=LOG_0;
	      //STATE bestH;
	      const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	      int numPredQ=predQ.size();
	      for(int ih=0 ; ih<numPredQ ; ++ih) {
		STATE h=predQ[ih];
		float inductiveP=diagValid ? diagPillar[h] : LOG_0;
		if(inductiveP==LOG_0) continue;
		float transP=getTransP(h,k);
		float logP=transP+emitP+inductiveP;
		logProbs.push_back(logP);
	      }
	    }
	      break;
	    case PHMM_INSERT: {
	      float emitP=getEmitP(k,IndexMap::UNDEFINED,j-1,PHMM_INSERT);
	      if(emitP==LOG_0) continue;
	      //float bestP=LOG_0;
	      //STATE bestH;
	      const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	      int numPredQ=predQ.size();
	      for(int ih=0 ; ih<numPredQ ; ++ih) {
		STATE h=predQ[ih];
		float inductiveP=leftValid ? leftPillar[h] : LOG_0;
		if(inductiveP==LOG_0) continue;
		float transP=getTransP(h,k);
		float logP=transP+emitP+inductiveP;
		logProbs.push_back(logP);
	      }
	    }
	      break;
	    case PHMM_DELETE: {
	      float emitP=getEmitP(k,y-1,IndexMap::UNDEFINED,PHMM_DELETE);
	      if(emitP==LOG_0) continue;
	      //float bestP=LOG_0;
	      //STATE bestH;//=0; //-1; //### DEBUGGING
	      const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	      int numPredQ=predQ.size();
	      for(int ih=0 ; ih<numPredQ ; ++ih) {
		STATE h=predQ[ih];
		float inductiveP=downPillar[h];
		if(inductiveP==LOG_0) continue;
		float transP=getTransP(h,k);
		float logP=transP+emitP+inductiveP;
		logProbs.push_back(logP);
	      }
	    } // case
	      break;
	      //default: continue;

	    } // switch
	  pillar[k]=sumLogProbs<float>(logProbs);
	} // foreach state
      } // foreach y
      if(j<maxX) frame.advanceWindow();
    } //foreach j
}



float HirschForwardSum::scorePath(StatePath &path)
{
  int L=path.length(), i=0, j=0;
  STATE prevState=0;
  double LL=0;
  for(int k=0 ; k<L ; ++k) {
    STATE q=path[k];
    hmm.updateColumnsFwd(q,i,j);
    double transP=getTransP(prevState,q);
    double emitP=emission(q,i-1,j-1,hmm.getStateType(q));
    LL+=transP+emitP;
    prevState=q;
  }
  LL+=getTransP(prevState,0);
  return LL;
}



double HirschForwardSum::computeLikelihood()
{
  WindowColumn &col0=frame.getThisCol();
  Array2D<float>::RowIn2DArray<float> fPillar=col0[0];
  fPillar.setAllTo(LOG_0);
  fPillar[0]=LOG_1;
  runPrescan();
  WindowColumn &col=frame.getThisCol();//getNextCol();
  int minY=col.minY, maxY=col.maxY;
  Vector<float> V;
  for(int y=minY ; y<=maxY ; ++y) {
    Array2D<float>::RowIn2DArray<float> pillar=col[y];
    for(STATE q=1 ; q<numStates ; ++q) {
      double inductive=pillar[q], trans=hmm.getTransP(q,0);
      double both=inductive+trans;
      if(isFinite(both)) V.push_back(both);
    }
  }
  double LL=sumLogProbs<float>(V);
  if(!isFinite(LL)) {
    cout<<"XXX: "<<V<<endl;
    //cout<<"THISCOL: "<<frame.getThisCol()<<"\nNEXTCOL: "<<frame.getNextCol()<<endl;
  }
  return LL;
}



void HirschForwardSum::debug()
{
  {
    WindowColumn &col0=frame.getThisCol();
    Array2D<float>::RowIn2DArray<float> fPillar=col0[0];
    fPillar.setAllTo(LOG_0);
    fPillar[0]=LOG_1;
  }

  WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();

  // FIRST INITIALIZATION STEP -- fill in remainder of thisCol:

  Vector<float> logProbs;
  Vector<STATE> states;
  int lowerY=max(minY,thisCol.minY), upperY=min(maxY,thisCol.maxY);
  cout<<"lowerY="<<lowerY<<" upperY="<<upperY<<endl;
  for(int y=lowerY+1 ; y<=upperY ; ++y) {
    Array2D<float>::RowIn2DArray<float> pillar=thisCol[y];
    Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[y-1];
    pillar.setAllTo(LOG_0);
    states.clear();
    getPrescanStates(y-1,minX-1,states);
    //cout<<"DDD "<<states<<endl;
    Vector<STATE>::iterator cur=states.begin(), end=states.end();
    for(; cur!=end ; ++cur) {
      logProbs.clear();
      STATE k=*cur;
      if(hmm.getStateType(k)!=PHMM_DELETE) continue;
      float emitP=getEmitP(k,y-1,IndexMap::UNDEFINED,PHMM_DELETE);
      if(k==3 && y==1) {
	cout<<"k=3 y=1 emit="<<emitP<<endl;
	FunctionalClass fc=hmm.getFunctionalClass(k,PARENT);
	cout<<"fc="<<fc<<endl;
	Array1D<double> &eqFreqs=fc.getEqFreqs();
	cout<<"eqFreqs="<<eqFreqs<<endl;
	Array3D<float>::IndexedTwice<float> row=
	  precomputedEmissions[PARENT][y-1][fc];
	Array1D<float> V(numAlpha);
	for(Symbol s=0 ; s<numAlpha ; ++s) {
	  double eq=getEQ(eqFreqs,s,firstTaxon,y-1,
			  hmm.getFunctionalClass(k,PARENT));
	  V[s]=safeAdd(/*consFactor**/row[s],/*eqFactor**/eq);
	  cout<<"  s="<<int(s)<<" row[s]="<<row[s]<<" eq="<<eq<<endl;
	}
	cout<<"V="<<V<<endl;
	double emit=sumLogProbs<float>(V);
	cout<<"emit="<<emit<<endl;
	INTERNAL_ERROR;
      }
      if(emitP==LOG_0) continue;
      const Array1D<STATE> &predQ=hmm.getPredecessors(k);
      int numPredQ=predQ.size();
      for(int ih=0 ; ih<numPredQ ; ++ih) {
	STATE h=predQ[ih];
	float inductiveP=prevPillar[h];
	if(k==3 && h==0 && y==1) cout<<"k=3 h=0 y=1 inductive="<<inductiveP<<endl;
	if(inductiveP==LOG_0) continue;
	float transP=getTransP(h,k);
	float logP=transP+emitP+inductiveP;
	logProbs.push_back(logP);
	if(k<4) cout<<"inductive="<<inductiveP<<" trans="<<transP<<" emit="<<emitP<<" k="<<k<<" h="<<h<<endl;
      }
      pillar[k]=sumLogProbs<float>(logProbs);
    }
  }

  for(int j=minX+1 ; j<=maxX ; ++j) 
    {
      // SECOND INITIALIZATION STEP -- fill in bottom pillar of nextCol:
      
      WindowColumn &thisCol=frame.getThisCol(), &nextCol=frame.getNextCol();
      int lowerY=max(minY,nextCol.minY), upperY=min(maxY,nextCol.maxY);
      Array2D<float>::RowIn2DArray<float> pillar=nextCol[lowerY];
      Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[lowerY];
      Array2D<STATE>::RowIn2DArray<STATE> preds=nextCol.indexPred(lowerY);
      pillar.setAllTo(LOG_0);
      states.clear();
      getPrescanStates(lowerY-1,j-1,states);
      if(j==200) cout<<"BBB y="<<lowerY<<" "<<states<<endl;
      Vector<STATE>::iterator cur=states.begin(), end=states.end();
      for(; cur!=end ; ++cur) {
	logProbs.clear();
	STATE k=*cur;
	if(hmm.getStateType(k)!=PHMM_INSERT) continue;
	float emitP=getEmitP(k,IndexMap::UNDEFINED,j-1,PHMM_INSERT);
	if(emitP==LOG_0) continue;
	const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	int numPredQ=predQ.size();
	for(int ih=0 ; ih<numPredQ ; ++ih) {
	  STATE h=predQ[ih];
	  float inductiveP=prevPillar[h];
	  if(inductiveP==LOG_0) continue;
	  float transP=getTransP(h,k);
	  float logP=transP+emitP+inductiveP;
	  logProbs.push_back(logP);
	}
	pillar[k]=sumLogProbs<float>(logProbs);
      }
      if(lowerY>minY && lowerY>thisCol.minY) {
	Array2D<float>::RowIn2DArray<float> prevPillar=thisCol[lowerY-1];
	cur=states.begin();
	for(; cur!=end ; ++cur) {
	  logProbs.clear();
	  STATE k=*cur;
	  if(hmm.getStateType(k)!=PHMM_MATCH) continue;
	  float emitP=getEmitP(k,lowerY-1,j-1,PHMM_MATCH);
	  if(emitP==LOG_0) continue;
	  const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	  int numPredQ=predQ.size();
	  for(int ih=0 ; ih<numPredQ ; ++ih) {
	    STATE h=predQ[ih];
	    float inductiveP=prevPillar[h];
	    if(inductiveP==LOG_0) continue;
	    float transP=getTransP(h,k);
	    float logP=transP+emitP+inductiveP;
	    logProbs.push_back(logP);
	  }
	  pillar[k]=sumLogProbs<float>(logProbs);
	}
      }
      
      // RECURRENCE -- fill in the rest of nextCol:

      for(int y=lowerY+1 ; y<=upperY ; ++y) {
	Array2D<float>::RowIn2DArray<float> pillar=nextCol[y];
	Array2D<float>::RowIn2DArray<float> leftPillar=thisCol[y];
	Array2D<float>::RowIn2DArray<float> downPillar=nextCol[y-1];
	Array2D<float>::RowIn2DArray<float> diagPillar=thisCol[y-1];
	Array2D<STATE>::RowIn2DArray<STATE> preds=nextCol.indexPred(y);
	bool leftValid=y<=thisCol.maxY, diagValid=y-1<=thisCol.maxY;
	pillar.setAllTo(LOG_0);//frame.zeroOut(pillar);
	states.clear();
	getPrescanStates(y-1,j-1,states);
	//cout<<states.size()<<" prescan states"<<endl;
	//if(j==200) cout<<"AAA y="<<y<<" "<<states<<endl;
	//if(states.size()==0) INTERNAL_ERROR; // ###
	Vector<STATE>::iterator cur=states.begin(), end=states.end();
	for(; cur!=end ; ++cur) {
	  logProbs.clear();
	  STATE k=*cur;
	  switch(hmm.getStateType(k)) 
	    {
	    case PHMM_MATCH: {
	      float emitP=getEmitP(k,y-1,j-1,PHMM_MATCH);
	      if(emitP==LOG_0) continue;
	      const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	      int numPredQ=predQ.size();
	      for(int ih=0 ; ih<numPredQ ; ++ih) {
		STATE h=predQ[ih];
		float inductiveP=diagValid ? diagPillar[h] : LOG_0;
		if(inductiveP==LOG_0) continue;
		float transP=getTransP(h,k);
		float logP=transP+emitP+inductiveP;
		logProbs.push_back(logP);
	      }
	    }
	      break;
	    case PHMM_INSERT: {
	      float emitP=getEmitP(k,IndexMap::UNDEFINED,j-1,PHMM_INSERT);
	      if(emitP==LOG_0) continue;
	      const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	      int numPredQ=predQ.size();
	      for(int ih=0 ; ih<numPredQ ; ++ih) {
		STATE h=predQ[ih];
		float inductiveP=leftValid ? leftPillar[h] : LOG_0;
		if(inductiveP==LOG_0) continue;
		float transP=getTransP(h,k);
		float logP=transP+emitP+inductiveP;
		logProbs.push_back(logP);
	      }
	    }
	      break;
	    case PHMM_DELETE: {
	      float emitP=getEmitP(k,y-1,IndexMap::UNDEFINED,PHMM_DELETE);
	      if(emitP==LOG_0) continue;
	      const Array1D<STATE> &predQ=hmm.getPredecessors(k);
	      int numPredQ=predQ.size();
	      for(int ih=0 ; ih<numPredQ ; ++ih) {
		STATE h=predQ[ih];
		float inductiveP=downPillar[h];
		if(inductiveP==LOG_0) continue;
		float transP=getTransP(h,k);
		float logP=transP+emitP+inductiveP;
		logProbs.push_back(logP);
	      }
	    } // case
	      break;
	    default: continue;

	    } // switch
	  pillar[k]=sumLogProbs<float>(logProbs);
	} // foreach state
      } // foreach y
      //if(j==199) cout<<nextCol<<endl;
      if(columnIsEmpty(j,lowerY,upperY)) {
	cout<<"EMPTY COLUMN: "<<j<<" lowerY="<<lowerY<<" upperY="<<upperY<<endl;
	INTERNAL_ERROR;
      }
      if(j<maxX) frame.advanceWindow();
    } //foreach j
}


