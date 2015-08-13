/****************************************************************
 LinkViterbi.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "LinkViterbi.H"
#include "BOOM/Constants.H"
#include "BOOM/Exceptions.H"
#include "BOOM/SumLogProbs.H"
using namespace std;
using namespace BOOM;


LinkViterbi::LinkViterbi(Taxon *parent,Taxon *child,LinkFelsenstein &F,
			 BandingType bandingType,int bandwidth)
  : parent(parent), child(child), bandingType(bandingType),
    bandwidth(bandwidth), F(F), posteriors(NULL), vc(VC_UNCONSTRAINED),
    alphabet(F.getAlphabet())
{
  // ctor

  numAlpha=alphabet.size();
  if(parent && child) {
    parentFP=&parent->getFunctionalParse();
    childFP=&child->getFunctionalParse();
    branch=parent->getBranchToChild(*child);
    hmm=branch->getHMM();
    numStates=hmm->getNumStates();
    L1=parent->getSeqLen();
    L2=child->getSeqLen();
    V.resize(L1+1,L2+1,numStates,bandwidth);
    T.resize(L1+1,L2+1,numStates,bandwidth);
    upMap=&branch->getUpMap();
    downMap=&branch->getDownMap();
    precomputeEmissions();
  }
  else {
    parentFP=childFP=NULL;
    branch=NULL;
    hmm=NULL;
    numStates=0;
    L1=L2=0;
    upMap=downMap=NULL;
  }
}



void LinkViterbi::usePosteriors(Posteriors *p)
{
  posteriors=p;
}



StatePath *LinkViterbi::decode(ViterbiConstraint constraint)
{
  vc=constraint;
  fillMatrix();
  return getPath(bestScore);
}



void LinkViterbi::fillMatrix()
{
  // Misc initialization:
  Vector<STATE> I, D, M;
  hmm->getStatesOfType(PHMM_INSERT,I);
  hmm->getStatesOfType(PHMM_DELETE,D);
  hmm->getStatesOfType(PHMM_MATCH,M);
  int nI=I.size(), nD=D.size(), nM=M.size();

  // Initialization of DP matrix:
  V.setAllTo(LOG_0);
  T.setAllTo(-1);
  V(0,0,0)=LOG_1;
  Array1D<float> *column, *prevCol;
  for(int i=1 ; i<=L1 ; ++i) {
    if(!V.getColumn(i,0,column)) break;
    V.getColumn(i-1,0,prevCol);
    for(int ik=0 ; ik<nD ; ++ik) { // deletion states
      STATE k=D[ik], bestH=0;
      float Pe=getEmitP(k,i-1,IndexMap::UNDEFINED,PHMM_DELETE);
      float bestScore=(*prevCol)[0]+getTransP(0,k);
      const Array1D<STATE> &pred=hmm->getPredecessors(k);
      STATE *p=&pred[0];
      int numPred=pred.size();
      for(int s=0 ; s<numPred ; ++s, ++p) {
	STATE h=*p;
	if(hmm->getStateType(h)!=PHMM_DELETE) continue;
	float score=(*prevCol)[h]+getTransP(h,k);
	if(score>bestScore) { bestH=h; bestScore=score; }
      }
      (*column)[k]=bestScore+Pe;
      T(i,0,k)=bestH;
    }
  }
  for(int j=1 ; j<=L2 ; ++j) {
    if(!V.getColumn(0,j,column)) break;
    V.getColumn(0,j-1,prevCol);
    for(int ik=0 ; ik<nI ; ++ik) { // insertion states
      STATE k=I[ik], bestH=0;
      float Pe=getEmitP(k,IndexMap::UNDEFINED,j-1,PHMM_INSERT);
      float bestScore=(*prevCol)[0]+getTransP(0,k);
      const Array1D<STATE> &pred=hmm->getPredecessors(k);
      STATE *p=&pred[0];
      int numPred=pred.size();
      for(int s=0 ; s<numPred ; ++s, ++p) {
	STATE h=*p;
	if(hmm->getStateType(h)!=PHMM_INSERT) continue;
	float score=(*prevCol)[h]+getTransP(h,k);
	if(score>bestScore) { bestH=h; bestScore=score; }
      }
      (*column)[k]=bestScore+Pe;
      T(0,j,k)=bestH;
    }
  }

  // Now for the recursion:

  Array1D<float> *Vij, *V_im1_jm1, *V_i_jm1, *V_im1_j;
  Array1D<short> *Tij;
  STATE *pk;
  for(int i=1 ; i<=L1 ; ++i) {
    for(int j=1 ; j<=L2 ; ++j) {
      if(!V.getColumn(i,j,Vij)) continue;
      V_im1_jm1=V_i_jm1=V_im1_j=NULL;
      V.getColumn(i-1,j-1,V_im1_jm1);
      V.getColumn(i,j-1,V_i_jm1); 
      V.getColumn(i-1,j,V_im1_j);
      T.getColumn(i,j,Tij);
      if(V_im1_jm1) {
	pk=&M[0];
	for(int ik=0 ; ik<nM ; ++ik, ++pk) { // match states
	  STATE k=*pk, bestH=-2;
	  float Pe=getEmitP(k,i-1,j-1,PHMM_MATCH);
	  float bestScore=NEGATIVE_INFINITY;
	  if(!isinf(Pe)) {
	    const Array1D<STATE> &pred=hmm->getPredecessors(k);
	    STATE *p=&pred[0];
	    int numPred=pred.size();
	    for(int s=0 ; s<numPred ; ++s, ++p) {
	      STATE h=*p;
	      float score=(*V_im1_jm1)[h]+getTransP(h,k);
	      if(score>bestScore) { bestH=h; bestScore=score; }
	    }
	  }
	  (*Vij)[k]=bestScore+Pe;
	  (*Tij)[k]=bestH;
	}
      }
      if(V_im1_j) {
	pk=&D[0];
	for(int ik=0 ; ik<nD ; ++ik, ++pk) { // deletion states
	  STATE k=*pk, bestH=-2;
	  float Pe=getEmitP(k,i-1,IndexMap::UNDEFINED,PHMM_DELETE);
	  float bestScore=NEGATIVE_INFINITY;
	  if(!isinf(Pe)) {
	    const Array1D<STATE> &pred=hmm->getPredecessors(k);
	    STATE *p=&pred[0];
	    int numPred=pred.size();
	    for(int s=0 ; s<numPred ; ++s, ++p) {
	      STATE h=*p;
	      float score=(*V_im1_j)[h]+getTransP(h,k);
	      if(score>bestScore) { bestH=h; bestScore=score; }
	    }
	  }
	  (*Vij)[k]=bestScore+Pe;
	  (*Tij)[k]=bestH;
	}
      }
      if(V_i_jm1) {
	pk=&I[0];
	for(int ik=0 ; ik<nI ; ++ik, ++pk) { // insertion states
	  STATE k=*pk, bestH=-2;
	  float Pe=getEmitP(k,IndexMap::UNDEFINED,j-1,PHMM_INSERT);
	  float bestScore=NEGATIVE_INFINITY;
	  if(!isinf(Pe)) {
	    const Array1D<STATE> &pred=hmm->getPredecessors(k);
	    STATE *p=&pred[0];
	    int numPred=pred.size();
	    for(int s=0 ; s<numPred ; ++s, ++p) {
	      STATE h=*p;
	      float score=(*V_i_jm1)[h]+getTransP(h,k);
	      if(score>bestScore) { bestH=h; bestScore=score; }
	    }
	  }
	  (*Vij)[k]=bestScore+Pe;
	  (*Tij)[k]=bestH;
	}
      }
    }
  }
}



StatePath *LinkViterbi::getPath(double &bestScore)
{
  // Initialization:

  StatePath *path=new StatePath(hmm);
  int bestH=-1;
  bestScore=NEGATIVE_INFINITY;
  Array1D<float> *V_L1_L2;
  V.getColumn(L1,L2,V_L1_L2);
  for(STATE h=1 ; h<numStates ; ++h) {
    float score=(*V_L1_L2)[h]+getTransP(h,0);
    if(score>bestScore) { bestH=h; bestScore=score; }
  }
  if(isinf(bestScore)) {
    for(STATE h=1 ; h<numStates ; ++h)
      cout<<(*V_L1_L2)[h]<<" "<<getTransP(h,0)<<endl;
    cout<<"L1="<<L1<<" L2="<<L2<<endl;
    cout<<"traceback failed, bestScore is -INF"<<endl;
    debugTraceback();
    throw "traceback failed, bestScore is -INF";
  }

  // Recursion:

  int i=L1, j=L2;
  STATE t=bestH;
  while(t>0) {
    path->push_back(t);
    STATE newT=T(i,j,t);
    hmm->updateColumnsRev(t,i,j);
    t=newT;
  }

  // Termination:
  StatePath *revPath=path->getReverse();
  delete path;
  revPath->setScore(bestScore);
  return revPath;
}



void LinkViterbi::debugTraceback()
{
  if(!parentFP) return;
  cout<<"&parentFP="<<parentFP<<endl;
  cout<<"parentFP="<<*parentFP<<endl;
  IndexMap &down=branch->getDownMap(), &up=branch->getUpMap();
  cout<<"downMap="<<down<<"\nupMap="<<up<<endl;
  int x=0;
  while(x<L1) {
    while(x<L1 && down[x]==IndexMap::UNDEFINED) ++x;
    int y=down[x];
    if(x==L1 || y==L2) break;
    bool f=false;
    for(STATE s=0 ; s<numStates ; ++s) {
      if(isFinite(V(x+1,y+1,s))) f=true;
    }
    cout<<"("<<x<<","<<y<<")="<<f<<endl;//<<" pFC="<<(*parentFP)[x]<<endl;
    if(!f) {
      for(int i=0 ; i<=x+1 ; ++i)
	for(int j=0 ; j<=y+1 ; ++j) {
	  cout<<"  ["<<i<<","<<j<<"]=";
	  for(STATE s=0 ; s<numStates ; ++s)
	    cout<<V(i,j,s)<<" ";
	  cout<<endl;
	}
      throw "LinkViterbi::debugTraceback";
    }
    ++x;
  }
}



float LinkViterbi::getTransP(STATE h,STATE k)
{
  float p=hmm->getTransP(h,k);
  if(posteriors) return finite(p) ? LOG_1 : LOG_0;
  else return p;
}



double LinkViterbi::getEmitP(STATE k,int p,int c,PairHMMStateType t)
{
  if(posteriors) return posteriors->compute(p,c,k);
  switch(vc)
    {
    case VC_UNCONSTRAINED: 
      break;
    case VC_PARENT_PARSE:
      if(emitsParent(t) && hmm->getFunctionalClass(k,PARENT)!=(*parentFP)[p])
	return LOG_0 ;
      break;
    case VC_CHILD_PARSE:
      if(emitsChild(t) && hmm->getFunctionalClass(k,CHILD)!=(*childFP)[c])
	return LOG_0;
      break;
    case VC_PARENT_AND_CHILD_PARSE:
      if(emitsParent(t) && hmm->getFunctionalClass(k,PARENT)!=(*parentFP)[p] ||
	 emitsChild(t) && hmm->getFunctionalClass(k,CHILD)!=(*childFP)[c])
	return LOG_0;
      break;
    case VC_INDEL_HISTORY:
      if(emitsParent(t) && (*downMap)[p]!=c ||
	 emitsChild(t) && (*upMap)[c]!=p)
	return LOG_0;
      break;
    case VC_INDEL_AND_PARENT_PARSE:
      if(emitsParent(t))
	if((*downMap)[p]!=c || 
	   hmm->getFunctionalClass(k,PARENT)!=(*parentFP)[p]) return LOG_0;
      if(emitsChild(t) && (*upMap)[c]!=p) return LOG_0 ;
      break;
    case VC_INDEL_AND_CHILD_PARSE:
      if(emitsParent(t) && (*downMap)[p]!=c) return LOG_0;
      if(emitsChild(t))
	if((*upMap)[c]!=p ||
	   hmm->getFunctionalClass(k,CHILD)!=(*childFP)[c]) return LOG_0;
      break;
    }
  return emission(k,p,c,t);
}



double LinkViterbi::emission(STATE k,int parentResidueIndex,
			     int childResidueIndex,
			     PairHMMStateType stateType)
{
  Array1D<float> V(numAlpha), V2(numAlpha);
  switch(stateType)
    {
    case PHMM_MATCH:{
      FunctionalClass fcParent=hmm->getFunctionalClass(k,PARENT);
      FunctionalClass fcChild=hmm->getFunctionalClass(k,CHILD);
      const Array1D<double> &eqFreqs=hmm->getEqFreqs(k);
      Array3D<float>::IndexedTwice<float> childRow=
	precomputedEmissions[CHILD][childResidueIndex][fcChild];
      Array3D<float>::IndexedTwice<float> parentRow=
	precomputedEmissions[PARENT][parentResidueIndex][fcParent];
      SubstitutionMatrix &Pt=*hmm->getSubstMatrix(k);
      for(Symbol sParent=0 ; sParent<numAlpha ; ++sParent) {
	double eq=log(eqFreqs[sParent]);
	for(Symbol sChild=0 ; sChild<numAlpha ; ++sChild)
	  V[sChild]=safeAdd(eq,parentRow[sParent],
			     safeAdd(Pt(sParent,sChild),childRow[sChild]));
	V2[sParent]=sumLogProbs<float>(V);
      }
      return sumLogProbs<float>(V2);
      break;}
    case PHMM_INSERT:{
      FunctionalClass fc=hmm->getFunctionalClass(k,CHILD);
      Array1D<double> &eq=fc.getEqFreqs();
      Array3D<float>::IndexedTwice<float> row=
	precomputedEmissions[CHILD][childResidueIndex][fc];
      for(Symbol s=0 ; s<numAlpha ; ++s)
	V[s]=safeAdd(row[s],log(eq[s]));
      return sumLogProbs<float>(V);
      break;}
    case PHMM_DELETE:{
      FunctionalClass fc=hmm->getFunctionalClass(k,PARENT);
      Array1D<double> &eq=fc.getEqFreqs();
      Array3D<float>::IndexedTwice<float> row=
	precomputedEmissions[PARENT][parentResidueIndex][fc];
      for(Symbol s=0 ; s<numAlpha ; ++s)
	V[s]=safeAdd(row[s],log(eq[s]));
      return sumLogProbs<float>(V);
      break;}
    default:
      INTERNAL_ERROR;
    }
}



float LinkViterbi::scorePath(StatePath &path)
{
  int L=path.length(), i=1, j=1;
  STATE prevState=0;
  float LL=0;
  for(int k=0 ; k<L ; ++k) {
    STATE q=path[k];
    float transP=getTransP(prevState,q);
    float emitP=emission(q,i-1,j-1,hmm->getStateType(q));
    LL+=transP+emitP;
    hmm->updateColumnsFwd(q,i,j);
    prevState=q;
  }
  LL+=getTransP(prevState,0);
  return LL;
}



void LinkViterbi::precomputeEmissions()
{
  precomputeEmissions(child,CHILD,INSIDE);
  precomputeEmissions(parent,PARENT,OUTSIDE,child);
}



void LinkViterbi::precomputeEmissions(Taxon *taxon,BranchEnd parentChild,
				      InsideOutside insideOutside,
				      Taxon *avoid)
{
  PrecomputedEmissions &E=precomputedEmissions[parentChild];

  // First, allocate the array
  int numAlpha=alphabet.size();
  int L=taxon->getSeqLen();
  int numFC=FunctionalClass::numClasses();
  E.resize(L,numFC,numAlpha);

  // Iterate over all sequence positions in the child taxon's sequence
  Array1D<float> temp(numAlpha);
  for(int pos=0 ; pos<L ; ++pos) {
    Array3D<float>::IndexedOnce<float> slice=E[pos];
    for(int fc=0 ; fc<numFC ; ++fc) {
      Array3D<float>::IndexedTwice<float> row=slice[fc];
      switch(insideOutside) {
      case INSIDE:  
	F.precomputeInside(temp,*taxon,pos,fc); 
	break;
      case OUTSIDE: 
	F.precomputeOutside(temp,*taxon,pos,*avoid,fc); 
	break;
      }
      for(Symbol s=0 ; s<numAlpha ; ++s) row[s]=temp[s];
    }
  }
}



