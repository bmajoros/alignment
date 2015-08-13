/****************************************************************
 PosteriorBackward.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "PosteriorBackward.H"
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
                      PosteriorBackward methods
 ****************************************************************/

PosteriorBackward::PosteriorBackward(const BranchHMM &hmm,Taxon &parent,
				     Taxon &child,int numTaxa,
				     Vector<int> &insideLeaves,
				     Vector<int> &outsideLeaves,
				     CollapsedOrthologyMatrix &childCOM,
				     CollapsedOrthologyMatrix &parentCOM,
				     Vector<bool> &parentReachability,
				     Vector<bool> &childReachability,
				     Array2D<SparseMatrix3D*> &matrices,
				     double logPseudocount,
				     Array1D<Taxon> &taxa,
				     Phylogeny *tree,
				     bool wantTransitivity,
				     float indelCoef)
  : hmm(hmm), parent(parent), child(child), numTaxa(numTaxa), 
    numStates(hmm.getNumStates()), insideLeaves(insideLeaves),
    outsideLeaves(outsideLeaves), childCOM(childCOM), parentCOM(parentCOM),
    parentReachability(parentReachability), matrices(matrices),
    childReachability(childReachability), logPseudocount(logPseudocount),
    taxa(taxa), tree(tree), wantTransitivity(wantTransitivity),
    indelCoef(indelCoef)
{
  int m=parent.getSeqLen(), n=child.getSeqLen();
  cachedEmitP.resize(numStates,m+1,n+1);
  cachedEmitP.setAllTo(log0);
  numInsideLeaves=insideLeaves.size();
  numOutsideLeaves=outsideLeaves.size();
  branch=child.getBranchToParent();
  upMap=&branch->getUpMap();
  downMap=&branch->getDownMap();
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  nI=QI.size(); nD=QD.size(); nM=QM.size();
  fillMatrix();
}



double PosteriorBackward::getLikelihood() const
{
  return B(0,0,0);
}



double PosteriorBackward::operator()(int i,int j,STATE k) const
{
  return B(i,j,k);
}



int PosteriorBackward::getFirstDim() const
{
  return B.getFirstDim();
}



int PosteriorBackward::getSecondDim() const
{
  return B.getSecondDim();
}



int PosteriorBackward::getThirdDim() const
{
  return B.getThirdDim();
}



float PosteriorBackward::getCachedEmitP(STATE k,int col1,int col2)
{
  return cachedEmitP(k,col1,col2);
}



double PosteriorBackward::getEmitP2(STATE q,int newI,int newJ)
{
  //cout<<"getEmitP"<<endl;
  float &cachedValue=cachedEmitP(q,newI,newJ);
  if(finite(cachedValue)) return cachedValue;
  bool justPseudo;
  PHMM_StateType stateType=hmm.getStateType(q);
  int parentL=parentCOM.getSeqLen(), childL=childCOM.getSeqLen();
  Vector<float> sumV;
  cout<<"a"<<endl;
  for(int il=0 ; il<numInsideLeaves ; ++il) {
    int insideLeaf=insideLeaves[il];
    CollapsedOrthologyMatrix::Entry ie=childCOM(newJ,insideLeaf);
    bool iDirect=ie.isDirectMatch();
    int iBegin=ie.begin, iEnd=ie.end;
    if(iEnd>ie.begin) --iEnd; // because for indel states it's noninclusive
    if(iEnd>taxa[insideLeaf].getSeqLen()) --iEnd;
    if(iBegin<0) iBegin=0;
    for(int ol=0 ; ol<numOutsideLeaves ; ++ol) {
      int outsideLeaf=outsideLeaves[ol];
      SparseMatrix3D &M=*matrices[outsideLeaf][insideLeaf]; 
      float pathLen=tree->distanceBetweenLeaves(taxa[insideLeaf].getID(),
						taxa[outsideLeaf].getID());
      CollapsedOrthologyMatrix::Entry oe=parentCOM(newI,outsideLeaf);
      bool oDirect=oe.isDirectMatch();
      Set<PHMM_StateType> stateTypes;
      bool DEL=(stateType==PHMM_DELETE || ie.end>ie.begin || 
		ie.begin<0 || ie.end>taxa[insideLeaf].getSeqLen());
      bool INS=(stateType==PHMM_INSERT || oe.end>oe.begin || 
		oe.begin<0 || oe.end>taxa[outsideLeaf].getSeqLen());
      if(!INS && !DEL) stateTypes.insert(PHMM_MATCH);
      else {
	if(INS) stateTypes.insert(PHMM_INSERT);
	if(DEL) stateTypes.insert(PHMM_DELETE);
      }
      int oBegin=oe.begin, oEnd=oe.end;
      if(oEnd>oe.begin) --oEnd; // because for indel states it's noninclusive
      if(oEnd>taxa[outsideLeaf].getSeqLen()) --oEnd;
      if(oBegin<0) oBegin=0;
      Vector<float> sumW;
      for(int opos=oBegin ; opos<=oEnd ; ++opos) {
	Set<PHMM_StateType>::iterator cur=stateTypes.begin(), 
	  end=stateTypes.end();
	cout<<"a"<<endl;
	for(; cur!=end ; ++cur) {
	  PHMM_StateType stateType=*cur;
	  EntryList &row=M(opos,stateType);
	  //EntryList::iterator cur=row.begin(), end=row.end();
	  cout<<"b"<<endl;
	  EntryList::iterator end=row.end();
	  EntryList::iterator &saved=M.getIter(opos,stateType);
	  for(; saved!=end ; ++saved) {
	    SparseMatrix3D::Entry &e=*saved;
	    if(e.y>=iBegin) break;
	  }
	  cout<<"c"<<endl;
	  EntryList::iterator cur=saved;
	  for(; cur!=end ; ++cur) {
	    SparseMatrix3D::Entry &e=*cur;
	    //cout<<"i "<<&e<<endl;
	    if(e.y>iEnd) break;
	    //cout<<"j"<<endl;
	    float value=e.value;
	    //cout<<"k"<<endl;
	    if(stateType!=PHMM_MATCH) value+=log(indelCoef); // ###
	    cout<<"l "<<sumW.size()<<endl;
	    //if(stateType!=PHMM_MATCH) value=logPseudocount;//###
	    sumW.push_back(value);
	    cout<<"m "<<e.y<<" "<<value<<endl;
	  }
	  cout<<"n"<<endl;
	}
	cout<<"o"<<endl;
      }
      float theSum=sumLogProbs(sumW);
      if(isFinite(theSum)) 
	sumV.push_back(theSum-log(pathLen));
      //else if(wantTransitivity)
      //if(stateType==PHMM_MATCH && iDirect && oDirect)
      //{justPseudo=true; return cachedValue=log0;}

    }
  }
  cout<<"c"<<endl;
  float ave=sumLogProbs(sumV)-log(sumV.size());
  if(!isFinite(ave)) return cachedValue=logPseudocount;
  return cachedValue=ave;
}



double PosteriorBackward::getEmitP(STATE q,int newI,int newJ)
{
  //cout<<"getEmitP"<<endl;
  float &cachedValue=cachedEmitP(q,newI,newJ);
  if(finite(cachedValue)) return cachedValue;
  bool justPseudo;
  PHMM_StateType stateType=hmm.getStateType(q);
  int parentL=parentCOM.getSeqLen(), childL=childCOM.getSeqLen();
  Vector<float> sumV;
  for(int il=0 ; il<numInsideLeaves ; ++il) {
    int insideLeaf=insideLeaves[il];
    CollapsedOrthologyMatrix::Entry ie=childCOM(newJ,insideLeaf);
    bool iDirect=ie.isDirectMatch();
    int iBegin=ie.begin, iEnd=ie.end;
    if(iEnd>ie.begin) --iEnd; // because for indel states it's noninclusive
    if(iEnd>taxa[insideLeaf].getSeqLen()) --iEnd;
    if(iBegin<0) iBegin=0;
    for(int ol=0 ; ol<numOutsideLeaves ; ++ol) {
      int outsideLeaf=outsideLeaves[ol];
      SparseMatrix3D &M=*matrices[outsideLeaf][insideLeaf]; 
      float pathLen=tree->distanceBetweenLeaves(taxa[insideLeaf].getID(),
						taxa[outsideLeaf].getID());
      CollapsedOrthologyMatrix::Entry oe=parentCOM(newI,outsideLeaf);
      bool oDirect=oe.isDirectMatch();
      Set<PHMM_StateType> stateTypes;
      bool DEL=(stateType==PHMM_DELETE || ie.end>ie.begin || 
		ie.begin<0 || ie.end>taxa[insideLeaf].getSeqLen());
      bool INS=(stateType==PHMM_INSERT || oe.end>oe.begin || 
		oe.begin<0 || oe.end>taxa[outsideLeaf].getSeqLen());
      if(!INS && !DEL) stateTypes.insert(PHMM_MATCH);
      else {
	if(INS) stateTypes.insert(PHMM_INSERT);
	if(DEL) stateTypes.insert(PHMM_DELETE);
      }
      int oBegin=oe.begin, oEnd=oe.end;
      if(oEnd>oe.begin) --oEnd; // because for indel states it's noninclusive
      if(oEnd>taxa[outsideLeaf].getSeqLen()) --oEnd;
      if(oBegin<0) oBegin=0;
      Vector<float> sumW;
      for(int opos=oBegin ; opos<=oEnd ; ++opos) {
	Set<PHMM_StateType>::iterator cur=stateTypes.begin(), 
	  end=stateTypes.end();
	for(; cur!=end ; ++cur) {
	  PHMM_StateType stateType=*cur;
	  EntryList &row=M(opos,stateType);
	  //EntryList::iterator cur=row.begin(), end=row.end();
	  EntryList::iterator cur=row.begin(), end=row.end();
	  for(; cur!=end ; ++cur) {
	    SparseMatrix3D::Entry e=*cur;
	    if(e.y<iBegin) continue;
	    if(e.y>iEnd) break;
	    float value=e.value;
	    if(stateType!=PHMM_MATCH) value+=log(indelCoef); // ###
	    //if(stateType!=PHMM_MATCH) value=logPseudocount;//###
	    sumW.push_back(value);
	  }
	}
      }
      float theSum=sumLogProbs(sumW);
      if(isFinite(theSum)) 
	sumV.push_back(theSum-log(pathLen));
      //else if(wantTransitivity)
	//if(stateType==PHMM_MATCH && iDirect && oDirect)
	  //{justPseudo=true; return cachedValue=log0;}

    }
  }
  float ave=sumLogProbs(sumV)-log(sumV.size());
  if(!isFinite(ave)) return cachedValue=logPseudocount;
  return cachedValue=ave;
}




void PosteriorBackward::resetSavedIterators()
{
  for(int il=0 ; il<numInsideLeaves ; ++il) {
    int insideLeaf=insideLeaves[il];
    for(int ol=0 ; ol<numOutsideLeaves ; ++ol) {
      int outsideLeaf=outsideLeaves[ol];
      SparseMatrix3D &M=*matrices[outsideLeaf][insideLeaf]; 
      int L=M.getFirstDim();
      for(int x=0 ; x<L ; ++x)
	for(int stateType=0 ; stateType<3 ; ++stateType)
	  M.getIter(x,stateType)=M.begin(x,stateType);
    }
  }
}



void PosteriorBackward::fillMatrix()
{
  // Initialization:
  //cout<<"B::init"<<endl;
  
  int m=parent.getSeqLen(), n=child.getSeqLen();
  B.resize(m+1,n+1,numStates);
  B.setAllTo(NEGATIVE_INFINITY);
  for(STATE k=1 ; k<numStates ; ++k) B(m,n,k)=LOG_1;
  Array1D<double> logProbs(nI);
  //resetSavedIterators(m);
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
    //resetSavedIterators(i);
    for(STATE k=1 ; k<numStates ; ++k) {
      for(int ih=0 ; ih<nD ; ++ih) {
	STATE h=QD[ih];
	logProbs[ih]=safeAdd(B(i+1,n,h),getEmitP(h,i,n));
      }
      B(i,n,k)=sumLogProbs<double>(logProbs);
    }
  }

  // Recurrence:
  //cout<<"B::recur"<<endl;
  logProbs.resize(numStates);
  logProbs[0]=NEGATIVE_INFINITY;
  for(int i=m-1 ; i>=0 ; --i) {
    //cout<<"resetting saved iterators"<<endl;
    //resetSavedIterators();
    //cout<<"done resetting"<<endl;
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
	      //cout<<"state: "<<h<<" statetype="<<hmm.getStateType(h)<<endl;
	      //throw "error in PosteriorBackward::fillMatrix (recurrence)";
	      break;
	    }
	  ++p;
	}
	B(i,j,k)=sumLogProbs<double>(logProbs);
      }
    }
  }

  // Termination:
  //cout<<"B::Term"<<endl;
  for(STATE h=1 ; h<numStates ; ++h) {
    switch(hmm.getStateType(h)) 
      {
      case PHMM_MATCH:
	//resetSavedIterators(0);
	logProbs[h]=safeAdd(getEmitP(h,0,0),B(1,1,h));
	break;
      case PHMM_INSERT:
	//resetSavedIterators(m);
	logProbs[h]=safeAdd(getEmitP(h,m,0),B(0,1,h));
	break;
      case PHMM_DELETE:
	//resetSavedIterators(0);
	logProbs[h]=safeAdd(getEmitP(h,0,n),B(1,0,h));
	break;
      default: throw "error in PosteriorBackward::fillMatrix (termination)";
      }
  }
  B(0,0,0)=sumLogProbs<double>(logProbs);
  //cout<<"/B"<<endl;
}

