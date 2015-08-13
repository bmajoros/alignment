/****************************************************************
 HirschPass.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "HirschPass.H"
#include "BOOM/Constants.H"
#include "BOOM/Exceptions.H"
#include "BOOM/SumLogProbs.H"
#include "Hirschberg.H"
using namespace std;
using namespace BOOM;

HirschPass::HirschPass(int minX,int maxX,int minY,int maxY,
		       const BandingPattern &bp,BranchHMM &hmm,
		       const Vector<STATE> &QI,int numI,
		       const Vector<STATE> &QD,int numD,
		       const Vector<STATE> &QM,int numM,
		       HirschbergDirection dir,ViterbiConstraint vc,
		       PrecomputedEmissions *precomputedEmissions,
		       FunctionalParse *parentFP, FunctionalParse *childFP,
		       IndexMap *upMap,IndexMap *downMap,
		       Hirschberg *hirschberg,ContentSensor *contentSensor)
  : hmm(hmm), bandingPattern(bp), minX(minX), maxX(maxX), 
    minY(minY), maxY(maxY), numStates(hmm.getNumStates()),
    frame(minX,maxX,minY,maxY,hmm.getNumStates(),bp,dir,hmm), vc(vc),
    QI(QI), QD(QD), QM(QM), numI(numI), numD(numD), numM(numM),
    precomputedEmissions(precomputedEmissions),
    parentFP(parentFP), childFP(childFP), upMap(upMap), downMap(downMap),
    numAlpha(4), hirschberg(hirschberg), bg(FunctionalClass::getBackground()),
    firstTaxon(hirschberg ? hirschberg->getFirstTaxon(): NULL),
    secondTaxon(hirschberg ? hirschberg->getSecondTaxon() : NULL),
    contentSensor(contentSensor), 
    constraintInterval(INT_MIN,INT_MAX)//minX-1,maxX+1)
{
  // ctor
}



void HirschPass::setConstraintInterval(int from,int to)
{
  constraintInterval.first=from;
  constraintInterval.second=to;
}



void HirschPass::resetConstraintInterval()
{
  constraintInterval.first=minX-1;
  constraintInterval.second=maxX+1;
}



double HirschPass::getEmitP(STATE k,int p,int c,PairHMMStateType t)
{
  //###if(posteriors) return posteriors->compute(p,c,k);
  //if(c<0 || c>=constraintInterval.first && c<=constraintInterval.second)
  switch(vc)
    {
    case VC_UNCONSTRAINED: 
      break;
    case VC_PARENT_PARSE:
      if(p<constraintInterval.first || p>constraintInterval.second) break;
      if(emitsParent(t) && hmm.getFunctionalClass(k,PARENT)!=(*parentFP)[p])
	return LOG_0 ;
      break;
    case VC_CHILD_PARSE:
      if(c<constraintInterval.first || c>constraintInterval.second) break;
      if(emitsChild(t) && hmm.getFunctionalClass(k,CHILD)!=(*childFP)[c])
	return LOG_0;
      break;
    case VC_PARENT_AND_CHILD_PARSE:
      if(emitsParent(t) && hmm.getFunctionalClass(k,PARENT)!=(*parentFP)[p] 
	 ||
	 emitsChild(t) && hmm.getFunctionalClass(k,CHILD)!=(*childFP)[c])
	return LOG_0;
      break;
    case VC_INDEL_HISTORY: // ### ?
      if(emitsParent(t) && (*downMap)[p]!=c ||
	 emitsChild(t) && (*upMap)[c]!=p)
	return LOG_0;
      break;
    case VC_INDEL_AND_PARENT_PARSE:
      if(emitsParent(t)) {
	if((*downMap)[p]!=c) {
	  return LOG_0;
	}
	if(hmm.getFunctionalClass(k,PARENT)!=(*parentFP)[p]) {
	  return LOG_0;
	}
      }
      if(emitsChild(t) && (*upMap)[c]!=p) {
	return LOG_0;
      }
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



double HirschPass::emission(STATE k,int parentResidueIndex,
			     int childResidueIndex,
			     PairHMMStateType stateType)
{
  //const double eqFactor=1.0;
  //const double consFactor=1.0;
  Array1D<float> V(numAlpha), V2(numAlpha);
  switch(stateType)
    {
    case PHMM_MATCH:{
      if(parentResidueIndex<0 || childResidueIndex<0) INTERNAL_ERROR; // ###
      FunctionalClass fcParent=hmm.getFunctionalClass(k,PARENT);
      FunctionalClass fcChild=hmm.getFunctionalClass(k,CHILD);
      const Array1D<double> &eqFreqs=hmm.getEqFreqs(k);
      Array3D<float>::IndexedTwice<float> childRow=
	precomputedEmissions[CHILD][childResidueIndex][fcChild];
      Array3D<float>::IndexedTwice<float> parentRow=
	precomputedEmissions[PARENT][parentResidueIndex][fcParent];
      SubstitutionMatrix &Pt=*hmm.getSubstMatrix(k);
      for(Symbol sParent=0 ; sParent<numAlpha ; ++sParent) {
	double eq=getEQ(eqFreqs,sParent,firstTaxon,parentResidueIndex,
			hmm.getFunctionalClass(k,PARENT));
	//eq*=eqFactor;//###
	for(Symbol sChild=0 ; sChild<numAlpha ; ++sChild) {
	  V[sChild]=safeAdd(eq,/*consFactor**/safeAdd(parentRow[sParent]
				    ,Pt(sParent,sChild),childRow[sChild]));
	}	
	V2[sParent]=sumLogProbs<float>(V);
      }
      return sumLogProbs<float>(V2); }
    case PHMM_INSERT:{
      if(childResidueIndex<0) INTERNAL_ERROR; // ###
      FunctionalClass fc=hmm.getFunctionalClass(k,CHILD);
      Array1D<double> &eqFreqs=fc.getEqFreqs(); // ###
      Array3D<float>::IndexedTwice<float> row=
	precomputedEmissions[CHILD][childResidueIndex][fc];
      for(Symbol s=0 ; s<numAlpha ; ++s) {
	double eq=getEQ(eqFreqs,s,secondTaxon,childResidueIndex,
			hmm.getFunctionalClass(k,CHILD));
	V[s]=safeAdd(/*consFactor**/row[s],/*eqFactor**/eq);
      }
      return sumLogProbs<float>(V); }
    case PHMM_DELETE:{
      if(parentResidueIndex<0) INTERNAL_ERROR; // ###
      FunctionalClass fc=hmm.getFunctionalClass(k,PARENT);
      Array1D<double> &eqFreqs=fc.getEqFreqs();
      Array3D<float>::IndexedTwice<float> row=
	precomputedEmissions[PARENT][parentResidueIndex][fc];
      for(Symbol s=0 ; s<numAlpha ; ++s) {
	double eq=getEQ(eqFreqs,s,firstTaxon,parentResidueIndex,
			hmm.getFunctionalClass(k,PARENT));
	V[s]=safeAdd(/*consFactor**/row[s],/*eqFactor**/eq);
      }
      return sumLogProbs<float>(V); }
    default:
      INTERNAL_ERROR;
    }
}



bool HirschPass::columnIsEmpty(int x,int minY,int maxY)
{
  /*
  {//###
    int numNonzero=0;
    int numStates=frame.getNumStates();
    WindowColumn &col=frame.getCol(x);
    for(int y=minY ; y<=maxY ; ++y) {
      Array2D<float>::RowIn2DArray<float> pillar=col[y];
      for(STATE q=0 ; q<numStates ; ++q) {
	if(isFinite(pillar[q])) ++numNonzero;
      }
    }
  }//###
  */

  int numStates=frame.getNumStates();
  WindowColumn &col=frame.getCol(x);
  for(int y=minY ; y<=maxY ; ++y) {
    Array2D<float>::RowIn2DArray<float> pillar=col[y];
    for(STATE q=0 ; q<numStates ; ++q) {
      if(isFinite(pillar[q])) return false;
    }
  }
  return true;
}



void HirschPass::checkColumn(int x,int minY,int maxY)
{
  if(columnIsEmpty(x,minY,maxY)) INTERNAL_ERROR;
}



void HirschPass::getPrescanStates(int firstIndex,int secondIndex,
				  Vector<STATE> &dest)
{
  // Append state 0 if applicable
  if(firstIndex<0 && secondIndex<0) {
    dest.push_back(0); // ### only applicable to Backward algorithm?
    return;
  }

  // First, add background-model states:
  dest.append(hmm.getBackgroundStates());
  
  // Now add states implied by parent & child foreground classes
  if(firstIndex<0) {
    //if(secondIndex<0) return;
    getCrossFuncStates(secondTaxon->getFcConstraints()[secondIndex],CHILD,
		       dest);
    return;
  }
  if(secondIndex<0) {
    firstTaxon->getFcConstraints();
    firstTaxon->getFcConstraints()[firstIndex];
    getCrossFuncStates(firstTaxon->getFcConstraints()[firstIndex],PARENT,
		       dest);
    return;
  }
  FuncClassSet &parentClasses=firstTaxon->getFcConstraints()[firstIndex];
  FuncClassSet &childClasses=secondTaxon->getFcConstraints()[secondIndex];

  FuncClassSet::iterator pCur=parentClasses.begin(), pEnd=parentClasses.end();
  FuncClassSet::iterator cCur=childClasses.begin(), cEnd=childClasses.end();
  FunctionalClass p, c;
  //FuncClassSet parentFCs, childFCs, retentionFCs;
  FuncClassSet retentionFCs;
  while(pCur!=pEnd && cCur!=cEnd) {
    p=*pCur;
    c=*cCur;
    if(p<c) {
      //parentFCs.push_back(p);
      ++pCur;
    }
    else if(c<p) {
      //childFCs.push_back(c);
      ++cCur;
    }
    else {
      retentionFCs.push_back(p);
      ++pCur;
      ++cCur;
    }
  }
  //parentFCs.append(pCur,pEnd);
  //childFCs.append(cCur,cEnd);

  //getCrossFuncStates(parentFCs,PARENT,dest);
  //getCrossFuncStates(childFCs,CHILD,dest);
  getCrossFuncStates(parentClasses,PARENT,dest); // ###
  getCrossFuncStates(childClasses,CHILD,dest);   // ###
  getRetentionStates(retentionFCs,dest);
}



void HirschPass::getCrossFuncStates(FuncClassSet &classes,BranchEnd branchEnd,
				    Vector<STATE> &dest)
{
  List<FunctionalClass>::iterator cur=classes.begin(), end=classes.end();
  for(; cur!=end ; ++cur) {
    FunctionalClass fc=*cur;
    for(PHMM_StateType t=0 ; t<3 ; ++t) {
      dest.append(hmm.getGainLossStates(fc,branchEnd,t));
    }	
  }
}



void HirschPass::getRetentionStates(FuncClassSet &classes,Vector<STATE> &dest)
{
  List<FunctionalClass>::iterator cur=classes.begin(), end=classes.end();
  for(; cur!=end ; ++cur) {
    FunctionalClass fc=*cur;
    for(PHMM_StateType t=0 ; t<3 ; ++t)
      dest.append(hmm.getRetentionStates(fc,t));
  }
}



