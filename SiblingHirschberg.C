/****************************************************************
 SiblingHirschberg.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "SiblingHirschberg.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Constants.H"
#include "BranchAttributes.H"
using namespace std;
using namespace BOOM;


SiblingHirschberg::SiblingHirschberg(BranchHMM *hmm,Taxon *left,Taxon *right,
				     LinkFelsenstein &F,int bandwidth,
				     bool usePrescan,int maxThreads,
				     ContentSensor *contentSensor)
  : Hirschberg(left->getParent(),NULL,F,bandwidth,left->getSeqLen(),
	       right->getSeqLen(),usePrescan,maxThreads,contentSensor),
    left(left), right(right)
{
  this->hmm=hmm;
  numStates=hmm->getNumStates();
  initStateSets();
  m=left->getSeqLen();
  n=right->getSeqLen();
  leftBranch=parent->getBranchToChild(*left);
  rightBranch=parent->getBranchToChild(*right);
  leftDownMap=&leftBranch->getDownMap();
  leftUpMap=&leftBranch->getUpMap();
  rightDownMap=&rightBranch->getDownMap();
  rightUpMap=&rightBranch->getUpMap();
  parent->setSeqLen(1);
  leftDownMap->resize(1);
  rightDownMap->resize(1);
  leftDownMapEntry=&(*leftDownMap)[0];
  rightDownMapEntry=&(*rightDownMap)[0];
  leftUpMap->resize(m);
  rightUpMap->resize(n);
  leftUpMap->setAllTo(0);
  rightUpMap->setAllTo(0);
  F.setFocalBranches(leftBranch,rightBranch); // ###
  precomputeEmissions();
}



double SiblingHirschberg::getEmitP(STATE k,int p,int c,PairHMMStateType t)
{
  //if(posteriors) return posteriors->compute(p,c,k);
  return emission(k,p,c,t);
}



double SiblingHirschberg::emission(STATE k,int leftIndex,int rightIndex,
				   PairHMMStateType stateType)
{
  INTERNAL_ERROR; // is this called?


  FunctionalClass fcLeft=hmm->getFunctionalClass(k,PARENT);
  FunctionalClass fcRight=hmm->getFunctionalClass(k,CHILD);
  bool useGainLossScores=hmm->shouldScoreGainLoss(k);
  Array1D<float> V(numAlpha), leftV(numAlpha), rightV(numAlpha);
  double P;
  switch(stateType)
    {
    case PHMM_MATCH: {
      PrecomputedEmissions &peLeft=precomputedEmissions[LEFT];
      PrecomputedEmissions &peRight=precomputedEmissions[RIGHT];
      Array3D<float>::IndexedTwice<float> leftRow=
	peLeft[leftIndex][fcLeft.getClassID()];
      Array3D<float>::IndexedTwice<float> rightRow=
	peRight[rightIndex][fcRight.getClassID()];
      FunctionalClass fcParent=pickForeground(fcLeft,fcRight);
      const Array1D<double> &eqFreqs=fcParent.getEqFreqs();
      SubstitutionMatrix &leftPt=
	*leftBranch->getHMM()->getSubstMatrix(fcParent,fcLeft);
      SubstitutionMatrix &rightPt=
	*rightBranch->getHMM()->getSubstMatrix(fcParent,fcRight);
      for(Symbol i=0 ; i<numAlpha; ++i) {
	for(Symbol j=0 ; j<numAlpha; ++j) {
	  leftV[j]=safeAdd(leftPt(i,j),leftRow[j]); // ### leftRow[i]
	  rightV[j]=safeAdd(rightPt(i,j),rightRow[j]);
	}
	double leftP=sumLogProbs<float>(leftV);
	double rightP=sumLogProbs<float>(rightV);
	V[i]=safeAdd(log(eqFreqs[i]),leftP,rightP);
      }
      P=sumLogProbs<float>(V);
    }
      break;
    case PHMM_DELETE: {
      Array3D<float>::IndexedTwice<float> row=
	precomputedEmissions[LEFT][leftIndex][fcLeft.getClassID()];
      const Array1D<double> &eqFreqs=fcLeft.getEqFreqs();
      for(Symbol i=0 ; i<numAlpha; ++i)
	V[i]=safeAdd(log(eqFreqs[i]),row[i]);
      P=sumLogProbs<float>(V);
    }
      break;
    case PHMM_INSERT:{
      Array3D<float>::IndexedTwice<float> row=
	precomputedEmissions[RIGHT][rightIndex][fcRight.getClassID()];
      const Array1D<double> &eqFreqs=fcRight.getEqFreqs();
      for(Symbol i=0 ; i<numAlpha; ++i)
	V[i]=safeAdd(log(eqFreqs[i]),row[i]);
      P=sumLogProbs<float>(V);
    }
      break;
    default: throw "error in SiblingViterbi::emission()";
    }
  return P;
}



void SiblingHirschberg::precomputeEmissions()
{
  precomputeEmissions(precomputedEmissions[LEFT],*left);
  precomputeEmissions(precomputedEmissions[RIGHT],*right);
}



void SiblingHirschberg::precomputeEmissions(PrecomputedEmissions &E,
					    Taxon &cladeRoot)
{
  int L=cladeRoot.getSeqLen();
  if(L==0) INTERNAL_ERROR;
  int numFC=FunctionalClass::numClasses();
  cladeRoot.getPrecomputedEmissions().resize(L,numFC,numAlpha);
  int numChildren=cladeRoot.getNodeType()==ROOT_NODE ? 1 : 
    (cladeRoot.getNodeType()==INTERNAL_NODE ? 2 : 0);
  PhylogenyNode *parent=cladeRoot.getNode();
  /*
  for(int i=0 ; i<numChildren ; ++i) {
    PhylogenyNode *child=parent->getChild(i);
    Taxon &taxon=static_cast<Taxon&>(*child->getDecoration());
    taxon.getPrecomputedEmissions().resize(taxon.getSeqLen(),numFC,numAlpha);
  }
  */
  F.precomputeInside(cladeRoot);
  E=cladeRoot.getPrecomputedEmissions();
}

