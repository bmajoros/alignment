/****************************************************************
 SiblingViterbi.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "SiblingViterbi.H"
#include "BOOM/SumLogProbs.H"
using namespace std;
using namespace BOOM;


SiblingViterbi::SiblingViterbi(BranchHMM *hmm,Taxon *left,Taxon *right,
			       LinkFelsenstein &F,BandingType bandingType,
			       int bandwidth)
  : LinkViterbi(left->getParent(),NULL,F,bandingType,bandwidth),
    left(left), right(right)
{
  this->hmm=hmm;
  numStates=hmm->getNumStates();
  L1=left->getSeqLen();
  L2=right->getSeqLen();
  V.resize(L1+1,L2+1,numStates,bandwidth);
  T.resize(L1+1,L2+1,numStates,bandwidth);
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
  leftUpMap->resize(L1);
  rightUpMap->resize(L2);
  leftUpMap->setAllTo(0);
  rightUpMap->setAllTo(0);
  F.setFocalBranches(leftBranch,rightBranch); // ###
  precomputeEmissions();
}



double SiblingViterbi::getEmitP(STATE k,int p,int c,PairHMMStateType t)
{
  if(posteriors) return posteriors->compute(p,c,k);
  return emission(k,p,c,t);
}



double SiblingViterbi::emission(STATE k,int leftIndex,
				int rightIndex,
				PairHMMStateType stateType)
{
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
	  leftV[j]=safeAdd(leftPt(i,j),leftRow[j]);//### fixed 10/2/09
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


/*
double SiblingViterbi::OLD_emission(STATE k,int leftIndex,
				int rightIndex,
				PairHMMStateType stateType)
{
  FunctionalClass fcLeft=hmm->getFunctionalClass(k,PARENT);
  FunctionalClass fcRight=hmm->getFunctionalClass(k,CHILD);
  bool useGainLossScores=hmm->shouldScoreGainLoss(k);
  double P;
  switch(stateType)
    {
    case PHMM_MATCH:
      *leftDownMapEntry=leftIndex;
      *rightDownMapEntry=rightIndex;
      left->setFcConstraint(fcLeft);
      right->setFcConstraint(fcRight);
      P=static_cast<LossyFelsenstein&>(F).logLikelihood(0,*parent,
					     pickForeground(fcLeft,fcRight),
							useGainLossScores);
      left->removeFcConstraint();
      right->removeFcConstraint();
      break;
    case PHMM_INSERT: 
      P=F.logLikelihood(rightIndex,*right,fcRight,useGainLossScores);
      break;
    case PHMM_DELETE:
      P=F.logLikelihood(leftIndex,*left,fcLeft,useGainLossScores);
      break;
    default: throw "error in LinkBackward::getEmitP";
    }
  return P;
}
*/


void SiblingViterbi::precomputeEmissions()
{
  precomputeEmissions(precomputedEmissions[LEFT],*left);
  precomputeEmissions(precomputedEmissions[RIGHT],*right);
}



void SiblingViterbi::precomputeEmissions(PrecomputedEmissions &E,
					 Taxon &cladeRoot)
{
  // First, allocate the array
  //int numAlpha=alphabet.size();
  int L=cladeRoot.getSeqLen();
  int numFC=FunctionalClass::numClasses();
  E.resize(L,numFC,numAlpha);

  // Iterate over all sequence positions in the clade root
  Array1D<float> temp(numAlpha);
  for(int pos=0 ; pos<L ; ++pos) {
    Array3D<float>::IndexedOnce<float> slice=E[pos];
    for(int fc=0 ; fc<numFC ; ++fc) {
      Array3D<float>::IndexedTwice<float> row=slice[fc];
      F.precomputeInside(temp,cladeRoot,pos,fc);
      for(Symbol s=0 ; s<numAlpha ; ++s)
	row[s]=temp[s];
    }
  }
}

