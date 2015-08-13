/****************************************************************
 BranchHMM.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BranchHMM.H"
#include "BOOM/AlphabetMap.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;



BranchHMM::BranchHMM(const PairHMM &hmm)
  : PairHMM(hmm), matrices(hmm.getNumStates()), eqFreqs(hmm.getNumStates()),
    parentFC(hmm.getNumStates()), childFC(hmm.getNumStates()),
    matrixByClass(FunctionalClass::numClasses(),FunctionalClass::numClasses()),
    gainFactor(0), lossFactor(0), retentionFactor(0), voidFactor(0),
    scoreGainAndLoss(hmm.getNumStates())
{
  scoreGainAndLoss.setAllTo(true);
  initMatrices();
}



BranchHMM::BranchHMM(const PairHMM &hmm,Array1D<SubstitutionMatrix*> &A)
  : PairHMM(hmm), matrices(A), eqFreqs(hmm.getNumStates()),
    parentFC(hmm.getNumStates()), childFC(hmm.getNumStates()),
    matrixByClass(FunctionalClass::numClasses(),FunctionalClass::numClasses()),
    gainFactor(LOG_1), lossFactor(LOG_1), retentionFactor(LOG_1), 
    voidFactor(LOG_0), scoreGainAndLoss(hmm.getNumStates())
{
  // ctor

  scoreGainAndLoss.setAllTo(true);
  matrixByClass.setAllTo(NULL);
  int n=getNumStates();
  bool bgInit=false;
  FunctionalClass bg=FunctionalClass::getBackground();
  for(STATE i=0 ; i<n ; ++i) // ### should skip state 0, for speed
    if(matrices[i]) {
      SubstitutionMatrix *Pt=matrices[i];
      Pt->getEqFreqs(eqFreqs[i]);
      if(!bgInit && getFunctionalClass(i,PARENT)==bg &&
	 getFunctionalClass(i,CHILD)==bg) {
	bgEqFreqs=eqFreqs[i];
	bgInit=true;
      }	
    }
}



BranchHMM::~BranchHMM()
{
  //cout<<"------"<<endl;
  int n=matrices.size();
  for(int i=0 ; i<n ; ++i) {
    //cout<<"deleting "<< &matrices[i]<<endl;
    delete matrices[i];
  }
}



void BranchHMM::initMatrices()
{
  Transducer T(*this);
  int n=getNumStates();
  for(STATE i=0 ; i<n ; ++i) { // ### should skip state 0, for speed
    matrices[i]=buildMatrix(T,i);
    matrices[i]->getEqFreqs(eqFreqs[i]);
  }
}


SubstitutionMatrix *BranchHMM::buildMatrix(Transducer &T,STATE q)
{
  const Alphabet &alphabet=getAlphabet();
  AlphabetMap alphabetMap(alphabet,PureDnaAlphabet::global());
  SubstitutionMatrix &M=*new SubstitutionMatrix(alphabet,alphabetMap);
  int nAlpha=alphabet.size();
  for(Symbol s=0 ; s<nAlpha ; ++s)
    for(Symbol t=0 ; t<nAlpha ; ++t)
      M(s,t)=T.getEmitP(q,s,t);
  return &M;
}



void BranchHMM::convertToLogs() 
{
  PairHMM::convertToLogs();
  int n=matrices.size();
  for(int i=0 ; i<n ; ++i) {
    SubstitutionMatrix *M=matrices[i];
    if(M) M->convertToLogs();
  }
}



BranchHMM *BranchHMM::eliminateSilentStates() const
{
  Array1D<STATE> stateMap;
  PairHMM *phmm=PairHMM::eliminateSilentStates(stateMap);
  int oldNumStates=getNumStates(), newNumStates=phmm->getNumStates();
  Array1D<SubstitutionMatrix*> newM(newNumStates);
  for(int i=0 ; i<newNumStates ; ++i) {
    SubstitutionMatrix *oldMatrix=matrices[stateMap[i]];
    newM[i]=oldMatrix ? new SubstitutionMatrix(*oldMatrix) : NULL;
  }
  BranchHMM *bhmm=new BranchHMM(*phmm,newM);
  delete phmm;
  return bhmm;
}



void BranchHMM::printOn(ostream &os) const 
{
  int n=getNumStates();
  os<<"functional classes: "<<n<<endl;
  for(int i=0 ; i<n ; ++i) {
    FunctionalClass fc=getFunctionalClass(i,PARENT);
    os<<"state #"<<i<<" : (PARENT) "<<fc.getLabel()<<" = "<<fc.getName()<<endl;
    fc=getFunctionalClass(i,CHILD);
    os<<"state #"<<i<<" : (CHILD) "<<fc.getLabel()<<" = "<<fc.getName()<<endl;
  }
  os<<"matrices:"<<endl;
  n=matrices.size();
  for(int i=0 ; i<n ; ++i)
    if(matrices[i])
      os<<"matrix "<<i<<":\n"<<*matrices[i]<<endl;
  n=eqFreqs.size();
  os<<"equilibrium frequencies:"<<endl;
  for(int i=0 ; i<n ; ++i)
    os<<"eq freqs for state "<<i<<":\n"<<eqFreqs[i]<<endl;
  static_cast<const PairHMM*>(this)->printOn(os);
}



ostream &operator<<(ostream &os,const BranchHMM &hmm)
{
  hmm.printOn(os);
  return os;
}



void BranchHMM::setFunctionalClass(STATE q,BranchEnd e,FunctionalClass fc)
{
  if(e==PARENT) parentFC[q]=fc;
  else childFC[q]=fc;
  //matrixByClass[fc]=matrices[q]; // ###
}



void BranchHMM::initFuncStateTypeMatrix()
{
  FunctionalClass bg=FunctionalClass::getBackground();
  const int numStateTypes=NUM_PHMM_STATETYPES;
  const int numFuncClasses=FunctionalClass::numClasses();
  classAndTypeToState.resize(numFuncClasses,numFuncClasses,numStateTypes);
  classAndTypeToState.setAllTo(INVALID_STATE);
  classEndTypeToState.resize(numFuncClasses,2,3);
  gainLossStates.resize(numFuncClasses,2,numStateTypes);
  retentionStates.resize(numFuncClasses,numStateTypes);
  const int numStates=getNumStates();
  for(STATE q=0 ; q<numStates ; ++q) {
    FunctionalClass fcParent=getFunctionalClass(q,PARENT);
    FunctionalClass fcChild=getFunctionalClass(q,CHILD);
    if(fcParent.isValid() && fcChild.isValid()) {
      PHMM_StateType stateType=getStateType(q);
      classAndTypeToState(fcParent,fcChild,stateType)=q;
      if(stateType<3) {
	classEndTypeToState(fcParent,PARENT,stateType).push_back(q);
	classEndTypeToState(fcChild,CHILD,stateType).push_back(q);
      }
      matrixByClass(fcParent,fcChild)=matrices[q];
      if(fcParent==fcChild) {
	if(fcChild==bg)
	  backgroundStates.push_back(q);
	else retentionStates(fcParent,stateType).push_back(q);
      }
      else {
	if(fcParent==bg) 
	  gainLossStates(fcChild,CHILD,stateType).push_back(q);
	else 
	  gainLossStates(fcParent,PARENT,stateType).push_back(q);
      }
    }
  }
}



void BranchHMM::initSubstGenerators()
{
  int numStates=getNumStates();
  for(int i=0 ; i<numStates ; ++i) {
    SubstitutionMatrix *Pt=getSubstMatrix(i);
    if(Pt) Pt->initGenerators();
  }
}



void BranchHMM::setGainLossFactor(GainLossType t,double f)
{
  switch(t)
    {
    case GLT_GAIN: gainFactor=log(f); voidFactor=log(1-f); break;
    case GLT_LOSS: lossFactor=log(f); break;
    case GLT_RETENTION: retentionFactor=log(f); break;
    case GLT_VOID: throw "BranchHMM::setGainLossFactor(GLT_VOID)";
    }

  // ### DEBUGGING
  //voidFactor=retentionFactor=0;
  // ### /DEBUGGING

}



void BranchHMM::setGainLossScorePolicy(STATE q,bool b)
{
  scoreGainAndLoss[q]=b;
}





