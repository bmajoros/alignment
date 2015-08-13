/****************************************************************
 TransducerTemplate.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "TransducerTemplate.H"
#include "BOOM/Exceptions.H"
using namespace std;
using namespace BOOM;



/****************************************************************
                     StateTemplate methods
 ****************************************************************/

StateTemplate::StateTemplate(int stateID,PairHMMStateType type,
			     FunctionalClass parentFC,FunctionalClass childFC)
  : stateID(stateID), stateType(type), parentFC(parentFC), childFC(childFC)
{
  switch(type) {
  case PHMM_MATCH:
  case PHMM_INSERT:
  case PHMM_DELETE:
    scoreGainLossFels=true;
    break;
  default:
    scoreGainLossFels=false;
    break;
  }
}



void StateTemplate::changeStateID(int newID)
{
  stateID=newID;
}



int StateTemplate::getStateID() const
{
  return stateID;
}



PairHMMStateType StateTemplate::getStateType() const
{
  return stateType;
}



RateMatrix *StateTemplate::getMatrix(BranchEnd e) const
{
  FunctionalClass fc=(e==PARENT ? parentFC : childFC);
  /*
  if(!fc.isValid()) {
    cout<<"Functional class not valid in StateTemplate: "
	<<fc.getClassID()<<endl;
    INTERNAL_ERROR;
  }
  */
  return fc.getMatrix();
}



Phylogeny *StateTemplate::getPhylogeny(BranchEnd e) const
{
  FunctionalClass fc=(e==PARENT ? parentFC : childFC);
  /*
  if(!fc.isValid()) {
    cout<<"Functional class not valid in StateTemplate: "
	<<fc.getClassID()<<endl;
    INTERNAL_ERROR;
  }
  */
  return fc.getPhylogeny();
}



void StateTemplate::setFunctionalClass(FunctionalClass fc,BranchEnd e)
{
  if(e==PARENT) parentFC=fc;
  else childFC=fc;
}



FunctionalClass StateTemplate::getFunctionalClass(BranchEnd e)
{
  return e==PARENT ? parentFC : childFC;
}



bool StateTemplate::shouldScoreGainLoss() const
{
  return scoreGainLossFels;
}



void StateTemplate::disableGainLossScoring()
{
  scoreGainLossFels=false;
}



GainLossType StateTemplate::crossFunctionalType()
{
  return FunctionalClass::classifyGainLoss(parentFC,childFC);
}




/****************************************************************
                     TransitionTemplate methods
 ****************************************************************/
TransitionTemplate::TransitionTemplate(StateTemplate *from,
				       StateTemplate *to,
				       Lambda::Closure *f)
  : from(from), to(to), transProbFunction(f)
{
  // ctor
}



StateTemplate *TransitionTemplate::getFrom() const
{
  return from;
}



StateTemplate *TransitionTemplate::getTo() const
{
  return to;
}



Lambda::Closure *TransitionTemplate::getTransProbFunc() const
{
  return transProbFunction;
}



/****************************************************************
                     TransducerTemplate methods
 ****************************************************************/

TransducerTemplate::TransducerTemplate()
  : gainFactor(NULL), lossFactor(NULL), retentionFactor(NULL)
{
  // ctor
}



TransducerTemplate::~TransducerTemplate()
{
  int n=states.size();
  Vector<StateTemplate*>::iterator sCur=states.begin(), sEnd=states.end();
  ++sCur; // Skip q0
  for(; sCur!=sEnd ; ++sCur) delete *sCur;

  Vector<TransitionTemplate*>::iterator tCur=transitions.begin(),
    tEnd=transitions.end();
  for(; tCur!=tEnd ; ++tCur) delete *tCur;

  /*
  delete gainFactor;
  delete lossFactor;
  delete retentionFactor;
  */
}



void TransducerTemplate::addState(StateTemplate *t)
{
  states.push_back(t);
}



void TransducerTemplate::addTransition(TransitionTemplate *t)
{
  transitions.push_back(t);
}



int TransducerTemplate::getNumStates() const
{
  return states.size();
}



int TransducerTemplate::getNumTransitions() const
{
  return transitions.size();
}



StateTemplate &TransducerTemplate::getIthState(int i)
{
  return *states[i];
}



TransitionTemplate &TransducerTemplate::getIthTransition(int i)
{
  return *transitions[i];
}



void TransducerTemplate::decomposeUnique(StateTemplate *&matchState,
					  StateTemplate *&insertState,
					  StateTemplate *&deleteState,
					  Array2D<Lambda::Closure*> &
					  transMatrix) const
{
  int numStates=states.size(), numTrans=transitions.size();
  for(int i=0 ; i<numStates ; ++i) {
    StateTemplate *state=states[i];
    switch(state->getStateType()) {
    case PHMM_MATCH:  matchState=state; break;
    case PHMM_INSERT: insertState=state; break;
    case PHMM_DELETE: deleteState=state; break;
    }
  }

  transMatrix.resize(NUM_PHMM_STATETYPES,NUM_PHMM_STATETYPES);
  transMatrix.setAllTo(NULL);
  for(int i=0 ; i<numTrans ; ++i) {
    TransitionTemplate *trans=transitions[i];
    PHMM_StateType from=trans->getFrom()->getStateType();
    PHMM_StateType to=trans->getTo()->getStateType();
    transMatrix[from][to]=trans->getTransProbFunc();
  }
}



void TransducerTemplate::addTransition(StateTemplate *from,StateTemplate *to,
				       Lambda::Closure *transProb)
{
  addTransition(new TransitionTemplate(from,to,transProb));
}



void TransducerTemplate::setGainLossFactor(Lambda::Closure *c,GainLossType t)
{
  switch(t)
    {
    case GLT_GAIN:      gainFactor=c;      break;
    case GLT_LOSS:      lossFactor=c;      break;
    case GLT_RETENTION: retentionFactor=c; break;
    case GLT_VOID: throw "TransducerTemplate::setGainLossFactor(GLT_VOID)";
    }
}



Lambda::Closure *TransducerTemplate::getGainLossFactor(GainLossType t)
{
  switch(t)
    {
    case GLT_GAIN:      return gainFactor;
    case GLT_LOSS:      return lossFactor;
    case GLT_RETENTION: return retentionFactor;
    case GLT_VOID:      return NULL;
    }
  throw "TransducerTemplate::getGainLossFactor";
}



/****************************************************************
                     MacrostateTemplate methods
 ****************************************************************/
MacrostateTemplate::MacrostateTemplate(int stateID,TransducerTemplate
				       *submodel)
  : stateID(stateID), submodel(submodel)
{
  // ctor
}



int MacrostateTemplate::getStateID() const
{
  return stateID;
}



TransducerTemplate *MacrostateTemplate::getSubmodel()
{
  return submodel;
}



void MacrostateTemplate::setSubmodel(TransducerTemplate *t)
{
  submodel=t;
}



/****************************************************************
                     MacrotransTemplate methods
 ****************************************************************/
MacrotransTemplate::MacrotransTemplate(MacrostateTemplate *from,
				       MacrostateTemplate *to,
				       Lambda::Closure *closure)
  : from(from), to(to), transProbFunction(closure)
{
  // ctor
}



MacrostateTemplate *MacrotransTemplate::getFrom() const
{
  return from;
}



MacrostateTemplate *MacrotransTemplate::getTo() const
{
  return to;
}



Lambda::Closure *MacrotransTemplate::getTransProbFunc() const
{
  return transProbFunction;
}



