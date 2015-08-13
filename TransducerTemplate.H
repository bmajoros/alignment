/****************************************************************
 TransducerTemplate.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_TransducerTemplate_H
#define INCL_TransducerTemplate_H
#include <iostream>
#include "PhyLib/Phylogeny.H"
#include "BOOM/String.H"
#include "PairHMM/PairHMM.H"
#include "PhyLib/RateMatrix.H"
#include "BOOM/lambda/Closure.H"
#include "FunctionalClass.H"
#include "GainLossType.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                         class StateTemplate
 ****************************************************************/
class StateTemplate {
  int stateID;
  PairHMMStateType stateType; // insert, match, delete, etc...
  FunctionalClass parentFC, childFC;
  bool scoreGainLossFels; // whether felsenstein should apply gain/loss probs
public:
  StateTemplate(int stateID,PairHMMStateType,FunctionalClass parentFC=
		FunctionalClass::NO_CLASS,FunctionalClass childFC=
		FunctionalClass::NO_CLASS);
  PairHMMStateType getStateType() const;
  RateMatrix *getMatrix(BranchEnd) const;
  Phylogeny *getPhylogeny(BranchEnd) const;
  int getStateID() const;
  void changeStateID(int newID);
  void setFunctionalClass(FunctionalClass,BranchEnd);
  FunctionalClass getFunctionalClass(BranchEnd);
  bool shouldScoreGainLoss() const;
  void disableGainLossScoring();
  GainLossType crossFunctionalType();
};



/****************************************************************
                      class TransitionTemplate
 ****************************************************************/
class TransitionTemplate {
  StateTemplate *from, *to;
  Lambda::Closure *transProbFunction; // evaluates to trans prob
public:
  TransitionTemplate(StateTemplate *from,StateTemplate *to,
		     Lambda::Closure *);
  ~TransitionTemplate() {/*delete transProbFunction;*/}
  StateTemplate *getFrom() const;
  StateTemplate *getTo() const;
  Lambda::Closure *getTransProbFunc() const;
};



/****************************************************************
                      class TransducerTemplate
 ****************************************************************/
class TransducerTemplate {
  Vector<StateTemplate*> states;
  Vector<TransitionTemplate*> transitions;
  Lambda::Closure *gainFactor, *lossFactor, *retentionFactor;
public:
  TransducerTemplate();
  virtual ~TransducerTemplate();
  void addState(StateTemplate *);
  void addTransition(TransitionTemplate *);
  void addTransition(StateTemplate *from,StateTemplate *to,
		     Lambda::Closure *transProb);
  int getNumStates() const;
  int getNumTransitions() const;
  StateTemplate &getIthState(int i);
  TransitionTemplate &getIthTransition(int i);
  void decomposeUnique(StateTemplate *&matchState,
		       StateTemplate *&insertState,
		       StateTemplate *&deleteState,
		       Array2D<Lambda::Closure*> &transMatrix) const;
  void setGainLossFactor(Lambda::Closure *,GainLossType);
  Lambda::Closure *getGainLossFactor(GainLossType);
};



/****************************************************************
                       class MacrostateTemplate
 ****************************************************************/
class MacrostateTemplate {
  int stateID;
  TransducerTemplate *submodel; // ### memory leak?
public:
  MacrostateTemplate(int stateID,TransducerTemplate *);
  int getStateID() const;
  TransducerTemplate *getSubmodel();
  void setSubmodel(TransducerTemplate *);
};



/****************************************************************
                      class MacrotransTemplate
 ****************************************************************/
class MacrotransTemplate {
  MacrostateTemplate *from, *to;
  Lambda::Closure *transProbFunction; // evaluates to trans prob
public:
  MacrotransTemplate(MacrostateTemplate *from,MacrostateTemplate *to,
		     Lambda::Closure *);
  MacrostateTemplate *getFrom() const;
  MacrostateTemplate *getTo() const;
  Lambda::Closure *getTransProbFunc() const;
};



#endif

