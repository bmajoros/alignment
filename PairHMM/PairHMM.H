/****************************************************************
 PairHMM.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PairHMM_H
#define INCL_PairHMM_H
#include <iostream>
#include "BOOM/Alphabet.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array3D.H"
#include "BOOM/String.H"
#include "BOOM/IndexMap.H"
#include "BOOM/Exceptions.H"
#include "State.H"
#include "StatePath.H"
using namespace std;
using namespace BOOM;

class PairHMM
{
public:
  PairHMM(const Alphabet &,Symbol gapSymbol,int numStates=0);
  PairHMM(const Alphabet &,Symbol gapSymbol,const String &filename);
  virtual ~PairHMM() {}
  void decode(const StatePath &,IndexMap &parentToChild,
	      IndexMap &childToParent,bool resize=true) const;
  bool save(const String &);
  static PairHMM *load(const Alphabet &,Symbol gapSymbol,
		       const String &filename);
  inline int getNumStates() const;
  inline Symbol getGapSymbol() const;
  inline const Alphabet &getAlphabet() const;
  virtual void convertToLogs();
  void convertFromLogs();
  inline double getEmitP(STATE from,Symbol,Symbol) const;
  inline double getTransP(STATE from,STATE to) const;
  void setNumStates(int);
  void setEmitP(STATE from,Symbol,Symbol,double);
  void setTransP(STATE from,STATE to,double);
  inline PairHMMStateType getStateType(STATE) const;
  inline bool isDeleteState(STATE) const;
  inline bool isInsertState(STATE) const;
  inline bool isMatchState(STATE) const;
  void setStateType(STATE,PairHMMStateType);
  inline void getStatesOfType(PairHMMStateType,Vector<STATE> &) const;
  inline bool isInLogSpace() const;
  PairHMM *eliminateSilentStates(Array1D<STATE> &stateMap) const;
    // stateMap maps the new HMM's states to this HMM's states
  void printOn(ostream &) const;
  double getTransProbs(const StatePath &) const;
  void initPredSuccLists(); // call AFTER converting to log space!
  inline const Array1D<STATE> &getPredecessors(STATE) const;
  inline const Array1D<STATE> &getSuccessors(STATE) const;
  inline void updateColumnsFwd(STATE,int &firstSeqIndex,int &secondSeqIndex) 
    const;
  inline void updateColumnsRev(STATE,int &firstSeqIndex,int &secondSeqIndex) 
    const;
  STATE sampleNextState(STATE fromState) const;
  void normalizeTransProbs();
protected:
  const Alphabet &alphabet;
  Symbol gapSymbol;
  Array2D<double> transP; // fromState x toState
  Array3D<double> emitP;  // fromState x symbol1 x symbol2
  Array1D< Array1D<STATE> > predecessors, successors;
  int numStates;
  Array1D<PHMM_StateType> stateTypes;
  bool isLog;
  bool load(const String &filename);
};

ostream &operator<<(ostream &,const PairHMM &);



Symbol PairHMM::getGapSymbol() const
{
  return gapSymbol;
}



int PairHMM::getNumStates() const
{
  return numStates;
}



double PairHMM::getEmitP(STATE fromState,Symbol x,Symbol y) const
{
  return emitP(fromState,x,y);
}



double PairHMM::getTransP(STATE fromState,STATE toState) const
{
  return transP[fromState][toState];
}



PairHMMStateType PairHMM::getStateType(STATE state) const
{
  return static_cast<PairHMMStateType>(stateTypes[state]);
}



void PairHMM::getStatesOfType(PairHMMStateType T,Vector<STATE> &V) const
{
  PHMM_StateType t=static_cast<PHMM_StateType>(T);
  for(STATE i=0 ; i<numStates ; ++i)
    if(stateTypes[i]==t)
      V.push_back(i);
}



bool PairHMM::isDeleteState(STATE s) const
{
  return getStateType(s)==PHMM_DELETE;
}



bool PairHMM::isInsertState(STATE s) const
{
  return getStateType(s)==PHMM_INSERT;
}



bool PairHMM::isMatchState(STATE s) const
{
  return getStateType(s)==PHMM_MATCH;
}



bool PairHMM::isInLogSpace() const
{
  return isLog;
}



const Alphabet &PairHMM::getAlphabet() const
{
  return alphabet;
}



const Array1D<STATE> &PairHMM::getPredecessors(STATE i) const
{
  return predecessors[i];
}



const Array1D<STATE> &PairHMM::getSuccessors(STATE i) const
{
  return successors[i];
}



void PairHMM::updateColumnsFwd(STATE q,int &x,int &y) const
{
  switch(getStateType(q))
    {
    case PHMM_MATCH:  ++x; ++y; break;
    case PHMM_INSERT:      ++y; break;
    case PHMM_DELETE: ++x;      break;
    default:  INTERNAL_ERROR;
    }
}



void PairHMM::updateColumnsRev(STATE q,int &x,int &y) const
{
  switch(getStateType(q))
    {
    case PHMM_MATCH:  --x; --y; break;
    case PHMM_INSERT:      --y; break;
    case PHMM_DELETE: --x;      break;
    default:  INTERNAL_ERROR;
    }
}



#endif

