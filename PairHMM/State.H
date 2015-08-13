/****************************************************************
 State.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_State_H
#define INCL_State_H
#include <iostream>
#include "BOOM/String.H"
using namespace std;
using namespace BOOM;

enum PairHMMStateType {
  PHMM_MATCH,
  PHMM_INSERT,
  PHMM_DELETE,
   // THE ABOVE STATES MUST REMAIN NUMBERED 0 THROUGH 2!
  PHMM_START_STOP,
  PHMM_SILENT,
  PHMM_OTHER,
  NUM_PHMM_STATETYPES=6
};

inline operator++(PairHMMStateType &t) {t=PairHMMStateType(int(t)+1);}
bool emitsParent(PairHMMStateType);
bool emitsChild(PairHMMStateType);
void updateCoords(PairHMMStateType,int &x,int &y);
void updateCoordsRev(PairHMMStateType,int &x,int &y);
PairHMMStateType swapIndels(PairHMMStateType);

ostream &operator<<(ostream &,PairHMMStateType);
PairHMMStateType getStateType(const String &stateType);
String toString(PairHMMStateType);

typedef PairHMMStateType PHMM_StateType;

typedef short int STATE; // # changed to short on 9/9/09

static const STATE INVALID_STATE=(short)-1; 

#endif

