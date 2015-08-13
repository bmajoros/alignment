/****************************************************************
 State.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "State.H"
using namespace std;


//const STATE INVALID_STATE=-1;

PairHMMStateType getStateType(const String &stateType)
{
  if(stateType=="MATCH") return PHMM_MATCH;
  if(stateType=="INSERT") return PHMM_INSERT;
  if(stateType=="DELETE") return PHMM_DELETE;
  if(stateType=="SILENT") return PHMM_SILENT;
  if(stateType=="START") return PHMM_START_STOP;
  if(stateType=="OTHER") return PHMM_OTHER;
  throw String("Unknown state type: ")+stateType;
}



ostream &operator<<(ostream &os,PairHMMStateType t)
{
  switch(t) 
    {
    case PHMM_MATCH:       os<<"MATCH";      break;
    case PHMM_INSERT:      os<<"INSERT";     break;
    case PHMM_DELETE:      os<<"DELETE";     break;
    case PHMM_SILENT:      os<<"SILENT";     break;
    case PHMM_START_STOP:  os<<"START/STOP"; break;
    case PHMM_OTHER:       os<<"OTHER";      break;
    default:               os<<"***unrecognized state type***"; break;
    }
  return os;
};



PairHMMStateType swapIndels(PairHMMStateType t)
{
  switch(t) 
    {
    case PHMM_INSERT: t=PHMM_DELETE; break;
    case PHMM_DELETE: t=PHMM_INSERT; break;
    }
  return t;
}



String toString(PairHMMStateType t)
{
  switch(t) 
    {
    case PHMM_MATCH:       return "MATCH";
    case PHMM_INSERT:      return "INSERT";
    case PHMM_DELETE:      return "DELETE";
    case PHMM_SILENT:      return "SILENT";
    case PHMM_START_STOP:  return "START/STOP";
    case PHMM_OTHER:       return "OTHER";
    default:               return "***unrecognized state type***";
    }
};



bool emitsParent(PairHMMStateType t)
{
  return t==PHMM_MATCH || t==PHMM_DELETE;
}



bool emitsChild(PairHMMStateType t)
{
  return t==PHMM_MATCH || t==PHMM_INSERT;
}



void updateCoords(PairHMMStateType t,int &x,int &y)
{
  switch(t) 
    {
    case PHMM_MATCH:  ++x; ++y; break;
    case PHMM_INSERT:      ++y; break;
    case PHMM_DELETE: ++x;      break;
    }
}



void updateCoordsRev(PairHMMStateType t,int &x,int &y)
{
  switch(t) 
    {
    case PHMM_MATCH:  --x; --y; break;
    case PHMM_INSERT:      --y; break;
    case PHMM_DELETE: --x;      break;
    }
}




