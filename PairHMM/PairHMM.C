/****************************************************************
 PairHMM.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include <math.h>
#include "PairHMM.H"
#include "BOOM/Map.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/RouletteWheel.H"
#include "BOOM/Exceptions.H"
using namespace std;
using namespace BOOM;


PairHMM::PairHMM(const Alphabet &alphabet,Symbol gapSymbol,int numStates)
  : alphabet(alphabet), numStates(numStates), stateTypes(numStates),
    gapSymbol(gapSymbol), isLog(false), transP(numStates,numStates),
    emitP(numStates,alphabet.size(),alphabet.size())
{
  // ctor

  if(numStates>0) {
    stateTypes.setAllTo(0); // ### debugging
    stateTypes[0]=PHMM_START_STOP;
  }
}



PairHMM::PairHMM(const Alphabet &alphabet,Symbol gapSymbol,
		 const String &filename)
  : alphabet(alphabet), gapSymbol(gapSymbol), isLog(false)
{
  // ctor

  load(filename);
}



bool PairHMM::save(const String &filename)
{
  ofstream os(filename.c_str());
  if(!os.good()) throw String("Error saving PairHMM to file: ")+filename;
  os<<numStates<<endl;
  os<<transP<<endl;
  os<<emitP<<endl;

  for(int i=0 ; i<numStates ; ++i)
    os<<static_cast<int>(stateTypes[i])<<" ";
  os<<endl;
}



bool PairHMM::load(const String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw String("Error reading PairHMM from file: ")+filename;
  is>>numStates;
  transP.resize(numStates,numStates);
  int nAlpha=alphabet.size();
  emitP.resize(numStates,nAlpha,nAlpha);
  stateTypes.resize(numStates);
    stateTypes.setAllTo(0); // ### debugging
  is>>transP>>emitP;
  int x;
  for(int i=0 ; i<numStates ; ++i) {
    is>>x;
    stateTypes[i]=static_cast<PHMM_StateType>(x);
  }

  for(STATE s=0 ; s<numStates ; ++s) {
    for(STATE t=0 ; t<numStates ; ++t) 
      if(transP(s,t)<0) isLog=true;
    for(int x=0 ; x<nAlpha ; ++x)
      for(int y=0 ; y<nAlpha ; ++y)
	if(emitP(s,x,y)<0) isLog=true;
  }
}



static PairHMM *PairHMM::load(const Alphabet &alphabet,Symbol gapSymbol,
			      const String &filename)
{
  return new PairHMM(alphabet,gapSymbol,filename);
}



void PairHMM::convertFromLogs()
{
  int nAlpha=alphabet.size();
  for(STATE s=0 ; s<numStates ; ++s) {
    for(STATE t=0 ; t<numStates ; ++t) 
      transP(s,t)=exp(transP(s,t));
    for(int x=0 ; x<nAlpha ; ++x)
      for(int y=0 ; y<nAlpha ; ++y)
	emitP(s,x,y)=exp(emitP(s,x,y));
  }
  isLog=false;
}



void PairHMM::convertToLogs()
{
  int nAlpha=alphabet.size();
  for(STATE s=0 ; s<numStates ; ++s) {
    for(STATE t=0 ; t<numStates ; ++t) 
      transP(s,t)=log(transP(s,t));
    for(int x=0 ; x<nAlpha ; ++x)
      for(int y=0 ; y<nAlpha ; ++y)
	emitP(s,x,y)=log(emitP(s,x,y));
  }
  isLog=true;
}



void PairHMM::setNumStates(int n)
{
  if(n==numStates) return;
  numStates=n;
  transP.resize(n,n);
  int nAlpha=alphabet.size();
  emitP.resize(n,nAlpha,nAlpha);
  stateTypes.resize(n);
  if(n>0) {
    stateTypes.setAllTo(0); // ### debugging
    stateTypes[0]=PHMM_START_STOP;
  }
}



void PairHMM::setEmitP(STATE fromState,Symbol x,Symbol y,double p)
{
  emitP(fromState,x,y)=p;
}



void PairHMM::setTransP(STATE fromState,STATE toState,double p)
{
  transP(fromState,toState)=p;
}



void PairHMM::setStateType(STATE state,PairHMMStateType stateType)
{
  stateTypes[state]=static_cast<PHMM_StateType>(stateType);
}



PairHMM *PairHMM::eliminateSilentStates(Array1D<STATE> &stateMap) const
{
  // ### We assume there are no loops among the silent states

  // First, compute stateMap and its inverse
  Vector<STATE> vStateMap;     // new -> old
  Map<STATE,STATE> inverseMap; // old -> new
  inverseMap[0]=0;
  vStateMap.push_back(0);
  for(STATE i=1 ; i<numStates ; ++i)
    if(stateTypes[i]!=PHMM_SILENT) {
      vStateMap.push_back(i);
      inverseMap[i]=vStateMap.size()-1;
    }
  stateMap=vStateMap;

  // Create a new HMM
  int newNumStates=stateMap.size();
  PairHMM *hmm=new PairHMM(alphabet,gapSymbol,newNumStates);
  hmm->transP.setAllTo(0.0);
  hmm->convertToLogs();
  int nAlpha=alphabet.size();
  for(STATE newState=0 ; newState<newNumStates ; ++newState) {
    STATE oldState=stateMap[newState];
    hmm->setStateType(newState,getStateType(oldState));
    for(Symbol a=0 ; a<nAlpha ; ++a)
      for(Symbol b=0 ; b<nAlpha ; ++b)
	hmm->setEmitP(newState,a,b,getEmitP(oldState,a,b));
  }

  // Now eliminate silent states using a recursive procedure:
  class Eliminator {
    const Array1D<STATE> &stateMap;
    const Map<STATE,STATE> &inverseMap;
    const PairHMM &oldHMM;
    PairHMM &newHMM;
    int oldN, newN;
    void addNewTrans(STATE oldFrom,STATE oldTo,double P) 
    {
      STATE newFrom=inverseMap[oldFrom], newTo=inverseMap[oldTo];
      double oldP=newHMM.getTransP(newFrom,newTo);
      double newP=sumLogProbs(oldP,P);
      newHMM.setTransP(newFrom,newTo,newP);
    }
    void recurs(STATE newDest,STATE oldDest,STATE oldCur,double P)
    {
      if(oldHMM.getStateType(oldCur)==PHMM_SILENT) {
	for(STATE q=0 ; q<oldN ; ++q) {
	  double p=oldHMM.getTransP(q,oldCur);
	  if(finite(p)) recurs(newDest,oldDest,q,P+p);
	}
      }
      else {
	addNewTrans(oldCur,oldDest,P);
      }
    }	
  public:
    Eliminator(const PairHMM &oldHMM,PairHMM &newHMM,
	       const Array1D<STATE> &stateMap,
	       const Map<STATE,STATE> &inverseMap)
      : oldHMM(oldHMM), newHMM(newHMM), stateMap(stateMap), 
	inverseMap(inverseMap), oldN(oldHMM.getNumStates()),
	newN(newHMM.getNumStates())
    {
      for(STATE newState=0 ; newState<newN ; ++newState) {
	STATE oldState=stateMap[newState];
	for(STATE oldPred=0 ; oldPred<oldN ; ++oldPred) {
	  double p=oldHMM.getTransP(oldPred,oldState);
	  if(finite(p)) recurs(newState,oldState,oldPred,p);
	}
      }
    }
    
  };
  Eliminator eliminator(*this,*hmm,stateMap,inverseMap);
  return hmm;
}



void PairHMM::printOn(ostream &os) const
{
  os<<"gapSymbol="<<gapSymbol<<endl;
  os<<numStates<<" states"<<endl;
  os<<"logspace="<<isLog<<endl;
  os<<"transition probabilities:\n"<<transP<<endl;
  os<<"emission probabilities:\n"<<emitP<<endl;
  os<<"stateTypes:\n"<<stateTypes<<endl;
}



ostream &operator<<(ostream &os,const PairHMM &hmm)
{
  hmm.printOn(os);
  return os;
}



void PairHMM::decode(const StatePath &phi,IndexMap &parentToChild,
		     IndexMap &childToParent,bool resize) const
{
  // Resize the maps, if necessary:
  int L=phi.length();
  if(resize) {
    int parentL=0, childL=0;
    phi.getSeqLengths(parentL,childL);
    parentToChild.resize(parentL);
    childToParent.resize(childL);
  }
  
  // Perform decoding:
  int parentI=0, childI=0;
  for(int i=0 ; i<L ; ++i)
    switch(getStateType(phi[i])) 
      {
      case PHMM_MATCH: 
	parentToChild[parentI]=childI;
	childToParent[childI]=parentI;
	++parentI;
	++childI;
	break;
      case PHMM_INSERT:
	childToParent[childI]=IndexMap::UNDEFINED;
	++childI;
	break;
      case PHMM_DELETE:
	parentToChild[parentI]=IndexMap::UNDEFINED;
	++parentI;
	break;
      }
}



double PairHMM::getTransProbs(const StatePath &path) const
{
  double logP=0.0;
  int n=path.size();
  STATE first=path[0];
  STATE prev=first;
  for(int i=1 ; i<n ; ++i) {
    STATE next=path[i];
    logP+=transP(prev,next);
    prev=next;
  }
  logP+=transP(0,first);
  logP+=transP(prev,0);
  return logP;
}



void PairHMM::initPredSuccLists()
{
  predecessors.resize(numStates);
  successors.resize(numStates);
  if(!isInLogSpace()) 
    throw "PairHMM::initPredSuccLists(): not in log space";
  for(STATE i=0 ; i<numStates ; ++i) {
    Vector<STATE> pred, succ;
    for(STATE j=0 ; j<numStates ; ++j) {
      if(finite(transP(i,j)))
	succ.push_back(j);
      if(finite(transP(j,i)))
	pred.push_back(j);
    }
    predecessors[i]=pred;
    successors[i]=succ;
  }
}



STATE PairHMM::sampleNextState(STATE fromState) const
{
  const Array1D<STATE> &succ=successors[fromState];
  int n=succ.size();
  RouletteWheel wheel;
  for(STATE k=0 ; k<n ; ++k) {
    const STATE toState=succ[k];
    double P=transP(fromState,toState);
    wheel.addSector(exp(P));
  }
  wheel.doneAddingSectors();
  int k=wheel.spin(); // could be made faster
  if(k>=n) k=n-1;
  return succ[k];
}



void PairHMM::normalizeTransProbs()
{
  for(int i=0 ; i<numStates ; ++i) {
    double sum=0.0;
    for(int j=0 ; j<numStates ; ++j) sum+=transP(i,j);
    if(sum<0) {
      cout<<"PairHMM::normalizeTransProbs() applied to log model"<<endl;
      throw 0;
    }
    for(int j=0 ; j<numStates ; ++j) transP(i,j)/=sum;
  }
}


