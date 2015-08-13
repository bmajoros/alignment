/****************************************************************
 Transducer.C
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "Transducer.H"
using namespace std;
using namespace BOOM;

Transducer::Transducer(const Alphabet &alphabet,Symbol gapSymbol,
		       int numStates)
  : PairHMM(alphabet,gapSymbol,numStates)
{
  // ctor
}



Transducer::Transducer(const Alphabet &,Symbol gapSymbol,
		       const String &filename)
  : PairHMM(alphabet,gapSymbol,filename)
{
  // ctor

  normalize();
}



Transducer::Transducer(const PairHMM &h)
  : PairHMM(h)
{
  // ctor

  normalize();
}



void Transducer::normalize() // normalize emissions to be conditional
{
  if(isInLogSpace()) 
    throw "Transducer::normalize(): take out of log space first!";
  int nAlpha=alphabet.size();
  for(STATE state=1 ; state<numStates ; ++state) // skip state 0
    for(Symbol a=0 ; a<nAlpha ; ++a) {
      double sum=0;
      for(Symbol b=0 ; b<nAlpha ; ++b)
	sum+=emitP(state,a,b);
      for(Symbol b=0 ; b<nAlpha ; ++b)
	emitP(state,a,b)/=sum;
    }
}







