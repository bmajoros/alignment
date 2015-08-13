/****************************************************************
 Transducer.H
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Transducer_H
#define INCL_Transducer_H
#include <iostream>
#include "PairHMM.H"
using namespace std;
using namespace BOOM;

class Transducer : public PairHMM
{
public:
  Transducer(const Alphabet &,Symbol gapSymbol,int numStates=0);
  Transducer(const Alphabet &,Symbol gapSymbol,const String &filename);
  Transducer(const PairHMM &);
  virtual ~Transducer() {}
private:
  void normalize();
};

#endif

