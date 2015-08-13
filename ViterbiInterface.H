/****************************************************************
 ViterbiInterface.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ViterbiInterface_H
#define INCL_ViterbiInterface_H
#include <iostream>
#include "ViterbiConstraint.H"
#include "Posteriors.H"
#include "PairHMM/StatePath.H"
using namespace std;
using namespace BOOM;


class ViterbiInterface {
public:
  virtual ~ViterbiInterface() {}
  virtual StatePath *decode(ViterbiConstraint=VC_UNCONSTRAINED)=0;
  virtual void usePosteriors(Posteriors *)=0;
  virtual float scorePath(StatePath &)=0;
};

#endif

