/****************************************************************
 HirschForwardSum.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_HirschForwardSum_H
#define INCL_HirschForwardSum_H
#include <iostream>
#include "HirschPass.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 The forward recurrence for the Hirschberg algorithm.
 ****************************************************************/
class HirschForwardSum : public HirschPass {
public:
  HirschForwardSum(int minX,int maxX,int minY,int maxY,
		const BandingPattern &,BranchHMM &,
		const Vector<STATE> &QI,int numI,
		const Vector<STATE> &QD,int numD,
		const Vector<STATE> &QM,int numM,
		ViterbiConstraint,PrecomputedEmissions *,
		FunctionalParse *parentFP, FunctionalParse *childFP,
		IndexMap *upMap,IndexMap *downMap,Hirschberg *,
		   ContentSensor *);
  virtual void run();
  virtual void runPrescan();
  float scorePath(StatePath &);
  double computeLikelihood();
  void debug();
};


#endif

