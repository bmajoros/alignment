/****************************************************************
 HirschBackwardSumPair.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_HirschBackwardSumPair_H
#define INCL_HirschBackwardSumPair_H
#include <iostream>
#include "HirschPass.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 The backward recurrence for the Hirschberg algorithm.
 ****************************************************************/
class HirschBackwardSumPair : public HirschPass {
  Taxon *parent, *child;
  Sequence &parentSeq, &childSeq;
public:
  HirschBackwardSumPair(int minX,int maxX,int minY,int maxY,
			const BandingPattern &,BranchHMM &,
			const Vector<STATE> &QI,int numI,
			const Vector<STATE> &QD,int numD,
			const Vector<STATE> &QM,int numM,
			ViterbiConstraint,
			Taxon *parent,Taxon *child,
			FunctionalParse *parentFP, FunctionalParse *childFP,
			IndexMap *upMap,IndexMap *downMap,Hirschberg *,
			ContentSensor *);
  virtual void run();
  virtual void runPrescan();
  virtual double getEmitP(STATE,int parentResidueIndex,
			  int childResidueIndex,PairHMMStateType);
  virtual double emission(STATE,int parentResidueIndex,
			  int childResidueIndex,PairHMMStateType);
};


#endif

