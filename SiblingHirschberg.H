/****************************************************************
 SiblingHirschberg.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_SiblingHirschberg_H
#define INCL_SiblingHirschberg_H
#include <iostream>
#include "Hirschberg.H"
using namespace std;
using namespace BOOM;

class BranchAttributes;

class SiblingHirschberg : public Hirschberg {
public:
  SiblingHirschberg(BranchHMM *,Taxon *left,Taxon *right,
		    LinkFelsenstein &F,int bandwidth,bool usePrescan,
		    int maxThreads,ContentSensor *);
  virtual ~SiblingHirschberg() {}
  virtual Taxon *getFirstTaxon() { return left; }
  virtual Taxon *getSecondTaxon() { return right; }
protected:
  Taxon *left, *right;
  BranchAttributes *leftBranch, *rightBranch;
  IndexMap *leftDownMap, *leftUpMap, *rightDownMap, *rightUpMap;
  int *leftDownMapEntry, *rightDownMapEntry;
  virtual double getEmitP(STATE,int parentResidueIndex,
			  int childResidueIndex,PairHMMStateType);
  virtual double emission(STATE,int parentResidueIndex,
			  int childResidueIndex,PairHMMStateType);
  virtual void precomputeEmissions();
  void precomputeEmissions(PrecomputedEmissions &,Taxon &cladeRoot);
};

#endif

