/****************************************************************
 SiblingViterbi.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_SiblingViterbi_H
#define INCL_SiblingViterbi_H
#include <iostream>
#include "Taxon.H"
#include "PairHMM/StatePath.H"
#include "BranchAttributes.H"
#include "BranchHMM.H"
#include "Banding.H"
#include "LinkFelsenstein.H"
#include "LinkViterbi.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 A version of LinkViterbi which aligns two sibling taxa.
 ****************************************************************/

class SiblingViterbi : public LinkViterbi {
public:
  SiblingViterbi(BranchHMM *hmm,Taxon *left,Taxon *right,
		 LinkFelsenstein &F,BandingType=NO_BANDING,
		 int bandwidth=100);
  virtual ~SiblingViterbi() {}
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

