/****************************************************************
 LinkViterbi.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_LinkViterbi_H
#define INCL_LinkViterbi_H
#include <iostream>
#include "Taxon.H"
#include "PairHMM/StatePath.H"
#include "BranchAttributes.H"
#include "BranchHMM.H"
#include "Banding.H"
#include "LinkFelsenstein.H"
#include "Posteriors.H"
#include "BOOM/BandedMatrix3D.H"
#include "PrecomputedEmissions.H"
#include "PhyLib/Phylogeny.H"
#include "ViterbiInterface.H"
using namespace std;
using namespace BOOM;

/****************************************************************
 A viterbi algorithm that operates on the "links" representation,
 as opposed to the ProfileViterbi class which operates on the 
 "profile" representation.  Whereas the profile viterbi operates
 on explicit alignments by aligning one subalignment to another,
 the LinkViterbi aligns a virtual parent sequence to a virtual
 child sequence, with Felsenstein's algorithm providing emission
 probabilities by traversing the "links" denoting indel histories
 across the tree.
 ****************************************************************/

class LinkViterbi : public ViterbiInterface {
public:
  LinkViterbi(Taxon *parent,Taxon *child,LinkFelsenstein &F,
	      BandingType=NO_BANDING,int bandwidth=100);
  virtual ~LinkViterbi() {}
  virtual StatePath *decode(ViterbiConstraint=VC_UNCONSTRAINED);
  virtual void usePosteriors(Posteriors *);
  virtual float scorePath(StatePath &);
protected:
  PrecomputedEmissions precomputedEmissions[2]; // index as PARENT/CHILD
                                                // or LEFT/RIGHT
  const Alphabet &alphabet;
  int numAlpha;
  Posteriors *posteriors;
  LinkFelsenstein &F;
  ViterbiConstraint vc;
  Taxon *parent, *child;
  FunctionalParse *parentFP, *childFP;
  BranchAttributes *branch;
  IndexMap *upMap, *downMap;
  BranchHMM *hmm;
  BandedMatrix3D<float> V; // DP matrix (ROW x COL x STATE)
  BandedMatrix3D<short> T; // traceback matrix; entries are states
  int L1, L2;       // lengths of S1 and S2
  int numStates;
  BandingType bandingType;
  int bandwidth;
  double bestScore;
  void fillMatrix();
  StatePath *getPath(double &bestScore);
  virtual double getEmitP(STATE,int parentResidueIndex,
			  int childResidueIndex,PairHMMStateType);
  virtual float getTransP(STATE h,STATE k);
  virtual double emission(STATE,int parentResidueIndex,
			  int childResidueIndex,PairHMMStateType);
  void debugTraceback();
  virtual void precomputeEmissions();
  void precomputeEmissions(Taxon *forTaxon,BranchEnd parentChild,
			   InsideOutside conditionalOn,Taxon *avoid=NULL);
};

#endif

