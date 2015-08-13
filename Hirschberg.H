/****************************************************************
 Hirschberg.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Hirschberg_H
#define INCL_Hirschberg_H
#include <iostream>
#include "BOOM/List.H"
#include "BOOM/Mutex.H"
#include "BOOM/Semaphore.H"
#include "EGGS/ContentSensor.H"
#include "ViterbiInterface.H"
#include "Taxon.H"
#include "Posteriors.H"
#include "PairHMM/StatePath.H"
#include "Banding.H"
#include "LinkFelsenstein.H"
#include "BandingPattern.H"
#include "PrecomputedEmissions.H"
using namespace std;
using namespace BOOM;
class HirschPass;


/****************************************************************
                        class HirschCell
 ****************************************************************/
struct HirschCell {
  HirschCell(int x,int y,STATE q,float forward=0.0,float backward=0.0);
  int x, y;
  STATE state;
  float forward, backward;
};



/****************************************************************
                    exception HIRSCHBERG_NO_PATH
 ****************************************************************/
struct HIRSCHBERG_NO_PATH {};



/****************************************************************
                         class Hirschberg

 Implements the Hirschberg algorithm for pair HMM's.  Uses classes
 HirschForward and HirschBackward.  This class deals specifically
 with the alignment of a parent taxon to a child taxon.  A subclass
 of this class (SiblingHirschberg) deals with alignment of sibling 
 taxa.
 ****************************************************************/
class Hirschberg : public ViterbiInterface {
public:
  Hirschberg(Taxon *parent,Taxon *child,LinkFelsenstein &F,
	     int bandwidth,int seq1Len,int seq2Len,bool usePrescan,
	     int maxThreads=8,ContentSensor * =NULL);
  virtual ~Hirschberg();
  virtual StatePath *decode(ViterbiConstraint=VC_UNCONSTRAINED);
  virtual void usePosteriors(Posteriors *);
  virtual float scorePath(StatePath &);
  float getTransP(STATE from,STATE to); // accounts for posteriors
  void recurse(List<HirschCell>::iterator fromIter,
	       List<HirschCell>::iterator toIter);
  virtual Taxon *getFirstTaxon() { return parent; }
  virtual Taxon *getSecondTaxon() { return child; }
  PrecomputedEmissions &getPrecomputedEmissions(int which);
protected:
  ContentSensor *contentSensor;
  Semaphore maxThreads; 
  Mutex computedCellsMutex; 
  BandingPattern bandingPattern;
  List<HirschCell> computedCells;
  Taxon *parent, *child;
  BranchAttributes *branch;
  BranchHMM *hmm;
  bool usePrescan;
  Vector<STATE> QD, QI, QM; // delete/insert/match states
  int numD, numI, numM; // number of delete/insert/match states
  int m, n; // lengths of the two sequences
  int numStates;
  const Alphabet &alphabet;
  int numAlpha;
  FunctionalParse *parentFP, *childFP;
  IndexMap *upMap, *downMap;
  PrecomputedEmissions precomputedEmissions[2]; // index as PARENT/CHILD
                                                // or LEFT/RIGHT
  LinkFelsenstein &F;
  Posteriors *posteriors;
  ViterbiConstraint vc;
  void endRecursion(List<HirschCell>::iterator fromIter,
		    List<HirschCell>::iterator toIter,
		    int minY,int maxY,int minX,int maxX);
  virtual void precomputeEmissions();
  void precomputeEmissions(Taxon *forTaxon,BranchEnd parentChild,
			   InsideOutside conditionalOn,Taxon *avoid=NULL);
  void initStateSets();
  void runFB(HirschPass &F,HirschPass &B);
};

#endif

