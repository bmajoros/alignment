/****************************************************************
 ConstrainedBackward.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ConstrainedBackward_H
#define INCL_ConstrainedBackward_H
#include <iostream>
#include "LinkFelsenstein.H"
#include "BOOM/GSL/Matrix.H"
#include "BOOM/GSL/Vector.H"
#include "BOOM/Symbol.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                  class ConstrainedBackward

  Implements backward algorithm for link representation, given
  a fixed alignment.  Useful for sampling annotations from an HMM
  with multiple match/indel states.
 ****************************************************************/
class ConstrainedBackward {
public:
  ConstrainedBackward(const BranchHMM &,Taxon &parent,Taxon &child,
		      const AlphabetMap &,int numTaxa,LinkFelsenstein &,
		      const FunctionalParse &parentParse);
  double getLikelihood() const; // P(S1,S2) summed over all paths
  double operator()(int i,STATE k) const;
  int getFirstDim() const; // largest possible index i
  int getSecondDim() const; // largest possible index k (states)
  double getCachedEmitP(STATE k,int alignmentColumn);
  Vector<PHMM_StateType> *getPath();
  bool ok();
  void printOn(ostream &) const;
private:
  bool aborted;
  const Alphabet &alphabet;
  const AlphabetMap &alphabetMap;
  BOOM::Symbol gapSymbol;
  const BranchHMM &hmm;
  const FunctionalParse &parentParse;
  bool constrainParent;
  Vector<STATE> QI, QD, QM; // insert/delete/match states, respectively
  int numStates, nI, nM, nD;
  Array2D<double> B; // the dynamic programming matrix
  Array2D<double> cachedEmitP; // indexed by state and column
  int numTaxa;
  Taxon &parent, &child;
  BranchAttributes *branch;
  Vector<PHMM_StateType> *path;
  int pathLen;
  IndexMap *upMap, *downMap;
  LinkFelsenstein &F;
  void fillMatrix();
  double getEmitP(STATE,int parentIndex,int childIndex,int alignmentColumn);
};

ostream &operator<<(ostream &,const ConstrainedBackward &);

#endif

