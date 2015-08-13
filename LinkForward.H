/****************************************************************
 LinkForward.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_LinkForward_H
#define INCL_LinkForward_H
#include <iostream>
#include "LinkFelsenstein.H"
#include "BOOM/GSL/Matrix.H"
#include "BOOM/GSL/Vector.H"
#include "BOOM/Symbol.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                        class LinkForward

  This class implements the "forward" algorithm for PhyloPairHMM's, 
  based on the "link" representation of an alignment (i.e., a 
  phylogeny with explicit indel histories).  A PhyloPairHMM is a 
  PairHMM with Felsenstein Emissions, so the emission into each of 
  the two channels of the PairHMM involves emitting into the 
  phylogeny network connected the root of that channel.  Here we 
  deal with a parent taxon and a child taxon.  The channel for the
  child taxon comprises the clade rooted at that child. For the 
  parent, it's the rest of the phylogeny accessible from the parent 
  taxon, not including any taxa in the child's clade.
 ****************************************************************/
class LinkForward
{
public:
  LinkForward(const BranchHMM &,Taxon &parent,Taxon &child,
	      const AlphabetMap &,int numTaxa,LinkFelsenstein &);
  double getLikelihood() const; // P(S1,S2) summed over all paths
  double operator()(int i,int j,STATE k) const;
  int getFirstDim() const; // largest possible index i
  int getSecondDim() const; // largest possible index j
  int getThirdDim() const; // largest possible index k
  double getCachedEmitP(STATE k,int col1,int col2);
private:
  const Alphabet &alphabet;
  const AlphabetMap &alphabetMap;
  BOOM::Symbol gapSymbol;
  const BranchHMM &hmm;
  Vector<STATE> QI, QD, QM; // insert/delete/match states, respectively
  int numStates;
  Array3D<double> F; // the dynamic programming matrix
  Array3D<double> cachedEmitP;
  int numTaxa;
  Taxon &parent, &child;
  LinkFelsenstein &felsenstein;
  double likelihood;
  void fillMatrix();
  double getEmitP(STATE,int col1,int col2);
};



#endif

