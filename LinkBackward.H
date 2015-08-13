/****************************************************************
 LinkBackward.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_LinkBackward_H
#define INCL_LinkBackward_H
#include <iostream>
#include "LinkFelsenstein.H"
#include "BOOM/GSL/Matrix.H"
#include "BOOM/GSL/Vector.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                        class LinkBackward

  This class implements the "backward" algorithm for PhyloPairHMM's, 
  based on the "link" representation of an alignment (i.e., a 
  phylogeny with explicit indel histories).  A PhyloPairHMM is a 
  PairHMM with Felsenstein Emissions, so the emission into each of 
  the to channels of the PairHMM involves emitting into the 
  phylogeny network connected the root of that channel.  Here we 
  deal with a parent taxon and a child taxon.  The channel for the
  child taxon comprises the clade rooted at that child. For the 
  parent, it's the rest of the phylogeny accessable from the parent 
  taxon, not including any taxa in the child's clade.
 ****************************************************************/
class LinkBackward
{
public:
  LinkBackward(const BranchHMM &,Taxon &parent,Taxon &child,
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
  Symbol gapSymbol;
  const BranchHMM &hmm;
  Vector<STATE> QI, QD, QM; // insert/delete/match states, respectively
  int numStates, nI, nM, nD;
  Array3D<double> B; // the dynamic programming matrix
  Array3D<double> cachedEmitP;
  int numTaxa;
  Taxon &parent, &child;
  BranchAttributes *branch;
  IndexMap *upMap, *downMap;
  LinkFelsenstein &F;
  void fillMatrix();
  double getEmitP(STATE,int col1,int col2);
};




/****************************************************************
                       class RootBackward

 A backward algorithm for sampling when the only thing known is
 the length of the parent sequence.  Useful for sampling a root
 sequence of a phylogeny when you know what length you want the
 sequence to be.  Emission probabilities are not used.
 ****************************************************************/
class RootBackward {
public:
  RootBackward(const BranchHMM &,Taxon &parent,Taxon &child,
	       const AlphabetMap &,int numTaxa,int parentLength);
  double getLikelihood() const; // P(parentLen) summed over all paths
  double operator()(int i,STATE k) const;
  int getFirstDim() const; // largest possible index i (parent pos)
  int getSecondDim() const; // largest possible index k (state)
  GSL::Vector &getTstarV();
private:
  int parentLength;
  GSL::Matrix T; // transition matrix for insertion states only
  GSL::Matrix Tstar; // sum of geometric series for T
  GSL::Vector V; // trans probs from insertion states to final state
  GSL::Vector TstarV; // Tstar * V
  const Alphabet &alphabet;
  const AlphabetMap &alphabetMap;
  Symbol gapSymbol;
  const BranchHMM &hmm;
  int numStates, nI, nM, nD, nMD;
  Array1D<int> h_to_ih;
  Array2D<double> B; // the dynamic programming matrix
  int numTaxa;
  Taxon &parent, &child;
  Vector<STATE> QI, QD, QM; // insert/delete/match states, respectively
  Vector<STATE> QMD; // match and delete states
  Array1D<double> sumW;
  void computeTstar();
  void computeW(GSL::Vector &W,GSL::Vector &TstarW,int i);
  void fillMatrix();
};



/****************************************************************
                       class ChildBackward

 A backward algorithm for sampling when the only thing known is
 the parent sequence (and the parent sequence's functional parse).  
 ****************************************************************/
class ChildBackward {
public:
  ChildBackward(const BranchHMM &,Taxon &parent,Taxon &child,
	       const AlphabetMap &,int numTaxa);
  double getLikelihood() const; // P(parentSeq) summed over all paths
  double operator()(int i,STATE k) const;
  int getFirstDim() const; // largest possible index i (parent pos)
  int getSecondDim() const; // largest possible index k (state)
  GSL::Vector &getTstarV();
  double getCachedEmitP(STATE k,int parentPosition);
private:
  GSL::Matrix T; // transition matrix for insertion states only
  GSL::Matrix Tstar; // sum of geometric series for T
  GSL::Vector V; // trans probs from insertion states to final state
  GSL::Vector TstarV; // Tstar * V
  const Alphabet &alphabet;
  const AlphabetMap &alphabetMap;
  Symbol gapSymbol;
  const BranchHMM &hmm;
  int numStates, nI, nM, nD, nMD;
  Array1D<int> h_to_ih;
  Array2D<double> B; // the dynamic programming matrix
  Array2D<double> cachedEmitP; // state x parent position
  int numTaxa;
  Taxon &parent, &child;
  FunctionalParse &parentFuncParse;
  Vector<STATE> QI, QD, QM; // insert/delete/match states, respectively
  Vector<STATE> QMD; // match and delete states
  Array1D<double> sumW;
  void computeTstar();
  void computeW(GSL::Vector &W,GSL::Vector &TstarW,int i);
  void fillMatrix();
  double getEmitP(STATE,int parentPosition);
};

#endif

