/****************************************************************
 LinkSampler.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_LinkSampler_H
#define INCL_LinkSampler_H
#include <iostream>
#include "PairHMM/StatePath.H"
#include "BranchHMM.H"
#include "LinkBackward.H"
#include "BranchAttributes.H"
#include "Taxon.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                          class LinkSampler

 Provides sampling of state paths conditional on fixed sequences.
 ****************************************************************/
class LinkSampler
{
public:
  LinkSampler(const BranchHMM &,const LinkBackward &,Taxon &parent,
	      Taxon &child,BranchAttributes *);
  StatePath *samplePath(double &pathScore,double &oldPathScore);
private:
  const BranchHMM &hmm;
  const LinkBackward &B;
  Taxon &parent, &child;
  BranchAttributes *branch;
  double computePathScore(StatePath &);
  void updateCols(STATE,int &i,int &j);
};




/****************************************************************
                        class RootSampler

 For sampling a state path conditional ONLY on the parent sequence
 length (but not the sequence's residues, nor on anything about the
 child, except the branch length to the child).
 ****************************************************************/
class RootSampler {
public:
  RootSampler(const BranchHMM &,const RootBackward &,Taxon &parent,
	      Taxon &child,BranchAttributes *);
  StatePath *samplePath(int parentLength,double &pathScore);
private:
  const BranchHMM &hmm;
  const RootBackward &B;
  Symbol gapSymbol;
  int numStates, nI, nM, nD;
  Vector<STATE> QI, QM, QD;
  Taxon &parent, &child;
  BranchAttributes *branch;
  double computePathScore(StatePath &);
  void updateCol(STATE,int &i);
};




/****************************************************************
                        class ChildSampler

 For sampling a state path conditional on a parent sequence (but
 not conditional on anything about the child, except the branch
 length to the child).
 ****************************************************************/
class ChildSampler {
public:
  ChildSampler(const BranchHMM &,const ChildBackward &,Taxon &parent,
	      Taxon &child,BranchAttributes *);
  StatePath *samplePath(double &pathScore);
private:
  const BranchHMM &hmm;
  const ChildBackward &B;
  Symbol gapSymbol;
  int numStates, nI, nM, nD;
  Vector<STATE> QI, QM, QD;
  Taxon &parent, &child;
  BranchAttributes *branch;
  double computePathScore(StatePath &);
  void updateCol(STATE,int &i);
};



#endif

