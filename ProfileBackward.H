/****************************************************************
 ProfileBackward.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ProfileBackward_H
#define INCL_ProfileBackward_H
#include <iostream>
#include "PairHMM/PairHMM.H"
#include "BOOM/Array3D.H"
#include "BOOM/Array1D.H"
#include "BOOM/Sequence.H"
#include "BOOM/Map.H"
#include "PhyLib/Phylogeny.H"
#include "AlignmentView.H"
using namespace std;
using namespace BOOM;


/*
  This implements the "backward" algorithm for PairHMM's with Felsenstein
  Emissions, meaning that each of the two channels of the PairHMM involves
  a subgraph of the phylogeny, and the probability of emission into a channel
  is the likelihood of observing the character states of the associated
  phylogeny subgraph.  

  #######
  WARNING: This code hasn't been updated to allow multiple insert/delete/match
           states in an HMM.
  #######

 */
class ProfileBackward
{
public:
  ProfileBackward(const PairHMM &,const AlignmentView &insideClade,
		  const AlignmentView &outsideClade,Phylogeny &,
		  const Array1D<int> &trackMap,const AlphabetMap &);
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
  const PairHMM &hmm;
  const AlignmentView &inside, &outside;
  const MultSeqAlignment &fullAlignment; // doesn't include unobserved taxa
  Phylogeny &phylogeny;
  MultSeqAlignment masterCol; // one column; both observed & unobserved taxa
  Array3D<double> B, cachedEmitP;
  PhylogenyNode &root;
  Map<String,int> nameToTrackId;
  const Array1D<int> &trackMap; // maps track ID's, fullAlignmt -> masterCol
  void initMasterCol();
  void initTrackMap();
  void resetMasterCol();
  void fillMatrix();
  double getEmitP(STATE,int col1,int col2);
  void installColumn(int col,const AlignmentView &);
};

#endif

