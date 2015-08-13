/****************************************************************
 ProfileViterbi.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ProfileViterbi_H
#define INCL_ProfileViterbi_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Array3D.H"
#include "PhyLib/Phylogeny.H"
#include "PairHMM/PairHMM.H"
#include "PairHMM/StatePath.H"
#include "BranchHMM.H"
#include "Taxon.H"
#include "ProfileFelsenstein.H"
#include "Banding.H"
using namespace std;
using namespace BOOM;


class ProfileViterbi
{
public:
  ProfileViterbi(const BranchHMM &,const MultSeqAlignment &A1,
		 PhylogenyNode *child1,const MultSeqAlignment &A2,
		 PhylogenyNode *child2,Phylogeny &,const AlphabetMap &,
		 BandingType,int bandwidth);
  StatePath *getPath(double &bestScore);
  MultSeqAlignment *decodeAlignment(const StatePath &);
  double operator()(int i,int j,STATE k) const;
  STATE getTracebackEntry(int i,int j,STATE k) const;
  int getFirstDim() const; // seq1 length
  int getSecondDim() const;// seq2 length
  int getThirdDim() const; // states
private:
  Array3D<double> V;  // dynamic-programming matrix
  Array3D<int> T;     // traceback matrix; entries are states
  const BranchHMM &hmm;
  const MultSeqAlignment &S1, &S2; // the profiles being aligned
  MultSeqAlignment masterCol; // combines one column from both profiles
  PhylogenyNode *child1, *child2; 
  Taxon *parentTaxon;
  Phylogeny &phylogeny;
  Array1D<int> trackMap1, trackMap2; // maps track ID's to masterCol
  Map<String,int> nameToTrackId;
  const int L1, L2;       // lengths of S1 and S2
  const int numStates;
  const Symbol gap;
  const Alphabet &alphabet;
  const AlphabetMap &alphabetMap;
  Array1D<double> eqFreqs;
  ProfileFelsenstein *F;
  BandingType bandingType;
  int bandwidth;
  void fillMatrix();
  void setupMasterCol();
  void setupMapping(const MultSeqAlignment &,Array1D<int> &trackMap);
  void copyColumn(const MultSeqAlignment &from,int fromCol,
		  MultSeqAlignment &to,int toCol,int numTracks,
		  const Array1D<int> &trackMap,bool useGaps=false);
  void installColumn(int col,const MultSeqAlignment &,Array1D<int> &trakMap);
  double getEmitP(STATE k,int col1,int col2);
};

#endif

