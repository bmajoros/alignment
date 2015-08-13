/****************************************************************
 PosteriorBackward.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PosteriorBackward_H
#define INCL_PosteriorBackward_H
#include <iostream>
#include "BOOM/GSL/Matrix.H"
#include "BOOM/GSL/Vector.H"
#include "CollapsedOrthologyMatrix.H"
#include "BranchHMM.H"
#include "SparseMatrix3D.H"
using namespace std;
using namespace BOOM;

// Performs the backward algorithm for ancestral taxa, using a mixture of
// pairwise leaf posteriors instead of emission & transition probabilities

class PosteriorBackward {
public:
  PosteriorBackward(const BranchHMM &hmm,Taxon &parent,Taxon &child,
		    int numTaxa,Vector<int> &insideLeaves,
		    Vector<int> &outsideLeaves,
		    CollapsedOrthologyMatrix &childCOM,
		    CollapsedOrthologyMatrix &parentCOM,
		    Vector<bool> &parentReachability,
		    Vector<bool> &childReachability,
		    Array2D<SparseMatrix3D*> &matrices,
		    double logPseudocount,Array1D<Taxon> &,
		    Phylogeny *,bool wantTransitivity,float indelCoef);
  double getLikelihood() const; // P(S1,S2) summed over all paths
  double operator()(int i,int j,STATE k) const;
  int getFirstDim() const; // largest possible index i
  int getSecondDim() const; // largest possible index j
  int getThirdDim() const; // largest possible index k
  double getEmitP(STATE,int col1,int col2);
  double getEmitP2(STATE,int col1,int col2);
  float getCachedEmitP(STATE k,int col1,int col2);
private:
  const BranchHMM &hmm;
  Vector<STATE> QI, QD, QM; // insert/delete/match states, respectively
  int numStates, nI, nM, nD;
  Array3D<double> B; // the dynamic programming matrix
  Array3D<float> cachedEmitP;
  int numTaxa;
  Array1D<Taxon> &taxa;
  Taxon &parent, &child;
  BranchAttributes *branch;
  IndexMap *upMap, *downMap;
  int numInsideLeaves, numOutsideLeaves;
  Vector<int> &insideLeaves, &outsideLeaves;
  CollapsedOrthologyMatrix &childCOM, &parentCOM;
  Vector<bool> &parentReachability, &childReachability;
  Array2D<SparseMatrix3D*> &matrices;
  double logPseudocount;
  float indelCoef;
  Phylogeny *tree;
  bool wantTransitivity;
  void fillMatrix();
  void resetSavedIterators();
  //double getEmitP2(STATE,int col1,int col2);
};


#endif

