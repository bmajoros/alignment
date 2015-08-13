/****************************************************************
 PostLeafBackward.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PostLeafBackward_H
#define INCL_PostLeafBackward_H
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

class PostLeafBackward {
public:
  PostLeafBackward(const BranchHMM &hmm,Taxon &parent,Taxon &child,
		   SparseMatrix3D &matrix,double logPseudocount,
		   float indelCoef);
  double getLikelihood() const; // P(S1,S2) summed over all paths
  double operator()(int i,int j,STATE k) const;
  int getFirstDim() const; // largest possible index i
  int getSecondDim() const; // largest possible index j
  int getThirdDim() const; // largest possible index k
  double getEmitP(STATE,int col1,int col2);
  float getCachedEmitP(STATE k,int col1,int col2);
private:
  const BranchHMM &hmm;
  Vector<STATE> QI, QD, QM; // insert/delete/match states, respectively
  int numStates, nI, nM, nD;
  Array3D<double> B; // the dynamic programming matrix
  Array3D<float> cachedEmitP;
  Taxon &parent, &child;
  int parentL, childL;
  SparseMatrix3D &M;
  double logPseudocount;
  float indelCoef;
  void fillMatrix();
};


#endif

