/****************************************************************
 BackwardAlgorithm.H
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BackwardAlgorithm_H
#define INCL_BackwardAlgorithm_H
#include <iostream>
#include "PairHMM.H"
#include "BOOM/Array3D.H"
#include "BOOM/Sequence.H"
using namespace std;
using namespace BOOM;


class BackwardAlgorithm
{
public:
  BackwardAlgorithm(const PairHMM &,const Sequence &,const Sequence &);
  double getLikelihood() const; // P(S1,S2) summed over all paths
  double operator()(int i,int j,STATE k) const;
  int getFirstDim() const; // largest possible index i
  int getSecondDim() const; // largest possible index j
  int getThirdDim() const; // largest possible index k
private:
  const PairHMM &hmm;
  const Sequence &S1, &S2;
  Array3D<double> B;
  void fillMatrix();
};

#endif

