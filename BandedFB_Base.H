/****************************************************************
 BandedFB_Base.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BandedFB_Base_H
#define INCL_BandedFB_Base_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/Array3D.H"
#include "BOOM/AlphabetMap.H"
#include "BOOM/Sequence.H"
#include "BranchHMM.H"
#include "BandingPattern.H"
using namespace std;
using namespace BOOM;

class BandedFB_Base {
public:
  BandedFB_Base(const Alphabet &,const AlphabetMap &,BOOM::Symbol gap,
		const BranchHMM &,BandingPattern &,Sequence &left,
		Sequence &right);
  int getFirstDim() const; // largest possible index i
  int getSecondDim() const; // largest possible index j
  int getThirdDim() const; // largest possible index k
  virtual double getLikelihood() const=0;
  double operator()(int i,int j,STATE k) const;
protected:
  const Alphabet &alphabet;
  const AlphabetMap &alphabetMap;
  BOOM::Symbol gapSymbol;
  const BranchHMM &hmm;
  BandingPattern &bandingPattern;
  Vector<STATE> QI, QD, QM; // insert/delete/match states, respectively
  int numStates, nI, nM, nD;
  Array3D<double> M; // the dynamic programming matrix
  Sequence &left, &right;
  double getEmitP(STATE,int col1,int col2);
  virtual void fillMatrix()=0;
};

#endif

