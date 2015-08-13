/****************************************************************
 Viterbi.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Viterbi_H
#define INCL_Viterbi_H
#include <iostream>
#include "PairHMM.H"
#include "StatePath.H"
#include "BOOM/Sequence.H"
#include "BOOM/Array3D.H"
using namespace std;
using namespace BOOM;


class Viterbi
{
public:
  Viterbi(const PairHMM &,const Sequence &,const Sequence &);
  StatePath *getPath();
  double operator()(int i,int j,STATE k) const;
  STATE getTracebackEntry(int i,int j,STATE k) const;
  int getFirstDim() const; // seq1 length
  int getSecondDim() const;// seq2 length
  int getThirdDim() const; // states
private:
  Array3D<double> V;  // dynamic-programming matrix
  Array3D<int> T;     // traceback matrix; entries are states
  const PairHMM &hmm;
  const Sequence &S1, &S2;
  const int L1, L2;       // lengths of S1 and S2
  const int numStates;
  const Symbol gap;
  void fillMatrix();
};

#endif

