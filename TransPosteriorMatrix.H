/****************************************************************
 TransPosteriorMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_TransPosteriorMatrix_H
#define INCL_TransPosteriorMatrix_H
#include <iostream>
#include "PosteriorMatrix.H"
#include "PairHMM/PairHMM.H"
using namespace std;
using namespace BOOM;


class TransPosteriorMatrix : public PosteriorMatrix {
public:
  TransPosteriorMatrix() {}
  TransPosteriorMatrix(BandedForward &,BandedBackward &,PairHMM &,
		       const Sequence &S1,const Sequence &S2);
  static TransPosteriorMatrix *loadBinary(const String &filename);
  const Vector<PHMM_Trans> &getTransitions() const;
protected:
  void compute(BandedForward &,BandedBackward &,PairHMM &);
  Vector<PHMM_Trans> transitions;
  const Sequence &S1, &S2;
};

#endif

