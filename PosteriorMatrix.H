/****************************************************************
 PosteriorMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PosteriorMatrix_H
#define INCL_PosteriorMatrix_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Array3D.H"
#include "BandedForward.H"
#include "BandedBackward.H"
using namespace std;
using namespace BOOM;


class PosteriorMatrix {
public:
  PosteriorMatrix() {}
  PosteriorMatrix(BandedForward &,BandedBackward &);
  double operator()(int i,int j,int k);
  bool save(double sparsenessThreshold,const String &filename);
  bool saveBinary(double sparsenessThreshold,const String &filename);
  static PosteriorMatrix *loadBinary(const String &filename);
protected:
  Array3D<double> M;
  void compute(BandedForward &,BandedBackward &);
  virtual void loadFromBinary(const String &filename);
};

#endif

