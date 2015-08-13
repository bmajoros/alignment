/****************************************************************
 BandingPattern.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BandingPattern_H
#define INCL_BandingPattern_H
#include <iostream>
#include "BOOM/Array1D.H"
using namespace std;
using namespace BOOM;

struct Bounds {
  int minY, maxY;
};

class BandingPattern {
public:
  BandingPattern(int yDim,int xDim,int bandWidth);// yDim = length of first seq
  void getBounds(int x,int &minY,int &maxY) const;
  int getXdim() const;
  int getYdim() const;
  int getBandwidth() const;
  int getMaxColumnHeight() const;
protected:
  Array1D<Bounds> bounds;
  int xDim, yDim, bandwidth, maxColHeight;
  void computeBounds();
};

#endif

