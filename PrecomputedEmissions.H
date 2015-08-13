/****************************************************************
 PrecomputedEmissions.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PrecomputedEmissions_H
#define INCL_PrecomputedEmissions_H
#include <iostream>
#include "BOOM/Array3D.H"
using namespace std;
using namespace BOOM;


class PrecomputedEmissions {
public:
  PrecomputedEmissions(int seqLen=0,int numFuncClasses=0,int alphabetSize=0);
  void resize(int seqLen,int numFuncClasses,int alphabetSize);
  inline Array3D<float>::IndexedOnce<float> operator[](int position);
  void printOn(ostream &os) {os<<A;}
protected:
  Array3D<float> A; // indexed as: A[position][funcClass][symbol]
};

inline ostream &operator<<(ostream &os,const PrecomputedEmissions &E) {E.printOn(os);return os;}

Array3D<float>::IndexedOnce<float> PrecomputedEmissions::operator[](int 
								    position)
{
  return A[position];
}

#endif

