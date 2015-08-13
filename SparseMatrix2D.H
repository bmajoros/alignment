/****************************************************************
 SparseMatrix2D.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_SparseMatrix2D_H
#define INCL_SparseMatrix2D_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array1D.H"
#include "BOOM/List.H"
using namespace std;
using namespace BOOM;


class SparseMatrix2D {
public:
  SparseMatrix2D(short L1,short L2);
  SparseMatrix2D(short L1); // initializes itself as identity matrix
  SparseMatrix2D(Array2D<float> &,float threshold);
  void initFrom(Array2D<float> &,float threshold);
  void becomeIdentity();
  static SparseMatrix2D *loadBinary(const String &filename);
  void saveBinary(const String &filename);
  struct Entry {
    short y;
    float value;
    Entry(short y,float v) : y(y) ,value(v) {}
  };
  typedef List<Entry> EntryList;
  EntryList &operator()(short x);
  void addEntry(short x,short y,float value);
  SparseMatrix2D *transpose();
  short getFirstDim() {return L1;}
  short getSecondDim() {return L2;}
  float find(short x,short y);
protected:
  Array1D<EntryList> M;
  short L1, L2;
};

#endif

