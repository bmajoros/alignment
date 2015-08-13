/****************************************************************
 SparseMatrix3D.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_SparseMatrix3D_H
#define INCL_SparseMatrix3D_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/Array3D.H"
#include "BOOM/List.H"
#include "PairHMM/State.H"
using namespace std;
using namespace BOOM;

class SparseMatrix3D {
public:
  SparseMatrix3D(short L1,short L2,short numStates);
  SparseMatrix3D(Array3D<float> &,float threshold);// L1 x L2 x states
  SparseMatrix3D(SparseMatrix3D &,float threshold);
  void initFrom(Array3D<float> &,float threshold); // L1 x L2 x states
  static SparseMatrix3D *loadBinary(const String &filename);
  void saveBinary(const String &filename);
  struct Entry {
    short y;
    float value;
    Entry(short y,float v) : y(y) ,value(v) {}
    Entry(const Entry &other) : y(other.y) ,value(other.value) {}
  };
  typedef List<Entry> EntryList;
  EntryList &operator()(short x,short state);
  float find(short x,short y,short state);
  void addEntry(short x,short y,short state,float value);
  SparseMatrix3D *transpose();
  short getFirstDim() {return L1;}
  short getSecondDim() {return L2;}
  short getThirdDim() {return numStates;}
  List<Entry>::iterator &getIter(int x,STATE);
  List<Entry>::iterator begin(int x,STATE);
  List<Entry>::iterator end(int x,STATE);
protected:
  Array2D<EntryList> M; // indexed as: (x,state)
  Array2D<EntryList::iterator> iter;
  short L1, L2, numStates;
};

#endif

