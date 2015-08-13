/****************************************************************
 CollapsedOrthologyMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_CollapsedOrthologyMatrix_H
#define INCL_CollapsedOrthologyMatrix_H
#include <iostream>
#include "ResidueOrthologyGraph.H"
#include "Taxon.H"
using namespace std;
using namespace BOOM;

/*
  NOTE: All COM indices (both domain and range) are 1-based!
 */
class CollapsedOrthologyMatrix {
public:
  CollapsedOrthologyMatrix(Taxon &,int numTaxa);
  static void compute(ResidueOrthologyGraph &);
  static const int UNDEFINED; // -1
  struct Entry {
    Entry() : begin(UNDEFINED), end(UNDEFINED) {}
    void set(int b,int e) {begin=b;end=e;}
    void set(int x) {begin=end=x;}
    void printOn(ostream &);
    bool isDirectMatch() {return begin==end;}
    Entry &operator+=(const Entry &);
    int begin;
    int end; // inclusive for match states, not inclusive for indel states
  };
  Entry &operator()(int pos,int taxonID);
  int getNumTaxa();
  int getSeqLen();
  inline Taxon &getTaxon() {return taxon;}
protected:
  Taxon &taxon;
  Array2D<Entry> M; // position X taxonID
  static void upPassInternal(Taxon &,ResidueOrthologyGraph &);
  static void downPass(Taxon &,ResidueOrthologyGraph &);
};

ostream &operator<<(ostream &,const CollapsedOrthologyMatrix::Entry &);

#endif

