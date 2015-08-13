/****************************************************************
 GapPattern.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_GapPattern_H
#define INCL_GapPattern_H
#include <iostream>
#include "BOOM/Sequence.H"
#include "BOOM/Array1D.H"
#include "GapPatternAlphabet.H"
using namespace std;
using namespace BOOM;

enum GapPatternElement {
  GPE_UNKNOWN,
  GPE_GAP,
  GPE_RESIDUE
};

class GapPattern {
public:
  GapPattern(const Sequence &); // alphabet must be GapPatternAlphabet
  GapPatternElement &operator[](int);
  int getLength() const;
  void printOn(ostream &) const;
  int countResidues() const;
private:
  Array1D<GapPatternElement> pattern;
};

ostream &operator<<(ostream &,const GapPattern &);

#endif

