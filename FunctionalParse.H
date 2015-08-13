/****************************************************************
 FunctionalParse.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_FunctionalParse_H
#define INCL_FunctionalParse_H
#include <iostream>
#include "FunctionalClass.H"
#include "PhyLib/Phylogeny.H"
#include "FunctionalElement.H"
using namespace std;
using namespace BOOM;

class Taxon;

/****************************************************************
  A "functional parse" is a sequence of labels.  It differs from a 
  StatePath in that distinct states may correspond to the same label.
  Thus, a unique functional parse can be inferred from a state path,
  but not vice-versa (though with additional information and constraints
  on the HMM, the reverse can be done---i.e., no functional class may
  have more than one match/insert/delete state, and the gap patterns
  of the two taxa must be available).
 ****************************************************************/
class FunctionalParse
{
  Vector<FunctionalClass> parse;
  Taxon *taxon;
public:
  FunctionalParse();
  FunctionalParse(const StatePath &,BranchEnd);
  int length() const;
  int size() const;
  void resize(int);
  void push_back(FunctionalClass);
  void clear(); // sets to zero length
  FunctionalClass operator[](int) const;
  FunctionalClass &operator[](int);
  void printOn(ostream &);
  Array1D<FunctionalElement> *getElements();
  inline Vector<FunctionalClass>::iterator begin() {return parse.begin();}
  inline Vector<FunctionalClass>::iterator end() {return parse.end();}
  Taxon *getTaxon() {return taxon;}
  void setTaxon(Taxon *t) {taxon=t;}
};

ostream &operator<<(ostream &,const FunctionalParse &);

#endif

