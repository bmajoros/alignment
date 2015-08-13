/****************************************************************
 GapPatternAlphabet.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_GapPatternAlphabet_H
#define INCL_GapPatternAlphabet_H
#include <iostream>
#include "BOOM/Alphabet.H"
using namespace std;
using namespace BOOM;

// Alphabet of 3 symbols: {*,-,?}. Asterisk means "any nucleotide",
// dash means "gap", and question mark means "one of the other two,
// but I don't know yet which one".  

class GapPatternAlphabet : public Alphabet
{
public:
  GapPatternAlphabet();
  static GapPatternAlphabet &global();
};


#endif

