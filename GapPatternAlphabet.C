/****************************************************************
 GapPatternAlphabet.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "GapPatternAlphabet.H"
using namespace std;
using namespace BOOM;

GapPatternAlphabet::GapPatternAlphabet()
  : Alphabet("?-*") // DO NOT CHANGE THE ORDERING OF LETTERS IN THIS STRING
{
  // ctor
}



GapPatternAlphabet &GapPatternAlphabet::global() {
  static GapPatternAlphabet alpha;
  return alpha;
}
