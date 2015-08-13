/****************************************************************
 GapPattern.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "GapPattern.H"
using namespace std;
using namespace BOOM;


GapPattern::GapPattern(const Sequence &seq)
{
  const Alphabet &alpha=GapPatternAlphabet::global();
  Symbol unknown=alpha.lookup('?');
  Symbol gap=alpha.lookup('-');
  Symbol residue=alpha.lookup('*');
  int L=seq.getLength();
  pattern.resize(L);
  for(int i=0 ; i<L ; ++i) {
    Symbol s=seq[i];
    if(s==unknown) pattern[i]=GPE_UNKNOWN;
    else if(s==gap) pattern[i]=GPE_GAP;
    else if(s==residue) pattern[i]=GPE_RESIDUE;
    else throw "Error in GapPattern::GapPattern()";
  }
}



GapPatternElement &GapPattern::operator[](int i)
{
  return pattern[i];
}



int GapPattern::getLength() const
{
  return pattern.size();
}



void GapPattern::printOn(ostream &os) const
{
  int L=getLength();
  for(int i=0 ; i<L ; ++i)
    switch(pattern[i]) {
    case GPE_UNKNOWN: os<<'?'; break;
    case GPE_GAP: os<<'-'; break;
    case GPE_RESIDUE: os<<'*'; break;
    }
}


ostream &operator<<(ostream &os,const GapPattern &g)
{
  g.printOn(os);
  return os;
}



int GapPattern::countResidues() const
{
  int L=getLength(), n=0;
  for(int i=0 ; i<L ; ++i)
    switch(pattern[i]) {
    case GPE_UNKNOWN: throw "GPE_UNKNOWN in GapPattern::countResidues()";
    case GPE_RESIDUE: ++n; break;
    }
  return n;
}


