/****************************************************************
 ProfileFelsenstein.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ProfileFelsenstein_H
#define INCL_ProfileFelsenstein_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Alphabet.H"
#include "BOOM/Vector.H"
#include "PhyLib/Phylogeny.H"
#include "BranchHMM.H"
using namespace std;
using namespace BOOM;


class ProfileFelsenstein
{
public:
  ProfileFelsenstein(const MultSeqAlignment &,const Alphabet &,
		     const AlphabetMap &);
  virtual double logLikelihood(int column,PhylogenyNode *,STATE);
private:
  const MultSeqAlignment &alignment;
  const Alphabet &alphabet;
  const int numAlpha;
  const AlphabetMap &alphabetMap;
};



#endif

