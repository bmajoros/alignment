/****************************************************************
 BandedForward.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BandedForward_H
#define INCL_BandedForward_H
#include <iostream>
#include "BandedFB_Base.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                        class BandedForward

  This class implements the "forward" algorithm for PairHMM's on
  leaf sequences.  The matrix is banded and constrained by prescans.
 ****************************************************************/
class BandedForward : public BandedFB_Base
{
public:
  BandedForward(const BranchHMM &,Sequence &left,Sequence &right,
		const AlphabetMap &,BandingPattern &);
  virtual double getLikelihood() const; // P(S1,S2) summed over all paths
private:
  Array3D<double> &F;
  double likelihood;
  virtual void fillMatrix();
};



#endif

