/****************************************************************
 BandedBackward.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BandedBackward_H
#define INCL_BandedBackward_H
#include <iostream>
#include "BandedFB_Base.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                        class BandedBackward

  This class implements the "backward" algorithm for PairHMM's on
  leaf sequences.  The matrix is banded and constrained by prescans.
 ****************************************************************/
class BandedBackward : public BandedFB_Base
{
public:
  BandedBackward(const BranchHMM &,Sequence &left,Sequence &right,
		 const AlphabetMap &,BandingPattern &);
  virtual double getLikelihood() const; // P(S1,S2) summed over all paths
private:
  Array3D<double> &B;
  virtual void fillMatrix();
};



#endif

