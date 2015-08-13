/****************************************************************
 ProfileSampler.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ProfileSampler_H
#define INCL_ProfileSampler_H
#include <iostream>
#include "PairHMM/PairHMM.H"
#include "PairHMM/StatePath.H"
#include "ProfileBackward.H"
#include "AlignmentView.H"
using namespace std;
using namespace BOOM;

class ProfileSampler
{
public:
  ProfileSampler(const PairHMM &,const ProfileBackward &,
		 const AlignmentView &,const AlignmentView &);
  StatePath *samplePath(double &pathScore,double &oldPathScore);
private:
  const PairHMM &hmm;
  const ProfileBackward &B;
  const AlignmentView &S1, &S2;
  double computeOldPathScore();
};

#endif

