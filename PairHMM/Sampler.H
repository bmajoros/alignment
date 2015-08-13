/****************************************************************
 Sampler.H
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Sampler_H
#define INCL_Sampler_H
#include <iostream>
#include "PairHMM.H"
#include "BackwardAlgorithm.H"
#include "StatePath.H"
using namespace std;
using namespace BOOM;

/****************************************************************
 Sampler -- Samples state paths conditional on pairs of sequences.
 ****************************************************************/
class Sampler {
public:
  Sampler(const PairHMM &,const BackwardAlgorithm &,
	  const Sequence &S1,const Sequence &S2);
  StatePath *samplePath();
protected:
  const PairHMM &hmm;
  const BackwardAlgorithm &B;
  const Sequence &S1, &S2;
};


/****************************************************************
 GenerativeSampler -- Samples state paths without conditioning on 
 anything.
 ****************************************************************/
class GenerativeSampler {
public:
  GenerativeSampler(const PairHMM &);
  StatePath *samplePath();
protected:
  const PairHMM &hmm;
};

#endif

