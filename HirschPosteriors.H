/****************************************************************
 HirschPosteriors.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_HirschPosteriors_H
#define INCL_HirschPosteriors_H
#include <iostream>
#include "Hirschberg.H"
#include "PhyLib/Phylogeny.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                         class HirschPosteriors
  Computes posterior probabilities for binding sites.
  ****************************************************************/
class HirschPosteriors : public Hirschberg {
public:
  HirschPosteriors(Taxon *parent,Taxon *child,LinkFelsenstein &F,
		   int bandwidth,int seq1Len,int seq2Len,bool usePrescan,
		   int maxThreads=8,ContentSensor * =NULL);
  double compute(int siteBegin,int siteEnd,BranchEnd=CHILD);
  double computeNumerator(int siteBegin,int siteEnd,BranchEnd=CHILD);
  double computeLikelihood(); // denominator
protected:
  //double likelihood;
};

#endif

