/****************************************************************
 TemplateInstantiator.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_TemplateInstantiator_H
#define INCL_TemplateInstantiator_H
#include <iostream>
#include "BranchHMM.H"
#include "TransducerTemplate.H"
#include "BOOM/Alphabet.H"
#include "BOOM/lambda/LambdaAPI.H"
#include "BOOM/AlphabetMap.H"
using namespace std;
using namespace BOOM;


struct TemplateInstantiator {
  TemplateInstantiator(LambdaAPI &);
  BranchHMM *instantiate(TransducerTemplate *,double branchLength,
			 ostream * =NULL);
private:
  LambdaAPI &lambda;
  AstCombination *caller;
  AstObject *closureWrapper;
  AstFloatLit *parmNode;
  const Alphabet &alphabet;
  BOOM::Symbol gapSymbol;
  DropGapMapping alphabetMap;
  void initCaller();
  void setGainLossFactor(BranchHMM *,TransducerTemplate *,GainLossType,
			 double branchLength);
  void setGainLossFactors(BranchHMM *,TransducerTemplate *,
			  double branchLength);
};

#endif

