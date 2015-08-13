/****************************************************************
 HirschPass.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_HirschPass_H
#define INCL_HirschPass_H
#include <iostream>
#include "BOOM/PureDnaAlphabet.H"
#include "BandingPattern.H"
#include "BranchHMM.H"
#include "HirschbergFrame.H"
#include "Hirschberg.H"
#include "Taxon.H"
#include "PairHMM/StatePath.H"
#include "BranchAttributes.H"
#include "BranchHMM.H"
#include "LinkFelsenstein.H"
#include "Posteriors.H"
#include "PrecomputedEmissions.H"
#include "PhyLib/Phylogeny.H"
#include "ViterbiConstraint.H"
using namespace std;
using namespace BOOM;

//class Hirschberg;

/****************************************************************
 Base class for the two Hirschberg algorithm "passes"
 ****************************************************************/
class HirschPass {
public:
  HirschPass(int minX,int maxX,int minY,int maxY,
	     const BandingPattern &bp,BranchHMM &hmm,
	     const Vector<STATE> &QI,int numI,
	     const Vector<STATE> &QD,int numD,
	     const Vector<STATE> &QM,int numM,
	     HirschbergDirection,ViterbiConstraint,
	     PrecomputedEmissions *precomputedEmissions,
	     FunctionalParse *parentFP, FunctionalParse *childFP,
	     IndexMap *upMap,IndexMap *downMap,Hirschberg *,
	     ContentSensor *);
  inline HirschbergFrame &getFrame();
  virtual void run()=0;
  virtual void runPrescan()=0;
  void setFirstTaxon(Taxon *t) {firstTaxon=t;}
  void setSecondTaxon(Taxon *t) {secondTaxon=t;}
  virtual double getEmitP(STATE,int parentResidueIndex,
			  int childResidueIndex,PairHMMStateType);
  void setConstraintInterval(int from,int to); // relative to child seq
  void resetConstraintInterval();
protected:
  PureDnaAlphabet alphabet;
  ContentSensor *contentSensor;
  Hirschberg *hirschberg;
  Taxon *firstTaxon, *secondTaxon;
  FunctionalClass bg;
  HirschbergFrame frame;
  ViterbiConstraint vc;
  pair<int,int> constraintInterval; // relative to child seq
  const BandingPattern &bandingPattern;
  BranchHMM &hmm;
  int numAlpha;
  int minX, maxX, minY, maxY, numStates;
  const Vector<STATE> &QD, &QI, &QM; // delete/insert/match states
  const int numD, numI, numM; // number of delete/insert/match states
  FunctionalParse *parentFP, *childFP;
  IndexMap *upMap, *downMap;
  PrecomputedEmissions *precomputedEmissions; // index as PARENT/CHILD
                                                // or LEFT/RIGHT
  virtual double emission(STATE,int parentResidueIndex,
			  int childResidueIndex,PairHMMStateType);
  inline double getEQ(const Array1D<double> &eqFreqs,Symbol,Taxon *,
		      int pos,FunctionalClass) ;
  inline float getTransP(STATE from,STATE to);
  void checkColumn(int x,int minY,int maxY);
  bool columnIsEmpty(int x,int minY,int maxY);
  void getPrescanStates(int firstIndex,int secondIndex,Vector<STATE> &into);
  void getCrossFuncStates(FuncClassSet &,BranchEnd branchEnd,
			  Vector<STATE> &into);
  void getRetentionStates(FuncClassSet &,Vector<STATE> &into);
};

inline HirschbergFrame &HirschPass::getFrame() {
  return frame;
}



double HirschPass::getEQ(const Array1D<double> &eqFreqs,Symbol s,
			 Taxon *taxon,int pos,FunctionalClass fc) 
{
  if(contentSensor && fc.isBackground()) {
    const Sequence &seq=taxon->getSeq();
    const String &seqStr=taxon->getSeqStr();
    const char c=alphabet.lookup(s);
    return contentSensor->scoreSingleBase(seq,seqStr,pos,s,c);
  }
  else return log(eqFreqs[s]);
}



float HirschPass::getTransP(STATE from,STATE to)
{
  return hirschberg ? hirschberg->getTransP(from,to) : hmm.getTransP(from,to);
}



#endif

