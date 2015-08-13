/****************************************************************
 BranchHMM.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BranchHMM_H
#define INCL_BranchHMM_H
#include <iostream>
#include "BOOM/Array1D.H"
#include "BOOM/Array3D.H"
#include "PairHMM/Transducer.H"
#include "PhyLib/SubstitutionMatrix.H"
#include "PhyLib/Phylogeny.H"
#include "FunctionalClass.H"
#include "FunctionalElement.H"
#include "GainLossType.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 A BranchHMM is a PairHMM which also maintains a substitution
 matrix.  Since the substitution matrix is a *conditional*
 distribution and since the PairHMM's emissions form a *joint*
 distribution, this class provides both PairHMM and transducer
 functionality --- i.e., emissions are available as either joint
 or conditional probabilities.
 ****************************************************************/
struct BranchHMM : public PairHMM {
  BranchHMM(const PairHMM &);
  BranchHMM(const PairHMM &,Array1D<SubstitutionMatrix*> &logarithmicMatrices);
  virtual ~BranchHMM();
  inline SubstitutionMatrix *getSubstMatrix(STATE) const;
  inline SubstitutionMatrix *getSubstMatrix(FunctionalClass,FunctionalClass) 
    const;
  BranchHMM *eliminateSilentStates() const;
  inline const Array1D<double> &getEqFreqs(STATE) const; // indexing~PureDnaAlphabet
  inline const Array1D<double> &getBgEqFreqs() const;
  inline void setBgEqFreqs(const Array1D<double> &);
  void printOn(ostream &) const;
  void setFunctionalClass(STATE,BranchEnd,FunctionalClass);
  inline FunctionalClass getFunctionalClass(STATE,BranchEnd) const;
  inline FunctionalElementType getFunctionalElementType(STATE) const;
  inline STATE getState(FunctionalClass parentFC,FunctionalClass childFC,
		 PairHMMStateType);
  void initFuncStateTypeMatrix();
  void initSubstGenerators();
  virtual void convertToLogs();
  void setGainLossFactor(GainLossType,double);
  inline double getGainLossFactor(GainLossType);
  inline bool shouldScoreGainLoss(STATE) const;
  void setGainLossScorePolicy(STATE,bool); // default is true
  inline GainLossType crossFunctionalType(STATE);
  inline const List<STATE> &getStates(FunctionalClass,BranchEnd,
				      PHMM_StateType) const;
  inline const List<STATE> &getBackgroundStates() const 
    { return backgroundStates; }
  inline const List<STATE> &getGainLossStates(FunctionalClass fc,BranchEnd be,
					      PHMM_StateType t) const
    { return gainLossStates(fc,be,t); }
  inline const List<STATE> &getRetentionStates(FunctionalClass fc,
					       PHMM_StateType t) const
    { return retentionStates(fc,t); }
protected:
  Array1D<bool> scoreGainAndLoss; // for each state; initialzed to all true
  Array1D<SubstitutionMatrix*> matrices; // in log-space
  Array2D<SubstitutionMatrix*> matrixByClass; // indexed by FunctionalClass's
  Array1D<Array1D<double> > eqFreqs; // same indexing as PureDnaAlphabet
  Array1D<double> bgEqFreqs;
  Array1D<FunctionalClass> parentFC, childFC;
  Array3D<STATE> classAndTypeToState; // class X class X type -> STATE
  Array3D<List<STATE> > classEndTypeToState; // class X BranchEnd X type -> state
  Array3D<List<STATE> > gainLossStates; // class X branchEnd X type -> state
  Array2D<List<STATE> > retentionStates; // class X type -> state

  List<STATE> backgroundStates; // both parent & child are background
  void initMatrices();
  SubstitutionMatrix *buildMatrix(Transducer &,STATE);
  double gainFactor, lossFactor, retentionFactor, voidFactor;
};

ostream &operator<<(ostream &,const BranchHMM &);


const List<STATE> &BranchHMM::getStates(FunctionalClass fc,BranchEnd end,
					PHMM_StateType t) const
{
  return classEndTypeToState(fc,end,t);
}



SubstitutionMatrix *BranchHMM::getSubstMatrix(STATE q) const
{
  return matrices[q];
}



SubstitutionMatrix *BranchHMM::getSubstMatrix(FunctionalClass fcParent,
					      FunctionalClass fcChild) const
{
  return matrixByClass(int(fcParent),int(fcChild));
}



const Array1D<double> &BranchHMM::getEqFreqs(STATE q) const
{
  return eqFreqs[q];
}



const Array1D<double> &BranchHMM::getBgEqFreqs() const
{
  return bgEqFreqs;
}



void BranchHMM::setBgEqFreqs(const Array1D<double> &f)
{
  bgEqFreqs=f;
  const List<STATE> &bgStates=getBackgroundStates();
  List<STATE>::const_iterator cur=bgStates.begin(), end=bgStates.end();
  for(; cur!=end ; ++cur) eqFreqs[*cur]=bgEqFreqs;
}



FunctionalClass BranchHMM::getFunctionalClass(STATE q,BranchEnd e) const
{
  return e==PARENT ? parentFC[q] : childFC[q];
}



FunctionalElementType BranchHMM::getFunctionalElementType(STATE q) const
{
  FunctionalClass fc=parentFC[q];
  if(fc.fg_or_bg()!=FOREGROUND) fc=childFC[q];
  if(fc.fg_or_bg()==FOREGROUND) return fc.getElementType();
  else return FunctionalElementType::NO_FUNC_ELEM;
}



STATE BranchHMM::getState(FunctionalClass parentFC,FunctionalClass childFC,
			  PairHMMStateType t)
{
  return classAndTypeToState(parentFC,childFC,t);
}



double BranchHMM::getGainLossFactor(GainLossType t)
{
  switch(t)
    {
    case GLT_GAIN:       return gainFactor;
    case GLT_LOSS:       return lossFactor;
    case GLT_RETENTION:  return retentionFactor;
    case GLT_VOID:       return voidFactor;
    }
}



bool BranchHMM::shouldScoreGainLoss(STATE q) const
{
  return scoreGainAndLoss[q];
}



GainLossType BranchHMM::crossFunctionalType(STATE q)
{
  return FunctionalClass::classifyGainLoss(parentFC[q],childFC[q]);
}



#endif

