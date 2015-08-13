/****************************************************************
 LossyFelsenstein.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/Constants.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Exceptions.H"
#include "LossyFelsenstein.H"
#include "Taxon.H"
#include "BranchAttributes.H"
using namespace std;
using namespace BOOM;



inline Taxon &nodeToTaxon(PhylogenyNode &node)
{
  return *static_cast<Taxon*>(node.getDecoration());    
}




/****************************************************************
                   LossyFelsenstein methods
 ****************************************************************/

LossyFelsenstein::LossyFelsenstein(const Alphabet &alphabet,
				   const AlphabetMap &alphabetMap,
				   int numTaxa)
  : LinkFelsenstein(alphabet,alphabetMap,numTaxa), 
    L(numTaxa,2,alphabet.size()), bgFC(FunctionalClass::getBackground()),
    V(alphabet.size())
{
  // ctor

  // deallocate base class' (unused) DP array:
  LinkFelsenstein::L.resize(0,0); 
  L.setAllTo(LOG_0);
}



double LossyFelsenstein::outsideLikelihood(int parentResidueIndex,
					   Taxon &parent,
					   Taxon &badChild,
					   FunctionalClass parentFC,
					   bool scoreGainLoss)
{
  BranchAttributes *branch=parent.getBranchToChild(badChild);
  IndexMap &downMap=branch->getDownMap();
  int &downM=downMap[parentResidueIndex];
  int childResidueIndex=downM;
  downM=IndexMap::UNDEFINED; // only temporary
  double p=
    ancestralLikelihood(parentResidueIndex,parent,parentFC,scoreGainLoss);
  downM=childResidueIndex;
  return p;
}



double LossyFelsenstein::ancestralLikelihood(int residueIndex,Taxon &taxon,
					     FunctionalClass fc,
					     bool scoreGainLoss)
{
  int rootResidueIndex;
  Taxon &root=taxon.findResidueRoot(residueIndex,rootResidueIndex);
  if(&root==&taxon) return logLikelihood(residueIndex,taxon,fc,scoreGainLoss);
  bool constrained=taxon.isFcConstrained();
  if(!constrained) taxon.setFcConstraint(fc);

  // ###DEBUGGING
  //taxon.removeFcConstraint();
  // ###

  double LL=logLikelihood(rootResidueIndex,root,fc,scoreGainLoss);
  if(!constrained) taxon.removeFcConstraint();
  return LL;
}



double LossyFelsenstein::logLikelihood(int residueIndex,Taxon &cladeRoot,
				       FunctionalClass fc,bool scoreGainLoss)
{
  useGainAndLossScores=scoreGainLoss;
  fgFC=fc;
  if(cladeRoot.isFcConstrained()) fc=cladeRoot.getFcConstraint();
  ForegroundBackground fb=fc.fg_or_bg();
  const Array1D<double> &eqFreqs=fc.getEqFreqs();
  postorder(cladeRoot,residueIndex);
  Array1D<double> V(numAlpha);
  int id=cladeRoot.getID();
  Array3D<double>::IndexedTwice<double> row=L[id][fb];
  for(Symbol i=0 ; i<numAlpha; ++i)
    V[i]=(i==gap ? LOG_0 : safeAdd(log(eqFreqs[i]),row[i]));
  return sumLogProbs<double>(V);
}



double LossyFelsenstein::logLikelihood(int residueIndex,Taxon &cladeRoot,
				       bool scoreGainLoss,
				       FunctionalClass foregroundFC)
{
  useGainAndLossScores=scoreGainLoss;
  fgFC=foregroundFC;
  ForegroundBackground fb=fgFC.fg_or_bg();
  const Array1D<double> &eqFreqs=fgFC.getEqFreqs();
  postorder(cladeRoot,residueIndex);
  Array1D<double> V(numAlpha);
  int id=cladeRoot.getID();
  Array3D<double>::IndexedTwice<double> row=L[id][fb];
  for(Symbol i=0 ; i<numAlpha; ++i)
    V[i]=(i==gap ? LOG_0 : safeAdd(log(eqFreqs[i]),row[i]));
  double fgLL=sumLogProbs<double>(V);
  if(fb==BACKGROUND) return fgLL;
  const Array1D<double> &bgEqFreqs=bgFC.getEqFreqs();
  row=L[id][BACKGROUND];
  for(Symbol i=0 ; i<numAlpha; ++i)
    V[i]=(i==gap ? LOG_0 : safeAdd(log(bgEqFreqs[i]),row[i]));
  double bgLL=sumLogProbs<double>(V);
  return max(fgLL,bgLL);
}



void LossyFelsenstein::postorderLeaf(Taxon &taxon,int residueIndex)
{
  int id=taxon.getID();
  Array3D<double>::IndexedOnce<double> row=L[id];
  Symbol a=taxon.getSeq()[residueIndex];//alphabetMap(taxon.getSeq()[residueIndex]); // ###?
  if(a==gap) throw "LossyFelsenstein::postorderLeaf() : gap in leaf sequence";
  Array3D<double>::IndexedTwice<double> fg=row[FOREGROUND];
  Array3D<double>::IndexedTwice<double> bg=row[BACKGROUND];
  if(!taxon.isFcConstrained())
    for(Symbol i=0 ; i<numAlpha; ++i)
      fg[i]=bg[i]=(i==a ? LOG_1 : LOG_0);
  else {
    FunctionalClass fc=taxon.getFcConstraint();
    if(fc.isBackground())
      for(Symbol i=0 ; i<numAlpha; ++i) {
	bg[i]=(i==a ? LOG_1 : LOG_0);
	fg[i]=LOG_0;
      }
    else
      for(Symbol i=0 ; i<numAlpha; ++i) {
	fg[i]=(i==a ? LOG_1 : LOG_0);
	bg[i]=LOG_0;
      }
  }
}



void LossyFelsenstein::postorderInternal(Taxon &taxon,int residueIndex)
{
  // First, recurse to the children (i.e., fill in their DP matrix
  // entries, and all those below them in the tree)
  BranchAttributes &leftBranch=*taxon.getBranch(LEFT);
  BranchAttributes &rightBranch=*taxon.getBranch(RIGHT);
  recurseToChild(leftBranch,residueIndex);
  recurseToChild(rightBranch,residueIndex);

  // Now process this node
  int id=taxon.getID();
  Taxon &leftChild=leftBranch.getChildTaxon();
  Taxon &rightChild=rightBranch.getChildTaxon();
  Array3D<double>::IndexedOnce<double> row=L[id];
  Array3D<double>::IndexedTwice<double> parentFG=row[FOREGROUND];
  Array3D<double>::IndexedTwice<double> parentBG=row[BACKGROUND];
  if(!taxon.isFcConstrained())
    for(Symbol a=0 ; a<numAlpha ; ++a) {
      if(a==gap) { parentFG[a]=parentBG[a]=LOG_0; continue; }
      double leftFgBg=processInternalChild(a,leftBranch,fgFC,bgFC);
      double leftFgFg=processInternalChild(a,leftBranch,fgFC,fgFC);
      double leftBgBg=processInternalChild(a,leftBranch,bgFC,bgFC);
      double rightFgBg=processInternalChild(a,rightBranch,fgFC,bgFC);
      double rightFgFg=processInternalChild(a,rightBranch,fgFC,fgFC);
      double rightBgBg=processInternalChild(a,rightBranch,bgFC,bgFC);
      parentFG[a]=safeAdd(max(leftFgFg,leftFgBg),max(rightFgFg,rightFgBg));
      parentBG[a]=safeAdd(leftBgBg,rightBgBg);
    }
  else {
    FunctionalClass fc=taxon.getFcConstraint();
    if(fc.isBackground())
      for(Symbol a=0 ; a<numAlpha ; ++a) {
	if(a==gap) { parentFG[a]=parentBG[a]=LOG_0; continue; }
	double leftBgBg=processInternalChild(a,leftBranch,bgFC,bgFC);
	double rightBgBg=processInternalChild(a,rightBranch,bgFC,bgFC);
	parentBG[a]=safeAdd(leftBgBg,rightBgBg);
	parentFG[a]=LOG_0;
      }
    else
      for(Symbol a=0 ; a<numAlpha ; ++a) {
	if(a==gap) { parentFG[a]=parentBG[a]=LOG_0; continue; }
	double leftFgBg=processInternalChild(a,leftBranch,fgFC,bgFC);
	double leftFgFg=processInternalChild(a,leftBranch,fgFC,fgFC);
	double rightFgBg=processInternalChild(a,rightBranch,fgFC,bgFC);
	double rightFgFg=processInternalChild(a,rightBranch,fgFC,fgFC);
	parentFG[a]=safeAdd(max(leftFgFg,leftFgBg),max(rightFgFg,rightFgBg));
	parentBG[a]=LOG_0;
      }
  }
}



int LossyFelsenstein::recurseToChild(BranchAttributes &branch,int parentIndex)
{
  const IndexMap &indexMap=branch.getDownMap();
  int childIndex=indexMap[parentIndex];
  Taxon &child=branch.getChildTaxon();
  if(childIndex==IndexMap::UNDEFINED) {
    Array3D<double>::IndexedOnce<double> row=L[child.getID()];
    Array3D<double>::IndexedTwice<double> fg=row[FOREGROUND], 
      bg=row[BACKGROUND];
    for(Symbol a=0 ; a<numAlpha ; ++a) 
      fg[a]=bg[a]=(a==gap ? LOG_0 : LOG_1);
  }
  else postorder(child,childIndex);
  return childIndex;
}



void LossyFelsenstein::postorderRoot(Taxon &taxon,int residueIndex)
{
  // First, visit the child subtree
  BranchAttributes *branch=taxon.getBranch(UNIQUE);
  recurseToChild(*branch,residueIndex);

  // Now process this node
  int id=taxon.getID();
  Array3D<double>::IndexedOnce<double> row=L[id];
  Array3D<double>::IndexedTwice<double> parentFG=row[FOREGROUND];
  Array3D<double>::IndexedTwice<double> parentBG=row[BACKGROUND];
  Taxon &child=branch->getChildTaxon();
  int childID=child.getID();
  if(!taxon.isFcConstrained())
    for(Symbol a=0 ; a<numAlpha ; ++a) {
      if(a==gap) { parentFG[a]=parentBG[a]=LOG_0; continue; }
      double childFgBg=processInternalChild(a,*branch,fgFC,bgFC);
      double childFgFg=processInternalChild(a,*branch,fgFC,fgFC);
      double childBgBg=processInternalChild(a,*branch,bgFC,bgFC);
      parentFG[a]=max(childFgFg,childFgBg);
      parentBG[a]=childBgBg;
    }
  else {
    FunctionalClass fc=taxon.getFcConstraint();
    if(fc.isBackground())
      for(Symbol a=0 ; a<numAlpha ; ++a) {
	if(a==gap) { parentFG[a]=parentBG[a]=LOG_0; continue; }
	parentBG[a]=processInternalChild(a,*branch,bgFC,bgFC);
	parentFG[a]=LOG_0;
      }
    else
      for(Symbol a=0 ; a<numAlpha ; ++a) {
	if(a==gap) { parentFG[a]=parentBG[a]=LOG_0; continue; }
	double childFgBg=processInternalChild(a,*branch,fgFC,bgFC);
	double childFgFg=processInternalChild(a,*branch,fgFC,fgFC);
	parentFG[a]=max(childFgFg,childFgBg);
	parentBG[a]=LOG_0;
      }
  }
}



double LossyFelsenstein::processInternalChild(Symbol parentSymbol,
					      BranchAttributes &branch,
					      FunctionalClass parentFC,
					      FunctionalClass childFC)
{
  BranchHMM *hmm=branch.getHMM();
  if(!hmm->getSubstMatrix(parentFC,childFC)) return LOG_0;
  Taxon &child=branch.getChildTaxon();
  if(child.isFcConstrained()) childFC=child.getFcConstraint();
  SubstitutionMatrix &Pt=*hmm->getSubstMatrix(parentFC,childFC);
  GainLossType gainLossType=
    FunctionalClass::classifyGainLoss(parentFC,childFC);
  double gainLossFactor=
    useGainAndLossScores && focalBranch1!=&branch && focalBranch2!=&branch
    ? hmm->getGainLossFactor(gainLossType) : 0;
  if(isinf(gainLossFactor)) {
    cout<<"gainLossFactor="<<gainLossFactor<<" use="<<useGainAndLossScores<<" lookup="<<hmm->getGainLossFactor(gainLossType)<<" type="<<gainLossType<<endl;
    throw 0;
  }
  ForegroundBackground FB=childFC.fg_or_bg();
  Array3D<double>::IndexedTwice<double> row=L[child.getID()][FB];
  Array1D<double> V(numAlpha);
  for(Symbol b=0 ; b<numAlpha ; ++b) {
    Symbol c=alphabetMap(b), p=alphabetMap(parentSymbol);
    V[b]=(b==gap) ? LOG_0 : 
      safeAdd(Pt(p,c),row[b],gainLossFactor);
  }
  return sumLogProbs<double>(V);
}



float LossyFelsenstein::precompInsideChild(Symbol parentSymbol,
					   Taxon &child,
					   IndexMap &indexMap,
					   BranchHMM *hmm,
					   FunctionalClass parentFC,
					   FunctionalClass childFC,
					   PrecomputedEmissions &peChild,
					   int childPos)
{
  if(child.isFcConstrained()) childFC=child.getFcConstraint();
  SubstitutionMatrix *m=hmm->getSubstMatrix(parentFC,childFC);
  if(!m) return LOG_0;
  SubstitutionMatrix &Pt=*m;
  if(childPos==IndexMap::UNDEFINED) return LOG_1;
  Array3D<float>::IndexedTwice<float> row=peChild[childPos][childFC];
  for(Symbol b=0 ; b<numAlpha ; ++b) {
    Symbol c=alphabetMap(b), p=alphabetMap(parentSymbol);
    V[b]=(b==gap) ? LOG_0 : safeAdd(Pt(p,c),row[b]);
  }
  return sumLogProbs<double>(V);
}



void LossyFelsenstein::precomputeInside(Taxon &taxon)
{
#define NEWLOSSYFELS 1
#ifdef NEWLOSSYFELS

  int L=taxon.getSeqLen();
  int numFC=FunctionalClass::numClasses();
  PrecomputedEmissions &peParent=taxon.getPrecomputedEmissions();
  if(taxon.isLeaf()) {
    for(int pos=0 ; pos<L ; ++pos) {
      Array3D<float>::IndexedOnce<float> slice=peParent[pos];
      Symbol a=alphabetMap(taxon.getSeq()[pos]);
      for(int fc=0 ; fc<numFC ; ++fc) {
	Array3D<float>::IndexedTwice<float> row=slice[fc];
	for(Symbol i=0 ; i<numAlpha; ++i) 
	  row[i]=(i==a ? LOG_1 : LOG_0);
      }
    }
  }
  //else if(taxon.isRoot()) INTERNAL_ERROR
  else {
    BranchAttributes &leftBranch=*taxon.getBranch(LEFT);
    BranchAttributes &rightBranch=*taxon.getBranch(RIGHT);
    Taxon &leftChild=leftBranch.getChildTaxon();
    Taxon &rightChild=rightBranch.getChildTaxon();
    IndexMap &leftIndexMap=leftBranch.getDownMap();
    IndexMap &rightIndexMap=rightBranch.getDownMap();
    BranchHMM *leftHMM=leftBranch.getHMM(), *rightHMM=rightBranch.getHMM();
    PrecomputedEmissions &peLeft=leftChild.getPrecomputedEmissions();
    PrecomputedEmissions &peRight=rightChild.getPrecomputedEmissions();
    for(int pos=0 ; pos<L ; ++pos) {
      int leftPos=leftIndexMap[pos];
      int rightPos=rightIndexMap[pos];
      Array3D<float>::IndexedOnce<float> slice=peParent[pos];
      for(int fc=0 ; fc<numFC ; ++fc) {
	const FunctionalClass fgFC=fc;
	ForegroundBackground fb=fgFC.fg_or_bg();
	Array3D<float>::IndexedTwice<float> row=slice[fc];
	if(!taxon.isFcConstrained())
	  for(Symbol a=0 ; a<numAlpha ; ++a) {
	    if(fb==FOREGROUND) {
	      float leftFgBg=precompInsideChild(a,leftChild,leftIndexMap,
						leftHMM,fgFC,bgFC,peLeft,
						leftPos);
	      float leftFgFg=precompInsideChild(a,leftChild,leftIndexMap,
						leftHMM,fgFC,fgFC,peLeft,
						leftPos);
	      float rightFgBg=precompInsideChild(a,rightChild,rightIndexMap,
						 rightHMM,fgFC,bgFC,peRight,
						 rightPos);
	      float rightFgFg=precompInsideChild(a,rightChild,rightIndexMap,
						 rightHMM,fgFC,fgFC,peRight,
						 rightPos);
	      row[a]=safeAdd(max(leftFgFg,leftFgBg),max(rightFgFg,rightFgBg));
	    }
	    else {
	      float leftBgBg=precompInsideChild(a,leftChild,leftIndexMap,
						leftHMM,bgFC,bgFC,peLeft,
						leftPos);
	      float rightBgBg=precompInsideChild(a,rightChild,rightIndexMap,
						 rightHMM,bgFC,bgFC,peRight,
						 rightPos);
	      row[a]=safeAdd(leftBgBg,rightBgBg);
	    }
	  }
	else {
	  FunctionalClass fc=taxon.getFcConstraint();
	  if(fc.isBackground())
	    for(Symbol a=0 ; a<numAlpha ; ++a) {
	      float leftBgBg=precompInsideChild(a,leftChild,leftIndexMap,
						leftHMM,bgFC,bgFC,peLeft,
						leftPos);
	      float rightBgBg=precompInsideChild(a,rightChild,rightIndexMap,
						 rightHMM,bgFC,bgFC,peRight,
						 rightPos);
	      row[a]=safeAdd(leftBgBg,rightBgBg);
	    }
	  else
	    for(Symbol a=0 ; a<numAlpha ; ++a) {
	      float leftFgBg=precompInsideChild(a,leftChild,leftIndexMap,
						leftHMM,fgFC,bgFC,peLeft,
						leftPos);
	      float leftFgFg=precompInsideChild(a,leftChild,leftIndexMap,
						leftHMM,fgFC,fgFC,peLeft,
						leftPos);
	      float rightFgBg=precompInsideChild(a,rightChild,rightIndexMap,
						 rightHMM,fgFC,bgFC,peRight,
						 rightPos);
	      float rightFgFg=precompInsideChild(a,rightChild,rightIndexMap,
						 rightHMM,fgFC,fgFC,peRight,
						 rightPos);
	      row[a]=safeAdd(max(leftFgFg,leftFgBg),max(rightFgFg,rightFgBg));
	    }
	}
      }
    }
  }
  
#else
  int len=taxon.getSeqLen();
  int numFC=FunctionalClass::numClasses();
  for(int pos=0 ; pos<len ; ++pos) {
    for(int fc=0 ; fc<numFC ; ++fc) {
      PrecomputedEmissions &peParent=taxon.getPrecomputedEmissions();
      Array3D<float>::IndexedTwice<float> peRow=peParent[pos][fc];
      fgFC=fc;
      ForegroundBackground fb=fgFC.fg_or_bg();
      postorder(taxon,pos);
      Array3D<double>::IndexedOnce<double> row=L[taxon.getID()];
      for(Symbol i=0 ; i<numAlpha ; ++i)
	peRow[i]=row[fb][i];
    }
  }
#endif
}



void LossyFelsenstein::precomputeOutside(Array1D<float> &dp,Taxon &taxon,
					 int pos,Taxon &badChild,
					 FunctionalClass fc)
{
  if(fc.isValid()) {
    fgFC=fc;
    ForegroundBackground fb=fc.fg_or_bg();
    outorder(taxon,pos,badChild);
    Array3D<double>::IndexedOnce<double> row=L[taxon.getID()];
    for(Symbol i=0 ; i<numAlpha ; ++i)
      dp[i]=row[fb][i];
  }
  else dp.setAllTo(NEGATIVE_INFINITY);
}



void LossyFelsenstein::outorderInternal(Taxon &taxon,int residueIndex,
					Taxon &badChild)
{
  int id=taxon.getID();
  Array3D<double>::IndexedOnce<double> row=L[id];
  Array3D<double>::IndexedTwice<double> rowFG=row[FOREGROUND];
  Array3D<double>::IndexedTwice<double> rowBG=row[BACKGROUND];
  BranchAttributes *leftBranch=taxon.getBranch(LEFT);
  BranchAttributes *rightBranch=taxon.getBranch(RIGHT);
  BranchAttributes *parentBranch=taxon.getBranchToParent();
  BranchAttributes *childBranch=
    &leftBranch->getChildTaxon()==&badChild ? rightBranch : leftBranch;
  Taxon &child=childBranch->getChildTaxon();
  Taxon *parent=parentBranch ? &parentBranch->getParentTaxon() : NULL;
  double parentFgBg, parentFgFg, parentBgBg;
  for(Symbol a=0 ; a<numAlpha ; ++a) {
    //if(a==gap) INTERNAL_ERROR; // ### DEBUGGING
    if(parent && parent->getSeqLen()>0) {
      //INTERNAL_ERROR; // ### DEBUGGING
      parentFgBg=processParent(a,*parentBranch,fgFC,bgFC);
      parentFgFg=processParent(a,*parentBranch,fgFC,fgFC);
      parentBgBg=processParent(a,*parentBranch,bgFC,bgFC);
    }
    else parentFgBg=parentFgFg=parentBgBg=LOG_1;
    if(!finite(parentFgBg) && !finite(parentFgFg) && !finite(parentBgBg))
      INTERNAL_ERROR;
    double childFgBg=processInternalChild(a,*childBranch,fgFC,bgFC);
    double childFgFg=processInternalChild(a,*childBranch,fgFC,fgFC);
    double childBgBg=processInternalChild(a,*childBranch,bgFC,bgFC);
    rowFG[a]=safeAdd(parentFgFg,max(childFgFg,childFgBg));
    rowBG[a]=safeAdd(max(parentFgBg,parentBgBg),childBgBg);
  }
}



double LossyFelsenstein::processParent(Symbol childSymbol,
					  BranchAttributes &branch,
					  FunctionalClass parentFC,
					  FunctionalClass childFC)
{
  Taxon &parent=branch.getParentTaxon();
  BranchHMM *hmm=branch.getHMM();
  SubstitutionMatrix &Pt=*hmm->getSubstMatrix(parentFC,childFC);
  if(!&Pt) return LOG_0;
  GainLossType gainLossType=
    FunctionalClass::classifyGainLoss(parentFC,childFC);
  double gainLossFactor=
    useGainAndLossScores && focalBranch1!=&branch && focalBranch2!=&branch
    ? hmm->getGainLossFactor(gainLossType) : 0;
  if(isinf(gainLossFactor)) {
    cout<<"gainLossFactor="<<gainLossFactor<<" use="<<useGainAndLossScores<<" lookup="<<hmm->getGainLossFactor(gainLossType)<<" type="<<gainLossType<<endl;
    INTERNAL_ERROR;
  }
  ForegroundBackground FB=parentFC.fg_or_bg();
  Array3D<double>::IndexedTwice<double> row=L[parent.getID()][FB];
  Array1D<double> V(numAlpha);
  for(Symbol b=0 ; b<numAlpha ; ++b) {
    Symbol c=alphabetMap(b), p=alphabetMap(childSymbol);
    V[b]=(b==gap) ? LOG_0 : 
      safeAdd(Pt(p,c),row[b],gainLossFactor);
  }
  return sumLogProbs<double>(V);
}

