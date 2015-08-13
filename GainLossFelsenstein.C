/****************************************************************
 GainLossFelsenstein.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/Constants.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Exceptions.H"
#include "GainLossFelsenstein.H"
#include "Taxon.H"
#include "BranchAttributes.H"
using namespace std;
using namespace BOOM;



inline Taxon &nodeToTaxon(PhylogenyNode &node)
{
  return *static_cast<Taxon*>(node.getDecoration());    
}




/****************************************************************
                   GainLossFelsenstein methods
 ****************************************************************/

GainLossFelsenstein::GainLossFelsenstein(const Alphabet &alphabet,
				   const AlphabetMap &alphabetMap,
				   int numTaxa)
  : LinkFelsenstein(alphabet,alphabetMap,numTaxa), 
    L(numTaxa,2,alphabet.size()), bgFC(FunctionalClass::getBackground())
{
  // ctor

  // deallocate base class' (unused) DP array:
  LinkFelsenstein::L.resize(0,0); 
  L.setAllTo(LOG_1);
}



double GainLossFelsenstein::outsideLikelihood(int parentResidueIndex,
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



double GainLossFelsenstein::ancestralLikelihood(int residueIndex,Taxon &taxon,
					     FunctionalClass fc,
					     bool scoreGainLoss)
{
  int rootResidueIndex;
  Taxon &root=taxon.findResidueRoot(residueIndex,rootResidueIndex);
  if(&root==&taxon) return logLikelihood(residueIndex,taxon,fc,scoreGainLoss);
  bool constrained=taxon.isFcConstrained();

  //constrained=scoreGainLoss=false; // ###

  if(!constrained) taxon.setFcConstraint(fc);

  // ### the fc in this call should be for the ancestral residue:
  double LL=logLikelihood(rootResidueIndex,root,fc,scoreGainLoss);
  if(!constrained) taxon.removeFcConstraint();

  return LL;
}



double GainLossFelsenstein::logLikelihood(int residueIndex,Taxon &cladeRoot,
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



double GainLossFelsenstein::logLikelihood(int residueIndex,Taxon &cladeRoot,
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



void GainLossFelsenstein::postorderLeaf(Taxon &taxon,int residueIndex)
{
  int id=taxon.getID();
  Array3D<double>::IndexedOnce<double> row=L[id];
  Symbol a=taxon.getSeq()[residueIndex];//alphabetMap(taxon.getSeq()[residueIndex]); // ###?
  if(a==gap) throw "GainLossFelsenstein::postorderLeaf() : gap in leaf sequence";
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



void GainLossFelsenstein::postorderInternal(Taxon &taxon,int residueIndex)
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
      //if(a==gap) { parentFG[a]=parentBG[a]=LOG_0; continue; }
      double leftFgBg=processInternalChild(a,leftBranch,fgFC,bgFC);
      double leftFgFg=processInternalChild(a,leftBranch,fgFC,fgFC);
      double leftBgBg=processInternalChild(a,leftBranch,bgFC,bgFC);
      double leftBgFg=processInternalChild(a,leftBranch,bgFC,fgFC);
      double rightFgBg=processInternalChild(a,rightBranch,fgFC,bgFC);
      double rightFgFg=processInternalChild(a,rightBranch,fgFC,fgFC);
      double rightBgBg=processInternalChild(a,rightBranch,bgFC,bgFC);
      double rightBgFg=processInternalChild(a,rightBranch,bgFC,fgFC);
      parentFG[a]=safeAdd(max(leftFgFg,leftFgBg),max(rightFgFg,rightFgBg));
      parentBG[a]=safeAdd(max(leftBgBg,leftBgFg),max(rightBgBg,rightBgFg));
    }
  else {
    FunctionalClass fc=taxon.getFcConstraint();
    if(fc.isBackground())
      for(Symbol a=0 ; a<numAlpha ; ++a) {
	//if(a==gap) { parentFG[a]=parentBG[a]=LOG_0; continue; }
	double leftBgBg=processInternalChild(a,leftBranch,bgFC,bgFC);
	double leftBgFg=processInternalChild(a,leftBranch,bgFC,fgFC);
	double rightBgBg=processInternalChild(a,rightBranch,bgFC,bgFC);
	double rightBgFg=processInternalChild(a,rightBranch,bgFC,fgFC);
	parentBG[a]=safeAdd(max(leftBgBg,leftBgFg),max(rightBgBg,rightBgFg));
	parentFG[a]=LOG_0;
      }
    else
      for(Symbol a=0 ; a<numAlpha ; ++a) {
	//if(a==gap) { parentFG[a]=parentBG[a]=LOG_0; continue; }
	double leftFgBg=processInternalChild(a,leftBranch,fgFC,bgFC);
	double leftFgFg=processInternalChild(a,leftBranch,fgFC,fgFC);
	double rightFgBg=processInternalChild(a,rightBranch,fgFC,bgFC);
	double rightFgFg=processInternalChild(a,rightBranch,fgFC,fgFC);
	parentFG[a]=safeAdd(max(leftFgFg,leftFgBg),max(rightFgFg,rightFgBg));
	parentBG[a]=LOG_0;
      }
  }
}



int GainLossFelsenstein::recurseToChild(BranchAttributes &branch,
					int parentIndex)
{
  const IndexMap &indexMap=branch.getDownMap();
  int childIndex=indexMap[parentIndex];
  Taxon &child=branch.getChildTaxon();
  if(childIndex==IndexMap::UNDEFINED) {
    Array3D<double>::IndexedOnce<double> row=L[child.getID()];
    Array3D<double>::IndexedTwice<double> fg=row[FOREGROUND], 
      bg=row[BACKGROUND];
    for(Symbol a=0 ; a<numAlpha ; ++a) 
      fg[a]=bg[a]=LOG_1;//(a==gap ? LOG_0 : LOG_1); // ###

    /*
    {//###
      bool found=false;
      for(Symbol a=0 ; a<numAlpha ; ++a) 
	if(isFinite(bg[a])) found=true;
      if(!found) {
	for(Symbol a=0 ; a<numAlpha ; ++a) 
	  cout<<"fg["<<int(a)<<"]="<<fg[a]<<" bg["<<int(a)<<"="<<bg[a]<<endl;
	INTERNAL_ERROR;
      }
    }//###
    */
  }
  else postorder(child,childIndex);
  /*
    {//###
    Array3D<double>::IndexedOnce<double> row=L[child.getID()];
    Array3D<double>::IndexedTwice<double> fg=row[FOREGROUND], 
      bg=row[BACKGROUND];
      bool found=false;
      for(Symbol a=0 ; a<numAlpha ; ++a) 
	if(isFinite(bg[a])) found=true;
      if(!found) {
	for(Symbol a=0 ; a<numAlpha ; ++a) 
	  cout<<"fg["<<int(a)<<"]="<<fg[a]<<" bg["<<int(a)<<"="<<bg[a]<<endl;
	INTERNAL_ERROR;
      }
      }//###
  */
  return childIndex;
}



void GainLossFelsenstein::postorderRoot(Taxon &taxon,int residueIndex)
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



double GainLossFelsenstein::processInternalChild(Symbol parentSymbol,
					      BranchAttributes &branch,
					      FunctionalClass parentFC,
					      FunctionalClass childFC)
{
  Taxon &child=branch.getChildTaxon();
  if(child.isFcConstrained()) childFC=child.getFcConstraint();
  BranchHMM *hmm=branch.getHMM();
  SubstitutionMatrix *pt=hmm->getSubstMatrix(parentFC,childFC);
  if(!pt) {INTERNAL_ERROR; return LOG_0;} // ###
  SubstitutionMatrix &Pt=*pt;
  GainLossType gainLossType=
    FunctionalClass::classifyGainLoss(parentFC,childFC);
  double gainLossFactor=
    useGainAndLossScores && focalBranch1!=&branch && focalBranch2!=&branch
    ? hmm->getGainLossFactor(gainLossType) : 0;
  if(!isFinite(gainLossFactor)) {
    cout<<"gainLossFactor="<<gainLossFactor<<" use="<<useGainAndLossScores<<" lookup="<<hmm->getGainLossFactor(gainLossType)<<" type="<<gainLossType<<endl;
    INTERNAL_ERROR;
  }
  ForegroundBackground FB=childFC.fg_or_bg();
  Array3D<double>::IndexedTwice<double> row=L[child.getID()][FB];
  Array1D<double> V(numAlpha);
  for(Symbol b=0 ; b<numAlpha ; ++b) {
    Symbol c=alphabetMap(b), p=alphabetMap(parentSymbol);
    V[b]=//(b==gap) ? LOG_0 : //###
      safeAdd(Pt(p,c),row[b],gainLossFactor);
  }
  return sumLogProbs<double>(V);
}



void GainLossFelsenstein::precomputeInside(Array1D<float> &dp,Taxon &taxon,
					   int pos,FunctionalClass fc)
{
  if(fc.isValid()) {
    fgFC=fc;
    ForegroundBackground fb=fc.fg_or_bg();
    postorder(taxon,pos);
    Array3D<double>::IndexedOnce<double> row=L[taxon.getID()];
    for(Symbol i=0 ; i<numAlpha ; ++i)
      dp[i]=row[fb][i];
  }
  else dp.setAllTo(NEGATIVE_INFINITY);
}



void GainLossFelsenstein::precomputeOutside(Array1D<float> &dp,Taxon &taxon,
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
    /*
    if(taxon.getName()=="A3" && pos==0 && fc.isBackground() && !isFinite(dp[0]) && !isFinite(dp[1]) && !isFinite(dp[2]) && !isFinite(dp[3])) {
      cout<<"FFF "<<"taxonID="<<taxon.getID()<<" "<<L<<endl;
      INTERNAL_ERROR;
    }
    */
  }
  else dp.setAllTo(NEGATIVE_INFINITY);
}



void GainLossFelsenstein::outorderInternal(Taxon &taxon,int residueIndex,
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
  double parentFgBg, parentFgFg, parentBgFg, parentBgBg;
  //bool found=false;//###DEBUGGING
  for(Symbol a=0 ; a<numAlpha ; ++a) {
    //if(a==gap) INTERNAL_ERROR; // ### DEBUGGING
    if(parentBranch) {
      parentFgBg=processParent(a,*parentBranch,fgFC,bgFC);
      parentFgFg=processParent(a,*parentBranch,fgFC,fgFC);
      parentBgFg=processParent(a,*parentBranch,bgFC,fgFC);
      parentBgBg=processParent(a,*parentBranch,bgFC,bgFC);
    }
    else parentFgBg=parentFgFg=parentBgFg=parentBgBg=LOG_1;
    double childFgBg=processInternalChild(a,*childBranch,fgFC,bgFC);
    double childFgFg=processInternalChild(a,*childBranch,fgFC,fgFC);
    double childBgBg=processInternalChild(a,*childBranch,bgFC,bgFC);
    double childBgFg=processInternalChild(a,*childBranch,bgFC,fgFC);
    rowFG[a]=safeAdd(max(parentFgFg,parentBgFg),max(childFgFg,childFgBg));
    rowBG[a]=safeAdd(max(parentFgBg,parentBgBg),max(childBgFg,childBgBg));
    //if(isFinite(rowFG[a]) || isFinite(rowBG[a])) found=true;//DEBUGGING
  }
  /*
  if(!found) {
    cout<<" taxon="<<taxon.getName()<<" taxonID="<<id<<" parent="<<(parentBranch ? parentBranch->getParentTaxon().getName() : "NONE")<<" child="<<(childBranch ? childBranch->getChildTaxon().getName() : "NONE")<<" badChild="<<badChild.getName()<<" residueIndex="<<residueIndex<<endl;
    for(Symbol a=0 ; a<numAlpha ; ++a) {
      if(parentBranch) {
	parentFgBg=processParent(a,*parentBranch,fgFC,bgFC);
	parentFgFg=processParent(a,*parentBranch,fgFC,fgFC);
	parentBgFg=processParent(a,*parentBranch,bgFC,fgFC);
	parentBgBg=processParent(a,*parentBranch,bgFC,bgFC);
      }
      else parentFgBg=parentFgFg=parentBgFg=parentBgBg=LOG_1;
      double childFgBg=processInternalChild(a,*childBranch,fgFC,bgFC);
      double childFgFg=processInternalChild(a,*childBranch,fgFC,fgFC);
      double childBgBg=processInternalChild(a,*childBranch,bgFC,bgFC);
      double childBgFg=processInternalChild(a,*childBranch,bgFC,fgFC);
      cout<<"childFgBg="<<childFgBg<<" childFgFg="<<childFgFg<<" childBgBg="<<childBgBg<<" childBgFg="<<childBgFg<<" parentFgBg="<<parentFgBg<<" parentFgFg="<<parentFgFg<<" parentBgBg="<<parentBgBg<<" parentBgFg="<<parentBgFg<<endl;
    }
    INTERNAL_ERROR;//DEBUGGING
  }
  */
}



double GainLossFelsenstein::processParent(Symbol childSymbol,
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
