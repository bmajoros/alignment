/****************************************************************
 LinkFelsenstein.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/Array2D.H"
#include "BOOM/Constants.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Exceptions.H"
#include "LinkFelsenstein.H"
#include "Taxon.H"
#include "BranchAttributes.H"
#include "ResidueAddress.H"
using namespace std;
using namespace BOOM;


const double log1=0.0;
const double log0=NEGATIVE_INFINITY;



inline Taxon &nodeToTaxon(PhylogenyNode &node)
{
  return *static_cast<Taxon*>(node.getDecoration());    
}




/****************************************************************
                   LinkFelsenstein methods
 ****************************************************************/

LinkFelsenstein::LinkFelsenstein(const Alphabet &alphabet,
				 const AlphabetMap &alphabetMap,
				 int numTaxa)
  : alphabet(alphabet), numAlpha(alphabet.size()),
    alphabetMap(alphabetMap), numTaxa(numTaxa), 
    gap(alphabet.lookup('-')), L(numTaxa,alphabet.size()),
    focalBranch1(NULL), focalBranch2(NULL),
    useGainAndLossScores(false) // ###
{
  // ctor
}



double LinkFelsenstein::outsideLikelihood(int parentResidueIndex,
					  Taxon &parent,Taxon &badChild,
					  bool scoreGainLoss)
{
  BranchAttributes *branch=parent.getBranchToChild(badChild);
  IndexMap &downMap=branch->getDownMap();
  int childResidueIndex=downMap[parentResidueIndex];
  downMap[parentResidueIndex]=IndexMap::UNDEFINED; // only temporary
  double p=ancestralLikelihood(parentResidueIndex,parent,scoreGainLoss);
  downMap[parentResidueIndex]=childResidueIndex;
  return p;
}



double LinkFelsenstein::outsideLikelihood(int parentResidueIndex,
					  Taxon &parent,Taxon &badChild,
					  FunctionalClass fc,
					  bool scoreGainLoss)
{
  BranchAttributes *branch=parent.getBranchToChild(badChild);
  IndexMap &downMap=branch->getDownMap();
  int childResidueIndex=downMap[parentResidueIndex];
  downMap[parentResidueIndex]=IndexMap::UNDEFINED; // only temporary
  double p=ancestralLikelihood(parentResidueIndex,parent,fc,scoreGainLoss);
  downMap[parentResidueIndex]=childResidueIndex;
  return p;
}



double LinkFelsenstein::ancestralLikelihood(int residueIndex,Taxon &taxon,
					    bool scoreGainLoss)
{
  int rootResidueIndex;
  Taxon &root=taxon.findResidueRoot(residueIndex,rootResidueIndex);
  return logLikelihood(rootResidueIndex,root,scoreGainLoss);
}



double LinkFelsenstein::ancestralLikelihood(int residueIndex,Taxon &taxon,
					    FunctionalClass fc,
					    bool scoreGainLoss)
{
  int rootResidueIndex;
  Taxon &root=taxon.findResidueRoot(residueIndex,rootResidueIndex);
  if(&root==&taxon)
    return logLikelihood(rootResidueIndex,root,fc,scoreGainLoss);
  else
    return logLikelihood(rootResidueIndex,root,scoreGainLoss);
}



double LinkFelsenstein::logLikelihood(int residueIndex,Taxon &cladeRoot,
				      bool scoreGainLoss)
{
  FunctionalClass fc=cladeRoot.getFunctionalParse()[residueIndex];
  return logLikelihood(residueIndex,cladeRoot,fc,scoreGainLoss);
}



double LinkFelsenstein::logLikelihood(int residueIndex,Taxon &cladeRoot,
				      FunctionalClass fc,
				      bool scoreGainLoss)
{
  useGainAndLossScores=scoreGainLoss;
  const Array1D<double> &eqFreqs=fc.getEqFreqs(); // ### fixed: 9/11/09
  postorder(cladeRoot,residueIndex);
  Array1D<double> V(numAlpha);
  int id=cladeRoot.getID();
  Array2D<double>::RowIn2DArray<double> row=L[id];
  for(Symbol i=0 ; i<numAlpha; ++i)
    V[i]=safeAdd(log(eqFreqs[i]),row[i]);
  return sumLogProbs<double>(V);
}



void LinkFelsenstein::postorder(Taxon &taxon,int residueIndex)
{
  switch(taxon.getNodeType()) 
    {
    case ROOT_NODE:     postorderRoot(taxon,residueIndex); break;
    case LEAF_NODE:     postorderLeaf(taxon,residueIndex); break;
    case INTERNAL_NODE: postorderInternal(taxon,residueIndex); break;
    }
}



void LinkFelsenstein::postorderLeaf(Taxon &taxon,int residueIndex)
{
  cout<<"postorderLeaf"<<endl;
  int id=taxon.getID();
  cout<<"a "<<id<<endl;
  Symbol a=alphabetMap(taxon.getSeq()[residueIndex]);
  cout<<"b "<<int(a)<<endl;
  Array2D<double>::RowIn2DArray<double> row=L[id];
  cout<<"c"<<endl;
  if(a==gap) throw "LinkFelsenstein::postorderLeaf() : gap in leaf sequence";
  for(Symbol i=0 ; i<numAlpha; ++i) {
    cout<<"d"<<endl;
    row[i]=(i==a ? log1 : log0);
    cout<<"e"<<endl;
  }
  cout<<"f"<<endl;
}



void LinkFelsenstein::postorderInternal(Taxon &taxon,int residueIndex)
{
  cout<<"postorderInternal"<<endl;
  // First, recurse to the children
  BranchAttributes &leftBranch=*taxon.getBranch(LEFT);
  BranchAttributes &rightBranch=*taxon.getBranch(RIGHT);
  cout<<"int1"<<endl;
  int leftIndex=recurseToChild(leftBranch,residueIndex);
  cout<<"int2"<<endl;
  int rightIndex=recurseToChild(rightBranch,residueIndex);
  cout<<"int3"<<endl;

  // Now process this node
  int id=taxon.getID();
  Taxon &leftChild=leftBranch.getChildTaxon();
  Taxon &rightChild=rightBranch.getChildTaxon();
  Array2D<double>::RowIn2DArray<double> row=L[id];
  FunctionalClass fcParent=taxon.getFunctionalParse()[residueIndex];
  cout<<"int4"<<endl;
  FunctionalClass fcLeft=
    (leftIndex==IndexMap::UNDEFINED ? fcParent : 
     leftChild.getFunctionalParse()[leftIndex]);
  cout<<"int4.1"<<endl;
  FunctionalClass fcRight=
    (rightIndex==IndexMap::UNDEFINED ? fcParent : 
     rightChild.getFunctionalParse()[rightIndex]);
  cout<<"int4.2 "<<fcParent<<" "<<fcLeft<<endl;
  SubstitutionMatrix &leftPt=
    *leftBranch.getHMM()->getSubstMatrix(fcParent,fcLeft);
  cout<<"int4.3"<<endl;
  SubstitutionMatrix &rightPt=
    *rightBranch.getHMM()->getSubstMatrix(fcParent,fcRight);
  cout<<"int5"<<endl;
  if(!&leftPt || !&rightPt) throw "you must resolve all functional classes before computing the likelihood"; // ### DEBUGGING
  for(Symbol a=0 ; a<numAlpha ; ++a)
    row[a]=
      (a==gap) ? log0
      : safeAdd(processInternalChild(a,leftChild,leftPt,leftBranch,
				     fcParent,fcLeft),
		processInternalChild(a,rightChild,rightPt,rightBranch,
				     fcParent,fcLeft));
  cout<<"int6"<<endl;
}



int LinkFelsenstein::recurseToChild(BranchAttributes &branch,int parentIndex)
{
  cout<<"rc1"<<endl;
  const IndexMap &indexMap=branch.getDownMap();
  int childIndex=indexMap[parentIndex];
  Taxon &child=branch.getChildTaxon();
  cout<<"rc2"<<endl;
  if(childIndex==IndexMap::UNDEFINED) {
    Array2D<double>::RowIn2DArray<double> row=L[child.getID()];
    for(Symbol a=0 ; a<numAlpha ; ++a) 
      row[a]=(a==gap ? log0 : log1);
  }
  else postorder(child,childIndex);
  cout<<"rc3"<<endl;
  return childIndex;
}



void LinkFelsenstein::postorderRoot(Taxon &taxon,int residueIndex)
{
  cout<<"postorderRoot"<<endl;

  // First, visit the child subtree
  BranchAttributes *branch=taxon.getBranch(UNIQUE);
  int childIndex=recurseToChild(*branch,residueIndex);

  // Now process this node
  int id=taxon.getID();
  Array2D<double>::RowIn2DArray<double> row=L[id];
  Taxon &child=taxon.getIthChild(0);
  int childID=child.getID();
  FunctionalClass fcParent=taxon.getFunctionalParse()[residueIndex];
  FunctionalClass fcChild=
    (childIndex==IndexMap::UNDEFINED ? fcParent
     : child.getFunctionalParse()[childIndex]);
  BranchHMM *hmm=branch->getHMM();
  SubstitutionMatrix &Pt=*hmm->getSubstMatrix(fcParent,fcChild);
  if(!&Pt) INTERNAL_ERROR;
  for(Symbol a=0 ; a<numAlpha ; ++a)
    row[a]=(a==gap) ? log0 
      : processInternalChild(a,child,Pt,*branch,fcParent,fcChild);
}



double LinkFelsenstein::processInternalChild(Symbol parentSymbol,
					     Taxon &child,
					     SubstitutionMatrix &Pt,
					     BranchAttributes &branch,
					     FunctionalClass parentFC,
					     FunctionalClass childFC)
{
  BranchHMM *hmm=branch.getHMM();
  GainLossType gainLossType=
    FunctionalClass::classifyGainLoss(parentFC,childFC);
  double gainLossFactor=
    useGainAndLossScores && focalBranch1!=&branch && focalBranch2!=&branch
    ? hmm->getGainLossFactor(gainLossType) : 0;
  Array2D<double>::RowIn2DArray<double> row=L[child.getID()];
  Array1D<double> V(numAlpha);
  for(Symbol b=0 ; b<numAlpha ; ++b)
    V[b]=(b==gap) ? log0 
      : safeAdd(Pt(parentSymbol,b),row[b],gainLossFactor);
  return sumLogProbs<double>(V);
}



void LinkFelsenstein::setFocalBranches(BranchAttributes *b1,
					BranchAttributes *b2)
{
  focalBranch1=b1;
  focalBranch2=b2;
}



const Alphabet &LinkFelsenstein::getAlphabet()
{
  return alphabet;
}



void LinkFelsenstein::precomputeInside(Array1D<float> &dp,Taxon &taxon,
				       int pos,FunctionalClass fc)
{
  postorder(taxon,pos);
  Array2D<double>::RowIn2DArray<double> row=L[taxon.getID()];
  for(Symbol i=0 ; i<numAlpha ; ++i)
    dp[i]=row[i];
}



void LinkFelsenstein::precomputeOutside(Array1D<float> &dp,Taxon &taxon,
					int pos,Taxon &avoidChild,
					FunctionalClass fc)
{
  outorder(taxon,pos,avoidChild);
  Array2D<double>::RowIn2DArray<double> row=L[taxon.getID()];
  for(Symbol i=0 ; i<numAlpha ; ++i)
    dp[i]=row[i];
}



void LinkFelsenstein::outorder(Taxon &taxon,int residueIndex,Taxon &avoid)
{
  switch(taxon.getNodeType()) 
    {
    case ROOT_NODE:     
      // not currently implemented (not needed by AVES)
      INTERNAL_ERROR;
      break;
    case LEAF_NODE:     
      postorderLeaf(taxon,residueIndex); 
      break;
    case INTERNAL_NODE: 
      {
	ResidueAddress ra(&taxon,residueIndex);
	ResidueAddress parent=ra.getParent();
	if(parent.isValid()) 
	  outorder(*parent.getTaxon(),parent.getIndex(),taxon);
	WhichChild oc=::otherChild(avoid.whichChild());
	ResidueAddress otherChild=ra.getChild(oc);
	
	//########
	if(otherChild.isValid())
	  postorder(*otherChild.getTaxon(),otherChild.getIndex());
	outorderInternal(taxon,residueIndex,avoid); 
      }
      break;
    }
}




