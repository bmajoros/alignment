/****************************************************************
 ProfileFelsenstein.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/Array2D.H"
#include "BOOM/Constants.H"
#include "BOOM/SumLogProbs.H"
#include "ProfileFelsenstein.H"
#include "Taxon.H"
#include "BranchAttributes.H"
using namespace std;
using namespace BOOM;


const double log1=0.0;
const double log0=NEGATIVE_INFINITY;


inline double safeAdd(double a,double b)
{
    if(isinf(a) || isinf(b)) return a;
    else return a+b;
}



/****************************************************************
                       class Likelihooder
 ****************************************************************/

struct Likelihooder : public TreeVisitor
{
public:
  Array2D<double> L;
  const MultSeqAlignment &A;
  int numNodes, numAlpha, column;
  const Alphabet &alpha;
  const AlphabetMap &alphabetMap;
  Symbol gap;
  STATE state; // fixed state: assumes "complete orthology"
               // (for progressive phase of aligner ONLY!)
  const Array1D<double> &eqFreqs;
    
  Likelihooder(int numNodes,int numAlpha,const Alphabet &alpha,
	       const AlphabetMap &alphabetMap,Symbol gapSymbol,
	       const MultSeqAlignment &A,int column,STATE state,
	       const Array1D<double> &eqFreqs)
    : A(A), numNodes(numNodes), numAlpha(numAlpha), alpha(alpha), 
      L(numNodes,numAlpha), column(column), gap(gapSymbol), state(state),
      eqFreqs(eqFreqs), alphabetMap(alphabetMap)
  {
    // ctor
  }
  virtual void processNode(LeafNode &u)
  {
    int id=u.getID();
    Symbol a=alphabetMap(A.getIthTrack(id)[column]);
    Array2D<double>::RowIn2DArray<double> row=L[id];
    if(a==gap)
      for(Symbol i=0 ; i<numAlpha; ++i)
	row[i]=log1;// missing data -- same as Seipel & Haussler
    else
      for(Symbol i=0 ; i<numAlpha; ++i) 
	row[i]=(i==a ? log1 : log0);
  }
  double processInternalChild(Symbol parentSymbol,int child,
			      SubstitutionMatrix &Pt)
  {
    Array2D<double>::RowIn2DArray<double> row=L[child];
    Array1D<double> V(numAlpha);
    for(Symbol b=0 ; b<numAlpha ; ++b)
      if(b==gap) V[b]=log0;
      else V[b]=row[b]+Pt(parentSymbol,b);
    return sumLogProbs<double>(V);
  }
  virtual void processNode(InternalNode &u) 
  {
    int id=u.getID();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    Taxon *taxon=static_cast<Taxon*>(u.getDecoration());
    SubstitutionMatrix &leftPt=getMatrix(*taxon,LEFT);
    SubstitutionMatrix &rightPt=getMatrix(*taxon,RIGHT);
    int left=u.getLeft()->getID(), right=u.getRight()->getID();
    for(Symbol a=0 ; a<numAlpha ; ++a)
	if(a!=gap) row[a]=
	  processInternalChild(a,left,leftPt)+
	  processInternalChild(a,right,rightPt);
  }
  SubstitutionMatrix &getMatrix(Taxon &taxon,WhichChild leftRight)
  {
    return *taxon.getIthBranch(leftRight)->getHMM()->getSubstMatrix(state);
  }
  virtual void processNode(RootNode &u)
  {
    // ### technically, the results of this function are never used in
    // the progressive aligner, so this is just for the sake of completness...

    Taxon &taxon=*static_cast<Taxon*>(u.getDecoration());
    int id=u.getID();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    Symbol a=alphabetMap(A.getIthTrack(id)[column]);
    int childID=u.getChild()->getID();
    SubstitutionMatrix &Pt=getMatrix(taxon,LEFT);
    if(a==gap) {
      for(Symbol a=0 ; a<numAlpha ; ++a)
	if(a==gap) continue;
	else row[a]=eqFreqs[a]+processInternalChild(a,childID,Pt);
    }
    else row[a]=eqFreqs[a]+processInternalChild(a,childID,Pt);
  }
  double getLikelihood(PhylogenyNode &node)
  {
    //Taxon &taxon=*static_cast<Taxon*>(node.getDecoration());
    Array1D<double> V(numAlpha);
    int id=node.getID();
    Array2D<double>::RowIn2DArray<double> row=L[id];
    for(Symbol i=0 ; i<numAlpha; ++i) V[i]=log(eqFreqs[i])+row[i];
    return sumLogProbs<double>(V);
  }
};




/****************************************************************
                   ProfileFelsenstein methods
 ****************************************************************/

ProfileFelsenstein::ProfileFelsenstein(const MultSeqAlignment &A,
				       const Alphabet &alphabet,
				       const AlphabetMap &alphabetMap)
  : alignment(A), alphabet(alphabet), numAlpha(alphabet.size()),
    alphabetMap(alphabetMap)
{
  // ctor
}




double ProfileFelsenstein::logLikelihood(int column,PhylogenyNode *node,
					 STATE state)
{
  const Array1D<double> &eqFreqs=
    static_cast<Taxon*>(node->getDecoration())->getEqFreqs(state);
  Likelihooder L(alignment.getNumTracks(),numAlpha,alphabet,alphabetMap,
		 alphabet.lookup('-'),alignment,column,state,
		 eqFreqs);
  node->postorderTraversal(L);
  return L.getLikelihood(*node);
}



