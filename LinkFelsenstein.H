/****************************************************************
 LinkFelsenstein.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_LinkFelsenstein_H
#define INCL_LinkFelsenstein_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Alphabet.H"
#include "BOOM/Vector.H"
#include "PhyLib/Phylogeny.H"
#include "BranchHMM.H"
#include "Taxon.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;


/************************************************************************
  This class implements Felsenstein's algorithm for the "link" framework,
  in which we have a network of residues linked to their ancestral and
  descendent orthologues, rather than a multi-sequence alignment.  The
  connected components of the link network take the place of the MSA 
  columns in the traditional version of Felsenstein's algorithm.

  PRECONDITION: putative functional classes have been installed into the
                appropriate FunctionalParse's of corresponding Taxon's
  PRECONDITION: for deletion states, the DownMap entry of the 
                BranchAttributes has been set to IndexMap::UNDEFINED
***********************************************************************/
class LinkFelsenstein
{
public:
  LinkFelsenstein(const Alphabet &,const AlphabetMap &,int numTaxa);

  // logLikelihood() computes the likelihood of the entire clade rooted
  // at the given taxon:
  virtual double logLikelihood(int residueIndex,Taxon &cladeRoot,
			       bool useGainAndLossScores);
  virtual double logLikelihood(int residueIndex,Taxon &cladeRoot,
			       FunctionalClass,bool useGainAndLossScores);

  // ancestralLikelihood() searches up the tree to find the highest ancestral
  // residue orthologous to the specified residue, and then computes the
  // log likelihood of the clade rooted at that ancestral residue:
  virtual double ancestralLikelihood(int residueIndex,Taxon &,FunctionalClass,
				     bool useGainAndLossScores);
  virtual double ancestralLikelihood(int residueIndex,Taxon &,
				     bool useGainAndLossScores);

  // outsideLikelihood() temporarily prunes the badChild from the parent,
  // then calls ancestralLikelihood(), then restores badChild, so the log
  // likelihood includes everything orthologous to the specified parent
  // residue EXCEPT those contained in the clade rooted at badChild:
  virtual double outsideLikelihood(int parentResidueIndex,Taxon &parent,
				   Taxon &badChild,FunctionalClass,
				   bool useGainAndLossScores);
  virtual double outsideLikelihood(int parentResidueIndex,Taxon &parent,
				   Taxon &badChild,bool useGainAndLossScores);

  // A "focal branch" is one for which gain/loss/retention scores are 
  // NOT to be assessed (since they are aready being assessed via 
  // transition probabilities in Viterbi or Forward/Backward):
  virtual void setFocalBranches(BranchAttributes*,BranchAttributes* =NULL);

  // precomputeInside() fills in a DP matrix row ("dp"), for a given Taxon at
  // a given residue position, in a given functional class, based only on
  // information *inside* the clade rooted at the given taxon:
  virtual void precomputeInside(Array1D<float> &dp,Taxon &,int pos,
				FunctionalClass);
  virtual void precomputeInside(Taxon &) {}

  // precomputeOutside() works like precomputeInside(), except that 
  // information is drawn from everywhere in the tree *except* from within
  // the clade rooted at "badChild":
  virtual void precomputeOutside(Array1D<float> &dp,Taxon &,int pos,
				 Taxon &badChild,FunctionalClass);

  const Alphabet &getAlphabet();
protected:
  const Alphabet &alphabet;
  const AlphabetMap &alphabetMap;
  const BOOM::Symbol gap;
  Array2D<double> L; // indexed as: L[taxon.ID][symbol]
  const int numTaxa, numAlpha;
  bool useGainAndLossScores;
  BranchAttributes *focalBranch1, *focalBranch2; // omit gain/loss scores

  void postorder(Taxon &,int residueIndex);
  virtual void postorderLeaf(Taxon &,int residueIndex);
  virtual void postorderInternal(Taxon &,int residueIndex);
  virtual int recurseToChild(BranchAttributes &branch,int parentIndex);
  virtual void postorderRoot(Taxon &,int residueIndex);
  double processInternalChild(Symbol parentSymbol,Taxon &child,
			      SubstitutionMatrix &,BranchAttributes &,
			      FunctionalClass parentFC,
			      FunctionalClass childFC);
  void outorder(Taxon &taxon,int residueIndex,Taxon &avoid);
  virtual void outorderInternal(Taxon &taxon,int residueIndex,Taxon &avoid) {};
};



#endif

