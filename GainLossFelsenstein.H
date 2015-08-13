/****************************************************************
 GainLossFelsenstein.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_GainLossFelsenstein_H
#define INCL_GainLossFelsenstein_H
#include <iostream>
#include "LinkFelsenstein.H"
#include "BOOM/Array3D.H"
#include "GainLossType.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;


/****************************************************************
 This class implements a variant of LinkFelsenstein, in which
 cross-functional states implementing gain or loss events are
 permitted.  Functional classes are not assumed known (except at
 the "root" of the current clade, and at all ancestors of this
 "root" taxon), so we maximize over all possible gain/loss/retention 
 histories along each descendent lineage.
 ****************************************************************/
class GainLossFelsenstein : public LinkFelsenstein
{
public:
  GainLossFelsenstein(const Alphabet &,const AlphabetMap &,int numTaxa);

  // logLikelihood() computes the likelihood of the entire clade rooted
  // at the given taxon, conditional on that "root" taxon being in a 
  // particular functional class.  Descendent taxa may be either in that
  // same class, or in the background FC (whichever is most probable).
  virtual double logLikelihood(int residueIndex,Taxon &cladeRoot,
			       FunctionalClass rootFC,
			       bool useGainAndLossScores);

  // This version of logLikelihood() maximizes over the possible functional
  // classes (foreground/background) of the cladeRoot (note thatforegroundFC 
  // might actually be NO_CLASS, and Taxon::fcConstraint may be set for any
  // taxon in the clade):
  virtual double logLikelihood(int residueIndex,Taxon &cladeRoot,
			       bool useGainAndLossScores,FunctionalClass
			       foregroundFC);

  // ancestralLikelihood() searches up the tree to find the highest ancestral
  // residue orthologous to the specified residue, and then computes the
  // log likelihood of the clade rooted at that ancestral residue.  Note
  // that the specified FunctionalClass is for the given taxon, NOT for the
  // ancestor!  For any node ancestral to the given taxon, the FC is assumed
  // already set in that ancestral taxon's functional parse.  For all others,
  // we choose the FC that maximizes the total probability.
  virtual double ancestralLikelihood(int residueIndex,Taxon &,
				     FunctionalClass,
				     bool useGainAndLossScores);

  // outsideLikelihood() temporarily prunes the badChild from the parent,
  // then calls ancestralLikelihood(), then restores badChild, so the log
  // likelihood includes everything orthologous to the specified parent
  // residue EXCEPT those contained in the clade rooted at badChild:
  virtual double outsideLikelihood(int parentResidueIndex,Taxon &parent,
				   Taxon &badChild,FunctionalClass parentFC,
				   bool useGainAndLossScores);

  // precomputeInside() fills in a DP matrix row ("dp"), for a given Taxon at
  // a given residue position, in a given functional class, based only on
  // information *inside* the clade rooted at the given taxon:
  virtual void precomputeInside(Array1D<float> &dp,Taxon &,int pos,
				FunctionalClass);

  // precomputeOutside() works like precomputeInside(), except that 
  // information is drawn from everywhere in the tree *except* from within
  // the clade rooted at "badChild":
  virtual void precomputeOutside(Array1D<float> &dp,Taxon &,int pos,
				 Taxon &badChild,FunctionalClass);

private:
  Array3D<double> L; // indexed as: L[taxonID][background/foreground][symbol]
  FunctionalClass fgFC, bgFC;

  virtual void postorderLeaf(Taxon &,int residueIndex);
  virtual void postorderInternal(Taxon &,int residueIndex);
  virtual int recurseToChild(BranchAttributes &branch,int parentIndex);
  virtual void postorderRoot(Taxon &,int residueIndex);
  double processInternalChild(BOOM::Symbol parentSymbol,BranchAttributes &,
			      FunctionalClass parentFC,
			      FunctionalClass childFC);
  virtual void outorderInternal(Taxon &taxon,int residueIndex,Taxon &avoid);
  double processParent(BOOM::Symbol childSymbol,BranchAttributes &,
		       FunctionalClass parentFC,FunctionalClass childFC);
};



#endif

