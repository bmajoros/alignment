/****************************************************************
 Taxon.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Taxon_H
#define INCL_Taxon_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/BitSet.H"
#include "BOOM/Poset.H"
#include "BOOM/String.H"
#include "PhyLib/Phylogeny.H"
#include "GapPattern.H"
#include "FunctionalParse.H"
#include "PrecomputedEmissions.H"
using namespace std;
using namespace BOOM;

class BranchAttributes;
class ResidueAddress;
class CollapsedOrthologyMatrix;

//typedef short int STATE;
typedef Poset<ResidueAddress>::Vertex<ResidueAddress> PosetNode;

class Taxon : public PhylogenyDecoration {
public:
  Taxon();
  void setNode(PhylogenyNode *);
  PhylogenyNode *getNode();
  PhylogenyNodeType getNodeType() const;
  bool isLeaf();
  bool isRoot();
  Sequence &getSeq();
  String &getSeqStr();
  const String &getName() const;
  int getID() const;
  MultSeqAlignment *getCladeAlignment();
  void setCladeAlignment(MultSeqAlignment *);
  BitSet &getCladeMembers();
  void setGapPattern(GapPattern *);
  GapPattern *getGapPattern();
  int getNumBranches();
  BranchAttributes *&getIthBranch(int i);
  BranchAttributes *&getBranch(WhichChild);
  BranchAttributes *getBranchToChild(Taxon &child);
  BranchAttributes *getBranchToParent();
  BranchAttributes *getBranchTo(Taxon &);
  void removeBranchToParent() {branchToParent=NULL;}
  void getNeighbors(Vector<Taxon*> &);
  void getNeighborsExcept(Taxon *except,Vector<Taxon*> &);
  Taxon &getIthChild(int i);
  const Array1D<double> &getEqFreqs(STATE) const;
  Taxon *getParent() const;
  inline int getSeqLen() const;
  inline void setSeqLen(int);
  FunctionalParse &getFunctionalParse();
  Taxon &findResidueRoot(int residueIndex,int &rootResidueIndex);
  void setPosetNode(int residueIndex,PosetNode *);
  PosetNode *getPosetNode(int residueIndex);
  bool isFcConstrained() const; // for focal species during Felsenstein
  FunctionalClass getFcConstraint() const; // ditto
  Array1D<FuncClassSet> &getFcConstraints(); // for pre-scanned binding sites
  void removeFcConstraint();
  void setFcConstraint(FunctionalClass);
  inline WhichChild whichChild();
  PrecomputedEmissions &getPrecomputedEmissions();
  CollapsedOrthologyMatrix *&getCOM();
  void initBranches();
private:
  PhylogenyNode *node;
  int seqLen;
  String seqStr;
  Sequence seq;
  Array1D<FuncClassSet> fcConstraints;
  Array1D<PosetNode*> posetNodes; // for building explicit alignment to emit
  MultSeqAlignment *cladeAlignment;
  BitSet cladeMembers;
  GapPattern *gapPattern;
  Array1D<BranchAttributes*> branches; // branches to children
  BranchAttributes *branchToParent;
  FunctionalParse functionalParse;
  FunctionalClass fcConstraint; // for current column during progressive phase
  PrecomputedEmissions precomputedEmissions;
  CollapsedOrthologyMatrix *com;
};


typedef pair<Taxon*,Taxon*> TaxonPair;


int Taxon::getSeqLen() const
{
  return seqLen;
}



void Taxon::setSeqLen(int l)
{
  if(l!=seqLen) {
    seqLen=l;
    //if(seqLen==627) INTERNAL_ERROR;
    posetNodes.resize(l);
    fcConstraints.resize(l);
  }
}



WhichChild Taxon::whichChild()
{
  return node->whichChild();
}



#endif

