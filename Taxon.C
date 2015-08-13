/****************************************************************
 Taxon.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "Taxon.H"
#include "BranchAttributes.H"
#include "ResidueAddress.H"
#include "BOOM/Exceptions.H"
#include "CollapsedOrthologyMatrix.H"
using namespace std;
using namespace BOOM;


Taxon::Taxon()
  : node(NULL), cladeAlignment(NULL), gapPattern(NULL), branchToParent(NULL),
    fcConstraint(FunctionalClass::NO_CLASS), seqLen(0), com(NULL)
{
  // ctor
}



void Taxon::setNode(PhylogenyNode *n)
{
  node=n;
  initBranches();
}



PhylogenyNode *Taxon::getNode()
{
  return node;
}



Sequence &Taxon::getSeq()
{
  return seq;
}



String &Taxon::getSeqStr()
{
  return seqStr;
}


const String &Taxon::getName() const
{
  return node->getName();
}



int Taxon::getID() const
{
  return node->getID();
}


MultSeqAlignment *Taxon::getCladeAlignment()
{
  return cladeAlignment;
}



void Taxon::setCladeAlignment(MultSeqAlignment *a)
{
  cladeAlignment=a;
}



BitSet &Taxon::getCladeMembers()
{
  return cladeMembers;
}



void Taxon::setGapPattern(GapPattern *p)
{
  delete gapPattern;
  gapPattern=p;
}



GapPattern *Taxon::getGapPattern()
{
  return gapPattern;
}



int Taxon::getNumBranches()
{
  return branches.size();
}



BranchAttributes *&Taxon::getBranch(WhichChild b)
{
  return branches[static_cast<int>(b)];
}



BranchAttributes *&Taxon::getIthBranch(int i)
{
  return branches[i];
}



void Taxon::initBranches()
{
  int numChildren;
  switch(node->getNodeType()) 
    {
    case ROOT_NODE:       numChildren=1; break;
    case LEAF_NODE:       numChildren=0; break;
    case INTERNAL_NODE:   numChildren=2; break;
    }
  branches.resize(numChildren);
  branches.setAllTo(NULL);
}



const Array1D<double> &Taxon::getEqFreqs(STATE q) const
{
  // ### DEBUGGING
  //q=FunctionalClass::getBackground(); // ### DEBUGGING
  // ### /DEBUGGING

  Taxon *parent=getParent();
  if(parent) {
    WhichChild whichChild=node->whichChild();
    BranchAttributes *branchAttr=parent->getBranch(whichChild);
    BranchHMM *hmm=branchAttr->getHMM();
    return hmm->getEqFreqs(q);
  }
  return branches[0]->getHMM()->getEqFreqs(q);
}



Taxon *Taxon::getParent() const
{
  PhylogenyNode *parent=node->getParent();
  if(!parent) return NULL;
  return static_cast<Taxon*>(parent->getDecoration());
}



bool Taxon::isLeaf()
{
  return node->getNodeType()==LEAF_NODE;
}



bool Taxon::isRoot()
{
  return node->getNodeType()==ROOT_NODE;
}



Taxon &Taxon::getIthChild(int i)
{
  return *static_cast<Taxon*>(getIthBranch(i)->getBranch()->getChild()->
			      getDecoration());
}



FunctionalParse &Taxon::getFunctionalParse()
{
  return functionalParse;
}



PhylogenyNodeType Taxon::getNodeType() const
{
  return node->getNodeType();
}



BranchAttributes *Taxon::getBranchToParent()
{
  if(branchToParent) return branchToParent;
  Taxon *parent=getParent();
  if(!parent) return NULL;
  int n=parent->getNumBranches();
  for(int i=0 ; i<n ; ++i) {
    BranchAttributes *branch=parent->getIthBranch(i);
    if(&branch->getChildTaxon()==this) 
      return branchToParent=branch; // yes, I mean assignment (=)
  }
  throw "Taxon::getBranchToParent()";
}



/*
  This function traces up the tree from the current taxon, to find the
  most ancestral taxon with an identifiable ortholog to the specified
  residue in this descendent taxon.
 */
Taxon &Taxon::findResidueRoot(int residueIndex,int &rootResidueIndex)
{
  rootResidueIndex=residueIndex;
  Taxon *taxon=this;
  BranchAttributes *branch=getBranchToParent();
  while(branch && branch->getUpMap().size()>0) {
    int mappedIndex=branch->getUpMap()[rootResidueIndex];
    if(mappedIndex==IndexMap::UNDEFINED) break;
    rootResidueIndex=mappedIndex;
    taxon=&branch->getParentTaxon();
    branch=taxon->getBranchToParent();
  }
  return *taxon;
}



BranchAttributes *Taxon::getBranchToChild(Taxon &child)
{
  int n=getNumBranches();
  for(int i=0 ; i<n ; ++i) {
    BranchAttributes *branch=getIthBranch(i);
    if(&branch->getChildTaxon()==&child) 
      return branch;
  }
  throw "Taxon::getBranchToChild() : child not found";//return NULL;
}



void Taxon::setPosetNode(int residueIndex,PosetNode *node)
{
  posetNodes[residueIndex]=node;
}



PosetNode *Taxon::getPosetNode(int residueIndex)
{
  return posetNodes[residueIndex];
}



void Taxon::getNeighborsExcept(Taxon *except,Vector<Taxon*> &v)
{
  v.clear();
  Taxon *parent=getParent();
  if(parent==this) INTERNAL_ERROR;
  if(parent && parent!=except) v.push_back(parent);
  int n=getNumBranches();
  for(int i=0 ; i<n ; ++i) {
    BranchAttributes *branch=getIthBranch(i);
    Taxon *child=&branch->getChildTaxon();
    if(child!=except) v.push_back(child);
  }
}



void Taxon::getNeighbors(Vector<Taxon*> &v)
{
  v.clear();
  Taxon *parent=getParent();
  if(parent) v.push_back(parent);
  int n=getNumBranches();
  for(int i=0 ; i<n ; ++i) {
    BranchAttributes *branch=getIthBranch(i);
    v.push_back(&branch->getChildTaxon());
  }
}



bool Taxon::isFcConstrained() const
{
  return fcConstraint.isValid();
}



FunctionalClass Taxon::getFcConstraint() const
{
  return fcConstraint;
}



void Taxon::removeFcConstraint()
{
  fcConstraint=FunctionalClass::NO_CLASS;
}



void Taxon::setFcConstraint(FunctionalClass fc)
{
  fcConstraint=fc;
}



Array1D<FuncClassSet> &Taxon::getFcConstraints()
{
  return fcConstraints;
}



PrecomputedEmissions &Taxon::getPrecomputedEmissions()
{
  return precomputedEmissions;
}



CollapsedOrthologyMatrix *&Taxon::getCOM()
{
  return com;
}



BranchAttributes *Taxon::getBranchTo(Taxon &other)
{
  if(&other==getParent()) return branchToParent;
  return getBranch(other.whichChild());
}


