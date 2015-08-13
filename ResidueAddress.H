/****************************************************************
 ResidueAddress.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ResidueAddress_H
#define INCL_ResidueAddress_H
#include <iostream>
#include "Taxon.H"
#include "BranchAttributes.H"
#include "BOOM/IndexMap.H"
#include "PhyLib/Phylogeny.H"
#include "FunctionalClass.H"
using namespace std;
using namespace BOOM;


class ResidueAddress
{
  Taxon *taxon;
  int index;
public:
  inline ResidueAddress();
  inline ResidueAddress(Taxon *,int);
  inline bool isValid() const;
  inline Taxon *getTaxon() const;
  inline int getIndex() const;
  inline FunctionalClass &getFunctionalClass(); // call only if isValid()==true
  inline ResidueAddress findRoot() const;
  inline bool isRoot() const;
  inline bool hasLeafDescendents() const;
  ResidueAddress getParent() const;
  ResidueAddress getLeftChild() const;
  ResidueAddress getRightChild() const;
  inline ResidueAddress getChild(WhichChild) const;
  inline void resetIndex(); // set index to 0
  inline void operator++(); // moves to next residue in sequence
  void printOn(ostream &) const;
};

ostream &operator<<(ostream &,const ResidueAddress &);


ResidueAddress::ResidueAddress(Taxon *t,int i) 
  : taxon(t), index(i) 
{
  // ctor
}



ResidueAddress::ResidueAddress() 
  : taxon(NULL), index(-1) 
{
  // default ctor
}



void ResidueAddress::operator++()
{
  ++index;
}



void ResidueAddress::resetIndex()
{
  index=0;
}



bool ResidueAddress::isValid() const
{
  return taxon!=NULL && index<taxon->getSeqLen() && index>=0;
}



Taxon *ResidueAddress::getTaxon() const
{
  return taxon;
}



int ResidueAddress::getIndex() const
{
  return index;
}


ResidueAddress ResidueAddress::findRoot() const
{
  int rootResidueIndex;
  Taxon &root=taxon->findResidueRoot(index,rootResidueIndex);
  return ResidueAddress(&root,rootResidueIndex);
}



bool ResidueAddress::hasLeafDescendents() const
{
  if(taxon->isLeaf()) return true;
  int n=taxon->getNumBranches();
  for(int i=0 ; i<n ; ++i) {
    BranchAttributes *branch=taxon->getIthBranch(i);
    int childIndex=branch->getDownMap()[index];
    if(childIndex!=IndexMap::UNDEFINED) {
      Taxon *child=&branch->getChildTaxon();
      if(ResidueAddress(child,childIndex).hasLeafDescendents()) 
	return true;
    }
  }
  return false;
}



ResidueAddress ResidueAddress::getChild(WhichChild c) const
{
  return c==LEFT ? getLeftChild() : getRightChild();
}



bool ResidueAddress::isRoot() const
{
  BranchAttributes *branch=taxon->getBranchToParent();
  if(!branch) return true;
  return branch->getUpMap()[index]==IndexMap::UNDEFINED;
}



FunctionalClass &ResidueAddress::getFunctionalClass()
{
  return taxon->getFunctionalParse()[index];
}




#endif

