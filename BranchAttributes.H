/****************************************************************
 BranchAttributes.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BranchAttributes_H
#define INCL_BranchAttributes_H
#include <iostream>
#include "BOOM/IndexMap.H"
#include "BranchHMM.H"
#include "PhyLib/Phylogeny.H"
#include "PairHMM/StatePath.H"
#include "Taxon.H"
using namespace std;
using namespace BOOM;


class BranchAttributes : public BranchDecoration {
public:
  BranchAttributes();
  BranchAttributes(BranchHMM *,PhylogenyBranch *);
  inline BranchHMM *getHMM();
  void changeHMM(BranchHMM *);
  PhylogenyBranch *getBranch();
  double getLength() const;
  int getParentID() const;
  int getChildID() const;
  inline IndexMap &getUpMap(); // maps child indices to parent indices
  inline IndexMap &getDownMap(); // maps parent indices to child indices
  ResidueAddress mapUp(int childIndex);
  ResidueAddress mapDown(int parentIndex);
  StatePath *getStatePath();
  Vector<PHMM_StateType> *inferStatePath(); // infers path from the IndexMaps
  void setStatePath(StatePath *); // also updates the upMap & downMap
  inline Taxon &getParentTaxon();
  inline Taxon &getChildTaxon();
private:
  BranchHMM *hmm;
  PhylogenyBranch *branch;
  IndexMap upMap; // maps child indices to parent indices
  IndexMap downMap; // maps parent indices to child indices
  StatePath *statePath;
};


BranchHMM *BranchAttributes::getHMM()
{
  return hmm;
}



IndexMap &BranchAttributes::getUpMap()
{
  return upMap;
}



IndexMap &BranchAttributes::getDownMap()
{
  return downMap;
}



Taxon &BranchAttributes::getParentTaxon()
{
  return *static_cast<Taxon*>(getBranch()->getParent()->getDecoration());
}



Taxon &BranchAttributes::getChildTaxon()
{
  return *static_cast<Taxon*>(getBranch()->getChild()->getDecoration());
}





#endif

