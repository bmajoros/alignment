/****************************************************************
 ResidueAddress.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "ResidueAddress.H"
using namespace std;
using namespace BOOM;



ResidueAddress ResidueAddress::getParent() const
{
  BranchAttributes *branch=taxon->getBranchToParent();
  return branch ? branch->mapUp(index) : ResidueAddress(NULL,-1);
}



ResidueAddress ResidueAddress::getLeftChild() const
{
  BranchAttributes *branch=taxon->getBranch(LEFT);
  return branch ? branch->mapDown(index) : ResidueAddress(NULL,-1);
}



ResidueAddress ResidueAddress::getRightChild() const
{
  BranchAttributes *branch=taxon->getBranch(RIGHT);
  return branch ? branch->mapDown(index) : ResidueAddress(NULL,-1);
}



void ResidueAddress::printOn(ostream &os) const
{
  if(isValid())
    os<<taxon->getName()<<":"<<index;
  else
    os<<"NULL:-1";
}



ostream &operator<<(ostream &os,const ResidueAddress &ra)
{
  ra.printOn(os);
  return os;
}


