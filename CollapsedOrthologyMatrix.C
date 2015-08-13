/****************************************************************
 CollapsedOrthologyMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "CollapsedOrthologyMatrix.H"
#include "BranchAttributes.H"
using namespace std;
using namespace BOOM;

static const int CollapsedOrthologyMatrix::UNDEFINED = -999;
typedef CollapsedOrthologyMatrix::Entry COM_Entry;

CollapsedOrthologyMatrix::CollapsedOrthologyMatrix(Taxon &taxon,int numTaxa)
  : taxon(taxon), M(taxon.getSeqLen()+1,numTaxa)
{
  // ctor
}



COM_Entry &CollapsedOrthologyMatrix::operator()(int pos,int taxonID)
{
  return M[pos][taxonID];
}



int CollapsedOrthologyMatrix::getNumTaxa()
{
  return M.getSecondDim();
}



int CollapsedOrthologyMatrix::getSeqLen()
{
  return M.getFirstDim()-1;
}



void CollapsedOrthologyMatrix::compute(ResidueOrthologyGraph &rog)
{
  Phylogeny *tree=rog.getTree();
  int numTaxa=rog.getNumTaxa();

  // Up-pass:
  Vector<PhylogenyNode*> nodes;
  tree->gatherNodes(nodes,POSTORDER);
  int numNodes=nodes.size();
  for(int i=0 ; i<numNodes ; ++i) {
    PhylogenyNode *node=nodes[i];
    Taxon &taxon=static_cast<Taxon&>(*node->getDecoration());
    CollapsedOrthologyMatrix *&com=taxon.getCOM();
    delete com;
    com=new CollapsedOrthologyMatrix(taxon,numTaxa);
    if(taxon.isLeaf()) {
      int taxonID=taxon.getID();
      int L=taxon.getSeqLen();
      for(int i=1 ; i<=L ; ++i) {
	CollapsedOrthologyMatrix::Entry &e=(*com)(i,taxonID);
	e.begin=e.end=i;
      }
      (*com)(0,taxonID).set(-1,0);
    }	
    else upPassInternal(taxon,rog);
  }

  // Down-pass:
  tree->gatherNodes(nodes,PREORDER);
  for(int i=0 ; i<numNodes ; ++i) {
    PhylogenyNode *node=nodes[i];
    Taxon &taxon=static_cast<Taxon&>(*node->getDecoration());
    downPass(taxon,rog);
  }
}



void CollapsedOrthologyMatrix::upPassInternal(Taxon &taxon,
					      ResidueOrthologyGraph &rog)
{
  Array1D<Taxon> &taxa=rog.getTaxa();
  int L=taxon.getSeqLen();
  CollapsedOrthologyMatrix &com=*taxon.getCOM();
  int numTaxa=com.getNumTaxa();
  int numBranches=taxon.getNumBranches();
  for(int i=0 ; i<numBranches ; ++i) {
    BranchAttributes &branch=*taxon.getIthBranch(i);
    Taxon &child=branch.getChildTaxon();
    CollapsedOrthologyMatrix &childCom=*child.getCOM();
    IndexMap &downMap=branch.getDownMap();
    BitSet &cladeMembers=child.getCladeMembers();
    Vector<int> leaves;
    for(int i=0 ; i<numTaxa ; ++i)
      if(taxa[i].isLeaf() && cladeMembers.isMember(i))
	leaves.push_back(i);
    int numLeaves=leaves.size();

    // Forward pass:
    for(int j=0 ; j<numLeaves ; ++j) 
      com(0,leaves[j]).begin=-1;
    for(int pos=0 ; pos<L ; ++pos) {
      int childPos=downMap[pos];
      if(childPos!=IndexMap::UNDEFINED)
	for(int j=0 ; j<numLeaves ; ++j) {
	  int leafID=leaves[j];
	  com(pos+1,leafID)=childCom(childPos+1,leafID);
	}	
      else for(int j=0 ; j<numLeaves ; ++j) {
	int leafID=leaves[j];
	com(pos+1,leafID).begin=com(pos,leafID).begin;
      }
      /*
      if(pos==0) 
	for(int j=0 ; j<numLeaves ; ++j)
	  cout<<"XXXf "<<taxon.getName()<<":"<<taxa[leaves[j]].getName()<<" "<<com(pos,leaves[j])<<endl;
      */
    }

    // Backward pass:
    for(int pos=L-1 ; pos>=0 ; --pos) {
      if(downMap[pos]==IndexMap::UNDEFINED)
	if(pos<L-1)
	  for(int j=0 ; j<numLeaves ; ++j) {
	    int leafID=leaves[j];
	    com(pos+1,leafID).end=com(pos+2,leafID).end;
	  }
	else
	  for(int j=0 ; j<numLeaves ; ++j) {
	    int leafID=leaves[j];
	    int leafLen=taxa[leafID].getSeqLen();
	    com(L,leafID).end=leafLen+1;
	  }
      /*
      if(pos==0) 
	for(int j=0 ; j<numLeaves ; ++j)
	  cout<<"XXXb "<<taxon.getName()<<":"<<taxa[leaves[j]].getName()<<" "<<com(pos,leaves[j])<<endl;
      */
    }
    for(int j=0 ; j<numLeaves ; ++j) {
      int leafID=leaves[j];
      com(0,leafID).end=com(1,leafID).end;
    }
  }
}



void CollapsedOrthologyMatrix::downPass(Taxon &taxon,
					ResidueOrthologyGraph &rog)
{
  Array1D<Taxon> &taxa=rog.getTaxa();
  int L=taxon.getSeqLen(), numTaxa=taxa.size();
  CollapsedOrthologyMatrix &com=*taxon.getCOM();
  BranchAttributes *branch=taxon.getBranchToParent();
  if(!branch) return;
  Taxon &parent=branch->getParentTaxon();
  CollapsedOrthologyMatrix &parentCom=*parent.getCOM();
  IndexMap &upMap=branch->getUpMap();
  BitSet &cladeMembers=taxon.getCladeMembers();
  Vector<int> leaves;
  for(int i=0 ; i<numTaxa ; ++i)
    if(taxa[i].isLeaf() && !cladeMembers.isMember(i)) 
      leaves.push_back(i);
  int numLeaves=leaves.size();

  // Forward pass:
  for(int j=0 ; j<numLeaves ; ++j)
    com(0,leaves[j]).begin=-1;
  for(int pos=0 ; pos<L ; ++pos) {
    int parentPos=upMap[pos];
    if(parentPos!=IndexMap::UNDEFINED)
      for(int j=0 ; j<numLeaves ; ++j) {
	int leafID=leaves[j];
	com(pos+1,leafID)=parentCom(parentPos+1,leafID);
      }	
    else for(int j=0 ; j<numLeaves ; ++j) {
      int leafID=leaves[j];
      com(pos+1,leafID).begin=com(pos,leafID).begin;
    }
    /*
    if(pos==0) 
      for(int j=0 ; j<numLeaves ; ++j)
	cout<<"XXXdf "<<taxon.getName()<<":"<<taxa[leaves[j]].getName()<<" "<<com(pos,leaves[j])<<endl;
    */
  }
  
  // Backward pass:
  for(int pos=L-1 ; pos>=0 ; --pos) {
    if(upMap[pos]==IndexMap::UNDEFINED)
      if(pos<L-1)
	for(int j=0 ; j<numLeaves ; ++j) {
	  int leafID=leaves[j];
	  com(pos+1,leafID).end=com(pos+2,leafID).end;
	}
      else for(int j=0 ; j<numLeaves ; ++j) {
	int leafID=leaves[j];
	int leafLen=taxa[leafID].getSeqLen();
	com(L,leafID).end=leafLen+1;
      }
    /*
    if(pos==0) 
      for(int j=0 ; j<numLeaves ; ++j)
	cout<<"XXXdb "<<taxon.getName()<<":"<<taxa[leaves[j]].getName()<<" "<<com(pos,leaves[j])<<endl;
    */
  }
  for(int j=0 ; j<numLeaves ; ++j) {
    int leafID=leaves[j];
    com(0,leafID).end=com(1,leafID).end;
  }
}



CollapsedOrthologyMatrix::Entry 
&CollapsedOrthologyMatrix::Entry::operator+=(const 
					     CollapsedOrthologyMatrix::Entry 
					     &other)
{
  begin=min(begin,other.begin);
  end=max(end,other.end);
  return *this;
}



void CollapsedOrthologyMatrix::Entry::printOn(ostream &os)
{
  os<<"("<<begin<<","<<end<<")";
}



ostream &operator<<(ostream &os,const CollapsedOrthologyMatrix::Entry &e)
{
  e.printOn(os);
  return os;
}


