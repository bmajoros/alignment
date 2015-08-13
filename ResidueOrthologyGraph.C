/****************************************************************
 ResidueOrthologyGraph.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/IndexMap.H"
#include "ResidueOrthologyGraph.H"
#include "BranchAttributes.H"
using namespace std;
using namespace BOOM;


ResidueOrthologyGraph::ResidueOrthologyGraph(Phylogeny *tree,
					     Array1D<Taxon> &taxa)
  : tree(tree), taxa(taxa), numTaxa(taxa.size())
{
  for(int i=0 ; i<numTaxa ; ++i)
    nameToTaxon[taxa[i].getName()]=i;
  tree->collectBranches(branches);//,PREORDER);
}



void ResidueOrthologyGraph::saveConnections(File &f)
{
  Vector<PhylogenyBranch*>::iterator cur=branches.begin(), end=branches.end();
  for(; cur!=end ; ++cur) {
    PhylogenyBranch *branch=*cur;
    BranchAttributes *attr=
      static_cast<BranchAttributes*>(branch->getDecoration());
    IndexMap &up=attr->getUpMap();
    IndexMap &down=attr->getDownMap();
    up.save(f);
    //cout<<"saving Downmap for branch "<<attr->getParentTaxon().getName()<<":"<<attr->getChildTaxon().getName()<<endl;
    down.save(f);
    //cout<<"done saving branch"<<endl;
  }
}



bool ResidueOrthologyGraph::loadConnections(File &f)
{
  Vector<PhylogenyBranch*>::iterator cur=branches.begin(), end=branches.end();
  for(; cur!=end ; ++cur) {
    PhylogenyBranch *branch=*cur;
    BranchAttributes *attr=
      static_cast<BranchAttributes*>(branch->getDecoration());
    IndexMap &up=attr->getUpMap();
    IndexMap &down=attr->getDownMap();
    if(!up.load(f)) return false;
    down.load(f);
    Taxon &parent=attr->getParentTaxon(), &child=attr->getChildTaxon();
    if(child.isLeaf() && up.size()!=child.getSeqLen() && child.getSeqLen()>0)
      up.resize(child.getSeqLen(),true);
    int parentLen=down.size(), childLen=up.size();
    //cout<<"loaded "<<parent.getName()<<"="<<parentLen<<":"<<parent.getSeqLen()<<" "<<child.getName()<<"="<<childLen<<":"<<child.getSeqLen()<<endl;
    parent.setSeqLen(parentLen);
    child.setSeqLen(childLen);
    parent.getFunctionalParse().resize(parentLen);
    child.getFunctionalParse().resize(childLen);
  }
  return true;
}



Vector<PhylogenyBranch*> &ResidueOrthologyGraph::getBranches()
{
  return branches;
}



int ResidueOrthologyGraph::getNumTaxa()
{
  return numTaxa;
}



Phylogeny *ResidueOrthologyGraph::getTree()
{
  return tree;
}



Taxon &ResidueOrthologyGraph::getTaxon(int taxonID)
{
  return taxa[taxonID];
}



int ResidueOrthologyGraph::getTaxonID(const String &name)
{
  return nameToTaxon[name];
}



void ResidueOrthologyGraph::getAlignment(int taxID1,int taxID2,
					 IndexMap &forward,
					 IndexMap &backward)
{
  Vector<PhylogenyNode*> path;
  Taxon &tax1=taxa[taxID1], &tax2=taxa[taxID2];
  tree->getPath(tax1.getNode(),tax2.getNode(),path);
  //forward.resize(tax1.getSeqLen());
  int pathLen=path.size();
  for(int i=0 ; i<pathLen-1 ; ++i) {
    Taxon *thisTax=static_cast<Taxon*>(path[i]->getDecoration());
    Taxon *nextTax=static_cast<Taxon*>(path[i+1]->getDecoration());
    IndexMap *im;
    if(thisTax->getParent()==nextTax) 
      im=&thisTax->getBranchToParent()->getUpMap();
    else
      im=&nextTax->getBranchToParent()->getDownMap();
    if(i==0) forward=*im;
    else forward.compose(*im,forward);
  }
  forward.invert(tax2.getSeqLen(),backward);
}


