/****************************************************************
 FunctionalDollo.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "FunctionalDollo.H"
using namespace std;
using namespace BOOM;



FunctionalDollo::FunctionalDollo(Phylogeny &phylogeny)
  : phylogeny(phylogeny), background(FunctionalClass::getBackground())
{
  // ctor
}



void FunctionalDollo::run() {
  changes=true;
  while(changes) {
    changes=false;
    Vector<PhylogenyNode*> nodes;
    phylogeny.gatherNodes(nodes);
    Vector<PhylogenyNode*>::iterator cur=nodes.begin(), end=nodes.end();
    for(; cur!=end ; ++cur) {
      PhylogenyNode *node=*cur;
      Taxon &taxon=static_cast<Taxon&>(*node->getDecoration());
      int L=taxon.getSeqLen();
      //cout<<"processing "<<taxon.getName()<<" L="<<L<<endl;
      for(int i=0 ; i<L ; ++i) {
	ResidueAddress ra(&taxon,i);
	if(ra.isRoot()) {
	  preorder(ra);
	  //cout<<ra.getIndex()<<" ";
	}
      }
      //cout<<endl;
    }
  }
}



void FunctionalDollo::preorder(ResidueAddress ra) {
  Taxon *taxon=ra.getTaxon();
  switch(taxon->getNodeType()) {
  case ROOT_NODE: 
    {
      ResidueAddress child=ra.getChild(0);
      if(!ra.getFunctionalClass().isBackground() &&
	 (!child.isValid() || child.getFunctionalClass().isBackground())) {
	//cout<<"changing "<<ra.getFunctionalClass()<<" to "<<background<<" at "<<ra<<endl;
	ra.getFunctionalClass()=background;
	changes=true;
      }
      if(child.isValid()) preorder(child);
    }
    break;
  case LEAF_NODE: break;
  case INTERNAL_NODE: 
    {
      ResidueAddress left=ra.getLeftChild(), right=ra.getRightChild();
      if(!ra.getFunctionalClass().isBackground()) {
	int numNeighbors=0;
	ResidueAddress parent=ra.getParent();
	if(parent.isValid() && !parent.getFunctionalClass().isBackground())
	  ++numNeighbors;
	if(left.isValid() && !left.getFunctionalClass().isBackground())
	  ++numNeighbors;
	if(right.isValid() && !right.getFunctionalClass().isBackground())
	  ++numNeighbors;
	if(numNeighbors<2) {
	  //cout<<"changing "<<ra.getFunctionalClass()<<" to "<<background<<" at "<<ra<<endl;
	  ra.getFunctionalClass()=background;
	  changes=true;
	}
      }	
      if(left.isValid()) preorder(left);
      if(right.isValid()) preorder(right);
    }
    break;
  }	
}



