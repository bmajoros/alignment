/****************************************************************
 PosetBuilder.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "PosetBuilder.H"
using namespace std;
using namespace BOOM;


Poset<ResidueAddress> *PosetBuilder::buildPoset(PhylogenyNode *root)
{
  Poset<ResidueAddress> *poset=new Poset<ResidueAddress>;
  struct Visitor : public TreeVisitor {
    Poset<ResidueAddress> &poset;
    PhylogenyNode *root;
    Visitor(Poset<ResidueAddress> &poset,PhylogenyNode *root) 
      : poset(poset), root(root) {}
    void processNode(InternalNode &v) {process(v);}
    void processNode(LeafNode &v) {process(v);}
    void processNode(RootNode &v) {process(v);}
    void process(PhylogenyNode &v) {
      Taxon *taxon=static_cast<Taxon*>(v.getDecoration());
      BranchAttributes *parentBranch=taxon->getBranchToParent();
      if(parentBranch && taxon->getNode()!=root) 
	processNonRoot(*taxon,*parentBranch);
      else processRoot(*taxon);
    }
    void processRoot(Taxon &taxon) {
      int L=taxon.getSeqLen();
      Node *prev=NULL;
      for(int i=0 ; i<L ; ++i) {
	Node *node=&poset.addVertex(ResidueAddress(&taxon,i));
	if(!node) INTERNAL_ERROR;
	taxon.setPosetNode(i,node);
	if(prev) prev->addSuccessor(node);
	prev=node;
      }
    }
    void processNonRoot(Taxon &taxon,BranchAttributes &branch) {
      const IndexMap &upMap=branch.getUpMap();
      Taxon *parent=&branch.getParentTaxon();
      int L=taxon.getSeqLen();
      Node *prev=NULL, *node;
      for(int i=0 ; i<L ; ++i) {
	int parentIndex=upMap[i];
	if(parentIndex==IndexMap::UNDEFINED) {
	  node=&poset.addVertex(ResidueAddress(&taxon,i));
	  if(prev) prev->addSuccessor(node);
	}
	else {
	  node=parent->getPosetNode(parentIndex);
	  if(i>0 && upMap[i-1]==IndexMap::UNDEFINED) prev->addSuccessor(node);
	}
	taxon.setPosetNode(i,node);
	prev=node;
      }
    }
  } V(*poset,root);
  root->preorderTraversal(V);
  return poset;
}



