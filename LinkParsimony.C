/****************************************************************
 LinkParsimony.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "LinkParsimony.H"
#include "BOOM/Random.H"
using namespace std;
using namespace BOOM;



LinkParsimony::LinkParsimony(const ResidueAddress &root,Alphabet &alphabet,
			     Symbol gap,int numTaxa,Symbol star)
  :alphabet(alphabet),
   numTaxa(numTaxa),
   alphabetSize(alphabet.size()),
   root(root),
   gapSymbol(gap),
   starSymbol(star==INVALID_SYMBOL ? gap : star)
{
  // ctor
  
  //if(star==INVALID_SYMBOL) INTERNAL_ERROR; // ### DEBUGGING

  //cout<<"parsimony "<<root<<endl;

  initialSets.resize(numTaxa);
  finalSets.resize(numTaxa);
  for(int i=0 ; i<numTaxa ; ++i) {
    initialSets[i].setSize(alphabetSize);
    finalSets[i].setSize(alphabetSize);
  }
}



void LinkParsimony::run() {
  for(int i=0 ; i<numTaxa ; ++i) {
    initialSets[i].purge();
    finalSets[i].purge();
  }
  LP_Up up(*this,root,gapSymbol,starSymbol);
  LP_Down down(*this,root,gapSymbol,starSymbol);
}



void LinkParsimony::runFullSeq()
{
  int L=root.getTaxon()->getSeqLen();
  root.resetIndex();
  for(int i=0 ; i<L ; ++i) {
    run();
    ++root;
  }
}



/****************************************************************
                          LP_Up methods
 ****************************************************************/
LP_Up::LP_Up(LinkParsimony &master,const ResidueAddress &root,
	     Symbol gapSymbol,Symbol star)
  : master(master),
    gapSymbol(gapSymbol),
    starSymbol(star)
{
  recurse(root.getTaxon(),root.getIndex());
}



void LP_Up::recurse(Taxon *taxon,int index)
{
  switch(taxon->getNodeType())
    {
    case ROOT_NODE:
      processNode(*static_cast<RootNode*>(taxon->getNode()),index);
      break;
    case LEAF_NODE:
      processNode(*static_cast<LeafNode*>(taxon->getNode()),index);
      break;
    case INTERNAL_NODE:
      processNode(*static_cast<InternalNode*>(taxon->getNode()),index);
      break;
    }
}



void LP_Up::recurseToChild(Taxon &taxon,WhichChild whichChild,
			   int parentIndex)
{
  BranchAttributes *branch=taxon.getIthBranch(whichChild);
  int childIndex=branch->getDownMap()[parentIndex];
  //cout<<"  rec "<<branch<<" "<< &branch->getChildTaxon()<<"="<<branch->getChildTaxon().getName()<<" "<<childIndex<<" (was "<<parentIndex<<") "<<branch->getDownMap().size()<<" "<<branch->getChildTaxon().getSeqLen()<<endl;
  if(childIndex!=IndexMap::UNDEFINED) 
    recurse(&branch->getChildTaxon(),childIndex);
}



void LP_Up::processNode(InternalNode &V,int index) 
{
  PhylogenyNode *leftNode=V.getLeft(), *rightNode=V.getRight();
  int ID=V.getID(), leftID=leftNode->getID(), rightID=rightNode->getID();
  BitSet &leftSet=master.initialSets[leftID];
  BitSet &rightSet=master.initialSets[rightID];
  BitSet &thisSet=master.initialSets[ID];
  Taxon &taxon=*static_cast<Taxon*>(V.getDecoration());
  recurseToChild(taxon,LEFT,index);
  recurseToChild(taxon,RIGHT,index);
  leftSet.intersect(rightSet,thisSet);
  if(thisSet.cardinality()==0) leftSet.unionWith(rightSet,thisSet);

  //###
  //if(thisSet.cardinality()==0) thisSet.addMember(starSymbol);
  //###

}



void LP_Up::processNode(LeafNode &V,int index) {
  int id=V.getID();
  Taxon &taxon=*static_cast<Taxon*>(V.getDecoration());
  Symbol s=taxon.getSeq()[index];
  BitSet &bitSet=master.initialSets[id];

  // ### DEBUGGING
  if(s>bitSet.getMaxSize()) { 
    cout<<index<<" "<<taxon.getSeq().getLength()<<" "<<taxon.getSeqLen()<<" "<<int(s)
	<<" "<<taxon.getName()<<" "<<taxon.getSeq()<<endl;
    INTERNAL_ERROR;
  }
  // ### DEBUGGING

  if(s!=gapSymbol) bitSet.addMember(s);
}



void LP_Up::processNode(RootNode &V,int index) {
  int id=V.getID();
  Taxon &taxon=*static_cast<Taxon*>(V.getDecoration());
  Symbol s=taxon.getSeq()[index];
  master.initialSets[id].addMember(s);
  master.finalSets[id].addMember(s);
}



/****************************************************************
                         LP_Down methods
 ****************************************************************/
LP_Down::LP_Down(LinkParsimony &master,const ResidueAddress &root,
		 Symbol gapSymbol,Symbol starSymbol)
  : master(master),
    gapSymbol(gapSymbol),
    starSymbol(starSymbol)
{
  recurse(root.getTaxon(),root.getIndex(),NULL,-1);
}



void LP_Down::recurse(Taxon *taxon,int index,Taxon *parent,int parentIndex)
{
  if(taxon->getNodeType()==INTERNAL_NODE)
    processNode(*static_cast<InternalNode*>(taxon->getNode()),index,
		parent,parentIndex);
}



void LP_Down::recurseToChild(Taxon &taxon,WhichChild whichChild,
			     int parentIndex)
{
  BranchAttributes *branch=taxon.getIthBranch(whichChild);
  int childIndex=branch->getDownMap()[parentIndex];
  if(childIndex!=IndexMap::UNDEFINED) 
    recurse(&branch->getChildTaxon(),childIndex,&taxon,parentIndex);
}



void LP_Down::processNode(InternalNode &V,int index,Taxon *parent,
			  int parentIndex) 
{
  int ID=V.getID();
  BitSet &Fp=master.finalSets[ID];
  BitSet &Sp=master.initialSets[ID];
  if(parent) {
    PhylogenyNode *left=V.getLeft(), *right=V.getRight();
    int parentID=parent->getID(), leftID=left->getID(), rightID=right->getID();
    BitSet &Fa=master.finalSets[parentID];
    Sp.intersect(Fa,Fp);
    if(Fp!=Fa) {
      BitSet &Sq=master.initialSets[leftID];
      BitSet &Sr=master.initialSets[rightID];
      BitSet Sq_I_Sr(master.alphabetSize);
      Sq.intersect(Sr,Sq_I_Sr);
      if(Sq_I_Sr.cardinality()>0) {
	BitSet Sq_U_Sr(master.alphabetSize);
	Sq.unionWith(Sr,Sq_U_Sr);
	BitSet Sq_U_Sr_I_Fa(master.alphabetSize);
	Sq_U_Sr.intersect(Fa,Sq_U_Sr_I_Fa);
	Sq_U_Sr_I_Fa.unionWith(Sp,Fp);
      }
      else {
	Sp.unionWith(Fa,Fp);
      }
    }
  }
  else Fp=Sp;
  
  // Resolve ties and store the inferred symbols into the alignment
  int n=Fp.cardinality();
  if(n>1) {
    int element=0;// randomness is bad for derivatives!
    Symbol s=Fp.getIthMember(element);
    Fp.purge();
    Fp.addMember(s);
  }
  Symbol s=n>0 ? Symbol(Fp.getIthMember(0)) : starSymbol;//gapSymbol;

  Taxon *taxon=static_cast<Taxon*>(V.getDecoration());
  taxon->getSeq()[index]=s;
  recurseToChild(*taxon,LEFT,index);
  recurseToChild(*taxon,RIGHT,index);
}

