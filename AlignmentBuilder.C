/****************************************************************
 AlignmentBuilder.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "AlignmentBuilder.H"
#include "BOOM/TopologicalSort.H"
#include "BOOM/ListQueue.H"
#include "BOOM/Exceptions.H"
#include "BranchAttributes.H"
#include "LinkParsimony.H"
using namespace std;
using namespace BOOM;


AlignmentBuilder::AlignmentBuilder(PhylogenyNode *root,Alphabet &alphabet,
				   Symbol gapSymbol,int numTaxa,
				   bool wantOrphans)
  : alphabet(alphabet), gapSymbol(gapSymbol), root(root), numTaxa(numTaxa),
    wantOrphans(wantOrphans)
{
  // ctor
}



void AlignmentBuilder::buildPoset(Poset<ResidueAddress> &poset)
{
  struct Visitor : public TreeVisitor {
    Poset<ResidueAddress> &poset;
    Set<Taxon*> &taxa;
    PhylogenyNode *root;
    Visitor(Poset<ResidueAddress> &poset,Set<Taxon*> &taxa,
	    PhylogenyNode *root) 
      : poset(poset), taxa(taxa), root(root) {}
    void processNode(InternalNode &v) {process(v);}
    void processNode(LeafNode &v) {process(v);}
    void processNode(RootNode &v) {process(v);}
    void process(PhylogenyNode &v) {
      Taxon *taxon=static_cast<Taxon*>(v.getDecoration());
      taxa.insert(taxon);
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
  } V(poset,taxa,root);
  root->preorderTraversal(V);
}



void AlignmentBuilder::dumpSeqLengths()
{
  struct Visitor : public TreeVisitor {
    void processNode(InternalNode &v) {
      Taxon &taxon=*static_cast<Taxon*>(v.getDecoration());
      cout<<"  LLL "<<taxon.getName()<<"="<<taxon.getSeqLen()<<endl;
    }
    void processNode(LeafNode &v) {
      Taxon &taxon=*static_cast<Taxon*>(v.getDecoration());
      cout<<"  LLL "<<taxon.getName()<<"="<<taxon.getSeqLen()<<endl;
    }
    void processNode(RootNode &v) {
      Taxon &taxon=*static_cast<Taxon*>(v.getDecoration());
      cout<<"  LLL "<<taxon.getName()<<"="<<taxon.getSeqLen()<<endl;
    }
  } V;
  root->preorderTraversal(V);
}



MultSeqAlignment *AlignmentBuilder::buildAlignment(bool includeAncestors,
						   bool useParsimony)
{
  // First, build a dependency graph indicating which connected 
  // components (represented by their residue roots) must come before 
  // others in the alignment:

  //if(useParsimony) {cout<<"L1"<<endl;dumpSeqLengths();}
  Poset<ResidueAddress> poset;
  //cout<<"build poset"<<endl;
  buildPoset(poset);

  //if(useParsimony) {cout<<"L2"<<endl;dumpSeqLengths();}
  //cout<<"sort dep graph"<<endl;
  // Sort the dependency graph to get a total ordering:
  TopologicalSort<ResidueAddress> sorter;
  typedef typename Poset<ResidueAddress>::Vertex<ResidueAddress> Node;
  //cout<<"sorter.sort"<<endl;
  Array1D<Node*> &sorted=*sorter.sort(poset);
  //cout<<"get align len"<<endl;
  int numColumns=getAlignmentLength(sorted);
  int numRoots=sorted.size();

  //if(useParsimony) {cout<<"L3"<<endl;dumpSeqLengths();}
  //cout<<"init empty align "<<numColumns<<Endl;
  // Initialize an empty alignment of the proper size
  MultSeqAlignment *msa=new MultSeqAlignment(alphabet,gapSymbol);
  Set<Taxon*>::iterator cur=taxa.begin(), end=taxa.end();
  for(; cur!=end ; ++cur) {
    Taxon *taxon=*cur;
    if(taxon->isLeaf() || includeAncestors) {
      AlignmentSeq &track=msa->findOrCreateTrack(taxon->getName());
      String &annoTrack=track.getAnnoTrack();
      annoTrack.padOrTruncate(numColumns);
    }
  }
  //if(useParsimony) {cout<<"L4"<<endl;dumpSeqLengths();}
  msa->extendToLength(numColumns);
  //if(useParsimony) {cout<<"L5"<<endl;dumpSeqLengths();}
  //cout<<"iterate total ord"<<Endl;
  //cout<<"BEFORE"<<endl; dumpSeqLengths();
  // Iterate across the total ordering, appending a new alignment column
  // for each element:
  if(useParsimony) {
    struct Visitor : public TreeVisitor {
      Symbol gapSymbol;
      Visitor(Symbol gap) : gapSymbol(gap) {}
      void processNode(InternalNode &v) {
	Taxon &taxon=*static_cast<Taxon*>(v.getDecoration());
	//cout<<"extending "<<taxon.getName()<<" "<<taxon.getSeqLen()<<endl;
	taxon.getSeq()=Sequence(gapSymbol,taxon.getSeqLen());
      }
      //void processNode(LeafNode &v) {}
      //void processNode(RootNode &v) {}
    } V(gapSymbol);
    root->preorderTraversal(V);
  }
  //cout<<"AFTER"<<endl;  dumpSeqLengths();
  //cout<<"last loop"<<endl;
  int col=0;
  for(int i=0 ; i<numRoots /*numColumns*/ ; ++i) {
    ResidueAddress addr=sorted[i]->getData();
    //cout<<"root "<<addr<<endl;
    if(!addr.hasLeafDescendents() && !wantOrphans) continue; // ### 
    //cout<<i<<" "<<useParsimony<<" "<<alphabet.size()<<endl;
    if(useParsimony) parsimony(addr);
    addColumn(col++,*msa,addr,includeAncestors);//NOT THE PROBLEM
    //cout<<"INDEX "<<i<<endl;
    //dumpSeqLengths();
  }
  //if(useParsimony) INTERNAL_ERROR; // DEBUGGING
  //cout<<"done last loop"<<endl;
  // Clean up
  delete &sorted;
  return msa;
}



int AlignmentBuilder::getAlignmentLength(
  Array1D<Poset<ResidueAddress>::Vertex<ResidueAddress>*> &roots)
{
  int n=roots.size(), r=0;
  for(int i=0 ; i<n ; ++i)
    if(roots[i]->getData().hasLeafDescendents() || wantOrphans) ++r;
  return r;
}



void AlignmentBuilder::addColumn(int col,MultSeqAlignment &msa,
				 const ResidueAddress &root,
				 bool includeAncestors)
{
  ListQueue<ResidueAddress> Q;
  Q.enqueue(root);
  while(!Q.isEmpty()) {
    ResidueAddress addr=Q.dequeue();
    Taxon *taxon=addr.getTaxon();
    int index=addr.getIndex();
    FunctionalParse &parse=taxon->getFunctionalParse();
    bool hasParse=(parse.length()>0);
    if(taxon->isLeaf() || includeAncestors) {
      AlignmentSeq &track=msa.getTrackByName(taxon->getName());
      Sequence &seq=taxon->getSeq();
      if(index<seq.getLength()) {
	track[col]=taxon->getSeq()[index];
	if(hasParse && index>=parse.size()) {
	  cout<<taxon->getName()<<" "<<parse.length()<<" "<<taxon->getSeqLen()<<" "<<seq.getLength()<<" "<<index<<endl;
	  INTERNAL_ERROR;//### DEBUGGING
	}
	track.getAnnoTrack()[col]=hasParse ? parse[index].getLabel() : ' ';
      }
      else {
	track[col]='*';
	track.getAnnoTrack()[col]=hasParse ? parse[index].getLabel() : '?';
      }	
    }
    if(!taxon->isLeaf()) {
      int numBranches=taxon->getNumBranches();
      for(int i=0 ; i<numBranches ; ++i) {
	BranchAttributes *branch=taxon->getIthBranch(i);
	Taxon *child=&branch->getChildTaxon();
	int childIndex=branch->getDownMap()[index];
	if(childIndex!=IndexMap::UNDEFINED)
	  Q.enqueue(ResidueAddress(child,childIndex));
      }
    }
  }
}



void AlignmentBuilder::parsimony(const ResidueAddress &root)
{
  //LinkParsimony P(root,alphabet,gapSymbol,numTaxa);
  Symbol s=alphabet.lookup('N');
  if(s==INVALID_SYMBOL) s=alphabet.size()-1;
  LinkParsimony P(root,alphabet,gapSymbol,numTaxa,s);
  P.run();
}



