/****************************************************************
 LinkParsimony.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_LinkParsimony_H
#define INCL_LinkParsimony_H
#include <iostream>
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Array1D.H"
#include "BOOM/BitSet.H"
#include "PhyLib/Phylogeny.H"
#include "ResidueAddress.H"
using namespace std;
using namespace BOOM;


struct LP_Up;
struct LP_Down;


/****************************************************************
 Implements Fitch's maximum parsimony algorithm for ancestral
 state reconstruction.
 ****************************************************************/
class LinkParsimony {
public:
  LinkParsimony(const ResidueAddress &root,Alphabet &,Symbol gapSymbol,
		int numTaxa,Symbol starSymbol=INVALID_SYMBOL);
  void run();
  void runFullSeq();
private:
  ResidueAddress root;
  int numTaxa;
  const Alphabet &alphabet;
  int alphabetSize;
  Array1D<BitSet> initialSets, finalSets;
  Symbol gapSymbol, starSymbol;
  
  friend class LP_Up;
  friend class LP_Down;
};



/****************************************************************
            Utility classes used by LinkParsimony
 ****************************************************************/

class LP_Up : public TreeVisitor {
public:
  LP_Up(LinkParsimony &,const ResidueAddress &,Symbol gap,Symbol starSymbol);
  virtual void processNode(InternalNode &,int index);
  virtual void processNode(LeafNode &,int index);
  virtual void processNode(RootNode &,int index);
private:
  LinkParsimony &master;
  Symbol gapSymbol, starSymbol;
  void recurse(Taxon *,int index);
  void recurseToChild(Taxon &,WhichChild,int parentIndex);
};


class LP_Down : public TreeVisitor {
public:
  LP_Down(LinkParsimony &,const ResidueAddress &,Symbol gap,Symbol star);
  virtual void processNode(InternalNode &,int index,Taxon *parent,
			  int parentIndex);
private:
  LinkParsimony &master;
  Symbol gapSymbol, starSymbol;
  void recurse(Taxon *,int index,Taxon *parent,int parentIndex);
  void recurseToChild(Taxon &,WhichChild,int parentIndex);
};


#endif

