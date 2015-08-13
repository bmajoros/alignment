/****************************************************************
 AlignmentBuilder.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_AlignmentBuilder_H
#define INCL_AlignmentBuilder_H
#include <iostream>
#include "BOOM/Array1D.H"
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Poset.H"
#include "BOOM/Alphabet.H"
#include "BOOM/Set.H"
#include "PhyLib/Phylogeny.H"
#include "Taxon.H"
#include "ResidueAddress.H"
using namespace std;
using namespace BOOM;


class AlignmentBuilder {
  typedef Poset<ResidueAddress>::Vertex<ResidueAddress> Node;
  PhylogenyNode *root;
  BOOM::Symbol gapSymbol;
  Alphabet &alphabet;
  Set<Taxon*> taxa;
  int numTaxa; // total, not how many we've seen here
  bool wantOrphans;
  void buildPoset(Poset<ResidueAddress> &);
  void addColumn(int index,MultSeqAlignment &,const ResidueAddress &root,
		 bool includeAncestors);
  int getAlignmentLength(
     Array1D<Poset<ResidueAddress>::Vertex<ResidueAddress>*> &roots);
  void parsimony(const ResidueAddress &root);
  void dumpSeqLengths();
public:
  AlignmentBuilder(PhylogenyNode *cladeRoot,Alphabet &,BOOM::Symbol gap,
		   int numTaxa,bool showOrphans=false);
  MultSeqAlignment *buildAlignment(bool includeAncestors,bool useParsimony);
};

#endif

