/****************************************************************
 PosetBuilder.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_PosetBuilder_H
#define INCL_PosetBuilder_H
#include <iostream>
#include "BOOM/Poset.H"
#include "PhyLib/Phylogeny.H"
#include "ResidueAddress.H"
using namespace std;
using namespace BOOM;

class PosetBuilder {
public:
  typedef Poset<ResidueAddress>::Vertex<ResidueAddress> Node;
  static Poset<ResidueAddress> *buildPoset(PhylogenyNode *root);
};

#endif

