/****************************************************************
 FunctionalDollo.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_FunctionalDollo_H
#define INCL_FunctionalDollo_H
#include <iostream>
#include "PhyLib/Phylogeny.H"
#include "FunctionalClass.H"
#include "ResidueAddress.H"
using namespace std;
using namespace BOOM;


/****************************************************************
 Applies Dollo parsimony to functional parses.
 ****************************************************************/

class FunctionalDollo {
public:
  FunctionalDollo(Phylogeny &);
  void run(); // all columns
private:
  Phylogeny &phylogeny;
  FunctionalClass background;
  bool changes;
  void preorder(ResidueAddress residueRoot);
};


#endif

