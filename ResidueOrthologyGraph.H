/****************************************************************
 ResidueOrthologyGraph.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ResidueOrthologyGraph_H
#define INCL_ResidueOrthologyGraph_H
#include <iostream>
#include "BOOM/Array1D.H"
#include "BOOM/File.H"
#include "BOOM/Vector.H"
#include "BOOM/Map.H"
#include "BOOM/IndexMap.H"
#include "PhyLib/Phylogeny.H"
#include "Taxon.H"
using namespace std;
using namespace BOOM;


class ResidueOrthologyGraph {
public:
  ResidueOrthologyGraph(Phylogeny *,Array1D<Taxon> &);
  void saveConnections(File &);
  bool loadConnections(File &); // false on EOF
  Vector<PhylogenyBranch*> &getBranches();
  int getNumTaxa();
  Phylogeny *getTree();
  Taxon &getTaxon(int taxonID);
  Array1D<Taxon> &getTaxa() {return taxa;}
  int getTaxonID(const String &name);
  void getAlignment(int taxID1,int taxID2,IndexMap &forward,
		    IndexMap &backward);
protected:
  Phylogeny *tree;
  Array1D<Taxon> &taxa; // indexed by node ID's from the phylogeny
  Map<String,int> nameToTaxon; // index into "taxa" array
  int numTaxa;
  Vector<PhylogenyBranch*> branches;
};

#endif

