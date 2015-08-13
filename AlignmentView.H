/****************************************************************
 AlignmentView.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_AlignmentView_H
#define INCL_AlignmentView_H
#include <iostream>
#include "BOOM/BitSet.H"
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Vector.H"

using namespace std;
using namespace BOOM;


class AlignmentView {
public:
  AlignmentView(const MultSeqAlignment &,const BitSet &taxonIDs,
		Symbol gapSymbol,const Array1D<int> &trackMap);
  int mapColumn(int) const; // maps virtual column indices to gapless columns
  const BitSet &getTaxonIDs() const;
  const MultSeqAlignment &getAlignment() const;
  int getLength() const;
private:
  const MultSeqAlignment &A;
  const BitSet &taxonIDs;
  const Symbol gapSymbol;
  Vector<int> colMap;
  const Array1D<int> trackMap; // maps taxon ID's to track ID's in A
  void computeColMap();
};

#endif

