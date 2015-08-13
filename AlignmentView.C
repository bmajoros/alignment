/****************************************************************
 AlignmentView.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "AlignmentView.H"
using namespace std;
using namespace BOOM;


AlignmentView::AlignmentView(const MultSeqAlignment &A,
			     const BitSet &taxonIDs,
			     Symbol gapSymbol,
			     const Array1D<int> &trackMap)
  : A(A), taxonIDs(taxonIDs), gapSymbol(gapSymbol), trackMap(trackMap)
{
  // ctor

  computeColMap();
}



int AlignmentView::mapColumn(int col) const
{
  return colMap[col];
}



const BitSet &AlignmentView::getTaxonIDs() const
{
  return taxonIDs;
}



void AlignmentView::computeColMap()
{
  int L=A.getLength(), n=taxonIDs.getMaxSize();//A.getNumTracks();
  for(int i=0 ; i<L ; ++i)
    for(int j=0 ; j<n ; ++j) {
      if(taxonIDs.isMember(j)) {
	int targetTrack=trackMap[j];
	if(targetTrack>=0 && A.getIthTrack(targetTrack)[i]!=gapSymbol) {
	  colMap.push_back(i);
	  break;
	}
      }
    }
}



const MultSeqAlignment &AlignmentView::getAlignment() const
{
  return A;
}



int AlignmentView::getLength() const
{
  return colMap.size();
}


