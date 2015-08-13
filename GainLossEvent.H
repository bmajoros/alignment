/****************************************************************
 GainLossEvent.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_GainLossEvent_H
#define INCL_GainLossEvent_H
#include <iostream>
#include "BOOM/Vector.H"
#include "PairHMM/StatePath.H"
#include "GainLossType.H"
#include "FunctionalElement.H"
#include "Taxon.H"
using namespace std;
using namespace BOOM;


class GainLossEvent {
public:
  GainLossEvent();
  GainLossEvent(GainLossType,FunctionalElement parentElem,
		FunctionalElement childElem);
  GainLossType getType() const;
  FunctionalElement getParentElem() const;
  FunctionalElement getChildElem() const;
  FunctionalElementType getElementType() const;
  Taxon &getParent() const;
  Taxon &getChild() const;
  int getParentBegin() const;
  int getChildBegin() const;
  int getParentEnd() const;
  int getChildEnd() const;
  static Vector<GainLossEvent> *getEvents(StatePath &,Taxon &parent,
					  Taxon &child);
protected:
  FunctionalElement parentElem, childElem;
  GainLossType type;
};

#endif

