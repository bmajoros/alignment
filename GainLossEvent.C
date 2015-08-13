/****************************************************************
 GainLossEvent.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "GainLossEvent.H"
#include "BOOM/Exceptions.H"
#include "BranchHMM.H"
using namespace std;
using namespace BOOM;

GainLossEvent::GainLossEvent()
  : type(GLT_VOID)
{
  // ctor
}

GainLossEvent::GainLossEvent(GainLossType t,FunctionalElement parentElem,
			     FunctionalElement childElem)
  : type(t), parentElem(parentElem), childElem(childElem)
{
  // ctor
}



GainLossType GainLossEvent::getType() const
{
  return type;
}



FunctionalElement GainLossEvent::getParentElem() const
{
  return parentElem;
}



FunctionalElement GainLossEvent::getChildElem() const
{
  return childElem;
}



FunctionalElementType GainLossEvent::getElementType() const
{
  switch(type) 
    {
    case GLT_GAIN:
    case GLT_RETENTION:
      return childElem.getType();
    case GLT_LOSS:
      return parentElem.getType();
    case GLT_VOID:
      return FunctionalElementType::NO_FUNC_ELEM;
    }
  INTERNAL_ERROR;
}



Taxon &GainLossEvent::getParent() const
{
  return *parentElem.getTaxon();
}



Taxon &GainLossEvent::getChild() const
{
  return *childElem.getTaxon();
}



int GainLossEvent::getParentBegin() const
{
  return parentElem.getBegin();
}



int GainLossEvent::getChildBegin() const
{
  return childElem.getBegin();
}



int GainLossEvent::getParentEnd() const
{
  return parentElem.getEnd();
}



int GainLossEvent::getChildEnd() const
{
  return childElem.getEnd();
}



Vector<GainLossEvent> *GainLossEvent::getEvents(StatePath &path,Taxon &parent,
						Taxon &child)
{
  const BranchHMM &hmm=*static_cast<const BranchHMM*>(path.getHMM());
  Vector<GainLossEvent> *events=new Vector<GainLossEvent>;
  int L=path.length(), parentPos=0, childPos=0;
  FunctionalClass prevParentC, prevChildC;
  for(int i=0 ; i<L ; ++i) {
    STATE q=path[i];
    FunctionalClass parentC=hmm.getFunctionalClass(q,PARENT);
    FunctionalClass childC=hmm.getFunctionalClass(q,CHILD);
    Strand parentStrand=parentC.getStrand();
    Strand childStrand=childC.getStrand();
    //if(i>0 && child.getName()=="moj3") cout<<i<<" "<<parentC<<" "<<childC<<" "<<bool(childC!=prevChildC)<<" "<<childC.beginsElement()<<" "<<childC.endsElement()<<endl;
    if(parentC!=prevParentC && parentC.beginsElement() || 
       childC!=prevChildC && childC.beginsElement()) {
      GainLossType gainLoss=FunctionalClass::classifyGainLoss(parentC,childC);
      Strand strand=parentC.isBackground() ? childStrand : parentStrand;
      FunctionalElementType parentType=parentC.getElementType();
      FunctionalElement parentElem(parentType,parentPos,
				   parentPos+parentType.getLength(),
				   strand,&parent);
      FunctionalElementType childType=childC.getElementType();
      FunctionalElement childElem(childType,childPos,
				   childPos+childType.getLength(),
				   strand,&child);
      events->push_back(GainLossEvent(gainLoss,parentElem,childElem));
      //cout<<"PUSHED"<<endl;
    }
    hmm.updateColumnsFwd(q,parentPos,childPos);
    prevParentC=parentC;
    prevChildC=childC;
  }
  return events;
}




