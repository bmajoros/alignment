/****************************************************************
 FunctionalElement.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "FunctionalElement.H"
#include "FunctionalClass.H"
#include "Taxon.H"
using namespace std;
using namespace BOOM;

Vector<String> FunctionalElementType::names;
Map<String,int> FunctionalElementType::nameToID;
const FunctionalElementType FunctionalElementType::NO_FUNC_ELEM(-1);
Map<int,Vector<FunctionalClass> > FunctionalElementType::fet2fc;



/****************************************************************
                  FunctionalElementType methods
 ****************************************************************/
FunctionalElementType::FunctionalElementType(int ID)
  : ID(ID)
{
  // ctor
}



FunctionalElementType::FunctionalElementType()
  : ID(-1)
{
  // ctor
}



FunctionalElementType::FunctionalElementType(const FunctionalElementType &t)
  : ID(t.ID)
{
  // ctor
}



void FunctionalElementType::resetAll()
{
  names.clear();
  nameToID.clear();
  fet2fc.clear();
}



Vector<FunctionalClass> &FunctionalElementType::getClasses()
{
  return fet2fc[ID];
}



void FunctionalElementType::appendClass(FunctionalClass fc)
{
  fet2fc[ID].push_back(fc);
}



FunctionalElementType FunctionalElementType::registerType(const String &name)
{
  int ID=numTypes();
  names.push_back(name);
  nameToID[name]=ID;
  return FunctionalElementType(ID);
}



int FunctionalElementType::numTypes()
{
  return names.size();
}



FunctionalElementType FunctionalElementType::getTypeByName(const String &name)
{
  if(nameToID.isDefined(name))
    return FunctionalElementType(nameToID[name]);
  return  NO_FUNC_ELEM;
}



bool FunctionalElementType::isValid() const
{
  return !(*this==NO_FUNC_ELEM);
}



bool FunctionalElementType::operator==(const FunctionalElementType &t) const
{
  return ID==t.ID;
}



bool FunctionalElementType::operator!=(const FunctionalElementType &t) const
{
  return ID!=t.ID;
}



const String &FunctionalElementType::getName()
{
  static String none("<NO_FUNC_ELEM>");
  return ID>=0 ? names[ID] : none;
}



int FunctionalElementType::getID() const
{
  return ID;
}



Strand FunctionalElementType::getStrand() 
{
  return getClasses()[0].getStrand();
}




/****************************************************************
                      FunctionalElement methods
 ****************************************************************/
FunctionalElement::FunctionalElement(FunctionalElementType t,int begin,
				     int end,Strand s,Taxon *taxon)
  : type(t), begin(begin), end(end), strand(s), taxon(taxon),
    parent(NULL)
{
  // ctor
}



FunctionalElement::FunctionalElement() 
  : type(FunctionalElementType::NO_FUNC_ELEM),
    begin(-1), end(-1), strand(NO_STRAND), taxon(NULL), parent(NULL)
{
  // ctor
}



int FunctionalElement::getBegin() const
{
  return begin;
}



int FunctionalElement::getEnd() const
{
  return end;
}



Strand FunctionalElement::getStrand() const
{
  return strand;
}



FunctionalElementType FunctionalElement::getType() const
{
  return type;
}



