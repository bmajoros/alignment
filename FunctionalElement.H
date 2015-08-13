/****************************************************************
 FunctionalElement.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_FunctionalElement_H
#define INCL_FunctionalElement_H
#include <iostream>
#include "EGGS/Strand.H"
#include "BOOM/String.H"
#include "BOOM/Vector.H"
#include "BOOM/Map.H"
using namespace std;
using namespace BOOM;

class FunctionalClass;
class Taxon;

/****************************************************************
                   class FunctionalElementType
 ****************************************************************/
class FunctionalElementType {
public:
  FunctionalElementType();
  FunctionalElementType(const FunctionalElementType &);
  const String &getName();
  int getID() const;
  bool isValid() const;
  bool operator==(const FunctionalElementType &) const;
  bool operator!=(const FunctionalElementType &) const;
  static const FunctionalElementType NO_FUNC_ELEM;
  static FunctionalElementType registerType(const String &name);
  static int numTypes();
  static FunctionalElementType getTypeByName(const String &name);
  Vector<FunctionalClass> &getClasses();
  void appendClass(FunctionalClass);
  static void resetAll();
  void printOn(ostream &os) {os<<ID;}
  Strand getStrand(); // because some element types are on "the other strand"
  int getLength() { return getClasses().size(); }
protected:
  FunctionalElementType(int ID);
  short int ID;
  static BOOM::Vector<String> names;
  static BOOM::Map<String,int> nameToID;
  static BOOM::Map<int,Vector<FunctionalClass> > fet2fc;
};

inline ostream &operator<<(ostream &os,const FunctionalElementType &fet)
  { fet.printOn(os); return os; }


/****************************************************************
                       class FunctionalElement
 ****************************************************************/
class FunctionalElement {
public:
  FunctionalElement();
  FunctionalElement(FunctionalElementType,int begin,int end,Strand,
		    Taxon * =NULL);
  int getBegin() const;
  int getEnd() const;
  Strand getStrand() const;
  FunctionalElementType getType() const;
  Taxon *getTaxon() {return taxon;}
  void setTaxon(Taxon *t) {taxon=t;}
  FunctionalElement *getParent() {return parent;}
  void setParent(FunctionalElement *p) {parent=p;}
protected:
  int begin, end; // begin is inclusive, end is not
  Strand strand;
  FunctionalElementType type;
  Taxon *taxon;
  FunctionalElement *parent;
};

#endif

