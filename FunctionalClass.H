/****************************************************************
 FunctionalClass.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_FunctionalClass_H
#define INCL_FunctionalClass_H
#include <iostream>
#include "BOOM/Vector.H"
#include "BOOM/Map.H"
#include "BOOM/String.H"
#include "BOOM/List.H"
#include "PairHMM/StatePath.H"
#include "PhyLib/RateMatrix.H"
#include "PhyLib/Phylogeny.H"
#include "GainLossType.H"
#include "FunctionalElement.H"
using namespace std;
using namespace BOOM;


enum ForegroundBackground {
  BACKGROUND=0,
  FOREGROUND=1,
  FB_NEITHER=2
};


/****************************************************************
 Instances of this class take only sizeof(int) bytes; static class
 members are used to store the mapping between these int codes and
 the textual descriptions ("exon","UTR",etc...)
 ****************************************************************/
class FunctionalClass
{
  short int classID;
  //static Vector<bool> background;
  static Vector<ForegroundBackground> foregroundBackground;
  static FunctionalClass bgFC;
  static Vector<String> names;
  static Vector<char> labels;
  static Map<String,int> nameToID;
  static Vector<RateMatrix*> matrices;
  static Vector< Array1D<double> > eqFreqs;
  static Vector<Phylogeny*> phylogenies;
  static Array1D<int> labelToClassID;
  static Vector<FunctionalElementType> elementTypes;
  static Vector<Strand> strands;
public:
  FunctionalClass();
  FunctionalClass(int classID);
  FunctionalClass(const FunctionalClass &);
  bool isValid() const;
  static void resetAll();
  static const FunctionalClass NO_CLASS;
  static FunctionalClass registerClass(const String &className,char label,
				       FunctionalElementType,
				       ForegroundBackground,
				       RateMatrix * =NULL,Phylogeny * =NULL);
  static int numClasses();
  static FunctionalClass getClassByName(const String &className);
  Strand getStrand() { return strands[classID]; }
  void setStrand(Strand s) { strands[classID]=s; }
  Array1D<double> &getEqFreqs();
  const String &getName() const;
  char getLabel() const;
  int getClassID() const;
  FunctionalElementType getElementType() const;
  bool isBackground() const;
  ForegroundBackground fg_or_bg() const;
  static FunctionalClass getBackground();
  RateMatrix *getMatrix() const;
  Phylogeny *getPhylogeny() const;
  operator int() const {return classID;}
  inline bool operator==(const FunctionalClass &other);
  inline bool operator!=(const FunctionalClass &other);
  inline bool operator<(const FunctionalClass &other);
  static GainLossType classifyGainLoss(FunctionalClass parentFC,
				       FunctionalClass childFC);
  static FunctionalClass classFromLabel(char label);
  bool beginsElement();
  bool endsElement();
};

FunctionalClass pickForeground(FunctionalClass,FunctionalClass);

ostream &operator<<(ostream &,const FunctionalClass &);


typedef List<FunctionalClass> FuncClassSet;


inline bool FunctionalClass::operator==(const FunctionalClass &other) 
{
  return classID==other.classID;
}



inline bool FunctionalClass::operator!=(const FunctionalClass &other) 
{
  return classID!=other.classID;
}



inline bool FunctionalClass::operator<(const FunctionalClass &other) 
{
  return classID<other.classID;
}



#endif

