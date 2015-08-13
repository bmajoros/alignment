/****************************************************************
 FunctionalClass.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "FunctionalClass.H"
#include "BOOM/Exceptions.H"
using namespace std;
using namespace BOOM;


// static members:
Vector<String>           FunctionalClass::names;
Vector<char>             FunctionalClass::labels;
Map<String,int>          FunctionalClass::nameToID;
Vector<RateMatrix*>      FunctionalClass::matrices;
Vector<Phylogeny*>       FunctionalClass::phylogenies;
const FunctionalClass    FunctionalClass::NO_CLASS(-1);
FunctionalClass          FunctionalClass::bgFC(-1);
//Vector<bool>             FunctionalClass::background;
Vector<ForegroundBackground> FunctionalClass::foregroundBackground;
Array1D<int>             FunctionalClass::labelToClassID(256);
Vector<Array1D<double> > FunctionalClass::eqFreqs;
Vector<FunctionalElementType> FunctionalClass::elementTypes;
Vector<Strand>           FunctionalClass::strands;



void FunctionalClass::resetAll()
{
  names.clear();
  labels.clear();
  nameToID.clear();

  Vector<RateMatrix*>::iterator rmCur=matrices.begin(), rmEnd=matrices.end();
  for(; rmCur!=rmEnd ; ++rmCur) delete *rmCur;
  matrices.clear();

  //Vector<Phylogeny*>::iterator pCur=phylogenies.begin(), 
  //  pEnd=phylogenies.end();
  //for(; pCur!=pEnd ; ++pCur) {cout<<"ra6.1 "<<*pCur<<endl;delete *pCur;}
  phylogenies.clear();

  foregroundBackground.clear();
  eqFreqs.clear();
  elementTypes.clear();
}



FunctionalClass::FunctionalClass()
  : classID(-1)
{
  // ctor
}



FunctionalClass::FunctionalClass(int classID)
  : classID(classID)
{
  // ctor
}



bool FunctionalClass::isValid() const
{
  return classID>=0 && fg_or_bg()!=FB_NEITHER;
}



FunctionalClass::FunctionalClass(const FunctionalClass &other)
  : classID(other.classID)
{
  // copy ctor
}



const String &FunctionalClass::getName() const
{
  return names[classID];
}



char FunctionalClass::getLabel() const
{
  return labels[classID];
}



RateMatrix *FunctionalClass::getMatrix() const
{
  return matrices[classID];
}



Phylogeny *FunctionalClass::getPhylogeny() const
{
  return phylogenies[classID];
}



bool FunctionalClass::isBackground() const
{
  return foregroundBackground[classID]==BACKGROUND;
}



FunctionalClass FunctionalClass::registerClass(const String &className,
					       char label,
					       FunctionalElementType T,
					       ForegroundBackground fb,
					       RateMatrix *Q,
					       Phylogeny *tree)
{
  int id=names.size();
  nameToID[className]=id;
  names.push_back(className);
  elementTypes.push_back(T);
  labels.push_back(label);
  labelToClassID[label]=id;
  matrices.push_back(Q);
  phylogenies.push_back(tree);
  foregroundBackground.push_back(fb);
  strands.push_back(PLUS_STRAND);
  if(fb==BACKGROUND) bgFC=FunctionalClass(id);

  // Get equilibrium frequencies
  Array1D<double> eq;
  if(Q) {
    SubstitutionMatrix *Pt=Q->instantiate(1.0);
    Pt->getEqFreqs(eq);
    delete Pt;
  }
  eqFreqs.push_back(eq);

  return FunctionalClass(id);
}



FunctionalElementType FunctionalClass::getElementType() const
{
  return elementTypes[classID];
}



FunctionalClass FunctionalClass::classFromLabel(char label)
{
  return labelToClassID[label];
}



int FunctionalClass::getClassID() const
{
  return classID;
}



int FunctionalClass::numClasses()
{
  return names.size();
}



FunctionalClass FunctionalClass::getClassByName(const String &className)
{
  if(nameToID.isDefined(className)) 
    return FunctionalClass(nameToID[className]);
  throw className+" undefined in FunctionClass::getClassByName()";
}



ostream &operator<<(ostream &os,const FunctionalClass &c)
{
  if(c.getClassID()>=0) os<<c.getName();
  else os<<"NO_CLASS";
  return os;
}



FunctionalClass FunctionalClass::getBackground()
{
  return bgFC;
}



ForegroundBackground FunctionalClass::fg_or_bg() const
{
  return foregroundBackground[classID];
}



GainLossType FunctionalClass::classifyGainLoss(FunctionalClass parentFC,
						FunctionalClass childFC)
{
  ForegroundBackground childFB=childFC.fg_or_bg();
  switch(parentFC.fg_or_bg()) {
  case FOREGROUND:
    switch(childFB) {
    case FOREGROUND: return GLT_RETENTION;
    case BACKGROUND: return GLT_LOSS;
    }
    break;
  case BACKGROUND:
    switch(childFB) {
    case FOREGROUND: return GLT_GAIN;
    case BACKGROUND: return GLT_VOID;
    }
    break;
  }
  return GLT_VOID; // ### 7/28/09
}



Array1D<double> &FunctionalClass::getEqFreqs()
{
  return eqFreqs[classID];
}



FunctionalClass pickForeground(FunctionalClass a,
			       FunctionalClass b)
{
  if(!a.isBackground()) return a;
  if(!b.isBackground()) return b;
  return FunctionalClass::getBackground();
}



bool FunctionalClass::beginsElement()
{
  if(!isValid()) return false;
  FunctionalElementType fet=getElementType();
  if(!fet.isValid()) return false;
  Vector<FunctionalClass> &classes=fet.getClasses();
  return classes.size()>0 && classes[0]==*this;
}



bool FunctionalClass::endsElement()
{
  if(!isValid()) return false;
  FunctionalElementType fet=getElementType();
  if(!fet.isValid()) return false;
  Vector<FunctionalClass> &classes=getElementType().getClasses();
  int n=classes.size();
  return n>0 && classes[n-1]==*this;
}



