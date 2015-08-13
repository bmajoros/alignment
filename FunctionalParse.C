/****************************************************************
 FunctionalParse.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "FunctionalParse.H"
#include "BranchHMM.H"
#include "Taxon.H"
using namespace std;
using namespace BOOM;

FunctionalParse::FunctionalParse()
  : taxon(NULL)
{
  // ctor
}




FunctionalParse::FunctionalParse(const StatePath &path,BranchEnd branchEnd)
  : taxon(NULL)
{
  PairHMM *hmm=path.getHMM();
  if(!hmm) throw "NULL hmm in StatePath";
  BranchHMM *bhmm=dynamic_cast<BranchHMM*>(hmm);
  if(!bhmm) throw "Attempt to get functional parse from non-BranchHMM";
  int L=path.length();
  for(int i=0 ; i<L ; ++i) {
    STATE s=path[i];
    bool keep;
    switch(hmm->getStateType(s))
      {
      case PHMM_MATCH:
	keep=true;
	break;
      case PHMM_INSERT:
	keep=(branchEnd==CHILD);
	break;
      case PHMM_DELETE:
	keep=(branchEnd==PARENT);
	break;
      default:
	keep=false;
	break;
      }
    if(keep) parse.push_back(bhmm->getFunctionalClass(s,branchEnd));
  }
}



int FunctionalParse::length() const
{
  return size();
}



int FunctionalParse::size() const
{
  return parse.size();
}



void FunctionalParse::push_back(FunctionalClass c)
{
  parse.push_back(c);
}



void FunctionalParse::clear()
{
  parse.clear();
}



FunctionalClass FunctionalParse::operator[](int i) const
{
  if(i>=parse.size()) 
    throw String("FunctionalParse::operator[] : out of bounds: ")+i;
  return parse[i];
}



FunctionalClass &FunctionalParse::operator[](int i)
{
  if(i>=parse.size()) 
    throw String("FunctionalParse::operator[] : out of bounds: ")+i;
  return parse[i];
}



void FunctionalParse::resize(int x)
{
  parse.resize(x);
}



void FunctionalParse::printOn(ostream &os)
{
  int L=size();
  //cout<<"L="<<L<<endl;  for(int i=0 ; i<L ; ++i) os<<i<<" "<<parse[i].getLabel()<<endl;//###
  for(int i=0 ; i<L ; ++i)
    os<<parse[i].getLabel();
}



ostream &operator<<(ostream &os,const FunctionalParse &phi)
{
  phi.printOn(os);
  return os;
}



Array1D<FunctionalElement> *FunctionalParse::getElements()
{
  Vector<FunctionalElement> elements;
  int L=parse.size();
  FunctionalElementType currentType=FunctionalElementType::NO_FUNC_ELEM;
  int begin;
  for(int i=0 ; i<L ; ++i) {
    FunctionalElementType t=parse[i].getElementType();
    if(t!=currentType) {
      if(t.isValid()) begin=i;
      else elements.push_back(FunctionalElement(currentType,begin,i,
						currentType.getStrand(),
						taxon));
      //else elements.push_back(FunctionalElement(currentType,begin,i,
      //					FORWARD_STRAND,taxon));
      currentType=t;
    }
  }
  return new Array1D<FunctionalElement>(elements);
}





