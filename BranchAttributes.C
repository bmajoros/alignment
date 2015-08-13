/****************************************************************
 BranchAttributes.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BranchAttributes.H"
#include "ResidueAddress.H"
using namespace std;
using namespace BOOM;


BranchAttributes::BranchAttributes()
  : branch(NULL), hmm(NULL), statePath(NULL)
{
  // default ctor
}



BranchAttributes::BranchAttributes(BranchHMM *hmm,PhylogenyBranch *branch)
  : branch(branch), hmm(hmm), statePath(NULL)
{
  // ctor
}



double BranchAttributes::getLength() const
{
  return branch->getLength();
}



PhylogenyBranch *BranchAttributes::getBranch()
{
  return branch;
}



int BranchAttributes::getParentID() const
{
  return branch->getParent()->getID();
}



int BranchAttributes::getChildID() const
{
  return branch->getChild()->getID();
}



ResidueAddress BranchAttributes::mapUp(int childIndex)
{
  int index=upMap[childIndex];
  if(index==IndexMap::UNDEFINED) return ResidueAddress(NULL,-1);
  return ResidueAddress(&getParentTaxon(),index);
}



ResidueAddress BranchAttributes::mapDown(int parentIndex)
{
  int index=downMap[parentIndex];
  if(index==IndexMap::UNDEFINED) return ResidueAddress(NULL,-1);
  return ResidueAddress(&getChildTaxon(),index);
}



StatePath *BranchAttributes::getStatePath()
{
  return statePath;
}



void BranchAttributes::setStatePath(StatePath *p)
{
  delete statePath;
  statePath=p;
  hmm->decode(*p,downMap,upMap,true);
  //if(!downMap.sanityCheck(getChildTaxon().getSeqLen())) INTERNAL_ERROR;
  //if(!upMap.sanityCheck(getParentTaxon().getSeqLen())) INTERNAL_ERROR;

  /*
  Taxon &parent=getParentTaxon(), &child=getChildTaxon();
  int parentLen=parent.getSeqLen(), childLen=child.getSeqLen();
  cout<<"parentLen="<<parentLen<<endl;
  downMap.resize(parentLen); upMap.resize(childLen);
  int len=p->length();
  cout<<"pathlen="<<len<<endl;
  int i=0, j=0;
  for(int k=0 ; k<len ; ++k) {
    cout<<"k="<<k<<endl;
    STATE q=(*statePath)[k];
    cout<<"q="<<q<<endl;
    cout<<"hmm="<<hmm<<endl;
    cout<<"statetype="<<hmm->getStateType(q)<<endl;
    switch(hmm->getStateType(q)) 
      {
      case PHMM_MATCH:
	downMap[i]=j; upMap[j]=i;
	++i; ++j;
	break;
      case PHMM_INSERT:
	upMap[j]=IndexMap::UNDEFINED;
	++j;
	break;
      case PHMM_DELETE:
	downMap[i]=IndexMap::UNDEFINED;
	++i;
	break;
      default: throw "BranchAtrributes::setStatePath";
      }
  }
  */
}



void BranchAttributes::changeHMM(BranchHMM *h)
{
  //cout<<"this="<<this<<endl;
  //cout<<"hmm="<<hmm<<" h="<<h<<endl;
  delete hmm;
  hmm=h;
  //cout<<"changed"<<endl;
}



Vector<PHMM_StateType> *BranchAttributes::inferStatePath()
{
  Vector<PHMM_StateType> *path=new Vector<PHMM_StateType>;
  int parentLen=downMap.size(), childLen=upMap.size();
  int i=0, j=0;
  //cout<<"parentLen="<<parentLen<<" childLen="<<childLen<<endl;
  for(; i<parentLen || j<childLen ; ) {
    PHMM_StateType t;
    if(i>=parentLen) {
      //cout<<"AAA j="<<j<<endl;
      path->push_back(PHMM_INSERT);
      ++j;
      continue;
    }
    int down=downMap[i];
    //cout<<"BBB i="<<i<<" down="<<down<<endl;
    if(down==j) {t=PHMM_MATCH; ++i; ++j;}
    else if(down==IndexMap::UNDEFINED) {t=PHMM_DELETE; ++i;}
    else {t=PHMM_INSERT; ++j;}
    //cout<<" ===> "<<t<<endl;
    path->push_back(t);
  }
  //cout<<" pathlen="<<path->size()<<endl;
  return path;
}



