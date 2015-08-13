/****************************************************************
 StatePath.C
 william.majoros@duke.edu

 This is open-source software, governed by the ARTISTIC LICENSE 
 (see www.opensource.org).
 ****************************************************************/
#include <fstream>
#include <iostream>
#include "StatePath.H"
#include "PairHMM.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


StatePath::StatePath(const PairHMM *hmm,int length)
  : path(length), hmm(hmm), score(NEGATIVE_INFINITY)
{
  // ctor
}



void StatePath::push_back(STATE q)
{
  path.push_back(q);
}



int StatePath::size() const
{
  return path.size();
}



int StatePath::length()
{
  return path.size();
}



STATE StatePath::operator[](int i) const
{
  return path[i];
}



STATE &StatePath::operator[](int i)
{
  return path[i];
}



StatePath *StatePath::getReverse() const
{
  StatePath *revPath=new StatePath(hmm);
  int n=size();
  for(int i=n-1 ; i>=0 ; --i) revPath->push_back(path[i]);
  return revPath;
}



void StatePath::setHMM(const PairHMM *h)
{
  hmm=h;
}



const PairHMM *StatePath::getHMM() const
{
  return hmm;
}



void StatePath::printOn(ostream &os) const
{
  int n=path.size();
  for(int i=0 ; i<n ; ++i) {
    os<<path[i];
    switch(hmm->getStateType(path[i])) {
    case PHMM_MATCH:       os<<'m'; break;
    case PHMM_INSERT:      os<<'i'; break;
    case PHMM_DELETE:      os<<'d'; break;
    case PHMM_START_STOP:  os<<'$'; break;
    case PHMM_SILENT:      os<<'_'; break;
    default:               os<<'?'; break;
    }
    os<<' ';
  }
}



ostream &operator<<(ostream &os,const StatePath &path)
{
  path.printOn(os);
  return os;
}



void StatePath::getSeqLengths(int &parentL,int &childL) const
{
  parentL=childL=0;
  int L=path.size();
  for(int i=0 ; i<L ; ++i)
    switch(hmm->getStateType(path[i])) 
      {
      case PHMM_MATCH:   ++parentL; ++childL; break;
      case PHMM_INSERT:  ++childL; break;
      case PHMM_DELETE:  ++parentL; break;
      }
}



void StatePath::setScore(double s)
{
  score=s;
}



double StatePath::getScore() const
{
  return score;
}



int StatePath::childPosToIndex(int childIndex)
{
  StatePath &self=*this;
  int childI=-1, pathIndex=0; 
  for( ; childI<childIndex ; ++pathIndex) {
    if(emitsChild(hmm->getStateType(self[pathIndex])))
       ++childI;
  }
  return pathIndex-1;
}



void StatePath::writeAsImage(ostream &os)
{
  // Vector<STATE> path
  int L=path.size();
  Vector<STATE>::iterator cur=path.begin(), end=path.end();
  int x=0, y=0;
  for(; cur!=end ; ++cur) {
    STATE q=*cur;
    hmm->updateColumnsFwd(q,x,y);
    os<<x<<" "<<y<<" 0 1"<<endl;
  }
}



void StatePath::writeAsImage(const String &filename)
{
  ofstream os(filename.c_str());
  writeAsImage(os);
}



void StatePath::clear()
{
  path.clear();
}


