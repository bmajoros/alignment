/****************************************************************
 StatePath.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_StatePath_H
#define INCL_StatePath_H
#include <iostream>
#include "State.H"
#include "BOOM/Vector.H"
#include "BOOM/String.H"
using namespace std;
using namespace BOOM;

class PairHMM;

class StatePath
{
  Vector<STATE> path; // usually doesn't include the start-stop state (q0)
  const PairHMM *hmm;
  double score;
public:
  StatePath(const PairHMM *hmm=NULL,int length=0);
  void push_back(STATE);
  int size() const;
  int length(); // same as size();
  void getSeqLengths(int &parentLen,int &childLen) const;
  STATE operator[](int) const;
  STATE &operator[](int);
  StatePath *getReverse() const;
  void setHMM(const PairHMM *);
  const PairHMM *getHMM() const;
  void printOn(ostream &) const;
  void setScore(double);
  double getScore() const;
  int childPosToIndex(int childIndex);
  void writeAsImage(ostream &);
  void writeAsImage(const String &filename);
  void clear();
};

ostream &operator<<(ostream &,const StatePath &);

#endif

