/****************************************************************
 Posteriors.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Posteriors_H
#define INCL_Posteriors_H
#include <iostream>
#include "LinkForward.H"
#include "LinkBackward.H"
using namespace std;
using namespace BOOM;


class Posteriors {
public:
  Posteriors(LinkForward &,LinkBackward &);
  double compute(int i,int j,STATE);
private:
  LinkForward &F;
  LinkBackward &B;
};

#endif

