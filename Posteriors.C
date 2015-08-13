/****************************************************************
 Posteriors.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "Posteriors.H"
using namespace std;
using namespace BOOM;

Posteriors::Posteriors(LinkForward &F,LinkBackward &B)
  : F(F), B(B)
{
  // ctor
}



double Posteriors::compute(int i,int j,STATE q)
{
  return F(i,j,q)+B(i,j,q)-B(0,0,0);
}



