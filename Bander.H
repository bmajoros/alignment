/****************************************************************
 Bander.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_Bander_H
#define INCL_Bander_H
#include <iostream>
#include "BOOM/Array2D.H"
using namespace std;
using namespace BOOM;


class Bander
{
public:
  Bander(int length1,int length2);
protected:
  Array2D<double> matrix;
};

#endif

