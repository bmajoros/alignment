/****************************************************************
 PrecomputedEmissions.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "PrecomputedEmissions.H"
using namespace std;
using namespace BOOM;


PrecomputedEmissions::PrecomputedEmissions(int seqLen,int numFuncClasses,
					   int alphabetSize)
  : A(seqLen,numFuncClasses,alphabetSize)
{
  // ctor
}



void PrecomputedEmissions::resize(int seqLen,int numFuncClasses,
				  int alphabetSize)
{
  A.resize(seqLen,numFuncClasses,alphabetSize);
}





