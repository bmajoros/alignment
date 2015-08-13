/****************************************************************
 HirschThread.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "HirschThread.H"
#include "Hirschberg.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                    HirschPassThread methods
 ****************************************************************/
HirschPassThread::HirschPassThread(HirschPass &F,bool usePrescan) 
  : F(F), usePrescan(usePrescan)
{ 
  start(); 
}



void HirschPassThread::f() 
{ 
  if(usePrescan)
    F.runPrescan();
  else
    F.run(); 
}




/****************************************************************
                    HirschRecurseThread methods
 ****************************************************************/

HirschRecurseThread::HirschRecurseThread(Hirschberg &H,
					 List<HirschCell>::iterator fromIter,
					 List<HirschCell>::iterator toIter) 
  : H(H), fromIter(fromIter), toIter(toIter) 
{ 
  start(); 
}



void HirschRecurseThread::f() 
{ 
  H.recurse(fromIter,toIter); 
}




