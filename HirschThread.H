/****************************************************************
 HirschThread.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_HirschThread_H
#define INCL_HirschThread_H
#include <iostream>
#include "BOOM/Thread.H"
#include "HirschPass.H"
#include "Hirschberg.H"
using namespace std;
using namespace BOOM;

class Hirschberg;
class HirschCell;

/****************************************************************
                       HirschPassThread
 ****************************************************************/
class HirschPassThread : public Thread {
  HirschPass &F;
  bool usePrescan;
public:
  HirschPassThread(HirschPass &,bool usePrescan);
  virtual void f();
};



/****************************************************************
                       HirschRecurseThread
 ****************************************************************/
class HirschRecurseThread : public Thread {
  Hirschberg &H;
  List<HirschCell>::iterator fromIter, toIter;
public:
  HirschRecurseThread(Hirschberg &H,List<HirschCell>::iterator fromIter,
		      List<HirschCell>::iterator toIter);
  virtual void f();
};

#endif

