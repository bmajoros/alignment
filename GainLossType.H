/****************************************************************
 GainLossType.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_GainLossType_H
#define INCL_GainLossType_H
#include <iostream>
#include "BOOM/String.H"
using namespace std;
using namespace BOOM;

enum GainLossType {
  GLT_VOID, // => no functional element in parent or child
  GLT_GAIN,
  GLT_LOSS,
  GLT_RETENTION
};

ostream &operator<<(ostream &,GainLossType);
String toString(GainLossType);

#endif
