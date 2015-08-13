/****************************************************************
 WigBinary.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_WigBinary_H
#define INCL_WigBinary_H
#include <iostream>
#include "BOOM/File.H"
#include "BOOM/String.H"
#include "BOOM/GffFeature.H"
using namespace std;
using namespace BOOM;

//typedef pair<int,int> WigInterval;
struct WigInterval {
  int begin, end;
  WigInterval(int b,int e) : begin(b), end(e) {}
  bool overlaps(GffFeature *g) 
    {return begin<g->getEnd() && g->getBegin()<end;}
  bool isEmpty() {return begin>end;}
};

class WigBinary {
public:
  WigBinary(const String &filename);
  float read(int pos); // pos is in residue units (not bytes!)
  void regionsAbove(float threshold,Vector<WigInterval> &into);
  void regionsAbove(float threshold,Vector<WigInterval> &into,
		    int from,int to);
  void regionsAbove(float threshold,Vector<WigInterval> &into,
		    const Vector<WigInterval> &lowerThreshRegions);
  int getLength(); // in residues, not bytes
  void getExtrema(float &min,float &max);
  void close();
protected:
  File F;
};

#endif

