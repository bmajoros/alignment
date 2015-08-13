/****************************************************************
 HirschbergFrame.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_HirschbergFrame_H
#define INCL_HirschbergFrame_H
#include <iostream>
#include "BOOM/Array2D.H"
#include "BOOM/Exceptions.H"
#include "BandingPattern.H"
#include "PairHMM/State.H"
#include "BranchHMM.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                      class WindowColumn
 ****************************************************************/
struct WindowColumn {
  WindowColumn(int maxHeight,int numStates) 
    : array(maxHeight,numStates), pred(maxHeight,numStates) {}
  Array2D<float>::RowIn2DArray<float> operator[](int globalY);
  Array2D<STATE>::RowIn2DArray<STATE> indexPred(int globalY);
  void setAllTo(float f) { array.setAllTo(f); }
  Array2D<float> array;
  Array2D<STATE> pred;
  int minY, maxY;
  void printOn(ostream &) const;
};
ostream &operator<<(ostream &,const WindowColumn &);



/****************************************************************
                     enum HirschbergDirection
 ****************************************************************/
enum HirschbergDirection {
  FORWARD,
  BACKWARD
};



/****************************************************************
                      class HirschbergFrame

 A three-dimensional dynamic-programming matrix for the banded
 Hirschberg algorithm (i.e., banded in the first two dimensions).
 ****************************************************************/

class HirschbergFrame {
public:
  HirschbergFrame(int minX,int maxX,int minY,int maxY,
		  int numStates,const BandingPattern &,
		  HirschbergDirection,BranchHMM &);
  virtual ~HirschbergFrame();
  int getMinX() const;
  int getMaxX() const;
  int getMinY() const;
  int getMaxY() const;
  const BandingPattern &getBandingPattern() const;
  void initWindow(int curColX);
  void advanceWindow();
  Array2D<float>::RowIn2DArray<float> operator()(int x,int y);
  inline WindowColumn &getThisCol();
  inline WindowColumn &getNextCol();
  inline WindowColumn &getCol(int x);
  inline int getNumStates() const;
  void backtrack(int &x,int &y,STATE &q);
  //inline void zeroOut(Array2D<float>::RowIn2DArray<float> &pillar);
  inline int getThisColPos() const {return winX;}
protected:
  //float *zeroPillar;
  int zeroPillarSize;
  WindowColumn *thisCol, *nextCol;
  int winX; // location of thisCol
  int minX, maxX, minY, maxY; // theoretical bounds of this submatrix
  int numStates;
  HirschbergDirection dir;
  int nextColOffset; // +1 if dir==FORWARD, -1 if dir==BACKWARD
  const BandingPattern &bandingPattern;
  BranchHMM &hmm;
  void swapColumns();
};



/*
inline void HirschbergFrame::zeroOut(Array2D<float>::RowIn2DArray<float> 
				     &pillar)
{
  memcpy(pillar.getRawData(),zeroPillar,zeroPillarSize);
}
*/



inline WindowColumn &HirschbergFrame::getThisCol()
{
  return *thisCol;
}



inline WindowColumn &HirschbergFrame::getNextCol()
{
  return *nextCol;
}



inline int HirschbergFrame::getNumStates() const
{
  return numStates;
}



inline WindowColumn &HirschbergFrame::getCol(int x)
{
  /*if(x!=winX && x!=winX+nextColOffset) {
    cout<<"INTERNAL_ERROR in HirschbergFrame.H"<<endl;
    INTERNAL_ERROR; 
    }*/

  //cout<<"winX="<<winX<<endl;
  return *(x==winX ? thisCol : nextCol);
}



#endif
