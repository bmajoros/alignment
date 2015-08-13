/****************************************************************
 HirschbergFrame.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "HirschbergFrame.H"
#include "BOOM/Constants.H"



/****************************************************************
                      WindowColumn methods
 ****************************************************************/

Array2D<float>::RowIn2DArray<float> WindowColumn::operator[](int globalY)
{
  return array[globalY-minY];
}



Array2D<STATE>::RowIn2DArray<STATE> WindowColumn::indexPred(int globalY)
{
  return pred[globalY-minY];
}



void WindowColumn::printOn(ostream &os) const
{
  int dim1=array.getFirstDim(), dim2=array.getSecondDim();
  int ySpan=maxY-minY+1;
  if(dim1>ySpan) dim1=ySpan;
  for(int j=0 ; j<dim2 ; ++j)
    {
      os<<"state "<<j<<": ("<<minY<<"-"<<maxY<<") ";
      for(int i=0 ; i<dim1 ; ++i)
	os << array.index(i,j) << ' ';
      os << endl;
    }

  //os<<array;
}



ostream &operator<<(ostream &os,const WindowColumn &col)
{
  col.printOn(os);
  return os;
}




/****************************************************************
                      HirschbergFrame methods
 ****************************************************************/

HirschbergFrame::HirschbergFrame(int minX,int maxX,int minY,int maxY,
				 int numStates,const BandingPattern &bp,
				 HirschbergDirection dir,BranchHMM &hmm)
  : minX(minX), maxX(maxX), minY(minY), maxY(maxY), numStates(numStates),
    bandingPattern(bp), dir(dir), nextColOffset(dir==FORWARD ? 1 : -1),
    hmm(hmm)
{
  int maxColHeight=bandingPattern.getMaxColumnHeight();
  thisCol=new WindowColumn(maxColHeight,numStates);
  nextCol=new WindowColumn(maxColHeight,numStates);
  
  // Initialize the zeroPillar
  /*
  zeroPillar=new float[numStates];
  for(int i=0 ; i<numStates ; ++i) zeroPillar[i]=LOG_0;
  zeroPillarSize=sizeof(float)*numStates;
  */
}



HirschbergFrame::~HirschbergFrame()
{
  delete thisCol;
  delete nextCol;
  //delete [] zeroPillar;
}



int HirschbergFrame::getMinX() const
{
  return minX;
}



int HirschbergFrame::getMaxX() const
{
  return maxX;
}



int HirschbergFrame::getMinY() const
{
  return minY;
}



int HirschbergFrame::getMaxY() const
{
  return maxY;
}



const BandingPattern &HirschbergFrame::getBandingPattern() const
{
  return bandingPattern;
}



void HirschbergFrame::initWindow(int curColX)
{
  winX=curColX;
  bandingPattern.getBounds(winX,thisCol->minY,thisCol->maxY);
  bandingPattern.getBounds(winX+nextColOffset,nextCol->minY,nextCol->maxY);
}



void HirschbergFrame::advanceWindow()
{
  swapColumns();
  winX+=nextColOffset;
  bandingPattern.getBounds(winX+nextColOffset,nextCol->minY,nextCol->maxY);
}



Array2D<float>::RowIn2DArray<float> HirschbergFrame::operator()(int x,int y)
{
  // ### DEBUGGING -- COMMENT THIS OUT TO IMPROVE SPEED
  if(x!=winX && x!=winX+nextColOffset) {
    cout<<"INTERNAL_ERROR in HirschbergFrame.C"<<endl;
    INTERNAL_ERROR;
    }
  // ###

  WindowColumn &col=*(x==winX ? thisCol : nextCol);
  return col[y];
}



void HirschbergFrame::swapColumns()
{
  WindowColumn *temp=thisCol;
  thisCol=nextCol;
  nextCol=temp;
}



void HirschbergFrame::backtrack(int &x,int &y,STATE &q)
{
  WindowColumn &col=*(x==winX ? thisCol : nextCol);
  STATE predQ=col.indexPred(y)[q];
  hmm.updateColumnsRev(q,y,x);
  q=predQ;
}



