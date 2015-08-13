/****************************************************************
 BandingPattern.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BandingPattern.H"
using namespace std;
using namespace BOOM;



BandingPattern::BandingPattern(int yDim,int xDim,int bandWidth)
  : xDim(xDim), yDim(yDim), bandwidth(bandWidth), bounds(xDim+1)
{
  computeBounds();
}



int BandingPattern::getMaxColumnHeight() const
{
  return maxColHeight;
}



void BandingPattern::getBounds(int x,int &minY,int &maxY) const
{
  const Bounds &b=bounds[x];
  minY=b.minY;
  maxY=b.maxY;
}



int BandingPattern::getXdim() const
{
  return xDim;
}



int BandingPattern::getYdim() const
{
  return yDim;
}



int BandingPattern::getBandwidth() const
{
  return bandwidth;
}



void BandingPattern::computeBounds()
{
  if(xDim==0 || yDim==0 || bandwidth==0) return;
  float m=yDim/float(xDim);
  if(bandwidth<m) bandwidth=int(m+1);
  maxColHeight=0;
  //bounds[0].minY=bounds[0].maxY=-1;
  for(int x=0 ; x<=xDim ; ++x) {
    int y=int(m*x);
    int beginY=y-bandwidth;
    int endY=y+bandwidth;
    if(beginY<0) beginY=0;
    if(endY>yDim) endY=yDim;
    Bounds &b=bounds[x];
    b.maxY=endY;
    b.minY=beginY;
    int colHeight=endY-beginY+1;
    if(colHeight>maxColHeight) maxColHeight=colHeight;
  }
}


