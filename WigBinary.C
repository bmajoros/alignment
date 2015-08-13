/****************************************************************
 WigBinary.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "WigBinary.H"
using namespace std;
using namespace BOOM;

WigBinary::WigBinary(const String &filename)
{
  F.open(filename);
}



float WigBinary::read(int pos)
{
  F.seek(pos*sizeof(float));
  return F.readFloat();
}



void WigBinary::close()
{
  F.close();
}



int WigBinary::getLength()
{
  return F.getSize()/sizeof(float);
}



void WigBinary::regionsAbove(float threshold,Vector<WigInterval> &into)
{
  regionsAbove(threshold,into,0,getLength()-1);
}



void WigBinary::regionsAbove(float threshold,Vector<WigInterval> &into,
			     const Vector<WigInterval> &lowerThreshRegions)
{
  Vector<WigInterval>::const_iterator cur=lowerThreshRegions.begin(),
    end=lowerThreshRegions.end();
  for(; cur!=end ; ++cur) {
    const WigInterval &I=*cur;
    regionsAbove(threshold,into,I.begin,I.end);
  }
}



void WigBinary::regionsAbove(float threshold,Vector<WigInterval> &into,
			     int from,int to)
{
  int begin=-1;
  F.seek(from*sizeof(float));
  for(int pos=from ; pos<=to ; ++pos) {
    float y=F.readFloat();
    //cout<<pos<<" "<<threshold<<" "<<y<<" "<<begin<<" "<<from<<" "<<to<<endl;
    if(begin<0) {
      if(y>=threshold) begin=pos;
    }
    else {
      if(y<threshold) {
	into.push_back(WigInterval(begin,pos));
	begin=-1;
      }
    }
  }
  if(begin>=0) into.push_back(WigInterval(begin,to));
}



void WigBinary::getExtrema(float &min,float &max)
{
  min=max=read(0);
  int L=getLength();
  F.seek(0);
  for(int i=0 ; i<L ; ++i){
    float y=F.readFloat();
    if(y<min) min=y;
    else if(y>max) max=y;
  }
}


