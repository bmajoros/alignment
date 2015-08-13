/****************************************************************
 PosteriorMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <fstream>
#include <iostream>
#include "PosteriorMatrix.H"
#include "BOOM/File.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


PosteriorMatrix::PosteriorMatrix(BandedForward &F,BandedBackward &B)
  : M(F.getFirstDim(),F.getSecondDim(),F.getThirdDim())
{
  compute(F,B);
}



double PosteriorMatrix::operator()(int i,int j,int k)
{
  return M(i,j,k);
}



void PosteriorMatrix::compute(BandedForward &F,BandedBackward &B)
{
  int L1=F.getFirstDim(), L2=F.getSecondDim(), numStates=F.getThirdDim();
  double LL=F.getLikelihood();
  for(int i=0 ; i<L1 ; ++i)
    for(int j=0 ; j<L2 ; ++j)
      for(int k=0 ; k<numStates ; ++k)
	M(i,j,k)=F(i,j,k)+B(i,j,k)-LL;
}



bool PosteriorMatrix::save(double cutoff,const String &filename)
{
  //double cutoff=log(0.001);
  int L1=M.getFirstDim(), L2=M.getSecondDim(), numStates=M.getThirdDim();
  ofstream os(filename.c_str());
  if(!os.good()) return false;
  for(int i=0 ; i<L1 ; ++i)
    for(int j=0 ; j<L2 ; ++j)
      for(int k=0 ; k<numStates ; ++k) {
	double x=M(i,j,k);
	if(x>cutoff) os<<i<<" "<<j<<" "<<k<<" "<<M(i,j,k)<<endl;
      }
  return true;
}



bool PosteriorMatrix::saveBinary(double cutoff,const String &filename)
{
  //double cutoff=log(0.001);
  short L1=M.getFirstDim(), L2=M.getSecondDim(), numStates=M.getThirdDim();
  File file(filename,"w");
  file.write(L1);
  file.write(L2);
  file.write(numStates);
  int numEntries=0;
  for(short i=0 ; i<L1 ; ++i)
    for(short j=0 ; j<L2 ; ++j)
      for(short k=0 ; k<numStates ; ++k)
	if(M(i,j,k)>cutoff) 
	  ++numEntries;
  file.write(numEntries);
  for(short i=0 ; i<L1 ; ++i)
    for(short j=0 ; j<L2 ; ++j)
      for(short k=0 ; k<numStates ; ++k) {
	float x=M(i,j,k);
	if(x>cutoff) {
	  file.write(i);
	  file.write(j);
	  file.write(k);
	  file.write(x);
	}
      }
  return true;
}



PosteriorMatrix *PosteriorMatrix::loadBinary(const String &filename)
{
  PosteriorMatrix *PM=new PosteriorMatrix;
  PM->loadFromBinary(filename);
  return PM;
}



void PosteriorMatrix::loadFromBinary(const String &filename)
{
  File file(filename);
  int L1=file.readShort();
  int L2=file.readShort();
  int numStates=file.readShort();
  int numEntries=file.readInt();
  M.resize(L1,L2,numStates);
  M.setAllTo(NEGATIVE_INFINITY);
  for(int e=0 ; e<numEntries ; ++e) {
    int i=file.readShort();
    int j=file.readShort();
    int k=file.readShort();
    M(i,j,k)=file.readFloat();
  }
}


