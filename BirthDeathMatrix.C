/****************************************************************
 BirthDeathMatrix.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BirthDeathMatrix.H"
#include "PhyLib/RateMatrix.H"
#include "BOOM/BinaryAlphabet.H"
#include "BOOM/Exceptions.H"
using namespace std;
using namespace BOOM;


BirthDeathMatrix::BirthDeathMatrix(float birthRate,float deathRate)
  : RateMatrix(DNA,BinaryAlphabet::global(),MT_GAINLOSS),
    birthRate(birthRate), deathRate(deathRate)
{
  init();
}



void BirthDeathMatrix::init()
{
  M(0,1)=birthRate;
  M(1,0)=deathRate;
  M(0,0)=-birthRate;
  M(1,1)=-deathRate;
}



void BirthDeathMatrix::save(ostream &os)
{
  os<<birthRate<<"\t"<<deathRate<<endl;
}



double BirthDeathMatrix::getIthParm(int i) const
{
  return i==0 ? birthRate : deathRate;
}



void BirthDeathMatrix::setIthParm(int i,double parm)
{
  if(i==0) birthRate=parm;
  else deathRate=parm;
}



RateMatrix *BirthDeathMatrix::average(const Array1D<RateMatrix*> &A) const
{
  int n=A.size();
  float sumB=0, sumD=0;
  for(int i=0 ; i<n ; ++i) {
    BirthDeathMatrix *M=dynamic_cast<BirthDeathMatrix*>(A[i]);
    if(!M) INTERNAL_ERROR;
    sumB+=M->birthRate;
    sumD+=M->deathRate;
  }
  float aveB=sumB/n, aveD=sumD/n;
  return new BirthDeathMatrix(aveB,aveD);
}



void BirthDeathMatrix::addNoise(GSL::ContinuousDistribution &d)
{
  birthRate+=d.random();
  deathRate+=d.random();
  init();
}



RateMatrix *BirthDeathMatrix::clone() const
{
  return new BirthDeathMatrix(birthRate,deathRate);
}


