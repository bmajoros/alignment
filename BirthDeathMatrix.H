/****************************************************************
 BirthDeathMatrix.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_BirthDeathMatrix_H
#define INCL_BirthDeathMatrix_H
#include <iostream>
#include "PhyLib/RateMatrix.H"
using namespace std;
namespace BOOM {


class BirthDeathMatrix : public RateMatrix {
public:
  BirthDeathMatrix(float birthRate,float deathRate);
  void exponentiate(float time,float &birthProb,float &deathProb);
  virtual void save(ostream &);
  virtual double getIthParm(int i) const;
  virtual void setIthParm(int i,double parm);
  virtual RateMatrix *average(const Array1D<RateMatrix*> &) const;
  virtual void addNoise(GSL::ContinuousDistribution &);
  virtual RateMatrix *clone() const;
  virtual void init();
protected:
  float birthRate, deathRate;
};


}
#endif

