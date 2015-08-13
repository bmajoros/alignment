/****************************************************************
 BandedFB_Base.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BandedFB_Base.H"
using namespace std;
using namespace BOOM;

BandedFB_Base::BandedFB_Base(const Alphabet &alphabet,
			     const AlphabetMap &alphabetMap,
			     BOOM::Symbol gap,
			     const BranchHMM &hmm,
			     BandingPattern &bandingPattern,
			     Sequence &left,
			     Sequence &right)
  : alphabet(alphabet), alphabetMap(alphabetMap), gapSymbol(gap),
    hmm(hmm), bandingPattern(bandingPattern), left(left), right(right),
    numStates(hmm.getNumStates()), M(1,1,1)
{
  hmm.getStatesOfType(PHMM_INSERT,QI);
  hmm.getStatesOfType(PHMM_DELETE,QD);
  hmm.getStatesOfType(PHMM_MATCH,QM);
  nI=QI.size(); nD=QD.size(); nM=QM.size();
}



double BandedFB_Base::getEmitP(STATE k,int leftIndex,int rightIndex)
{
  /*
  if(left[leftIndex]>4) throw String("left[leftIndex] ")+leftIndex+" "+int(left[leftIndex])+" "+left.getLength();
  if(right[rightIndex]>4) throw String("right[rightIndex] ")+rightIndex+" "+int(right[rightIndex])+" "+right.getLength();
  */

  double P;
  const Array1D<double> &eqFreqs=hmm.getEqFreqs(k);
  switch(hmm.getStateType(k)) 
    {
    case PHMM_MATCH: {
      Symbol a=left[leftIndex], b=right[rightIndex];
      SubstitutionMatrix &Pt=*hmm.getSubstMatrix(k);
      double eq=log(eqFreqs[a]);
      double subst=Pt(a,b);
      //cout<<"eq="<<eq<<" subst="<<subst<<endl;
      P=eq+subst;
      break;
    }
    case PHMM_INSERT: {
      Symbol a=right[rightIndex];//left[leftIndex];
      double eq=log(eqFreqs[a]);
      //cout<<"eq="<<eq<<endl;
      P=eq;
      break;
    }
    case PHMM_DELETE: {
      Symbol b=left[leftIndex];//right[rightIndex];
      double eq=log(eqFreqs[b]);
      //cout<<"eq="<<eq<<endl;
      P=eq;
      break;
    }
    default: throw "error in BandedFB_Base::getEmitP";
    }
  return P;
}



int BandedFB_Base::getFirstDim() const
{
  return M.getFirstDim();
}



int BandedFB_Base::getSecondDim() const
{
  return M.getSecondDim();
}



int BandedFB_Base::getThirdDim() const
{
  return M.getThirdDim();
}



double BandedFB_Base::operator()(int i,int j,STATE k) const
{
  return M(i,j,k);
}



