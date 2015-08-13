/****************************************************************
 HirschPosteriors.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "HirschPosteriors.H"
#include "HirschForwardSum.H"
#include "HirschBackwardSum.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                       HirschPosteriors methods
 ****************************************************************/
HirschPosteriors::HirschPosteriors(Taxon *parent,Taxon *child,
				   LinkFelsenstein &F,int bandwidth,
				   int parentLen,int childLen,
				   bool usePrescan,int maxThreads,
				   ContentSensor *contentSensor)
  : Hirschberg(parent,child,F,bandwidth,parentLen,childLen,usePrescan,
	       maxThreads,contentSensor)
{
  vc=VC_UNCONSTRAINED;
}



double HirschPosteriors::computeLikelihood()
{
  if(!parent || !child) INTERNAL_ERROR;
  int a=0, b=m, c=0, d=n;
  HirschForwardSum F(c,d,a,b,bandingPattern,*hmm,QI,numI,QD,
		     numD,QM,numM,vc,precomputedEmissions,parentFP,
		     childFP,upMap,downMap,this,contentSensor);
  double LL=F.computeLikelihood();
  if(!isFinite(LL)) {
    cout<<"debugging"<<endl;
    HirschForwardSum F(c,d,a,b,bandingPattern,*hmm,QI,numI,QD,
		       numD,QM,numM,vc,precomputedEmissions,parentFP,
		       childFP,upMap,downMap,this,contentSensor);
    F.debug();
  }
  return LL;
}



double HirschPosteriors::compute(int siteBegin,int siteEnd,BranchEnd
				 parentChild)
{
  double LL=computeLikelihood();
  double numerator=computeNumerator(siteBegin,siteEnd,parentChild);
  double posterior=numerator-LL;
  return posterior;
}



double HirschPosteriors::computeNumerator(int siteBegin,int siteEnd,BranchEnd
					  parentChild)
{
  ViterbiConstraint vc=parentChild==PARENT ? VC_PARENT_PARSE : VC_CHILD_PARSE;
  HirschForwardSum F2(0,n,0,m,bandingPattern,*hmm,QI,numI,QD,
		     numD,QM,numM,vc,precomputedEmissions,parentFP,
		     childFP,upMap,downMap,this,contentSensor);
  F2.setConstraintInterval(siteBegin,siteEnd);
  double numerator=F2.computeLikelihood();
  return numerator;
}

