/****************************************************************
 Hirschberg.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "Hirschberg.H"
#include "BOOM/Constants.H"
#include "BOOM/Exceptions.H"
#include "HirschForwardMax.H"
#include "HirschBackwardMax.H"
#include "HirschForwardSum.H"
#include "HirschBackwardSum.H"
#include "HirschThread.H"
using namespace std;
using namespace BOOM;

// change Max to Sum in these next two lines to use posterior decoding:
#define HirschForward HirschForwardMax
#define HirschBackward HirschBackwardMax


// DO NOT CHANGE THIS LINE:
#define HFHFHF HirschForwardMax


/****************************************************************
                       HirschCell methods
 ****************************************************************/
HirschCell::HirschCell(int x,int y,STATE q,float forward,float backward)
  : x(x), y(y), forward(forward), backward(backward), state(q)
{
  // ctor
}



/****************************************************************
                       Hirschberg methods
 ****************************************************************/
Hirschberg::Hirschberg(Taxon *parent,Taxon *child,LinkFelsenstein &F,
		       int bandwidth,int parentLen,int childLen,
		       bool usePrescan,int maxThreads,ContentSensor *
		       contentSensor)
  : parent(parent), child(child), F(F), vc(VC_UNCONSTRAINED), 
    alphabet(F.getAlphabet()), posteriors(NULL),
    bandingPattern(parentLen+1,childLen+1,bandwidth),
    maxThreads(maxThreads), usePrescan(usePrescan),
    contentSensor(contentSensor)
{
  numAlpha=alphabet.size();
  if(parent && child) {
    m=parent->getSeqLen();
    n=child->getSeqLen();
    parentFP=&parent->getFunctionalParse();
    childFP=&child->getFunctionalParse();
    branch=parent->getBranchToChild(*child);
    hmm=branch->getHMM();
    initStateSets();
    numStates=hmm->getNumStates();
    upMap=&branch->getUpMap();
    downMap=&branch->getDownMap();
    precomputeEmissions();
  }
  else {
    parentFP=childFP=NULL;
    branch=NULL;
    hmm=NULL;
    numStates=0;
    m=n=0;
    upMap=downMap=NULL;
  }
}



Hirschberg::~Hirschberg()
{
}



PrecomputedEmissions &Hirschberg::getPrecomputedEmissions(int which)
{
  return precomputedEmissions[which];
}



void Hirschberg::initStateSets()
{
  hmm->getStatesOfType(PHMM_MATCH,QM);
  hmm->getStatesOfType(PHMM_INSERT,QI);
  hmm->getStatesOfType(PHMM_DELETE,QD);
  numM=QM.size();
  numI=QI.size();
  numD=QD.size();
}



StatePath *Hirschberg::decode(ViterbiConstraint vc)
{
  // INITIALIZATION:

  this->vc=vc;
  HirschCell firstCell(0,0,0,LOG_1);
  computedCells.push_back(firstCell);

  // RECURSION:

  recurse(computedCells.begin(),computedCells.end());

  // TERMINATION:
  StatePath *path=new StatePath(hmm);
  List<HirschCell>::iterator cur=computedCells.begin(), 
    end=computedCells.end();
  float score=LOG_0;
  for(; cur!=end ; ++cur) {
    const HirschCell &cell=*cur;
    if(cell.state>0) path->push_back(cell.state);
    if(cell.state>0) score=cell.forward+cell.backward;
  }
  path->setScore(score);
  return path;
}



void Hirschberg::recurse(List<HirschCell>::iterator fromIter,
			 List<HirschCell>::iterator toIter)
{
  // Get bounds a, b, c, and d
  HirschCell &fromCell=*fromIter;
  int a=fromCell.y, b, c=fromCell.x, d;
  STATE toState;
  float toValue;
  computedCellsMutex.wait();
  bool endException=(toIter==computedCells.end());
  computedCellsMutex.signal();
  if(endException) {
    b=m;
    d=n;
  }
  else {
    HirschCell &toCell=*toIter;
    b=toCell.y;
    d=toCell.x;
    toState=toCell.state;
    toValue=toCell.backward;
  }

  // Detect end of recursion
  if(c>d-2) {
    endRecursion(fromIter,toIter,a,b,c,d);
    return;
  }

  // Partition and run forward & backward viterbi algorithms
  int jStar=(c+d)/2;
  HirschForward F(c,jStar,a,b,bandingPattern,*hmm,QI,numI,QD,
		  numD,QM,numM,vc,precomputedEmissions,parentFP,
		  childFP,upMap,downMap,this,contentSensor);
  HirschBackward B(jStar,d,a,b,bandingPattern,*hmm,QI,numI,QD,
		   numD,QM,numM,vc,precomputedEmissions,parentFP,
		   childFP,upMap,downMap,this,contentSensor);
  HirschbergFrame &fFrame=F.getFrame(), &bFrame=B.getFrame();
  {
    WindowColumn &fCol=fFrame.getThisCol(), &bCol=bFrame.getThisCol();
    Array2D<float>::RowIn2DArray<float> bPillar=bCol[b];
    if(endException) {
      //bCol.setAllTo(LOG_0); // ### DEBUGGING
      //bFrame.getNextCol().setAllTo(LOG_0); // ### DEBUGGING
      for(STATE q=0 ; q<numStates ; ++q) bPillar[q]=getTransP(q,0);
    }
    else {
      bPillar.setAllTo(LOG_0);
      bPillar[toState]=toValue;
    }
    Array2D<float>::RowIn2DArray<float> fPillar=fCol[fromCell.y];
    //fCol.setAllTo(LOG_0); // ### DEBUGGING
    //fFrame.getNextCol().setAllTo(LOG_0); // ### DEBUGGING
    fPillar.setAllTo(LOG_0);
    fPillar[fromCell.state]=fromCell.forward;
    runFB(F,B);
  }

  // Find optimal crossing point
  WindowColumn &fCol=fFrame.getCol(jStar), &bCol=bFrame.getCol(jStar);
  int minY=max(a,fCol.minY), maxY=min(b,fCol.maxY);
  float bestP=LOG_0;
  int iStar;
  STATE qStar;
  for(int y=minY ; y<=maxY ; ++y) {
    Array2D<float>::RowIn2DArray<float> fPillar=fCol[y], bPillar=bCol[y];
    for(STATE q=0 ; q<numStates ; ++q) {
      float logP=fPillar[q]+bPillar[q];
      if(logP>bestP) {
	bestP=logP;
	iStar=y;
	qStar=q;
      }
    }
  }
  if(isInfinity(bestP)) {

    //if(allowZeroPath) throw HIRSCHBERG_NO_PATH;

    cout<<"optimal crossing point in Hirschberg has zero probability"
	<<endl; 
    cout<<"minY="<<minY<<", maxY="<<maxY<<endl;
    for(int y=minY ; y<=maxY ; ++y) {
      Array2D<float>::RowIn2DArray<float> fPillar=fCol[y], bPillar=bCol[y];
      for(STATE q=0 ; q<numStates ; ++q) {
	if(isFinite(fPillar[q]) || isFinite(bPillar[q])) 
	  cout<<"  j*="<<jStar<<"  y="<<y<<" q="<<q<<" fw="<<fPillar[q]<<" bw="<<bPillar[q]<<endl;
      }
    }
    throw HIRSCHBERG_NO_PATH(); // ### DEBUGGING
    INTERNAL_ERROR;
  }	
  //cout<<"(i*,j*,q*)=("<<iStar<<","<<jStar<<","<<qStar<<") logP="<<bestP<<" f="<<fCol[iStar][qStar]<<" b="<<bCol[iStar][qStar]<<" "<<hmm->getStateType(qStar)<<endl;

  // Insert computed cell
  HirschCell newCell(jStar,iStar,qStar,fCol[iStar][qStar],bCol[iStar][qStar]);
  computedCellsMutex.wait();
  list<HirschCell>::iterator newIter=
    computedCells.insertBefore(toIter,newCell);
  computedCellsMutex.signal();

  // Recurse to left and right subproblems
  if(maxThreads.tryWait()) {
    HirschRecurseThread hrt(*this,fromIter,newIter);
    recurse(newIter,toIter);
    hrt.join();
    maxThreads.post();
  }
  else {
    recurse(fromIter,newIter);
    recurse(newIter,toIter);
  }
  return;
}



void Hirschberg::runFB(HirschPass &F,HirschPass &B)
{
  if(maxThreads.tryWait()) {
    HirschPassThread hft(F,usePrescan);
    if(usePrescan) B.runPrescan(); else B.run();
    hft.join();
    maxThreads.post();
  }
  else {
    if(usePrescan) {
      F.runPrescan(); 
      B.runPrescan(); 
    }
    else { 
      F.run(); 
      B.run(); 
    }
  }
}



void Hirschberg::endRecursion(List<HirschCell>::iterator fromIter,
			      List<HirschCell>::iterator toIter,
			      int a,int b,int c,int d)
{
  // PRECONDITION: only one or two columns remain

  // First, run the forward viterbi on the remaining column(s)
  HFHFHF F(c,d,a,b,bandingPattern,*hmm,QI,numI,QD,
		  numD,QM,numM,vc,precomputedEmissions,parentFP,
		  childFP,upMap,downMap,this,contentSensor);
  HirschbergFrame &fFrame=F.getFrame();
  //if(toIter==computedCells.end()) INTERNAL_ERROR; // ### debugging
  HirschCell &fromCell=*fromIter; 
  WindowColumn &fCol=fFrame.getThisCol();
  fCol[a].setAllTo(LOG_0);
  fCol[a][fromCell.state]=fromCell.forward;
  if(usePrescan) F.runPrescan(); else F.run();
  
  // Now backtrack and allocate completed cells
  int x=d, y=b;
  STATE q=-1;
  if(toIter==computedCells.end()) {
    HirschbergFrame &frame=F.getFrame();
    Array2D<float>::RowIn2DArray<float> pillar=frame(x,y);
    float bestP=LOG_0;
    for(STATE qStar=1 ; qStar<numStates ; ++qStar) {
      float logP=pillar[qStar]+getTransP(qStar,0);
      if(logP>bestP) { bestP=logP; q=qStar; }
    }
    HirschCell newCell(x,y,q,bestP);
    toIter=computedCells.insertBefore(toIter,newCell);
  }
  else q=(*toIter).state;
  List<HirschCell>::iterator iter=toIter;
  while(x>c || y>a) {
    if(x<0 || y<0) INTERNAL_ERROR;
    fFrame.backtrack(x,y,q);
    if(x==c && y==a) break;
    float logP=fFrame(x,y)[q];
    HirschCell newCell(x,y,q,logP);
    computedCellsMutex.wait(); // ###
    iter=computedCells.insertBefore(iter,newCell);
    computedCellsMutex.signal(); // ###
  }
}



void Hirschberg::precomputeEmissions()
{
  precomputeEmissions(child,CHILD,INSIDE);
  precomputeEmissions(parent,PARENT,OUTSIDE,child);
}



void Hirschberg::precomputeEmissions(Taxon *taxon,BranchEnd parentChild,
				      InsideOutside insideOutside,
				      Taxon *avoid)
{
  PrecomputedEmissions &E=precomputedEmissions[parentChild];

  // First, allocate the array
  int L=taxon->getSeqLen();
  int numFC=FunctionalClass::numClasses();
  E.resize(L,numFC,numAlpha);

  // Iterate over all sequence positions
  Array1D<float> temp(numAlpha);
  for(int pos=0 ; pos<L ; ++pos) {
    Array3D<float>::IndexedOnce<float> slice=E[pos];
    for(int fc=0 ; fc<numFC ; ++fc) {
      Array3D<float>::IndexedTwice<float> row=slice[fc];
      switch(insideOutside) {
      case INSIDE:  
	F.precomputeInside(temp,*taxon,pos,fc); 
	break;
      case OUTSIDE: 
	F.precomputeOutside(temp,*taxon,pos,*avoid,fc); 
	break;
      }
      for(Symbol s=0 ; s<numAlpha ; ++s) row[s]=temp[s];
    }
  }
  taxon->getPrecomputedEmissions()=E;
}



void Hirschberg::usePosteriors(Posteriors *p)
{
  posteriors=p;
}



float Hirschberg::getTransP(STATE h,STATE k)
{
  float p=hmm->getTransP(h,k);
  if(posteriors) return finite(p) ? LOG_1 : LOG_0;
  else return p;
}



float Hirschberg::scorePath(StatePath &path)
{
  int a=0, b=m, c=0, d=n;
  HirschForward F(c,d,a,b,bandingPattern,*hmm,QI,numI,QD,
		  numD,QM,numM,vc,precomputedEmissions,parentFP,
		  childFP,upMap,downMap,this,contentSensor);
  return F.scorePath(path);
}

