/****************************************************************
 ProfileBackward.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "ProfileBackward.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Constants.H"
#include "BOOM/PureDnaAlphabet.H"
#include "ProfileFelsenstein.H"
using namespace std;
using namespace BOOM;


const double log0=NEGATIVE_INFINITY;
const double log1=log(1.0);


ProfileBackward::ProfileBackward(const PairHMM &hmm,
				 const AlignmentView &insideClade,
				 const AlignmentView &outsideClade,
				 Phylogeny &phylogeny,
				 const Array1D<int> &trackMap,
				 const AlphabetMap &alphabetMap)
  : hmm(hmm), inside(insideClade), outside(outsideClade),
    phylogeny(phylogeny), alphabet(PureDnaAlphabet::global()),
    fullAlignment(insideClade.getAlignment()), gapSymbol(hmm.getGapSymbol()),
    masterCol(hmm.getAlphabet(),hmm.getGapSymbol()), trackMap(trackMap),
    root(*phylogeny.getRoot()), alphabetMap(alphabetMap)
{
  // ctor

  initMasterCol();
  fillMatrix();
}



double ProfileBackward::getLikelihood() const 
{
  return B(0,0,0);
}



double ProfileBackward::operator()(int i,int j,STATE k) const
{
  return B(i,j,k);
}



int ProfileBackward::getFirstDim() const
{
  return B.getFirstDim();
}



int ProfileBackward::getSecondDim() const
{
  return B.getSecondDim();
}



int ProfileBackward::getThirdDim() const
{
  return B.getThirdDim();
}



void ProfileBackward::fillMatrix()
{
  if(!hmm.isInLogSpace()) throw "HMM is not in log space";
  
  // Initialization:
  
  Vector<STATE> I, D, M;
  hmm.getStatesOfType(PHMM_INSERT,I);
  hmm.getStatesOfType(PHMM_DELETE,D);
  hmm.getStatesOfType(PHMM_MATCH,M);
  int nI=I.size(), nD=D.size(), nM=M.size();
  int m=outside.getLength(), n=inside.getLength();
  int numStates=hmm.getNumStates();
  cachedEmitP.resize(numStates,m+1,n+1);
  cachedEmitP.setAllTo(NEGATIVE_INFINITY);
  B.resize(m+1,n+1,numStates);
  B.setAllTo(NEGATIVE_INFINITY);
  for(STATE k=1 ; k<numStates ; ++k)
    B(m,n,k)=hmm.getTransP(k,0);
  Array1D<double> logProbs(nI);
  for(int j=n-1 ; j>=0 ; --j) {
    for(STATE k=1 ; k<numStates ; ++k) {
      for(int ih=0 ; ih<nI ; ++ih) {
	STATE h=I[ih];
	logProbs[ih]=
	  hmm.getTransP(k,h)+getEmitP(h,m,j)+B(m,j+1,h);
      }
      B(m,j,k)=sumLogProbs<double>(logProbs);
    }
  }
  logProbs.resize(nD);
  for(int i=m-1 ; i>=0 ; --i) {
    for(STATE k=1 ; k<numStates ; ++k) {
      for(int ih=0 ; ih<nD ; ++ih) {
	STATE h=D[ih];
	logProbs[ih]=
	  hmm.getTransP(k,h)+B(i+1,n,h)+getEmitP(h,i,n);//###joint
	  //hmm.getTransP(k,h)+B(i+1,n,h); //###conditional
      }
      B(i,n,k)=sumLogProbs<double>(logProbs);
    }
  }

  // Recurrence:

  logProbs.resize(numStates);
  for(int i=m-1 ; i>=0 ; --i) {
    for(int j=n-1 ; j>=0 ; --j) {
      for(STATE k=1 ; k<numStates ; ++k) {
	logProbs.setAllTo(log0);
	for(STATE h=1 ; h<numStates ; ++h) {
	  double trans=hmm.getTransP(k,h);
	  switch(hmm.getStateType(h)) 
	    {
	    case PHMM_MATCH:
	      logProbs[h]=trans+getEmitP(h,i,j)+B(i+1,j+1,h);
	      break;
	    case PHMM_INSERT:
	      logProbs[h]=trans+getEmitP(h,m,j)+B(i,j+1,h);
	      break;
	    case PHMM_DELETE:
	      logProbs[h]=trans+getEmitP(h,i,n)+B(i+1,j,h);
	      break;
	    default: throw "error in ProfileBackward.C";
	    }
	}
	B(i,j,k)=sumLogProbs<double>(logProbs);
      }
    }
  }

  // Termination:
  for(STATE h=1 ; h<numStates ; ++h) {
    double trans=hmm.getTransP(0,h);
    switch(hmm.getStateType(h)) 
      {
      case PHMM_MATCH:
	logProbs[h]=trans+getEmitP(h,0,0)+B(1,1,h);
	break;
      case PHMM_INSERT:
	logProbs[h]=trans+getEmitP(h,m,0)+B(0,1,h);
	break;
      case PHMM_DELETE:
	logProbs[h]=trans;+getEmitP(h,0,n)+B(1,0,h);//###joint
	//logProbs[h]=trans;//###conditional
	break;
      default: throw "error in ProfileBackward.C";
      }
  }
  B(0,0,0)=sumLogProbs<double>(logProbs);
}



void ProfileBackward::initMasterCol()
{
  struct TV : public TreeVisitor {
    MultSeqAlignment &A;
    Symbol gap;
    Map<String,int> &nameToTrackId;
    TV(MultSeqAlignment &A,Symbol gap,Map<String,int> &m) 
      : A(A), gap(gap), nameToTrackId(m) {}
    void f(PhylogenyNode &node) {
      A.findOrCreateTrack(node.getName()).getSeq().append(gap);
      nameToTrackId[node.getName()]=node.getID();
    }
    void processNode(InternalNode &node) {f(node);}
    void processNode(RootNode &node)     {f(node);}
    void processNode(LeafNode &node)     {f(node);}
  } visitor(masterCol,gapSymbol,nameToTrackId);
  phylogeny.postorderTraversal(visitor);
  //initTrackMap();
}



/*
void ProfileBackward::initTrackMap()
{
  int n=phylogeny.getNumNodes();//fullAlignment.getNumTracks();
  trackMap.resize(n);
  trackMap.setAllTo(-1);
  for(int i=0 ; i<n ; ++i) {
    AlignmentSeq &track=fullAlignment.getIthTrack(i);
    trackMap[nameToTrackId[track.getName()]]=i;
  }
}
*/



void ProfileBackward::resetMasterCol()
{
  int numTracks=masterCol.getNumTracks();
  for(int i=0 ; i<numTracks ; ++i)
    masterCol.getIthTrack(i).getSeq()[0]=gapSymbol;
}



double ProfileBackward::getEmitP(STATE k,int col1,int col2)
{
  double cachedValue=cachedEmitP(k,col1,col2);
  if(finite(cachedValue)) return cachedValue;
  resetMasterCol();
  ProfileFelsenstein F(masterCol,alphabet,alphabetMap);
  switch(hmm.getStateType(k)) 
    {
    case PHMM_MATCH:
      installColumn(col1,outside);
      installColumn(col2,inside);
      break;
    case PHMM_INSERT: // ### would be faster to run on just the clade
      installColumn(col2,inside);
      break;
    case PHMM_DELETE:
      installColumn(col1,outside);
      break;
    default: throw "error in ProfileBackward";
    }
  return cachedEmitP(k,col1,col2)=F.logLikelihood(0,&root,k);
}



void ProfileBackward::installColumn(int col,const AlignmentView &view)
{
  col=view.mapColumn(col);
  const BitSet &taxonIDs=view.getTaxonIDs();
  int numTaxa=phylogeny.getNumNodes();
  for(int i=0 ; i<numTaxa ; ++i)
    if(taxonIDs.isMember(i)) {
      int targetTrack=trackMap[i];
      if(targetTrack<0) continue;
      masterCol.getIthTrack(i).getSeq()[0]=
	fullAlignment.getIthTrack(targetTrack).getSeq()[col];
    }
}



double ProfileBackward::getCachedEmitP(STATE k,int col1,int col2)
{
  return cachedEmitP(k,col1,col2);
}

