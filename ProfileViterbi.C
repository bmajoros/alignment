/****************************************************************
 ProfileViterbi.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "ProfileViterbi.H"
#include "BOOM/Constants.H"
#include "BOOM/DnaDashAlphabet.H"
#include "BOOM/PureDnaAlphabet.H"
#include "ProfileFelsenstein.H"
using namespace std;
using namespace BOOM;

const int GAP=-1;


ProfileViterbi::ProfileViterbi(const BranchHMM &hmm,
			       const MultSeqAlignment &S1,
			       PhylogenyNode *child1,
			       const MultSeqAlignment &S2,
			       PhylogenyNode *child2,
			       Phylogeny &phylogeny,
			       const AlphabetMap &alphabetMap,
			       BandingType bandingType,
			       int bandwidth)
  : hmm(hmm), S1(S1), S2(S2), child1(child1), child2(child2),
    V(S1.getLength()+1,S2.getLength()+1,hmm.getNumStates()),
    T(S1.getLength()+1,S2.getLength()+1,hmm.getNumStates()),
    L1(S1.getLength()), L2(S2.getLength()), phylogeny(phylogeny),
    gap(hmm.getGapSymbol()), numStates(hmm.getNumStates()),
    masterCol(hmm.getAlphabet(),hmm.getGapSymbol()),
    alphabet(PureDnaAlphabet::global()), alphabetMap(alphabetMap),
    bandingType(bandingType), bandwidth(bandwidth)
{
  parentTaxon=static_cast<Taxon*>(child1->getParent()->getDecoration());
  setupMasterCol();
  if(!hmm.isInLogSpace()) throw "HMM not in log space";
  F=new ProfileFelsenstein(masterCol,alphabet,alphabetMap);
  fillMatrix();
}



double ProfileViterbi::operator()(int i,int j,STATE k) const
{
  return V(i,j,k);
}



int ProfileViterbi::getFirstDim() const
{
  return V.getFirstDim();
}



int ProfileViterbi::getSecondDim() const
{
  return V.getSecondDim();
}



int ProfileViterbi::getThirdDim() const
{
  return V.getThirdDim();
}



void ProfileViterbi::fillMatrix()
{
  // Misc initialization:
  Vector<STATE> I, D, M;
  hmm.getStatesOfType(PHMM_INSERT,I);
  hmm.getStatesOfType(PHMM_DELETE,D);
  hmm.getStatesOfType(PHMM_MATCH,M);
  int nI=I.size(), nD=D.size(), nM=M.size();

  // Initialization of DP matrix:
  cout<<"initialization"<<endl;
  V.setAllTo(LOG_0);
  T.setAllTo(-1);
  V(0,0,0)=LOG_1;
  for(int i=1 ; i<=L1 ; ++i)
    for(int ik=0 ; ik<nD ; ++ik) { // deletion states
      STATE k=D[ik], bestH=0;
      double Pe=getEmitP(k,i-1,GAP);
      double bestScore=V(i-1,0,0)+hmm.getTransP(0,k);
      const Array1D<STATE> &pred=hmm.getPredecessors(k);
      STATE *p=&pred[0];
      int numPred=pred.size();
      //cout<<k<<" has "<<numPred<<" D preds: "<<pred<<endl;
      for(int s=0 ; s<numPred ; ++s, ++p) {
	STATE h=*p;
	if(hmm.getStateType(h)!=PHMM_DELETE) continue;
	double score=V(i-1,0,h)+hmm.getTransP(h,k);
	if(score>bestScore) { bestH=h; bestScore=score; }
      }
      /*
      for(int ih=0 ; ih<nD ; ++ih) {
	STATE h=D[ih];
	double score=V(i-1,0,h)+hmm.getTransP(h,k);
	if(score>bestScore) { bestH=h; bestScore=score; }
      }
      */
      V(i,0,k)=bestScore+Pe;
      T(i,0,k)=bestH;
    }
  for(int j=1 ; j<=L2 ; ++j) 
    for(int ik=0 ; ik<nI ; ++ik) { // insertion states
      STATE k=I[ik], bestH=0;
      double Pe=getEmitP(k,GAP,j-1);
      double bestScore=V(0,j-1,0)+hmm.getTransP(0,k);
      const Array1D<STATE> &pred=hmm.getPredecessors(k);
      STATE *p=&pred[0];
      int numPred=pred.size();
      //cout<<k<<" has "<<numPred<<" I preds: "<<pred<<endl;
      for(int s=0 ; s<numPred ; ++s, ++p) {
	STATE h=*p;
	if(hmm.getStateType(h)!=PHMM_INSERT) continue;
	double score=V(0,j-1,h)+hmm.getTransP(h,k);
	if(score>bestScore) { bestH=h; bestScore=score; }
      }
      /*
      for(int ih=0 ; ih<nI ; ++ih) {
	STATE h=I[ih];
	double score=V(0,j-1,h)+hmm.getTransP(h,k);
	if(score>bestScore) { bestH=h; bestScore=score; }
      }
      */
      V(0,j,k)=bestScore+Pe;
      T(0,j,k)=bestH;
    }
  
  // Now for the recursion:
  cout<<"recursion"<<endl;
  for(int i=1 ; i<=L1 ; ++i) {
    Array3D<double>::IndexedOnce<double> Vi=V[i], Vim1=V[i-1];
    Array3D<int>::IndexedOnce<int> Ti=T[i];
    for(int j=1 ; j<=L2 ; ++j) {
      if(bandingType==FIXED_WIDTH_BANDING && abs(i-j)>=bandwidth) continue;
      Array3D<double>::IndexedTwice<double> Vij=Vi[j],
	V_im1_jm1=Vim1[j-1], V_i_jm1=Vi[j-1], V_im1_j=Vim1[j];
      Array3D<int>::IndexedTwice<int> Tij=Ti[j];
      for(int ik=0 ; ik<nM ; ++ik) { // match states
	STATE k=M[ik], bestH=-2;
	double Pe=getEmitP(k,i-1,j-1);
	double bestScore=NEGATIVE_INFINITY;

	const Array1D<STATE> &pred=hmm.getPredecessors(k);
	STATE *p=&pred[0];
	int numPred=pred.size();
	for(int s=0 ; s<numPred ; ++s, ++p) {
	  STATE h=*p;
	  double score=V_im1_jm1[h]+hmm.getTransP(h,k);
	  if(score>bestScore) { bestH=h; bestScore=score; }
	}
	/*
	for(STATE h=0 ; h<numStates ; ++h) {
	  double score=V_im1_jm1[h]+hmm.getTransP(h,k);
	  if(score>bestScore) { bestH=h; bestScore=score; }
	}
	*/
	Vij[k]=bestScore+Pe;
	Tij[k]=bestH;
      }
      for(int ik=0 ; ik<nD ; ++ik) { // deletion states
	STATE k=D[ik], bestH=-2;
	double Pe=getEmitP(k,i-1,GAP);
	double bestScore=NEGATIVE_INFINITY;

	const Array1D<STATE> &pred=hmm.getPredecessors(k);
	STATE *p=&pred[0];
	int numPred=pred.size();
	for(int s=0 ; s<numPred ; ++s, ++p) {
	  STATE h=*p;
	  double score=V_im1_j[h]+hmm.getTransP(h,k);
	  if(score>bestScore) { bestH=h; bestScore=score; }
	}
	/*
	for(STATE h=0 ; h<numStates ; ++h) {
	  double score=V_im1_j[h]+hmm.getTransP(h,k);
	  if(score>bestScore) { bestH=h; bestScore=score; }
	}
	*/
	Vij[k]=bestScore+Pe;
	Tij[k]=bestH;
      }
      for(int ik=0 ; ik<nI ; ++ik) { // insertion states
	STATE k=I[ik], bestH=-2;
	double Pe=getEmitP(k,GAP,j-1);
	double bestScore=NEGATIVE_INFINITY;

	const Array1D<STATE> &pred=hmm.getPredecessors(k);
	STATE *p=&pred[0];
	int numPred=pred.size();
	for(int s=0 ; s<numPred ; ++s, ++p) {
	  STATE h=*p;
	  double score=V_i_jm1[h]+hmm.getTransP(h,k);
	  if(score>bestScore) { bestH=h; bestScore=score; }
	}
	/*
	for(STATE h=0 ; h<numStates ; ++h) {
	  double score=V_i_jm1[h]+hmm.getTransP(h,k);
	  if(score>bestScore) { bestH=h; bestScore=score; }
	}
	*/
	Vij[k]=bestScore+Pe;
	Tij[k]=bestH;
      }
    }
  }
}



StatePath *ProfileViterbi::getPath(double &bestScore)
{
  // Initialization:

  StatePath *path=new StatePath(&hmm);
  int bestH=-1;
  bestScore=NEGATIVE_INFINITY;
  Array3D<double>::IndexedTwice<double> V_L1_L2=V[L1][L2];
  for(STATE h=1 ; h<numStates ; ++h) {
    double score=V_L1_L2[h]+hmm.getTransP(h,0);
    if(score>bestScore) { bestH=h; bestScore=score; }
  }

  // Recursion:

  int i=L1, j=L2;
  STATE t=bestH;
  while(t>0) {
    path->push_back(t);
    STATE newT=T(i,j,t);
    switch(hmm.getStateType(t)) 
      {
      case PHMM_MATCH:  --i; --j; break;
      case PHMM_INSERT:      --j; break;
      case PHMM_DELETE: --i;      break;
      default:  throw "error";
      }
    t=newT;
  }

  // Termination:
  StatePath *revPath=path->getReverse();
  delete path;
  return revPath;
}



void ProfileViterbi::setupMasterCol()
{
  // First, populate the masterCol with a single column across all
  // tracks implied by the phylogeny:

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
  } visitor(masterCol,gap,nameToTrackId);
  phylogeny.postorderTraversal(visitor);

  // Now prepare the mapping arrays that map track ID's from the two
  // alignments to the masterCol aligment:

  setupMapping(S1,trackMap1);
  setupMapping(S2,trackMap2);
}



void ProfileViterbi::setupMapping(const MultSeqAlignment &A,
				  Array1D<int> &trackMap)
{
  int n=A.getNumTracks();
  trackMap.resize(n);
  for(int i=0 ; i<n ; ++i) {
    AlignmentSeq &track=A.getIthTrack(i);
    trackMap[i]=nameToTrackId[track.getName()];
  }
}



void ProfileViterbi::installColumn(int col,const MultSeqAlignment &A,
				   Array1D<int> &trackMap)
{
  int n=trackMap.size();
  for(int i=0 ; i<n ; ++i) {
    int ID=trackMap[i];
    masterCol.getIthTrack(ID).getSeq()[0]=A.getIthTrack(i).getSeq()[col];
  }
}



MultSeqAlignment *ProfileViterbi::decodeAlignment(const StatePath &path)
{
  // Create an empty alignment with appropriate tracks
  MultSeqAlignment &A=*new MultSeqAlignment(DnaDashAlphabet::global(),gap);
  int n1=S1.getNumTracks(), n2=S2.getNumTracks();
  Array1D<int> trackMap1(n1), trackMap2(n2);
  for(int i=0 ; i<n1 ; ++i) 
    trackMap1[i]=A.findOrCreateTrack(S1.getIthTrack(i).getName()).getID();
  for(int j=0 ; j<n2 ; ++j) 
    trackMap2[j]=A.findOrCreateTrack(S2.getIthTrack(j).getName()).getID();
  A.extendToLength(path.size());

  // Step along the path
  int i=0, j=0, L=path.size();
  A.enableAnnotation();
  Array1D<char> &anno=A.getAnnotationTrack();
  anno.resize(L);
  for(int k=0 ; k<L ; ++k) {
    FunctionalClass fc=hmm.getFunctionalClass(path[k],PARENT);
    anno[k]=fc.getLabel();
    switch(hmm.getStateType(path[k])) 
      {
      case PHMM_MATCH:
	copyColumn(S1,i,A,k,n1,trackMap1);
	copyColumn(S2,j,A,k,n2,trackMap2);
	++i;
	++j;
	break;
      case PHMM_INSERT:
	copyColumn(S1,i,A,k,n1,trackMap1,true);
	copyColumn(S2,j,A,k,n2,trackMap2);
	++j; 
	break;
      case PHMM_DELETE:
	copyColumn(S1,i,A,k,n1,trackMap1);
	copyColumn(S2,j,A,k,n2,trackMap2,true);
	++i;
	break;
      }
  }
  return &A;
}



void ProfileViterbi::copyColumn(const MultSeqAlignment &from,int fromCol,
				MultSeqAlignment &to,int toCol,
				int numTracks,const Array1D<int> &trackMap,
				bool useGaps)
{
  if(useGaps)
    for(int i=0 ; i<numTracks ; ++i)
      to.getIthTrack(trackMap[i])[toCol]=gap;
  else
    for(int i=0 ; i<numTracks ; ++i)
      to.getIthTrack(trackMap[i])[toCol]=from.getIthTrack(i)[fromCol];
}



double ProfileViterbi::getEmitP(STATE k,int col1,int col2)
{
  // Run Felsenstein's algorithm on the masterCol single-column alignment
  // (state k determines whether to run Fels on the combined clades or
  // just one clade -- i.e., MATCH versus INSERT/DELETE)
  switch(hmm.getStateType(k)) 
    {
    case PHMM_MATCH: {
      installColumn(col1,S1,trackMap1);
      installColumn(col2,S2,trackMap2);
      return F->logLikelihood(0,child1->getParent(),k);
    }
    case PHMM_INSERT: {
      installColumn(col2,S2,trackMap2);
      return F->logLikelihood(0,child2,k);
    }
    case PHMM_DELETE: {
      installColumn(col1,S1,trackMap1);
      return F->logLikelihood(0,child1,k);
    }
    }
  throw "error in ProfileViterbi";
}
