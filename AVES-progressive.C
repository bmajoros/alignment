/****************************************************************
 AVES-progressive.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <fstream>
#include "AVES.H"
#include "BOOM/Stack.H"
#include "SiblingHirschberg.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;



/****************************************************************
                         AVES::progressive()
 ****************************************************************/
void AVES::progressive()
{
  // Phase I : The UP-PASS
  LossyFelsenstein F1(pureDnaAlphabet,identityMap,numTaxa);
  //GainLossFelsenstein F1(pureDnaAlphabet,identityMap,numTaxa);//###DEBUGGING
  cout<<"uppass()"<<endl;
  uppass(F1);
  cout<<"dollo()"<<endl;
  dollo();
  if(baselineMode) {
    transducerTemplate=modelCompiler->getIthModel(1);
    if(!transducerTemplate) //throw "-B requires two models";
      transducerTemplate=modelCompiler->getIthModel(0);
    reinitBranches();
  }

 // Phase II : The DOWN-PASS
  //usePrescan=false; // ###
  cout<<"inferring root's functional parse"<<endl;
  GainLossFelsenstein F2(pureDnaAlphabet,identityMap,numTaxa);
  //LossyFelsenstein F2(pureDnaAlphabet,identityMap,numTaxa); //###
  inferRootFuncParse(F2);
  reconstructAlignment(true);
  if(baselineMode) {
    completeOrthology(fullAlignment);
    return;
  }
  fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',
			    60,true,true);
  cout<<"downpass()"<<endl;
  if(noPrescanOnDownpass) {usePrescan=false;cout<<"no prescans"<<endl;}
  downpass(F2);

  // Phase III : Iterative Refinement
  cout<<"iterativeRefinement()"<<endl;
  iterativeRefinement(F2);
}



/****************************************************************
                     AVES::completeOrthology()
 ****************************************************************/
void AVES::completeOrthology(MultSeqAlignment *A)
{
  String rootName=taxa[tree->getRoot()->getID()].getName();
  AlignmentSeq &rootTrack=A->getTrackByName(rootName);
  const String &rootAnno=rootTrack.getAnnoTrack();
  int numTracks=A->getNumTracks();
  for(int i=0 ; i<numTracks ; ++i) {
    AlignmentSeq &track=A->getIthTrack(i);
    track.getAnnoTrack()=rootAnno;
  }
  Vector<PhylogenyNode*> nodes;
  tree->gatherNodes(nodes,PREORDER);
  int N=nodes.size();
  FunctionalClass bg=FunctionalClass::getBackground();
  for(int i=0 ; i<N ; ++i) {
    PhylogenyNode *node=nodes[i];
    Taxon &taxon=static_cast<Taxon&>(*node->getDecoration());
    BranchAttributes *branch=taxon.getBranchToParent();
    if(!branch) continue;
    Taxon &parent=branch->getParentTaxon();
    IndexMap &up=branch->getUpMap();
    FunctionalParse &fp=taxon.getFunctionalParse();
    int L=taxon.getSeqLen();
    fp.resize(L);
    FunctionalParse &parentFP=parent.getFunctionalParse();
    for(int j=0 ; j<L ; ++j) {
      int parentPos=up[j];
      fp[j]=(parentPos==IndexMap::UNDEFINED) ? bg : parentFP[parentPos];
    }
  }
}



/****************************************************************
                         AVES::uppass()
 ****************************************************************/
void AVES::uppass(LinkFelsenstein &F)
{
  Vector<PhylogenyNode*> nodes;
  tree->gatherNodes(nodes,POSTORDER);
  int n=nodes.size();
  InternalNode *iNode;
  for(int i=0 ; i<n ; ++i) {
    PhylogenyNode *node=nodes[i];
    switch(node->getNodeType()) {
    case ROOT_NODE: 
      progressiveRoot(node);
      break;
    case LEAF_NODE: break;
    case INTERNAL_NODE:
      iNode=static_cast<InternalNode*>(node);
      align(iNode->getLeft(),iNode->getRight(),F);
      break;
    }
  }
}



/****************************************************************
                      AVES::progressiveRoot()
 ****************************************************************/
void AVES::progressiveRoot(PhylogenyNode *node)
{
  RootNode *root=static_cast<RootNode*>(node);
  Taxon *rootTaxon=static_cast<Taxon*>(root->getDecoration());
  BranchAttributes *branch=rootTaxon->getBranch(UNIQUE);
  Taxon *childTaxon=&branch->getChildTaxon();
  IndexMap &downMap=branch->getDownMap();
  IndexMap &upMap=branch->getUpMap();
  int childLen=childTaxon->getSeqLen();
  rootTaxon->setSeqLen(childLen);
  downMap.resize(childLen);
  upMap.resize(childLen);
  for(int i=0 ; i<childLen ; ++i) connect(downMap,i,upMap,i);
}



/****************************************************************
                           AVES::align()
 ****************************************************************/
void AVES::align(PhylogenyNode *leftNode,PhylogenyNode *rightNode,
		 LinkFelsenstein &F)
{
  cout<<"align()"<<endl;

  // Get the two taxa
  Taxon &left=taxa[leftNode->getID()], &right=taxa[rightNode->getID()];
  Taxon &parent=*left.getParent();
  cout<<"left="<<left.getName()<<" right="<<right.getName()
      <<" parent="<<parent.getName()<<endl;

  // Instantiate a new PairHMM for the combined branch lengths connecting
  // the two sibling taxa to be aligned
  BranchAttributes *leftBranch=left.getBranchToParent();
  BranchAttributes *rightBranch=right.getBranchToParent();
  double combinedBranchLengths=
    leftBranch->getLength()+rightBranch->getLength();
  BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
						   combinedBranchLengths,
						   NULL);

  // Perform alignment
  ViterbiInterface &viterbi=
    *(useHirschberg ?
      (ViterbiInterface*) new SiblingHirschberg(hmm,&left,&right,F,bandwidth,
						usePrescan,numThreads,
						contentSensor) :
      (ViterbiInterface*) new SiblingViterbi(hmm,&left,&right,F,
					     bandingType,bandwidth));
  LinkForward *forward;
  LinkBackward *backward;
  if(usePosteriors) {
    cout<<"<<usePosteriors"<<endl;
    delete posteriors;
    forward=new LinkForward(*hmm,left,right,identityMap,numTaxa,F);
    backward=new LinkBackward(*hmm,left,right,identityMap,numTaxa,F);
    posteriors=new Posteriors(*forward,*backward);
    viterbi.usePosteriors(posteriors);
  }
  cout<<"decode()"<<endl;
  StatePath &path=*viterbi.decode();
  left.getFunctionalParse()=FunctionalParse(path,PARENT);
  right.getFunctionalParse()=FunctionalParse(path,CHILD);
  cout<<"path score="<<path.getScore()<<endl;

  cout<<"install indel histories"<<endl;
  installSiblingHistories(path,parent,leftBranch,rightBranch,hmm);
  delete &path;
  delete hmm;
  if(usePosteriors) { delete forward; delete backward; }

  if(usePrescan) {
    cout<<"update prescan constraints"<<endl;
    updatePrescans(parent);
  }

  if(useContentSensor) {
    cout<<"reconstructing ancestral sequence for content sensing"<<endl;
    parent.getSeq().resize(parent.getSeqLen());
    ResidueAddress ra(&parent,0);
    LinkParsimony parsimony(ra,pureDnaAlphabet,gapSymbol,numTaxa);
    parsimony.runFullSeq();
    parent.getSeqStr()=parent.getSeq()(alphabet);
  }

  // Print out alignment so far
  /*
  cout<<"infer parses"<<endl;
  //###inferParentChildParse(F,parent,left);
  //###inferParentChildParse(F,parent,right);
  cout<<"print alignments"<<endl;
  AlignmentBuilder AB(parent.getNode(),alphabet,gapSymbol,numTaxa);
  MultSeqAlignment *A=AB.buildAlignment(true,true);
  A->printSlice(cout,0,A->getLength(),'+',60,true,true);
  delete A;
  */
}



/****************************************************************
                     AVES::updatePrescans()
 ****************************************************************/
void AVES::updatePrescans(Taxon &parent)
{
  BranchAttributes *leftBranch=parent.getBranch(LEFT);
  BranchAttributes *rightBranch=parent.getBranch(RIGHT);
  Taxon &left=leftBranch->getChildTaxon();
  Taxon &right=rightBranch->getChildTaxon();
  IndexMap &leftMap=leftBranch->getDownMap();
  IndexMap &rightMap=rightBranch->getDownMap();
  int parentLen=parent.getSeqLen();
  Array1D<FuncClassSet> &parentConstraints=parent.getFcConstraints();
  Array1D<FuncClassSet> &leftConstraints=left.getFcConstraints();
  Array1D<FuncClassSet> &rightConstraints=right.getFcConstraints();
  parentConstraints.resize(parent.getSeqLen());
  for(int pos=0 ; pos<parentLen ; ++pos) {
    int leftPos=leftMap[pos], rightPos=rightMap[pos];
    FuncClassSet &parentSet=parentConstraints[pos];
    if(leftPos==IndexMap::UNDEFINED) {
      //if(rightPos==IndexMap::UNDEFINED) INTERNAL_ERROR;//### DEBUGGING
      parentSet=rightConstraints[rightPos];
    }
    else if(rightPos==IndexMap::UNDEFINED)
      parentSet=leftConstraints[leftPos];
    else {
      parentSet.clear();
      leftConstraints[leftPos].unionWith(rightConstraints[rightPos],
					 parentSet);
    }
  }
}



/****************************************************************
                     AVES::rebuildPrescans()
 ****************************************************************/
void AVES::rebuildPrescans()
{
  Vector<PhylogenyNode*> nodes;
  tree->gatherNodes(nodes,POSTORDER);
  Vector<PhylogenyNode*>::iterator cur=nodes.begin(), end=nodes.end();
  for(; cur!=end ; ++cur) {
    PhylogenyNode *node=*cur;
    if(node->getNodeType()!=INTERNAL_NODE) continue;
    Taxon *taxon=static_cast<Taxon*>(node->getDecoration());
    updatePrescans(*taxon);
  }
}



/****************************************************************
                  AVES::installSiblingHistories()
 ****************************************************************/
void AVES::installSiblingHistories(StatePath &path,Taxon &parent,
				   BranchAttributes *leftBranch,
				   BranchAttributes *rightBranch,
				   BranchHMM *hmm)
{
  int pathLen=path.size();
  if(pathLen==0) INTERNAL_ERROR;
  parent.setSeqLen(pathLen);
  IndexMap &leftUpMap=leftBranch->getUpMap();
  IndexMap &leftDownMap=leftBranch->getDownMap();
  IndexMap &rightUpMap=rightBranch->getUpMap();
  IndexMap &rightDownMap=rightBranch->getDownMap();
  leftDownMap.resize(pathLen); rightDownMap.resize(pathLen);
  leftUpMap.resize(leftBranch->getChildTaxon().getSeqLen());
  rightUpMap.resize(rightBranch->getChildTaxon().getSeqLen());
  //if(rightBranch->getChildTaxon().getSeqLen()<1) INTERNAL_ERROR;//###DEBUGGING
  int i=0, j=0;
  for(int k=0 ; k<pathLen ; ++k) {
    STATE q=path[k];
    switch(hmm->getStateType(q)) {
    case PHMM_MATCH:
      connect(leftDownMap,k,leftUpMap,i);
      connect(rightDownMap,k,rightUpMap,j);
      ++i; ++j;
      break;
    case PHMM_INSERT:
      connect(rightDownMap,k,rightUpMap,j);
      leftDownMap[k]=IndexMap::UNDEFINED; 
      ++j;
      break;
    case PHMM_DELETE:
      connect(leftDownMap,k,leftUpMap,i);
      rightDownMap[k]=IndexMap::UNDEFINED;
      ++i;
      break;
    default: throw String("AVES::align() ")+hmm->getStateType(q)+" "+q;
    }
  }
}



/****************************************************************
                         AVES::connect()
 ****************************************************************/
void AVES::connect(IndexMap &parent,int parentIndex,IndexMap &child,
		    int childIndex)
{
  parent[parentIndex]=childIndex;
  child[childIndex]=parentIndex;
}



/****************************************************************
                           AVES::emit()
 ****************************************************************/
void AVES::emit(const String &outfile)
{
  // Emit the alignment
  fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',
			    60,true,true);
  fullAlignment->save(outfile,false);

  // Emit the annotation
  Array1D<Array1D<FunctionalElement>*> taxonElems(numTaxa);
  ofstream gff(gffFile.c_str());
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    FunctionalParse &phi=taxon.getFunctionalParse();
    Array1D<FunctionalElement> &elems=*phi.getElements();
    int n=elems.size();
    for(int i=0 ; i<n ; ++i) {
      FunctionalElement E=elems[i];
      String substrate=taxon.getName(), source="AVES";
      String type=E.getType().getName();
      int begin=E.getBegin(), end=E.getEnd(), score=0, phase=0;
      Strand strand=E.getStrand();
      if(type[0]=='-') {
	strand='-';
	type=type.substring(1);
      }
      --end; // ###
      gff<<substrate<<"\t"<<source<<"\t"<<type<<"\t"<<begin<<"\t"<<end
	 <<"\t"<<score<<"\t"<<strand<<"\t"<<phase<<endl;
    }
    taxonElems[i]=&elems;
    //delete &elems;
  }

  // Construct & emit the functional element histories
  if(historiesFile!="") {
    ofstream os(historiesFile.c_str());
    for(int i=0 ; i<numTaxa ; ++i) {
      Taxon &taxon=taxa[i];
      BranchAttributes *branch=taxon.getBranchToParent();
      if(!branch) {
	Array1D<FunctionalElement> &elems=*taxonElems[i];
	int numElems=elems.size();
	for(int j=0 ; j<numElems ; ++j) {
	  FunctionalElement *elem=&elems[j];
	  int b=elem->getBegin(), e=elem->getEnd();
	  os<<elem->getType().getName()<<"@("<<b<<","<<e<<","
	    <<taxon.getName()<<") is a root"<<endl;
	}
      }
      IndexMap &upMap=branch->getUpMap();
      int parentID=taxon.getParent()->getID();
      Taxon *parentTaxon=taxon.getParent();
      Array1D<FunctionalElement> &parentElems=*taxonElems[parentID];
      Array1D<FunctionalElement> &elems=*taxonElems[i];
      int numElems=elems.size(), numParentElems=parentElems.size();
      int parentJ=0;
      for(int j=0 ; j<numElems ; ++j) {
	FunctionalElement *elem=&elems[j];
	int b=elem->getBegin(), e=elem->getEnd();
	int parentBegin=upMap[b], parentEnd=upMap[e];
	bool isRoot=true;
	for(int k=parentJ ; k<numParentElems ; ++k) {
	  FunctionalElement *parentElem=&parentElems[k];
	  if(parentBegin<parentElem->getEnd() && 
	     parentElem->getBegin()<parentEnd) {
	    elem->setParent(parentElem);
	    parentJ=k+1;
	    os<<elem->getType().getName()<<"@("<<b<<","<<e<<","
	      <<taxon.getName()<<") descends from "
	      <<parentElem->getType().getName()<<"@("<<parentElem->getBegin()
	      <<","<<parentElem->getEnd()<<","<<parentTaxon->getName()<<")"
	      <<endl;
	    isRoot=false;
	    break;
	  }
	}
	if(isRoot)
	  os<<elem->getType().getName()<<"@("<<b<<","<<e<<","
	    <<taxon.getName()<<") is a root"<<endl;
      }
    }
  }

  // Clean up
  for(int i=0 ; i<numTaxa ; ++i) delete taxonElems[i];
}



/****************************************************************
                  AVES::installAnnotationTrack()
 ****************************************************************/
void AVES::installAnnotationTrack(MultSeqAlignment &A,
					 FunctionalParse &f)
{
  Array1D<char> &annotationTrack=A.getAnnotationTrack();
  int L=f.size();
  annotationTrack.resize(L);
  for(int i=0 ; i<L ; ++i) 
    annotationTrack[i]=f[i].getLabel();
  A.enableAnnotation();
}


/****************************************************************
                      AVES::gatherCladeMembers()
 ****************************************************************/
void AVES::gatherCladeMembers()
{
  struct TV : public TreeVisitor {
    Array1D<Taxon> &taxa;
    TV(Array1D<Taxon> &taxa) : taxa(taxa) {}
    void processNode(InternalNode &node) {
      int id=node.getID();
      BitSet &M=taxa[id].getCladeMembers();
      M=taxa[node.getLeft()->getID()].getCladeMembers();
      M+=taxa[node.getRight()->getID()].getCladeMembers();
      M.addMember(id);
    }
    void processNode(LeafNode &node) {
      int id=node.getID();
      BitSet &M=taxa[id].getCladeMembers();
      M.setSize(taxa.size());
      M.addMember(id);
    }
    void processNode(RootNode &node) {
      int id=node.getID();
      BitSet &M=taxa[id].getCladeMembers();
      M=taxa[node.getChild()->getID()].getCladeMembers();
      M.addMember(id);
    }
  } visitor(taxa);
  tree->postorderTraversal(visitor);
}



/****************************************************************
  This function is performed after the initial, greedy-progressive 
  phase, and uses parsimony to infer a single gap pattern for each 
  unobservable taxon.  Since the progressive phase is greedy anyway, 
  this should be OK.
 ****************************************************************/
void AVES::inferGapPatterns(const MultSeqAlignment &A)
{
  // First, construct an alignment on 3 symbols: *, -, and ?, representing
  // a nucleotide, gap, or unknown element
  const Alphabet &gapAlphabet=GapPatternAlphabet::global();
  BOOM::Symbol unknown=gapAlphabet.lookup('?');
  BOOM::Symbol gap=gapAlphabet.lookup('-');
  BOOM::Symbol residue=gapAlphabet.lookup('*');
  int L=A.getLength();
  MultSeqAlignment augmentedAlignment(gapAlphabet,gap);
  int numTaxa=taxa.size();
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    const String taxonName=taxon.getName();
    AlignmentSeq &newTrack=augmentedAlignment.findOrCreateTrack(taxonName);
    newTrack.extendToLength(L,unknown);
    switch(taxon.getNode()->getNodeType())
      {
      case ROOT_NODE:
	if(!A.trackExists(taxonName)) break;
	// fall through...
      case LEAF_NODE: {
	AlignmentSeq &oldTrack=A.getTrackByName(taxonName);
	for(int j=0 ; j<L ; ++j) 
	  newTrack[j]=(oldTrack[j]==gapSymbol ? gap : residue);
	break; }
      case INTERNAL_NODE: break;
      }
  }	

  // Next, compute the most parsimonious 3-symbol strings for the 
  // unobserved sequences
  SingleGainParsimony parsimony(*tree,gapAlphabet,augmentedAlignment,
  		    augmentedAlignment.getGapSymbols());
  //FitchParsimony parsimony(*tree,gapAlphabet,augmentedAlignment,
  //			    augmentedAlignment.getGapSymbols());
  parsimony.run();

  // Finally, convert the 3-symbol strings to GapPattern objects
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    GapPattern *gp=new GapPattern(augmentedAlignment.getIthTrack(i).getSeq());
    taxon.setGapPattern(gp);
  }
  augmentedAlignment.printSlice(cout,0,augmentedAlignment.getLength(),
				'+',60,true);
}



/****************************************************************
                         AVES::initBranches()
 ****************************************************************/
void AVES::initBranches()
{
  ostream *os=NULL;
  if(shouldDumpModel) {
    os=new ofstream(dumpFile.c_str());
    if(!os->good()) throw String("Error writing to file: ")+dumpFile;
  }
  tree->gatherBranches(branches);
  int n=branches.size();
  branchAttributes.resize(n);
  for(int i=0 ; i<n ; ++i) {
    PhylogenyBranch *branch=&branches[i];
    PhylogenyNode *parentNode=branch->getParent();
    PhylogenyNode *childNode=branch->getChild();
    if(os) (*os)<<"BranchHMM for branch "<<parentNode->getName()<<":"
		<<childNode->getName()<<endl;
    BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
						     branch->getLength(),os);
    branchAttributes[i]=BranchAttributes(hmm,branch);
    Taxon *parent=&taxa[parentNode->getID()];
    parentNode->getDecoration()=parent;
    switch(parentNode->getNodeType()) 
      {
      case ROOT_NODE: parent->getIthBranch(0)=&branchAttributes[i]; break;
      case LEAF_NODE: break;
      case INTERNAL_NODE: {
	InternalNode *in=static_cast<InternalNode*>(parentNode);
	WhichChild c=in->whichChildIsThis(childNode);
	parent->getIthBranch(c)=&branchAttributes[i];
	}
	break;
      }
  }
  if(shouldDumpModel) delete os;
}



void AVES::reinitBranches()
{
  int n=branches.size();
  for(int i=0 ; i<n ; ++i) {
    PhylogenyBranch *branch=&branches[i];
    BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
						     branch->getLength(),NULL);
    branchAttributes[i].changeHMM(hmm);
  }
}



/****************************************************************
                     AVES::loadSeedAlignment()
 ****************************************************************/
void AVES::loadSeedAlignment(const String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw String("Error opening file: ")+filename;
  MultiAlignment *maf=MultiAlignment::nextAlignmentFromMAF(is);
  fullAlignment=new MultSeqAlignment(*maf,alphabet,gapSymbol);
  inferGapPatterns(*fullAlignment);
  delete maf;
  is.close();
  taxa[tree->getRoot()->getID()].setCladeAlignment(fullAlignment);
}



/****************************************************************
  This function initializes the indel histories (the IndexMap's) 
  on the branches of the tree, based on the gap patterns which were 
  inferred using parsimony.  Since this is only done once, after 
  the initial, greedy progressive phase, the use of max-parsimony 
  gap patterns should be OK.  The reason it needs to be done at all 
  is because the progressive phase aligns sibling taxa, not parent-
  child pairs of taxa, so each pairwise alignment between sibling 
  taxa leaves some parent-child indel information undetermined.  
  During sampling this problem doesn't arise.
 ****************************************************************/
void AVES::reconstructHistories()
{
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &parent=taxa[i];
    int numBranches=parent.getNumBranches();
    for(int j=0 ; j<numBranches ; ++j) {
      BranchAttributes *branch=parent.getIthBranch(j);
      Taxon &child=branch->getChildTaxon();
      reconstructHistory(parent,child,branch->getDownMap(),
			 branch->getUpMap());
    }
  }
}



/****************************************************************
                   AVES::reconstructHistory()
 ****************************************************************/
void AVES::reconstructHistory(Taxon &parentTaxon,Taxon &childTaxon,
			      IndexMap &downMap,IndexMap &upMap)
{
  // First, resize the maps
  const GapPattern &parent=*parentTaxon.getGapPattern();
  const GapPattern &child=*childTaxon.getGapPattern();
  int parentL=parent.countResidues(), childL=child.countResidues();
  int L=parent.getLength();
  downMap.resize(parentL);
  upMap.resize(childL);
  parentTaxon.setSeqLen(parentL);
  parentTaxon.getSeq()=Sequence(alphabet.lookup('A'),parentL);
  childTaxon.setSeqLen(childL);

  // Now reconstruct the maps
  int parentI=0, childI=0;
  for(int i=0 ; i<L ; ++i) {
    GapPatternElement a=parent[i], b=child[i];
    if(a==GPE_RESIDUE)
      if(b==GPE_RESIDUE) {
	downMap[parentI]=childI;
	upMap[childI++]=parentI++;
      }
      else downMap[parentI++]=IndexMap::UNDEFINED;
    else if(b==GPE_RESIDUE) upMap[childI++]=IndexMap::UNDEFINED;
  }
}



/****************************************************************
                       reconstructHistory()
 ****************************************************************/
void AVES::reconstructHistory(StatePath &statePath,PairHMM &pairHMM,
			      IndexMap &upMap,IndexMap &downMap)
{
  int len=statePath.length();
  int x=0, y=0;
  for(int i=0 ; i<len ; ++i) {
    STATE q=statePath[i];
    switch(pairHMM.getStateType(q))
      {
      case PHMM_MATCH:  
	downMap[x]=y; 
	upMap[y]=x; 
	++x;
	++y;
	break;
      case PHMM_INSERT: 
	upMap[y]=IndexMap::UNDEFINED;
	++y;
	break;
      case PHMM_DELETE:
	downMap[x]=IndexMap::UNDEFINED;
	++x;
	break;
      default: throw "AVES::reconstructHistory() : unknown state type";
      }
  }
}



/****************************************************************
                         AVES::dollo()
 ****************************************************************/
void AVES::dollo()
{
  cout<<"reconstructing alignment"<<endl;
  reconstructAlignment();
  cout<<"inferring gap patterns"<<endl;
  inferGapPatterns(*fullAlignment);
  cout<<"reconstructHistories"<<endl;
  reconstructHistories();
  if(usePrescan) {
    cout<<"rebuilding prescan constraints"<<endl;
    rebuildPrescans();
  }
}



/****************************************************************
                         AVES::downpass()
 ****************************************************************/
void AVES::downpass(LinkFelsenstein &F)
{
  cout<<"downpass..."<<endl;
  Vector<PhylogenyNode*> nodes;
  tree->gatherNodes(nodes,PREORDER);
  int n=nodes.size();
  for(int i=0 ; i<n ; ++i) {
    PhylogenyNode *node=nodes[i];
    Taxon &taxon=*static_cast<Taxon*>(node->getDecoration());
    switch(node->getNodeType()) {
    case ROOT_NODE: {
      Taxon &child=taxon.getBranch(UNIQUE)->getChildTaxon();
      child.getFunctionalParse()=taxon.getFunctionalParse();
    }
      break;
    case LEAF_NODE:
      // nothing to do
      break;
    case INTERNAL_NODE:
      downpass(taxon,F,taxon.getBranch(LEFT));
      downpass(taxon,F,taxon.getBranch(RIGHT));
      break;
    }
  }
}



/****************************************************************
                         AVES::downpass()
 ****************************************************************/
void AVES::downpass(Taxon &parent,LinkFelsenstein &F,
		    BranchAttributes *branch)
{
  Taxon &child=branch->getChildTaxon();
  ViterbiInterface &V=
    *(useHirschberg ?
      (ViterbiInterface*) new Hirschberg(&parent,&child,F,bandwidth,
					 parent.getSeqLen(),
					 child.getSeqLen(),usePrescan,
					 numThreads) :
      (ViterbiInterface*) new LinkViterbi(&parent,&child,F,bandingType,
					  bandwidth));
  LinkForward *forward;
  LinkBackward *backward;
  if(usePosteriors) {
    cout<<"using posteriors"<<endl;
    delete posteriors;
    BranchHMM *hmm=branch->getHMM();
    forward=new LinkForward(*hmm,parent,child,alphabetMap,numTaxa,F);
    backward=new LinkBackward(*hmm,parent,child,alphabetMap,numTaxa,F);
    posteriors=new Posteriors(*forward,*backward);
    V.usePosteriors(posteriors);
  }
  cout<<"decode() : parent="<<parent.getName()<<" child="<<child.getName()<<endl;
  cout<<"parentLen="<<parent.getSeqLen()<<" childLen="<<child.getSeqLen()<<endl;
  ViterbiConstraint constraint=
    baselineMode ? VC_INDEL_AND_PARENT_PARSE : VC_PARENT_PARSE;
  StatePath *path=V.decode(constraint);
  if(usePosteriors) {delete forward; delete backward;}
  //###if(path->length()==0) throw "empty path in AVES::downpass()";
  cout<<"assign functional parse of length "<<path->length()<<endl;

  updateFcConstraints(*path,child);
  
  cout<<child.getName()<<" parse has length "<<child.getFunctionalParse().length()<<endl;
  cout<<"NEW path score: "<<V.scorePath(*path)<<" ("<<path->getScore()<<")"<<endl;
  cout<<"update indel history"<<endl;
  branch->setStatePath(path); // also updates the indel history

  reconstructAlignment(true);
  cout<<"NEW ALIGNMENT FOR "<<parent.getName()<<"->"<<child.getName()
      <<":"<<endl;
  fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',
			    60,true,true);
}



void AVES::updateFcConstraints(StatePath &path,Taxon &taxon)
{
  FunctionalParse fp(path,CHILD);
  taxon.getFunctionalParse()=fp;
  Array1D<FuncClassSet> &fcConstraints=taxon.getFcConstraints();
  FuncClassSet *fcs=&fcConstraints[0];
  Vector<FunctionalClass>::iterator cur=fp.begin();
  int L=fp.size();
  for(int i=L ; i ; --i, ++fcs, ++cur) {
    FuncClassSet &s=*fcs;
    s.clear();
    s.push_back(*cur);
  }
}



/****************************************************************
                     AVES::inferRootFuncParse()
 ****************************************************************/
void AVES::inferRootFuncParse(LinkFelsenstein &F)
{
  PhylogenyNode *rootNode=tree->getRoot();
  Taxon &root=*static_cast<Taxon*>(rootNode->getDecoration());
  Taxon &child=root.getIthChild(0);
  cout<<root.getName()<<"->"<<child.getName()<<endl;
  ViterbiInterface &V=
    *(useHirschberg ?
      (ViterbiInterface*) new Hirschberg(&root,&child,F,bandwidth,
					 root.getSeqLen(),
					 child.getSeqLen(),usePrescan,
					 numThreads) :
      (ViterbiInterface*) new LinkViterbi(&root,&child,F,bandingType,
					  bandwidth));
  StatePath *path=V.decode(VC_INDEL_HISTORY); // posteriors not necessary!
  root.getFunctionalParse()=FunctionalParse(*path,PARENT);
  child.getBranchToParent()->setStatePath(path);
}



/****************************************************************
                     AVES::inferParentChildParse()
 ****************************************************************/
void AVES::inferParentChildParse(LinkFelsenstein &F,Taxon &parent,
				 Taxon &child)
{
  ViterbiInterface &V=
    *(useHirschberg ?
      (ViterbiInterface*) new Hirschberg(&parent,&child,F,bandwidth,
					 parent.getSeqLen(),
					 child.getSeqLen(),usePrescan,
					 numThreads) :
      (ViterbiInterface*) new LinkViterbi(&parent,&child,F,bandingType,
					  bandwidth));
  StatePath *path=V.decode(VC_INDEL_HISTORY);
  child.getBranchToParent()->setStatePath(path);
}



/****************************************************************
                         AVES::iterativeRefinement()
 ****************************************************************/
void AVES::iterativeRefinement(LinkFelsenstein &F)
{
  int numBranches=branchAttributes.size();
  for(int i=0 ; i<refinementIterations ; ++i) {
    int j=RandomNumber(numBranches);
    BranchAttributes &branch=branchAttributes[j];
    refine(branch,F);
  }
}



void AVES::updateIndels(BranchAttributes &branch,Taxon &parent,
			Taxon &child,StatePath &path)
{
  BranchHMM *hmm=branch.getHMM();
  int pathLen=path.size();
  IndexMap &upMap=branch.getUpMap();
  IndexMap &downMap=branch.getDownMap();
  int i=0, j=0;
  for(int k=0 ; k<pathLen ; ++k) {
    STATE q=path[k];
    switch(hmm->getStateType(q)) {
    case PHMM_MATCH:
      connect(downMap,i,upMap,j);
      ++i; ++j;
      break;
    case PHMM_INSERT:
      upMap[j]=IndexMap::UNDEFINED; 
      ++j;
      break;
    case PHMM_DELETE:
      downMap[i]=IndexMap::UNDEFINED;
      ++i;
      break;
    default: throw "AVES::updateIndels()";
    }
  }
}



void AVES::refine(BranchAttributes &branch,LinkFelsenstein &F)
{
  // First, realign along this branch
  Taxon &parent=branch.getParentTaxon(), &child=branch.getChildTaxon();
  ViterbiInterface &viterbi=
    *(useHirschberg ?
      (ViterbiInterface*) new Hirschberg(&parent,&child,F,bandwidth,
					 parent.getSeqLen(),
					 child.getSeqLen(),usePrescan,
					 numThreads) :
      (ViterbiInterface*) new LinkViterbi(&parent,&child,F,
					  bandingType,bandwidth));
  LinkForward *forward;
  LinkBackward *backward;
  if(usePosteriors) {
    delete posteriors;
    BranchHMM *hmm=branch.getHMM();
    forward=new LinkForward(*hmm,parent,child,alphabetMap,numTaxa,F);
    backward=new LinkBackward(*hmm,parent,child,alphabetMap,numTaxa,F);
    posteriors=new Posteriors(*forward,*backward);
    viterbi.usePosteriors(posteriors);
  }
  StatePath *path=viterbi.decode(VC_UNCONSTRAINED);
  if(usePosteriors) {delete forward; delete backward;}
  parent.getFunctionalParse()=FunctionalParse(*path,PARENT);
  child.getFunctionalParse()=FunctionalParse(*path,CHILD);
  updateIndels(branch,parent,child,*path);

  // Now propagate changes in all directions radiating away from this
  // branch
  Stack<TaxonPair> S;
  Vector<Taxon*> V;
  parent.getNeighborsExcept(&child,V);
  int n=V.size();
  for(int i=0 ; i<n ; ++i) {
    if(&parent==V[i]) INTERNAL_ERROR;
    S.push(TaxonPair(&parent,V[i]));
  }
  V.clear();
  child.getNeighborsExcept(&parent,V);
  n=V.size();
  for(int i=0 ; i<n ; ++i) {
    if(&child==V[1]) INTERNAL_ERROR;
    S.push(TaxonPair(&child,V[i]));
  }
  while(!S.isEmpty()) {
    TaxonPair taxa=S.pop();
    propagate(taxa.first,taxa.second,F);
    V.clear();
    taxa.second->getNeighborsExcept(taxa.first,V);
    n=V.size();
    for(int i=0 ; i<n ; ++i) {
      if(taxa.second==V[i]) INTERNAL_ERROR;
      S.push(TaxonPair(taxa.second,V[i]));
    }
  }
}



void AVES::propagate(Taxon *from,Taxon *to,LinkFelsenstein &F)
{
  bool fromIsParent=to->getParent()==from;
  Taxon *parent, *child;
  ViterbiConstraint vc;
  if(fromIsParent) {
    parent=from;
    child=to;
    vc=VC_PARENT_PARSE;
  }
  else {
    parent=to;
    child=from;
    vc=VC_CHILD_PARSE;
  }
  BranchAttributes *branch=child->getBranchToParent();
  ViterbiInterface &V=
    *(useHirschberg ?
      (ViterbiInterface*) new Hirschberg(parent,child,F,bandwidth,
					 parent->getSeqLen(),
					 child->getSeqLen(),usePrescan,
					 numThreads) :
      (ViterbiInterface*) new LinkViterbi(parent,child,F,bandingType,
					  bandwidth));
  LinkForward *forward;
  LinkBackward *backward;
  if(usePosteriors) {
    delete posteriors;
    BranchHMM *hmm=branch->getHMM();
    forward=new LinkForward(*hmm,*parent,*child,alphabetMap,numTaxa,F);
    backward=new LinkBackward(*hmm,*parent,*child,alphabetMap,numTaxa,F);
    posteriors=new Posteriors(*forward,*backward);
    V.usePosteriors(posteriors);
  }
  StatePath *path=V.decode(vc);
  if(fromIsParent)
    child->getFunctionalParse()=FunctionalParse(*path,CHILD);
  else
    parent->getFunctionalParse()=FunctionalParse(*path,PARENT);
  updateIndels(*branch,*parent,*child,*path);
}



