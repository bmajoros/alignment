/****************************************************************
 AVES-progressive.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <fstream>
#include "AVES.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;



/****************************************************************
                         AVES::progressive()
 ****************************************************************/
void AVES::progressive()
{
  struct TV : public TreeVisitor 
  {
    Array1D<Taxon> &taxa;
    Alphabet &alphabet;
    BOOM::Symbol gap;
    AVES &app;
    StatePath *statePath;
    TV(Array1D<Taxon>&t,Alphabet &a,BOOM::Symbol g,AVES &app) 
      : taxa(t), alphabet(a), gap(g), app(app) {}
    void processNode(InternalNode &node) 
    {
      MultSeqAlignment *A=
	app.align(node.getLeft(),node.getRight(),statePath);
      Taxon &taxon=app.taxa[node.getID()];
      taxon.setCladeAlignment(A);
      //if(!taxon.getParent())
      //taxon.getFunctionalParse()=FunctionalParse(*statePath,PARENT);//###WRONG
    }
    void processNode(LeafNode &node) 
    {
      Taxon &taxon=app.taxa[node.getID()];
      MultSeqAlignment *A=
	new MultSeqAlignment(taxon.getSeq(),taxon.getName(),alphabet,gap);
      A->setScore(0.0);
      taxon.setCladeAlignment(A);
    }
    void processNode(RootNode &node) 
    {
      Taxon &rootTaxon=app.taxa[node.getID()];
      Taxon &childTaxon=app.taxa[node.getChild()->getID()];
      MultSeqAlignment *A=childTaxon.getCladeAlignment();
      rootTaxon.setCladeAlignment(A);
    }
  } visitor(taxa,alphabet,gapSymbol,*this);
  tree->postorderTraversal(visitor);
  fullAlignment=taxa[tree->getRoot()->getID()].getCladeAlignment();
  functionalParse=new FunctionalParse(*visitor.statePath,PARENT);
}



/****************************************************************
                           AVES::emit()
 ****************************************************************/
void AVES::emit(const String &outfile)
{
  fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,true);
  fullAlignment->save(outfile);
}



/****************************************************************
                           AVES::align()
 ****************************************************************/
MultSeqAlignment *AVES::align(PhylogenyNode *left,PhylogenyNode *right,
				     StatePath *&path)
{
  // Get the two taxa and their subalignments
  Taxon &taxon1=taxa[left->getID()], &taxon2=taxa[right->getID()];
  MultSeqAlignment *A1=taxon1.getCladeAlignment();
  MultSeqAlignment *A2=taxon2.getCladeAlignment();
  int parentID=left->getParent()->getID();
  Taxon &parent=taxa[parentID];

  // Instantiate a new PairHMM for the combined branch lengths connecting
  // the two sibling taxa to be aligned
  double combinedBranchLengths=parent.getIthBranch(0)->getLength()+
    parent.getIthBranch(1)->getLength();
  BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
						   combinedBranchLengths,
						   NULL);

  // Perform alignment
  ProfileViterbi viterbi(*hmm,*A1,left,*A2,right,*tree,alphabetMap,
			 bandingType,bandwidth);
  path=viterbi.getPath(proposalScore);
  delete functionalParse;
  MultSeqAlignment *A=viterbi.decodeAlignment(*path);
  cout<<"alignment:"<<endl;
  A->printSlice(cout,0,A->getLength(),"+",60,true);
  return A;
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
                        AVES::decodeAlignment()
 ****************************************************************/
MultSeqAlignment *AVES::decodeAlignment(const StatePath &path,
					const AlignmentView &inside,
					const AlignmentView &outside,
					BranchHMM *hmm)
{
  // Create an empty alignment with appropriate tracks
  MultSeqAlignment &A=*new MultSeqAlignment(alphabet,gapSymbol);
  int n=fullAlignment->getNumTracks();
  for(int i=0 ; i<n ; ++i) 
    A.findOrCreateTrack(fullAlignment->getIthTrack(i).getName());
  int L=path.size();
  A.extendToLength(L);

  // Step along the path
  int i=0, j=0;
  for(int k=0 ; k<L ; ++k) {
    switch(hmm->getStateType(path[k])) 
      {
      case PHMM_MATCH:
	copyColumn(outside,i,A,k);
	copyColumn(inside,j,A,k);
	++i;
	++j;
	break;
      case PHMM_INSERT:
	copyColumn(inside,j,A,k);
	++j; 
	break;
      case PHMM_DELETE:
	copyColumn(outside,i,A,k);
	++i;
	break;
      default: throw "Error in AVES::decodeAlignment()";
      }
  }
  return &A;
}



/****************************************************************
                       AVES::copyColumn()
 ****************************************************************/
void AVES::copyColumn(const AlignmentView &view,int fromCol,
			     MultSeqAlignment &to,int toCol)
{
  const BitSet &taxonIDs=view.getTaxonIDs();
  int numTaxa=taxa.size();
  for(int i=0 ; i<numTaxa ; ++i)
    if(taxonIDs.isMember(i) && trackMap[i]>=0)
      to.getIthTrack(trackMap[i])[toCol]=
	fullAlignment->getIthTrack(trackMap[i])[view.mapColumn(fromCol)];
}



/****************************************************************
                         AVES::initTrackMap()
 ****************************************************************/
void AVES::initTrackMap()
{
  int n=tree->getNumNodes();
  trackMap.resize(n);
  trackMap.setAllTo(-1);
  for(int i=0 ; i<n ; ++i) {
    PhylogenyNode *node=taxa[i].getNode();
    if(node->getNodeType()==LEAF_NODE)
      trackMap[i]=
	fullAlignment->getTrackByName(node->getName()).getID();
  }
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
  MultSeqAlignment augmentedAlignment(gapAlphabet,gap);//unknown);
  //augmentedAlignment.extendToLength(L);
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
  During sampling this isn't an issue.
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
  This function deterministically infers a StatePath from the indel
  history (i.e., IndexMap's) of a branch.  It's only useful after the 
  initial greedy-progressive phase.
 ****************************************************************/
void AVES::reconstructStatePaths()
{
  // First, we need to reconstruct some (crude) functional parses
  // for all taxa, because LinkFelsenstein requires this in order
  // to assess substitution rates when computing emissions:

  int L=functionalParse->length();
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    GapPattern *gaps=taxon.getGapPattern();
    if(!gaps) throw "missing gap pattern";
    if(gaps->getLength()!=L) 
      throw "gap pattern length doesn't match functional parse length";
    int numResidues=gaps->countResidues();
    taxon.setSeqLen(numResidues);
    int residueIndex=0;
    FunctionalParse &taxonFP=taxon.getFunctionalParse();
    for(int j=0 ; j<L ; ++j)
      if((*gaps)[j]==GPE_RESIDUE)
	taxonFP[residueIndex++]=(*functionalParse)[j];
  }
  
  // Now we can reconstruct the actual state paths:

  struct TV : public TreeVisitor {
    void processNode(InternalNode &node) {
      reconstructStatePath(node,*node->getLeft());
      reconstructStatePath(node,*node->getRight());
    }
    void processNode(RootNode &node) {
      reconstructStatePath(node,*node.getChild());
    }
  } visitor;
  tree->preorderTraversal(visitor);

  /*
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &parent=taxa[i];
    int numBranches=parent.getNumBranches();
    for(int j=0 ; j<numBranches ; ++j) {
      BranchAttributes *branch=parent.getIthBranch(j);
      Taxon &child=branch->getChildTaxon();
      StatePath *path=reconstructStatePath(parent,child,branch->getHMM());
      branch->setStatePath(path);
    }
  } 
  */
}



/****************************************************************
                   AVES::reconstructStatePath()
 ****************************************************************/
StatePath *AVES::reconstructStatePath(PhylogenyNode &parentNode,
				      PhylogenyNode &childNode)
{
  throw "AVES::reconstructStatPath() needs to be re-implemented";

  Taxon &parent=taxa[parentNode->getID()], &child=taxa[childNode->getID()];
  if(parent.isRoot() && 


  LinkViterbi V(parent,child,bandingType,bandwidth);
  V.decode();


  /*
  StatePath *path=new StatePath(hmm);
  FunctionalParse &parentParse=parent.getFunctionalParse();
  FunctionalParse &childParse=child.getFunctionalParse()p;
  GapPattern &parentGaps=*parent.getGapPattern();
  GapPattern &childGaps=*child.getGapPattern();
  int L=parentGaps.getLength(), childPos=0;
  for(int i=0 ; i<L ; ++i) {
    GapPatternElement p=parentGaps[i], c=childGaps[i];
    FunctionalClass fc=fParse[i];
    PairHMMStateType stateType;
    if(p==GPE_RESIDUE)
      stateType=(c==GPE_RESIDUE ? PHMM_MATCH : PHMM_DELETE);
    else 
      stateType=(c==GPE_RESIDUE ? PHMM_INSERT : PHMM_SILENT);
    if(stateType!=PHMM_SILENT) {
      STATE q=hmm->getState(fc,stateType);
      path->push_back(q);
    }
  }
  return path;
  */
}



