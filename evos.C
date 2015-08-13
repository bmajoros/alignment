/****************************************************************
 evos.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "evos.H"
#include "BOOM/SummaryStats.H"
#include "BOOM/Map.H"
#include "BOOM/Time.H"
#include "BOOM/GSL/Random.H"
#include "GainLossEvent.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;

int main(int argc,char *argv[])
  {
    try
      {
	Evos app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(const BOOM::RootException &e)
      {
	cerr<<e.getMessage()<<endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Evos::Evos()
  : alphabetMap(DnaDashAlphabet::global(),PureDnaAlphabet::global()),
    functionalParse(NULL), alignmentBuilder(NULL), gainsLosses(4)
{
  gapSymbol=alphabet.lookup('-');
  modelCompiler=new ModelCompiler;
  modelCompiler->setGCthreshold(1000000);//10000000);
  templateInstantiator=modelCompiler->getInstantiator();
  gainsLosses.setAllTo(0);
}



int Evos::main(int argc,char *argv[])
{
  Time timer;
  timer.startCounting();

  // Process command line
  CommandLine cmd(argc,argv,"Ld:s:l:");
  if(cmd.numArgs()!=5)
    throw String(
"evos [opts] <in.lambda> <root-length> <out.fasta> <out.maf> <out.anno>\n\
     where: -L = generate leaf sequences only\n\ 
            -d <outfile> : dump model parameters to file\n\
            -s ### : seed randomizer with number\n\
            -l <X> : rescale branch lengths by factor X\n\
");
  String lambdaFile=cmd.arg(0);
  rootLength=cmd.arg(1).asInt();
  String fastaFile=cmd.arg(2);
  String mafFile=cmd.arg(3);
  String annoFile=cmd.arg(4);
  gff.open(annoFile.c_str());
  unsigned seed=
    cmd.option('s') ? cmd.optParm('s').asUnsigned() : unsigned(time(0));
  //SeedRandomizer(unsigned(seed));
  GSL::Random::seed(seed);
  cout<<"randomization seed: "<<seed<<endl;
  //Lambda::GarbageCollector::setThreshold(10); // ###
  leavesOnly=cmd.option('L');
  shouldDumpModel=cmd.option('d');
  if(shouldDumpModel) dumpFile=cmd.optParm('d');
  
  // Load HMM
  modelCompiler->parse(lambdaFile);
  if(modelCompiler->getNumModels()<1)
    throw "use (register-transducer ...) to register a transducer";
  transducerTemplate=modelCompiler->getIthModel(0);

  // Load phylogeny
  tree=modelCompiler->getGuideTree();
  tree->gatherNodes(phylogenyNodes);
  if(cmd.option('l')) 
    tree->scaleBranches(cmd.optParm('l').asFloat());

  // Link taxa to phylogeny nodes
  numTaxa=phylogenyNodes.size();
  taxa.resize(numTaxa);
  for(int i=0 ; i<numTaxa ; ++i) {
    PhylogenyNode *node=phylogenyNodes[i];
    int ID=node->getID();
    nameToTaxon[node->getName()]=ID;
    Taxon &taxon=taxa[ID];
    taxon.setNode(node);
    node->getDecoration()=&taxon;
  }
  alignmentBuilder=
    new AlignmentBuilder(tree->getRoot(),alphabet,gapSymbol,numTaxa);

  // Instantiate BranchHMM's and substitution matrices on tree
  initBranches();
  gatherCladeMembers();

  // Generate root sequence
  cout<<"generating root"<<endl;
  generateRootSeq();

  // Generate all other sequences
  cout<<"evolving"<<endl;
  evolve();

  // Emit results
  cout<<"writing output"<<Endl;
  emit(fastaFile,mafFile,annoFile);

  timer.stopCounting();
  cout<<"Total time: "<<timer.elapsedSeconds()<<" sec"<<endl;
  return 0;
}



/****************************************************************
                      gatherCladeMembers()
 ****************************************************************/
void Evos::gatherCladeMembers()
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
                         Evos::initBranches()
 ****************************************************************/
void Evos::initBranches()
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
                         Evos::generateRootSeq()
 ****************************************************************/
/*
void Evos::generateRootSeq()
{
  // First, sample a state path:
  int rootID=tree->getRoot()->getID();
  Taxon &root=taxa[rootID];
  BranchAttributes *branch=root.getIthBranch(0);
  Taxon &child=branch->getChildTaxon();
  BranchHMM *hmm=branch->getHMM();
  //hmm->normalizeTransProbs(); // ### BAD!  ALREADY IN LOG SPACE!
  RootBackward B(*hmm,root,child,alphabetMap,numTaxa,rootLength);
  RootSampler sampler(*hmm,B,root,child,branch);
  double pathScore;
  StatePath &path=*sampler.samplePath(rootLength,pathScore);

  // Extract a functional parse from the state path:
  root.getFunctionalParse()=FunctionalParse(path,PARENT);
  
  // Now use each state's nucleotide equilibrium frequencies to sample
  // actual residues for the root sequence:
  root.setSeqLen(rootLength);
  Sequence &seq=root.getSeq();
  seq.resize(rootLength);
  int pathLen=path.length(), pos=0;
  for(int i=0 ; i<pathLen ; ++i) {
    STATE q=path[i];
    if(hmm->getStateType(q)!=PHMM_INSERT) {
      const Array1D<double> &eq=hmm->getEqFreqs(q);
      BOOM::Symbol s=sampleResidue(eq);
      seq[pos]=s;
      ++pos;
    }
  }

  // =============
  BranchAttributes *branch=parent.getIthBranch(i);
  Taxon &child=branch->getChildTaxon();
  BranchHMM *hmm=branch->getHMM();
  StatePath *path=NULL;
  ChildBackward B(*hmm,parent,child,alphabetMap,numTaxa);
  ChildSampler sampler(*hmm,B,parent,child,branch);
  double pathScore;
  path=sampler.samplePath(pathScore);
  child.getFunctionalParse()=FunctionalParse(*path,CHILD);
  sampleChildSeq(*branch,parent,child,*path,*hmm);

}
*/

void Evos::generateRootSeq()
{
  // First, sample a state path:
  int rootID=tree->getRoot()->getID();
  Taxon &root=taxa[rootID];
  BranchAttributes *branch=root.getIthBranch(0);
  Taxon &child=branch->getChildTaxon();
  BranchHMM *hmm=branch->getHMM();
  //hmm->normalizeTransProbs(); // ### BAD!  ALREADY IN LOG SPACE!
  //cout<<"running backward"<<endl;
  RootBackward B(*hmm,root,child,alphabetMap,numTaxa,rootLength);
  //cout<<"allocating sampler"<<endl;
  RootSampler sampler(*hmm,B,root,child,branch);
  double pathScore;
  //cout<<"running sampler"<<endl;
  StatePath &path=*sampler.samplePath(rootLength,pathScore);
  //cout<<"sampler finished"<<endl;

  // Extract a functional parse from the state path:
  root.getFunctionalParse()=FunctionalParse(path,PARENT);
  
  // Now use each state's nucleotide equilibrium frequencies to sample
  // actual residues for the root sequence:
  root.setSeqLen(rootLength);
  Sequence &seq=root.getSeq();
  seq.resize(rootLength);
  int pathLen=path.length(), pos=0;
  for(int i=0 ; i<pathLen ; ++i) {
    STATE q=path[i];
    if(hmm->getStateType(q)!=PHMM_INSERT) {
      const Array1D<double> &eq=hmm->getEqFreqs(q);
      BOOM::Symbol s=sampleResidue(eq);
      seq[pos]=s;
      ++pos;
    }
  }
  //seq.printOn(cout,pureDnaAlphabet);cout<<endl;
}




/****************************************************************
                         Evos::sampleResidue()
 ****************************************************************/
BOOM::Symbol Evos::sampleResidue(const Array1D<double> &eqFreqs) 
{
  double r=GSL::Random::randomDouble(0.0,1.0), cumulative=0.0;
  BOOM::Symbol i;
  for(i=0 ; i<4 ; ++i) {
    cumulative+=eqFreqs[i];
    if(cumulative>=r) break;
  }
  if(i>3) i=3;
  return i;
}



/****************************************************************
                         Evos::evolve()
 ****************************************************************/
void Evos::evolve()
{
  int rootID=tree->getRoot()->getID();
  Taxon &root=taxa[rootID];
  evolve(root);
}



/****************************************************************
                         Evos::evolve()
 ****************************************************************/
void Evos::evolve(Taxon &parent)
{
  int numBranches=parent.getNumBranches();
  for(int i=0 ; i<numBranches ; ++i) {
    BranchAttributes *branch=parent.getIthBranch(i);
    Taxon &child=branch->getChildTaxon();
    BranchHMM *hmm=branch->getHMM();

    // Sample a state path
    StatePath *path=NULL; //branch->getStatePath(); // ###
    if(!path) {
      ChildBackward B(*hmm,parent,child,alphabetMap,numTaxa);
      ChildSampler sampler(*hmm,B,parent,child,branch);
      double pathScore;
      path=sampler.samplePath(pathScore);
      branch->setStatePath(path);
    }
    
    // Set the child's functional parse
    child.getFunctionalParse()=FunctionalParse(*path,CHILD);

    // Sample sequence for child taxon
    sampleChildSeq(*branch,parent,child,*path,*hmm);

    // Convert state path to explicit indel history
    IndexMap &upMap=branch->getUpMap(), &downMap=branch->getDownMap();
    int parentLen=parent.getSeq().getLength();
    int childLen=child.getSeq().getLength();
    parent.setSeqLen(parentLen);
    child.setSeqLen(childLen);
    upMap.resize(childLen);
    downMap.resize(parentLen);
    reconstructHistory(*path,*hmm,branch->getUpMap(),branch->getDownMap());

    // Tabulate gains and losses
    countGainLoss(*path,*hmm,parent,child,branch->getLength());

    // Recurse to child
    evolve(child);
  }
}



/****************************************************************
                   Evos::countGainLoss()
 ****************************************************************/
void Evos::countGainLoss(StatePath &path,BranchHMM &hmm,Taxon &parent,
			 Taxon &child,float branchLen)
{
  Array1D<int> localCounts(4);
  localCounts.setAllTo(0);
  Vector<GainLossEvent> &events=
    *GainLossEvent::getEvents(path,parent,child);
  //cout<<events.size()<<" events found"<<endl;
  Vector<GainLossEvent>::iterator cur=events.begin(), end=events.end();
  for(; cur!=end ; ++cur) {
    GainLossEvent &event=*cur;
    GainLossType glt=event.getType();
    ++gainsLosses[glt];
    ++localCounts[glt];
    FunctionalElementType elemType=event.getElementType();
    gff<<parent.getName()<<":"<<child.getName()<<"\tevos\t"
       <<toString(glt)<<"\t"<<event.getParentBegin()<<"\t"
       <<event.getChildBegin()
       <<"\t.\t"<<elemType.getStrand()<<"\t.\t/factor="
       <<elemType.getName().substitute("-","")<<endl;
  }
  delete &events;
  FunctionalParse parentParse(path,PARENT), childParse(path,CHILD);
  Array1D<FunctionalElement> *parentElements=parentParse.getElements();
  Array1D<FunctionalElement> *childElements=childParse.getElements();
  int numParentElem=parentElements->size();
  int numChildElem=childElements->size();
  delete parentElements;
  delete childElements;
  cout<<"Branch "<<parent.getName()<<"->"<<child.getName()<<" gains="
      <<localCounts[GLT_GAIN]<<" losses="<<localCounts[GLT_LOSS]
      <<" retentions="<<localCounts[GLT_RETENTION]<<" turnovers="
      <<localCounts[GLT_GAIN]+localCounts[GLT_LOSS]<<" branchlen="
      <<branchLen<<" parentN="<<numParentElem<<" childN="
      <<numChildElem<<endl;
}



void Evos::countGainLoss_OLD(StatePath &path,BranchHMM &hmm,Taxon &parent,
			 Taxon &child,float branchLen)
{
  Array1D<int> localCounts(4);
  localCounts.setAllTo(0);
  int len=path.length();
  GainLossType prevGLT=GLT_VOID;
  for(int i=0 ; i<len ; ++i) {
    STATE q=path[i];
    GainLossType glt=hmm.crossFunctionalType(q);
    FunctionalElementType fet=hmm.getFunctionalElementType(q);
    if(i==0 || fet!=hmm.getFunctionalElementType(path[i-1])) {
      //if(glt!=prevGLT && glt!=GLT_VOID) {
      ++gainsLosses[glt];
      ++localCounts[glt];

      /*
      String substrate=parent.getName()+":"+child.getName();
      String factor=fet.getName();
      String type=;
      gff<<substrate<<"\tevos\t"<<glt<<"\t"<<
      */
      /*
      String substrate=taxon.getName(), source="evos";
      string type=E.getType().getName();
      int begin=E.getBegin(), end=E.getEnd(), score=0, phase=0;
      Strand strand=E.getStrand();
      gff<<substrate<<"\t"<<glt<<"\t"<<type<<"\t"<<begin<<"\t"<<end
	 <<"\t"<<score<<"\t"<<strand<<"\t"<<phase<<endl;
       */
    }
    prevGLT=glt;
  }
  FunctionalParse parentParse(path,PARENT), childParse(path,CHILD);
  Array1D<FunctionalElement> *parentElements=parentParse.getElements();
  Array1D<FunctionalElement> *childElements=childParse.getElements();
  int numParentElem=parentElements->size();
  int numChildElem=childElements->size();
  delete parentElements;
  delete childElements;
  cout<<"Branch "<<parent.getName()<<"->"<<child.getName()<<" gains="
      <<localCounts[GLT_GAIN]<<" losses="<<localCounts[GLT_LOSS]
      <<" retentions="<<localCounts[GLT_RETENTION]<<" turnovers="
      <<localCounts[GLT_GAIN]+localCounts[GLT_LOSS]<<" branchlen="
      <<branchLen<<" parentN="<<numParentElem<<" childN="
      <<numChildElem<<endl;
}



/****************************************************************
                   Evos::sampleChildSeq()
 ****************************************************************/
void Evos::sampleChildSeq(BranchAttributes &branch,Taxon &parent,
			      Taxon &child,StatePath &path,BranchHMM &hmm)
{
  //hmm.initSubstGenerators();
  //cout<<child.getName()<<endl;
  //cout<<path<<endl;
  Sequence &parentSeq=parent.getSeq(), &childSeq=child.getSeq();
  int parentLen=parentSeq.getLength();
  int pathLen=path.length();
  int parentPos=0;
  SubstitutionMatrix *Pt;
  for(int i=0 ; i<pathLen ; ++i) {
    STATE q=path[i];
    switch(hmm.getStateType(q)) {
    case PHMM_MATCH:
      Pt=hmm.getSubstMatrix(q);
      Pt->setAlphabetMap(AlphabetIdentityMap(Pt->getAlphabet()));
      childSeq+=Pt->mutate(parentSeq[parentPos]);
      ++parentPos;
      break;
    case PHMM_INSERT:
      childSeq+=sampleResidue(hmm.getEqFreqs(q));
      break;
    case PHMM_DELETE: 
      ++parentPos;
      break;
    default: throw "Evos::sampleChildSeq() : unknown state type";
    }
  }
}



/****************************************************************
                   Evos::reconstructHistory()
 ****************************************************************/
void Evos::reconstructHistory(StatePath &statePath,PairHMM &pairHMM,
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
                    Evos::alignmentStats()
 ****************************************************************/
void Evos::alignmentStats(MultSeqAlignment &A,ostream &os)
{
  // First, collect stats related to just the alignment (no annotations)
  Vector<int> gapLengths;
  int residues=0, gaps=0, N=A.getNumTracks(), gapBegin=-1;
  BOOM::Symbol GAP=A.getGapSymbol();
  for(int i=0 ; i<N ; ++i) {
    AlignmentSeq &track=A.getIthTrack(i);
    const Sequence &seq=track.getSeq();
    const String &anno=track.getAnnoTrack();
    int L=seq.getLength();
    Vector<BOOM::Symbol>::iterator cur=seq.begin(), end=seq.end();
    for(int pos=0; cur!=end ; ++cur, ++pos) {
      if(*cur==GAP) {
	++gaps; 
	if(gapBegin<0) gapBegin=pos;
      }
      else {
	++residues;
	if(gapBegin>=0) {
	  gapLengths.push_back(pos-gapBegin);
	  gapBegin=-1;
	}
      }
    }
    if(gapBegin>=0) {
      gapLengths.push_back(L-gapBegin);
      gapBegin=-1;
    }
  }

  // Now collect stats related to annotations
  PhylogenyNode *root=tree->getRoot();
  String rootName=root->getName();
  int L=A.getLength(), totalSites=0, rootSites=0, siteColumns=0;
  Vector<float> overallCons, fgCons, bgCons;
  for(int pos=0 ; pos<L ; ++pos) {
    int numPresent=0, nonGap=0;
    bool foreground=false;
    Map<int,float> P;
    for(int i=0 ; i<4 ; ++i) P[i]=0.0;
    for(int i=0 ; i<N ; ++i) {
      AlignmentSeq &track=A.getIthTrack(i);
      BOOM::Symbol a=track.getSeq()[pos];
      if(a==GAP) continue;
      P[a]+=1.0;
      ++nonGap;
      bool isRoot=track.getName()==rootName;
      const String &anno=track.getAnnoTrack();
      char label=anno[pos];
      FunctionalClass fc=FunctionalClass::classFromLabel(label);
      if(!fc.isBackground()) {
	foreground=true;
	FunctionalClass prevFC=
	  pos>0 ? FunctionalClass::classFromLabel(anno[pos-1])
	  : FunctionalClass::getBackground();
	if(prevFC.isBackground()) {
	  ++numPresent;
	  if(isRoot) ++rootSites;
	}
      }
    }
    if(numPresent>0) {
      degreesOfOrthology.push_back(numPresent/float(N));
      totalSites+=numPresent;
      ++siteColumns;
    }
    double H=0.0;
    for(int i=0 ; i<4 ; ++i) {
      double p=P[i]/nonGap;
      if(p>0) H-=p*log2(p);
    }
    float cons=(log2(4)-H)/log2(4);
    overallCons.push_back(cons);
    if(foreground) fgCons.push_back(cons);
    else bgCons.push_back(cons);
  }

  // Gather some basic sequence stats
  int totalDNA=0;
  for(int i=0 ; i<numTaxa ; ++i) totalDNA+=taxa[i].getSeqLen();
  int rootLen=static_cast<Taxon*>(root->getDecoration())->getSeqLen();
   
  // Report results
  int totalPositions=gaps+residues;
  float gapDensity=int(1000*gaps/float(totalPositions)+5/9.0)/10.0;
  os<<gaps<<" gapped positions, out of "<<totalPositions
    <<" total positions ("<<gapDensity<<"%)"<<endl;
  SummaryStats stats(gapLengths);
  double aveLen=stats.getMean(), minLen=stats.getMin(), maxLen=stats.getMax(),
    sdLen=stats.getStdDev();
  int numGaps=stats.getN();
  os<<"ave gap length: "<<aveLen<<" +/- "<<sdLen<<" ("<<minLen<<"-"
    <<maxLen<<", N="<<numGaps<<")"<<endl;
  SummaryStats annoStats(degreesOfOrthology);
  double aveDOO=annoStats.getMean();
  os<<"ave DOO="<<aveDOO<<", totalSites="<<totalSites
    <<" rootSites="<<rootSites<<" siteColumns="<<siteColumns<<endl;
  double SPK=totalSites/float(totalDNA)*1000;
  double rootSPK=rootSites/float(rootLen)*1000;
  os<<"sites/kb="<<SPK<<" ("<<rootSPK<<" in root)"<<endl;
  SummaryStats consStats(overallCons);
  cout<<"overall conservation: "<<consStats.getMean()<<endl;
  SummaryStats fgConsStats(fgCons);
  cout<<"foreground conservation: "<<fgConsStats.getMean()<<endl;
  SummaryStats bgConsStats(bgCons);
  cout<<"background conservation: "<<bgConsStats.getMean()<<endl;
  cout<<"TOTAL GAINS: "<<gainsLosses[GLT_GAIN]<<endl;
  cout<<"TOTAL LOSSES: "<<gainsLosses[GLT_LOSS]<<endl;
  cout<<"TOTAL RETENTIONS: "<<gainsLosses[GLT_RETENTION]<<endl;
}




/****************************************************************
                        Evos::emit()
 ****************************************************************/
void Evos::emit(const String &outfile,const String &mafFile,
		    const String &annoFile)
{
  // Emit alignment
  Taxon *targetSpecies=NULL;
  for(int i=0 ; i<numTaxa ; ++i)
    if(taxa[i].isLeaf()) { targetSpecies=&taxa[i]; break; }
  //cout<<"building alignment"<<endl;
  MultSeqAlignment *A=alignmentBuilder->buildAlignment(true,false);
  //cout<<"saving alignment"<<endl;
  A->save(mafFile,false); // true);

  // Report some stats about the alignment to stdout
  alignmentStats(*A,cout);

  // Emit individual sequences
  FastaWriter fastaWriter;
  ofstream os(outfile.c_str());
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    if(taxon.isLeaf()) {
      Sequence &seq=taxon.getSeq();
      String *seqStr=seq.toString(pureDnaAlphabet);
      String defline=String(">")+taxon.getName();
      fastaWriter.addToFasta(defline,*seqStr,os);
      delete seqStr;
    }
  }

  // Emit annotations
  //ofstream gff(annoFile.c_str());
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    FunctionalParse &phi=taxon.getFunctionalParse();
    Array1D<FunctionalElement> &elems=*phi.getElements();
    int n=elems.size();
    //cout<<n<<" elements"<<endl;
    for(int i=0 ; i<n ; ++i) {
      FunctionalElement E=elems[i];
      String substrate=taxon.getName(), source="evos";
      string type=E.getType().getName().substitute("-","");
      int begin=E.getBegin(), end=E.getEnd(), score=0, phase=0;
      Strand strand=E.getStrand();
      gff<<substrate<<"\t"<<source<<"\t"<<type<<"\t"<<begin<<"\t"<<end
	 <<"\t"<<score<<"\t"<<strand<<"\t"<<phase<<endl;
    }
    delete &elems;
  }
}


