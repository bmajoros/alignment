/****************************************************************
 marge.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "marge.H"
#include "BOOM/PowerMean.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/List.H"
#include "PostLeafBackward.H"

#define VERBOSE 1
void verbose(const String &x) {if(VERBOSE)cout<<x<<endl;}


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
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
    catch(const RootException &e)
      {
	cerr << "Exception: "<<e.getMessage()<<endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



ostream &operator<<(ostream &os,const EditType &e)
{
  switch(e)
    {
    case INSERT: os<<"INSERT"; break;
    case DELETE: os<<"DELETE"; break;
    }
  return os;
}



void EditEvent::printOn(ostream &os) const
{
  os<<type<<":"<<pos<<"->"<<mapTo;
}



ostream &operator<<(ostream &os,const EditEvent &e)
{
  e.printOn(os);
  return os;
}



Application::Application()
  : pathRegex("(.*\/)([^\/]+).post2"),
    filenameRegex("([^\/]+)-([^\/\\.]+)\\.post2"),
    transProbs(4,4)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  MPI *mpi=new MPI(argc,&argv);
  slaveDriver=new MpiSlaveDriver(*mpi);
  if(slaveDriver->amItheMaster()) return master(argc,argv);
  else return slave(argc,argv);
}




void Application::processCmdLine(int argc,char *argv[]) 
{
  // Process command line
  CommandLine cmd(argc,argv,"b:Gr:c:ps:i:");
  if(cmd.numArgs()!=5)
    throw String("\n\
MARGE = Multiple Alignment via Random Graph Estimation\n\
\n\
usage: marge <model.lambda> <library> <*.fasta> <#samples> <outfile>\n\
            -b <strategy> = band using strategy:\n\
                f### = fixed-width banding of width 2*###\n\
            -G = don't perform garbage collection\n\
            -r N = use randomization seed N\n\
            -c N = use pseudocount of N (default=0.01)\n\
            -p = print each alignment after sampling it\n\
            -s N = skip first N samples (default=0)\n\
            -i c = indel coef\n\
");
  String lambdaFile=cmd.arg(0);
  String libraryDir=cmd.arg(1);
  String fastaFile=cmd.arg(2);
  numSamples=cmd.arg(3).asInt();
  outfile=cmd.arg(4);
  if(cmd.option('r')) seed=cmd.optParm('r').asUnsigned();
  else seed=GetRandomSeed();
  logPseudoCount=log(cmd.option('c') ? cmd.optParm('c').asFloat() : 0.01);
  SeedRandomizer(seed);
  indelCoef=cmd.option('i') ? cmd.optParm('i').asFloat() : 1.0;

  if(cmd.option('b')) {
    String strategy=cmd.optParm('b');
    if(strategy[0]=='f' && strategy.length()>1) {
      bandwidth=strategy.substring(1,strategy.length()-1).asInt();
      bandingType=FIXED_WIDTH_BANDING;
    }
    else throw "invalid banding strategy specified";
  }
  else bandingType=NO_BANDING;
  skipSamples=cmd.option('s') ? cmd.optParm('s').asInt() : 0;
  wantPrint=cmd.option('p');
  numThreads=cmd.option('t') ? cmd.optParm('t').asInt()-1 : 0;
  cout<<(numThreads+1)<<" threads"<<endl;
  disableGC=cmd.option('G');
  if(disableGC) modelCompiler->setGCthreshold(LARGEST_INTEGER);
  shouldSeed=false;
  baselineMode=false;
  usePosteriors=false;
  shouldHillClimb=false;
  wantFuncDollo=false;
  explicitHistories=true;
  shouldDumpModel=false;
  useHirschberg=true;
  useContentSensor=false;
  usePrescan=true;
  noPrescanOnDownpass=true;
  
  // Load HMM
  cout<<"executing "<<lambdaFile<<"..."<<endl;
  LambdaAPI &lambda=modelCompiler->getLambda();
  lambda.getGC().setSilence(true);
  modelCompiler->parse(lambdaFile);
  lambda.checkGarbageLevel();
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  transducerTemplate=modelCompiler->instantiate(ctor);
  cout<<transducerTemplate->getNumStates()<<" states"<<endl;
  hmm=templateInstantiator->instantiate(transducerTemplate,1.0,NULL);
  numStates=hmm->getNumStates();
  getStateTypes(stateTypes,dI,dJ);
  dot=alphabet.lookup('.');

  // Get phylogeny
  cout<<"loading tree"<<endl;
  tree=modelCompiler->getGuideTree(transducerTemplate);
  tree->gatherNodes(phylogenyNodes);
  tree->constructBranches();
  if(cmd.option('l'))
    tree->scaleBranches(cmd.optParm('l').asFloat());
  tree->computePairwiseDistances();

  // Link taxa to phylogeny nodes
  cout<<"linking taxa to phylogeny nodes"<<endl;
  numTaxa=phylogenyNodes.size();
  taxa.resize(numTaxa);
  linkSeqs.resize(numTaxa);
  for(int i=0 ; i<numTaxa ; ++i) {
    PhylogenyNode *node=phylogenyNodes[i];
    int ID=node->getID();
    nameToTaxon[node->getName()]=ID;
    Taxon &taxon=taxa[ID];
    taxon.setNode(node);
    node->getDecoration()=&taxon;
  }
  bool showOrphans=true; // ###
  alignmentBuilder=
    new AlignmentBuilder(tree->getRoot(),DnaDashDotAlphabet::global(),
			 DnaDashDotAlphabet::global().lookup('-'),
			 numTaxa,showOrphans);
  /*
    new AlignmentBuilder(tree->getRoot(),alphabet,
			 alphabet.lookup('-'),numTaxa,showOrphans);
  */
  gatherCladeMembers(); // ### ?

  // Instantiate BranchHMM's and substitution matrices on tree
  cout<<"instantiating HMM's"<<endl;
  initBranches();

  // Load posterior matrices
  loadMatrices(libraryDir);

  if(shouldSeed) {
    cout<<"loading seed alignment"<<endl;
    loadSeedAlignment(cmd.optParm('i'));
    // ### other stuff needs to be done here...otherwise probably segfault...

    cout<<"loading sequences"<<endl;
    FastaReader fastaReader(fastaFile,alphabet);
    String def, seq, name, junk;
    while(fastaReader.nextSequence(def,seq)) {
      FastaReader::parseDefline(def,name,junk);
      if(!nameToTaxon.isDefined(name)) continue;
      Sequence S(seq,alphabet);
      if(S.getLength()==0) S.append(alphabet.lookup('A'));
      S.replaceAll(INVALID_SYMBOL,gapSymbol);
      S.useOnlyTheseSymbols(isNucleotide);
      int ID=nameToTaxon[name];
      Taxon &taxon=taxa[ID];
      taxon.getSeq()=S;
      taxon.setSeqLen(S.getLength());
      taxon.getSeqStr()=S(alphabet);
    }
    
    cout<<"reconstructHistories"<<endl;
    reconstructHistories();
  }
  else 
    {
      // Load sequences
      cout<<"loading sequences"<<endl;
      FastaReader fastaReader(fastaFile,alphabet);
      String def, seq, name, junk;
      while(fastaReader.nextSequence(def,seq)) {
	FastaReader::parseDefline(def,name,junk);
	if(!nameToTaxon.isDefined(name)) continue;
	Sequence S(seq,alphabet);
	if(S.getLength()==0) S.append(alphabet.lookup('A'));
	S.replaceAll(INVALID_SYMBOL,gapSymbol);
	S.useOnlyTheseSymbols(isNucleotide);
	int ID=nameToTaxon[name];
	Taxon &taxon=taxa[ID];
	taxon.getSeq()=S;
	taxon.setSeqLen(S.getLength());
	taxon.getSeqStr()=S(alphabet);
      }
      initLinkSeqs();
      
      for(int i=0 ; i<numTaxa ; ++i) {
	Taxon &taxon=taxa[i];
	Array1D<FuncClassSet> &fcConstraints=taxon.getFcConstraints();
	fcConstraints.resize(taxon.getSeqLen());
      }

      // Perform alignment
      cout<<"performing alignment"<<endl;
      progressive();
      MemoryProfiler::report("Total memory used:");
    }
  
  for(int i=0;i<numTaxa;++i) taxa[i].getFunctionalParse().clear();
  
  cout<<"reconstructing initial alignment"<<endl;
  reconstructAlignment(true);
  if(slaveDriver->amItheMaster()) {
    cout<<"initial alignment:"<<endl;
    fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,true);
  }	
}



void Application::initLinkSeqs()
{
  LinkedResidue r;
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &tax=taxa[i];
    if(!tax.isLeaf()) continue;
    LinkSeq &linkSeq=linkSeqs[i];
    int L=tax.getSeqLen();
    for(int j=0 ; j<L ; ++j) linkSeq.append(r);
    Sequence &seq=tax.getSeq();
    LinkSeq::iterator cur=linkSeq.begin(), end=linkSeq.end();
    for(int j=0 ; j<L ; ++j, ++cur) {
      LinkedResidue &r=*cur;
      r.s=seq[j];
    }	
  }
}



int Application::master(int argc,char *argv[])
{
  timer.startCounting();
  processCmdLine(argc,argv);
  cout.flush();

  int numSlaves=slaveDriver->getNumSlaves();
  int samplesPerSlave=numSamples/numSlaves;
  int soFar=0;
  for(int i=0 ; i<numSlaves ; ++i) {
    MpiVariableMsg &msg=*new MpiVariableMsg(SAMPLE);
    int slaveSamples=samplesPerSlave;
    if(i==numSlaves-1) slaveSamples=numSamples-soFar;
    msg<<slaveSamples;
    msg.close();
    slaveDriver->addWork(&msg);
    soFar+=slaveSamples;
  }

  Vector<MpiFixedMsg*> results;
  slaveDriver->waitForResults(results);
  Vector<MpiFixedMsg*>::iterator cur=results.begin(), end=results.end();
  for(; cur!=end ; ++cur) delete *cur;
  slaveDriver->terminateSlaves();

  // Consolidate alignment files from slaves
  MPI &mpi=slaveDriver->getMPI();
  File out(outfile,"w");
  int n=numSamples-numSlaves*skipSamples;
  out.write(n);
  for(int i=0 ; i<numSlaves ; ++i) {
    int slaveID=mpi.getIthSlaveID(i);
    String slaveFile=outfile+"-"+slaveID;
    File in(slaveFile);
    long size=in.getSize();
    for(int j=0 ; j<size ; ++j) {
      unsigned char c=in.readByte();
      out.write(c);
    }
  }

  cout<<"master done"<<endl;
  timer.stopCounting();
  cout<<"MASTER ELAPSED TIME: "<<timer.elapsedTime()<<endl;
  return 0;
}



int Application::slave(int argc,char *argv[])
{
  processCmdLine(argc,argv);
  SeedRandomizer(seed+slaveDriver->getMPI().getProcessID());
  
  // Start the timer
  timer.startCounting();

  // Perform gradient ascent
  slave_serveTheMaster();

  // Clean up
  timer.stopCounting();
  MemoryProfiler::report("Total memory used by slave:");
  return 0;
}



void Application::slave_serveTheMaster()
{
  while(1) {
    MpiFixedMsg *msg=slaveDriver->acceptWork();
    switch(msg->getTag())
      {
      case SAMPLE:
	slave_sample(*msg);
	break;
      case TERMINATE_SLAVE:
	return;
      }
    delete msg;
  }
}



void Application::slave_sample(MpiFixedMsg &msg)
{
  msg>>numSamples;

  samplingLoop();

  MpiVariableMsg *newMsg=new MpiFixedMsg(DONE_SAMPLING);
  slaveDriver->replyToMaster(newMsg);
}



/****************************************************************
                     Application::loadMatrices()
 ****************************************************************/
void Application::loadMatrices(const String &inDir)
{  
  matrices.resize(numTaxa,numTaxa);
  Array1D<int> lengths(numTaxa);
  Vector<String> files;
  if(!File::getFileList(inDir,files)) 
    throw String("Can't get dir listing for ")+inDir;
  Vector<String>::iterator cur=files.begin(), end=files.end();
  for(; cur!=end ; ++cur) {
    const String &filename=*cur;
    String baseName, path;
    if(pathRegex.match(filename)) {
      path=pathRegex[1];
      baseName=pathRegex[2];
    }
    else baseName=filename;
    if(!filenameRegex.match(baseName)) continue;
    String firstSpecies=filenameRegex[1], secondSpecies=filenameRegex[2];
    if(!nameToTaxon.isDefined(firstSpecies) ||
       !nameToTaxon.isDefined(secondSpecies)) continue;
    int tax1=nameToTaxon[firstSpecies], tax2=nameToTaxon[secondSpecies];
    SparseMatrix3D *M=SparseMatrix3D::loadBinary(inDir+"/"+filename);
    matrices[tax1][tax2]=M;
    matrices[tax2][tax1]=M->transpose();
    lengths[tax1]=M->getFirstDim();
    lengths[tax2]=M->getSecondDim();
  }
  for(int tax=0 ; tax<numTaxa ; ++tax) {
    if(!taxa[tax].isLeaf()) continue;
    int L=lengths[tax];//taxa[tax].getSeqLen();
    SparseMatrix3D &M=*new SparseMatrix3D(L+1,L+1,3);
    matrices[tax][tax]=&M;
    for(int i=1 ; i<=L ; ++i) 
      for(int j=1 ; j<=L ; ++j)
	M.addEntry(i,j,PHMM_MATCH,LOG_1);
  }
}



/****************************************************************
                     Application::progressive()
 ****************************************************************/
void Application::progressive()
{
  // Phase I : The UP-PASS
  LossyFelsenstein F1(pureDnaAlphabet,identityMap,numTaxa);
  cout<<"uppass()"<<endl;
  uppass(F1);
  cout<<"dollo()"<<endl;
  dollo();
}



/****************************************************************
                     Application::samplingLoop()
 ****************************************************************/
void Application::samplingLoop()
{
  outfile=outfile+"-"+slaveDriver->getMPI().getProcessID();
  File f(outfile,"w");
  ResidueOrthologyGraph rog(tree,taxa);
  for(int i=0 ; i<numSamples ; ++i) {
    cout<<"sample #"<<i<<endl;
    sample(f,rog,i);
    /*
    if(i>=skipSamples) {
      rog.saveConnections(f);
      if(wantPrint) {
	reconstructAlignment(true);
	fullAlignment->printSlice(cout,0,fullAlignment->getLength(),
				  '+',60,true);
      }
    }
    */
  }
}



/****************************************************************
                     Application::getStateTypes()
 ****************************************************************/
void Application::getStateTypes(Array1D<PHMM_StateType> &stateTypes,
				Array1D<int> &dI,Array1D<int> &dJ)
{
  int numStates=hmm->getNumStates();
  stateTypes.resize(numStates); dI.resize(numStates); dJ.resize(numStates);
  for(int i=0 ; i<numStates ; ++i) {
    STATE t=stateTypes[i]=hmm->getStateType(i);
    switch(t) {
    case PHMM_MATCH:  dI[i]=1; dJ[i]=1; matchState=i;  break;
    case PHMM_INSERT: dI[i]=0; dJ[i]=1; insertState=i; break;
    case PHMM_DELETE: dI[i]=1; dJ[i]=0; deleteState=i; break;
    }
  }
}



/****************************************************************
                   Application::initBranches()
 ****************************************************************/
void Application::initBranches()
{
  //tree->gatherBranches(branches);
  Vector<PhylogenyBranch*> tmp;
  tree->collectBranches(tmp);
  int n=tmp.size();
  branches.resize(n);
  for(int i=0 ; i<n ; ++i) branches[i]=*tmp[i];
  branchAttributes.resize(n);
  for(int i=0 ; i<n ; ++i) {
    PhylogenyBranch *branch=&branches[i];
    PhylogenyNode *parentNode=branch->getParent();
    PhylogenyNode *childNode=branch->getChild();
    BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
						     branch->getLength(),NULL);
    branchAttributes[i]=BranchAttributes(hmm,branch);
    branch->getDecoration()=&branchAttributes[i];
    tmp[i]->getDecoration()=&branchAttributes[i];
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
}



/****************************************************************
                   Application::classifyState()
 ****************************************************************/
PHMM_StateType Application::classifyState(IndexMap &upMap,IndexMap &downMap,
					  int parentPos,int childPos)
{
  if(childPos>=upMap.size()) return PHMM_DELETE;
  if(parentPos>=downMap.size()) return PHMM_INSERT;
  int upPos=upMap[childPos];
  //  childPos<upMap.size() ? upMap[childPos] : IndexMap::UNDEFINED;
  int downPos=downMap[parentPos];
    //  parentPos<downMap.size() ? downMap[parentPos] : IndexMap::UNDEFINED;
  int numDefined=0;
  if(upPos!=IndexMap::UNDEFINED) ++numDefined;
  if(downPos!=IndexMap::UNDEFINED) ++numDefined;
  switch(numDefined)
    {
    case 0: 
      return Random0to1()<0.5 ? PHMM_INSERT : PHMM_DELETE;
    case 1: return upPos==IndexMap::UNDEFINED ? PHMM_INSERT : PHMM_DELETE;
    case 2: return PHMM_MATCH;
    }
}



/****************************************************************
                     Application::sample()
 ****************************************************************/
void Application::sample(File &f,ResidueOrthologyGraph &rog,int sampleNum)
{
  for(int i=0 ; i<linkSeqs.size() ; ++i)
    if(taxa[i].isLeaf()) linkSeqs[i].removeLinks();
    else {
      linkSeqs[i].clear();
      taxa[i].setSeqLen(0);
    }
  int numNodes=tree->getNumNodes();
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &tax=taxa[i];
    if(!tax.isLeaf()) tax.setSeqLen(0);
  }
  int numLeaves=0;
  for(int i=0 ; i<numNodes ; ++i) if(taxa[i].isLeaf()) ++numLeaves;
  int numPairs=numLeaves-1; // num pairs in a spanning tree
  Array2D<bool> equivClasses(numNodes,numNodes);
  equivClasses.setAllTo(false);
  //BitSet nodesUsed(numNodes);
  List<PairwiseAlignment*> alignments;
  for(int i=0 ; i<numPairs ; ++i) {
    //CollapsedOrthologyMatrix::compute(rog);
    int tax1, tax2;
    samplePair(tax1,tax2,numNodes,equivClasses);
    equivClasses[tax1][tax2]=true;
    for(int j=tax1+1 ; j<numNodes ; ++j)
      if(equivClasses[tax1][j]) equivClasses[tax2][j]=true;
    for(int j=tax2+1 ; j<numNodes ; ++j)
      if(equivClasses[tax2][j]) equivClasses[tax1][j]=true;
    Taxon &leaf1=taxa[tax1], &leaf2=taxa[tax2];
    IndexMap fwdMap, backMap;
    PairwiseAlignment *aln=alignLeaves(tax1,tax2);
    alignments.append(aln);
    /*
    if(wantPrint) {
      reconstructAlignment(true); //false);
      fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,
				true);
    }
    */
  }

  // Visit the edges in a manner that grows a connected component, 
  // incrementally inferring the ancestral relations along the way
  BitSet touched(numTaxa);
  bool first=true;
  while(!alignments.isEmpty()) {
    //cout<<"alignments size="<<alignments.size()<<endl;
    //cout<<"touched="<<touched<<endl;
    List<PairwiseAlignment*>::iterator cur=alignments.begin(), 
      end=alignments.end();
    for(; cur!=end ; ) {
      PairwiseAlignment *aln=*cur;
      if(first || touched.isMember(aln->tax1) || 
	 touched.isMember(aln->tax2)) {
	if(!first && !touched.isMember(aln->tax1)) aln->invert();
	updateAncestors(*aln);
	touched.addMember(aln->tax1);
	touched.addMember(aln->tax2);
	List<PairwiseAlignment*>::iterator victim=cur;
	++cur;
	alignments.erase(victim);
	delete aln;
	first=false;
      }
    }
  }
}



/****************************************************************
                     Application::alignLeaves()
 ****************************************************************/
PairwiseAlignment *Application::alignLeaves(int tax1,int tax2)
{
  Taxon &leaf1=taxa[tax1], &leaf2=taxa[tax2];
  PairwiseAlignment *aln=new PairwiseAlignment;
  aln->tax1=tax1; aln->tax2=tax2;
  PostLeafBackward B(*hmm,leaf1,leaf2,*matrices[tax1][tax2],logPseudoCount,
		     indelCoef);
  int parentL=leaf1.getSeqLen(), childL=leaf2.getSeqLen();
  StatePath path(hmm);
  int i=0, j=0; // 1-based
  STATE prevQ=0;
  while(i<parentL || j<childL) {
    Vector<STATE> states;
    Vector<float> probs;
    for(int q=1 ; q<numStates ; ++q) {
      PHMM_StateType stateType=stateTypes[q];
      int newI=i+dI[q], newJ=j+dJ[q];
      if(newI>parentL || newJ>childL) continue;
      float cellP=B.getEmitP(q,newI,newJ);
      cellP=safeAdd(cellP,B(newI,newJ,q),-B(i,j,prevQ));
      if(isFinite(cellP)) {
	probs.push_back(exp(cellP));
	states.push_back(q);
      }
    }
    if(probs.size()==0) INTERNAL_ERROR;
    float sum=0;
    Vector<float>::iterator cur=probs.begin(), end=probs.end();
    for(; cur!=end ; ++cur) sum+=*cur;
    RouletteWheel wheel;
    cur=probs.begin();
    for(int i=0 ; cur!=end ; ++cur, ++i) {
      wheel.addSector(probs[i]/sum);
    }
    wheel.doneAddingSectors();
    int x=wheel.spin();
    STATE q=states[x];
    path.push_back(q);
    updateCoords(hmm->getStateType(q),i,j);
    prevQ=q;
  }
  hmm->decode(path,aln->fwdMap,aln->backMap,true);
  return aln;
}



/****************************************************************
                  Application::samplePair()
 ****************************************************************/
struct Edge {
  Edge() {}
  Edge(int i,int j,float d) : tax1(i), tax2(j), distance(d), p(1/d) {}
  int tax1, tax2; // INVARIANT: tax1<tax2
  float distance, p;
};
void Application::samplePair(int &tax1,int &tax2,int numNodes,
			     Array2D<bool> &equivClasses)
{
  Vector<Edge> freeEdges;
  for(int i=0 ; i<numNodes ; ++i) {
    if(!taxa[i].isLeaf()) continue;
    for(int j=i+1 ; j<numNodes ; ++j) {
      if(!taxa[j].isLeaf() || equivClasses[i][j]) continue;
      float dist=tree->distanceBetweenLeaves(i,j);
      freeEdges.push_back(Edge(i,j,dist));
    }
  }
  int n=freeEdges.size();
  float sum=0;
  for(int i=0 ; i<n ; ++i) {
    Edge &e=freeEdges[i];
    sum+=e.p;
  }
  RouletteWheel wheel;
  for(int i=0 ; i<n ; ++i) {
    Edge &e=freeEdges[i];
    wheel.addSector(e.p/sum);
  }
  wheel.doneAddingSectors();
  Edge &e=freeEdges[wheel.spin()];
  tax1=e.tax1;
  tax2=e.tax2;
}



/****************************************************************
                 Application::updateAncestors()
 ****************************************************************/
void Application::updateAncestors(PairwiseAlignment &aln)
{
  Taxon &taxon1=taxa[aln.tax1], &taxon2=taxa[aln.tax2];
  Vector<PhylogenyNode*> path;
  tree->getPath(taxon1.getNode(),taxon2.getNode(),path);
  IndexMap fwd=aln.fwdMap;
  PhylogenyNode *fromNode=path[0];
  path.cut(0);
  int n=path.size();
  for(int i=0 ; i<n ; ++i) {
    PhylogenyNode *toNode=path[0];
    path.cut(0);
    Taxon &from=taxa[fromNode->getID()], &to=taxa[toNode->getID()];
    if(to.getSeqLen()==0) extendSimple(from,to,fwd);
    else if(to.isLeaf()) extendIntoLeaf(from,to,fwd);
    else extend(from,to,fwd);
    fromNode=toNode;
  }
  for(int i=0 ; i<linkSeqs.size() ; ++i)
    linkSeqs[i].assignPositions();
  int numBranches=branchAttributes.size();
  for(int i=0 ; i<numBranches ; ++i) {
    BranchAttributes &branch=branchAttributes[i];
    resolveIndexMaps(branch);
  }
}



/****************************************************************
                  Application::allocateLinkSeq()
 ****************************************************************/
void Application::allocateLinkSeq(LinkSeq &S,int len)
{
  BOOM::Symbol star=alphabet.lookup('.');
  LinkedResidue r(star);
  for(int i=0 ; i<len ; ++i) S.append(r);
}



/****************************************************************
                  Application::getDir()
 ****************************************************************/
BranchDirection Application::getDir(Taxon &from,Taxon &to)
{
  BranchAttributes *branch=from.getBranchTo(to);
  if(&branch->getParentTaxon()==&to) return BRANCH_UP;
  return to.whichChild()==LEFT ? BRANCH_LEFT : BRANCH_RIGHT;
}



/****************************************************************
                  Application::extendSimple()
 ****************************************************************/
void Application::extendSimple(Taxon &from,Taxon &to,IndexMap &fwd)
{
  BranchDirection dir=getDir(from,to), revDir=getDir(to,from);
  int newL=fwd.countNonzero();
  LinkSeq &fromSeq=linkSeqs[from.getID()], &toSeq=linkSeqs[to.getID()];
  to.setSeqLen(newL);
  allocateLinkSeq(toSeq,newL);
  IndexMap newFwd(newL);
  LinkSeq::iterator fCur=fromSeq.begin(), fEnd=fromSeq.end();
  LinkSeq::iterator tCur=toSeq.begin(), tEnd=toSeq.end();
  int nextTarget=0, fwdL=fwd.size();
  for(int i=0 ; i<fwdL ; ++i, ++fCur) {
    if(fwd[i]!=IndexMap::UNDEFINED) {
      LinkedResidue &fromRes=*fCur, &toRes=*tCur;
      fromRes.link[dir]=&toRes;
      toRes.link[revDir]=&fromRes;
      newFwd[nextTarget]=fwd[i];
      ++nextTarget;
      ++tCur;
    }
  }
  fwd=newFwd;
}



/****************************************************************
                  Application::extend()
 ****************************************************************/
void Application::extend(Taxon &from,Taxon &to,IndexMap &fwd)
{
  BranchDirection dir=getDir(from,to), revDir=getDir(to,from);
  LinkSeq &fromSeq=linkSeqs[from.getID()], &toSeq=linkSeqs[to.getID()];
  LinkSeq::iterator fCur=fromSeq.begin(), fEnd=fromSeq.end();
  LinkSeq::iterator tCur=toSeq.begin(), tEnd=toSeq.end();
  int fromL=fwd.size(), targetPos=0;
  for(int i=0 ; i<fromL ; ++i, ++fCur) {
    if(fwd[i]==IndexMap::UNDEFINED) continue;
    LinkedResidue &fRes=*fCur;
    if(fRes.link[dir]) {
      while(&*tCur!=fRes.link[dir]) { ++tCur; ++targetPos; }
      //newFwd[targetPos]=fwd[i];
      continue;
    }
    LinkedResidue r(dot); r.pos=666;
    toSeq.insertAfter(r,tCur);
    ++tCur;
    ++targetPos;
    fRes.link[dir]=&*tCur;
    (*tCur).link[revDir]=&fRes;
    //cout<<"p "<<targetPos<<" "<<i<<" "<<toL<<" "<<fwd.size()<<endl;
    //newFwd[targetPos]=fwd[i];
  }
  int toL=toSeq.size();
  to.setSeqLen(toL);
  IndexMap newFwd(toL); newFwd.clear();
  fCur=fromSeq.begin();
  toSeq.assignPositions();
  //cout<<"targetPos="<<targetPos<<" fwd.size="<<fwd.size()<<" fwd.nonzero="<<fwd.countNonzero()<<endl;
  //cout<<"dir="<<dir<<" revdir="<<revDir<<" from="<<from.getName()<<" to="<<to.getName()<<endl;
  //cout<<"XXX "<<toSeq<<endl;
  for(int i=0; fCur!=fEnd ; ++fCur, ++i) {
    LinkedResidue &fRes=*fCur;
    if(fRes.link[dir]) {
      //cout<<"d "<<*fRes.link[dir]<<endl;
      //cout<<"e "<<i<<" "<<fwd[i]<<endl;
      newFwd[fRes.link[dir]->pos]=fwd[i];
    }
  }
  fwd=newFwd;
}



/****************************************************************
                  Application::extendIntoLeaf()
 ****************************************************************/
void Application::extendIntoLeaf(Taxon &from,Taxon &to,IndexMap &fwd)
{
  BranchDirection dir=getDir(from,to), revDir=getDir(to,from);
  LinkSeq &fromSeq=linkSeqs[from.getID()], &toSeq=linkSeqs[to.getID()];
  LinkSeq::iterator fCur=fromSeq.begin(), fEnd=fromSeq.end();
  LinkSeq::iterator tCur=toSeq.begin(), tEnd=toSeq.end();
  int fromL=fwd.size(), toL=to.getSeqLen(), targetPos=0;
  allocateLinkSeq(toSeq,toL);
  toSeq.assignPositions();
  for(int i=0 ; i<fromL ; ++i, ++fCur) {
    if(fwd[i]==IndexMap::UNDEFINED) continue;
    LinkedResidue &fRes=*fCur;
    while((*tCur).pos<fwd[i]) ++tCur;
    LinkedResidue &tRes=*tCur;
    fRes.link[dir]=&tRes;
    tRes.link[revDir]=&fRes;
  }
}



/****************************************************************
                  Application::resolveIndexMaps()
 ****************************************************************/
void Application::resolveIndexMaps(BranchAttributes &branch)
{
  IndexMap &upMap=branch.getUpMap(), &downMap=branch.getDownMap();
  Taxon &parent=branch.getParentTaxon(), &child=branch.getChildTaxon();
  int parentID=parent.getID(), childID=child.getID();
  BranchDirection dir=getDir(parent,child), revDir=getDir(child,parent);
  LinkSeq &parentSeq=linkSeqs[parentID], &childSeq=linkSeqs[childID];
  LinkSeq::iterator pCur=parentSeq.begin(), pEnd=parentSeq.end();
  LinkSeq::iterator cCur=childSeq.begin(), cEnd=childSeq.end();
  downMap.resize(parent.getSeqLen()); upMap.resize(child.getSeqLen());
  downMap.clear(); upMap.clear();
  for(; pCur!=pEnd ; ++pCur) {
    LinkedResidue r=*pCur;
    LinkedResidue *p=r.link[dir];
    if(p) downMap[r.pos]=p->pos;
  }
  for(; cCur!=cEnd ; ++cCur) {
    LinkedResidue r=*cCur;
    LinkedResidue *p=r.link[revDir];
    if(p) upMap[r.pos]=p->pos;
  }
}



/****************************************************************
                     struct LinkedResidue
 ****************************************************************/
LinkedResidue::LinkedResidue() 
{
  for(int i=0 ; i<3 ; ++i) link[i]=NULL;
}



LinkedResidue::LinkedResidue(BOOM::Symbol s)
  : s(s)
{
  for(int i=0 ; i<3 ; ++i) link[i]=NULL;
}



LinkedResidue::LinkedResidue(const LinkedResidue &other) 
{
  s=other.s;
  pos=other.pos;
  for(int i=0 ; i<3 ; ++i) link[i]=other.link[i];
};



void LinkedResidue::removeLinks()
{
  for(int i=0 ; i<3 ; ++i) link[i]=NULL;
}



void LinkedResidue::printOn(ostream &os) const
{
  os<<int(s)<<"@"<<pos<<"["<<link[0]<<","<<link[1]<<","<<link[2]<<"]";
}



ostream &operator<<(ostream &os,const LinkedResidue &r)
{
  r.printOn(os);
  return os;
}



ostream &operator<<(ostream &os,BranchDirection d)
{
  switch(d) 
    {
    case BRANCH_UP: os<<"UP"; break;
    case BRANCH_LEFT: os<<"LEFT"; break;
    case BRANCH_RIGHT: os<<"RIGHT"; break;
    }
  return os;
}




/****************************************************************
                         class LinkSeq
 ****************************************************************/
void LinkSeq::assignPositions()
{
  iterator cur=begin(), end=this->end();
  for(int i=0 ; cur!=end ; ++cur, ++i)
    (*cur).pos=i;
}



void LinkSeq::removeLinks()
{
  iterator cur=begin(), end=this->end();
  for(int i=0 ; cur!=end ; ++cur, ++i) 
    (*cur).removeLinks();
}



void LinkSeq::printOn(ostream &os) const
{
  const_iterator cur=begin(), end=this->end();
  for(int i=0 ; cur!=end ; ++cur, ++i)
    os<<(*cur)<<" ";
}



ostream &operator<<(ostream &os,const LinkSeq &s) 
{
  s.printOn(os);
  return os;
}




/****************************************************************
                     struct pairwiseAlignment
 ****************************************************************/
void PairwiseAlignment::invert()
{
  int tmp=tax1;
  tax1=tax2;
  tax2=tmp;
  IndexMap t=fwdMap;
  fwdMap=backMap;
  backMap=t;
}


