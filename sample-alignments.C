/****************************************************************
 sample-alignments.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "sample-alignments.H"
#include "BOOM/PowerMean.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "PosteriorBackward.H"

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
  CommandLine cmd(argc,argv,"b:Gr:c:pMti:s:");
  if(cmd.numArgs()!=5)
    throw String("\n\
sample-alignments <model.lambda> <library> <*.fasta> <#samples> <outfile>\n\
            -b <strategy> = band using strategy:\n\
                f### = fixed-width banding of width 2*###\n\
            -G = don't perform garbage collection\n\
            -r N = use randomization seed N\n\
            -c N = use pseudocount of N (default=0.01)\n\
            -p = print each alignment after sampling it\n\
            -M = use semi-Markov emissions instead of pseudocounts\n\
            -t = enforce transitivity in match probabilities\n\
            -i X = use indel coefficient X (default=1.0)\n\
            -s N = skip first N samples (default=0)\n\
");
  String lambdaFile=cmd.arg(0);
  String libraryDir=cmd.arg(1);
  String fastaFile=cmd.arg(2);
  numSamples=cmd.arg(3).asInt();
  outfile=cmd.arg(4);
  samplerType=(cmd.option('g') ? GIBBS : METROPOLIS_HASTINGS);
  wantSemiMarkov=cmd.option('M');
  wantTransitivity=cmd.option('t');
  if(wantTransitivity) cout<<"transitivity enabled"<<endl;
  //if(samplerType==GIBBS) throw "option -g is currently disabled";
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
  getHMMprobs();

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
  Vector<PhylogenyBranch*> branches;
  tree->collectBranches(branches,PREORDER);
  Vector<int> taxonIndices;
  for(int i=0 ; i<numTaxa ; ++i) taxonIndices.push_back(i);
  for(int i=0 ; i<numSamples ; ++i) {
    cout<<"sample #"<<i<<endl;
    VectorSorter<PhylogenyBranch*>::shuffle(branches);
    VectorSorter<int>::shuffle(taxonIndices);
    sample(f,branches,taxonIndices,rog,i);
    if(i>=skipSamples) {
      rog.saveConnections(f);
      if(wantPrint) {
	reconstructAlignment(true);
	fullAlignment->printSlice(cout,0,fullAlignment->getLength(),
				  '+',60,true);
      }
    }
  }
}



/****************************************************************
                     Application::sample()
 ****************************************************************/
void Application::sample(File &f,Vector<PhylogenyBranch*> &branches,
			 Vector<int> &taxonIndices,ResidueOrthologyGraph &rog,
			 int sampleNum)
{
  cout<<"sample branches (#"<<sampleNum<<")"<<endl;
  int numBranches=branches.size();
  for(int branchIndex=0 ; branchIndex<numBranches ; ++branchIndex) {
    CollapsedOrthologyMatrix::compute(rog);
    PhylogenyBranch *branch=branches[branchIndex];
    cout<<"sampling branch "<<branch->getParent()->getName()<<":"
	<<branch->getChild()->getName()<<endl;
    
    { //### DEBUGGING
      //if(!(branch->getChild()->getName()=="dm3")) continue;
    } //###

    sampleBranch2(branch);
    //sampleBranch(branch);
    /*
    if(wantPrint) {
      reconstructAlignment(true); //false);
      fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,
				true);
    }
    */
  }
}



/****************************************************************
                   Application::getLeafVector()
 ****************************************************************/
void Application::getLeafVector(BitSet &s,Vector<int> &v) 
{
  v.clear();
  for(int i=0 ; i<numTaxa ; ++i)
    if(s.isMember(i) && taxa[i].isLeaf())
      v.push_back(i);
}



/****************************************************************
                     Application::getLeaves()
 ****************************************************************/
void Application::getLeaves(Vector<int> &insideLeaves,
			    Vector<int> &outsideLeaves,
			    Taxon &parent,Taxon &child)
{
  BitSet insideTaxa=child.getCladeMembers();
  BitSet outsideTaxa(numTaxa);
  outsideTaxa.addAll();
  outsideTaxa-=insideTaxa;
  getLeafVector(insideTaxa,insideLeaves);
  getLeafVector(outsideTaxa,outsideLeaves);
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
                    Application::updateIndels()
 ****************************************************************/
void Application::updateIndels(BranchAttributes &branch,Taxon &parent,
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
    case PHMM_DELETE:
      downMap[i]=IndexMap::UNDEFINED; 
      ++i;
      break;
    case PHMM_INSERT:
      upMap[j]=IndexMap::UNDEFINED;
      ++j;
      break;
    default: throw "REAPR::updateIndels()";
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



CollapsedOrthologyMatrix::Entry 
  Application::indexCOM(CollapsedOrthologyMatrix &COM,int pos,int leaf)
{
  CollapsedOrthologyMatrix::Entry e=COM(pos,leaf);
  if(e.end>e.begin) --e.end; // noninclusive for indel states
  return e;
}



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



void Application::computeReachabilityVector(CollapsedOrthologyMatrix &com,
					    Vector<int> &leaves,
					    Vector<bool> &V)
{
  int L=com.getSeqLen();
  V.resize(L+1);
  for(int i=0 ; i<=L ; ++i)
    V[i]=reachesLeaf(com,leaves,i);
}



bool Application::reachesLeaf(CollapsedOrthologyMatrix &com,
			      Vector<int> &leaves,int pos)
{
  Vector<int>::iterator cur=leaves.begin(), end=leaves.end();
  for(; cur!=end ; ++cur) {
    CollapsedOrthologyMatrix::Entry e=com(pos,*cur);
    if(e.begin==e.end) return true;
  }
  return false;
}



bool Application::isDeadEnd(Taxon &taxon,int pos,InsideOutside io,
			    Taxon *exclude)
{
  int numBranches=taxon.getNumBranches();
  switch(io) 
    {
    case INSIDE:
      for(int i=0 ; i<numBranches ; ++i) {
	BranchAttributes *branch=taxon.getIthBranch(i);
	if(branch->getDownMap()[pos]!=IndexMap::UNDEFINED) return false;
      }
      break;
    case OUTSIDE: {
      BranchAttributes *parentBranch=taxon.getBranchToParent();
      if(parentBranch && 
	 parentBranch->getUpMap()[pos]!=IndexMap::UNDEFINED) return false;
      for(int i=0 ; i<numBranches ; ++i) {
	BranchAttributes *branch=taxon.getIthBranch(i);
	if(&branch->getChildTaxon()==exclude) continue;
	if(branch->getDownMap()[pos]!=IndexMap::UNDEFINED) return false;
      }
    } break; }
  return true;
}

    

/****************************************************************
                   Application::sandwichProb()
 ****************************************************************/
float Application::sandwichProb(InsideOutside io,PHMM_StateType stateType,
				Vector<int> &insideLeaves,
				Vector<int> &outsideLeaves,
				CollapsedOrthologyMatrix &childCOM,
				CollapsedOrthologyMatrix &parentCOM,
				int childPos,int parentPos,bool &justPseudo)
{
  int numInsideLeaves=insideLeaves.size();
  int numOutsideLeaves=outsideLeaves.size();
  Vector<float> sumV;
  for(int il=0 ; il<numInsideLeaves ; ++il) {
    int insideLeaf=insideLeaves[il];
    CollapsedOrthologyMatrix::Entry ie=indexCOM(childCOM,childPos,insideLeaf);
    if(io==INSIDE && childPos+1<childCOM.getSeqLen())
      ie+=indexCOM(childCOM,childPos+1,insideLeaf);
    int iBegin=ie.begin, iEnd=ie.end;
    if(iEnd>taxa[insideLeaf].getSeqLen()) --iEnd;
    if(iBegin<0) iBegin=0;
    for(int ol=0 ; ol<numOutsideLeaves ; ++ol) {
      int outsideLeaf=outsideLeaves[ol];
      SparseMatrix3D &M=*matrices[outsideLeaf][insideLeaf];
      float pathLen=tree->distanceBetweenLeaves(taxa[insideLeaf].getID(),
						taxa[outsideLeaf].getID());
      CollapsedOrthologyMatrix::Entry oe=indexCOM(parentCOM,parentPos,
						  outsideLeaf);
      if(io==OUTSIDE && parentPos+1<parentCOM.getSeqLen())
	oe+=indexCOM(parentCOM,parentPos+1,outsideLeaf);
      int oBegin=oe.begin, oEnd=oe.end;
      if(oEnd>taxa[outsideLeaf].getSeqLen()) --oEnd;
      if(oBegin<0) oBegin=0;
      for(int opos=oBegin ; opos<=oEnd ; ++opos) {
	EntryList &row=M(opos,stateType);
	EntryList::iterator cur=row.begin(), end=row.end();
	for(; cur!=end ; ++cur) {
	  SparseMatrix3D::Entry e=*cur;
	  if(e.y<iBegin) continue;
	  if(e.y>iEnd) break;
	  sumV.push_back(e.value-log(pathLen));
	}
      }
    }
  }
  if(sumV.size()==0) {
    sumV.push_back(logPseudoCount);
    justPseudo=true;
    //return 0; //###
  }
  else justPseudo=false;
  float ave=sumLogProbs(sumV)-log(sumV.size());
  return exp(ave);
}



void Application::updateTentacles(BranchAttributes &branch,Taxon &parent,
				  Taxon &child,
				  CollapsedOrthologyMatrix &parentCOM,
				  CollapsedOrthologyMatrix &childCOM,
				  int parentL,int childL,
				  Vector<int> &outsideLeaves,
				  Vector<int> &insideLeaves,
				  int numOutsideLeaves,int numInsideLeaves,
				  StatePath &path)
{

  double FUDGE=30.0;//###DEBUGGING
  double exponent=1.0;//5.0;

  /*
   First, sample events that change the parent & child sequence lengths
  */
  int i=0, j=0;
  IndexMap &upMap=branch.getUpMap(), &downMap=branch.getDownMap();
  PHMM_StateType q=0;
  Vector<EditEvent> parentEvents, childEvents;
  while(i<parentL || j<childL) {
    PHMM_StateType newQ=classifyState(upMap,downMap,i,j);
    switch(newQ) 
      {
      case PHMM_MATCH: {
	Vector<InsideOutside> DE;
	if(isDeadEnd(parent,i,OUTSIDE,&child)) DE.push_back(OUTSIDE);
	if(!child.isLeaf() && isDeadEnd(child,j,INSIDE)) DE.push_back(INSIDE);
	if(DE.size()==0) break;

	if(DE.size()==2) { // REQUIRES CHANGES TO accommodateEdits()
	  if(Random0to1()<0.5)
	    childEvents.push_back(EditEvent(DELETE,j));
	  else
	    parentEvents.push_back(EditEvent(DELETE,i));
	  break;
	}

	int which=RandomNumber(DE.size());
	switch(DE[which]) {
	case INSIDE: 
	  {
	    /*
	      float delP=transProbs[q][PHMM_DELETE];
	      float matchP=transProbs[q][PHMM_MATCH];
	      float denom=delP+matchP;
	      if(Random0to1()<delP/denom) 
	      childEvents.push_back(EditEvent(DELETE,j));
	    */
	    bool justPseudo1, justPseudo2;
	    float delP=sandwichProb(INSIDE,PHMM_DELETE,insideLeaves,
				    outsideLeaves,childCOM,parentCOM,j,i,
				    justPseudo1);
	    float matchP=sandwichProb(INSIDE,PHMM_MATCH,insideLeaves,
				      outsideLeaves,childCOM,parentCOM,j,i,
				      justPseudo2);
	    if(justPseudo1 || justPseudo2) {
	      delP=transProbs[q][PHMM_DELETE];
	      matchP=transProbs[q][PHMM_MATCH];
	    }
	    delP=powf(delP,exponent);
	    matchP=powf(matchP,exponent);
	    float denom=delP+matchP;
	    if(Random0to1()<delP/denom/FUDGE) 
	      childEvents.push_back(EditEvent(DELETE,j));
	  } break;
	case OUTSIDE: 
	  {
	    /*
	      float insP=transProbs[q][PHMM_INSERT];
	      float matchP=transProbs[q][PHMM_MATCH];
	      float denom=insP+matchP;
	      if(Random0to1()<insP/denom)   
	      parentEvents.push_back(EditEvent(DELETE,i));
	    */
	    bool justPseudo1, justPseudo2;
	    float insP=sandwichProb(OUTSIDE,PHMM_INSERT,insideLeaves,
				    outsideLeaves,childCOM,parentCOM,j,i,
				    justPseudo1);
	    float matchP=sandwichProb(OUTSIDE,PHMM_MATCH,insideLeaves,
				      outsideLeaves,childCOM,parentCOM,j,i,
				      justPseudo2);
	    if(justPseudo1 || justPseudo2) {
	      insP=transProbs[q][PHMM_INSERT];
	      matchP=transProbs[q][PHMM_MATCH];
	    }
	    insP=powf(insP,exponent);
	    matchP=powf(matchP,exponent);
	    float denom=insP+matchP;
	    if(Random0to1()<insP/denom/FUDGE) 
	      parentEvents.push_back(EditEvent(DELETE,i));
	  } break; }
      } break;
      case PHMM_INSERT: {
	/*
	if(!reachesLeaf(childCOM,insideLeaves,j)) {
	  childEvents.push_back(EditEvent(DELETE,j));
	  break;
	}
	*/
	if(!reachesLeaf(childCOM,insideLeaves,j)) break;//###
	bool justPseudo1, justPseudo2;
	float insP=sandwichProb(OUTSIDE,PHMM_INSERT,insideLeaves,outsideLeaves,
				childCOM,parentCOM,j,i,justPseudo1);
	float matchP=sandwichProb(OUTSIDE,PHMM_MATCH,insideLeaves,
				  outsideLeaves,childCOM,parentCOM,j,i,
				  justPseudo2);
	if(justPseudo1 || justPseudo2) {
	  insP=transProbs[q][PHMM_INSERT];
	  matchP=transProbs[q][PHMM_MATCH];
	}
	insP=powf(insP,exponent);
	matchP=powf(matchP,exponent);
	float denom=insP+matchP;
	if(Random0to1()>insP/denom/FUDGE) 
	  parentEvents.push_back(EditEvent(INSERT,i,j));//grow the tentacle
      } break;
      case PHMM_DELETE:
	if(!reachesLeaf(parentCOM,outsideLeaves,i)) {
	  if(isDeadEnd(parent,i,OUTSIDE,&child)) 
	    parentEvents.push_back(EditEvent(DELETE,i));//### ???
	  break;//###
	}
	if(child.isLeaf()) break;
	bool justPseudo1, justPseudo2;
	float delP=sandwichProb(INSIDE,PHMM_DELETE,insideLeaves,outsideLeaves,
				childCOM,parentCOM,j,i,justPseudo1);
	float matchP=sandwichProb(INSIDE,PHMM_MATCH,insideLeaves,
				  outsideLeaves,childCOM,parentCOM,j,i,
				  justPseudo2);
	delP=powf(delP,exponent);
	matchP=powf(matchP,exponent);
	if(justPseudo1 || justPseudo2) {
	  delP=transProbs[q][PHMM_DELETE];
	  matchP=transProbs[q][PHMM_MATCH];
	}
	float denom=delP+matchP;
	if(Random0to1()>delP/denom/FUDGE) 
	  childEvents.push_back(EditEvent(INSERT,j,i));//grow the tentacle
	break;
      }
    updateCoords(newQ,i,j);
    q=newQ;
  }

  /* Now update all of the IndexMaps on the incident edges, based on the
     sampled changes to the sequence lengths */
  IndexMap childSeqMap;
  accommodateChildEdits(child,branch,childEvents,upMap,downMap,childSeqMap);
  accommodateParentEdits(parent,child,branch,parentEvents,upMap,downMap,
  			 childSeqMap);
}



void Application::getHMMprobs()
{
  transProbs.setAllTo(0.0);
  for(STATE i=0 ; i<numStates ; ++i) {
    PHMM_StateType iType=hmm->getStateType(i);
    if(iType>3) continue;
    for(STATE j=0 ; j<numStates ; ++j) {
      PHMM_StateType jType=hmm->getStateType(j);
      if(jType>3) continue;
      float P=exp(hmm->getTransP(i,j));
      transProbs[iType][jType]+=P;
    }
  }
}



void Application::getHMMprobs(PairHMM *hmm)
{
  transProbs.setAllTo(0.0);
  for(STATE i=0 ; i<numStates ; ++i) {
    PHMM_StateType iType=hmm->getStateType(i);
    if(iType>3) continue;
    for(STATE j=0 ; j<numStates ; ++j) {
      PHMM_StateType jType=hmm->getStateType(j);
      if(jType>3) continue;
      float P=exp(hmm->getTransP(i,j));
      transProbs[iType][jType]+=P;
    }
  }
}



int Application::getNewLength(int oldLen,Vector<EditEvent> &edits)
{
  int newLen=oldLen;
  Vector<EditEvent>::iterator cur=edits.begin(), end=edits.end();
  for(; cur!=end ; ++cur)
    if((*cur).type==INSERT) ++newLen; else --newLen;
  return newLen;
}



void Application::computeSeqMap(int oldLen,Vector<EditEvent> &events,
				IndexMap &seqMap)
{
  seqMap.resize(oldLen);
  Vector<EditEvent>::iterator cur=events.begin(), end=events.end();
  for(int pos=0, image=0 ; pos<oldLen ;) {
    if(cur!=end) {
      const EditEvent &e=*cur;
      while(pos<e.pos) seqMap[pos++]=image++;
      if(e.type==INSERT) ++image;
      else seqMap[pos++]=IndexMap::UNDEFINED;
      ++cur;
    }
    else while(pos<oldLen) seqMap[pos++]=image++;
  }
}



void Application::accommodateChildEdits(Taxon &child,
					BranchAttributes &branch,
					Vector<EditEvent> &childEvents,
					IndexMap &upMap,IndexMap &downMap,
					IndexMap &childSeqMap)
{
  // Update the sequence map, which maps old positions to new positions
  int oldLen=child.getSeqLen();
  int newLen=getNewLength(oldLen,childEvents);
  computeSeqMap(oldLen,childEvents,childSeqMap);
  if(!childSeqMap.sanityCheck(newLen)) INTERNAL_ERROR;
  if(child.isLeaf()) return;
  IndexMap invSeqMap;
  childSeqMap.invert(newLen,invSeqMap);
  if(!invSeqMap.sanityCheck(oldLen)) INTERNAL_ERROR;

  // Update the parent's downmap
  downMap.compose(childSeqMap,downMap);
  if(!downMap.sanityCheck(newLen)) INTERNAL_ERROR;

  // Update the grandchildrens' upmaps
  BranchAttributes *childBranch1=child.getIthBranch(0); 
  BranchAttributes *childBranch2=child.getIthBranch(1);
  IndexMap &um1=childBranch1->getUpMap();
  IndexMap &um2=childBranch2->getUpMap();
  um1.compose(childSeqMap,um1);
  um2.compose(childSeqMap,um2);
  if(!um1.sanityCheck(newLen)) INTERNAL_ERROR;
  if(!um2.sanityCheck(newLen)) INTERNAL_ERROR;

  // Update the maps pointing away from this child
  IndexMap newUpMap(newLen);
  newUpMap.clear();
  invSeqMap.compose(upMap,newUpMap);
  if(!newUpMap.sanityCheck(branch.getParentTaxon().getSeqLen())) INTERNAL_ERROR;
  IndexMap temp=invSeqMap;
  Vector<EditEvent>::iterator cur=childEvents.begin(), end=childEvents.end();
  for(; cur!=end ; ++cur) {
    const EditEvent &e=*cur;
    if(e.type!=INSERT) continue;
    int pos=findStartingPos(e.pos,childSeqMap);
    while(temp[pos]!=IndexMap::UNDEFINED) ++pos;
    temp[pos]=666;
    newUpMap[pos]=e.mapTo;
    downMap[e.mapTo]=pos;
  }
  if(!newUpMap.sanityCheck(branch.getParentTaxon().getSeqLen())) INTERNAL_ERROR;
  if(!downMap.sanityCheck(newLen)) INTERNAL_ERROR;
  upMap=newUpMap;
  IndexMap &dm1=childBranch1->getDownMap(), newDM1(newLen);
  IndexMap &dm2=childBranch2->getDownMap(), newDM2(newLen);
  newDM1.setAllTo(IndexMap::UNDEFINED);
  newDM2.setAllTo(IndexMap::UNDEFINED);
  invSeqMap.compose(dm1,newDM1);
  invSeqMap.compose(dm2,newDM2);
  dm1=newDM1;
  dm2=newDM2;
  if(!dm1.sanityCheck(childBranch1->getChildTaxon().getSeqLen())) INTERNAL_ERROR;
  if(!dm2.sanityCheck(childBranch2->getChildTaxon().getSeqLen())) INTERNAL_ERROR;

  // Update the taxon's seqlen
  child.setSeqLen(newLen);
}



void Application::accommodateParentEdits(Taxon &parent,Taxon &child,
					 BranchAttributes &branch,
					 Vector<EditEvent> &parentEvents, 
					 IndexMap &upMap,IndexMap &downMap,
					 const IndexMap &childSeqMap)
{
  // Compute a sequence map for the parent
  int oldLen=parent.getSeqLen();
  int newLen=getNewLength(oldLen,parentEvents);
  IndexMap parentSeqMap(oldLen);
  computeSeqMap(oldLen,parentEvents,parentSeqMap);
  if(!parentSeqMap.sanityCheck(newLen)) INTERNAL_ERROR;
  IndexMap invSeqMap;
  parentSeqMap.invert(newLen,invSeqMap);
  if(!invSeqMap.sanityCheck(oldLen)) INTERNAL_ERROR;

  // Update all IndexMaps that point to this taxon
  BranchAttributes *parentBranchUp=parent.getBranchToParent(); 
  BranchAttributes *parentBranchDown=
    parent.getIthBranch(otherChild(child.whichChild()));
  upMap.compose(parentSeqMap,upMap);
  if(!upMap.sanityCheck(newLen)) INTERNAL_ERROR;

  if(parentBranchUp) {
    IndexMap &dm=parentBranchUp->getDownMap();
    dm.compose(parentSeqMap,dm);
    if(!dm.sanityCheck(newLen)) INTERNAL_ERROR;
  }
  IndexMap &um=parentBranchDown->getUpMap();
  if(!um.sanityCheck(oldLen)) INTERNAL_ERROR;
  um.compose(parentSeqMap,um);
  if(!um.sanityCheck(newLen)) INTERNAL_ERROR;

  // Update the maps pointing away from this taxon
  IndexMap newDownMap;
  invSeqMap.compose(downMap,newDownMap);
  if(!newDownMap.sanityCheck(child.getSeqLen())) INTERNAL_ERROR;
  IndexMap temp=invSeqMap;
  Vector<EditEvent>::iterator cur=parentEvents.begin(), end=parentEvents.end();
  for(; cur!=end ; ++cur) {
    const EditEvent &e=*cur;
    if(e.type!=INSERT) continue;
    int pos=findStartingPos(e.pos,parentSeqMap);
    while(temp[pos]!=IndexMap::UNDEFINED) ++pos;
    temp[pos]=666;
    if(pos>=invSeqMap.size() || pos>=newDownMap.size()) INTERNAL_ERROR;//###
    int newTo=childSeqMap[e.mapTo];
    newDownMap[pos]=newTo;
    upMap[newTo]=pos;
  }

  // ### THIS IS IT
  if(!newDownMap.sanityCheck(child.getSeqLen())) //INTERNAL_ERROR;//<<<<<<<<XXX
    {
      //cout<<parentSeqMap<<endl;
      IndexMap newDownMap;
      invSeqMap.compose(downMap,newDownMap);
      //newDownMap.compose(childSeqMap,newDownMap);
      cout<<newDownMap<<endl;
      IndexMap temp=invSeqMap;
      Vector<EditEvent>::iterator cur=parentEvents.begin(), 
	end=parentEvents.end();
      for(; cur!=end ; ++cur) {
	const EditEvent &e=*cur;
	if(e.type!=INSERT) continue;
	cout<<e<<endl;
	int pos=findStartingPos(e.pos,parentSeqMap);
	cout<<"initial pos="<<pos<<endl;
	//while(invSeqMap[pos]!=IndexMap::UNDEFINED) ++pos;
	//while(newDownMap[pos]!=IndexMap::UNDEFINED) ++pos;
	while(temp[pos]!=IndexMap::UNDEFINED) ++pos;
	temp[pos]=666;
	if(pos>=invSeqMap.size() || pos>=newDownMap.size()) INTERNAL_ERROR;//###
	int newTo=childSeqMap[e.mapTo];
	cout<<"setting "<<pos<<"->"<<newTo<<endl;
	newDownMap[pos]=newTo;
	upMap[newTo]=pos;
	if(!newDownMap.sanityCheck(child.getSeqLen())) {
	  cout<<"childSeqMap: "<<childSeqMap<<endl;
	  cout<<"temp: "<<temp<<endl;
	  cout<<"invSeqMap: "<<invSeqMap<<endl;
	  cout<<"downMap: "<<downMap<<endl;
	  cout<<"newDownMap: "<<newDownMap<<endl;
	  //cout<<temp<<endl;
	  INTERNAL_ERROR;
	}
      }
    }
  // ###


  if(!upMap.sanityCheck(newLen)) INTERNAL_ERROR;
  downMap=newDownMap;
  IndexMap &dm=parentBranchDown->getDownMap();
  if(!dm.sanityCheck(parentBranchDown->getChildTaxon().getSeqLen())) INTERNAL_ERROR;
  newDownMap.setAllTo(IndexMap::UNDEFINED);
  invSeqMap.compose(dm,newDownMap);
  dm=newDownMap;
  if(!dm.sanityCheck(parentBranchDown->getChildTaxon().getSeqLen())) INTERNAL_ERROR;
  if(parentBranchUp) {
    IndexMap &um=parentBranchUp->getUpMap();
    IndexMap newUpMap(newLen);
    newUpMap.setAllTo(IndexMap::UNDEFINED);
    invSeqMap.compose(um,newUpMap);
    um=newUpMap;
    if(!um.sanityCheck(parentBranchUp->getParentTaxon().getSeqLen())) INTERNAL_ERROR;
  }

  // Update the taxon's seqlen
  parent.setSeqLen(newLen);
  //cout<<"setting parent len="<<newLen<<endl;
}



int Application::findStartingPos(int epos,IndexMap &seqMap)
{
  int i;
  for(i=epos-1 ; i>0 && seqMap[i]==IndexMap::UNDEFINED ; --i);
  return i>=0 && seqMap[i]!=IndexMap::UNDEFINED ? seqMap[i]+1 : 0;
}



/****************************************************************
                     Application::sampleBranch()
 ****************************************************************/
void Application::sampleBranch(PhylogenyBranch *br)
{
  BranchAttributes *branch=
    static_cast<BranchAttributes*>(br->getDecoration());
  getHMMprobs(branch->getHMM());
  Taxon &parent=branch->getParentTaxon();
  Taxon &child=branch->getChildTaxon();
  CollapsedOrthologyMatrix &parentCOM=*parent.getCOM();
  CollapsedOrthologyMatrix &childCOM=*child.getCOM();
  int parentL=parent.getSeqLen(), childL=child.getSeqLen();
  Vector<int> outsideLeaves, insideLeaves;
  getLeaves(insideLeaves,outsideLeaves,parent,child);
  int numInsideLeaves=insideLeaves.size();
  int numOutsideLeaves=outsideLeaves.size();
  int i=0, j=0; // 1-based
  StatePath &path=*new StatePath(hmm);
  Vector<bool> parentReachability, childReachability;
  computeReachabilityVector(parentCOM,outsideLeaves,parentReachability);
  computeReachabilityVector(childCOM,insideLeaves,childReachability);
  while(i<parentL || j<childL) {
    PHMM_StateType skipState=PHMM_SILENT;
    if(i<parentL && j<childL) {
      bool parentDead=!parentReachability[i+1];
      bool childDead=!childReachability[j+1];
      if(parentDead && childDead) {
	path.push_back(matchState);
	++i; ++j;
	continue;
      }
      else if(parentDead) skipState=PHMM_INSERT;
      else if(childDead) skipState=PHMM_DELETE;
    }
    Vector<STATE> states;
    Vector<float> probs;
    bool allPseudo=true;
    for(int q=1 ; q<numStates ; ++q) {
      PHMM_StateType stateType=stateTypes[q];
      int newI=i+dI[q], newJ=j+dJ[q];
      if(newI>parentL || newJ>childL) continue;
      bool justPseudo;
      //cout<<"computeCellP()"<<endl;
      float cellP=computeCellP(insideLeaves,outsideLeaves,numInsideLeaves,
			       numOutsideLeaves,childCOM,parentCOM,q,newI,
			       newJ,justPseudo,parentReachability,
			       childReachability);
      //cout<<"/computeCellP"<<endl;
      if(justPseudo) continue;
      else allPseudo=false;
      if(isFinite(cellP)) {
	probs.push_back(exp(cellP));
	states.push_back(q);
      }
    }
    if(allPseudo) {
      bool skipped=false;
      if(wantSemiMarkov) {
	/*
	skipped=semiMarkov(path,i,j,insideLeaves,outsideLeaves,numInsideLeaves,
			   numOutsideLeaves,childCOM,parentCOM,parentL,childL,
			   parentReachability,childReachability);*/
	//cout<<"calling manhattan"<<endl;
	int oldI=i, oldJ=j;
	//cout<<"parentL="<<parentL<<" childL="<<childL<<endl;
	skipped=manhattan(path,i,j,insideLeaves,outsideLeaves,numInsideLeaves,
			   numOutsideLeaves,childCOM,parentCOM,parentL,childL,
			   parentReachability,childReachability);
	if(i==oldI && j==oldJ) INTERNAL_ERROR;
	//cout<<"manhattan returned "<<skipped<<" "<<i<<" "<<j<<endl;
	if(skipped) continue;
      }
      lookahead(probs,states,insideLeaves,outsideLeaves,numInsideLeaves,
		numOutsideLeaves,childCOM,parentCOM,parentL,childL,i,j,
		skipState);
    }
    float sum=0;
    Vector<float>::iterator cur=probs.begin(), end=probs.end();
    for(; cur!=end ; ++cur) sum+=*cur;
    RouletteWheel wheel;
    cur=probs.begin();
    for(int i=0 ; cur!=end ; ++cur, ++i) wheel.addSector(probs[i]/sum);
    wheel.doneAddingSectors();
    STATE q=states[wheel.spin()];
    if(q==deleteState && skipState==PHMM_DELETE) q=matchState;
    else if(q==insertState && skipState==PHMM_INSERT) q=matchState;
    else if(q==matchState) {
      while(i<parentL-2 && !parentReachability[i+1]) {// semi-markov move
	path.push_back(deleteState);
	++i;
      }
      while(j<childL-2 && !childReachability[j+1]) { // semi-markov move
	path.push_back(insertState);
	++j;
      }
    }
    else if(q==deleteState && skipState==PHMM_INSERT) {
      while(i<parentL-1 && !parentReachability[i+1]) {// semi-markov move
	path.push_back(q);
	++i;
      }
    }
    else if(q==insertState && skipState==PHMM_DELETE) {
      while(j<childL-1 && !childReachability[j+1]) { // semi-markov move
	path.push_back(q);
	++j;
      }
    }
    path.push_back(q);
    updateCoords(hmm->getStateType(q),i,j);
  }
  //cout<<"updateIndels()"<<endl;
  updateIndels(*branch,parent,child,path);
  delete &path;//branch->setStatePath(&path);
  //cout<<"updateTentacles()"<<endl;
  updateTentacles(*branch,parent,child,parentCOM,childCOM,parentL,childL,
		  outsideLeaves,insideLeaves,numOutsideLeaves,numInsideLeaves,
		  path);
  //cout<<"/updateTentacles"<<endl;
  cout.flush();
}



/****************************************************************
                    Application::computeCellP()
 ****************************************************************/
float Application::computeCellP(Vector<int> &insideLeaves,
				Vector<int> &outsideLeaves,
				int numInsideLeaves,int numOutsideLeaves,
				CollapsedOrthologyMatrix &childCOM,
				CollapsedOrthologyMatrix &parentCOM,STATE q,
				int newI,int newJ,bool &justPseudo,
				Vector<bool> &parentReachability,
				Vector<bool> &childReachability)
{
  PHMM_StateType stateType=hmm->getStateType(q);

  int parentL=parentCOM.getSeqLen(), childL=childCOM.getSeqLen();
  switch(stateType)
    {
    case PHMM_MATCH:
      while(newJ<childL-1 && !childReachability[newJ]) ++newJ;
      while(newI<parentL-1 && !parentReachability[newI]) ++newI;
      break;
    case PHMM_INSERT:
      while(newJ<childL-1 && !childReachability[newJ]) ++newJ;
      break;
    case PHMM_DELETE:
      while(newI<parentL-1 && !parentReachability[newI]) ++newI;
      break;
    }
  if(!childReachability[newJ] && !parentReachability[newI]) {
    justPseudo=true;
    return NEGATIVE_INFINITY;
  }

  Vector<float> sumV;
  for(int il=0 ; il<numInsideLeaves ; ++il) {
    int insideLeaf=insideLeaves[il];
    CollapsedOrthologyMatrix::Entry ie=childCOM(newJ,insideLeaf);
    bool iDirect=ie.isDirectMatch();
    int iBegin=ie.begin, iEnd=ie.end;
    if(iEnd>ie.begin) --iEnd; // because for indel states it's noninclusive
    if(iEnd>taxa[insideLeaf].getSeqLen()) --iEnd;
    if(iBegin<0) iBegin=0;
    for(int ol=0 ; ol<numOutsideLeaves ; ++ol) {
      int outsideLeaf=outsideLeaves[ol];
      SparseMatrix3D &M=*matrices[outsideLeaf][insideLeaf];
      float pathLen=tree->distanceBetweenLeaves(taxa[insideLeaf].getID(),
						taxa[outsideLeaf].getID());
      CollapsedOrthologyMatrix::Entry oe=parentCOM(newI,outsideLeaf);
      bool oDirect=oe.isDirectMatch();
      Set<PHMM_StateType> stateTypes;
      bool DEL=(stateType==PHMM_DELETE || ie.end>ie.begin || 
		ie.begin<0 || ie.end>taxa[insideLeaf].getSeqLen());
      bool INS=(stateType==PHMM_INSERT || oe.end>oe.begin || 
		oe.begin<0 || oe.end>taxa[outsideLeaf].getSeqLen());
      if(!INS && !DEL) stateTypes.insert(PHMM_MATCH);
      else {
	if(INS) stateTypes.insert(PHMM_INSERT);
	if(DEL) stateTypes.insert(PHMM_DELETE);
      }
      int oBegin=oe.begin, oEnd=oe.end;
      if(oEnd>oe.begin) --oEnd; // because for indel states it's noninclusive
      if(oEnd>taxa[outsideLeaf].getSeqLen()) --oEnd;
      if(oBegin<0) oBegin=0;
      Vector<float> sumW;
      for(int opos=oBegin ; opos<=oEnd ; ++opos) {
	Set<PHMM_StateType>::iterator cur=stateTypes.begin(), 
	  end=stateTypes.end();
	for(; cur!=end ; ++cur) {
	  PHMM_StateType stateType=*cur;
	  EntryList &row=M(opos,stateType);
	  EntryList::iterator cur=row.begin(), end=row.end();
	  for(; cur!=end ; ++cur) {
	    SparseMatrix3D::Entry e=*cur;
	    if(e.y<iBegin) continue;
	    if(e.y>iEnd) break;
	    sumW.push_back(e.value);
	  }
	}
      }
      float theSum=sumLogProbs(sumW);
      if(isFinite(theSum)) sumV.push_back(theSum-log(pathLen));
      else if(wantTransitivity)
	if(stateType==PHMM_MATCH && iDirect && oDirect)
	  {justPseudo=true; return NEGATIVE_INFINITY;}

    }
  }
  if(sumV.size()==0) {
    justPseudo=true;
    return NEGATIVE_INFINITY;
  }
  else justPseudo=false;
  float ave=sumLogProbs(sumV)-log(sumV.size());
  //float ave=PowerMean::compute_log(sumV,5);//1.5);
  return ave;
}



/****************************************************************
                     Application::lookahead()
 ****************************************************************/
void Application::lookahead(Vector<float> &probs,Vector<STATE> &states,
			    Vector<int> &insideLeaves,
			    Vector<int> &outsideLeaves,int numInsideLeaves,
			    int numOutsideLeaves,
			    CollapsedOrthologyMatrix &childCOM,
			    CollapsedOrthologyMatrix &parentCOM,
			    int parentL,int childL,int i,int j,
			    PHMM_StateType skip)
{
  // THIS FUNCTION IS OBSOLETE

  // ###########
  INTERNAL_ERROR;
  // ###########


  Vector<float> oldProbs=probs;
  Vector<STATE> oldStates=states;
  probs.clear();
  states.clear();
  for(int q=1 ; q<numStates ; ++q) {
    PHMM_StateType stateType=stateTypes[q];
    if(stateType==skip) continue;
    int newI=i, newJ=j;
    Vector<float> sumV;
    while(1) {
      newI+=dI[q];
      newJ+=dJ[q];
      if(stateType==PHMM_MATCH) break;
      if(newI>parentL || newJ>childL) break;
      PHMM_StateType stateType=hmm->getStateType(q);
      for(int il=0 ; il<numInsideLeaves ; ++il) {
	int insideLeaf=insideLeaves[il];
	CollapsedOrthologyMatrix::Entry ie=childCOM(newJ,insideLeaf);
	for(int ol=0 ; ol<numOutsideLeaves ; ++ol) {
	  int outsideLeaf=outsideLeaves[ol];
	  SparseMatrix3D &M=*matrices[insideLeaf][outsideLeaf];
	  float pathLen=tree->distanceBetweenLeaves(taxa[insideLeaf].getID(),
						    taxa[outsideLeaf].getID());
	  CollapsedOrthologyMatrix::Entry oe=parentCOM(newI,outsideLeaf);
	  Set<PHMM_StateType> stateTypes;
	  bool INS=(stateType==PHMM_INSERT || ie.end>ie.begin || 
		    ie.begin<0 || ie.end>taxa[insideLeaf].getSeqLen());
	  bool DEL=(stateType==PHMM_DELETE || oe.end>oe.begin || 
		    oe.begin<0 || oe.end>taxa[outsideLeaf].getSeqLen());
	  if(!INS && !DEL) stateTypes.insert(PHMM_MATCH);
	  else {
	    if(INS) stateTypes.insert(PHMM_INSERT);
	    if(DEL) stateTypes.insert(PHMM_DELETE);
	  }
	  int iBegin=ie.begin, iEnd=ie.end, oBegin=oe.begin, oEnd=oe.end;
	  if(iEnd>ie.begin) --iEnd; // indel states: it's noninclusive
	  if(iEnd>taxa[insideLeaf].getSeqLen()) --iEnd;
	  if(iBegin<0) iBegin=0;
	  if(oEnd>oe.begin) --oEnd; // indel states: it's noninclusive
	  if(oEnd>taxa[outsideLeaf].getSeqLen()) --oEnd;
	  if(oBegin<0) oBegin=0;
	  for(int ipos=iBegin ; ipos<=iEnd ; ++ipos) {
	    Set<PHMM_StateType>::iterator cur=stateTypes.begin(), 
	      end=stateTypes.end();
	    for(; cur!=end ; ++cur) {
	      PHMM_StateType stateType=*cur;
	      EntryList &row=M(ipos,stateType);
	      EntryList::iterator cur=row.begin(), end=row.end();
	      for(; cur!=end ; ++cur) {
		SparseMatrix3D::Entry e=*cur;
		if(e.y<oBegin) continue;
		if(e.y>oEnd) break;
		sumV.push_back(e.value-log(pathLen));
	      }
	    }
	  }
	}
      }
      if(sumV.size()>0) break;
    }
    if(newI>parentL || newJ>childL) continue;
    if(sumV.size()==0) continue;
    float ave=sumLogProbs(sumV)-log(sumV.size());
    float cellP=ave;
    if(isFinite(cellP)) {
      probs.push_back(exp(cellP));
      states.push_back(q);
    }
  }
  if(states.size()==0) {
    states=oldStates;
    probs=oldProbs;
  }
}



/****************************************************************
                     Application::semiMarkov()
 ****************************************************************/
bool Application::semiMarkov(StatePath &path,int &i,int &j,
			     Vector<int> &insideLeaves,
			     Vector<int> &outsideLeaves,int numInsideLeaves,
			     int numOutsideLeaves,
			     CollapsedOrthologyMatrix &childCOM,
			     CollapsedOrthologyMatrix &parentCOM,
			     int parentL,int childL,
			     Vector<bool> &parentReachability,
			     Vector<bool> &childReachability)
{
  bool justPseudo;
  int newI=i+1, newJ=j+1;
  float insP=NEGATIVE_INFINITY, delP=NEGATIVE_INFINITY;
  for(; newI<=parentL || newJ<=childL ; ++newI, ++newJ) {
    if(newI<parentL && j<childL) {
      insP=computeCellP(insideLeaves,outsideLeaves,numInsideLeaves, 
			numOutsideLeaves,childCOM,parentCOM,matchState,
			i+1,newJ+1,justPseudo,parentReachability,
			childReachability);
      if(justPseudo) insP=NEGATIVE_INFINITY;
      if(isFinite(insP)) break;
    }
    if(newJ<childL && i<parentL) {
      delP=computeCellP(insideLeaves,outsideLeaves,numInsideLeaves, 
			numOutsideLeaves,childCOM,parentCOM,matchState,
			newI+1,j+1,justPseudo,parentReachability,
			childReachability);
      if(justPseudo) delP=NEGATIVE_INFINITY;
      if(isFinite(delP)) break;
    }
  }
  int offset=newI-i;
  if(isFinite(insP) && isFinite(delP)) {
    float ratio=exp(insP-sumLogProbs(insP,delP));
    if(Random0to1()>ratio) insP=NEGATIVE_INFINITY;
  }
  if(isFinite(insP)) {
    for(int k=0 ; k<offset ; ++k) path.push_back(insertState);
    j=newJ;
    return true;
  }
  else if(isFinite(delP)) {
    for(int k=0 ; k<offset ; ++k) path.push_back(deleteState);
    i=newI;
    return true;
  }
  if(Random0to1()<0.5) {
    for(; i<parentL ; ++i) path.push_back(deleteState);
    for(; j<childL ; ++j) path.push_back(insertState);
    return true;
  }
  else {
    for(; j<childL ; ++j) path.push_back(insertState);
    for(; i<parentL ; ++i) path.push_back(deleteState);
    return true;
  }

  return false;
}




bool Application::cellIsZero(Vector<int> &insideLeaves,
			      Vector<int> &outsideLeaves,
			      int numInsideLeaves,int numOutsideLeaves,
			      CollapsedOrthologyMatrix &childCOM,
			      CollapsedOrthologyMatrix &parentCOM,STATE q,
			      int newI,int newJ,
			      Vector<bool> &parentReachability,
			      Vector<bool> &childReachability)
{
  PHMM_StateType stateType=hmm->getStateType(q);
  int parentL=parentCOM.getSeqLen(), childL=childCOM.getSeqLen();
  Vector<float> sumV;
  for(int il=0 ; il<numInsideLeaves ; ++il) {
    int insideLeaf=insideLeaves[il];
    CollapsedOrthologyMatrix::Entry ie=childCOM(newJ,insideLeaf);
    int iBegin=ie.begin, iEnd=ie.end;
    if(iEnd>ie.begin) --iEnd; // because for indel states it's noninclusive
    if(iEnd>taxa[insideLeaf].getSeqLen()) --iEnd;
    if(iBegin<0) iBegin=0;
    for(int ol=0 ; ol<numOutsideLeaves ; ++ol) {
      int outsideLeaf=outsideLeaves[ol];
      SparseMatrix3D &M=*matrices[insideLeaf][outsideLeaf];
      CollapsedOrthologyMatrix::Entry oe=parentCOM(newI,outsideLeaf);
      Set<PHMM_StateType> stateTypes;
      bool INS=(stateType==PHMM_INSERT || ie.end>ie.begin || 
		ie.begin<0 || ie.end>taxa[insideLeaf].getSeqLen());
      bool DEL=(stateType==PHMM_DELETE || oe.end>oe.begin || 
		oe.begin<0 || oe.end>taxa[outsideLeaf].getSeqLen());
      if(!INS && !DEL) stateTypes.insert(PHMM_MATCH);
      else {
	if(INS) stateTypes.insert(PHMM_INSERT);
	if(DEL) stateTypes.insert(PHMM_DELETE);
      }
      int oBegin=oe.begin, oEnd=oe.end;
      if(oEnd>oe.begin) --oEnd; // because for indel states it's noninclusive
      if(oEnd>taxa[outsideLeaf].getSeqLen()) --oEnd;
      if(oBegin<0) oBegin=0;
      Vector<float> sumW;
      for(int ipos=iBegin ; ipos<=iEnd ; ++ipos) {
	Set<PHMM_StateType>::iterator cur=stateTypes.begin(), 
	  end=stateTypes.end();
	for(; cur!=end ; ++cur) {
	  PHMM_StateType stateType=*cur;
	  EntryList &row=M(ipos,stateType);
	  EntryList::iterator cur=row.begin(), end=row.end();
	  for(; cur!=end ; ++cur) {
	    SparseMatrix3D::Entry e=*cur;
	    if(e.y<oBegin) continue;
	    if(e.y>oEnd) break;
	    sumW.push_back(e.value);
	  }
	}
      }
      float theSum=sumLogProbs(sumW);
      if(isFinite(theSum)) sumV.push_back(theSum);
    }
  }
  if(sumV.size()==0) return true;
  float ave=sumLogProbs(sumV)-log(sumV.size());
  return !isFinite(ave);
}



/****************************************************************
                    Application::debugCellP()
 ****************************************************************/
float Application::debugCellP(Vector<int> &insideLeaves,
			      Vector<int> &outsideLeaves,
			      int numInsideLeaves,int numOutsideLeaves,
			      CollapsedOrthologyMatrix &childCOM,
			      CollapsedOrthologyMatrix &parentCOM,STATE q,
			      int newI,int newJ,bool &justPseudo,
			      Vector<bool> &parentReachability,
			      Vector<bool> &childReachability)
{
  PHMM_StateType stateType=hmm->getStateType(q);

  int parentL=parentCOM.getSeqLen(), childL=childCOM.getSeqLen();
  switch(stateType)
    {
    case PHMM_MATCH:
      while(newJ<childL-1 && !childReachability[newJ]) ++newJ;
      while(newI<parentL-1 && !parentReachability[newI]) ++newI;
      break;
    case PHMM_DELETE:
      while(newJ<childL-1 && !childReachability[newJ]) ++newJ;
      break;
    case PHMM_INSERT:
      while(newI<parentL-1 && !parentReachability[newI]) ++newI;
      break;
    }
  if(!childReachability[newJ] && !parentReachability[newI]) {
    justPseudo=true;
    return NEGATIVE_INFINITY;
  }

  Vector<float> sumV;
  for(int il=0 ; il<numInsideLeaves ; ++il) {
    int insideLeaf=insideLeaves[il];
    CollapsedOrthologyMatrix::Entry ie=childCOM(newJ,insideLeaf);
    int iBegin=ie.begin, iEnd=ie.end;
    if(iEnd>ie.begin) --iEnd; // because for indel states it's noninclusive
    if(iEnd>taxa[insideLeaf].getSeqLen()) --iEnd;
    if(iBegin<0) iBegin=0;
    cout<<"d "<<taxa[insideLeaf].getName()<<" "<<taxa[insideLeaf].getSeqLen()<<" "<<iBegin<<" "<<iEnd<<endl;
    for(int ol=0 ; ol<numOutsideLeaves ; ++ol) {
      int outsideLeaf=outsideLeaves[ol];
      SparseMatrix3D &M=*matrices[insideLeaf][outsideLeaf];
      float pathLen=tree->distanceBetweenLeaves(taxa[insideLeaf].getID(),
						taxa[outsideLeaf].getID());
      CollapsedOrthologyMatrix::Entry oe=parentCOM(newI,outsideLeaf);
      Set<PHMM_StateType> stateTypes;
      bool INS=(stateType==PHMM_INSERT || ie.end>ie.begin || 
		ie.begin<0 || ie.end>taxa[insideLeaf].getSeqLen());
      bool DEL=(stateType==PHMM_DELETE || oe.end>oe.begin || 
		oe.begin<0 || oe.end>taxa[outsideLeaf].getSeqLen());
      if(!INS && !DEL) stateTypes.insert(PHMM_MATCH);
      else {
	if(INS) stateTypes.insert(PHMM_INSERT);
	if(DEL) stateTypes.insert(PHMM_DELETE);
      }
      int oBegin=oe.begin, oEnd=oe.end;
      if(oEnd>oe.begin) --oEnd; // because for indel states it's noninclusive
      if(oEnd>taxa[outsideLeaf].getSeqLen()) --oEnd;
      if(oBegin<0) oBegin=0;
      cout<<"e "<<taxa[outsideLeaf].getName()<<" "<<taxa[outsideLeaf].getSeqLen()<<" "<<oBegin<<" "<<oEnd<<" "<<INS<<" "<<DEL<<endl;
      Vector<float> sumW;
      for(int ipos=iBegin ; ipos<=iEnd ; ++ipos) {
	cout<<"f "<<ipos<<endl;
	Set<PHMM_StateType>::iterator cur=stateTypes.begin(), 
	  end=stateTypes.end();
	for(; cur!=end ; ++cur) {
	  PHMM_StateType stateType=*cur;
	  EntryList &row=M(ipos,stateType);
	  EntryList::iterator cur=row.begin(), end=row.end();
	  for(; cur!=end ; ++cur) {
	    SparseMatrix3D::Entry e=*cur;
	    if(e.y<oBegin) continue;
	    if(e.y>oEnd) break;
	    sumW.push_back(e.value);
	    cout<<"push "<<e.value<<endl;
	  }
	}
      }
      float theSum=sumLogProbs(sumW);
      cout<<"sum="<<theSum<<" pathLen="<<pathLen<<endl;
      if(isFinite(theSum)) sumV.push_back(theSum-log(pathLen));
    }
  }
  if(sumV.size()==0) {
    justPseudo=true;
    return NEGATIVE_INFINITY;
  }
  else justPseudo=false;
  float ave=sumLogProbs(sumV)-log(sumV.size());
  cout<<"ave="<<ave<<endl;
  //float ave=PowerMean::compute_log(sumV,5);//1.5);
  return ave;
}



/****************************************************************
                     Application::manhattan()
 ****************************************************************/
bool Application::manhattan(StatePath &path,int &i,int &j,
			     Vector<int> &insideLeaves,
			     Vector<int> &outsideLeaves,int numInsideLeaves,
			     int numOutsideLeaves,
			     CollapsedOrthologyMatrix &childCOM,
			     CollapsedOrthologyMatrix &parentCOM,
			     int parentL,int childL,
			     Vector<bool> &parentReachability,
			     Vector<bool> &childReachability)
{
  bool justPseudo;
  int newI=i+2, newJ=j+2;
  float insP=NEGATIVE_INFINITY, delP=NEGATIVE_INFINITY;
  for(; newI<=parentL || newJ<=childL ; ++newI, ++newJ) {
    for(STATE q=1 ; q<4 ; ++q) {
      if(newI<parentL && j<childL) {
	insP=computeCellP(insideLeaves,outsideLeaves,numInsideLeaves, 
			  numOutsideLeaves,childCOM,parentCOM,q,
			  i,newJ,justPseudo,parentReachability,
			  childReachability);
	if(justPseudo) insP=NEGATIVE_INFINITY;
	if(isFinite(insP)) goto BREAK;
      }
      if(newJ<childL && i<parentL) {
	delP=computeCellP(insideLeaves,outsideLeaves,numInsideLeaves, 
			  numOutsideLeaves,childCOM,parentCOM,q,
			  newI,j,justPseudo,parentReachability,
			  childReachability);
	if(justPseudo) delP=NEGATIVE_INFINITY;
	if(isFinite(delP)) goto BREAK;
      }
    }
  }
 BREAK:
  int offset=newI-i-1;//### added -1
  if(isFinite(insP) && isFinite(delP)) { // RESOLVE TIES
    float ratio=exp(insP-sumLogProbs(insP,delP));
    if(Random0to1()>ratio) insP=NEGATIVE_INFINITY;
  }
  if(isFinite(insP)) { // HORIZONTAL
    for(int k=0 ; k<offset ; ++k) path.push_back(insertState);
    j=newJ-1;
    return true;
  }
  else if(isFinite(delP)) { // VERTICAL
    for(int k=0 ; k<offset ; ++k) path.push_back(deleteState);
    i=newI-1;
    return true;
  }
  cout<<"BAILING OUT "<<i<<" "<<j<<endl;
  if(Random0to1()<0.5) { // ABORT -- INSERT/DELETE ALL IN BOTH SEQS
    for(; i<parentL ; ++i) path.push_back(deleteState);
    for(; j<childL ; ++j) path.push_back(insertState);
    return true;
  }
  else {
    for(; j<childL ; ++j) path.push_back(insertState);
    for(; i<parentL ; ++i) path.push_back(deleteState);
    return true;
  }

  return false;
}
/*
  bool justPseudo;
  for(int offset=1 ; i+offset<=parentL || j+offset<=childL ; ++offset) {
    for(int il=0 ; il<numInsideLeaves ; ++il) {
      int insideLeaf=insideLeaves[il];
      CollapsedOrthologyMatrix::Entry ie=childCOM(j+offset,insideLeaf);
      int iBegin=ie.begin, iEnd=ie.end;
      if(iEnd>ie.begin) --iEnd; // because for indel states it's noninclusive
      if(iEnd>taxa[insideLeaf].getSeqLen()) --iEnd;
      if(iBegin<0) iBegin=0;
      for(int ol=0 ; ol<numOutsideLeaves ; ++ol) {
	int outsideLeaf=outsideLeaves[ol];
	SparseMatrix3D &M1=*matrices[insideLeaf][outsideLeaf];
	SparseMatrix3D &M2=*matrices[outsideLeaf][insideLeaf];
	      CollapsedOrthologyMatrix::Entry oe=parentCOM(newI,outsideLeaf);

      }
    }
  }
======
*/



/****************************************************************
                    Application::jumpToCell()
 ****************************************************************/
void Application::jumpToCell(StatePath &path,int &i,int &j,
			     int newI,int newJ) {
  if(Random0to1()<0.5) {
    for(int k=i+1 ; k<newI ; ++k) path.push_back(deleteState);
    for(int k=j+1 ; k<newJ ; ++k) path.push_back(insertState);
  }
  else {
    for(int k=j+1 ; k<newJ ; ++k) path.push_back(insertState);
    for(int k=i+1 ; k<newI ; ++k) path.push_back(deleteState);
  }
}



/****************************************************************
                     Application::sampleBranch2()
 ****************************************************************/
void Application::sampleBranch2(PhylogenyBranch *br)
{
  //verbose("samapleBranch2");
  BranchAttributes *branch=
    static_cast<BranchAttributes*>(br->getDecoration());
  getHMMprobs(branch->getHMM());
  Taxon &parent=branch->getParentTaxon();
  Taxon &child=branch->getChildTaxon();
  CollapsedOrthologyMatrix &parentCOM=*parent.getCOM();
  CollapsedOrthologyMatrix &childCOM=*child.getCOM();
  int parentL=parent.getSeqLen(), childL=child.getSeqLen();
  Vector<int> outsideLeaves, insideLeaves;
  getLeaves(insideLeaves,outsideLeaves,parent,child);
  int numInsideLeaves=insideLeaves.size();
  int numOutsideLeaves=outsideLeaves.size();
  StatePath path(hmm);
  Vector<bool> parentReachability, childReachability;
  //verbose("a");
  computeReachabilityVector(parentCOM,outsideLeaves,parentReachability);
  computeReachabilityVector(childCOM,insideLeaves,childReachability);
  //verbose("b");
  PosteriorBackward B(*hmm,parent,child,numTaxa,insideLeaves,outsideLeaves,
		      childCOM,parentCOM,parentReachability,
		      childReachability,matrices,logPseudoCount,
		      taxa,tree,wantTransitivity,indelCoef);
  //verbose("c");
  int i=0, j=0; // 1-based
  STATE prevQ=0;
  while(i<parentL || j<childL) {
    Vector<STATE> states;
    Vector<float> probs;
    for(int q=1 ; q<numStates ; ++q) {
      PHMM_StateType stateType=stateTypes[q];
      int newI=i+dI[q], newJ=j+dJ[q];
      if(newI>parentL || newJ>childL) continue;
      bool justPseudo;

      /*
      float cellP=computeCellP(insideLeaves,outsideLeaves,numInsideLeaves,
			       numOutsideLeaves,childCOM,parentCOM,q,newI,
			       newJ,justPseudo,parentReachability,
			       childReachability);
      */
      //float cellP=B.getCachedEmitP(q,newI,newJ);
      float cellP=B.getEmitP(q,newI,newJ);

      //if(justPseudo) continue;// ### 2/8
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
  //verbose("updateIndels");
  updateIndels(*branch,parent,child,path);
  //verbose("updateTentacles");
  updateTentacles(*branch,parent,child,parentCOM,childCOM,parentL,childL,
		  outsideLeaves,insideLeaves,numOutsideLeaves,numInsideLeaves,
		  path);
  cout.flush();
  //verbose("/sampleBranch2");
}



