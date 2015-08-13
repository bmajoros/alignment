/****************************************************************
 AVES.C : Alignment Via Evolutionary Sampling
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "AVES.H"
#include "BOOM/Stack.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/MPI.H"
#include "SiblingHirschberg.H"
#include "FunctionalDollo.H"
#include "HirschPosteriors.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;

extern Alphabet alphabet;
Alphabet alphabet;


/****************************************************************
                          AVES::AVES()
 ****************************************************************/
AVES::AVES()
  : alphabetMap(DnaDashAlphabet::global(),PureDnaAlphabet::global()),
    functionalParse(NULL), targetSpecies(NULL), alignmentBuilder(NULL),
    fullAlignment(NULL), refinementIterations(0), posteriors(NULL),
    identityMap(PureDnaAlphabet::global()), shouldDumpModel(false),
    contentSensor(NULL), shouldSeed(false), baselineMode(false),
    usePosteriors(false), shouldHillClimb(false), wantFuncDollo(false),
    explicitHistories(true), useHirschberg(true), useContentSensor(false),
    usePrescan(true), noPrescanOnDownpass(true)
{
  ::alphabet=PureDnaAlphabet::global();
  gapSymbol=alphabet.lookup('-');
  modelCompiler=new ModelCompiler;
  modelCompiler->setGCthreshold(100);//10000000);
  templateInstantiator=modelCompiler->getInstantiator();
  isNucleotide.resize(alphabet.size());
  isNucleotide.setAllTo(false);
  isNucleotide[alphabet.lookup('A')]=true;
  isNucleotide[alphabet.lookup('C')]=true;
  isNucleotide[alphabet.lookup('G')]=true;
  isNucleotide[alphabet.lookup('T')]=true;
}



/****************************************************************
                          AVES::main()
 ****************************************************************/
int AVES::main(int argc,char *argv[])
{
  mpi=new MPI(argc,&argv);
  slaveDriver=new MpiSlaveDriver(*mpi);
  return slaveDriver->amItheMaster() ?
    master(argc,argv) : slave(argc,argv);
}




int AVES::master(int argc,char *argv[]) 
{
  cout<<"numSlaves="<<slaveDriver->getNumSlaves()<<endl;
  processCmdLine(argc,argv);

  if(!shouldSeed) {
    // Perform alignment
    cout<<"performing alignment"<<endl;
    timer.startCounting();
    progressive();
    timer.stopCounting();
    cout<<"Elapsed time: "<<timer.elapsedTime()<<endl;
      MemoryProfiler::report("Total memory used:");
  }
  
  if(!baselineMode) reconstructAlignment(true);
  fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,true);

  // Apply functional Dollo, to localize binding sites to lower clades
  if(wantFuncDollo) {
    cout<<"Localizing sites via functional Dollo..."<<endl;
    localize();
    if(!baselineMode) reconstructAlignment(true);
  }

  // Emit results
  cout<<"\n\nFINAL ALIGNMENT:\n"<<endl;
  cout<<"length="<<fullAlignment->getLength()<<endl;
  emit(outfile);

  cout<<"done"<<endl;
  cout<<"terminating slaves"<<endl;
  slaveDriver->terminateSlaves();
  return 0;
}



int AVES::slave(int argc,char *argv[]) 
{
  processCmdLine(argc,argv);
  LossyFelsenstein F(pureDnaAlphabet,identityMap,numTaxa);
  bool askedForWork=false;
  while(1) {
    MpiFixedMsg *msg;
    if(askedForWork) msg=mpi->waitForMessage();
    else {
      msg=slaveDriver->acceptWork();
      askedForWork=true;
    }
    //cout<<"slave "<<mpi->getProcessID()<<" recd msg="<<msg<<endl;
    switch(msg->getTag())  
      {
      case TERMINATE_SLAVE: return 0;
      case MSG_ALIGN_SIBS: {
	//cout<<"slave "<<mpi->getProcessID()<<" recd: MSG_ALIGN_SIBS"<<endl;
	int parentID;
	(*msg)>>parentID;
	Taxon &parent=taxa[parentID];
	align(parent.getIthChild(LEFT).getNode(),
	      parent.getIthChild(RIGHT).getNode(),F);
	askedForWork=false;
      }	break;
      case MSG_STATE_PATH: {
	int leftID, rightID;
	StatePath *path=unpackPath(*msg,leftID,rightID);
	//cout<<"slave "<<mpi->getProcessID()<<" recd: MSG_STATE_PATH "<<msg->getBufferSize()<<" "<<taxa[leftID].getName()<<" "<<taxa[rightID].getName()<<endl;
	installStatePath(*path,leftID,rightID);
	Taxon &left=taxa[leftID];
	Taxon &parent=*left.getParent();
	//cout<<"slave "<<mpi->getProcessID()<<" now "<<parent.getName()<<" length="<<parent.getSeqLen()<<endl;
	delete path;
      } break;
      }
    delete msg;
  }
}



int AVES::processCmdLine(int argc,char *argv[]) 
{
  // Process command line
  CommandLine cmd(argc,argv,"S:Ds:hgi:md:b:r:p:PBl:Avt:e:ZGH:");
  if(cmd.numArgs()!=5)
    throw String(
"AVES [options] <in.fasta> <in.lambda> <target-track> <out.maf> <out.gff>\n\
     where: -A = just align, don't try to annotate\n\
            -b <strategy> = band using strategy:\n\
                f### = fixed-width banding of width 2*###\n\
            -B = baseline (background align, then annotate w/PhyloHMM)\n\
            -d <outfile> = dump model parameters from merged model\n\
            -D = apply Dollo parsimony to functional classes\n\
            -e <*.model> = use given content sensor for EQ frequencies\n\
            -g = use Gibbs sampling instead of Metropolis-Hastings\n\
            -G = don't perform garbage collection\n\
            -h = perform stochastic hill-climbing instead of sampling\n\
            -H <outfile> = emit site histories into outfile\n\
            -i <*.maf> = initialize with the given alignment\n\
            -l <X> = rescale branch lengths by factor X\n\
            -m = marginalize over indel histories\n\
            -p <in.gff> = use pre-scanned sites to limit search space\n\
            -P = use Posterior Viterbi for progressive stage\n\
            -r <N> = use N iterations of iterative refinement\n\
            -s <N> = sample N alignments\n\
            -S LLR|POSTERIOR|NONE = scoring method (default=POSTERIOR)\n\
            -t N = use N threads\n\
            -v = use standard Viterbi, rather than Hirschberg\n\
            -Z = ignore pre-scans during downpass\n\
");
  String fastaFile=cmd.arg(0);
  String lambdaFile=cmd.arg(1);
  String targetName=cmd.arg(2);
  outfile=cmd.arg(3);
  if(cmd.option('S')) {
    String strategy=cmd.optParm('S');
    if(strategy=="LLR") scoringStrategy=SS_LLR;
    else if(strategy=="POSTERIOR") scoringStrategy=SS_POSTERIOR;
    else if(strategy=="NONE") scoringStrategy=SS_NONE;
    else throw String("Unknown scoring strategy: ")+strategy;
  }
  else scoringStrategy=SS_POSTERIOR;
  gffFile=cmd.arg(4);
  baselineMode=cmd.option('B');
  if(cmd.option('r'))
    refinementIterations=cmd.optParm('r').asInt();
  usePosteriors=cmd.option('P');
  shouldSample=cmd.option('s');
  numSamples=cmd.optParm('s');
  shouldHillClimb=cmd.option('h');
  shouldSeed=cmd.option('i');
  wantFuncDollo=cmd.option('D');
  samplerType=(cmd.option('g') ? GIBBS : METROPOLIS_HASTINGS);
  //if(samplerType==GIBBS) throw "option -g is currently disabled";
  explicitHistories=!cmd.option('m');
  randomize();
  //Lambda::GarbageCollector::setThreshold(10000000); // ###
  shouldDumpModel=cmd.option('d');
  useHirschberg=!cmd.option('v');
  dumpFile=cmd.optParm('d');
  useContentSensor=cmd.option('e');
  if(cmd.option('b')) {
    String strategy=cmd.optParm('b');
    if(strategy[0]=='f' && strategy.length()>1) {
      bandwidth=strategy.substring(1,strategy.length()-1).asInt();
      bandingType=FIXED_WIDTH_BANDING;
    }
    else throw "invalid banding strategy specified";
  }
  else bandingType=NO_BANDING;
  usePrescan=cmd.option('p');
  String prescanFile=(usePrescan ? cmd.optParm('p') : "");
  numThreads=cmd.option('t') ? cmd.optParm('t').asInt()-1 : 0;
  cout<<(numThreads+1)<<" threads"<<endl;
  noPrescanOnDownpass=cmd.option('Z');
  disableGC=cmd.option('G');
  if(disableGC) modelCompiler->setGCthreshold(LARGEST_INTEGER);
  if(cmd.option('H')) historiesFile=cmd.optParm('H');

  // Load content sensor for higher-order equilibrium frequencies
  contentSensor=
    useContentSensor ?
    ContentSensor::load(cmd.optParm('e')) :
    NULL;
  //cout<<"contentSensor="<<contentSensor<<endl;
  
  // Load HMM
  cout<<"executing "<<lambdaFile<<"..."<<endl;
  LambdaAPI &lambda=modelCompiler->getLambda();
  modelCompiler->parse(lambdaFile);
  lambda.checkGarbageLevel();

  cout<<modelCompiler->getNumModels()<<" models loaded"<<endl;
  //if(modelCompiler->getNumModels()<1)
  //  throw "use (register-transducer ...) to register a transducer";
  //transducerTemplate=modelCompiler->getIthModel(0);
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  transducerTemplate=modelCompiler->instantiate(ctor);
  cout<<transducerTemplate->getNumStates()<<" states"<<endl;
  defaultHMM=templateInstantiator->instantiate(transducerTemplate,1.0,NULL);

  // Get phylogeny
  cout<<"loading tree"<<endl;
  tree=modelCompiler->getGuideTree(transducerTemplate);
  tree->gatherNodes(phylogenyNodes);
  if(cmd.option('l'))
    tree->scaleBranches(cmd.optParm('l').asFloat());

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
    if(node->getName()==targetName) targetSpecies=&taxon;
  }
  if(!targetSpecies) throw "Target species not found";
  alignmentBuilder=
    new AlignmentBuilder(tree->getRoot(),alphabet,gapSymbol,numTaxa);

  cout<<"gatherCladeMembers"<<endl;
  gatherCladeMembers();

  // Instantiate BranchHMM's and substitution matrices on tree
  cout<<"instantiating HMM's"<<endl;
  initBranches();

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
    
    // Load pre-scanned annotations to limit search space
    cout<<"loading pre-scans"<<endl;
    for(int i=0 ; i<numTaxa ; ++i) {
      Taxon &taxon=taxa[i];
      Array1D<FuncClassSet> &fcConstraints=taxon.getFcConstraints();
      fcConstraints.resize(taxon.getSeqLen());
    }
    if(usePrescan) loadFcConstraints(prescanFile);

    cout<<"reconstructHistories"<<endl;
    reconstructHistories();
    if(usePrescan) {
      cout<<"rebuilding prescan constraints"<<endl;
      rebuildPrescans(tree);
    }
    baselineMode=true;//###
    cout<<"inferring root's functional parse"<<endl;
    GainLossFelsenstein F2(pureDnaAlphabet,identityMap,numTaxa);
    inferRootFuncParse(F2);

    //###
    //reconstructAlignment(true);//###???
    //###

    cout<<"downpass()"<<endl;
    downpass(F2);
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

      // Load pre-scanned annotations to limit search space
      cout<<"loading pre-scans"<<endl;
      for(int i=0 ; i<numTaxa ; ++i) {
	Taxon &taxon=taxa[i];
	Array1D<FuncClassSet> &fcConstraints=taxon.getFcConstraints();
	fcConstraints.resize(taxon.getSeqLen());
      }
      if(usePrescan) loadFcConstraints(prescanFile);
    }
  return 0;
}



/****************************************************************
                             AVES::sample()
 ****************************************************************/
void AVES::sample()
{
  // First, select a suitable branch
  int numNodes=phylogenyNodes.size();
  int childID=RandomNumber(numNodes);
  PhylogenyNode *parent=taxa[childID].getNode()->getParent();
  while(!parent || parent->getNodeType()!=INTERNAL_NODE) {
    childID=(childID+1)%numNodes;
    parent=taxa[childID].getNode()->getParent();
  }
  Taxon &childTaxon=taxa[childID], &parentTaxon=taxa[parent->getID()];
  cout<<"re-sampling above clade: "<<childTaxon.getName()<<" ... ";

  // Get the Branch HMM for that branch
  InternalNode *internalParent=static_cast<InternalNode*>(parent);
  WhichChild whichChild=
    internalParent->whichChildIsThis(childTaxon.getNode());
  BranchAttributes *branch=parentTaxon.getIthBranch(whichChild);
  BranchHMM *hmm=branch->getHMM();

  // Resample the pairwise alignment along that branch
  cout<<"Running the backward-algorithm"<<endl;
  LinkFelsenstein F(alphabet,alphabetMap,numTaxa);
  LinkBackward B(*hmm,parentTaxon,childTaxon,alphabetMap,numTaxa,F);
  cout<<"Sampling pairwise alignment..."<<endl;
  LinkSampler sampler(*hmm,B,parentTaxon,childTaxon,branch);
  double P_old_given_new, P_new_given_old;
  StatePath *newPath=sampler.samplePath(P_new_given_old,P_old_given_new);
  StatePath *oldPath=branch->getStatePath();
  branch->setStatePath(newPath);
  IndexMap oldUpMap=branch->getUpMap(), oldDownMap=branch->getDownMap();
  reconstructHistory(*newPath,*hmm,branch->getUpMap(),branch->getDownMap());
  double P_old=alignmentScore, P_new=likelihood();

  // Decide whether to accept or reject
  bool accept;
  if(shouldHillClimb) accept=(P_new>P_old);
  else {
    double logHastings=P_new - P_old + P_old_given_new - P_new_given_old;
    double hastings=exp(logHastings);
    if(samplerType==METROPOLIS_HASTINGS) 
      cout<<"hastings ratio="<<hastings<<" ("<<P_new<<"-"<<P_old
	  <<"+"<<P_old_given_new<<"-"<<P_new_given_old<<" = "
	  <<logHastings<<")"<<endl;
    double r=Random0to1();
    accept=(r<=hastings || samplerType==GIBBS);
    cout<<(accept ? "ACCEPT" : "REJECT")<<endl;
  }

  // If accepted, install the new alignment
  if(accept) {
    alignmentScore=P_new;
    delete oldPath;
    parentTaxon.getFunctionalParse()=FunctionalParse(*newPath,PARENT);
    childTaxon.getFunctionalParse()=FunctionalParse(*newPath,CHILD);

   //reconstructHistory(*newPath,*hmm,branch->getUpMap(),branch->getDownMap());

      // ### debugging:
      reconstructAlignment();
      fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,true);
  }
  else {
    branch->setStatePath(oldPath);
    branch->getUpMap()=oldUpMap;
    branch->getDownMap()=oldDownMap;
    delete newPath;
  }
}



/****************************************************************
                         AVES::likelihood()
 ****************************************************************/
double AVES::likelihood()
{
  // First, compute the contribution of the emission terms by applying
  // Felsenstein's algorithm to all connected components in the graph
  double logP=0;
  LinkFelsenstein F(alphabet,alphabetMap,numTaxa);
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    int L=taxon.getSeqLen();
    BranchAttributes *parentBranch=taxon.getBranchToParent();
    if(parentBranch) {
      IndexMap &upMap=parentBranch->getUpMap();
      for(int j=0 ; j<L ; ++j) {
	if(upMap[j]==IndexMap::UNDEFINED)// residue is a root
	  {
	    logP+=F.logLikelihood(j,taxon,false);
	    if(!finite(logP)) 
	      {cout<<taxon.getName()<<" "<<j<<endl; INTERNAL_ERROR;}
	  }
      }
    }
    else for(int j=0 ; j<L ; ++j) {
      logP+=F.logLikelihood(j,taxon,false);
      if(!finite(logP)) 
	{cout<<taxon.getName()<<" "<<j<<endl; INTERNAL_ERROR;}
    }
  }

  // Next, process each branch individually to assess the transition 
  // probabilities for the transducers
  logP+=sumTransLogProbs();

  return logP;
}



/****************************************************************
                      AVES::sumTransProbs()
 ****************************************************************/
double AVES::sumTransLogProbs()
{
  double logP=0;
  int n=branches.size();
  for(int i=0 ; i<n ; ++i) {
    BranchAttributes &branch=branchAttributes[i];
    BranchHMM *hmm=branch.getHMM();
    StatePath *path=branch.getStatePath();
    if(!path) INTERNAL_ERROR;
    logP+=hmm->getTransProbs(*path);
  }
  return logP;
}


/****************************************************************
                    AVES::reconstructAlignment()

 Reconstructs an explicit alignment from the "linked-residues"
 representation on the phylogeny.
 ****************************************************************/
void AVES::reconstructAlignment(bool includeAncestors)
{
  delete fullAlignment;
  fullAlignment=alignmentBuilder->buildAlignment(includeAncestors,
						 includeAncestors);
  //SingleGainParsimony P(*tree,alphabet,*fullAlignment,
  //		   fullAlignment->getGapSymbols());
  //P.run();
}

 

void AVES::loadFcConstraints(const String &filename)
{
  // Read features
  GffReader reader(filename);
  GffFeature *feature;
  while(feature=reader.nextFeature()) {
    const String &substrate=feature->getSubstrate();
    String typeName=feature->getFeatureType();
    int begin=feature->getBegin()+1; // +1 is because GffReader subtracts 1
    char strand=feature->getStrand();
    if(strand=='-') typeName=String("-")+typeName;
    if(!nameToTaxon.isDefined(substrate)) continue;
    int taxonID=nameToTaxon[substrate];
    Taxon &taxon=taxa[taxonID];
    //cout<<taxon.getName()<<"="<<substrate<<endl;//###
    Array1D<FuncClassSet> &fcConstraints=taxon.getFcConstraints();
    FunctionalElementType type=
      FunctionalElementType::getTypeByName(typeName);
    const Vector<FunctionalClass> &classes=type.getClasses();
    int n=classes.size();
    if(n==0) {
      //cout<<"\""<<typeName<<"\""<<" -> "<<type<<endl;
      //INTERNAL_ERROR; // ### DEBUGGING
    }
    for(int i=0 ; i<n ; ++i) 
      fcConstraints[begin+i].push_back(classes[i]);
    delete feature;
  }
   
  // Sort arrays (so we can later use fast intersect/union algorithms)
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    Array1D<FuncClassSet> &fcConstraints=taxon.getFcConstraints();
    int n=fcConstraints.size();
    for(int i=0 ; i<n ; ++i) fcConstraints[i].sortInPlace();
  }
}



/****************************************************************
                         AVES::progressive()
 ****************************************************************/
void AVES::progressive()
{
  // Phase I : The UP-PASS
  LossyFelsenstein F1(pureDnaAlphabet,identityMap,numTaxa);
  //GainLossFelsenstein F1(pureDnaAlphabet,identityMap,numTaxa);//###DEBUGGING
  cout<<"uppass()"<<endl;
  uppass2(F1);
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
  if(refinementIterations>0) {
    cout<<"iterativeRefinement()"<<endl;
    iterativeRefinement(F2);
  }

  //dollo();
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
void AVES::uppass2(LinkFelsenstein &F)
{
  BitSet resolved(numTaxa);
  for(int i=0 ; i<numTaxa ; ++i) 
    if(taxa[i].isLeaf()) resolved.addMember(i);
  while(resolved.cardinality()<numTaxa) {
    for(int i=0 ; i<numTaxa ; ++i) {
      if(resolved.isMember(i)) continue;
      Taxon &parent=taxa[i];
      Taxon &left=parent.getIthChild(LEFT), &right=parent.getIthChild(RIGHT);
      if(!resolved.isMember(left.getID()) ||
	 !resolved.isMember(right.getID())) continue;
      MpiVariableMsg *msg=new MpiVariableMsg(MSG_ALIGN_SIBS);
      (*msg)<<i;
      msg->close();
      slaveDriver->addWork(msg);
    }
    Vector<MpiFixedMsg*> results;
    //cout<<"master waiting for results"<<endl;
    slaveDriver->waitForResults(results);
    //cout<<"master install state paths"<<endl;
    installStatePaths(results,resolved);
  }
}
void AVES::installStatePaths(Vector<MpiFixedMsg*> &msgs,BitSet &resolved)
{
  int n=msgs.size();
  for(int i=0 ; i<n ; ++i) {
    MpiFixedMsg *msg=msgs[i];
    //cout<<"master msg="<<msg<<" "<<" tag="<<msg->getTag()<<endl;
    int leftID, rightID;
    StatePath *path=unpackPath(*msg,leftID,rightID);
    //cout<<"master path="<<path<<endl;
    delete msg;
    //cout<<"master install"<<endl;
    installStatePath(*path,leftID,rightID);
    //cout<<"master delete path"<<endl;
    delete path;
    //cout<<"master end of loop"<<endl;
    resolved.addMember(taxa[leftID].getParent()->getID());
  }
}
void AVES::installStatePath(StatePath &path,int leftID,int rightID)
{
  Taxon &left=taxa[leftID], &right=taxa[rightID];
  //cout<<"process "<<mpi->getProcessID()<<" : "<<left.getName()<<"-"<<right.getName()<<endl;
  Taxon &parent=*left.getParent();
  //cout<<"process "<<mpi->getProcessID()<<" constructing functional parses"<<endl;

  LossyFelsenstein F(pureDnaAlphabet,identityMap,numTaxa);
  int numFC=FunctionalClass::numClasses();
  int numAlpha=pureDnaAlphabet.size();
  left.getPrecomputedEmissions().resize(left.getSeqLen(),numFC,numAlpha);
  F.precomputeInside(left);
  right.getPrecomputedEmissions().resize(right.getSeqLen(),numFC,numAlpha);
  F.precomputeInside(right);

  left.getFunctionalParse()=FunctionalParse(path,PARENT);
  right.getFunctionalParse()=FunctionalParse(path,CHILD);
  //cout<<"process "<<mpi->getProcessID()<<" install indel histories"<<endl;
  BranchAttributes *leftBranch=left.getBranchToParent();
  BranchAttributes *rightBranch=right.getBranchToParent();
  installSiblingHistories(path,parent,leftBranch,rightBranch,defaultHMM);
  if(usePrescan) updatePrescans(parent);
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
                         AVES::packPath()
 ****************************************************************/
MpiVariableMsg *AVES::packPath(StatePath &path,MPI_MESSAGE_TAG tag,
			       int leftID,int rightID)
{
  MpiVariableMsg *msg=new MpiVariableMsg(tag);
  int n=path.length();
  (*msg)<<leftID<<rightID<<n;
  for(int i=0 ; i<n ; ++i) {
    //if(path[i]>872) INTERNAL_ERROR;
    (*msg)<<path[i];
  }
  msg->close();
  //cout<<"sending path of length "<<n<<" : "<<path<<endl;
  return msg;
}



/****************************************************************
                         AVES::unpackPath()
 ****************************************************************/
StatePath *AVES::unpackPath(MpiFixedMsg &msg,int &leftID,int &rightID)
{
  StatePath *path=new StatePath;
  msg.beginExtraction();
  int n;
  msg>>leftID>>rightID>>n;
  //cout<<"receiving path: "<<&msg<<" "<<leftID<<" "<<rightID<<" "<<n<<" "<<msg.getBufferSize()<<" pid="<<mpi->getProcessID()<<endl;
  STATE q;
  for(int i=0 ; i<n ; ++i) {
    msg>>q;
    path->push_back(q);
    //cout<<q<<endl;
    //if(q>872) INTERNAL_ERROR;
  }
  path->setHMM(defaultHMM);
  return path;
}



/****************************************************************
                           AVES::align()
 ****************************************************************/
void AVES::align(PhylogenyNode *leftNode,PhylogenyNode *rightNode,
		 LinkFelsenstein &F)
{
  //cout<<"align()"<<endl;

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
  //cout<<"instantiate hmm"<<endl;
  BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
						   combinedBranchLengths,
						   NULL);
  // Perform alignment
  //cout<<"instantiating viterbi"<<endl;
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
  //cout<<"decode() "<<left.getSeqLen()<<" "<<right.getSeqLen()<<endl;
  StatePath &path=*viterbi.decode();
  //cout<<"constructing functional parses"<<endl;
  left.getFunctionalParse()=FunctionalParse(path,PARENT);
  right.getFunctionalParse()=FunctionalParse(path,CHILD);
  cout<<"path score="<<path.getScore()<<endl;

  //cout<<"install indel histories"<<endl;
  installSiblingHistories(path,parent,leftBranch,rightBranch,hmm);

  int myID=mpi->getProcessID(), numProcesses=mpi->getNumProcesses();
  MpiVariableMsg *msg=
    packPath(path,MSG_STATE_PATH,left.getID(),right.getID());
  for(int i=0 ; i<numProcesses ; ++i) {
    if(i==myID) continue;
    //cout<<"slave "<<mpi->getProcessID()<<" sending path to "<<i<<endl;
    mpi->send(*msg,i);
  }
  delete msg;
  delete &path;
  delete hmm;
  if(usePosteriors) { delete forward; delete backward; }

  if(usePrescan) {
    //cout<<"update prescan constraints"<<endl;
    updatePrescans(parent);
  }

  if(useContentSensor) {
    //cout<<"reconstructing ancestral sequence for content sensing"<<endl;
    parent.getSeq().resize(parent.getSeqLen());
    ResidueAddress ra(&parent,0);
    LinkParsimony parsimony(ra,pureDnaAlphabet,gapSymbol,numTaxa);
    parsimony.runFullSeq();
    parent.getSeqStr()=parent.getSeq()(alphabet);
  }
  //delete &viterbi;
}



/****************************************************************
                     AVES::updatePrescans()
 ****************************************************************/
void AVES::updatePrescans(Taxon &parent)
{
  //cout<<"updatePrescans "<<parent.getName()<<" "<<parent.getSeqLen()<<endl;
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

        //###
      /*
	cout<<left.getName()<<" "<<leftPos<<": ";
	FuncClassSet &leftSet=leftConstraints[leftPos];
	FuncClassSet::iterator cur=leftSet.begin(), end=leftSet.end();
	for(; cur!=end ; ++cur) cout<<*cur<<" ";
	cout<<right.getName()<<" "<<rightPos<<": ";
	FuncClassSet &rightSet=rightConstraints[rightPos];
	cur=rightSet.begin(), end=rightSet.end();
	for(; cur!=end ; ++cur) cout<<*cur<<" ";
      */
	//cout<<endl;
	//###

      parentSet.clear();
      leftConstraints[leftPos].unionWith(rightConstraints[rightPos],
					 parentSet);
    }

    // ### DEBUGGING
    /*
    if(left.isLeaf()) {
      //cout<<"TAXON="<<parent.getName()<<endl;
      cout<<pos<<": ";
      FuncClassSet::iterator cur=parentSet.begin(), end=parentSet.end();
      for(; cur!=end ; ++cur) cout<<*cur<<" ";
      cout<<endl;
    }
    */
    // ###

  }
}



/****************************************************************
                     AVES::rebuildPrescans()
 ****************************************************************/
void AVES::rebuildPrescans(Phylogeny *tree)
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
                       AVES::functionalDollo()
 ****************************************************************/
void AVES::localize()
{
  FunctionalDollo dollo(*tree);
  dollo.run();
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
    double LL;
    if(scoringStrategy==SS_POSTERIOR) LL=computePostDenominator(taxon);
    for(int i=0 ; i<n ; ++i) {
      FunctionalElement E=elems[i];
      String substrate=taxon.getName(), source="AVES";
      String type=E.getType().getName();
      int begin=E.getBegin(), end=E.getEnd(), phase=0;
      double score;//=taxon.isLeaf() ? computeScore(taxon,E) : 0.0;
      //exp(computePostNumerator(taxon,E)-LL);
      if(scoringStrategy==SS_POSTERIOR) 
	score=exp(computePostNumerator(taxon,E)-LL);
      else score=computeScore(taxon,E);
      Strand strand=E.getStrand();
      if(type[0]=='-') {
	strand='-';
	type=type.substring(1);
      }
      --end; // ###
      if(taxon.isLeaf()) {
	if(scoringStrategy==SS_POSTERIOR) {
	  double LLR=computeScoreLLR(taxon,E);
	  gff<<substrate<<"\t"<<source<<"\t"<<type<<"\t"<<begin<<"\t"<<end
	     <<"\t"<<score<<"\t"<<strand<<"\t"<<phase<<" /LLR="<<LLR<<endl;
	}
	else 
	  gff<<substrate<<"\t"<<source<<"\t"<<type<<"\t"<<begin<<"\t"<<end
	     <<"\t"<<score<<"\t"<<strand<<"\t"<<phase<<endl;
      }
      else
	gff<<substrate<<"\t"<<source<<"\t"<<type<<"\t"<<begin<<"\t"<<end
	   <<"\t"<<score<<"\t"<<strand<<"\t"<<phase<<endl;
    }
    taxonElems[i]=&elems;
    //delete &elems;
  }
  gff.close();

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
	continue;
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
  //for(int i=0 ; i<numTaxa ; ++i) delete taxonElems[i];
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
  tree->gatherCladeMembers();
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
	  //newTrack[j]=(oldTrack[j]==gapSymbol ? gap : residue);
	  newTrack[j]=(oldTrack[j]>3 ? gap : residue); //###
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

  /*
  augmentedAlignment.printSlice(cout,0,augmentedAlignment.getLength(),
				'+',60,true);
  INTERNAL_ERROR;
  */
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
    //cout<<"i="<<i<<endl;
    PhylogenyBranch *branch=&branches[i];
    PhylogenyNode *parentNode=branch->getParent();
    PhylogenyNode *childNode=branch->getChild();
    if(os) (*os)<<"BranchHMM for branch "<<parentNode->getName()<<":"
		<<childNode->getName()<<endl;
    BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
						     branch->getLength(),os);
    branchAttributes[i]=BranchAttributes(hmm,branch);
    branch->getDecoration()=&branchAttributes[i];
    //cout<<"CCC "<<branch<<" -> "<<&branchAttributes[i]<<endl;
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
			 branch->getUpMap(),true);
    }
  }
}



/****************************************************************
                   AVES::reconstructHistory()
 ****************************************************************/
void AVES::reconstructHistory(Taxon &parentTaxon,Taxon &childTaxon,
			      IndexMap &downMap,IndexMap &upMap,
			      bool noParentSeqs)
{
  // First, resize the maps
  const GapPattern &parent=*parentTaxon.getGapPattern();
  const GapPattern &child=*childTaxon.getGapPattern();
  int parentL=parent.countResidues(), childL=child.countResidues();
  int L=parent.getLength();
  downMap.resize(parentL);
  upMap.resize(childL);
  if(noParentSeqs) {
    parentTaxon.setSeqLen(parentL);
    parentTaxon.getSeq()=Sequence(alphabet.lookup('A'),parentL);
  }
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
      else {
	downMap[parentI++]=IndexMap::UNDEFINED;
      }
    else {
      if(b==GPE_RESIDUE) upMap[childI++]=IndexMap::UNDEFINED;
    }
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
  cout<<"RUNNING DOLLO"<<endl;
  //cout<<"reconstructing alignment"<<endl;
  reconstructAlignment(true);
  //cout<<"inferring gap patterns"<<endl;
  inferGapPatterns(*fullAlignment);
  //cout<<"reconstructHistories"<<endl;
  reconstructHistories();
  if(usePrescan) {
    //cout<<"rebuilding prescan constraints"<<endl;
    rebuildPrescans(tree);
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
  //cout<<"decode() : parent="<<parent.getName()<<" child="<<child.getName()<<endl;
  //cout<<"parentLen="<<parent.getSeqLen()<<" childLen="<<child.getSeqLen()<<endl;
  ViterbiConstraint constraint=
    baselineMode ? VC_INDEL_AND_PARENT_PARSE : VC_PARENT_PARSE;
  StatePath *path=V.decode(constraint);
  //cout<<"delete V"<<endl;
  //delete &V;
  if(usePosteriors) {delete forward; delete backward;}

  //cout<<"updateFcConstraints"<<endl;
  updateFcConstraints(*path,child);
  
  //cout<<child.getName()<<" parse has length "<<child.getFunctionalParse().length()<<endl;
  cout<<"NEW path score: "<<V.scorePath(*path)<<" ("<<path->getScore()<<")"<<endl;
  //cout<<"update indel history"<<endl;
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
  Taxon &child=root.getIthChild(1); // 0);//###
  //cout<<root.getName()<<"->"<<child.getName()<<endl;
  ViterbiInterface &V=
    *(useHirschberg ?
      (ViterbiInterface*) new Hirschberg(&root,&child,F,bandwidth,
					 root.getSeqLen(),
					 child.getSeqLen(),usePrescan,
					 numThreads) :
      (ViterbiInterface*) new LinkViterbi(&root,&child,F,bandingType,
					  bandwidth));
  StatePath *path=V.decode(VC_INDEL_HISTORY); // posteriors not necessary!
  //delete &V;
  root.getFunctionalParse()=FunctionalParse(*path,PARENT);
  child.getBranchToParent()->setStatePath(path);
  //cout<<*path<<"\n"<<root.getFunctionalParse()<<endl;
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
  //delete &V;
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
  //delete &viterbi;
  if(usePosteriors) {delete forward; delete backward;}
  parent.getFunctionalParse()=FunctionalParse(*path,PARENT);
  child.getFunctionalParse()=FunctionalParse(*path,CHILD);
  updateIndels(branch,parent,child,*path);
  branch.setStatePath(path);

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
  //delete &V;
  if(fromIsParent)
    child->getFunctionalParse()=FunctionalParse(*path,CHILD);
  else
    parent->getFunctionalParse()=FunctionalParse(*path,PARENT);
  updateIndels(*branch,*parent,*child,*path);
  branch->setStatePath(path);
}



/*
int AVES::getPathIndex(StatePath &path,int childIndex,BranchHMM &hmm)
{
  int childI=-1, pathIndex=0; 
  for( ; childI<childIndex ; ++pathIndex) {
    if(emitsChild(hmm.getStateType(path[pathIndex])))
       ++childI;
  }
  return pathIndex-1;
}
*/



double AVES::computeScore(Taxon &taxon,FunctionalElement &elem)
{
  switch(scoringStrategy)
    {
    case SS_LLR: return computeScoreLLR(taxon,elem);
    case SS_POSTERIOR: return computeScorePosterior(taxon,elem);
    case SS_NONE: return 0.0;
    }
  INTERNAL_ERROR;
}



double AVES::computeScoreLLR(Taxon &taxon,FunctionalElement &elem)
{
  int begin=elem.getBegin();
  FunctionalClass bg=FunctionalClass::getBackground();
  const Array1D<double> &bgEqFreqs=bg.getEqFreqs();
  double score=0.0;
  BranchAttributes *branch=taxon.getBranchToParent();
  if(branch) {
    StatePath *path=branch->getStatePath();
    if(!path) INTERNAL_ERROR;
    BranchHMM *hmm=branch->getHMM();
    int pathIndex=path->childPosToIndex(begin);
    SubstitutionMatrix &bgPt=*hmm->getSubstMatrix(bg,bg);
    IndexMap &upmap=branch->getUpMap();
    Taxon &parent=*taxon.getParent();
    PrecomputedEmissions &childE=taxon.getPrecomputedEmissions(); 
    PrecomputedEmissions &parentE=parent.getPrecomputedEmissions();
    Sequence &seq=taxon.getSeq();
    FunctionalElementType T=elem.getType();
    Vector<FunctionalClass> &classes=T.getClasses();
    int L=classes.size();
    for(int i=0, pos=begin ; i<L ; ++i, ++pos) {
      STATE state=(*path)[pathIndex];
      FunctionalClass fcParent=hmm->getFunctionalClass(state,PARENT);
      FunctionalClass fcChild=hmm->getFunctionalClass(state,CHILD);
      const Array1D<double> &eqFreqs=hmm->getEqFreqs(state);
      PHMM_StateType stateType=hmm->getStateType(state);
      BOOM::Symbol s=seq[pos];
      SubstitutionMatrix &Pt=*hmm->getSubstMatrix(state);
      int parentPos=upmap[pos];
      double fgScore=
	computeScore(stateType,Pt,eqFreqs,parentPos,pos,parentE,childE,
		     fcParent,fcChild,parent,taxon);
      double bgScore=
	computeScore(stateType,bgPt,bgEqFreqs,parentPos,pos,parentE,childE,
		     bg,bg,parent,taxon);
      score+=fgScore-bgScore;
      ++pathIndex;
      while(!emitsChild(hmm->getStateType((*path)[pathIndex]))) ++pathIndex;
    }
  }
  return score;
}



double AVES::computeScore(PHMM_StateType stateType,SubstitutionMatrix &Pt,
			  const Array1D<double> &eqFreqs,int parentPos,
			  int childPos,PrecomputedEmissions &parentE,
			  PrecomputedEmissions &childE,
			  FunctionalClass fcParent,FunctionalClass fcChild,
			  Taxon &parent,Taxon &child)
{
  const int numAlpha=PureDnaAlphabet::global().size();
  Array1D<float> V(numAlpha), V2(numAlpha);
  switch(stateType)
    {
    case PHMM_MATCH:{
      Array3D<float>::IndexedTwice<float> childRow=
	childE[childPos][fcChild];
      Array3D<float>::IndexedTwice<float> parentRow=
	parentE[parentPos][fcParent];
      for(BOOM::Symbol sParent=0 ; sParent<numAlpha ; ++sParent) {
	double eq=getEQ(eqFreqs,sParent,&parent,parentPos,fcParent);
	for(BOOM::Symbol sChild=0 ; sChild<numAlpha ; ++sChild) {
	  V[sChild]=safeAdd(eq,safeAdd(parentRow[sParent],
				       Pt(sParent,sChild),childRow[sChild]));
	}	
	V2[sParent]=sumLogProbs<float>(V);
      }
      return sumLogProbs<float>(V2); }
    case PHMM_INSERT:{
      Array3D<float>::IndexedTwice<float> row=childE[childPos][fcChild];
      for(BOOM::Symbol s=0 ; s<numAlpha ; ++s) {
	double eq=getEQ(eqFreqs,s,&child,childPos,fcChild);
	V[s]=safeAdd(row[s],eq);
      }
      return sumLogProbs<float>(V); }
    }
  INTERNAL_ERROR;
}



double AVES::getEQ(const Array1D<double> &eqFreqs,BOOM::Symbol s,
		   Taxon *taxon,int pos,FunctionalClass fc) 
{
  if(contentSensor && fc.isBackground()) {
    const Sequence &seq=taxon->getSeq();
    const String &seqStr=taxon->getSeqStr();
    const char c=alphabet.lookup(s);
    return contentSensor->scoreSingleBase(seq,seqStr,pos,s,c);
  }
  else return log(eqFreqs[s]);
}



double AVES::computeScorePosterior(Taxon &taxon,FunctionalElement &elem)
{
  GainLossFelsenstein F(pureDnaAlphabet,identityMap,numTaxa);
  Taxon *parent=taxon.getParent();
  double logP;
  if(parent) {
    HirschPosteriors H(parent,&taxon,F,bandwidth,parent->getSeqLen(),
		       taxon.getSeqLen(),usePrescan);
    logP=H.compute(elem.getBegin(),elem.getEnd()-1,CHILD);
  }
  else {
    Taxon *child=&taxon.getIthChild(0);
    HirschPosteriors H(&taxon,child,F,bandwidth,taxon.getSeqLen(),
		       child->getSeqLen(),usePrescan);
    logP=H.compute(elem.getBegin(),elem.getEnd()-1,PARENT);
  }
  return exp(logP);
}



double AVES::computePostNumerator(Taxon &taxon,FunctionalElement &elem)
{
  GainLossFelsenstein F(pureDnaAlphabet,identityMap,numTaxa);
  Taxon *parent=taxon.getParent();
  double logP;
  if(parent) {
    HirschPosteriors H(parent,&taxon,F,bandwidth,parent->getSeqLen(),
		       taxon.getSeqLen(),usePrescan);
    logP=H.computeNumerator(elem.getBegin(),elem.getEnd()-1,CHILD);
  }
  else {
    Taxon *child=&taxon.getIthChild(0);
    HirschPosteriors H(&taxon,child,F,bandwidth,taxon.getSeqLen(),
		       child->getSeqLen(),usePrescan);
    logP=H.computeNumerator(elem.getBegin(),elem.getEnd()-1,PARENT);
  }
  return logP;
}



double AVES::computePostDenominator(Taxon &taxon)
{
  GainLossFelsenstein F(pureDnaAlphabet,identityMap,numTaxa);
  Taxon *parent=taxon.getParent();
  double logP;
  if(parent) {
    HirschPosteriors H(parent,&taxon,F,bandwidth,parent->getSeqLen(),
		       taxon.getSeqLen(),usePrescan);
    logP=H.computeLikelihood();
  }
  else {
    Taxon *child=&taxon.getIthChild(0);
    HirschPosteriors H(&taxon,child,F,bandwidth,taxon.getSeqLen(),
		       child->getSeqLen(),usePrescan);
    logP=H.computeLikelihood();
  }
  return logP;
}


