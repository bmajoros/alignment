/****************************************************************
 sample-annotations.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <unistd.h>
#include "sample-annotations.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "ConstrainedBackward.H"
#include "GainLossFelsenstein.H"

#define VERBOSE 0
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



Application::Application()
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
  CommandLine cmd(argc,argv,"b:Gr:c:pMti:s:g");
  if(cmd.numArgs()!=5)
    throw String("\n\
sample-annotations <model.lambda> <*.fasta> <*.aln> <#samples-per-alignment> <outfile>\n\
            -G = don't perform garbage collection\n\
            -r N = use randomization seed N\n\
            -p = print each alignment after sampling it\n\
            -s N = skip first N samples (default=0)\n\
            -g = emit GFF format instead of BART format\n\
");
  String lambdaFile=cmd.arg(0);
  String fastaFile=cmd.arg(1);
  alnFile=cmd.arg(2);
  numSamples=cmd.arg(3).asInt();
  outfile=cmd.arg(4);
  useGFF=cmd.option('g');
  //if(samplerType==GIBBS) throw "option -g is currently disabled";
  if(cmd.option('r')) seed=cmd.optParm('r').asUnsigned();
  else seed=GetRandomSeed();
  cout<<"seed="<<seed<<endl;
  SeedRandomizer(seed);
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
  //cout<<"executing "<<lambdaFile<<"..."<<endl;
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

  // Get phylogeny
  //cout<<"loading tree"<<endl;
  tree=modelCompiler->getGuideTree(transducerTemplate);
  tree->gatherNodes(phylogenyNodes);
  tree->constructBranches();
  if(cmd.option('l'))
    tree->scaleBranches(cmd.optParm('l').asFloat());
  tree->computePairwiseDistances();

  // Link taxa to phylogeny nodes
  //cout<<"linking taxa to phylogeny nodes"<<endl;
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
  //cout<<"instantiating HMM's"<<endl;
  initBranches();

  // Load sequences
  //cout<<"loading sequences"<<endl;
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

  for(int i=0;i<numTaxa;++i) taxa[i].getFunctionalParse().clear();
  cout.flush();
}



int Application::master(int argc,char *argv[])
{
  timer.startCounting();
  processCmdLine(argc,argv);

  int numSlaves=slaveDriver->getNumSlaves();
  int samplesPerSlave=numSamples/numSlaves;
  int leftovers=numSamples-samplesPerSlave*numSlaves;
  int nextID=1;
  for(int i=0 ; i<numSlaves ; ++i) {
    MpiVariableMsg &msg=*new MpiVariableMsg(SAMPLE);
    int num=samplesPerSlave+(i==numSlaves-1 ? leftovers : 0);
    msg<<num<<nextID;
    msg.close();
    slaveDriver->addWork(&msg);
    //cout<<"XXX "<<nextID<<" "<<num<<endl;
    nextID+=num;
  }

  Vector<MpiFixedMsg*> results;
  slaveDriver->waitForResults(results);
  Vector<MpiFixedMsg*>::iterator cur=results.begin(), end=results.end();
  for(; cur!=end ; ++cur) delete *cur;
  slaveDriver->terminateSlaves();

  // Consolidate files from slaves
  //cout<<"master combining files"<<endl;
  MPI &mpi=slaveDriver->getMPI();
  File out(outfile,"w");
  int n=numSamples*numSlaves;
  for(int i=0 ; i<numSlaves ; ++i) {
    int slaveID=mpi.getIthSlaveID(i);
    String slaveFile=outfile+"-"+slaveID;
    File in(slaveFile);
    long size=in.getSize();
    char *buffer=new char[size];
    in.read(size,(void*)buffer);
    out.write(size,(void*)buffer);
    delete [] buffer;
    unlink(slaveFile.c_str());
  }
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
  msg>>numSamples>>sampleID;

  samplingLoop();

  MpiVariableMsg *newMsg=new MpiFixedMsg(DONE_SAMPLING);
  slaveDriver->replyToMaster(newMsg);
}



/****************************************************************
                     Application::samplingLoop()
 ****************************************************************/
void Application::samplingLoop()
{
  int slaveID=slaveDriver->getMPI().getProcessID();
  outfile=outfile+"-"+slaveID;
  ofstream f(outfile.c_str());
  ResidueOrthologyGraph rog(tree,taxa);
  Vector<PhylogenyBranch*> branches;//=rog.getBranches();
  tree->collectBranches(branches,PREORDER);
  Vector<int> taxonIndices;
  for(int i=0 ; i<numTaxa ; ++i) taxonIndices.push_back(i);
  File aln(alnFile.c_str(),"r");
  while(rog.loadConnections(aln)) {
    if(aln.eof()) break;
    //dollo();
    for(int i=0 ; i<numSamples ; ++i) {
      //cout<<"sample #"<<i<<" slaveID="<<slaveID<<endl;
      if(!sample(f,branches,rog))
	cout<<"slave "<<slaveID<<" rejected sample "<<i<<endl;
      //else cout<<"ACCEPT"<<endl;
      if(wantPrint) {
	//dollo();
	fullAlignment->capByAnno();
	fullAlignment->printSlice(cout,0,fullAlignment->getLength(),
				  '+',60,true);
      }
    }
  }
  //cout<<"NO MORE ALIGNMENTS IN FILE"<<endl;
}



/****************************************************************
                     Application::sample()
 ****************************************************************/
bool Application::sample(ostream &f,Vector<PhylogenyBranch*> &branches,
			 ResidueOrthologyGraph &rog)
{
  f<<"sample "<<sampleID<<endl;
  int numBranches=branches.size();
  Vector<StatePath*> paths;
  for(int branchIndex=0 ; branchIndex<numBranches ; ++branchIndex) {
    //CollapsedOrthologyMatrix::compute(rog);
    PhylogenyBranch *branch=branches[branchIndex];
    /*
    cout<<"sampling branch "<<branch->getParent()->getName()<<":"
	<<branch->getChild()->getName()<<endl;
    */
    if(!sampleBranch(branch,f,paths)) return false;
    /*
    if(wantPrint) {
      reconstructAlignment(true); //false);
      fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,
				true);
    }
    */
  }

  for(int branchIndex=0 ; branchIndex<numBranches ; ++branchIndex) {
    PhylogenyBranch *br=branches[branchIndex];
    BranchAttributes *branch=
      static_cast<BranchAttributes*>(br->getDecoration());
    savePath(*paths[branchIndex],f,branch->getParentTaxon(),
	     branch->getChildTaxon(),*hmm);
  }
  f<<endl;
  ++sampleID;

  return true;
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
  hmm->getStatesOfType(PHMM_INSERT,QI);
  hmm->getStatesOfType(PHMM_DELETE,QD);
  hmm->getStatesOfType(PHMM_MATCH,QM);
  int numStates=hmm->getNumStates();
  stateTypes.resize(numStates); dI.resize(numStates); dJ.resize(numStates);
  for(int i=0 ; i<numStates ; ++i) {
    STATE t=stateTypes[i]=hmm->getStateType(i);
    switch(t) {
    case PHMM_MATCH:  dI[i]=1; dJ[i]=1; break;
    case PHMM_INSERT: dI[i]=0; dJ[i]=1; break;
    case PHMM_DELETE: dI[i]=1; dJ[i]=0; break;
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
                     Application::sampleBranch()
 ****************************************************************/
bool Application::sampleBranch(PhylogenyBranch *br,ostream &f,
			       Vector<StatePath*> &paths)
{
  BranchAttributes *branch=
    static_cast<BranchAttributes*>(br->getDecoration());
  BranchHMM *hmm=branch->getHMM();
  Taxon &parent=branch->getParentTaxon();
  Taxon &child=branch->getChildTaxon();
  int parentL=parent.getSeqLen(), childL=child.getSeqLen();
  StatePath &path=*new StatePath(hmm);
  FunctionalParse parentParse=parent.getFunctionalParse();
  GainLossFelsenstein F(pureDnaAlphabet,identityMap,numTaxa);
  ConstrainedBackward B(*hmm,parent,child,alphabetMap,numTaxa,F,parentParse);
  if(!B.ok()) return false; // very rare
  int i=0, j=0; // 1-based
  STATE prevQ=0;
  Vector<PHMM_StateType> &alignmentPath=*B.getPath();
  int pathLen=alignmentPath.size();
  for(int c=0 ; c<pathLen ; ++c) {
    Vector<STATE> states;
    Vector<float> probs;
    PHMM_StateType t=alignmentPath[c];
    const Array1D<STATE> &succ=hmm->getSuccessors(prevQ);
    STATE *p=&succ[0];
    int numSucc=succ.size();
    for(int s=0 ; s<numSucc ; ++s, ++p) {
      STATE q=*p;
      double trans=hmm->getTransP(prevQ,q);
      float cellP=safeAdd(trans,B.getCachedEmitP(q,c));
      cellP=safeAdd(cellP,B(c+1,q),-B(c,prevQ));
      if(isFinite(cellP)) {
	probs.push_back(exp(cellP));
	states.push_back(q);
      }
    }
    if(probs.size()==0) {
      //INTERNAL_ERROR;
      cout<<"XXX probs.size()==0 in sampleBranch"<<endl;
      cout<<"c="<<c<<" L="<<pathLen<<" prevQ="<<prevQ<<endl;
      cout<<br->getParent()->getName()<<":"<<br->getChild()->getName()<<endl;
      cout<<parentParse<<endl;
      cout<<alignmentPath<<endl;
      cout<<B<<endl;
      return false;
    }	
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
  installFuncParse(path,f,parent,child,*hmm);
  paths.push_back(&path);
  //cout<<path<<endl;
  cout.flush(); 
  return true;
}



void Application::savePath(StatePath &path,ostream &f,Taxon &parent,
			   Taxon &child,BranchHMM &hmm)
{
  FunctionalParse &parentParse=parent.getFunctionalParse();
  FunctionalParse &childParse=child.getFunctionalParse();
  if(parentParse.size()==0) parentParse.resize(parent.getSeqLen());
  if(childParse.size()==0) childParse.resize(child.getSeqLen());
  int PID=slaveDriver->getMPI().getProcessID();
  const String &parentName=parent.getName(), &childName=child.getName();
  int L=path.length(), parentL=parent.getSeqLen(), childL=child.getSeqLen();
  int parentPos=0, childPos=0, parentBegin, childBegin;
  FunctionalClass prevParentC, prevChildC;
  for(int i=0 ; i<L ; ++i) {
    STATE q=path[i];
    PHMM_StateType stateType=hmm.getStateType(q);
    FunctionalClass parentC=hmm.getFunctionalClass(q,PARENT);
    FunctionalClass childC=hmm.getFunctionalClass(q,CHILD);
    Strand parentStrand=parentC.getStrand();
    Strand childStrand=childC.getStrand();
    if(emitsParent(stateType)) parentParse[parentPos]=parentC;
    if(emitsChild(stateType)) childParse[childPos]=childC;

    bool beginsHere=false;
    if(parentC!=prevParentC) {
      if(parentC.beginsElement()) {
	parentBegin=parentPos;
	childBegin=childPos;
	beginsHere=true;
      }
      if(parentC.endsElement() && parent.getBranchToParent()==NULL &&
	 child.whichChild()==LEFT)
	if(useGFF) 
	  f<<parentName<<"\tBART\t"<<parentC.getElementType().getName()
	   <<"\t"<<parentBegin<<"\t"<<parentPos+1<<"\t.\t"<<parentStrand
	   <<"\t.\tsample="<<sampleID<<"\n";
	else
	  f<<parentName<<"\t"<<parentC.getElementType().getName()
	   <<"\t"<<parentBegin<<"\t"<<parentPos+1<<"\t"<<parentStrand
	   <<"\t"<<endl;
    }
    if(childC!=prevChildC) {
      if(childC.beginsElement()) {
	parentBegin=parentPos;
	childBegin=childPos;
	beginsHere=true;
      }
      if(childC.endsElement())
	if(useGFF)
	  f<<childName<<"\tBART\t"<<childC.getElementType().getName()<<"\t"
	   <<childBegin<<"\t"<<childPos+1<<"\t.\t"<<childStrand
	   <<"\t.\tsample="<<sampleID<<"\n";
	else 
	  f<<childName<<"\t"<<childC.getElementType().getName()<<"\t"
	   <<childBegin<<"\t"<<childPos+1<<"\t"<<childStrand<<"\t"<<endl;
    }
    if(beginsHere) {
      GainLossType gainLoss=FunctionalClass::classifyGainLoss(parentC,childC);
      String factor=parentC.isBackground() ? 
	childC.getElementType().getName() : parentC.getElementType().getName();
      Strand strand=parentC.isBackground() ? childStrand : parentStrand;
      if(useGFF)
	f<<parentName<<":"<<childName<<"\tBART\t"<<toString(gainLoss)
	 <<"\t"<<parentBegin<<"\t"<<childBegin<<"\t.\t"<<strand<<
	  "\t.\tsample="<<sampleID<<";factor="<<factor<<"\n";
      else
	f<<parentName<<":"<<childName<<"\t"<<toString(gainLoss)
	 <<"\t"<<parentBegin<<"\t"<<childBegin<<"\t"<<strand
	 <<"\t"<<factor<<endl;
    }
    hmm.updateColumnsFwd(q,parentPos,childPos);
    prevParentC=parentC;
    prevChildC=childC;
  }
}




void Application::installFuncParse(StatePath &path,ostream &f,Taxon &parent,
			   Taxon &child,BranchHMM &hmm)
{
  FunctionalParse &parentParse=parent.getFunctionalParse();
  FunctionalParse &childParse=child.getFunctionalParse();
  if(parentParse.size()==0) parentParse.resize(parent.getSeqLen());
  if(childParse.size()==0) childParse.resize(child.getSeqLen());
  int L=path.length();
  int parentPos=0, childPos=0, parentBegin, childBegin;
  FunctionalClass prevParentC, prevChildC;
  for(int i=0 ; i<L ; ++i) {
    STATE q=path[i];
    PHMM_StateType stateType=hmm.getStateType(q);
    FunctionalClass parentC=hmm.getFunctionalClass(q,PARENT);
    FunctionalClass childC=hmm.getFunctionalClass(q,CHILD);
    if(emitsParent(stateType)) parentParse[parentPos]=parentC;
    if(emitsChild(stateType)) childParse[childPos]=childC;
    hmm.updateColumnsFwd(q,parentPos,childPos);
    prevParentC=parentC;
    prevChildC=childC;
  }
}


