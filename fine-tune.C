/****************************************************************
 fine-tune.C : re-estimates parameters from unannotated sequences
 bmajoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "BOOM/MPI.H"
#include "BOOM/Vector.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Regex.H"
#include "BOOM/GSL/Optimizer.H"
#include "BOOM/Map.H"
#include "BandingPattern.H"
#include "HirschBackwardSumPair.H"
#include "AVES.H"
#include "SiblingHirschberg.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;
#include "FineTuner.H"


const double STEP_SIZE=0.1;//0.1;  ### BEST VALUE IS 0.1 ###
const double EPSILON=0.1;//2.2e-5; ### BEST VALUE IS 0.1 ###
const double TOLERANCE=0.001;//10;//0.1;// 0.01; // ### NO EFFECT ###
const double GRADIENT_THRESHOLD=1.0;//0.1;
const int MAX_ITERATIONS=50;
const int MAX_PAIRS=10;

/****************************************************************
                           main()
 ****************************************************************/
int main(int argc,char *argv[])
  {
    try
      {
	FineTuner app;
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



/****************************************************************
                          FineTuner::FineTuner()
 ****************************************************************/
FineTuner::FineTuner()
{
  pairType=
    //PT_UNIFORM_LENGTHS;
    PT_ONE_PAIR_PER_CLADE;
}



/****************************************************************
                          FineTuner::main()
 ****************************************************************/
int FineTuner::main(int argc,char *argv[])
{
  MPI *mpi=new MPI(argc,&argv);
  slaveDriver=new MpiSlaveDriver(*mpi);
  if(slaveDriver->amItheMaster()) return master(argc,argv);
  else return slave(argc,argv);
}



int FineTuner::master(int argc,char *argv[])
{
  cout<<"master"<<endl;

  // Process command line
  CommandLine cmd(argc,argv,"b:vp:t:e:GT:");
  if(cmd.numArgs()!=4)
    throw String(
"fine-tune [options] <in.fasta> <model.lambda> <target-track> <parms.lambda>\n\
     where: -b <strategy> = band using strategy:\n\
                f### = fixed-width banding of width 2*###\n\
            -v = use standard Viterbi, rather than Hirschberg\n\
            -p <in.gff> = use pre-scanned sites to limit search space\n\
            -t N = use N threads\n\
            -e <*.model> = use given content sensor for EQ frequencies\n\
            -G = don't perform garbage collection\n\
            -T N = use triplets instead of pairs\n\
\n\
     WARNING: <parms.lambda> will be overwritten!\n\
\n\
");
  String fastaFile=cmd.arg(0);
  lambdaFile=cmd.arg(1);
  String targetName=cmd.arg(2);
  outfile=cmd.arg(3);
  shouldSeed=shouldHillClimb=baselineMode=usePosteriors=shouldSample=false;
  noPrescanOnDownpass=shouldDumpModel=explicitHistories=false;
  useTriplets=cmd.option('T');
  numSamples=0;
  randomize();
  useHirschberg=!cmd.option('v');
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
  disableGC=cmd.option('G');
  if(disableGC) modelCompiler->setGCthreshold(LARGEST_INTEGER);
  numAlpha=alphabet.size();

  // Load content sensor for higher-order equilibrium frequencies
  contentSensor=
    useContentSensor ?
    ContentSensor::load(cmd.optParm('e')) :
    NULL;
  
  // Load HMM
  cout<<"executing "<<lambdaFile<<"..."<<endl;
  LambdaAPI &lambda=modelCompiler->getLambda();
  modelCompiler->parse(lambdaFile);
  cout<<"deleting foreign objs"<<endl;
  modelCompiler->deleteForeignObjects();
  cout<<"checking garbage"<<endl;
  lambda.checkGarbageLevel();

  cout<<"getting ctor"<<endl;
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  cout<<"instantiating"<<endl;
  transducerTemplate=modelCompiler->instantiate(ctor);
  cout<<"getting num states"<<endl;
  cout<<transducerTemplate->getNumStates()<<" states"<<endl;
  collectModelParms();

  // Get phylogeny
  tree=modelCompiler->getGuideTree(transducerTemplate);
  tree->gatherNodes(phylogenyNodes);

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
    if(node->getName()==targetName) targetSpecies=&taxon;
  }
  if(!targetSpecies) throw "Target species not found";
  alignmentBuilder=
    new AlignmentBuilder(tree->getRoot(),alphabet,gapSymbol,numTaxa);
  gatherCladeMembers();

  // Instantiate BranchHMM's and substitution matrices on tree
  cout<<"instantiating HMM's"<<endl;
  initBranches();

  // Load sequences
  cout<<"loading sequences"<<endl;
  FastaReader fastaReader(fastaFile,alphabet);
  String def, seq, name, junk;
  while(fastaReader.nextSequence(def,seq)) {
    FastaReader::parseDefline(def,name,junk);
    if(!nameToTaxon.isDefined(name)) continue;
    Sequence S(seq,alphabet);
    BOOM::Symbol A=alphabet.lookup('A');
    if(S.getLength()==0) S.append(A);
    //S.replaceAll(INVALID_SYMBOL,gapSymbol);
    //S.replaceAll(alphabet.lookup('N'),gapSymbol);
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
  cout<<"constructing pairs/triples"<<endl;
  if(useTriplets) collectTriplets(cmd.optParm('T').asInt());
  else constructPairs();

  // Start the timer
  timer.startCounting();

  // Perform gradient ascent
  cout<<"optimizing"<<endl;
  fineTune();

  // Clean up
  timer.stopCounting();
  cout<<"Elapsed time: "<<timer.elapsedTime()<<endl;
  MemoryProfiler::report("Total memory used by master:");
  cout<<"done"<<endl;
  return 0;
}



/****************************************************************
                       mapping probabilities
 ****************************************************************/

const double pi=2.0*acos(0.0);

inline double mapToProbability(double x)
{
  return atan(x)/pi+0.5;
}

inline double mapFromProbability(double x) 
{
  return tan(pi*(x-0.5));
}

GSL::Vector mapToProbabilities(const GSL::Vector &v,int numBranches)
{
  GSL::Vector r=v;
  int dimensionality=r.getDim();
  for(int i=0 ; i<dimensionality ; ++i)
    r[i]=mapToProbability(r[i]);
  return r;
}

GSL::Vector mapFromProbabilities(const GSL::Vector &v,int numBranches)
{
  GSL::Vector r=v;
  int dimensionality=r.getDim();
  for(int i=0 ; i<dimensionality ; ++i)
    r[i]=mapFromProbability(r[i]);
  return r;
}



/****************************************************************
                       FineTuner methods
 ****************************************************************/

void FineTuner::initBranches()
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
    //BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
    //					     branch->getLength(),os);
    BranchHMM *hmm=NULL;
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



Taxon *FineTuner::findLeaf(Taxon *parent)
{
  while(!parent->isLeaf()) parent=&parent->getIthChild(0);
  return parent;
}



void FineTuner::uniformLengths()
{
  Vector<SpeciesPair> pairs;
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxonI=taxa[i];
    if(taxonI.getNodeType()!=LEAF_NODE) continue;
    PhylogenyNode *left=taxonI.getNode();
    for(int j=i+1 ; j<numTaxa ; ++j) {
      Taxon &taxonJ=taxa[j];
      if(taxonJ.getNodeType()!=LEAF_NODE) continue;
      PhylogenyNode *right=taxonJ.getNode();
      double distance=tree->distanceBetween(left,right);
      pairs.push_back(SpeciesPair(&taxonI,&taxonJ,NULL));
      pairs[pairs.size()-1].distance=distance;
    }
  }

  struct PairCmp : public Comparator<SpeciesPair> {
    bool equal(SpeciesPair &a,SpeciesPair &b) {return a.distance==b.distance;}
    bool greater(SpeciesPair &a,SpeciesPair &b){return a.distance>b.distance;}
    bool less(SpeciesPair &a,SpeciesPair &b)  {return a.distance<b.distance;}
  } cmp;
  VectorSorter<SpeciesPair> sorter(pairs,cmp);
  sorter.sortAscendInPlace();

  // MAX_PAIRS
  int n=pairs.size();
  float delta=n/float(MAX_PAIRS-1);
  if(MAX_PAIRS>=n || delta<=0) {
    this->pairs=pairs;
  }
  else {
    this->pairs.push_back(pairs[0]);
    this->pairs.push_back(pairs[n-1]);
    float f=delta;
    for(int i=0 ; i<MAX_PAIRS-2 ; ++i) {
      this->pairs.push_back(pairs[int(f)]);
      f+=delta;
    }
  }
  n=this->pairs.size();
  for(int i=0 ; i<n ; ++i) {
    SpeciesPair &p=this->pairs[i];
    p.hmm=
	templateInstantiator->instantiate(transducerTemplate,p.distance,NULL);
  }
  BranchHMM *hmm=this->pairs[0].hmm;
  hmm->getStatesOfType(PHMM_INSERT,QI);
  hmm->getStatesOfType(PHMM_DELETE,QD);
  hmm->getStatesOfType(PHMM_MATCH,QM);
  numI=QI.size(), numD=QD.size(), numM=QM.size();
}



void FineTuner::onePairPerClade()
{
  bool first=true;
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    if(taxon.getNodeType()==INTERNAL_NODE) {
      Taxon *left=findLeaf(&taxon.getIthChild(0));
      Taxon *right=findLeaf(&taxon.getIthChild(1));
      double branchLen=
	tree->distanceBetween(left->getNode(),right->getNode());
      BranchHMM *hmm=
	templateInstantiator->instantiate(transducerTemplate,branchLen,NULL);
      pairs.push_back(SpeciesPair(left,right,hmm));
      pairs[pairs.size()-1].distance=branchLen;
      //pairs.push_back(SpeciesPair(right,left,hmm));
      //pairs[pairs.size()-1].distance=branchLen;
      if(first) {
	hmm->getStatesOfType(PHMM_INSERT,QI);
	hmm->getStatesOfType(PHMM_DELETE,QD);
	hmm->getStatesOfType(PHMM_MATCH,QM);
	numI=QI.size(), numD=QD.size(), numM=QM.size();
	first=false;
      }
    }
  }
}



void FineTuner::constructPairs()
{
  switch(pairType) 
    {
    case PT_ONE_PAIR: throw "not implemented";
    case PT_ALL_PAIRS: throw "not implemented";
    case PT_ONE_PAIR_PER_CLADE:
      onePairPerClade();
      break;
    case PT_N_OVER_2: throw "not implemented";
    case PT_UNIFORM_LENGTHS:
      uniformLengths();
      break;
    }
}



double FineTuner::pairwiseLikelihood(int whichPair)
{
  double sum=0.0;
  Vector<SpeciesPair>::iterator cur=pairs.begin(), end=pairs.end();
  int pairIndex=0;
  for(; cur!=end ; ++cur, ++pairIndex) {
    if(whichPair>=0 && whichPair!=pairIndex) continue;
    SpeciesPair &p=*cur;
    Taxon *left=p.first, *right=p.second;
    //cout<<"left="<<left->getName()<<" right="<<right->getName()<<endl;
    BranchHMM *hmm=p.hmm;
    int leftLen=left->getSeqLen(), rightLen=right->getSeqLen();
    BandingPattern bandingPattern(leftLen,rightLen,bandwidth);//###
    HirschBackwardSumPair B(//0,leftLen,0,rightLen, //###
			0,rightLen,0,leftLen,
			bandingPattern,*hmm,QI,numI,QD,numD,QM,numM,
			VC_UNCONSTRAINED,left,right,
			NULL, // parentFP
			NULL, // childFP
			NULL, // upMap
			NULL, // downMap
			NULL, // Hirschberg *
			contentSensor);
    B.setFirstTaxon(left);
    B.setSecondTaxon(right);

    // Init Hirschberg frame
    HirschbergFrame &frame=B.getFrame();
    WindowColumn &bCol=frame.getThisCol();
    Array2D<float>::RowIn2DArray<float> bPillar=bCol[leftLen];//###rightLen];
    int numStates=frame.getNumStates();
    for(STATE q=0 ; q<numStates ; ++q) {
      bPillar[q]=hmm->getTransP(q,0);
    }

    B.runPrescan();
    double LL=frame.getNextCol()[0][0];
    sum+=LL;
  }
  return sum;
}



void FineTuner::fineTune()
{
  const double tolerance=TOLERANCE;
  const double gradientThreshold=GRADIENT_THRESHOLD;
  const double stepSize=STEP_SIZE;////0.001;
  const int maxIterations=MAX_ITERATIONS;
  dimensionality=parms.size();

  // Compute an initial likelihood (summed over the pairs of sequences)
  ObjectiveFunction f(this);
  GSL::Vector initialPoint, optimalPoint;
  initialPoint.resize(dimensionality);
  copyIn(initialPoint);
  optimalPoint=initialPoint;
  cout<<"INITIAL POINT: "<<initialPoint<<endl;
  double initialLL=distributedLikelihood(initialPoint);
  cout<<"initialLL="<<initialLL<<endl;
  float optimalValue=initialLL;

  // Perform gradient ascent on the likelihood => find optimal model parms
  optimizer=new GSL::Optimizer(GSL::BFGS,f,initialPoint,stepSize,
                               GSL::BY_EITHER,tolerance,gradientThreshold,
                               maxIterations);

  cout<<"running optimizer"<<endl;
  optimizer->run();
  cout<<"optimization complete."<<endl;
  slaveDriver->terminateSlaves();

  // Collect the optimal point
  //cout<<"Collect the optimal point"<<endl;
  GSL::Vector bestPoint;
  bestPoint=f.getBestPoint();
  double bestValue=f.getBestScore();
  if(!isFinite(optimalValue) || bestValue>optimalValue 
     || bestPoint.getDim()!=optimalPoint.getDim())
    {
      optimalValue=bestValue;
      optimalPoint=bestPoint;
    }
  delete optimizer;
  cout<<"writing output"<<endl;
  writeFile(optimalPoint);
  cout<<"BEST POINT for this iteration: "<<bestPoint<<endl;
  cout<<"logL increase: "
      <<f.getFirstScore()
      <<" to "
      <<bestValue
      <<" in "<<f.getNumEvaluations()<<" evals + "
      <<f.getNumCacheHits()<<" cache hits"
      <<endl;
}


/*
void FineTuner::parseParmsFile(const String &filename)
{
  Regex regex("define\\s+'(\\S+)\\s+([^)\\s]+)");
  File f(filename);
  while(!f.eof()) {
    String line=f.getline();
    if(regex.search(line)) {
      String name=regex[0], value=regex[1];
      parms.push_back(ParmPair(name,value.asFloat()));
    }
  }
  f.close();
}
*/


void FineTuner::copyIn(GSL::Vector &pt)
{
  //cout<<"copyIn"<<endl;
  int n=parms.size();
  for(int i=0 ; i<n ; ++i) {
    ModelParm &p=parms[i];
    pt[i]=p.unconstrain(p.getValue()); // ### DEBUGGING
  }
  //cout<<"/copyIn"<<endl;
}



void FineTuner::copyOut(const GSL::Vector &pt)
{
  //cout<<"copyOut"<<endl;
  LambdaAPI &lambda=modelCompiler->getLambda();
  RunTimeEnvironment &env=lambda.getEnvironment();

  int n=parms.size();
  //cout<<"n="<<n<<endl;
  for(int i=0 ; i<n ; ++i) {
    //cout<<"i="<<i<<endl;
    ModelParm &p=parms[i];
    //cout<<"constrain"<<endl;
    double value=p.constrain(pt[i]);
    //cout<<"value="<<value<<endl;
    //double a=p.constrain(pt[i]);
    //double b=p.unconstrain(a);
    LambdaFloat *f=lambda.makeFloat(value,true);
    //cout<<"f="<<f<<endl;
    env.defineGlobal(p.getName(),f);
    //cout<<"copyOut: "<<p.getName()<<"="<<value<<endl;
  }
  //cout<<"reparse"<<endl;
  reparse();
  //cout<<"/copyOut"<<endl;
  //cout<<"Evaluating point: "<<pt<<endl;
}



void FineTuner::writeFile(const GSL::Vector &pt)
{
  //cout<<"writeFile"<<endl;
  //LambdaAPI &lambda=modelCompiler->getLambda();
  //RunTimeEnvironment &env=lambda.getEnvironment();

  ofstream os(outfile.c_str());
  int n=parms.size();
  for(int i=0 ; i<n ; ++i) {
    ModelParm &p=parms[i];
    double value=p.constrain(pt[i]);
    os<<"(define '"<<p.getName()<<" "<<value<<")"<<endl;
  }
  //cout<<"/writeFile"<<endl;
}



void FineTuner::collectModelParms()
{
  LambdaAPI &lambda=modelCompiler->getLambda();
  parms=modelCompiler->getModelParms();
  Vector<ModelParm>::iterator cur=parms.begin(), end=parms.end();
  for(; cur!=end ; ++cur) {
    ModelParm &parm=*cur;
    LambdaObject *obj=lambda.lookupGlobal(parm.getName());
    if(!obj) throw String("Undefined parm: ")+parm.getName();
    LambdaFloat *f=dynamic_cast<LambdaFloat*>(obj);
    if(!f) throw String("Parm has wrong type: ")+parm.getName();
    parm.setValue(f->getValue());
  }
}



void FineTuner::reparse()
{
  //cout<<"reparse"<<endl;
  delete transducerTemplate;
  //modelCompiler->deleteClosures();
  LambdaAPI &lambda=modelCompiler->getLambda();

  //FunctionalClass::resetAll();
  //FunctionalElementType::resetAll();

  //modelCompiler->recompile();
  //cout<<"DEL FOREIGN OBJ"<<endl;
  modelCompiler->deleteForeignObjects();
  //delete modelCompiler->getGuideTree();

  //lambda.checkGarbageLevel();
  //int numModels=modelCompiler->getNumModels();
  //cout<<numModels<<" models are in memory"<<endl;
  //transducerTemplate=modelCompiler->getIthModel(numModels-1);
  //cout<<"getting ctor"<<endl;
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  //cout<<"instantiate ctor "<<ctor<<endl;
  transducerTemplate=modelCompiler->instantiate(ctor);

  //cout<<"instantiation loop"<<endl;
  int n=useTriplets ? triplets.size() : pairs.size();
  for(int i=0 ; i<n ; ++i) {
    if(useTriplets) {
      SpeciesTriplet *triplet=triplets[i];
      Vector<PhylogenyBranch*> branches;
      triplet->tree->collectBranches(branches);
      int numBranches=branches.size();
      for(int j=0 ; j<numBranches ; ++j) {
	PhylogenyBranch *branch=branches[j];
	BranchAttributes *attr=
	  static_cast<BranchAttributes*>(branch->getDecoration());
	BranchHMM *hmm=
	  templateInstantiator->instantiate(transducerTemplate,
					    branch->getLength(),NULL);
	attr->changeHMM(hmm);
      }

    }
    else {
      SpeciesPair &p=pairs[i];
      //cout<<"del old hmm"<<endl;
      delete p.hmm;
      //cout<<"instantiate"<<endl;
      p.hmm=
	templateInstantiator->instantiate(transducerTemplate,p.distance,NULL);
    }
  }
  //cout<<"/reparse"<<endl;
}



void FineTuner::computeGradient(const GSL::Vector &currentPt,
				GSL::Vector &gradient)
{
  //cout<<"computeGradient"<<endl;
  // First, send all the jobs to the slaves

  writeFile(currentPt);

  int numPairs=useTriplets ? triplets.size() : pairs.size();
  const double epsilon=EPSILON;
  GSL::Vector perturbed=currentPt;
  int n=gradient.getDim();
  for(int i=0 ; i<n ; ++i)
    {
      double &x=perturbed[i];
      const double x0=x;
      double dx=epsilon*x;
      double temp=x+dx; 
      dx=temp-x; // make sure it's an "exact value" in the machine
      double twoDX=2*dx;
      x-=dx;
      for(int j=0 ; j<numPairs ; ++j) {
	MpiVariableMsg &msg0=*new MpiVariableMsg(COMPUTE_LIKELIHOOD);
	msg0<<i<<j<<0<<perturbed;
	msg0.close();
	slaveDriver->addWork(&msg0);
      }
      x+=twoDX;
      for(int j=0 ; j<numPairs ; ++j) {
	MpiVariableMsg &msg1=*new MpiVariableMsg(COMPUTE_LIKELIHOOD);
	msg1<<i<<j<<1<<perturbed;
	msg1.close();
	slaveDriver->addWork(&msg1);
      }	
      x=x0;
    }

  // Wait for the slaves to finish
  Vector<MpiFixedMsg*> results;
  slaveDriver->waitForResults(results);

  // Assemble results into a gradient vector
  gradient.resize(dimensionality);
  gradient.setAllTo(0.0);
  int i, positive;
  double LL;
  Vector<MpiFixedMsg*>::iterator cur=results.begin(), end=results.end();
  for(; cur!=end ; ++cur) {
    MpiFixedMsg &msg=**cur;
    msg>>i>>positive>>LL;
    double x=currentPt[i];
    double dx=epsilon*x;
    double temp=x+dx; 
    dx=temp-x; // make sure it's an "exact value" in the machine
    double twoDX=2*dx;
    if(!positive) LL=-LL;
    gradient[i]+=LL/twoDX;
    delete &msg;
  }
  cout<<"gradient = "<<gradient<<endl;
  //cout<<"/computeGradient"<<endl;
}



double FineTuner::distributedLikelihood(const GSL::Vector &currentPoint)
{
  int n=useTriplets ? triplets.size() : pairs.size();
  for(int i=0 ; i<n ; ++i) {
    MpiVariableMsg &msg=*new MpiVariableMsg(PAIR_LL);
    msg<<i<<currentPoint;
    msg.close();
    slaveDriver->addWork(&msg);
  }
  Vector<MpiFixedMsg*> results;
  slaveDriver->waitForResults(results);
  double LL=0.0;
  Vector<MpiFixedMsg*>::iterator cur=results.begin(), end=results.end();
  for(; cur!=end ; ++cur) {
    MpiFixedMsg &msg=**cur;
    double ll;
    msg>>ll;
    LL+=ll;
    delete &msg;
  }
  return LL;
}




/****************************************************************
                   ObjectiveFunction methods
 ****************************************************************/

ObjectiveFunction::ObjectiveFunction(FineTuner *app)
  : app(app), parms(app->getParms()), bestScore(NEGATIVE_INFINITY),
    numEvaluations(0), cacheHits(0), firstScore(NEGATIVE_INFINITY)
{
}



void ObjectiveFunction::gradient(const GSL::Vector &currentPt,
				   GSL::Vector &gradient)
{
  //cout<<"OF::gradient"<<endl;
  // PRECONDITION: currentPt has *not* been constrained

  static int iterations=0;
  twoPointSymmetric(currentPt,gradient);
  cout<<"#"<<(++iterations)<<" logL="<<likelihood<<"\tGRAD="
      <<gradient.norm()<<endl;
  //cout<<"/OF::gradient"<<endl;
}



void ObjectiveFunction::twoPointSymmetric(const GSL::Vector &currentPt,
					  GSL::Vector &gradient)
{
  app->computeGradient(currentPt,gradient);
}



void ObjectiveFunction::copyIn(GSL::Vector &pt)
{
  app->copyIn(pt);
}



void ObjectiveFunction::copyOut(const GSL::Vector &pt)
{
  app->copyOut(pt);
}




double ObjectiveFunction::f(const GSL::Vector &currentPoint)
{
  if(cache.isDefined(currentPoint)) {
    double value=cache[currentPoint];
    ++cacheHits;
    return value;
  }

  ++numEvaluations;

  //copyOut(currentPoint);
  //likelihood=app->pairwiseLikelihood();
  likelihood=app->distributedLikelihood(currentPoint);

  cout<<"newest likelihood: "<<currentPoint<<" LL="<<likelihood<<endl;
  if(!isFinite(bestScore) || likelihood>bestScore) {
    bestPoint=currentPoint;
    bestScore=likelihood;
  }
  cache[currentPoint]=-likelihood;
  if(!isFinite(firstScore)) firstScore=likelihood;
  return -likelihood;
}



bool ObjectiveFunction::pointIsValid(const GSL::Vector &point)
{
  INTERNAL_ERROR;
}



/****************************************************************
                         slave methods
 ****************************************************************/

int FineTuner::slave(int argc,char *argv[])
{
  //cout<<"slave"<<endl;
  // Process command line
  CommandLine cmd(argc,argv,"b:vp:t:e:GT:");
  if(cmd.numArgs()!=4)
    throw String(
"fine-tune [options] <in.fasta> <model.lambda> <target-track> <parms.lambda>\n\
     where: -b <strategy> = band using strategy:\n\
                f### = fixed-width banding of width 2*###\n\
            -v = use standard Viterbi, rather than Hirschberg\n\
            -p <in.gff> = use pre-scanned sites to limit search space\n\
            -t N = use N threads\n\
            -e <*.model> = use given content sensor for EQ frequencies\n\
            -G = don't perform garbage collection\n\
\n\
     WARNING: <parms.lambda> will be overwritten!\n\
\n\
");
  String fastaFile=cmd.arg(0);
  lambdaFile=cmd.arg(1);
  String targetName=cmd.arg(2);
  outfile=cmd.arg(3);
  shouldSeed=shouldHillClimb=baselineMode=usePosteriors=shouldSample=false;
  noPrescanOnDownpass=shouldDumpModel=explicitHistories=false;
  numSamples=0;
  randomize();
  useTriplets=cmd.option('T');
  useHirschberg=!cmd.option('v');
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
  //cout<<(numThreads+1)<<" threads"<<endl;
  disableGC=cmd.option('G');
  if(disableGC) modelCompiler->setGCthreshold(LARGEST_INTEGER);
  numAlpha=alphabet.size();

  // Load content sensor for higher-order equilibrium frequencies
  contentSensor=
    useContentSensor ?
    ContentSensor::load(cmd.optParm('e')) :
    NULL;
  
  // Load HMM
  //cout<<"executing "<<lambdaFile<<"..."<<endl;
  LambdaAPI &lambda=modelCompiler->getLambda();
  modelCompiler->parse(lambdaFile);
  //cout<<"slave done executing"<<endl;
  modelCompiler->deleteForeignObjects();
  //cout<<"checking garbage level"<<endl;
  lambda.checkGarbageLevel();

  //cout<<modelCompiler->getNumModels()<<" models loaded"<<endl;
  //if(modelCompiler->getNumModels()<1)
  // throw "use (register-transducer ...) to register a transducer";
  //transducerTemplate=modelCompiler->getIthModel(0);
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  transducerTemplate=modelCompiler->instantiate(ctor);
  //cout<<transducerTemplate->getNumStates()<<" states"<<endl;
  collectModelParms();

  // Get phylogeny
  //cout<<"loading tree"<<endl;
  tree=modelCompiler->getGuideTree(transducerTemplate);
  tree->gatherNodes(phylogenyNodes);

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
    if(node->getName()==targetName) targetSpecies=&taxon;
  }
  if(!targetSpecies) throw "Target species not found";
  alignmentBuilder=
    new AlignmentBuilder(tree->getRoot(),alphabet,gapSymbol,numTaxa);
  //cout<<"gatherCladeMembers"<<endl;
  gatherCladeMembers();

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
    //S.replaceAll(INVALID_SYMBOL,gapSymbol);
    //S.replaceAll(alphabet.lookup('N'),gapSymbol);
    S.useOnlyTheseSymbols(isNucleotide);
    int ID=nameToTaxon[name];
    Taxon &taxon=taxa[ID];
    taxon.getSeq()=S;
    taxon.setSeqLen(S.getLength());
    taxon.getSeqStr()=S(alphabet);
  }
  
  // Load pre-scanned annotations to limit search space
  //cout<<"loading pre-scans"<<endl;
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    Array1D<FuncClassSet> &fcConstraints=taxon.getFcConstraints();
    fcConstraints.resize(taxon.getSeqLen());
  }
  if(usePrescan) loadFcConstraints(prescanFile);
  if(useTriplets) collectTriplets(cmd.optParm('T').asInt());
  else constructPairs();

  // Start the timer
  timer.startCounting();

  // Perform gradient ascent
  //cout<<"calling serveTheMaster()"<<endl;
  serveTheMaster();
  //cout<<"DONE SERVING THE MASTER"<<endl;

  // Clean up
  MemoryProfiler::report("Total memory used by slave:");
  /*
  timer.stopCounting();
  cout<<"Elapsed time: "<<timer.elapsedTime()<<endl;
  MemoryProfiler::report("Total memory used:");
  */
  //cout<<"done"<<endl;
  return 0;
}



void FineTuner::serveTheMaster()
{
  while(1) {
    MpiFixedMsg *msg=slaveDriver->acceptWork();
    switch(msg->getTag())
      {
      case COMPUTE_LIKELIHOOD:
	//cout<<"recv'd COMP_LIK"<<endl;
	slaveLikelihood(*msg);
	break;
      case TERMINATE_SLAVE:
	//cout<<"recv'd TERMINATE"<<endl;
	return;
      case PAIR_LL:
	//cout<<"recv'd PAIR_LL"<<endl;
	slavePairLL(*msg);
	break;
      }
    delete msg;
  }
}



void FineTuner::slaveLikelihood(MpiFixedMsg &msg)
{
  GSL::Vector pt;
  int dimension, firstOrSecond, whichPair;
  msg>>dimension>>whichPair>>firstOrSecond>>pt;
  copyOut(pt);
  double LL;
  try {
    LL=useTriplets ? 
      -tripletLikelihood(whichPair) : -pairwiseLikelihood(whichPair);
  }
  catch(...) { LL=NEGATIVE_INFINITY; }

  MpiVariableMsg *newMsg=new MpiVariableMsg(LIKELIHOOD);
  *newMsg<<dimension<<firstOrSecond<<LL;
  newMsg->close();
  slaveDriver->replyToMaster(newMsg);
}



void FineTuner::slavePairLL(MpiFixedMsg &msg)
{
  GSL::Vector pt;
  int whichPair;
  //cout<<"parse message"<<endl;
  msg>>whichPair>>pt;
  //cout<<"copy out"<<endl;
  copyOut(pt);
  double LL;
  //cout<<"compute"<<endl;
  try {
    LL=useTriplets ? 
      tripletLikelihood(whichPair) : pairwiseLikelihood(whichPair);
  }
  catch(...) { LL=NEGATIVE_INFINITY; }
  //cout<<"send result"<<endl;
  MpiVariableMsg *newMsg=new MpiVariableMsg(LIKELIHOOD);
  *newMsg<<LL;
  newMsg->close();
  slaveDriver->replyToMaster(newMsg);
}



int SpeciesTriplet::find(const String &name,Array1D<Taxon> &taxa)
{
  int n=taxa.size();
  for(int i=0 ; i<n ; ++i)
    if(taxa[i].getName()==name)
      return i;
  INTERNAL_ERROR;
}


void SpeciesTriplet::initTaxa(Array1D<Taxon> &T)
{
  Vector<PhylogenyNode*> nodes;
  tree->gatherNodes(nodes);
  int numNodes=nodes.size();
  taxa.resize(numNodes);
  for(int i=0 ; i<numNodes ; ++i) {
    PhylogenyNode *node=nodes[i];
    node->setID(i);
    int taxonIndex=find(node->getName(),T);
    taxa[i]=T[taxonIndex];
    //cout<<"XXX "<<taxa[i].getFcConstraints().size()<<" vs. "<<T[taxonIndex].getFcConstraints().size()<<" "<<T[taxonIndex].getName()<<endl;
    Taxon *taxonI=&taxa[i];
    delete node->getDecoration();
    node->getDecoration()=taxonI;
    taxonI->setNode(node);
    //taxonI->removeFcConstraint();
    taxonI->setCladeAlignment(NULL);
    taxonI->setGapPattern(NULL);
    taxonI->removeBranchToParent();
    switch(node->getNodeType())
      {
      case ROOT_NODE: {
	PhylogenyBranch *branch=&node->getBranch(0);
	delete branch->getDecoration();
	BranchAttributes *attr=new BranchAttributes(NULL,branch);
	branch->getDecoration()=attr;
	taxonI->getBranch(0)=attr;
        }
	break;
      case INTERNAL_NODE: {
	PhylogenyBranch *leftBranch=&node->getBranch(0);
	PhylogenyBranch *rightBranch=&node->getBranch(1);
	delete leftBranch->getDecoration();
	delete rightBranch->getDecoration();
	BranchAttributes *leftAttr=new BranchAttributes(NULL,leftBranch);
	BranchAttributes *rightAttr=new BranchAttributes(NULL,rightBranch);
	leftBranch->getDecoration()=leftAttr;
	rightBranch->getDecoration()=rightAttr;
	taxonI->getBranch(LEFT)=leftAttr;
	taxonI->getBranch(RIGHT)=rightAttr;
        }
	break;
      }
  }
}



void FineTuner::collectTriplets(int numTriplets)
{
  Vector<int> leafIDs;
  for(int i=0 ; i<numTaxa ; ++i)
    if(taxa[i].isLeaf()) leafIDs.push_back(i);
  int numLeaves=leafIDs.size();
  Vector<SpeciesTriplet*> triplets;
  for(int i=0 ; i<numLeaves-2 ; ++i)
    for(int j=i+1 ; j<numLeaves-1 ; ++j)
      for(int k=j+1 ; k<numLeaves ; ++k) {
	float d=tree->getSpannedDistance(i,j,k);
	triplets.push_back(new SpeciesTriplet(leafIDs[i],leafIDs[j],
					      leafIDs[k],d));
      }
  struct Cmp : public BOOM::Comparator<SpeciesTriplet*> {
    bool equal(SpeciesTriplet* &a,SpeciesTriplet* &b)
      { return a->distance==b->distance; }
    bool greater(SpeciesTriplet* &a,SpeciesTriplet* &b)
      { return a->distance>b->distance; }
    bool less(SpeciesTriplet* &a,SpeciesTriplet* &b)
      { return a->distance<b->distance; }
  } cmp;
  VectorSorter<SpeciesTriplet*> sorter(triplets,cmp);
  sorter.sortAscendInPlace();
  if(numTriplets>triplets.size()) numTriplets=triplets.size();//INTERNAL_ERROR;
  float inc=float(triplets.size())/numTriplets;
  int next=0;
  for(int i=0 ; i<numTriplets ; ++i) {
    SpeciesTriplet *triplet=triplets[next];
    this->triplets.push_back(triplet);
    Phylogeny *tripletTree=triplet->tree=tree->clone();
    Set<String> keepTaxa;
    for(int j=0 ; j<3 ; ++j)
      keepTaxa.insert(taxa[triplet->ID[j]].getName());
    tripletTree->prune(keepTaxa);
    Vector<PhylogenyNode*> nodes;
    tripletTree->gatherNodes(nodes);
    tripletTree->constructBranches();
    triplet->initTaxa(taxa);
    triplets[next]=NULL;
    Vector<PhylogenyBranch*> branches;
    tripletTree->collectBranches(branches);
    int numBranches=branches.size();
    for(int j=0 ; j<numBranches ; ++j) {
      PhylogenyBranch *branch=branches[j];
      BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
      					     branch->getLength(),NULL);
      static_cast<BranchAttributes*>(branch->getDecoration())->changeHMM(hmm);
    }
    next=int(next+inc);
    if(next>=triplets.size()) next=triplets.size()-1;
    //cout<<"TREE: "<<*tripletTree<<endl;
  }
  int n=triplets.size();
  for(int i=0 ; i<n ; ++i) delete triplets[i];
}



double FineTuner::tripletLikelihood(int whichTriplet)
{
  double LL=0.0;
 
  int numTriplets=triplets.size();
  for(int i=0 ; i<numTriplets ; ++i) {
    if(whichTriplet>=0 && i!=whichTriplet) continue;
    SpeciesTriplet *triplet=triplets[i];
    int numNodes=triplet->tree->getNumNodes();
    
    // Phase I : The UP-PASS
    LossyFelsenstein F1(pureDnaAlphabet,identityMap,numNodes);
    uppass(F1,triplet);
    dollo(triplet);
    
    // Phase II : The DOWN-PASS
    GainLossFelsenstein F2(pureDnaAlphabet,identityMap,numNodes);
    inferRootFuncParse(F2,triplet);
    downpass(F2,triplet);
    
    LL+=likelihood(triplet);
  }
  return LL;
}


//OK
void FineTuner::uppass(LinkFelsenstein &F,SpeciesTriplet *triplet)
{
  Vector<PhylogenyNode*> nodes;
  triplet->tree->gatherNodes(nodes,POSTORDER);
  int n=nodes.size();
  InternalNode *iNode;
  for(int i=0 ; i<n ; ++i) {
    PhylogenyNode *node=nodes[i];
    switch(node->getNodeType()) {
    case ROOT_NODE: 
      progressiveRoot(node,triplet);
      break;
    case LEAF_NODE: break;
    case INTERNAL_NODE:
      iNode=static_cast<InternalNode*>(node);
      align(iNode->getLeft(),iNode->getRight(),F,triplet);
      break;
    }
  }
}


//OK
void FineTuner::progressiveRoot(PhylogenyNode *node,SpeciesTriplet *triplet)
{
  RootNode *root=static_cast<RootNode*>(node);
  Taxon *rootTaxon=static_cast<Taxon*>(root->getDecoration());
  BranchAttributes *branch=
    static_cast<BranchAttributes*>(root->getBranch(0).getDecoration());
  Taxon *childTaxon=&branch->getChildTaxon();
  IndexMap &downMap=branch->getDownMap();
  IndexMap &upMap=branch->getUpMap();
  int childLen=childTaxon->getSeqLen();
  rootTaxon->setSeqLen(childLen);
  downMap.resize(childLen);
  upMap.resize(childLen);
  for(int i=0 ; i<childLen ; ++i) connect(downMap,i,upMap,i);
}



//OK
void FineTuner::align(PhylogenyNode *leftNode,PhylogenyNode *rightNode,
		      LinkFelsenstein &F,SpeciesTriplet *triplet)
{
  //cout<<"align()"<<endl;

  // Get the two taxa
  Taxon &left=static_cast<Taxon&>(*leftNode->getDecoration());
  Taxon &right=static_cast<Taxon&>(*rightNode->getDecoration());
  PhylogenyNode *parentNode=leftNode->getParent();
  Taxon &parent=static_cast<Taxon&>(*parentNode->getDecoration());
  //cout<<"left="<<left.getName()<<" right="<<right.getName()
  //  <<" parent="<<parent.getName()<<endl;

  // Instantiate a new PairHMM for the combined branch lengths connecting
  // the two sibling taxa to be aligned
  BranchAttributes *leftBranch=
    static_cast<BranchAttributes*>(parentNode->getBranch(LEFT).
				   getDecoration());
  BranchAttributes *rightBranch=
    static_cast<BranchAttributes*>(parentNode->getBranch(RIGHT).
				   getDecoration());
  //cout<<*triplet->tree<<endl;
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
  StatePath *path=NULL;
  try { path=viterbi.decode(); }
  catch(...) {
    delete &viterbi;
    delete path;
    delete hmm;
    throw;
  }
  delete &viterbi;
  left.getFunctionalParse()=FunctionalParse(*path,PARENT);
  right.getFunctionalParse()=FunctionalParse(*path,CHILD);

  installSiblingHistories(*path,parent,leftBranch,rightBranch,hmm);
  delete path;
  delete hmm;

  if(usePrescan) updatePrescans(parent);

  if(useContentSensor) {
    parent.getSeq().resize(parent.getSeqLen());
    ResidueAddress ra(&parent,0);
    LinkParsimony parsimony(ra,pureDnaAlphabet,gapSymbol,numTaxa);
    parsimony.runFullSeq();
    parent.getSeqStr()=parent.getSeq()(alphabet);
  }
}


//OK
void FineTuner::dollo(SpeciesTriplet *triplet)
{
  Phylogeny *tree=triplet->tree;
  int numNodes=tree->getNumNodes();
  Vector<PhylogenyNode*> nodes;
  tree->gatherNodes(nodes);
  for(int i=0 ; i<numNodes ; ++i) {
    PhylogenyNode *node=nodes[i];
    Taxon &taxon=static_cast<Taxon&>(*node->getDecoration());
    taxon.getFunctionalParse().clear();
  }

  //reconstructAlignment();
  AlignmentBuilder alignmentBuilder(tree->getRoot(),alphabet,gapSymbol,
				    numNodes);
  delete fullAlignment;
  fullAlignment=alignmentBuilder.buildAlignment(true,true);
  inferGapPatterns(*fullAlignment,triplet);
  reconstructHistories(triplet);
  if(usePrescan) rebuildPrescans(tree);
}



// OK
void FineTuner::reconstructHistories(SpeciesTriplet *triplet)
{
  Phylogeny *tree=triplet->tree;
  Vector<PhylogenyBranch*> branches;
  tree->collectBranches(branches);
  int n=branches.size();
  for(int i=0 ; i<n ; ++i) {
    BranchAttributes *branch=
      static_cast<BranchAttributes*>(branches[i]->getDecoration());
    Taxon &child=branch->getChildTaxon();
    Taxon &parent=branch->getParentTaxon();
    reconstructHistory(parent,child,branch->getDownMap(),
		       branch->getUpMap());
  }
}


//OK
void FineTuner::inferRootFuncParse(LinkFelsenstein &F,SpeciesTriplet *triplet)
{
  Phylogeny *tree=triplet->tree;
  PhylogenyNode *rootNode=tree->getRoot();
  PhylogenyNode *childNode=rootNode->getChild(0);
  Taxon &root=*static_cast<Taxon*>(rootNode->getDecoration());
  Taxon &child=*static_cast<Taxon*>(childNode->getDecoration());
  //cout<<root.getName()<<"->"<<child.getName()<<endl;
  ViterbiInterface &V=
    *(useHirschberg ?
      (ViterbiInterface*) new Hirschberg(&root,&child,F,bandwidth,
					 root.getSeqLen(),
					 child.getSeqLen(),usePrescan,
					 numThreads) :
      (ViterbiInterface*) new LinkViterbi(&root,&child,F,bandingType,
					  bandwidth));
  StatePath *path=NULL;
  try { path=V.decode(VC_INDEL_HISTORY); }
  catch(...) {
    delete &V;
    throw;
  }
  delete &V;
  root.getFunctionalParse()=FunctionalParse(*path,PARENT);

  BranchAttributes *branch=
    static_cast<BranchAttributes*>(rootNode->
       getBranch(childNode->whichChild()).getDecoration());
  branch->setStatePath(path);
}


//OK
void FineTuner::downpass(LinkFelsenstein &F,SpeciesTriplet *triplet)
{
  //cout<<"downpass.."<<endl;
  Vector<PhylogenyNode*> nodes;
  Phylogeny *tree=triplet->tree;
  tree->gatherNodes(nodes,PREORDER);
  int n=nodes.size();
  for(int i=0 ; i<n ; ++i) {
    PhylogenyNode *node=nodes[i];
    Taxon &taxon=*static_cast<Taxon*>(node->getDecoration());
    switch(node->getNodeType()) {
    case ROOT_NODE: {
      RootNode *root=static_cast<RootNode*>(node);
      Taxon &child=static_cast<Taxon&>(*root->getChild()->getDecoration());
      child.getFunctionalParse()=taxon.getFunctionalParse();
    }
      break;
    case LEAF_NODE:
      // nothing to do
      break;
    case INTERNAL_NODE:
      downpass(taxon,F,taxon.getBranch(LEFT),triplet);
      downpass(taxon,F,taxon.getBranch(RIGHT),triplet);
      break;
    }
  }
}



void FineTuner::downpass(Taxon &parent,LinkFelsenstein &F,
			 BranchAttributes *branch,SpeciesTriplet *triplet)
{
  Phylogeny *tree=triplet->tree;
  int numNodes=tree->getNumNodes();
  Taxon &child=branch->getChildTaxon();
  ViterbiInterface &V=
    *(useHirschberg ?
      (ViterbiInterface*) new Hirschberg(&parent,&child,F,bandwidth,
					 parent.getSeqLen(),
					 child.getSeqLen(),usePrescan,
					 numThreads) :
      (ViterbiInterface*) new LinkViterbi(&parent,&child,F,bandingType,
					  bandwidth));
  ViterbiConstraint constraint=
    baselineMode ? VC_INDEL_AND_PARENT_PARSE : VC_PARENT_PARSE;
  StatePath *path=NULL;
  try { path=V.decode(constraint); }
  catch(...) {
    delete &V;
    throw;
  }
  delete &V;
  //cout<<"assign functional parse of length "<<path->length()<<endl;

  updateFcConstraints(*path,child);
  
  //cout<<child.getName()<<" parse has length "<<child.getFunctionalParse().length()<<endl;
  //cout<<"NEW path score: "<<V.scorePath(*path)<<" ("<<path->getScore()<<")"<<endl;
  //cout<<"update indel history"<<endl;
  branch->setStatePath(path); // also updates the indel history

  //reconstructAlignment(true);

  /* ###
  AlignmentBuilder alignmentBuilder(tree->getRoot(),alphabet,gapSymbol,
				    numNodes);
  delete fullAlignment;
  fullAlignment=alignmentBuilder.buildAlignment(true,true);
  */

  /*
  cout<<"NEW ALIGNMENT FOR "<<parent.getName()<<"->"<<child.getName()
      <<":"<<endl;
  fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',
			    60,true,true);
  */
}



double FineTuner::likelihood(SpeciesTriplet *triplet)
{
  // First, compute the contribution of the emission terms by applying
  // Felsenstein's algorithm to all connected components in the graph
  
  double logP=0;
  LinkFelsenstein F(alphabet,alphabetMap,numTaxa);
  Vector<PhylogenyBranch*> branches;
  Phylogeny *tree=triplet->tree;
  tree->collectBranches(branches);
  int n=branches.size();
  for(int i=0 ; i<n ; ++i) {
    PhylogenyBranch *pBranch=branches[i];
    BranchAttributes *branch=
      static_cast<BranchAttributes*>(pBranch->getDecoration());
    Taxon &taxon=branch->getChildTaxon();
    IndexMap &upMap=branch->getUpMap();
    int L=taxon.getSeqLen();
    for(int j=0 ; j<L ; ++j) {
      if(upMap[j]==IndexMap::UNDEFINED) {// residue is a root
	logP+=F.logLikelihood(j,taxon,false);
	if(!finite(logP)) 
	  {cout<<taxon.getName()<<" "<<j<<endl; INTERNAL_ERROR;}
      }
    }
  }
  Taxon &taxon=static_cast<Taxon&>(*tree->getRoot()->getDecoration());
  int L=taxon.getSeqLen();
  for(int j=0 ; j<L ; ++j) {
    logP+=F.logLikelihood(j,taxon,false);
    //if(!finite(logP)) {cout<<taxon.getName()<<" "<<j<<endl; INTERNAL_ERROR;}
  }

  // Next, process each branch individually to assess the transition 
  // probabilities for the transducers
  logP+=sumTransLogProbs(triplet);

  return logP;
}



void FineTuner::inferGapPatterns(const MultSeqAlignment &A,
				 SpeciesTriplet *triplet)
{
  // First, construct an alignment on 3 symbols: *, -, and ?, representing
  // a nucleotide, gap, or unknown element
  const Alphabet &gapAlphabet=GapPatternAlphabet::global();
  BOOM::Symbol unknown=gapAlphabet.lookup('?');
  BOOM::Symbol gap=gapAlphabet.lookup('-');
  BOOM::Symbol residue=gapAlphabet.lookup('*');
  int L=A.getLength();
  MultSeqAlignment augmentedAlignment(gapAlphabet,gap);
  Phylogeny *tree=triplet->tree;
  Vector<PhylogenyNode*> nodes;
  tree->gatherNodes(nodes);
  int numNodes=nodes.size();
  for(int i=0 ; i<numNodes ; ++i) {
    PhylogenyNode *node=nodes[i];
    Taxon &taxon=static_cast<Taxon&>(*node->getDecoration());
    const String taxonName=taxon.getName();
    AlignmentSeq &newTrack=augmentedAlignment.findOrCreateTrack(taxonName);
    newTrack.extendToLength(L,unknown);
    switch(node->getNodeType())
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
  for(int i=0 ; i<numNodes ; ++i) {
    PhylogenyNode *node=nodes[i];
    Taxon &taxon=static_cast<Taxon&>(*node->getDecoration());
    GapPattern *gp=new GapPattern(augmentedAlignment.getIthTrack(i).getSeq());
    taxon.setGapPattern(gp);
  }
  //augmentedAlignment.printSlice(cout,0,augmentedAlignment.getLength(),
  //			'+',60,true);
}



double FineTuner::sumTransLogProbs(SpeciesTriplet *triplet)
{
  double logP=0;
  Vector<PhylogenyBranch*> branches;
  triplet->tree->collectBranches(branches);
  int n=branches.size();
  for(int i=0 ; i<n ; ++i) {
    BranchAttributes &branch=
      static_cast<BranchAttributes&>(*branches[i]->getDecoration());
    BranchHMM *hmm=branch.getHMM();
    StatePath *path=branch.getStatePath();
    if(!path) INTERNAL_ERROR;
    logP+=hmm->getTransProbs(*path);
  }
  return logP;
}
