/****************************************************************
 train-background.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaDashAlphabet.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/Map.H"
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Random.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Alphabet.H"
#include "BOOM/FastaWriter.H"
#include "BOOM/Time.H"
#include "BOOM/Exceptions.H"
#include "BOOM/MPI.H"
#include "BOOM/GSL/Optimizer.H"
#include "BOOM/Environment.H"
#include "PairHMM/PairHMM.H"
#include "PairHMM/Transducer.H"
#include "PhyLib/Phylogeny.H"
#include "PhyLib/SingleGainParsimony.H"
#include "Taxon.H"
#include "GapPattern.H"
#include "BranchAttributes.H"
#include "FunctionalParse.H"
#include "AlignmentBuilder.H"
#include "LinkBackward.H"
#include "LinkSampler.H"
#include "ModelCompiler.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;


/****************************************************************
                         enum MessageType
 ****************************************************************/
enum MessageTag {
  MSG_COMPUTE_LIKELIHOOD,
  LIKELIHOOD_IS
};



/****************************************************************
                         class Application
 ****************************************************************/
class Application
{
  MPI *mpi;
  MpiSlaveDriver *slaveDriver;
  ModelCompiler *modelCompiler; // parses lambda program -> TransducerTemplate
  TransducerTemplate *transducerTemplate; // transducer with free variable t
  TemplateInstantiator *templateInstantiator; // makes transducer for fixed t
  Phylogeny *tree;
  Array1D<Taxon> taxa; // indexed by node ID's from the phylogeny
  MultSeqAlignment *fullAlignment;
  FunctionalParse *functionalParse;
  Vector<PhylogenyNode*> phylogenyNodes;
  Vector<PhylogenyBranch> branches;
  Vector<BranchAttributes> branchAttributes;
  DnaDashAlphabet alphabet;
  DropGapMapping alphabetMap;
  PureDnaAlphabet pureDnaAlphabet;
  BOOM::Symbol gapSymbol;
  BitSet gapSymbols;
  bool leavesOnly, shouldDumpModel; // command-line options
  Map<String,int> nameToTaxon; // index into "taxa" array
  int numTaxa;
  AlignmentBuilder *alignmentBuilder;
  int master(int argc,char *argv[]);
  int slave(int argc,char *argv[]);
  void gatherCladeMembers();
  void initBranches();
  void optimize();
public:
  Application();
  int main(int argc,char *argv[]);
};



/****************************************************************
                     class ObjectiveFunction
 ****************************************************************/
class ObjectiveFunction : public GSL::ObjectiveFunction
{
public:
  ObjectiveFunction();
  virtual double f(const GSL::Vector &currentPoint);
  virtual void gradient(const GSL::Vector &currentPoint,
			GSL::Vector &gradient);
};



/****************************************************************
                            main()
 ****************************************************************/
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
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



/****************************************************************
                    Application::Application()
 ****************************************************************/
Application::Application()
  : alphabetMap(DnaDashAlphabet::global(),PureDnaAlphabet::global()),
    functionalParse(NULL), alignmentBuilder(NULL)
{
  gapSymbol=alphabet.lookup('-');
  modelCompiler=new ModelCompiler;
  templateInstantiator=modelCompiler->getInstantiator();
  gapSymbols.setSize(alphabet.size());
  gapSymbols.addMember((unsigned long)(int)gapSymbol);
}



/****************************************************************
                      Application::main()
 ****************************************************************/
int Application::main(int argc,char *argv[])
{
  // Initialize MPI
  mpi=new MPI(argc,argv);
  slaveDriver=new MpiSlaveDriver(*mpi);
  int exitCode;

  // Determine whether this process is a MASTER or a SLAVE
  if(mpi->getProcessID()==0) exitCode=master(argc,argv);
  else exitCode=slave(argc,argv);

  return exitCode;
}



/****************************************************************
                   Application::gatherCladeMembers()
 ****************************************************************/
void Application::gatherCladeMembers()
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
                     Application::initBranches()
 ****************************************************************/
void Application::initBranches()
{
  ostream *os=NULL;
  String dumpFile="train.dump";
  if(shouldDumpModel) {
    os=new ofstream(dumpFile.c_str());
    if(!os->good()) throw String("Error writing to file: ")+dumpFile;
  }
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
                             MASTER
 ****************************************************************/
int Application::master(int argc,char *argv[])
{
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw String("train-background <model.lambda> <training-alignment.maf> <outfile.lambda>");
  String lambdaFile=cmd.arg(0);
  String trainAlignFile=cmd.arg(1);
  String outfile=cmd.arg(2);
  Lambda::GarbageCollector::setThreshold(10); // ###



  cout<<"terminating slaves..."<<endl;
  slaveDriver->terminateSlaves();
  cout<<"MASTER exiting..."<<endl;
  return 0;
}



/****************************************************************
                              SLAVE
 ****************************************************************/
int Application::slave(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  String lambdaFile=cmd.arg(0);
  String trainAlignFile=cmd.arg(1);
  String outfile=cmd.arg(2);
  Lambda::GarbageCollector::setThreshold(10); // ###

  // Serve the master...
  bool done=false;
  while(!done) {
    cout<<mpi->getProcessID()<<" waiting for work"<<endl;
    MpiFixedMsg &work=*slaveDriver->acceptWork();
    cout<<mpi->getProcessID()<<" received work"<<endl;
    switch(work.getTag()) 
      {
      case MSG_COMPUTE_LIKELIHOOD:
	{

	  // ... perform alignment and evaluate likelihood ...

	  MpiVariableMsg &reply=*new MpiVariableMsg(MSG_LIKELIHOOD_IS);
	  reply<<logLikelihood;
	  reply<<endmsg;
	  cout<<"replying to master "<<mpi->getProcessID()<<endl;
	  slaveDriver->replyToMaster(&reply);
	  cout<<"done replying "<<mpi->getProcessID()<<endl;
	}
	break;
      case TERMINATE_SLAVE:
	cout<<"slave "<<mpi->getProcessID()<<" quitting"<<endl;
	done=true;
	break;
      default:
	cout<<"SLAVE: unknown message!"<<endl;
	break;
      }
    delete &work;
  }
  return 0;
}



/****************************************************************
                     Application::optimize()
 ****************************************************************/
void Application::optimize()
{
  ObjectiveFunction f();
  GSL::Vector initialPoint;
  initialPoint.resize(dimensionality);
  f.copyIn(initialPoint);
  cout<<"INITIAL POINT: "<<initialPoint<<endl;
  GSL::OptimizerType optimizerType=
      GSL::stringToOptimizerType(optimizerTypeStr);
  GSL::StoppingCriterion stoppingCriterion=
      GSL::stringToStoppingCriterion(stoppingCriterionStr);

  cout<<"init bounds"<<endl;
  initBounds();
  
  cout<<"mapping initial point"<<endl;
  initialPoint=mapFromProbabilities(initialPoint,numBranches);
  cout<<"allocating optimizer"<<endl;
  optimizer=new GSL::Optimizer(optimizerType,f,initialPoint,stepSize,
                               stoppingCriterion,tolerance,gradientThreshold,
                               maxIterations);
  cout<<"running optimizer"<<endl;
  optimizer->run();
  cout<<"optimization complete."<<endl;

  // Collect the optimal point
  GSL::Vector bestPoint;
  bestPoint=f.getBestPoint();
  double bestValue=f.getBestScore();
  if(!isFinite(optimalValue) || bestValue>optimalValue 
     || bestPoint.getDim()!=optimalPoint.getDim()) // ### 7/23/07
    {
      optimalValue=bestValue;
      optimalPoint=bestPoint;
    }
  delete optimizer;
  copyOut(optimalPoint);
  cout<<"BEST POINT for this iteration: "<<bestPoint<<endl;
  cout<<"logL increase: "
      <<f.getFirstScore()
      <<" to "
      <<bestValue
      <<" in "<<f.getNumEvaluations()<<" evals + "
      <<f.getNumCacheHits()<<" cache hits"
      <<endl;
}



/****************************************************************
                   ObjectiveFunction methods
****************************************************************/
ObjectiveFunction::ObjectiveFunction()
{
}


double ObjectiveFunction::f(const GSL::Vector &currentPoint)
{
  for(int i=0 ; i<100 ; ++i) {
    MpiVariableMsg &message=*new MpiVariableMsg(MSG_WORK);
    message<<i;
    message<<endmsg;
    slaveDriver->addWork(&message);
  }
  Vector<MpiFixedMsg*> results;
  slaveDriver->waitForResults(results);
  int n=results.size();
  for(int i=0 ; i<n ; ++i) {
    MpiFixedMsg *msg=results[i];
    int x, x2;
    (*msg)>>x>>x2;
    cout<<"RESULT: "<<x<<" "<<x2<<endl;
    delete msg;
  }

}



void ObjectiveFunction::gradient(const GSL::Vector &currentPoint,
				 GSL::Vector &gradient)
{
}



/*
int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw String("train-background <model.lambda> <training-alignment.maf> <outfile.lambda>");
  String lambdaFile=cmd.arg(0);
  String trainAlignFile=cmd.arg(1);
  String outfile=cmd.arg(2);
  Lambda::GarbageCollector::setThreshold(10); // ###

  // Load HMM
  modelCompiler->parse(lambdaFile);
  if(modelCompiler->getNumModels()<1)
    throw "use (register-transducer ...) to register a transducer";
  transducerTemplate=modelCompiler->getIthModel(0);

  // Load phylogeny
  tree=modelCompiler->getGuideTree();
  tree->gatherNodes(phylogenyNodes);
  tree->gatherBranches(branches);

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

  // Load training alignment
  ifstream is(trainAlignFile.c_str());
  if(!is.good()) throw String("Error opening file: ")+trainAlignFile;
  MultiAlignment *maf=MultiAlignment::nextAlignmentFromMAF(is);
  fullAlignment=new MultSeqAlignment(*maf,alphabet,gapSymbol);
  delete maf;
  is.close();
  int alignmentLength=fullAlignment->getLength();

  // Infer ancestral sequences via Dollo parsimony
  cout<<"inferring ancestral sequences"<<endl;
  SingleGainParsimony parsimony(*tree,alphabet,*fullAlignment,gapSymbols);
  cout<<"dollo.run()"<<endl;
  parsimony.run();

  // Estimate initial values for parameters
  double sumT=0, sumAlpha=0, sumBeta=0;
  int numBranches=branches.size();
  for(int i=0 ; i<numBranches ; ++i) {
    PhylogenyBranch *branch=&branches[i];
    Taxon *parent=static_cast<Taxon*>(branch->getParent()->getDecoration());
    Taxon *child=static_cast<Taxon*>(branch->getChild()->getDecoration());
    String parentName=parent->getName(), childName=child->getName();
    AlignmentSeq &parentTrack=fullAlignment->getTrackByName(parentName);
    AlignmentSeq &childTrack=fullAlignment->getTrackByName(childName);
    const Sequence &parentSeq=parentTrack.getSeq(), 
      &childSeq=childTrack.getSeq();
    Vector<PHMM_StateType> path;
    for(int pos=0 ; pos<alignmentLength ; ++pos) {
      BOOM::Symbol a=parentSeq[pos], b=childSeq[pos];
      PHMM_StateType stateType;
      if(a==gapSymbol)
	if(b==gapSymbol)
	  continue;
	else
	  stateType=PHMM_INSERT;
      else 
	if(b==gapSymbol)
	  stateType=PHMM_DELETE;
	else
	  stateType=PHMM_MATCH;
      path.push_back(stateType);
    }
    double t=branch->getLength();
    sumT+=t;
    int pathLen=path.size();
    int matchSelfTrans=0, indelSelfTrans=0, outMatch=0, outIndel=0;
    for(int i=1 ; i<pathLen ; ++i) {
      PHMM_StateType prev=path[i-1], cur=path[i];
      if(prev==PHMM_MATCH) ++outMatch;
      else ++outIndel;
      if(prev!=cur) continue;
      if(prev==PHMM_MATCH) ++matchSelfTrans;
      else ++indelSelfTrans;
    }
    cout<<parentName<<":"<<childName<<" match="<<matchSelfTrans<<" indel="<<indelSelfTrans<<endl;
    float x=outMatch>0 ? matchSelfTrans/float(outMatch) : 0;
    float y=outIndel>0 ? indelSelfTrans/float(outIndel) : 0;
    float alpha=-log(x)/t;
    float beta=-log(1-y)/t;
    sumAlpha+=alpha*t;
    sumBeta+=beta*t;
    cout<<"t="<<t<<" alpha="<<alpha<<" beta="<<beta<<" x="<<x<<" y="<<y<<endl;
  }

  double aveAlpha=sumAlpha/sumT, aveBeta=sumBeta/sumT;
  cout<<"ave(alpha)="<<aveAlpha<<" ave(beta)="<<aveBeta<<endl;
  
  return 0;
  }
*/


