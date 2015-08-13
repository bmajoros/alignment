/****************************************************************
 build-library.C : compute posterior library for REAPR
 bmajoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include "BOOM/MPI.H"
#include "BOOM/Vector.H"
#include "BOOM/Map.H"
#include "BOOM/Exceptions.H"
#include "BandingPattern.H"
#include "BandedForward.H"
#include "BandedBackward.H"
#include "PosteriorMatrix.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;
#include "LibraryBuilder.H"


/****************************************************************
                           main()
 ****************************************************************/
int main(int argc,char *argv[])
  {
    try
      {
	LibraryBuilder app;
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
                 LibraryBuilder::LibraryBuilder()
 ****************************************************************/
LibraryBuilder::LibraryBuilder()
{
}



/****************************************************************
                          LibraryBuilder::main()
 ****************************************************************/
int LibraryBuilder::main(int argc,char *argv[])
{
  MPI *mpi=new MPI(argc,&argv);
  slaveDriver=new MpiSlaveDriver(*mpi);
  if(slaveDriver->amItheMaster()) return master(argc,argv);
  else return slave(argc,argv);
}



void LibraryBuilder::processCmdLine(int argc,char *argv[]) {
  // Process command line
  CommandLine cmd(argc,argv,"b:Gt:T:e");
  if(cmd.numArgs()!=4) throw String(
"\n\
build-library [options] <in.fasta> <model.lambda> <output-dir> <factor>\n\
     where: -b <strategy> = band using strategy:\n\
                f### = fixed-width banding of width 2*###\n\
            -G = don't perform garbage collection\n\
            -t ### = use given sparseness threshold (default: 0.001)\n\
            -T X = only build tables that pair any taxon with X\n\
            -e = re-estimate equilibrium frequecies from input seqs\n\
\n\
");
  fastaFile=cmd.arg(0);
  lambdaFile=cmd.arg(1);
  outDir=cmd.arg(2);
  factor=cmd.arg(3);
  sparsenessThreshold=
    cmd.option('t') ? cmd.optParm('t').asFloat() : 0.001;
  if(cmd.option('b')) {
    String strategy=cmd.optParm('b');
    if(strategy[0]=='f' && strategy.length()>1) {
      bandwidth=strategy.substring(1,strategy.length()-1).asInt();
      bandingType=FIXED_WIDTH_BANDING;
    }
    else throw "invalid banding strategy specified";
  }
  else bandingType=NO_BANDING;
  reestimateEq=cmd.option('e');
  if(cmd.option('T')) targetTaxon=cmd.optParm('T');
  disableGC=cmd.option('G');
  if(disableGC) modelCompiler->setGCthreshold(LARGEST_INTEGER);
  numAlpha=alphabet.size();
  
  // Load HMM
  //cout<<"executing "<<lambdaFile<<"..."<<endl;
  LambdaAPI &lambda=modelCompiler->getLambda();
  lambda.getGC().setSilence(true);
  modelCompiler->parse(lambdaFile);
  modelCompiler->deleteForeignObjects();
  lambda.checkGarbageLevel();

  //cout<<"instantiating model"<<endl;
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  transducerTemplate=modelCompiler->instantiate(ctor);
  //cout<<transducerTemplate->getNumStates()<<" states"<<endl;

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
  }
  alignmentBuilder=NULL;
  gatherCladeMembers();
  if(!targetTaxon.isEmpty()) {
    if(!nameToTaxon.isDefined(targetTaxon))
      throw String("Can't find taxon: ")+targetTaxon;
    targetID=nameToTaxon[targetTaxon];
  }
  else targetID=-1;

  // Instantiate BranchHMM's and substitution matrices on tree
  initBranches();

  // Load sequences
  //cout<<"loading sequences"<<endl;
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

  if(reestimateEq) estimateBgEq();
}



int LibraryBuilder::master(int argc,char *argv[])
{
  processCmdLine(argc,argv);
  
  // Start the timer
  timer.startCounting();

  // Perform gradient ascent
  //cout<<"building library"<<endl;
  buildLibrary();

  // Clean up
  timer.stopCounting();
  cout<<"Elapsed time: "<<timer.elapsedTime()<<endl;
  MemoryProfiler::report("Total memory used by master:");
  cout<<"done"<<endl;
  return 0;
}



/****************************************************************
                       LibraryBuilder methods
 ****************************************************************/

void LibraryBuilder::buildLibrary()
{
  // Send messages to slaves
  for(int i=0 ; i<numTaxa ; ++i) {
    if(!taxa[i].isLeaf()) continue;
    //for(int j=i+1 ; j<numTaxa ; ++j) {
    for(int j=0 ; j<numTaxa ; ++j) {
      if(j==i) continue;
      if(targetID>-1 && i!=targetID && j!=targetID) continue;
      if(!taxa[j].isLeaf()) continue;
      MpiVariableMsg &msg=*new MpiVariableMsg(COMPUTE_TABLE);
      msg<<i<<j;
      msg.close();
      slaveDriver->addWork(&msg);
    }
  }

  // Wait for slaves to finish
  Vector<MpiFixedMsg*> results;
  slaveDriver->waitForResults(results);
  Vector<MpiFixedMsg*>::iterator cur=results.begin(), end=results.end();
  for(; cur!=end ; ++cur) delete *cur;
  slaveDriver->terminateSlaves();
}



void LibraryBuilder::initBranches()
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



/****************************************************************
                         slave methods
 ****************************************************************/

int LibraryBuilder::slave(int argc,char *argv[])
{
  processCmdLine(argc,argv);
  
  // Start the timer
  timer.startCounting();

  // Perform gradient ascent
  slave_serveTheMaster();

  // Clean up
  timer.stopCounting();
  MemoryProfiler::report("Total memory used by slave:");
  return 0;
}



void LibraryBuilder::slave_serveTheMaster()
{
  while(1) {
    MpiFixedMsg *msg=slaveDriver->acceptWork();
    switch(msg->getTag())
      {
      case COMPUTE_TABLE:
	slave_computeTable(*msg);
	break;
      case TERMINATE_SLAVE:
	return;
      }
    delete msg;
  }
}



void LibraryBuilder::slave_computeTable(MpiFixedMsg &msg)
{
  int i, j;
  msg>>i>>j;
  Taxon &left=taxa[i], &right=taxa[j];
  //cout<<i<<"x"<<j<<" = "<<left.getName()<<" "<<right.getName()<<endl;
  double dist=tree->distanceBetween(left.getNode(),right.getNode());
  Sequence &seq1=left.getSeq(), &seq2=right.getSeq();
  int L1=seq1.getLength(), L2=seq2.getLength();
  BandingPattern bandingPattern(L2,L1,bandwidth);
  BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
						   dist,NULL);
  if(reestimateEq) hmm->setBgEqFreqs(bgFreqs);
  BandedForward F(*hmm,seq1,seq2,alphabetMap,bandingPattern);
  BandedBackward B(*hmm,seq1,seq2,alphabetMap,bandingPattern);
  delete hmm;
  double LL1=F.getLikelihood(), LL2=B.getLikelihood();
  if(!isFinite(LL1) || !isFinite(LL2)) 
    throw "Can't compute likelihood: banding parameter might be too small";
  if(fabs(LL1-LL2)>1.0e-8) throw "Likelihoods not equal";
  PosteriorMatrix post(F,B);
  //TransPosteriorMatrix transPost(F,B,*hmm,seq1,seq2,post);
  cout<<factor<<" "<<left.getName()<<" "<<right.getName()<<" "<<LL1<<endl;
  float logThreshold=log(sparsenessThreshold);
  String outfile=
    outDir+"/"+factor+"-"+left.getName()+"-"+right.getName()+".post";
  if(!post.saveBinary(logThreshold,outfile)) 
    throw String("Error writing to file: ")+outfile;
  /*
  String transOutfile=
    outDir+"/"+factor+"-"+left.getName()+"-"+right.getName()+".tpost";
  if(!transPost.saveBinary(logThreshold,transOutfile)) 
    throw String("Error writing to file: ")+transOutfile;
  */
  MpiVariableMsg *newMsg=new MpiFixedMsg(DONE_COMPUTING);
  slaveDriver->replyToMaster(newMsg);
}



void LibraryBuilder::estimateBgEq()
{
  Array1D<double> eq(4);
  eq.setAllTo(0.0);
  int N=0;
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    if(!taxon.isLeaf()) continue;
    Sequence &seq=taxon.getSeq();
    int L=seq.getLength();
    for(int j=0 ; j<L ; ++j) {
      eq[seq[j]]+=1.0;
      ++N;
    }
  }
  for(int i=0 ; i<eq.size() ; ++i) eq[i]/=N;
  bgFreqs=eq;
}

