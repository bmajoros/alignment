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
                         class Application
 ****************************************************************/
class Application
{
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
  bool pairwise; // whether to estimate initial parms via pairwise leaves
  Map<String,int> nameToTaxon; // index into "taxa" array
  int numTaxa, alignmentLength;
  AlignmentBuilder *alignmentBuilder;
  void gatherCladeMembers();
  void initBranches();
  void sumBranches(double &sumAlpha,double &sumBeta,double &sumT);
  void sumPairwise(double &sumAlpha,double &sumBeta,double &sumT);
  void sumPair(Taxon *parent,Taxon *child,double &sumAlpha,double &sumBeta,
	       float branchLength);
public:
  Application();
  int main(int argc,char *argv[]);
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
    functionalParse(NULL), alignmentBuilder(NULL), pairwise(false)
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
  // Process command line
  CommandLine cmd(argc,argv,"p");
  if(cmd.numArgs()!=3)
    throw String("\n\
train-background <model.lambda> <training-alignment.maf> <outfile.lambda>\n\
    where:  -p = use pairs of leaves to estimate initial parm values\n\
");
  String lambdaFile=cmd.arg(0);
  String trainAlignFile=cmd.arg(1);
  String outfile=cmd.arg(2);
  //Lambda::GarbageCollector::setThreshold(10); // ###
  if(cmd.option('p')) pairwise=true;

  // Load HMM
  cout<<"loading HMM"<<endl;
  modelCompiler->parse(lambdaFile);
  //if(modelCompiler->getNumModels()<1)
  //throw "use (register-transducer ...) to register a transducer";
  //transducerTemplate=modelCompiler->getIthModel(0);
  Lambda::Closure *ctor=modelCompiler->getCtor(0);
  transducerTemplate=modelCompiler->instantiate(ctor);

  // Load phylogeny
  cout<<"loading phylogeny"<<endl;
  tree=modelCompiler->getGuideTree(transducerTemplate);
  tree->gatherNodes(phylogenyNodes);
  tree->gatherBranches(branches);

  // Link taxa to phylogeny nodes
  cout<<"linking taxa to phylogeny"<<endl;
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
  cout<<"loading alignment"<<endl;
  ifstream is(trainAlignFile.c_str());
  if(!is.good()) throw String("Error opening file: ")+trainAlignFile;
  fullAlignment=NULL;
  while(!is.eof()) {
    MultiAlignment *maf=MultiAlignment::nextAlignmentFromMAF(is);
    if(!maf) break;
    if(!fullAlignment) 
      fullAlignment=new MultSeqAlignment(*maf,alphabet,gapSymbol);
    else {
      MultSeqAlignment temp(*maf,alphabet,gapSymbol);
      fullAlignment->append(temp);
    }
    delete maf;
  }
  is.close();
  alignmentLength=fullAlignment->getLength();

  // Infer ancestral sequences via Dollo parsimony
  cout<<"inferring ancestral sequences"<<endl;
  if(!pairwise) {
    cout<<"inferring ancestral sequences"<<endl;
    SingleGainParsimony parsimony(*tree,alphabet,*fullAlignment,gapSymbols);
    cout<<"dollo.run()"<<endl;
    parsimony.run();
  }
  //fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,true);

  // Estimate initial values for parameters
  cout<<"estimating initial parm values"<<endl;
  double sumT=0, sumAlpha=0, sumBeta=0;
  if(pairwise) sumPairwise(sumAlpha,sumBeta,sumT);
  else sumBranches(sumAlpha,sumBeta,sumT);
  double aveAlpha=sumAlpha/sumT, aveBeta=sumBeta/sumT;
  cout<<"ave(alpha)="<<aveAlpha<<" ave(beta)="<<aveBeta<<endl;
  float lambda=modelCompiler->lookupFloat("lambda");
  float mu=modelCompiler->lookupFloat("mu");

  // Write output to file
  ofstream os(outfile.c_str());
  os<<"(define 'alpha "<<aveAlpha<<")"<<endl;
  os<<"(define 'beta "<<aveBeta<<")"<<endl;
  os<<"(define 'lambda "<<lambda<<")"<<endl;
  os<<"(define 'mu "<<mu<<")"<<endl;
  
  return 0;
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



void Application::sumPair(Taxon *parent,Taxon *child,double &sumAlpha,
			  double &sumBeta,float branchLength)
{
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
  float alpha=-log(x)/branchLength;
  float beta=-log(1-y)/branchLength;
  sumAlpha+=alpha*branchLength; // ### weighted average
  sumBeta+=beta*branchLength; // ### weighted average
  cout<<"t="<<branchLength<<" alpha="<<alpha<<" beta="<<beta<<" x="<<x<<" y="<<y<<endl;
}



void Application::sumBranches(double &sumAlpha,double &sumBeta,double &sumT)
{
  int numBranches=branches.size();
  for(int i=0 ; i<numBranches ; ++i) {
    PhylogenyBranch *branch=&branches[i];
    float branchLength=branch->getLength();
    sumT+=branchLength;
    Taxon *parent=static_cast<Taxon*>(branch->getParent()->getDecoration());
    Taxon *child=static_cast<Taxon*>(branch->getChild()->getDecoration());
    sumPair(parent,child,sumAlpha,sumBeta,branchLength);
  }
}



void Application::sumPairwise(double &sumAlpha,double &sumBeta,double &sumT)
{
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon *a=&taxa[i];
    if(!a->isLeaf()) continue;
    for(int j=i+1 ; j<numTaxa ; ++j) {
      Taxon *b=&taxa[j];
      if(!b->isLeaf()) continue;
      float branchLength=tree->distanceBetween(a->getNode(),b->getNode());
      sumT+=branchLength;
      sumPair(a,b,sumAlpha,sumBeta,branchLength);
    }
  }
}





