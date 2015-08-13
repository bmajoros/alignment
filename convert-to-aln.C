#include <iostream>
#include <fstream>
#include "AVES.H"
#include "ResidueOrthologyGraph.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;

class Application : public AVES {
  virtual void initBranches();
  virtual void inferGapPatterns(const MultSeqAlignment &);
  virtual void loadSeedAlignment(const String &filename);
  virtual void reconstructHistories();
public:
  Application() {}
  virtual int main(int argc,char *argv[]);
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




int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4)
    throw String(
"\nconvert-to-aln <in.maf> <in.fasta> <tree.phy> <out.aln>\n\n");
  String alignFile=cmd.arg(0);
  String fastaFile=cmd.arg(1);
  String lambdaFile=cmd.arg(2);
  outfile=cmd.arg(3);
  baselineMode=false;
  disableGC=cmd.option('G');
  if(disableGC) modelCompiler->setGCthreshold(LARGEST_INTEGER);
  else modelCompiler->getLambda().getGC().setSilence(true);
  contentSensor=NULL;
  
  // Load phylogeny
  tree=new Phylogeny(lambdaFile);
  tree->constructBranches();
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
  alignmentBuilder=
    new AlignmentBuilder(tree->getRoot(),alphabet,gapSymbol,numTaxa);
  gatherCladeMembers();
  initBranches();
  loadSeedAlignment(alignFile);

  // ### other stuff needs to be done here...otherwise probably segfault...
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
  
  reconstructHistories();
  ResidueOrthologyGraph rog(tree,taxa);
  File f(outfile,"w");
  rog.saveConnections(f);

  //for(int i=0 ; i<numTaxa ; ++i) cout<<taxa[i].getName()<<" => "<<taxa[i].getSeq()(alphabet)<<endl;

  /*
  delete fullAlignment;
  fullAlignment=alignmentBuilder->buildAlignment(true,false);
  fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',
			    60,true,true);
  */
}



void Application::reconstructHistories()
{
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &parent=taxa[i];
    int numBranches=parent.getNumBranches();
    for(int j=0 ; j<numBranches ; ++j) {
      BranchAttributes *branch=parent.getIthBranch(j);
      Taxon &child=branch->getChildTaxon();
      reconstructHistory(parent,child,branch->getDownMap(),
			 branch->getUpMap(),false);
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
    /*
    BranchHMM *hmm=templateInstantiator->instantiate(transducerTemplate,
						     branch->getLength(),NULL);
    branchAttributes[i]=BranchAttributes(hmm,branch);
    */
    branchAttributes[i]=BranchAttributes(NULL,branch);
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
                     AVES::loadSeedAlignment()
 ****************************************************************/
void Application::loadSeedAlignment(const String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw String("Error opening file: ")+filename;
  MultiAlignment *maf=MultiAlignment::nextAlignmentFromMAF(is);
  fullAlignment=new MultSeqAlignment(*maf,alphabet,gapSymbol);
  inferGapPatterns(*fullAlignment);
  for(int i=0 ; i<numTaxa ; ++i) {
    Taxon &taxon=taxa[i];
    AlignmentSeq &track=fullAlignment->getTrackByName(taxon.getName());
    const Sequence &S=track.getSeq();
    taxon.getSeq()=S;
    taxon.setSeqLen(S.getLength());
    taxon.getSeqStr()=S(alphabet);
    //cout<<taxon.getName()<<" "<<taxon.getSeqStr()<<endl;
  }
  delete maf;
  is.close();
  taxa[tree->getRoot()->getID()].setCladeAlignment(fullAlignment);
}



/****************************************************************
                     AVES::inferGapPatterns()
 ****************************************************************/
void Application::inferGapPatterns(const MultSeqAlignment &A)
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
      case INTERNAL_NODE:
	if(!A.trackExists(taxonName)) break;
	// fall through...
      case LEAF_NODE: {
	AlignmentSeq &oldTrack=A.getTrackByName(taxonName);
	for(int j=0 ; j<L ; ++j) 
	  //newTrack[j]=(oldTrack[j]==gapSymbol ? gap : residue);
	  newTrack[j]=(oldTrack[j]>3 ? gap : residue); //###
	break; }
      }
  }	

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
