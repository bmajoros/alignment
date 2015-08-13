/****************************************************************
 AVES.H : Alignment Via Evolutionary Sampling
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/

#ifndef INCL_AVES_H
#define INCL_AVES_H

#include <iostream>
#include <fstream>
#include "LinkParsimony.H"
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/DnaDashAlphabet.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/Map.H"
#include "BOOM/MPI.H"
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/Random.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Time.H"
#include "BOOM/GffReader.H"
#include "BOOM/MemoryProfiler.H"
#include "BOOM/Constants.H"
#include "PairHMM/Viterbi.H"
#include "PairHMM/PairHMM.H"
#include "PairHMM/Transducer.H"
#include "PhyLib/Phylogeny.H"
#include "PhyLib/FitchParsimony.H"
#include "PhyLib/SingleGainParsimony.H"
#include "EGGS/ContentSensor.H"
#include "Taxon.H"
#include "ProfileViterbi.H"
#include "AlignmentView.H"
#include "ProfileSampler.H"
#include "ProfileFelsenstein.H"
#include "LinkFelsenstein.H"
#include "LinkBackward.H"
#include "LinkSampler.H"
#include "GapPattern.H"
#include "BranchAttributes.H"
#include "ModelCompiler.H"
#include "FunctionalParse.H"
#include "AlignmentBuilder.H"
#include "Banding.H"
#include "Posteriors.H"
#include "LossyFelsenstein.H"
#include "SiblingViterbi.H"
#include "GainLossFelsenstein.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;

enum SamplerType { GIBBS, METROPOLIS_HASTINGS };
enum ScoringStrategy { SS_LLR, SS_POSTERIOR, SS_NONE };
enum AvesMessageTag { MSG_ALIGN_SIBS, MSG_STATE_PATH };

class AVES {
protected:
  SamplerType samplerType;
  int bandwidth; // for constant-width banding
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
  //DnaDashAlphabet alphabet;
  DnaDashDotAlphabet alphabet;
  PureDnaAlphabet pureDnaAlphabet;
  DropGapMapping alphabetMap;
  AlphabetIdentityMap identityMap;
  BOOM::Symbol gapSymbol;
  ContentSensor *contentSensor;
  ScoringStrategy scoringStrategy;
  String outfile;
  bool shouldSample, shouldHillClimb, shouldSeed; // command-line flags
  bool explicitHistories, shouldDumpModel; // more command-line flags
  bool baselineMode, useHirschberg;
  bool usePrescan; // whether to limit search space to pre-scanned TFBS's
  bool noPrescanOnDownpass; 
  bool useContentSensor; // use a Markov chain instead of simple EQ freqs
  bool disableGC;
  bool wantFuncDollo;
  BandingType bandingType;
  int numSamples; // from command line
  Map<String,int> nameToTaxon; // index into "taxa" array
  int numTaxa;
  int refinementIterations;
  int numThreads;
  MPI *mpi;
  MpiSlaveDriver *slaveDriver;
  double alignmentScore, proposalScore; // for current sampling epoch
  AlignmentBuilder *alignmentBuilder;
  Taxon *targetSpecies;
  String dumpFile; // for dumping the model parameters
  String gffFile, historiesFile;
  bool usePosteriors;
  Posteriors *posteriors;
  Time timer;
  BranchHMM *defaultHMM;
  Array1D<bool> isNucleotide;
  virtual void progressive();
  void uppass(LinkFelsenstein &);
  void align(PhylogenyNode *left,PhylogenyNode *right,LinkFelsenstein &);
  void dollo();
  void localize();
  void downpass(LinkFelsenstein &);
  void downpass(Taxon &parent,LinkFelsenstein &,BranchAttributes *);
  virtual void updateIndels(BranchAttributes &,Taxon &parent,
			    Taxon &child,StatePath&);
  void inferRootFuncParse(LinkFelsenstein &);
  void completeOrthology(MultSeqAlignment *);
  void iterativeRefinement(LinkFelsenstein &);
  void refine(BranchAttributes &,LinkFelsenstein &);
  void propagate(Taxon *from,Taxon *to,LinkFelsenstein &F);
  int slave(int argc,char *argv[]);
  int master(int argc,char *argv[]);
  void emit(const String &outfile);
  void installSiblingHistories(StatePath &path,Taxon &parent,
			       BranchAttributes *leftBranch,
			       BranchAttributes *rightBranch,
			       BranchHMM *hmm);
  void inferParentChildParse(LinkFelsenstein &,Taxon &parent,
			     Taxon &child);
  void connect(IndexMap &parent,int parentIdx,IndexMap &child,int childIdx);
  MpiVariableMsg *packPath(StatePath &,MPI_MESSAGE_TAG,int leftID,int rightID);
  StatePath *unpackPath(MpiFixedMsg &,int &leftID,int &rightID);
  void uppass2(LinkFelsenstein &);
  void installStatePaths(Vector<MpiFixedMsg*> &,BitSet &resolved);
  void installStatePath(StatePath &,int leftID,int rightID);
  void progressiveRoot(PhylogenyNode *);
  void gatherCladeMembers();
  virtual void initBranches();
  void reinitBranches();
  void sample();
  virtual void inferGapPatterns(const MultSeqAlignment &);
  double likelihood();
  virtual void loadSeedAlignment(const String &filename);
  void loadFcConstraints(const String &filename);
  int processCmdLine(int argc,char *argv[]);
  virtual void reconstructHistories();
  void reconstructHistory(Taxon &parent,Taxon &child,IndexMap &downMap,
			  IndexMap &upMap,bool noParentSeqs=true);
  void reconstructHistory(StatePath &,PairHMM &,IndexMap &,IndexMap &);
  void rebuildPrescans(Phylogeny *);
  void installAnnotationTrack(MultSeqAlignment &,FunctionalParse &);
  void reconstructAlignment(bool includeAncestors=false);
  void updatePrescans(Taxon &);
  double sumTransLogProbs();
  void updateFcConstraints(StatePath &,Taxon &);
  double computeScore(Taxon &,FunctionalElement &);
  double computeScoreLLR(Taxon &,FunctionalElement &);
  double computeScorePosterior(Taxon &,FunctionalElement &);
  double computeScore(PHMM_StateType,SubstitutionMatrix &,
		      const Array1D<double> &eqFreqs,int parentPos,
		      int childPos,PrecomputedEmissions &parentE,
		      PrecomputedEmissions &childE,
		      FunctionalClass fcParent,FunctionalClass fcChild,
		      Taxon &parent,Taxon &child);
  double computePostNumerator(Taxon &,FunctionalElement &);
  double computePostDenominator(Taxon &);
  double getEQ(const Array1D<double> &eqFreqs,BOOM::Symbol,Taxon *,int pos,
	       FunctionalClass);
public:
  AVES();
  virtual int main(int argc,char *argv[]);
};



#endif

