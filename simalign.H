#ifndef INCL_simalign_H
#define INCL_simalign_H
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
#include "PairHMM/PairHMM.H"
#include "PairHMM/Transducer.H"
#include "PhyLib/Phylogeny.H"
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



class Simalign
{
  int rootLength;
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
  bool leavesOnly, shouldDumpModel; // command-line options
  Map<String,int> nameToTaxon; // index into "taxa" array
  int numTaxa;
  AlignmentBuilder *alignmentBuilder;
  String dumpFile;
  int totalSites, lossEvents, gainEvents;
  Vector<float> degreesOfOrthology;
  Array1D<int> gainsLosses;
  void emit(const String &outfile,const String &mafFile,
	    const String &annoFile);
  void gatherCladeMembers();
  void initBranches();
  void sample();
  void installAnnotationTrack(MultSeqAlignment &,FunctionalParse &);
  void reconstructAlignment();
  void generateRootSeq();
  void evolve(Taxon &);
  void evolve();
  BOOM::Symbol sampleResidue(const Array1D<double> &eqFreqs);
  void sampleChildSeq(BranchAttributes &,Taxon &parent,Taxon &child,
		      StatePath &,BranchHMM &hmm);
  void reconstructHistory(StatePath &,PairHMM &,IndexMap &upMap,
			  IndexMap &downMap);
  void alignmentStats(MultSeqAlignment &,ostream &);
  void countGainLoss(StatePath &path,BranchHMM &);

public:
  Simalign();
  int main(int argc,char *argv[]);
};



#endif
