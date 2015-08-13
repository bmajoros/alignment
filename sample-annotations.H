#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Regex.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/MPI.H"
#include "AVES.H"
#include "SparseMatrix3D.H"
#include "BranchHMM.H"
#include "ResidueOrthologyGraph.H"
#include "CollapsedOrthologyMatrix.H"
using namespace std;
using namespace BOOM;

typedef SparseMatrix3D::EntryList EntryList;

enum MessageTag {
  SAMPLE,
  DONE_SAMPLING,
  PRINTOUT
};

class Application : public AVES {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  String outfile;
  BranchHMM *hmm;
  int numStates;
  unsigned seed;
  bool wantPrint;
  bool useGFF;
  int skipSamples, sampleID;
  MpiSlaveDriver *slaveDriver;
  String alnFile;
  Array1D<PHMM_StateType> stateTypes;
  Array1D<int> dI, dJ;
  Array2D<float> transProbs;
  Vector<STATE> QI, QD, QM; // insert/delete/match states, respectively
  void processCmdLine(int argc,char *argv[]);
  int master(int argc,char *argv[]);
  int slave(int argc,char *argv[]);
  void slave_serveTheMaster();
  void slave_sample(MpiFixedMsg &msg);
  void samplingLoop();
  bool sample(ostream &f,Vector<PhylogenyBranch*> &,ResidueOrthologyGraph &);
  bool sampleBranch(PhylogenyBranch *,ostream &,Vector<StatePath*> &paths);
  void getLeaves(Vector<int> &insideLeaves,Vector<int> &outsideLeaves,
		 Taxon &parent,Taxon &child);
  void getStateTypes(Array1D<PHMM_StateType> &stateTypes,
		     Array1D<int> &dI,Array1D<int> &dJ);
  void getLeafVector(BitSet &cladeMembers,Vector<int> &leafTaxa);
  virtual void initBranches();
  PHMM_StateType classifyState(IndexMap &upMap,IndexMap &downMap,
			       int parentPos,int childPos);
  void savePath(StatePath &,ostream &,Taxon &parent,Taxon &child,
		BranchHMM &);
  void installFuncParse(StatePath &,ostream &,Taxon &parent,Taxon &child,
		BranchHMM &);
};

