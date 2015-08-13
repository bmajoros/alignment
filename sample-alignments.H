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

enum EditType {
  INSERT, DELETE
};
ostream &operator<<(ostream &,const EditType &);

struct EditEvent {
  EditEvent(EditType t,int p,int m=0) : type(t), pos(p), mapTo(m) {}
  EditType type;
  int pos;
  int mapTo;
  void printOn(ostream &) const;
};
ostream &operator<<(ostream &,const EditEvent &);


class Application : public AVES {
public:
  Application();
  virtual void progressive();
  int main(int argc,char *argv[]);
protected:
  String outfile;
  Array2D<SparseMatrix3D*> matrices;
  Regex filenameRegex, pathRegex;
  float logPseudoCount;
  float indelCoef;
  BranchHMM *hmm;
  int numStates;
  unsigned seed;
  STATE matchState, insertState, deleteState;
  bool wantPrint, wantSemiMarkov, wantTransitivity;
  int skipSamples;
  MpiSlaveDriver *slaveDriver;
  Array1D<PHMM_StateType> stateTypes;
  Array1D<int> dI, dJ;
  Array2D<float> transProbs;
  void processCmdLine(int argc,char *argv[]);
  int master(int argc,char *argv[]);
  int slave(int argc,char *argv[]);
  void slave_serveTheMaster();
  void slave_sample(MpiFixedMsg &msg);
  void samplingLoop();
  void sample(File &f,Vector<PhylogenyBranch*> &,
	      Vector<int> &taxonIndices,ResidueOrthologyGraph &,
	      int sampleNum);
  void sampleBranch(PhylogenyBranch *);
  void sampleBranch2(PhylogenyBranch *);
  void loadMatrices(const String &dirname);
  void getLeaves(Vector<int> &insideLeaves,Vector<int> &outsideLeaves,
		 Taxon &parent,Taxon &child);
  void getStateTypes(Array1D<PHMM_StateType> &stateTypes,
		     Array1D<int> &dI,Array1D<int> &dJ);
  void getHMMprobs();
  void getHMMprobs(PairHMM *);
  void getLeafVector(BitSet &cladeMembers,Vector<int> &leafTaxa);
  float computeCellP(Vector<int> &insideLeaves,Vector<int> &outsideLeaves,
		     int numInsideLeaves,int numOutsideLeaves,
		     CollapsedOrthologyMatrix &childCOM,
		     CollapsedOrthologyMatrix &parentCOM,STATE q,
		     int newI,int newJ,bool &justPseudo,
		     Vector<bool> &parentReachability,
		     Vector<bool> &childReachability);
  float debugCellP(Vector<int> &insideLeaves,Vector<int> &outsideLeaves,
		   int numInsideLeaves,int numOutsideLeaves,
		   CollapsedOrthologyMatrix &childCOM,
		   CollapsedOrthologyMatrix &parentCOM,STATE q,
		   int newI,int newJ,bool &justPseudo,
		   Vector<bool> &parentReachability,
		   Vector<bool> &childReachability);
  virtual void initBranches();
  virtual void updateIndels(BranchAttributes &,Taxon &parent,
			    Taxon &child,StatePath&);
  void lookahead(Vector<float> &probs,Vector<STATE> &states,
		 Vector<int> &insideLeaves,Vector<int> &outsideLeaves,
		 int numInsideLeaves,int numOutsideLeaves,
		 CollapsedOrthologyMatrix &childCOM,
		 CollapsedOrthologyMatrix &parentCOM,
		 int parentL,int childL,int i,int j,PHMM_StateType skip);
  void updateTentacles(BranchAttributes &branch,Taxon &parent,
		       Taxon &child,
		       CollapsedOrthologyMatrix &parentCOM,
		       CollapsedOrthologyMatrix &childCOM,
		       int parentL,int childL,
		       Vector<int> &outsideLeaves,
		       Vector<int> &insideLeaves,
		       int numOutsideLeaves,int numInsideLeaves,
		       StatePath &path);
  bool reachesLeaf(CollapsedOrthologyMatrix &,Vector<int> &leaves,int pos);
  void computeReachabilityVector(CollapsedOrthologyMatrix &,
				 Vector<int> &leaves,Vector<bool> &into);
  bool isDeadEnd(Taxon &,int pos,InsideOutside,Taxon *exclude=NULL);
  PHMM_StateType classifyState(IndexMap &upMap,IndexMap &downMap,
			       int parentPos,int childPos);
  int getNewLength(int oldLen,Vector<EditEvent> &);
  void computeSeqMap(int oldLen,Vector<EditEvent> &,IndexMap &);
  void accommodateChildEdits(Taxon &child,BranchAttributes &,
			     Vector<EditEvent> &childEvents,
			     IndexMap &upMap,IndexMap &downMap,
			     IndexMap &childSeqMap);
  void accommodateParentEdits(Taxon &parent,Taxon &child,BranchAttributes &,
			      Vector<EditEvent> &parentEvents, 
			      IndexMap &upMap,IndexMap &downMap,
			      const IndexMap &childSeqMap);
  int findStartingPos(int epos,IndexMap &seqMap);
  CollapsedOrthologyMatrix::Entry indexCOM(CollapsedOrthologyMatrix &COM,
					   int pos,int leaf);
  float sandwichProb(InsideOutside io,PHMM_StateType stateType,
		     Vector<int> &insideLeaves,Vector<int> &outsideLeaves,
		     CollapsedOrthologyMatrix &childCOM,
		     CollapsedOrthologyMatrix &parentCOM,
		     int childPos,int parentPos,bool &justPseudo);  
  bool semiMarkov(StatePath &path,int &i,int &j,Vector<int> &insideLeaves,
		  Vector<int> &outsideLeaves,int numInsideLeaves,
		  int numOutsideLeaves,CollapsedOrthologyMatrix &childCOM,
		  CollapsedOrthologyMatrix &parentCOM,int parentL,int childL,
		  Vector<bool> &parentReachability,
		  Vector<bool> &childReachability);
  bool manhattan(StatePath &path,int &i,int &j,Vector<int> &insideLeaves,
		  Vector<int> &outsideLeaves,int numInsideLeaves,
		  int numOutsideLeaves,CollapsedOrthologyMatrix &childCOM,
		  CollapsedOrthologyMatrix &parentCOM,int parentL,int childL,
		  Vector<bool> &parentReachability,
		  Vector<bool> &childReachability);
  bool cellIsZero(Vector<int> &insideLeaves,
		   Vector<int> &outsideLeaves,
		   int numInsideLeaves,int numOutsideLeaves,
		   CollapsedOrthologyMatrix &childCOM,
		   CollapsedOrthologyMatrix &parentCOM,STATE q,
		   int newI,int newJ,
		   Vector<bool> &parentReachability,
		   Vector<bool> &childReachability);
  void jumpToCell(StatePath &,int &i,int &j,int newI,int newJ);
};

