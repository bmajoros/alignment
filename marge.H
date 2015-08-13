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

/****************************************************************
                          enum MessageTag
 ****************************************************************/
enum MessageTag {
  SAMPLE,
  DONE_SAMPLING,
  PRINTOUT
};

/****************************************************************
                          enum EditType
 ****************************************************************/
enum EditType {
  INSERT, DELETE
};
ostream &operator<<(ostream &,const EditType &);

/****************************************************************
                         struct EditEvent
 ****************************************************************/
struct EditEvent {
  EditEvent(EditType t,int p,int m=0) : type(t), pos(p), mapTo(m) {}
  EditType type;
  int pos;
  int mapTo;
  void printOn(ostream &) const;
};
ostream &operator<<(ostream &,const EditEvent &);



/****************************************************************
                     struct pairwiseAlignment
 ****************************************************************/
struct PairwiseAlignment {
  int tax1, tax2;
  IndexMap fwdMap, backMap;
  void invert();
};



/****************************************************************
                     struct LinkedResidue
 ****************************************************************/
enum BranchDirection {
  BRANCH_UP,
  BRANCH_LEFT,
  BRANCH_RIGHT
};
ostream &operator<<(ostream &,BranchDirection);
struct LinkedResidue {
  LinkedResidue();
  LinkedResidue(BOOM::Symbol s);
  LinkedResidue(const LinkedResidue &other);
  void printOn(ostream &) const;
  void removeLinks();
  int pos;
  BOOM::Symbol s;
  LinkedResidue *link[3];
};
ostream &operator<<(ostream &,const LinkedResidue &);
class LinkSeq : public List<LinkedResidue> {
public:
  void assignPositions();
  void removeLinks();
  void printOn(ostream &) const;
};
ostream &operator<<(ostream &,const LinkSeq &);



/****************************************************************
                           Application
 ****************************************************************/
class Application : public AVES {
public:
  Application();
  virtual void progressive();
  int main(int argc,char *argv[]);
protected:
  BOOM::Symbol dot;
  Array1D<LinkSeq> linkSeqs; // indexed by taxon ID
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
  void sample(File &f,ResidueOrthologyGraph &,int sampleNum);
  void loadMatrices(const String &dirname);
  void getStateTypes(Array1D<PHMM_StateType> &stateTypes,
		     Array1D<int> &dI,Array1D<int> &dJ);
  virtual void initBranches();
  PHMM_StateType classifyState(IndexMap &upMap,IndexMap &downMap,
			       int parentPos,int childPos);
  PairwiseAlignment *alignLeaves(int tax1,int tax);
  void samplePair(int &tax1,int &tax2,int numNodes,
		  Array2D<bool> &equivClasses);
  void updateAncestors(PairwiseAlignment &);
  void extendSimple(Taxon &from,Taxon &to,IndexMap &);
  void extend(Taxon &from,Taxon &to,IndexMap &);
  void extendIntoLeaf(Taxon &from,Taxon &to,IndexMap &);
  void initLinkSeqs();
  void allocateLinkSeq(LinkSeq &,int len);
  BranchDirection getDir(Taxon &from,Taxon &to);
  void resolveIndexMaps(BranchAttributes &);
};

