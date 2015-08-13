#ifndef INCL_FineTuner_H
#define INCL_FineTuner_H

inline double mapToProbability(double x);
inline double mapFromProbability(double x) ;
GSL::Vector mapToProbabilities(const GSL::Vector &v,int numBranches);
GSL::Vector mapFromProbabilities(const GSL::Vector &v,int numBranches);


/****************************************************************
                          enum MessageTag
 ****************************************************************/
enum MessageTag {
  COMPUTE_LIKELIHOOD, // master -> slave  (all pairs)
  LIKELIHOOD,         // slave -> master
  PAIR_LL             // master -> slave  (just one pair)
};


/****************************************************************
                          enum PairType
 ****************************************************************/
enum PairType {
  PT_ONE_PAIR, // not implemented
  PT_ALL_PAIRS,// not implemented
  PT_ONE_PAIR_PER_CLADE,
  PT_N_OVER_2, // not implemented
  PT_UNIFORM_LENGTHS
};



/****************************************************************
                        struct SpeciesPair
 ****************************************************************/
struct SpeciesPair {
  Taxon *first, *second;
  double distance;
  BranchHMM *hmm;
  SpeciesPair(Taxon *f,Taxon *s,BranchHMM *h) 
    : first(f), second(s), hmm(h) {}
};



/****************************************************************
                        struct SpeciesTriplet
 ****************************************************************/
struct SpeciesTriplet {
  int ID[3];
  float distance;
  Phylogeny *tree;
  Array1D<Taxon> taxa;
  SpeciesTriplet(int id0,int id1,int id2,float dist)
    : distance(dist), tree(NULL) { ID[0]=id0; ID[1]=id1; ID[2]=id2; }
  void initTaxa(Array1D<Taxon> &);
  int find(const String &name,Array1D<Taxon> &);
};



/****************************************************************
                         class FineTuner
 ****************************************************************/
class FineTuner : public AVES {
  Vector<ModelParm> parms;
  String lambdaFile;
  GSL::Optimizer *optimizer;
  Vector<STATE> QI, QD, QM;
  int numI, numD, numM;
  String outfile;
  String optimizerTypeStr, stoppingCriterionStr;
  double tolerance, gradientThreshold, gradientEpsilon, stepSize;
  int numAlpha, dimensionality;
  PairType pairType;
  Vector<SpeciesPair> pairs;
  MpiSlaveDriver *slaveDriver;
  Vector<SpeciesTriplet*> triplets;
  bool useTriplets;
  void fineTune();
  void onePairPerClade();
  void uniformLengths();
  Taxon *findLeaf(Taxon *);
  void constructPairs();
  void collectModelParms();
  virtual void initBranches();
  int master(int argc,char *argv[]);
  int slave(int argc,char *argv[]);
  void serveTheMaster();
  void slaveLikelihood(MpiFixedMsg &);
  void slavePairLL(MpiFixedMsg &);
  void writeFile(const GSL::Vector &);
  void collectTriplets(int numTriplets);
  void uppass(LinkFelsenstein &F,SpeciesTriplet *triplet);
  void progressiveRoot(PhylogenyNode *node,SpeciesTriplet *triplet);
  void align(PhylogenyNode *leftNode,PhylogenyNode *rightNode,
	     LinkFelsenstein &F,SpeciesTriplet *triplet);
  void dollo(SpeciesTriplet *triplet);
  void reconstructHistories(SpeciesTriplet *triplet);
  void inferRootFuncParse(LinkFelsenstein &F,SpeciesTriplet *triplet);
  void downpass(LinkFelsenstein &F,SpeciesTriplet *triplet);
  void downpass(Taxon &parent,LinkFelsenstein &,BranchAttributes *,
		SpeciesTriplet *);
  double likelihood(SpeciesTriplet *triplet);
  void inferGapPatterns(const MultSeqAlignment &,SpeciesTriplet *);
  double sumTransLogProbs(SpeciesTriplet *);
public:
  FineTuner();
  virtual int main(int argc,char *argv[]);
  void copyIn(GSL::Vector &);
  void copyOut(const GSL::Vector &);
  double pairwiseLikelihood(int whichPair=-1);
  double tripletLikelihood(int whichTriplet=-1);
  double distributedLikelihood(const GSL::Vector &);
  Vector<ModelParm> &getParms() {return parms;}
  void reparse();
  void computeGradient(const GSL::Vector &currentPt,
		       GSL::Vector &gradient);
};



/****************************************************************
                     class ObjectiveFunction
 ****************************************************************/
class ObjectiveFunction : public GSL::ObjectiveFunction
{
  FineTuner *app;
  double worstLikelihood, likelihood;
  GSL::Optimizer **optimizer;
  BOOM::Map<const GSL::Vector,double> cache;
  int cacheHits;
  int numEvaluations;
  double firstScore, bestScore;
  GSL::Vector bestPoint;
  Vector<ModelParm> &parms;
  
  void twoPointSymmetric(const GSL::Vector &currentPoint,
			  GSL::Vector &gradient);
  bool pointIsValid(const GSL::Vector &point);
public:
  ObjectiveFunction(FineTuner *app);
  virtual double f(const GSL::Vector &currentPoint);
  virtual void gradient(const GSL::Vector &currentPoint,
			GSL::Vector &gradient);
  int getNumEvaluations() const {return numEvaluations;}
  int getNumCacheHits() const {return cacheHits;}
  void copyIn(GSL::Vector &);
  void copyOut(const GSL::Vector &);
  double getFirstScore() {return firstScore;}
  double getBestScore() {return bestScore;}
  GSL::Vector &getBestPoint() {return bestPoint;}
};


#endif

