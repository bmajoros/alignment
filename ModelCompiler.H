/****************************************************************
 ModelCompiler.H
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#ifndef INCL_ModelCompiler_H
#define INCL_ModelCompiler_H
#include <iostream>
#include "PhyLib/Phylogeny.H"
#include "BOOM/String.H"
#include "BOOM/Vector.H"
#include "TransducerTemplate.H"
#include "TemplateInstantiator.H"
#include "BOOM/lambda/LambdaAPI.H"
#include "GainLossType.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                   lambda "foreign objects"
 ****************************************************************/

enum ForeignObjectType {
  FOT_RATE_MATRIX,
  FOT_RATE_MATRICES,
  FOT_PHYLOGENY,
  FOT_STATES,
  FOT_STATE,
  FOT_TRANSITION,
  FOT_TRANSDUCER,
  FOT_STATE_TYPE,
  FOT_MACROSTATE,
  FOT_MACROSTATES,
  FOT_MACROTRANS,
  FOT_FUNC_CLASS,
  FOT_FUNC_ELEM_TYPE,
  FOT_GAINLOSS_FACTOR
};

ostream &operator<<(ostream &,ForeignObjectType);

struct FO_RateMatrix : public ForeignObject {
  FO_RateMatrix(RateMatrix *m) 
    : ForeignObject(reinterpret_cast<void*>(m),FOT_RATE_MATRIX) {}
  RateMatrix *getMatrix() const 
    { return reinterpret_cast<RateMatrix*>(getObject()); }
};

struct FO_RateMatrices : public ForeignObject {
  FO_RateMatrices(int numMatrices)
    : ForeignObject(NULL,FOT_RATE_MATRICES), matrices(numMatrices) {}
  void setMatrix(int i,FO_RateMatrix *m) { matrices[i]=m; }
  int getNumMatrices() const {return matrices.size();}
  FO_RateMatrix *getMatrix(int i) const { return matrices[i]; }
  virtual void pushAccessibles(MarkStack &s);
protected:
  Array1D<FO_RateMatrix*> matrices;
};

struct FO_Phylogeny : public ForeignObject {
  FO_Phylogeny(Phylogeny *p) 
    : ForeignObject(reinterpret_cast<void*>(p),FOT_PHYLOGENY) {}
  Phylogeny *getPhylogeny() const 
    { return reinterpret_cast<Phylogeny*>(getObject()); }
};

struct FO_State : public ForeignObject {
  FO_State(StateTemplate *t)
    : ForeignObject(reinterpret_cast<void*>(t),FOT_STATE) {}
  StateTemplate *getState()
    { return reinterpret_cast<StateTemplate*>(getObject()); }
};

struct FO_States : public ForeignObject {
  FO_States(Array1D<FO_State*> &a)
    : ForeignObject(reinterpret_cast<void*>(&a),FOT_STATES) {}
  Array1D<FO_State*> &getStates() 
  { return *reinterpret_cast<Array1D<FO_State*>*>(getObject()); }
  virtual ~FO_States() {}//delete &getStates();}
  virtual void pushAccessibles(MarkStack &s);
};

class FO_StateType : public ForeignObject {
  PairHMMStateType stateType;
public:
  FO_StateType(PairHMMStateType t) 
    : ForeignObject(NULL,FOT_STATE_TYPE), stateType(t) {}
  PairHMMStateType getStateType() {return stateType;}
};

struct FO_Transition : public ForeignObject {
  FO_Transition(TransitionTemplate *t)
    : ForeignObject(reinterpret_cast<void*>(t),FOT_TRANSITION) {}
  TransitionTemplate *getTransition()
    { return reinterpret_cast<TransitionTemplate*>(getObject()); }
};

struct FO_Transducer : public ForeignObject {
  FO_Transducer(TransducerTemplate *t)
    : ForeignObject(reinterpret_cast<void*>(t),FOT_TRANSDUCER) {}
  TransducerTemplate *getTemplate()
    { return reinterpret_cast<TransducerTemplate*>(getObject());}
};

struct FO_Macrostate : public ForeignObject {
  FO_Macrostate(MacrostateTemplate *t)
    : ForeignObject(reinterpret_cast<void*>(t),FOT_MACROSTATE) {}
  MacrostateTemplate *getState()
    { return reinterpret_cast<MacrostateTemplate*>(getObject()); }
  virtual ~FO_Macrostate() {}//delete getState();}
};

struct FO_Macrostates : public ForeignObject {
  FO_Macrostates(Array1D<FO_Macrostate*> &a)
    : ForeignObject(reinterpret_cast<void*>(&a),FOT_MACROSTATES) {}
  Array1D<FO_Macrostate*> &getStates() 
  { return *reinterpret_cast<Array1D<FO_Macrostate*>*>(getObject()); }
  virtual ~FO_Macrostates() {}//delete &getStates();}
  virtual void pushAccessibles(MarkStack &s);
};

struct FO_Macrotrans : public ForeignObject {
  FO_Macrotrans(MacrotransTemplate *t)
    : ForeignObject(reinterpret_cast<void*>(t),FOT_MACROTRANS) {}
  MacrotransTemplate *getTransition()
    { return reinterpret_cast<MacrotransTemplate*>(getObject()); }
  virtual ~FO_Macrotrans() {}//delete getTransition();}
};

class FO_FuncClass : public ForeignObject {
  RateMatrix *Q;
  Phylogeny *phy;
  char label;
  String name;
  FunctionalClass fc;
public:
  FO_FuncClass(RateMatrix *Q,Phylogeny *phy,char label,String name,
	       FunctionalClass fc) 
    : ForeignObject(NULL,FOT_FUNC_CLASS), Q(Q), phy(phy), label(label), 
      name(name), fc(fc) {}
  RateMatrix *getMatrix() const {return Q;}
  Phylogeny *getPhylogeny() const {return phy;}
  char getLabel() const {return label;}
  const String &getName() const {return name;}
  FunctionalClass getClass() const {return fc;}
};

class FO_FuncElemType : public ForeignObject {
  Array1D<FO_FuncClass*> classes;
public:
  FO_FuncElemType(int numClasses) : ForeignObject(NULL,FOT_FUNC_ELEM_TYPE), 
				    classes(numClasses) {}
  void setClass(int i,FO_FuncClass *c) {classes[i]=c;}
  int getNumClasses() const {return classes.size();}
  FO_FuncClass *getIthClass(int i) {return classes[i];}
  virtual void pushAccessibles(MarkStack &s);
};

class FO_GainLossFactor : public ForeignObject {
  Lambda::Closure *factor;
  GainLossType t;
public:
  FO_GainLossFactor(Lambda::Closure *factor,GainLossType t) 
    : ForeignObject(NULL,FOT_GAINLOSS_FACTOR), factor(factor), t(t) {}
  Lambda::Closure *getFactor() const {return factor;}
  GainLossType getType() const {return t;}
  virtual void pushAccessibles(MarkStack &s) {s.push(factor);}
};



/****************************************************************
                         struct ModelParm
****************************************************************/
class ModelParm {
  String name;
  float value;
  float min, max, range;
public:
  ModelParm() {}
  ModelParm(String n,float min,float max);
  float constrain(float);
  float unconstrain(float);
  float getValue() {return value;}
  void setValue(float x) {value=x;}
  float getMin() {return min;}
  float getMax() {return max;}
  const String &getName() {return name;}
  bool isValid(float x) {return x>=min && x<=max;}
};



/****************************************************************
                       class ModelCompiler
 ****************************************************************/
class ModelCompiler : public MarkHook
{
public:
  ModelCompiler();
  virtual ~ModelCompiler();
  void parse(const String &filename);
  int getNumModels() const;
  TransducerTemplate *getIthModel(int i);
  TemplateInstantiator *getInstantiator();
  Phylogeny *getGuideTree(); // returns any phylogeny it can find
  Phylogeny *getGuideTree(TransducerTemplate *);
  static ModelCompiler *getGlobal();
  float lookupFloat(const String &variableName);
  void setGCthreshold(int);
  LambdaAPI &getLambda() {return lambda;}
  Vector<ModelParm> &getModelParms() {return parms;}
  AstForest *getSourceCode();
  void recompile();
  void deleteForeignObjects();
  void deleteClosures();
  Lambda::Closure *getCtor(int i);
  TransducerTemplate *instantiate(Lambda::Closure *ctor);
  virtual void pushAccessibles(MarkStack &);

protected:
  static ModelCompiler *global;
  AstForest *sourceCode;
  LambdaAPI &lambda;
  StateTemplate *q0;
  MacrostateTemplate *mq0;
  FunctionalClass fc0;
  Lambda::Closure *prob1; // function that evaluates to 1
  Vector<ModelParm> parms;
  static Vector<TransducerTemplate*> templates; // ### DEPRECATED
  Vector<Lambda::Closure*> ctors;
  Vector<Lambda::ForeignObject*> foreignObjects;
  Vector<Lambda::Closure*> closures;
  AstCombination *caller;
  AstObject *closureWrapper;
  double prevTreeScale;
  void registerFO(Lambda::ForeignObject*fo);
  void defineGlobals(bool redefine=false);
  void initCaller();
  static String declareSymbol(LambdaObject *,const String &funcName,
			      RunTimeEnvironment &);
  TransducerTemplate *compose(Vector<MacrostateTemplate*> &,
			      Vector<MacrotransTemplate*> &,
			      Vector<FO_GainLossFactor*> &gainLossFactors);
  void compose(TransducerTemplate &fromSub,TransducerTemplate &toSub,
	       Lambda::Closure &macroTrans,TransducerTemplate &,
	       bool allowFromQ0,bool allowToQ0);
  Lambda::Closure *tripleProduct(Lambda::Closure *,Lambda::Closure *,
				 Lambda::Closure *);
  Lambda::Closure *AoverBplusC(Lambda::Closure *,Lambda::Closure *,
				 Lambda::Closure *);
  Lambda::Closure *multiplyClosures(Lambda::Closure *,Lambda::Closure *);
  Lambda::Closure *sumClosures(Lambda::Closure *,Lambda::Closure *);
  Lambda::Closure *oneOverOneMinus(Lambda::Closure *x);
  Lambda::Closure *oneOverSum(Lambda::Closure *,Lambda::Closure *);
  TransducerTemplate *eliminateSilentStates(TransducerTemplate *in,
					    Array1D<STATE> &stateMap);
  TransducerTemplate *composeLossModel(TransducerTemplate *ttBackground,
				       Array1D<FunctionalClass> &foreground);
  TransducerTemplate *composeGainModel(TransducerTemplate *ttBackground,
				       Array1D<FunctionalClass> &foreground);
  TransducerTemplate *composeSimpleLossModel(TransducerTemplate *ttBackground,
				       Array1D<FunctionalClass> &foreground);
  TransducerTemplate *composeSimpleGainModel(TransducerTemplate *ttBackground,
				       Array1D<FunctionalClass> &foreground);
  TransducerTemplate *compose2tierLoss(TransducerTemplate *ttBackground,
				       Array1D<FunctionalClass> &foreground);
  TransducerTemplate *compose2tierGain(TransducerTemplate *ttBackground,
				       Array1D<FunctionalClass> &foreground);
  void decomposeBackground(TransducerTemplate *bg,
			   StateTemplate *&insertState,
			   StateTemplate *&deleteState,
			   StateTemplate *&matchState,
			   Lambda::Closure *&insertNorm,
			   Lambda::Closure *&deleteNorm,
			   Lambda::Closure *&matchNorm,
			   Array2D<Lambda::Closure*> &transMatrix,
			   FunctionalClass &);
  void normalize(PHMM_StateType from,PHMM_StateType to,Lambda::Closure *
		 normTerm,Array2D<Lambda::Closure*> &transMatrix);
  static LambdaObject *functionalClass(RunTimeEnvironment &,
				       ForegroundBackground,String func);

  // Callback functions (lambda "foreign functions"):
  static LambdaObject *ff_LoadRateMatrix(RunTimeEnvironment &,
					 ContinuationType &);
  static LambdaObject *ff_LoadRateMatrices(RunTimeEnvironment &,
					   ContinuationType &);
  static LambdaObject *ff_LoadPhylogeny(RunTimeEnvironment &,
					ContinuationType &);
  static LambdaObject *ff_RegisterTransducer(RunTimeEnvironment &,
					     ContinuationType &);
  static LambdaObject *ff_state(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_states(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_transition(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_transducer(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_macrostate(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_macrostates(RunTimeEnvironment &,ContinuationType&);
  static LambdaObject *ff_macrotrans(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_compose(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_functionalClass(RunTimeEnvironment &,
					  ContinuationType &);
  static LambdaObject *ff_funcElemType(RunTimeEnvironment &,
				       ContinuationType &);
  static LambdaObject *ff_background(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_lossOfFunction(RunTimeEnvironment &,
					 ContinuationType &);
  static LambdaObject *ff_gainOfFunction(RunTimeEnvironment &,
					 ContinuationType &);
  static LambdaObject *ff_simpleLossModel(RunTimeEnvironment &,
					 ContinuationType &);
  static LambdaObject *ff_simpleGainModel(RunTimeEnvironment &,
					 ContinuationType &);
  static LambdaObject *ff_2tierLossModel(RunTimeEnvironment &,
					 ContinuationType &);
  static LambdaObject *ff_2tierGainModel(RunTimeEnvironment &,
					 ContinuationType &);
  static LambdaObject *ff_bindingSite(RunTimeEnvironment &,ContinuationType&);
  static LambdaObject *ff_gainFactor(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_lossFactor(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_retentionFactor(RunTimeEnvironment &,
					  ContinuationType &);
  static LambdaObject *ff_scaleTree(RunTimeEnvironment &,ContinuationType &);
  static LambdaObject *ff_registerParm(RunTimeEnvironment &,
				       ContinuationType &);
  static LambdaObject *ff_registerConstructor(RunTimeEnvironment &,
					      ContinuationType &);
  static LambdaObject *ff_reverseStrand(RunTimeEnvironment &,
					ContinuationType &);
};

#endif

