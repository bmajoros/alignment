/****************************************************************
 ModelCompiler.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "TransducerTemplate.H"
#include "PhyLib/RateMatrix.H"
#include "ModelCompiler.H"
#include "FunctionalElement.H"
#include "BOOM/Random.H"
#include "BOOM/lambda/Parser.H"
#include "BOOM/lambda/AstPrinter.H"
#include "BOOM/lambda/LambdaException.H"
#include "BOOM/Array2D.H"
using namespace std;
using namespace BOOM;


ModelCompiler *ModelCompiler::global=NULL;
Vector<TransducerTemplate*> ModelCompiler::templates; // static


ostream &operator<<(ostream &os,ForeignObjectType t)
{
  switch(t) 
    {
    case FOT_RATE_MATRIX: os<<"FOT_RATE_MATRIX"; break;
    case FOT_PHYLOGENY:   os<<"FOT_PHYLOGENY";   break;
    case FOT_STATES:      os<<"FOT_STATES";      break;
    case FOT_STATE:       os<<"FOT_STATE";       break;
    case FOT_TRANSITION:  os<<"FOT_TRANSITION";  break;
    case FOT_TRANSDUCER:  os<<"FOT_TRANSDUCER";  break;
    case FOT_STATE_TYPE:  os<<"FOT_STATE_TYPE";  break;
    default: throw "operator<<(ForeignObjectType)";
    }
  return os;
}



/****************************************************************
                     ModelCompiler methods
 ****************************************************************/

ModelCompiler::ModelCompiler()
  : //lambda(GetRandomSeed()), 
  q0(NULL), lambda(*new LambdaAPI(GetRandomSeed())), prevTreeScale(1.0)
{
  // ctor

  /*if(!global)*/ global=this;
  lambda.getGC().setSilence(false);
  AstForest *prob1src=lambda.parseSource("[t|1]");
  lambda.makeImmortal(prob1src); // ###
  prob1=dynamic_cast<Lambda::Closure*>(lambda.evaluate(prob1src));
  lambda.makeImmortal(prob1); // ###

  defineGlobals();
  initCaller();
  lambda.registerMarkHook(this);
}



ModelCompiler::~ModelCompiler()
{
  templates.clear();
  delete &lambda;
}



void ModelCompiler::deleteClosures()
{
  /*
  int n=closures.size();
  Lambda::Closure **p=&closures[0];
  for(int i=0 ; i<n ; ++i) {
    lambda.registerWithGC(*p);
    ++p;
  }
  */

  // They're already registered with the GC, so we just have to stop
  // pushing them on the MarkStack:
  closures.clear();
}



void ModelCompiler::pushAccessibles(MarkStack &s)
{
  int n=closures.size();
  Lambda::Closure **p=&closures[0];
  for(int i=0 ; i<n ; ++i) {
    s.push(*p);
    ++p;
  }
  n=ctors.size();
  p=&ctors[0];
  for(int i=0 ; i<n ; ++i) {
    s.push(*p);
    ++p;
  }
}



void ModelCompiler::initCaller()
{
  AstForest *args=new AstForest;
  caller=lambda.makeCombination(NULL,args,false);
  closureWrapper=new AstObject(NULL);
  caller->changeFunction(closureWrapper);
  
  //###
  //lambda.makeImmortal(args);
  delete args;
  lambda.makeImmortal(caller);
  lambda.makeImmortal(closureWrapper);
  //###
}



Lambda::Closure *ModelCompiler::getCtor(int i)
{
  return ctors[i];
}



TransducerTemplate *ModelCompiler::instantiate(Lambda::Closure *ctor)
{
  closureWrapper->changeObject(ctor);
  LambdaObject *result=lambda.evaluate(caller); // ### SEG FAULT ####
  FO_Transducer *fo=dynamic_cast<FO_Transducer*>(result);
  if(!fo) throw LambdaException("contructor does not evaluate to transducer");
  TransducerTemplate *t=fo->getTemplate();
  //lambda.registerWithGC(fo);
  return t;
}



void ModelCompiler::deleteForeignObjects()
{
  Vector<Lambda::ForeignObject*>::iterator cur=foreignObjects.begin(),
    end=foreignObjects.end();
  for(; cur!=end ; ++cur) lambda.registerWithGC(*cur);
  foreignObjects.clear();
}



ModelCompiler *ModelCompiler::getGlobal()
{
  return global;
}



String ModelCompiler::declareSymbol(LambdaObject *obj,const String &funcName,
				    RunTimeEnvironment &env)
{
  if(!obj || obj->getType()!=OBTYPE_SYMBOL)
    throw LambdaException(String("function \"")+funcName+
			  "\" requires symbol");
  Lambda::Symbol *symbol=static_cast<Lambda::Symbol*>(obj);
  env.declareGlobalAndBackpatch(symbol->getLexeme());
  return symbol->getLexeme();
}



void ModelCompiler::defineGlobals(bool redefine)
{
  if(!redefine) {
    lambda.declareGlobal("MATCH");  // a state type
    lambda.declareGlobal("INSERT"); // a state type
    lambda.declareGlobal("DELETE"); // a state type
    lambda.declareGlobal("SILENT"); // a state type
    lambda.declareGlobal("start-state"); // start state for transducers
    lambda.declareGlobal("end-state");   // end state for transducers
    lambda.declareGlobal("start"); // start state for metamodels
    lambda.declareGlobal("end");   // end state for metamodels
  }

  fc0=FunctionalClass::registerClass("start/stop",'*',
				     FunctionalElementType::NO_FUNC_ELEM,
				     FB_NEITHER);
  delete q0;
  q0=new StateTemplate(0,PHMM_START_STOP,fc0,fc0);
  static FO_State *fo_q0=new FO_State(q0);
  static TransducerTemplate *emptySubmodel=new TransducerTemplate;
  if(redefine) {
    delete fo_q0;
    delete emptySubmodel;
    fo_q0=new FO_State(q0);
    emptySubmodel=new TransducerTemplate;
  }
  emptySubmodel->addState(q0);
  emptySubmodel->addTransition(new TransitionTemplate(q0,q0,prob1));
  mq0=new MacrostateTemplate(0,emptySubmodel);
  static FO_Macrostate *fo_mq0=new FO_Macrostate(mq0);
  lambda.assignToGlobal("MATCH",new FO_StateType(PHMM_MATCH));
  lambda.assignToGlobal("INSERT",new FO_StateType(PHMM_INSERT));
  lambda.assignToGlobal("DELETE",new FO_StateType(PHMM_DELETE));
  lambda.assignToGlobal("SILENT",new FO_StateType(PHMM_SILENT));
  lambda.assignToGlobal("start-state",fo_q0);
  lambda.assignToGlobal("end-state",fo_q0);
  lambda.assignToGlobal("start",fo_mq0);
  lambda.assignToGlobal("end",fo_mq0);

  lambda.defineFunction(&ff_LoadRateMatrix,"load-rate-matrix",2);
  lambda.defineFunction(&ff_LoadRateMatrices,"load-rate-matrices",2);
  lambda.defineFunction(&ff_LoadPhylogeny,"load-phylogeny",1);
  lambda.defineFunction(&ff_RegisterTransducer,"register-transducer",1);
  lambda.defineFunction(&ff_states,"states",VARIADIC);
  lambda.defineFunction(&ff_state,"state",VARIADIC); // 3 or 4 parms
  lambda.defineFunction(&ff_transition,"transition",3);
  lambda.defineFunction(&ff_transducer,"transducer",VARIADIC);
  lambda.defineFunction(&ff_macrostates,"macrostates",VARIADIC);
  lambda.defineFunction(&ff_macrostate,"macrostate",2);
  lambda.defineFunction(&ff_macrotrans,"macrotrans",3);
  lambda.defineFunction(&ff_compose,"compose",VARIADIC);
  lambda.defineFunction(&ff_functionalClass,"functional-class",5);
  lambda.defineFunction(&ff_funcElemType,"functional-element-type",4);
  lambda.defineFunction(&ff_background,"background",4);
  lambda.defineFunction(&ff_lossOfFunction,"full-loss-model",VARIADIC);
  lambda.defineFunction(&ff_gainOfFunction,"full-gain-model",VARIADIC);
  lambda.defineFunction(&ff_simpleLossModel,"loss-of-function",VARIADIC);
  lambda.defineFunction(&ff_simpleGainModel,"gain-of-function",VARIADIC);
  lambda.defineFunction(&ff_2tierLossModel,"two-tier-loss",VARIADIC);
  lambda.defineFunction(&ff_2tierGainModel,"two-tier-gain",VARIADIC);
  lambda.defineFunction(&ff_bindingSite,"binding-site",VARIADIC);
  lambda.defineFunction(&ff_gainFactor,"gain-factor",1);
  lambda.defineFunction(&ff_lossFactor,"loss-factor",1);
  lambda.defineFunction(&ff_retentionFactor,"retention-factor",1);
  lambda.defineFunction(&ff_scaleTree,"scale-tree",2);
  lambda.defineFunction(&ff_registerParm,"register-parm",3);
  lambda.defineFunction(&ff_registerConstructor,"register-constructor",1);
  lambda.defineFunction(&ff_reverseStrand,"reverse-strand",1);
}



void ModelCompiler::parse(const String &filename)
{
  AstForest *forest=lambda.parseFile(filename);
  lambda.makeImmortal(forest);
  sourceCode=forest;
  LambdaObject *result=lambda.evaluate(forest);
  //if(result) lambda.makeImmortal(result);
}



AstForest *ModelCompiler::getSourceCode()
{
  return sourceCode;
}



void ModelCompiler::recompile()
{
  parms.clear();
  templates.clear();
  defineGlobals(true);

  LambdaObject *result=lambda.evaluate(sourceCode); // returns NULL
  //lambda.makeImmortal(result);
}



int ModelCompiler::getNumModels() const
{
  return templates.size();
}



TransducerTemplate *ModelCompiler::getIthModel(int i)
{
  //if(i>=templates.size()) INTERNAL_ERROR; // ### DEBUGGING
  return i<templates.size() ? templates[i] : NULL;
}



TemplateInstantiator *ModelCompiler::getInstantiator()
{
  return new TemplateInstantiator(lambda);
}



Phylogeny *ModelCompiler::getGuideTree()
{
  Vector<TransducerTemplate*>::iterator cur=templates.begin(),
    end=templates.end();
  for(; cur!=end ; ++cur) {
    TransducerTemplate *TT=*cur;
    Phylogeny *phy=getGuideTree(TT);
    if(phy) return phy;
  }
  throw "ModelCompiler::getGuideTree";
}

Phylogeny *ModelCompiler::getGuideTree(TransducerTemplate *TT)
{
  int n=TT->getNumStates();
  for(int i=0 ; i<n ; ++i) {
    StateTemplate &st=TT->getIthState(i);
    Phylogeny *phy=st.getPhylogeny(PARENT);
    if(phy) return phy;
    phy=st.getPhylogeny(CHILD);
    if(phy) return phy;
  }
  return NULL;
}




/****************************************************************
                   lambda "foreign functions"
 ****************************************************************/

LambdaObject *ModelCompiler::ff_registerParm(RunTimeEnvironment &env,
					     ContinuationType &)
{
  // usage: (register-parm 'name <min-value> <max-value>)

  LambdaObject *parm0=env.getParm(0);
  if(parm0->getType()!=OBTYPE_SYMBOL)
    throw LambdaException("register-parm requires symbol parameter");
  const String &symbol=static_cast<Lambda::Symbol*>(parm0)->getLexeme();

  LambdaObject *parm1=env.getParm(1);
  if(parm1->getType()!=OBTYPE_FLOAT)
    throw LambdaException("register-parm requires two float parameters");
  float minVal=static_cast<LambdaFloat*>(parm1)->getValue();

  LambdaObject *parm2=env.getParm(2);
  if(parm2->getType()!=OBTYPE_FLOAT)
    throw LambdaException("register-parm requires two float parameters");
  float maxVal=static_cast<LambdaFloat*>(parm2)->getValue();

  global->parms.push_back(ModelParm(symbol,minVal,maxVal));
}



LambdaObject *ModelCompiler::ff_registerConstructor(RunTimeEnvironment &env,
						    ContinuationType &)
{
  // usage: (register-constructor closure)

  LambdaObject *parm0=env.getParm(0);
  if(parm0->getType()!=OBTYPE_CLOSURE)
    throw LambdaException("register-parm requires closure parameter");
  Lambda::Closure *ctor=static_cast<Lambda::Closure*>(parm0);
  global->closures.push_back(ctor);//env.makeImmortal(ctor);
  global->ctors.push_back(ctor);
}



LambdaObject *ModelCompiler::ff_LoadRateMatrix(RunTimeEnvironment &env,
					       ContinuationType &) 
{
  LambdaObject *parm=env.getParm(0);
  if(parm->getType()!=OBTYPE_STRING)
    throw LambdaException("load-rate-matrix requires string parameter");
  const String &filename=static_cast<LambdaString*>(parm)->getValue();

  LambdaFloat *parm1=dynamic_cast<LambdaFloat*>(env.getParm(1));
  if(!parm1)
    throw LambdaException("load-rate-matrix requires scale parameter");
  float scale=parm1->getValue();

  RateMatrix *matrix=RateMatrix::load(filename);
  matrix->rescale(scale);
  FO_RateMatrix *fo=new FO_RateMatrix(matrix);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_LoadRateMatrices(RunTimeEnvironment &env,
						 ContinuationType &) 
{
  LambdaObject *parm=env.getParm(0);
  if(parm->getType()!=OBTYPE_STRING)
    throw LambdaException("load-rate-matrices requires filename parameter");
  const String &filename=static_cast<LambdaString*>(parm)->getValue();

  LambdaFloat *parm1=dynamic_cast<LambdaFloat*>(env.getParm(1));
  if(!parm1)
    throw LambdaException("load-rate-matrices requires scale parameter");
  float scale=parm1->getValue();

  ifstream is(filename.c_str());
  if(!is.good()) throw LambdaException(String("can't open file ")+filename);
  int numMatrices;
  is>>numMatrices;
  FO_RateMatrices *matrices=new FO_RateMatrices(numMatrices);
  global->registerFO(matrices);
  for(int i=0 ; i<numMatrices ; ++i) {
    RateMatrix *matrix=RateMatrix::load(is);
    matrix->rescale(scale);
    FO_RateMatrix *fo=new FO_RateMatrix(matrix);
    global->registerFO(fo);
    matrices->setMatrix(i,fo);
  }
  return matrices;
}



LambdaObject *ModelCompiler::ff_LoadPhylogeny(RunTimeEnvironment &env,
					       ContinuationType &) 
{
  LambdaObject *parm=env.getParm(0);
  if(parm->getType()!=OBTYPE_STRING)
    throw LambdaException("load-phylogeny requires string parameter");
  const String &filename=static_cast<LambdaString*>(parm)->getValue();
  Phylogeny *phy=new Phylogeny(filename);
  FO_Phylogeny *fo=new FO_Phylogeny(phy);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_RegisterTransducer(RunTimeEnvironment &env,
						   ContinuationType &) 
{
  LambdaObject *parm=env.getParm(0);
  FO_Transducer *fo=dynamic_cast<FO_Transducer*>(parm);
  if(!fo) throw LambdaException("register-transducer : bad argument");
  // ########  global->registerFO(fo);//env.makeImmortal(fo); // ###
  ModelCompiler::templates.push_back(fo->getTemplate());
  return NULL;
}



LambdaObject *ModelCompiler::ff_state(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (state num type parentFuncClass childFuncClass)
  //     or: (state num type funcClass)
  //     or: (state num type) for SILENT states only

  // Process parms
  int numParms=env.getNumParms();
  if(numParms<2 || numParms>4) 
    throw LambdaException("(state) takes 2, 3, or 4 parameters");
  LambdaObject *parm0=env.getParm(0), *parm1=env.getParm(1);

  // Get the state ID
  LambdaInt *parm0int=dynamic_cast<LambdaInt*>(parm0);
  if(!parm0int) throw LambdaException("invalid state ID");
  int stateID=parm0int->getValue();
  if(stateID<1) throw LambdaException("state ID must be >0");

  // Get the state type
  FO_StateType *fo_stateType=dynamic_cast<FO_StateType*>(parm1);
  if(!fo_stateType) throw LambdaException("invalid state type");
  PairHMMStateType stateType=fo_stateType->getStateType();
  
  // Get the functional classes
  FunctionalClass parentFC, childFC;
  if(numParms==2)
    parentFC=childFC=FunctionalClass::NO_CLASS;
  else {
    FO_FuncClass *fo_parentFC=dynamic_cast<FO_FuncClass*>(env.getParm(2));
    if(!fo_parentFC) throw LambdaException("invalid parent functional class");
    parentFC=fo_parentFC->getClass();
    if(numParms==3) childFC=parentFC;
    else {
      FO_FuncClass *fo_childFC=dynamic_cast<FO_FuncClass*>(env.getParm(3));
      if(!fo_childFC) throw LambdaException("invalid child functional class");
      childFC=fo_childFC->getClass();
    }
  }

  // Allocate a state template and return it
  StateTemplate *stateTemplate=
     new StateTemplate(stateID,stateType,parentFC,childFC);
  FO_State *fo_State=new FO_State(stateTemplate);
  global->registerFO(fo_State);
  return fo_State;
}



LambdaObject *ModelCompiler::ff_states(RunTimeEnvironment &env,
				       ContinuationType &) 
{
  int numParms=env.getNumParms();
  Array1D<FO_State*> &states=*new Array1D<FO_State*>(numParms);
  for(int i=0 ; i<numParms ; ++i) states[i]=env.getParm(i);
  FO_States *fo=new FO_States(states);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_transition(RunTimeEnvironment &env,
					   ContinuationType &) 
{
  LambdaObject *parm0=env.getParm(0), *parm1=env.getParm(1),
    *parm2=env.getParm(2);
  FO_State *from=dynamic_cast<FO_State*>(parm0);
  FO_State *to=dynamic_cast<FO_State*>(parm1);
  if(!from || !to) 
    throw LambdaException("invalid state in transition construct");
  Lambda::Closure *func=dynamic_cast<Lambda::Closure*>(parm2);
  global->closures.push_back(func);//env.makeImmortal(func);//###
  if(!func) throw LambdaException("invalid transition function");
  TransitionTemplate *transTemplate=
    new TransitionTemplate(from->getState(),to->getState(),func);
  FO_Transition *fo=new FO_Transition(transTemplate);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_transducer(RunTimeEnvironment &env,
					   ContinuationType &)
{
  TransducerTemplate *t=new TransducerTemplate;
  t->addState(getGlobal()->q0);
  int numParms=env.getNumParms();
  for(int i=0 ; i<numParms ; ++i) {
    LambdaObject *parm=env.getParm(i);
    if(!parm) continue;
    ForeignObject *obj=dynamic_cast<ForeignObject*>(parm);
    switch(obj->getTag()) 
      {
      case FOT_STATES: 
	{
	  FO_States *states=static_cast<FO_States*>(obj);
	  Array1D<FO_State*> &statesArray=states->getStates();
	  int numStates=statesArray.size();
	  for(int i=0 ; i<numStates ; ++i) 
	    t->addState(statesArray[i]->getState());
	}
	break;
      case FOT_TRANSITION:
	t->addTransition(static_cast<FO_Transition*>(obj)->getTransition());
	break;
      }
  }
  Array1D<STATE> stateMap;
  t=global->eliminateSilentStates(t,stateMap);
  FO_Transducer *fo=new FO_Transducer(t);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_macrostate(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (macrostate num submodel)

  // Process parms
  LambdaObject *parm0=env.getParm(0), *parm1=env.getParm(1);

  // Get the state ID
  LambdaInt *parm0int=dynamic_cast<LambdaInt*>(parm0);
  if(!parm0int) throw LambdaException("invalid state ID in (macrostate)");
  int stateID=parm0int->getValue();
  if(stateID<1) throw LambdaException("state ID must be >0 in (macrostate)");

  // Get the submodel
  FO_Transducer *fo_submodel=dynamic_cast<FO_Transducer*>(parm1);
  if(!fo_submodel) throw LambdaException("invalid submodel");

  // Allocate a macrostate template and store it in the global namespace
  MacrostateTemplate *stateTemplate=
    new MacrostateTemplate(stateID,fo_submodel->getTemplate());
  FO_Macrostate *fo_macrostate=new FO_Macrostate(stateTemplate);
  global->registerFO(fo_macrostate);
  return fo_macrostate;
}



LambdaObject *ModelCompiler::ff_macrostates(RunTimeEnvironment &env,
				       ContinuationType &) 
{
  int numParms=env.getNumParms();
  Array1D<FO_Macrostate*> &states=*new Array1D<FO_Macrostate*>(numParms);
  for(int i=0 ; i<numParms ; ++i) states[i]=env.getParm(i);
  FO_Macrostates *fo=new FO_Macrostates(states);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_macrotrans(RunTimeEnvironment &env,
					   ContinuationType &) 
{
  LambdaObject *parm0=env.getParm(0), *parm1=env.getParm(1),
    *parm2=env.getParm(2);
  FO_Macrostate *from=dynamic_cast<FO_Macrostate*>(parm0);
  FO_Macrostate *to=dynamic_cast<FO_Macrostate*>(parm1);
  if(!from || !to) 
    throw LambdaException("invalid macrostate in macrotrans construct");
  Lambda::Closure *func=dynamic_cast<Lambda::Closure*>(parm2);
  global->closures.push_back(func);//env.makeImmortal(func);
  if(!func) throw LambdaException("invalid transition function in macrotrans");
  MacrotransTemplate *transTemplate=
    new MacrotransTemplate(from->getState(),to->getState(),func);
  FO_Macrotrans *fo=new FO_Macrotrans(transTemplate);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_functionalClass(RunTimeEnvironment &env,
						ContinuationType &) 
{
  // usage: (functional-class matrix phylogeny label name elementType)

  return functionalClass(env,FOREGROUND,"functional-class");
}


LambdaObject *ModelCompiler::ff_background(RunTimeEnvironment &env,
					ContinuationType &) 
{
  // usage: (background matrix phylogeny label name)

  return functionalClass(env,BACKGROUND,"background");
}



LambdaObject *ModelCompiler::functionalClass(RunTimeEnvironment &env,
					     ForegroundBackground fb,
					     String func)
{
  LambdaObject *parm0=env.getParm(0), *parm1=env.getParm(1),
    *parm2=env.getParm(2), *parm3=env.getParm(3);
  FO_RateMatrix *matrix=dynamic_cast<FO_RateMatrix*>(parm0);
  if(!matrix)
    throw LambdaException(func+" requires matrix parameter");
  FO_Phylogeny *phylogeny=dynamic_cast<FO_Phylogeny*>(parm1);
  if(!phylogeny) 
    throw LambdaException(func+" requires phylogeny parameter");
  
  // Get the state label
  LambdaString *labelStr=dynamic_cast<LambdaString*>(parm2);
  if(!labelStr)
    throw LambdaException(func+" requires label parameter");
  const String &labStr=labelStr->getValue();
  if(labStr.length()!=1) throw LambdaException("label must have length 1");
  char label=labStr[0];

  // Get the state name
  LambdaString *nameStr=dynamic_cast<LambdaString*>(parm3);
  if(!nameStr)
    throw LambdaException(func+" requires name parameter");
  const String &name=nameStr->getValue();

  // Get the element-type
  FunctionalElementType type;
  if(fb==BACKGROUND) type=FunctionalElementType::NO_FUNC_ELEM;
  else {
    LambdaString *typeStr=dynamic_cast<LambdaString*>(env.getParm(4));
    if(!typeStr)
      throw LambdaException(func+" requires element-type parameter");
    const String &typeName=typeStr->getValue();
    type=FunctionalElementType::getTypeByName(typeName);
    if(!type.isValid()) type=FunctionalElementType::registerType(typeName);
  }
  
  // Allocate a functional class and store it in the global namespace
  RateMatrix *Q=matrix->getMatrix();
  Phylogeny *tree=phylogeny->getPhylogeny();
  FunctionalClass fc=
    FunctionalClass::registerClass(name,label,type,fb,Q,tree);
  FO_FuncClass *fo_funcClass=new FO_FuncClass(Q,tree,label,name,fc);
  global->registerFO(fo_funcClass);
  return fo_funcClass;
}



LambdaObject *ModelCompiler::ff_funcElemType(RunTimeEnvironment &env,
					     ContinuationType &) 
{
  static String func="functional-element-type";
  LambdaObject *parm0=env.getParm(0), *parm1=env.getParm(1),
    *parm2=env.getParm(2), *parm3=env.getParm(3);

  // Get the element-type
  LambdaString *typeStr=dynamic_cast<LambdaString*>(parm0);
  if(!typeStr)
    throw LambdaException(func+" requires element-type parameter");
  const String &typeName=typeStr->getValue();
  FunctionalElementType type=FunctionalElementType::getTypeByName(typeName);
  if(!type.isValid()) type=FunctionalElementType::registerType(typeName);

  // Get the single-letter label
  LambdaString *labelStr=dynamic_cast<LambdaString*>(parm1);
  if(!labelStr) throw LambdaException(func+" requires label parameter");
  const String &labStr=labelStr->getValue();
  if(labStr.length()!=1) throw LambdaException("label must have length 1");
  char label=labStr[0];

  FO_RateMatrices *matrices=dynamic_cast<FO_RateMatrices*>(parm2);
  if(!matrices) throw LambdaException(func+" requires matrix parameter");
  int numClasses=matrices->getNumMatrices();

  FO_Phylogeny *phylogeny=dynamic_cast<FO_Phylogeny*>(parm3);
  if(!phylogeny) throw LambdaException(func+" requires phylogeny parameter");
  Phylogeny *tree=phylogeny->getPhylogeny();
  
  // Allocate functional classes
  FO_FuncElemType *fo_funcElemType=new FO_FuncElemType(numClasses);
  for(int i=0 ; i<numClasses ; ++i) {
    RateMatrix *Q=matrices->getMatrix(i)->getMatrix();
    String name=labStr+"_"+i;
    FunctionalClass fc=
      FunctionalClass::registerClass(name,label,type,FOREGROUND,Q,tree);
    type.appendClass(fc);
    FO_FuncClass *fo_funcClass=new FO_FuncClass(Q,tree,label,name,fc);
    global->registerFO(fo_funcClass);
    fo_funcElemType->setClass(i,fo_funcClass);
  }
  global->registerFO(fo_funcElemType);
  return fo_funcElemType;
}



LambdaObject *ModelCompiler::ff_reverseStrand(RunTimeEnvironment &env,
					     ContinuationType &) 
{
  static String func="reverse-strand";
  LambdaObject *parm0=env.getParm(0);
  FO_FuncElemType *type=dynamic_cast<FO_FuncElemType*>(parm0);
  if(!type)
    throw LambdaException(func+" requires element-type parameter");
  //FunctionalElementType &fet=type->getIthClass(0)->getClass().getElementType();
  int n=type->getNumClasses();
  for(int i=0 ; i<n ; ++i) {
    FO_FuncClass *fo=type->getIthClass(i);
    FunctionalClass c=fo->getClass();
    c.setStrand(complement(c.getStrand()));
  }
  return parm0;
}



LambdaObject *ModelCompiler::ff_compose(RunTimeEnvironment &env,
					   ContinuationType &)
{
  Vector<MacrostateTemplate*> macrostates;
  Vector<MacrotransTemplate*> macrotrans;
  Vector<FO_GainLossFactor*> gainLossFactors;
  macrostates.push_back(getGlobal()->mq0);
  int numParms=env.getNumParms();
  for(int i=0 ; i<numParms ; ++i) {
    ForeignObject *obj=dynamic_cast<ForeignObject*>(env.getParm(i));
    switch(obj->getTag()) 
      {
      case FOT_MACROSTATES: 
	{
	  FO_Macrostates *states=static_cast<FO_Macrostates*>(obj);
	  Array1D<FO_Macrostate*> &statesArray=states->getStates();
	  int numStates=statesArray.size();
	  for(int i=0 ; i<numStates ; ++i) {
	    MacrostateTemplate *oldState=statesArray[i]->getState();
	    macrostates.push_back(oldState);
	  }
	}
	break;
      case FOT_MACROTRANS:
	macrotrans.push_back(static_cast<FO_Macrotrans*>(obj)->
			     getTransition()); 
	break;
      case FOT_GAINLOSS_FACTOR:
	gainLossFactors.push_back(static_cast<FO_GainLossFactor*>(obj));
	break;
      }
  }
  TransducerTemplate *tt=
    getGlobal()->compose(macrostates,macrotrans,gainLossFactors);
  FO_Transducer *fo=new FO_Transducer(tt);
  global->registerFO(fo);
  return fo;
}



TransducerTemplate *ModelCompiler::compose(Vector<MacrostateTemplate*> 
					   &macrostates,
					   Vector<MacrotransTemplate*> 
					   &macrotransitions,
					   Vector<FO_GainLossFactor*> 
					   &gainLossFactors)
{
  TransducerTemplate *t=new TransducerTemplate;
  t->addState(q0);

  // Process states and within-submodel transitions
  int n=macrostates.size();
  for(int i=0 ; i<n ; ++i) {
    MacrostateTemplate *mst=macrostates[i];
    TransducerTemplate *submodel=mst->getSubmodel();
    int m=submodel->getNumStates();
    for(int j=0 ; j<m ; ++j) {
      StateTemplate &st=submodel->getIthState(j);
      if(st.getStateType()==PHMM_START_STOP) continue;
      t->addState(&st);
      st.changeStateID(t->getNumStates()-1);
    }
    m=submodel->getNumTransitions();
    for(int j=0 ; j<m ; ++j) {
      TransitionTemplate &tt=submodel->getIthTransition(j);
      bool fromStart=(tt.getFrom()->getStateType()==PHMM_START_STOP);
      bool toStop=(tt.getTo()->getStateType()==PHMM_START_STOP);
      if(!fromStart && !toStop) t->addTransition(&tt);
    }
  }

  // Process between-submodel transitions
  n=macrotransitions.size();
  for(int i=0 ; i<n ; ++i) {
    MacrotransTemplate &mtt=*macrotransitions[i];
    MacrostateTemplate *fromMacro=mtt.getFrom(), *toMacro=mtt.getTo();
    TransducerTemplate *fromSub=fromMacro->getSubmodel();
    TransducerTemplate *toSub=toMacro->getSubmodel();
    Lambda::Closure *probFunc=mtt.getTransProbFunc();
    bool allowFromQ0=(fromMacro->getStateID()==0);
    bool allowToQ0=(toMacro->getStateID()==0);
    compose(*fromSub,*toSub,*probFunc,*t,allowFromQ0,allowToQ0);
  }

  // Process gain/loss/retention factors
  n=gainLossFactors.size();
  for(int i=0 ; i<n ; ++i) {
    FO_GainLossFactor *factor=gainLossFactors[i];
    t->setGainLossFactor(factor->getFactor(),factor->getType());
  }

  return t;
}



void ModelCompiler::compose(TransducerTemplate &fromSub,
			    TransducerTemplate &toSub,
			    Lambda::Closure &macroTrans,
			    TransducerTemplate &t,
			    bool allowFromQ0,bool allowToQ0)
{
  const int nFrom=fromSub.getNumTransitions(), nTo=toSub.getNumTransitions();
  for(int i=0 ; i<nFrom ; ++i) {
    TransitionTemplate &ttFrom=fromSub.getIthTransition(i);
    if(ttFrom.getTo()->getStateType()!=PHMM_START_STOP) continue;
    StateTemplate *fromState=ttFrom.getFrom();
    if(fromState==q0 && !allowFromQ0) continue;
    Lambda::Closure *a=ttFrom.getTransProbFunc();
    for(int j=0 ; j<nTo ; ++j) {
      TransitionTemplate &ttTo=toSub.getIthTransition(j);
      if(ttTo.getFrom()->getStateType()!=PHMM_START_STOP) continue;
      StateTemplate *toState=ttTo.getTo();
      if(toState==q0 && !allowToQ0) continue;
      Lambda::Closure *b=ttTo.getTransProbFunc();
      Lambda::Closure *composedProb=tripleProduct(a,b,&macroTrans);
      TransitionTemplate *newTrans=
	new TransitionTemplate(fromState,toState,composedProb);
      t.addTransition(newTrans);
    }
  }
}



Lambda::Closure *ModelCompiler::AoverBplusC(Lambda::Closure *a,
					    Lambda::Closure *b,
					    Lambda::Closure *c)
{
  if(!a || !b || !c) throw "ModelCompiler::AoverBplusC()";
  Lambda::AstForest *forest=
    lambda.parseSource("[t|(/ (exp t) (+ (exp t) (exp t)))]");
  Lambda::Closure *z=
    dynamic_cast<Lambda::Closure*>(lambda.getEvaluator().evaluate(forest));
  if(!z) throw "Error composing three closures in ModelCompiler";
  z->setStaticChain(lambda.getEnvironment().getGlobalAR());
  struct V : public Lambda::AstPostorderTraversal {
    Array1D<Lambda::Closure*> closures;
    int nextClosure;
    V(Lambda::Closure *a,Lambda::Closure *b,Lambda::Closure *c)
      : closures(3), nextClosure(0)
    { closures[0]=a; closures[1]=b; closures[2]=c; }
    void process(AstCombination &combo) {
      if(combo.getParms().size()==1)
	combo.changeFunction(new AstObject(closures[nextClosure++]));
    }
  } v(a,b,c);
  v.traverse(z->getFunction());
  global->closures.push_back(z);//lambda.makeImmortal(z);//###
  return z;
}



Lambda::Closure *ModelCompiler::tripleProduct(Lambda::Closure *a,
					      Lambda::Closure *b,
					      Lambda::Closure *c)
{
  if(!a || !b || !c) throw "ModelCompiler::tripleProduct()";
  Lambda::AstForest *forest=
    lambda.parseSource("[t|(* (* (exp t) (exp t)) (exp t))]");
  Lambda::Closure *z=
    dynamic_cast<Lambda::Closure*>(lambda.getEvaluator().evaluate(forest));
  if(!z) throw "Error composing three closures in ModelCompiler";
  z->setStaticChain(lambda.getEnvironment().getGlobalAR());
  struct V : public Lambda::AstPostorderTraversal {
    Array1D<Lambda::Closure*> closures;
    int nextClosure;
    V(Lambda::Closure *a,Lambda::Closure *b,Lambda::Closure *c)
      : closures(3), nextClosure(0)
    { closures[0]=a; closures[1]=b; closures[2]=c; }
    void process(AstCombination &combo) {
      if(combo.getParms().size()==1)
	combo.changeFunction(new AstObject(closures[nextClosure++]));
    }
  } v(a,b,c);
  v.traverse(z->getFunction());
  global->closures.push_back(z);//lambda.makeImmortal(z);
  return z;
}



Lambda::Closure *ModelCompiler::multiplyClosures(Lambda::Closure *a,
					      Lambda::Closure *b)
{
  if(!a || !b) throw "ModelCompiler::multiplyClosures()";
  Lambda::AstForest *forest=lambda.parseSource("[t|(* (exp t) (exp t))]");
  Lambda::Closure *z=
    dynamic_cast<Lambda::Closure*>(lambda.getEvaluator().evaluate(forest));
  if(!z) throw "Error composing two closures in ModelCompiler";
  z->setStaticChain(lambda.getEnvironment().getGlobalAR());
  struct V : public Lambda::AstPostorderTraversal {
    Array1D<Lambda::Closure*> closures;
    int nextClosure;
    V(Lambda::Closure *a,Lambda::Closure *b)
      : closures(2), nextClosure(0) { closures[0]=a; closures[1]=b; }
    void process(AstCombination &combo) {
      if(combo.getParms().size()==1)
	combo.changeFunction(new AstObject(closures[nextClosure++]));
    }
  } v(a,b);
  v.traverse(z->getFunction());
  global->closures.push_back(z);//lambda.makeImmortal(z);
  return z;
}



Lambda::Closure *ModelCompiler::sumClosures(Lambda::Closure *a,
					    Lambda::Closure *b)
{
  if(!a || !b) throw "ModelCompiler::sumClosures()";
  Lambda::AstForest *forest=lambda.parseSource("[t|(+ (exp t) (exp t))]");
  Lambda::Closure *z=
    dynamic_cast<Lambda::Closure*>(lambda.getEvaluator().evaluate(forest));
  if(!z) throw "Error adding two closures in ModelCompiler";
  z->setStaticChain(lambda.getEnvironment().getGlobalAR());
  struct V : public Lambda::AstPostorderTraversal {
    Array1D<Lambda::Closure*> closures;
    int nextClosure;
    V(Lambda::Closure *a,Lambda::Closure *b)
      : closures(2), nextClosure(0) { closures[0]=a; closures[1]=b; }
    void process(AstCombination &combo) {
      if(combo.getParms().size()==1)
	combo.changeFunction(new AstObject(closures[nextClosure++]));
    }
  } v(a,b);
  v.traverse(z->getFunction());
  global->closures.push_back(z);//lambda.makeImmortal(z);
  return z;
}



TransducerTemplate *
ModelCompiler::eliminateSilentStates(TransducerTemplate *oldTT,
					  Array1D<STATE> &stateMap)
{
  // ### We assume there are no loops among the silent states

  TransducerTemplate *newTT=new TransducerTemplate;

  // First, compute a matrix representation of the transitions
  int numStates=oldTT->getNumStates(), numTrans=oldTT->getNumTransitions();
  Array2D<Lambda::Closure*> transMatrix(numStates,numStates);// for oldTT
  transMatrix.setAllTo(NULL);
  for(int i=0 ; i<numTrans ; ++i) {
    TransitionTemplate *tt=&oldTT->getIthTransition(i);
    STATE fromID=tt->getFrom()->getStateID();
    STATE toID=tt->getTo()->getStateID();
    transMatrix[fromID][toID]=tt->getTransProbFunc();
  }

  // Compute stateMap and its inverse
  Vector<STATE> vStateMap;     // new -> old
  Map<STATE,STATE> inverseMap; // old -> new
  for(STATE i=0 ; i<numStates ; ++i) {
    StateTemplate *st=&oldTT->getIthState(i);
    if(st->getStateType()!=PHMM_SILENT) {
      STATE newStateID=vStateMap.size();
      vStateMap.push_back(i);
      inverseMap[i]=newStateID;
      newTT->addState(st);
      st->changeStateID(newStateID);
    }
  }
  stateMap=vStateMap;
  int newNumStates=stateMap.size();

  // Now eliminate silent state transitions using a recursive procedure:
  class Eliminator {
    ModelCompiler &compiler;
    const Array1D<STATE> &stateMap;
    const Map<STATE,STATE> &inverseMap;
    Array2D<Lambda::Closure*> &oldTransMatrix, newTransMatrix;
    const TransducerTemplate &oldTT;
    TransducerTemplate &newTT;
    int oldN, newN;
    void addNewTrans(STATE oldFrom,STATE oldTo,Lambda::Closure *P) 
    {
      STATE newFrom=inverseMap[oldFrom], newTo=inverseMap[oldTo];
      Lambda::Closure *p=newTransMatrix[newFrom][newTo];
      if(p) P=compiler.sumClosures(p,P);
      newTransMatrix[newFrom][newTo]=P;
    }
    void recurs(STATE newDest,STATE oldDest,STATE oldCur,Lambda::Closure *P,
		int depth=0)
    {
      if(oldTT.getIthState(oldCur).getStateType()==PHMM_SILENT)
	for(STATE q=0 ; q<oldN ; ++q) {
	  Lambda::Closure *p=oldTransMatrix[q][oldCur];
	  if(p) {
	    if(depth>2*oldN) throw "Infinite recursion detected in "
	      "ModelCompiler::eliminateSilentStates()";
	    recurs(newDest,oldDest,q,compiler.multiplyClosures(P,p),depth+1);
	  }
	}
      else addNewTrans(oldCur,oldDest,P);
    }	
  public:
    Eliminator(const TransducerTemplate &oldTT,TransducerTemplate &newTT,
	       Array2D<Lambda::Closure*> &transMatrix,ModelCompiler &cmp,
	       const Array1D<STATE> &stateMap,
	       const Map<STATE,STATE> &inverseMap)
      : oldTT(oldTT), newTT(newTT), stateMap(stateMap), compiler(cmp),
	inverseMap(inverseMap), oldN(oldTT.getNumStates()),
	newN(newTT.getNumStates()), oldTransMatrix(transMatrix)
    {
      newTransMatrix.resize(newN,newN);
      newTransMatrix.setAllTo(NULL);
      for(STATE newState=0 ; newState<newN ; ++newState) {
	STATE oldState=stateMap[newState];
	for(STATE oldPred=0 ; oldPred<oldN ; ++oldPred) {
	  Lambda::Closure *p=oldTransMatrix[oldPred][oldState];
	  if(p) recurs(newState,oldState,oldPred,p);
	}
      }
      for(int i=0 ; i<newN ; ++i) for(int j=0 ; j<newN ; ++j) {
	Lambda::Closure *p=newTransMatrix[i][j];
	if(p) {
	  StateTemplate &a=newTT.getIthState(i), &b=newTT.getIthState(j);
	  TransitionTemplate *tt=new TransitionTemplate(&a,&b,p);
	  newTT.addTransition(tt);
	}
      }
    }
  };
  Eliminator eliminator(*oldTT,*newTT,transMatrix,*this,stateMap,inverseMap);
  return newTT;
}



LambdaObject *ModelCompiler::ff_lossOfFunction(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (full-loss-model background-model func-elem-type)

  // Process parms
  int numParms=env.getNumParms();
  if(numParms!=2)
    throw LambdaException("'full-loss-model' requires two parameters");

  // Get the background model
  FO_Transducer *fo_bg=dynamic_cast<FO_Transducer*>(env.getParm(0));
  if(!fo_bg) throw LambdaException("full-loss-model : bad argument");
  TransducerTemplate *ttBackground=fo_bg->getTemplate();

  // Get the functional classes for the foreground model
  FO_FuncElemType *fo_funcElemType=
    dynamic_cast<FO_FuncElemType*>(env.getParm(1));
  if(!fo_funcElemType) 
    throw LambdaException("wrong parameter type in 'full-loss-model'");
  int numColumns=fo_funcElemType->getNumClasses();
  Array1D<FunctionalClass> foregroundClasses(numColumns);
  for(int i=0 ; i<numColumns ; ++i)
    foregroundClasses[i]=fo_funcElemType->getIthClass(i)->getClass();

  // Compose foreground and background into a loss model
  TransducerTemplate *lossModel=
    global->composeLossModel(ttBackground,foregroundClasses);
  
  // Package it inside a Lambda object and return it
  FO_Transducer *fo=new FO_Transducer(lossModel);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_gainOfFunction(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (full-gain-model background-model func-elem-type)

  // Process parms
  int numParms=env.getNumParms();
  if(numParms!=2)
    throw LambdaException("'full-gain-model' requires two parameters");

  // Get the background model
  FO_Transducer *fo_bg=dynamic_cast<FO_Transducer*>(env.getParm(0));
  if(!fo_bg) throw LambdaException("full-gain-model : bad argument");
  TransducerTemplate *ttBackground=fo_bg->getTemplate();

  // Get the functional classes for the foreground model
  FO_FuncElemType *fo_funcElemType=
    dynamic_cast<FO_FuncElemType*>(env.getParm(1));
  if(!fo_funcElemType) 
    throw LambdaException("wrong parameter type in 'full-gain-model'");
  int numColumns=fo_funcElemType->getNumClasses();
  Array1D<FunctionalClass> foregroundClasses(numColumns);
  for(int i=0 ; i<numColumns ; ++i)
    foregroundClasses[i]=fo_funcElemType->getIthClass(i)->getClass();
  // Compose foreground and background into a loss model
  TransducerTemplate *gainModel=
    global->composeGainModel(ttBackground,foregroundClasses);
  // Package it inside a Lambda object and return it
  FO_Transducer *fo=new FO_Transducer(gainModel);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_simpleLossModel(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (loss-of-function background-model func-elem-type)

  // Process parms
  int numParms=env.getNumParms();
  if(numParms!=2)
    throw LambdaException("'loss-of-function' requires two parameters");

  // Get the background model
  FO_Transducer *fo_bg=dynamic_cast<FO_Transducer*>(env.getParm(0));
  if(!fo_bg) throw LambdaException("loss-of-function : bad argument");
  TransducerTemplate *ttBackground=fo_bg->getTemplate();

  // Get the functional classes for the foreground model
  FO_FuncElemType *fo_funcElemType=
    dynamic_cast<FO_FuncElemType*>(env.getParm(1));
  if(!fo_funcElemType) 
    throw LambdaException("wrong parameter type in 'loss-of-function'");
  int numColumns=fo_funcElemType->getNumClasses();
  Array1D<FunctionalClass> foregroundClasses(numColumns);
  for(int i=0 ; i<numColumns ; ++i)
    foregroundClasses[i]=fo_funcElemType->getIthClass(i)->getClass();

  // Compose foreground and background into a loss model
  TransducerTemplate *lossModel=
    global->composeSimpleLossModel(ttBackground,foregroundClasses);
  
  // Package it inside a Lambda object and return it
  FO_Transducer *fo=new FO_Transducer(lossModel);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_simpleGainModel(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (gain-of-function background-model func-elem-type)

  // Process parms
  int numParms=env.getNumParms();
  if(numParms!=2)
    throw LambdaException("'gain-of-function' requires two parameters");

  // Get the background model
  FO_Transducer *fo_bg=dynamic_cast<FO_Transducer*>(env.getParm(0));
  if(!fo_bg) throw LambdaException("gain-of-function : bad argument");
  TransducerTemplate *ttBackground=fo_bg->getTemplate();

  // Get the functional classes for the foreground model
  FO_FuncElemType *fo_funcElemType=
    dynamic_cast<FO_FuncElemType*>(env.getParm(1));
  if(!fo_funcElemType) 
    throw LambdaException("wrong parameter type in 'gain-of-function'");
  int numColumns=fo_funcElemType->getNumClasses();
  Array1D<FunctionalClass> foregroundClasses(numColumns);
  for(int i=0 ; i<numColumns ; ++i)
    foregroundClasses[i]=fo_funcElemType->getIthClass(i)->getClass();

  // Compose foreground and background into a loss model
  TransducerTemplate *gainModel=
    global->composeSimpleGainModel(ttBackground,foregroundClasses);
  
  // Package it inside a Lambda object and return it
  FO_Transducer *fo=new FO_Transducer(gainModel);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_2tierLossModel(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (two-tier-loss background-model func-elem-type)

  // Process parms
  int numParms=env.getNumParms();
  if(numParms!=2)
    throw LambdaException("'two-tier-loss' requires two parameters");

  // Get the background model
  FO_Transducer *fo_bg=dynamic_cast<FO_Transducer*>(env.getParm(0));
  if(!fo_bg) throw LambdaException("two-tier-loss : bad argument");
  TransducerTemplate *ttBackground=fo_bg->getTemplate();

  // Get the functional classes for the foreground model
  FO_FuncElemType *fo_funcElemType=
    dynamic_cast<FO_FuncElemType*>(env.getParm(1));
  if(!fo_funcElemType) 
    throw LambdaException("wrong parameter type in 'two-tier-loss'");
  int numColumns=fo_funcElemType->getNumClasses();
  Array1D<FunctionalClass> foregroundClasses(numColumns);
  for(int i=0 ; i<numColumns ; ++i)
    foregroundClasses[i]=fo_funcElemType->getIthClass(i)->getClass();

  // Compose foreground and background into a loss model
  TransducerTemplate *lossModel=
    global->compose2tierLoss(ttBackground,foregroundClasses);
  
  // Package it inside a Lambda object and return it
  FO_Transducer *fo=new FO_Transducer(lossModel);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_2tierGainModel(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (two-tier-gain background-model func-elem-type)

  // Process parms
  int numParms=env.getNumParms();
  if(numParms!=2)
    throw LambdaException("'two-tier-gain' requires two parameters");

  // Get the background model
  FO_Transducer *fo_bg=dynamic_cast<FO_Transducer*>(env.getParm(0));
  if(!fo_bg) throw LambdaException("two-tier-gain : bad argument");
  TransducerTemplate *ttBackground=fo_bg->getTemplate();

  // Get the functional classes for the foreground model
  FO_FuncElemType *fo_funcElemType=
    dynamic_cast<FO_FuncElemType*>(env.getParm(1));
  if(!fo_funcElemType) 
    throw LambdaException("wrong parameter type in 'two-tier-gain'");
  int numColumns=fo_funcElemType->getNumClasses();
  Array1D<FunctionalClass> foregroundClasses(numColumns);
  for(int i=0 ; i<numColumns ; ++i)
    foregroundClasses[i]=fo_funcElemType->getIthClass(i)->getClass();

  // Compose foreground and background into a loss model
  TransducerTemplate *gainModel=
    global->compose2tierGain(ttBackground,foregroundClasses);
  
  // Package it inside a Lambda object and return it
  FO_Transducer *fo=new FO_Transducer(gainModel);
  global->registerFO(fo);
  return fo;
}



/*
LambdaObject *ModelCompiler::ff_lossOfFunction(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (loss-of-function background-model fg1 fg2 fg3 ... fgN)

  // Process parms
  int numParms=env.getNumParms();
  if(numParms<2)
    throw LambdaException("(loss-of-function) requires >1 parameters");
  LambdaObject *parm0=env.getParm(0);

  // Get the background model
  FO_Transducer *fo_bg=dynamic_cast<FO_Transducer*>(env.getParm(0));
  if(!fo_bg) throw LambdaException("loss-of-function : bad argument");
  TransducerTemplate *ttBackground=fo_bg->getTemplate();

  // Get the functional classes for the foreground model
  int numColumns=numParms-1;
  Array1D<FunctionalClass> foregroundClasses(numColumns);
  for(int i=0 ; i<numColumns ; ++i) {
    FO_FuncClass *fo_fg=dynamic_cast<FO_FuncClass*>(env.getParm(i+1));
    if(!fo_fg) throw 
      LambdaException("invalid functional class in call to loss-of-function");
    foregroundClasses[i]=fo_fg->getClass();
  }

  // Compose foreground and background into a loss model
  TransducerTemplate *lossModel=
    global->composeLossModel(ttBackground,foregroundClasses);
  
  // Package it inside a Lambda object and return it
  return new FO_Transducer(lossModel);
}



LambdaObject *ModelCompiler::ff_gainOfFunction(RunTimeEnvironment &env,
				      ContinuationType &) 
{
  // syntax: (gain-of-function background-model fg1 fg2 fg3 ... fgN)

  // Process parms
  int numParms=env.getNumParms();
  if(numParms<2)
    throw LambdaException("(gain-of-function) requires >1 parameters");
  LambdaObject *parm0=env.getParm(0);

  // Get the background model
  FO_Transducer *fo_bg=dynamic_cast<FO_Transducer*>(env.getParm(0));
  if(!fo_bg) throw LambdaException("gain-of-function : bad argument");
  TransducerTemplate *ttBackground=fo_bg->getTemplate();

  // Get the functional classes for the foreground model
  int numColumns=numParms-1;
  Array1D<FunctionalClass> foregroundClasses(numColumns);
  for(int i=0 ; i<numColumns ; ++i) {
    FO_FuncClass *fo_fg=dynamic_cast<FO_FuncClass*>(env.getParm(i+1));
    if(!fo_fg) throw 
      LambdaException("invalid functional class in call to gain-of-function");
    foregroundClasses[i]=fo_fg->getClass();
  }

  // Compose foreground and background into a gain model
  TransducerTemplate *gainModel=
    global->composeGainModel(ttBackground,foregroundClasses);
  
  // Package it inside a Lambda object and return it
  return new FO_Transducer(gainModel);
}

 */


void ModelCompiler::decomposeBackground(TransducerTemplate *bg,
					StateTemplate *&insertState,
					StateTemplate *&deleteState,
					StateTemplate *&matchState,
					Lambda::Closure *&insertNorm,
					Lambda::Closure *&deleteNorm,
					Lambda::Closure *&matchNorm,
					Array2D<Lambda::Closure*> 
					 &transMatrix,
					FunctionalClass &fc)
{
  insertState=deleteState=matchState=NULL;
  bg->decomposeUnique(matchState,insertState,deleteState,transMatrix);
  if(!insertState || !deleteState || !matchState)
    throw "missing state in ModelCompiler::decomposeBackground()";
  matchNorm=transMatrix[PHMM_MATCH][PHMM_START_STOP];
  insertNorm=transMatrix[PHMM_INSERT][PHMM_START_STOP];
  deleteNorm=transMatrix[PHMM_DELETE][PHMM_START_STOP];
  StateTemplate *state=insertState;
  if(!state) state=(matchState ? matchState : deleteState);
  fc=state->getFunctionalClass(PARENT);
}



void ModelCompiler::normalize(PHMM_StateType from,PHMM_StateType to,
			      Lambda::Closure *normTerm,
			      Array2D<Lambda::Closure*> &transMatrix)
{
  return;//###
  Lambda::Closure *&entry=transMatrix[from][to];
  if(entry) entry=multiplyClosures(entry,normTerm->deepCopy());
}



Lambda::Closure *ModelCompiler::oneOverSum(Lambda::Closure *a,
					   Lambda::Closure *b)
{
  if(!a || !b) INTERNAL_ERROR;
  Lambda::AstForest *forest=
    lambda.parseSource("[t|(/ 1 (+ (exp t) (exp t)))]");
  Lambda::Closure *z=
    dynamic_cast<Lambda::Closure*>(lambda.getEvaluator().evaluate(forest));
  if(!z) throw "Error composing two closures in ModelCompiler";
  z->setStaticChain(lambda.getEnvironment().getGlobalAR());
  struct V : public Lambda::AstPostorderTraversal {
    Array1D<Lambda::Closure*> closures;
    int nextClosure;
    V(Lambda::Closure *a,Lambda::Closure *b)
      : closures(2), nextClosure(0) { closures[0]=a; closures[1]=b; }
    void process(AstCombination &combo) {
      if(combo.getParms().size()==1)
	combo.changeFunction(new AstObject(closures[nextClosure++]));
    }
  } v(a,b);
  v.traverse(z->getFunction());
  global->closures.push_back(z);//lambda.makeImmortal(z);
  return z; 
}



Lambda::Closure *ModelCompiler::oneOverOneMinus(Lambda::Closure *x)
{
  if(!x) INTERNAL_ERROR;
  Lambda::AstForest *forest=lambda.parseSource("[t|(/ 1 (- 1 (exp t)))]");
  Lambda::Closure *z=
    dynamic_cast<Lambda::Closure*>(lambda.getEvaluator().evaluate(forest));
  if(!z) throw "Error composing two closures in ModelCompiler";
  z->setStaticChain(lambda.getEnvironment().getGlobalAR());
  struct V : public Lambda::AstPostorderTraversal {
    Lambda::Closure *x;
    V(Lambda::Closure *x) : x(x) {}
    void process(AstCombination &combo) {
      if(combo.getParms().size()==1)
	combo.changeFunction(new AstObject(x));
    }
  } v(x);
  v.traverse(z->getFunction());
  global->closures.push_back(z);//lambda.makeImmortal(z);
  return z;
}



TransducerTemplate *ModelCompiler::composeLossModel(TransducerTemplate *bg,
						   Array1D<FunctionalClass>
						   &fg)
{
  // PRECONDITION: The background model has nonzero-probability transitions
  //               from q0 to both the match and deletion states.

  // Decompose background model
  Array2D<Lambda::Closure*> transMatrix;
  Lambda::Closure *insertNorm, *deleteNorm, *matchNorm; // normalizing terms
  StateTemplate *insertState, *deleteState, *matchState;
  FunctionalClass bgFC;
  decomposeBackground(bg,insertState,deleteState,matchState,insertNorm,
		      deleteNorm,matchNorm,transMatrix,bgFC);
  Lambda::Closure *q0_Mat=transMatrix[PHMM_START_STOP][PHMM_MATCH];
  Lambda::Closure *q0_Del=transMatrix[PHMM_START_STOP][PHMM_DELETE];
  if(!q0_Mat || !q0_Del)
    throw "Background lacks transition from q0 to match or delete state";

  // Construct Lambda expressions for normalization terms
  Lambda::Closure *iNorm=(insertNorm ? oneOverOneMinus(insertNorm) : NULL);
  Lambda::Closure *dNorm=(deleteNorm ? oneOverOneMinus(deleteNorm) : NULL);
  Lambda::Closure *mNorm=(matchNorm ? oneOverOneMinus(matchNorm) : NULL);
  Lambda::Closure *q0toMatch=AoverBplusC(q0_Mat,q0_Mat,q0_Del);
  Lambda::Closure *q0toDelete=AoverBplusC(q0_Del,q0_Mat,q0_Del);

  // Renormalize transition probabilities to account for the lack of a
  // transition to the STOP state
  if(iNorm) {
    normalize(PHMM_INSERT,PHMM_INSERT,iNorm,transMatrix);
    normalize(PHMM_INSERT,PHMM_MATCH,iNorm,transMatrix);
    normalize(PHMM_INSERT,PHMM_DELETE,iNorm,transMatrix);
  }
  if(dNorm) {
    normalize(PHMM_DELETE,PHMM_INSERT,dNorm,transMatrix);
    normalize(PHMM_DELETE,PHMM_MATCH,dNorm,transMatrix);
    normalize(PHMM_DELETE,PHMM_DELETE,dNorm,transMatrix);
  }
  if(mNorm) {
    normalize(PHMM_MATCH,PHMM_INSERT,mNorm,transMatrix);
    normalize(PHMM_MATCH,PHMM_MATCH,mNorm,transMatrix);
    normalize(PHMM_MATCH,PHMM_DELETE,mNorm,transMatrix);
  }

  // Build model
  TransducerTemplate *tt=new TransducerTemplate();
  tt->addState(q0);
  int numColumns=fg.size(), nextID=1;
  StateTemplate *pMat, *pDel;
  for(int i=0 ; i<numColumns ; ++i) {
    FunctionalClass fgFC=fg[i];

    // Add a new column of states
    StateTemplate *sMat=new StateTemplate(nextID++,PHMM_MATCH,fgFC,bgFC);
    StateTemplate *sDel=new StateTemplate(nextID++,PHMM_DELETE,fgFC,bgFC);
    tt->addState(sMat); tt->addState(sDel);
    StateTemplate *sIns;
    if(i>0) {
      sIns=new StateTemplate(nextID++,PHMM_INSERT,bgFC,bgFC);
      tt->addState(sIns);
      sMat->disableGainLossScoring();
      sDel->disableGainLossScoring();
      sIns->disableGainLossScoring();
    
      // Add within-column transitions
      tt->addTransition(sIns,sMat,transMatrix[PHMM_INSERT][PHMM_MATCH]);
      tt->addTransition(sIns,sDel,transMatrix[PHMM_INSERT][PHMM_DELETE]);
      tt->addTransition(sIns,sIns,transMatrix[PHMM_INSERT][PHMM_INSERT]);
    }

    // Add transitions from previous column to this column
    if(i==0) {
      tt->addTransition(q0,sMat,q0toMatch);
      tt->addTransition(q0,sDel,q0toDelete);
    }
    else {
      tt->addTransition(pDel,sDel,transMatrix[PHMM_DELETE][PHMM_DELETE]);
      tt->addTransition(pDel,sMat,transMatrix[PHMM_DELETE][PHMM_MATCH]);
      tt->addTransition(pDel,sIns,transMatrix[PHMM_DELETE][PHMM_INSERT]);
      tt->addTransition(pMat,sDel,transMatrix[PHMM_MATCH][PHMM_DELETE]);
      tt->addTransition(pMat,sMat,transMatrix[PHMM_MATCH][PHMM_MATCH]);
      tt->addTransition(pMat,sIns,transMatrix[PHMM_MATCH][PHMM_INSERT]);
    }

    // Get ready for next column
    pMat=sMat; pDel=sDel;
  }

  // Attach the final column to the STOP state
  tt->addTransition(pMat,q0,prob1);
  tt->addTransition(pDel,q0,prob1);

  // Return the new model
  return tt;
}



TransducerTemplate *ModelCompiler::composeGainModel(TransducerTemplate *bg,
						   Array1D<FunctionalClass>
						   &fg)
{
  // PRECONDITION: The background model has nonzero-probability transitions
  //               from q0 to both the match and insertion states.

  // Decompose background model
  Array2D<Lambda::Closure*> transMatrix;
  Lambda::Closure *insertNorm, *deleteNorm, *matchNorm; // normalizing terms
  StateTemplate *insertState, *deleteState, *matchState;
  FunctionalClass bgFC;
  decomposeBackground(bg,insertState,deleteState,matchState,insertNorm,
		      deleteNorm,matchNorm,transMatrix,bgFC);
  Lambda::Closure *q0_Mat=transMatrix[PHMM_START_STOP][PHMM_MATCH];
  Lambda::Closure *q0_Ins=transMatrix[PHMM_START_STOP][PHMM_INSERT];
  if(!q0_Mat || !q0_Ins)
    throw "Background lacks transition from q0 to match or insert state";

  // Construct Lambda expressions for normalization terms
  Lambda::Closure *iNorm=(insertNorm ? oneOverOneMinus(insertNorm) : NULL);
  Lambda::Closure *dNorm=(deleteNorm ? oneOverOneMinus(deleteNorm) : NULL);
  Lambda::Closure *mNorm=(matchNorm ? oneOverOneMinus(matchNorm) : NULL);
  Lambda::Closure *q0toMatch=AoverBplusC(q0_Mat,q0_Mat,q0_Ins);
  Lambda::Closure *q0toInsert=AoverBplusC(q0_Ins,q0_Mat,q0_Ins);

  // Renormalize transition probabilities to account for the lack of a
  // transition to the STOP state
  if(iNorm) {
    normalize(PHMM_INSERT,PHMM_INSERT,iNorm,transMatrix);
    normalize(PHMM_INSERT,PHMM_MATCH,iNorm,transMatrix);
    normalize(PHMM_INSERT,PHMM_DELETE,iNorm,transMatrix);
  }
  if(dNorm) {
    normalize(PHMM_DELETE,PHMM_INSERT,dNorm,transMatrix);
    normalize(PHMM_DELETE,PHMM_MATCH,dNorm,transMatrix);
    normalize(PHMM_DELETE,PHMM_DELETE,dNorm,transMatrix);
  }
  if(mNorm) {
    normalize(PHMM_MATCH,PHMM_INSERT,mNorm,transMatrix);
    normalize(PHMM_MATCH,PHMM_MATCH,mNorm,transMatrix);
    normalize(PHMM_MATCH,PHMM_DELETE,mNorm,transMatrix);
  }

  // Build model
  TransducerTemplate *tt=new TransducerTemplate();
  tt->addState(q0);
  int numColumns=fg.size(), nextID=1;
  StateTemplate *pMat, *pIns;
  for(int i=0 ; i<numColumns ; ++i) {
    FunctionalClass fgFC=fg[i];

    // Add a new column of states
    StateTemplate *sMat=new StateTemplate(nextID++,PHMM_MATCH,bgFC,fgFC);
    StateTemplate *sIns=new StateTemplate(nextID++,PHMM_INSERT,bgFC,fgFC);
    tt->addState(sMat); tt->addState(sIns); 
    StateTemplate *sDel;
    if(i>0) {
      sDel=new StateTemplate(nextID++,PHMM_DELETE,bgFC,bgFC);
      tt->addState(sDel);
      sMat->disableGainLossScoring();
      sDel->disableGainLossScoring();
      sIns->disableGainLossScoring();

      // Add within-column transitions
      tt->addTransition(sDel,sMat,transMatrix[PHMM_DELETE][PHMM_MATCH]);
      tt->addTransition(sDel,sIns,transMatrix[PHMM_DELETE][PHMM_INSERT]);
      tt->addTransition(sDel,sDel,transMatrix[PHMM_DELETE][PHMM_DELETE]);
    }

    // Add transitions from previous column to this column
    if(i==0) {
      tt->addTransition(q0,sMat,q0toMatch);
      tt->addTransition(q0,sIns,q0toInsert);
    }
    else {
      tt->addTransition(pIns,sIns,transMatrix[PHMM_INSERT][PHMM_INSERT]);
      tt->addTransition(pIns,sMat,transMatrix[PHMM_INSERT][PHMM_MATCH]);
      tt->addTransition(pIns,sDel,transMatrix[PHMM_INSERT][PHMM_DELETE]);
      tt->addTransition(pMat,sIns,transMatrix[PHMM_MATCH][PHMM_INSERT]);
      tt->addTransition(pMat,sMat,transMatrix[PHMM_MATCH][PHMM_MATCH]);
      tt->addTransition(pMat,sDel,transMatrix[PHMM_MATCH][PHMM_DELETE]);
    }

    // Get ready for next column
    pMat=sMat; pIns=sIns;
  }

  // Attach the final column to the STOP state
  tt->addTransition(pMat,q0,prob1);
  tt->addTransition(pIns,q0,prob1);

  // Return the new model
  return tt;
}



TransducerTemplate *ModelCompiler::composeSimpleLossModel(
    TransducerTemplate *bg,Array1D<FunctionalClass> &fg)
{
  TransducerTemplate *tt=new TransducerTemplate();
  tt->addState(q0);
  int numColumns=fg.size(), nextID=1;
  StateTemplate *pMat;
  for(int i=0 ; i<numColumns ; ++i) {
    FunctionalClass fgFC=fg[i];
    StateTemplate *sMat=new StateTemplate(nextID++,PHMM_DELETE,fgFC,
					  FunctionalClass::getBackground());
    tt->addState(sMat);
    if(i==0) tt->addTransition(q0,sMat,prob1);
    else {
      sMat->disableGainLossScoring();
      tt->addTransition(pMat,sMat,prob1);
    }
    pMat=sMat;
  }
  tt->addTransition(pMat,q0,prob1);
  return tt;
}



TransducerTemplate *ModelCompiler::composeSimpleGainModel(
    TransducerTemplate *bg,Array1D<FunctionalClass> &fg)
{
  TransducerTemplate *tt=new TransducerTemplate();
  tt->addState(q0);
  int numColumns=fg.size(), nextID=1;
  StateTemplate *pMat;
  FunctionalClass bgFC=FunctionalClass::getBackground();
  for(int i=0 ; i<numColumns ; ++i) {
    FunctionalClass fgFC=fg[i];
    StateTemplate *sMat=new StateTemplate(nextID++,PHMM_INSERT,bgFC,fgFC);
    tt->addState(sMat);
    if(i==0) tt->addTransition(q0,sMat,prob1);
    else {
      sMat->disableGainLossScoring();
      tt->addTransition(pMat,sMat,prob1);
    }
    pMat=sMat;
  }
  tt->addTransition(pMat,q0,prob1);
  return tt;
}



TransducerTemplate *ModelCompiler::compose2tierLoss(TransducerTemplate *bg,
						   Array1D<FunctionalClass>
						   &fg)
{
  // PRECONDITION: The background model has nonzero-probability transitions
  //               from q0 to both the match and deletion states.

  // Decompose background model
  Array2D<Lambda::Closure*> transMatrix;
  Lambda::Closure *insertNorm, *deleteNorm, *matchNorm; // normalizing terms
  StateTemplate *insertState, *deleteState, *matchState;
  FunctionalClass bgFC;
  decomposeBackground(bg,insertState,deleteState,matchState,insertNorm,
		      deleteNorm,matchNorm,transMatrix,bgFC);
  Lambda::Closure *q0_Mat=transMatrix[PHMM_START_STOP][PHMM_MATCH];
  Lambda::Closure *q0_Del=transMatrix[PHMM_START_STOP][PHMM_DELETE];
  if(!q0_Mat || !q0_Del)
    throw "Background lacks transition from q0 to match or delete state";

  // Build model
  TransducerTemplate *tt=new TransducerTemplate();
  tt->addState(q0);
  int numColumns=fg.size(), nextID=1;
  StateTemplate *pMat, *pDel;
  for(int i=0 ; i<numColumns ; ++i) {
    FunctionalClass fgFC=fg[i];

    // Add a new column of states
    StateTemplate *sMat=new StateTemplate(nextID++,PHMM_MATCH,fgFC,bgFC);
    StateTemplate *sDel=new StateTemplate(nextID++,PHMM_DELETE,fgFC,bgFC);
    tt->addState(sMat); tt->addState(sDel);
    if(i>0) {
      sMat->disableGainLossScoring();
      sDel->disableGainLossScoring();
    }

    // Add transitions from previous column to this column
    if(i==0) {
      tt->addTransition(q0,sMat,q0_Mat);
      tt->addTransition(q0,sDel,q0_Del);
    }
    else {
      tt->addTransition(pDel,sDel,transMatrix[PHMM_DELETE][PHMM_DELETE]);
      tt->addTransition(pDel,sMat,transMatrix[PHMM_DELETE][PHMM_MATCH]);
      tt->addTransition(pMat,sDel,transMatrix[PHMM_MATCH][PHMM_DELETE]);
      tt->addTransition(pMat,sMat,transMatrix[PHMM_MATCH][PHMM_MATCH]);
    }

    // Get ready for next column
    pMat=sMat; pDel=sDel;
  }

  // Attach the final column to the STOP state
  tt->addTransition(pMat,q0,prob1);
  tt->addTransition(pDel,q0,prob1);

  // Return the new model
  return tt;
}



TransducerTemplate *ModelCompiler::compose2tierGain(TransducerTemplate *bg,
						   Array1D<FunctionalClass>
						   &fg)
{
  // PRECONDITION: The background model has nonzero-probability transitions
  //               from q0 to both the match and insertion states.

  // Decompose background model
  Array2D<Lambda::Closure*> transMatrix;
  Lambda::Closure *insertNorm, *deleteNorm, *matchNorm; // normalizing terms
  StateTemplate *insertState, *deleteState, *matchState;
  FunctionalClass bgFC;
  decomposeBackground(bg,insertState,deleteState,matchState,insertNorm,
		      deleteNorm,matchNorm,transMatrix,bgFC);
  Lambda::Closure *q0_Mat=transMatrix[PHMM_START_STOP][PHMM_MATCH];
  Lambda::Closure *q0_Ins=transMatrix[PHMM_START_STOP][PHMM_INSERT];
  if(!q0_Mat || !q0_Ins)
    throw "Background lacks transition from q0 to match or insert state";

  // Build model
  TransducerTemplate *tt=new TransducerTemplate();
  tt->addState(q0);
  int numColumns=fg.size(), nextID=1;
  StateTemplate *pMat, *pIns;
  for(int i=0 ; i<numColumns ; ++i) {
    FunctionalClass fgFC=fg[i];

    // Add a new column of states
    StateTemplate *sMat=new StateTemplate(nextID++,PHMM_MATCH,bgFC,fgFC);
    StateTemplate *sIns=new StateTemplate(nextID++,PHMM_INSERT,fgFC,fgFC);
    tt->addState(sMat); tt->addState(sIns); 
    if(i>0) {
      sMat->disableGainLossScoring();
      sIns->disableGainLossScoring();
    }

    // Add transitions from previous column to this column
    if(i==0) {
      tt->addTransition(q0,sMat,q0_Mat);
      tt->addTransition(q0,sIns,q0_Ins);
    }
    else {
      tt->addTransition(pIns,sIns,transMatrix[PHMM_INSERT][PHMM_INSERT]);
      tt->addTransition(pIns,sMat,transMatrix[PHMM_INSERT][PHMM_MATCH]);
      tt->addTransition(pMat,sIns,transMatrix[PHMM_MATCH][PHMM_INSERT]);
      tt->addTransition(pMat,sMat,transMatrix[PHMM_MATCH][PHMM_MATCH]);
    }

    // Get ready for next column
    pMat=sMat; pIns=sIns;
  }

  // Attach the final column to the STOP state
  tt->addTransition(pMat,q0,prob1);
  tt->addTransition(pIns,q0,prob1);

  // Return the new model
  return tt;
}



LambdaObject *ModelCompiler::ff_bindingSite(RunTimeEnvironment &env,
					    ContinuationType &) 
{
  // syntax: (binding-site functionalElementType)

  // Process parms
  if(env.getNumParms()!=1) 
    throw LambdaException("'binding-site' requires a parameter");
  StateTemplate * const q0=global->q0;
  Lambda::Closure * const prob1=global->prob1;

  // Get the functional classes for the foreground model
  FO_FuncElemType *fo_funcElemType=
    dynamic_cast<FO_FuncElemType*>(env.getParm(0));
  if(!fo_funcElemType) 
    throw LambdaException("'binding-site' requires functional element type");

  // Compose model
  int numColumns=fo_funcElemType->getNumClasses();
  TransducerTemplate *tt=new TransducerTemplate();
  tt->addState(q0);
  StateTemplate *prevState=q0;
  for(int i=0 ; i<numColumns ; ++i) {
    FunctionalClass f=fo_funcElemType->getIthClass(i)->getClass();
    StateTemplate *state=new StateTemplate(i+1,PHMM_MATCH,f,f);
    tt->addState(state);
    tt->addTransition(prevState,state,prob1);
    prevState=state;
  }
  tt->addTransition(prevState,q0,prob1);
  
  // Package it inside a Lambda object and return it
  FO_Transducer *fo=new FO_Transducer(tt);
  global->registerFO(fo);
  return fo;
}


/*
LambdaObject *ModelCompiler::ff_bindingSite(RunTimeEnvironment &env,
					    ContinuationType &) 
{
  // syntax: (binding-site fc1 fc2 fc3 ... fcN)

  // Process parms
  int numParms=env.getNumParms();
  if(numParms<1)
    throw LambdaException("(binding-site) requires >0 parameters");
  StateTemplate * const q0=global->q0;
  Lambda::Closure * const prob1=global->prob1;

  // Get the functional classes for the foreground model
  int numColumns=numParms;
  Array1D<FunctionalClass> fc(numColumns);
  for(int i=0 ; i<numColumns ; ++i) {
    FO_FuncClass *fo_fg=dynamic_cast<FO_FuncClass*>(env.getParm(i));
    if(!fo_fg) throw 
      LambdaException("invalid functional class in call to binding-site");
    fc[i]=fo_fg->getClass();
  }

  // Compose model
  TransducerTemplate *tt=new TransducerTemplate();
  tt->addState(q0);
  StateTemplate *prevState=q0;
  for(int i=0 ; i<numColumns ; ++i) {
    FunctionalClass f=fc[i];
    StateTemplate *state=new StateTemplate(i+1,PHMM_MATCH,f,f);
    tt->addState(state);
    tt->addTransition(prevState,state,prob1);
    prevState=state;
  }
  tt->addTransition(prevState,q0,prob1);
  
  // Package it inside a Lambda object and return it
  return new FO_Transducer(tt);
}
 */


LambdaObject *ModelCompiler::ff_gainFactor(RunTimeEnvironment &env,
					   ContinuationType &)
{
  // syntax: (gain-factor [t|...])

  // Process parms
  LambdaObject *parm0=env.getParm(0);

  // Get the value
  Lambda::Closure *closure=dynamic_cast<Lambda::Closure*>(parm0);
  if(!closure) throw LambdaException("gain-factor requires closure parm");
  global->closures.push_back(closure);//env.makeImmortal(closure);

  // Allocate a wrapper object
  LambdaObject *fo=new FO_GainLossFactor(closure,GLT_GAIN);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_scaleTree(RunTimeEnvironment &env,
					  ContinuationType &)
{
  // syntax: (scale-tree tree factor)

  // Get the tree
  FO_Phylogeny *phylogeny=dynamic_cast<FO_Phylogeny*>(env.getParm(0));
  if(!phylogeny) throw 
    LambdaException("scale-tree requires phylogeny parameter");
  Phylogeny *tree=phylogeny->getPhylogeny();

  // Get the value
  LambdaFloat *parm1=dynamic_cast<LambdaFloat*>(env.getParm(1));
  if(!parm1)
    throw LambdaException("scale-tree requires scale parameter");
  float scale=parm1->getValue();

  // Scale the tree
  tree->scaleBranches(scale/global->prevTreeScale);
  global->prevTreeScale=scale;

  // Return the tree
  return phylogeny;
}



LambdaObject *ModelCompiler::ff_lossFactor(RunTimeEnvironment &env,
					   ContinuationType &)
{
  // syntax: (loss-factor [t|...])

  // Process parms
  LambdaObject *parm0=env.getParm(0);

  // Get the value
  Lambda::Closure *closure=dynamic_cast<Lambda::Closure*>(parm0);
  if(!closure) throw LambdaException("loss-factor requires closure parm");
  global->closures.push_back(closure);//env.makeImmortal(closure);

  // Allocate a wrapper object
  LambdaObject *fo=new FO_GainLossFactor(closure,GLT_LOSS);
  global->registerFO(fo);
  return fo;
}



LambdaObject *ModelCompiler::ff_retentionFactor(RunTimeEnvironment &env,
						ContinuationType &)
{
  // syntax: (retention-factor [t|...])

  // Process parms
  LambdaObject *parm0=env.getParm(0);

  // Get the value
  Lambda::Closure *closure=dynamic_cast<Lambda::Closure*>(parm0);
  if(!closure) 
    throw LambdaException("retention-factor requires closure parm");
  global->closures.push_back(closure);//env.makeImmortal(closure);

  // Allocate a wrapper object
  LambdaObject *fo=new FO_GainLossFactor(closure,GLT_RETENTION);
  global->registerFO(fo);
  return fo;
}



float ModelCompiler::lookupFloat(const String &variableName)
{
  LambdaObject *obj=lambda.lookupGlobal(variableName);
  if(!obj) throw String("Undefined variable: ")+variableName;
  LambdaFloat *flt=dynamic_cast<LambdaFloat*>(obj);
  if(!flt) throw variableName+" not a float";
  return flt->getValue();
}



void ModelCompiler::setGCthreshold(int t)
{
  lambda.setGCthreshold(t);
}



void ModelCompiler::registerFO(Lambda::ForeignObject*fo) 
{
  /*
  switch(fo->getTag())
    {
    case FOT_RATE_MATRIX:
    case FOT_RATE_MATRICES:
    case FOT_PHYLOGENY:
    case FOT_STATES:
    case FOT_STATE:
    case FOT_TRANSITION:
      return;//
    case FOT_TRANSDUCER: // BAD
      //return;
    case FOT_STATE_TYPE:
      //return;
    case FOT_MACROSTATE:
    case FOT_MACROSTATES:
      //return;
    case FOT_MACROTRANS:
    case FOT_FUNC_CLASS:
    case FOT_GAINLOSS_FACTOR:
      break;
    }
  */
  foreignObjects.push_back(fo);
}



/****************************************************************
                         struct ModelParm
****************************************************************/

const double pi=2.0*acos(0.0);

ModelParm::ModelParm(String n,float min,float max)
  : name(n), min(min), max(max), range(max-min)
{
  // ctor
}



float ModelParm::constrain(float x)
{
  // This function maps an unconstrained value produced by the optimizer
  // into the range acceptable for this particular parameter

  return min+range*(0.5+atan(x)/pi); 
}



float ModelParm::unconstrain(float y)
{
  // This function maps this parameter into the wider range that the
  // optimizer works with (-inf to +inf)

  return tan(pi*((y-min)/range-0.5));
}



/****************************************************************
                     FO_FuncElemType methods
****************************************************************/

void FO_FuncElemType::pushAccessibles(MarkStack &s)
{
  int n=classes.size();
  FO_FuncClass **p=&classes[0];
  for(int i=0 ; i<n ; ++i) {
    s.push(*p);
    ++p;
  }
}


/****************************************************************
                     FO_RateMatrices methods
****************************************************************/
void FO_RateMatrices::pushAccessibles(MarkStack &s)
{
  int n=matrices.size();
  FO_RateMatrix **p=&matrices[0];
  for(int i=0 ; i<n ; ++i) {
    s.push(*p);
    ++p;
  }
}


/****************************************************************
                       FO_States methods
****************************************************************/
void FO_States::pushAccessibles(MarkStack &s)
{
  Array1D<FO_State*> &states=getStates() ;
  int n=states.size();
  FO_State **p=&states[0];
  for(int i=0 ; i<n ; ++i) {
    s.push(*p);
    ++p;
  }
}

/****************************************************************
                     FO_Macrostates methods
****************************************************************/
void FO_Macrostates::pushAccessibles(MarkStack &s)
{
  Array1D<FO_Macrostate*> &states=getStates() ;
  int n=states.size();
  FO_Macrostate **p=&states[0];
  for(int i=0 ; i<n ; ++i) {
    s.push(*p);
    ++p;
  }
}


