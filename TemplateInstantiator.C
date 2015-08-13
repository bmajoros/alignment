/****************************************************************
 TemplateInstantiator.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "TemplateInstantiator.H"
#include "BOOM/DnaDashAlphabet.H"
#include "BOOM/PureDnaAlphabet.H"
#include "BOOM/lambda/LambdaException.H"
#include "PhyLib/SubstitutionMixtureModel.H"
using namespace std;
using namespace BOOM;



TemplateInstantiator::TemplateInstantiator(LambdaAPI &lambda)
  : alphabet(DnaDashAlphabet::global()), lambda(lambda),
    alphabetMap(PureDnaAlphabet::global(),DnaDashAlphabet::global())
{
  gapSymbol=alphabet.lookup('-');
  initCaller();
}



void TemplateInstantiator::initCaller()
{
  parmNode=new AstFloatLit(0.0);
  AstForest *args=new AstForest;
  args->appendTree(parmNode);
  caller=lambda.makeCombination(NULL,args,false);
  closureWrapper=new AstObject(NULL);
  caller->changeFunction(closureWrapper);
  
  //###
  lambda.makeImmortal(parmNode);
  lambda.makeImmortal(args);
  //delete args;
  lambda.makeImmortal(caller);
  lambda.makeImmortal(closureWrapper);
  //###
}



BranchHMM *TemplateInstantiator::instantiate(TransducerTemplate *T,
					     double branchLength,
					     ostream *os)
{
  parmNode->changeValue(branchLength);
  int numStates=T->getNumStates(), numTransitions=T->getNumTransitions();
  PairHMM hmm(alphabet,gapSymbol,numStates);
  Array1D<SubstitutionMatrix*> matrices(numStates);
  matrices.setAllTo(NULL);
  hmm.setStateType(0,PHMM_START_STOP);
  Array1D<double> eqFreqs;
  int nAlpha=PureDnaAlphabet::global().size();
  for(int i=0 ; i<numStates ; ++i) {
    StateTemplate &state=T->getIthState(i);
    STATE id=state.getStateID();
    FunctionalClass fcParent=state.getFunctionalClass(PARENT);
    FunctionalClass fcChild=state.getFunctionalClass(CHILD);
    SubstitutionMatrix *Pt;
    if(fcParent==fcChild) {
      RateMatrix *Q=state.getMatrix(PARENT);
      Pt=(Q ? Q->instantiate(branchLength) : NULL);
    }
    else {
      RateMatrix *Q1=state.getMatrix(PARENT), *Q2=state.getMatrix(CHILD);
      Pt=new SubstitutionMixtureModel(*Q1,*Q2,branchLength);
    }
    matrices[id]=Pt;
    hmm.setStateType(id,state.getStateType());
    if(Pt) {
      Pt->getEqFreqs(eqFreqs); // resizes eqFreqs (OK)
      for(BOOM::Symbol j=0 ; j<nAlpha ; ++j)
	for(BOOM::Symbol k=0 ; k<nAlpha ; ++k) {
	  double emitP=(*Pt)(j,k)*eqFreqs[j];
	  hmm.setEmitP(id,alphabetMap(j),alphabetMap(k),emitP); 
	}
    }
    if(os) {
      *os<<"#"<<id<<": "<<state.getStateType()<<" "
	 <<state.crossFunctionalType()<<" "
	 <<state.getFunctionalClass(PARENT).getName()<<"/"
	 <<state.getFunctionalClass(CHILD).getName()<<" "
	 <<endl
	 <<"subst matrix for state "<<id<<":\n";
      if(Pt) {
	*os<<*Pt<<endl;
	Array1D<double> eq;
	Pt->getEqFreqs(eq);
	*os<<"  EQ: "<<eq<<endl;
      }
      else *os<<"NULL"<<endl;
    }
  }
  for(STATE i=0 ; i<numStates ; ++i)
    for(STATE j=0 ; j<numStates ; ++j)
      hmm.setTransP(i,j,0.0);
  for(int i=0 ; i<numTransitions ; ++i) {
    TransitionTemplate &trans=T->getIthTransition(i);
    StateTemplate *from=trans.getFrom(), *to=trans.getTo();
    STATE fromID=from->getStateID(), toID=to->getStateID();
    Lambda::Closure *f=trans.getTransProbFunc();
    closureWrapper->changeObject(f);

    //cout<<"crashes here:"<<endl;
    LambdaObject *result=lambda.evaluate(caller); // ### SEG FAULT ####
    //cout<<"did it crash?"<<endl;

    if(!lambda.isNumeric(result))
      throw LambdaException("transition function does not "
			    "evaluate to a real number");
    float transP=lambda.asFloat(result);
    if(transP<0.0) {
      cout<<transP<<" : Negative transition probability between states "
	  <<fromID<<" and "<<toID<<" numStates="<<numStates<<" from="
	  <<from->getStateType()<<" to="<<to->getStateType()<<endl;
      cout<<"lambda="<<lambda.asFloat(lambda.lookupGlobal("lambda"))<<endl;
      cout<<"mu="<<lambda.asFloat(lambda.lookupGlobal("mu"))<<endl;
      cout<<"ratio-bcd="<<lambda.asFloat(lambda.lookupGlobal("ratio-bcd"))<<endl;
      cout<<"ratio-cad="<<lambda.asFloat(lambda.lookupGlobal("ratio-cad"))<<endl;
      cout<<"ratio-gt="<<lambda.asFloat(lambda.lookupGlobal("ratio-gt"))<<endl;
      cout<<"ratio-hb="<<lambda.asFloat(lambda.lookupGlobal("ratio-hb"))<<endl;
      cout<<"ratio-kni="<<lambda.asFloat(lambda.lookupGlobal("ratio-kni"))<<endl;
      cout<<"ratio-Kr="<<lambda.asFloat(lambda.lookupGlobal("ratio-Kr"))<<endl;
      cout<<"ratio-tll="<<lambda.asFloat(lambda.lookupGlobal("ratio-tll"))<<endl;
      cout<<"density="<<lambda.asFloat(lambda.lookupGlobal("density"))<<endl;
      lambda.print(f->getFunction(),cout);
      throw String("Negative transition probability between states ")+
	fromID+" and "+toID;
    }
    else if(transP>1.0) 
      throw String("Transition probability >1 between states ")+
	fromID+" and "+toID;
    hmm.setTransP(fromID,toID,transP);
    if(os) *os<<from->getFunctionalClass(PARENT).getName()<<"/"
	      <<from->getFunctionalClass(CHILD).getName()<<":"
	      <<from->getStateType()<<"("<<fromID<<")"<<"->"
	      <<to->getFunctionalClass(PARENT).getName()<<"/"
	      <<to->getFunctionalClass(CHILD).getName()<<":"
	      <<to->getStateType()
	      <<"("<<toID<<")="<<transP<<endl;
  }
  BranchHMM *branchHMM=new BranchHMM(hmm,matrices);
  for(int i=0 ; i<numStates ; ++i) {
    StateTemplate &state=T->getIthState(i);
    STATE id=state.getStateID();
    branchHMM->setFunctionalClass(id,PARENT,state.getFunctionalClass(PARENT));
    branchHMM->setFunctionalClass(id,CHILD,state.getFunctionalClass(CHILD));
    branchHMM->setGainLossScorePolicy(id,state.shouldScoreGainLoss());
  }
  setGainLossFactors(branchHMM,T,branchLength);
  branchHMM->initSubstGenerators();
  branchHMM->normalizeTransProbs();
  branchHMM->convertToLogs();
  branchHMM->initFuncStateTypeMatrix();
  branchHMM->initPredSuccLists();
  lambda.checkGarbageLevel();
  return branchHMM;
}



void TemplateInstantiator::setGainLossFactor(BranchHMM *hmm,
					     TransducerTemplate *t,
					     GainLossType glt,
					     double branchLength)
{
  Lambda::Closure *c=t->getGainLossFactor(glt);
  if(c) {
    //lambda.makeImmortal(c);//###
    closureWrapper->changeObject(c);
    LambdaObject *result=lambda.evaluate(caller);
    if(!lambda.isNumeric(result))
      throw LambdaException("gain/loss factor does not evaluate to number");
    float f=lambda.asFloat(result);
    hmm->setGainLossFactor(glt,f);
  }
}



void TemplateInstantiator::setGainLossFactors(BranchHMM *hmm,
					      TransducerTemplate *t,
					      double branchLength)
{
  setGainLossFactor(hmm,t,GLT_GAIN,branchLength);
  setGainLossFactor(hmm,t,GLT_LOSS,branchLength);
  setGainLossFactor(hmm,t,GLT_RETENTION,branchLength);
}


