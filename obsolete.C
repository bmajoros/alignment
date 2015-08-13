  /*
  // Resample the alignment path along that branch
  int rootID=tree->getRoot()->getID();
  BitSet &cladeMembers=childTaxon.getCladeMembers();
  BitSet nonCladeMembers=cladeMembers;
  nonCladeMembers.complement();
  AlignmentView cladeView(*fullAlignment,cladeMembers,gapSymbol,trackMap);
  AlignmentView nonCladeView(*fullAlignment,nonCladeMembers,gapSymbol,
			     trackMap);
  ProfileBackward backward(*hmm,cladeView,nonCladeView,*tree,trackMap,
			   alphabetMap);
  ProfileSampler sampler(*hmm,backward,nonCladeView,cladeView);
  double newScore, oldScore;
  StatePath *path=sampler.samplePath(newScore,oldScore);
  
  // Convert the HMM path to an alignment
  MultSeqAlignment *newAlignment=
    decodeAlignment(*path,cladeView,nonCladeView,hmm);
  delete path;
  inferGapPatterns(*newAlignment);
  double logP=likelihood();
  newAlignment->setScore(logP);
  cout<<"logP="<<logP<<endl;

  if(shouldHillClimb) {
    if(logP>alignmentScore) {
      alignmentScore=logP;
      delete fullAlignment;
      fullAlignment=newAlignment;
      fullAlignment->setScore(logP);
      fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,true);
    }
    else delete newAlignment;
  }
  else {
    proposalScore=oldScore;
    double P_old_given_new=oldScore; // P(old|new)
    double P_new_given_old=newScore; // P(new|old)
    double P_new=logP; // P(new)
    double P_old=alignmentScore; // P(old)
    double logHastings=P_new - P_old + P_old_given_new - P_new_given_old;
    double hastings=exp(logHastings);
    if(samplerType==METROPOLIS_HASTINGS) 
      cout<<"hastings ratio="<<hastings<<" ("<<proposalScore<<"-"<<newScore
	  <<"+"<<logP<<"-"<<alignmentScore<<" = "<<logHastings<<")"<<endl;
    double r=Random0to1();
    if(r<=hastings || samplerType==GIBBS) {
      if(samplerType==METROPOLIS_HASTINGS) cout<<"ACCEPT"<<endl;
      alignmentScore=logP;
      proposalScore=newScore;
      delete fullAlignment;
      fullAlignment=newAlignment;
      fullAlignment->setScore(logP);
      fullAlignment->printSlice(cout,0,fullAlignment->getLength(),'+',60,true);
    }
    else {
      cout<<"REJECT (r="<<r<<")"<<endl;
      delete newAlignment;
    }
  }
}
  */



/****************************************************************
                  APES::getTransProb_nonLinked()
 ****************************************************************/
double APES::getTransProbs_nonLinked(BranchHMM &hmm,
				     GapPattern &parentGaps,
				     GapPattern &childGaps)
{
  // ### This function will have to change when functional states are
  // incorporated into the model

  Vector<STATE> matchStates, insertionStates, deletionStates;
  hmm.getStatesOfType(PHMM_MATCH,matchStates);
  hmm.getStatesOfType(PHMM_INSERT,insertionStates);
  hmm.getStatesOfType(PHMM_DELETE,deletionStates);
  STATE M=matchStates[0], I=insertionStates[0], D=deletionStates[0];//### hack
  double logP=0;
  int L=parentGaps.getLength();
  STATE prevState=0, nextState;
  for(int j=0 ; j<L ; ++j) {
    GapPatternElement p=parentGaps[j], c=childGaps[j];
    if(p==GPE_RESIDUE)
      if(c==GPE_RESIDUE) nextState=M;
      else nextState=D;
    else
      if(c==GPE_RESIDUE) nextState=I;
      else continue; // both are gaps
    logP+=hmm.getTransP(prevState,nextState);
    prevState=nextState;
  }
  logP+=hmm.getTransP(prevState,0);
  return logP;
}



/****************************************************************
                  APES::likelihood_nonLinked()
 ****************************************************************/
double APES::likelihood_nonLinked(const MultSeqAlignment &A)
{
  // First, compute the contribution of the emission terms by applying
  // Felsenstein's algorithm to all columns of the alignment
  MultSeqAlignment augmented(alphabet,gapSymbol);
  int numTaxa=taxa.size(), L=A.getLength();
  for(int i=0 ; i<numTaxa ; ++i) {
    const String name=taxa[i].getName();
    if(taxa[i].getNode()->getNodeType()==LEAF_NODE) {
      AlignmentSeq *newTrack=new AlignmentSeq(name,i,A.getGapSymbols());
      newTrack->getSeq()=A.getTrackByName(name).getSeq();
      augmented.addTrack(newTrack);
    }
    else augmented.findOrCreateTrack(name).extendToLength(L,gapSymbol);
  }
  double logP=0;
  //  int L=augmented.getLength();
  PhylogenyNode *root=tree->getRoot();
  ProfileFelsenstein F(augmented,PureDnaAlphabet::global(),alphabetMap);


  //#########
  STATE state=1; // ### NEED TO CHANGE THIS WHEN I ADD FUNCTIONAL STATES
  //#########



  for(int i=0 ; i<L ; ++i)
    logP+=F.logLikelihood(i,root,state);

  // Next, process each branch individually to assess the transition 
  // probabilities for the transducers
  int n=branches.size();
  for(int i=0 ; i<n ; ++i) {
    BranchAttributes &attr=branchAttributes[i];
    BranchHMM *hmm=attr.getHMM();
    int parentID=attr.getParentID(), childID=attr.getChildID();
    GapPattern &parentGaps=*taxa[parentID].getGapPattern();
    GapPattern &childGaps=*taxa[childID].getGapPattern();
    logP+=getTransProbs_nonLinked(*hmm,parentGaps,childGaps);
  }
  return logP;
}


