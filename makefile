CC		= g++
MPICC		= mpiCC
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -w
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
LIBS		= -lpthread  -lgsl -lm -lgslcblas -LPhyLib -lphylib -LBOOM/lambda -llambda -LEGGS -lEGGS -LBOOM -lBOOM 
PHYLIB		= PhyLib
PHMM		= PairHMM

EGGS:
	ln -s ../EGGS

$(OBJ):
	mkdir $(OBJ)

clean:
	rm obj/*.o

#---------------------------------------------------------
CLASSES = \
		$(OBJ)/PrecomputedEmissions.o \
		$(OBJ)/GainLossEvent.o \
		$(OBJ)/BirthDeathMatrix.o \
		$(OBJ)/PostLeafBackward.o \
		$(OBJ)/ConstrainedBackward.o \
		$(OBJ)/PosteriorBackward.o \
		$(OBJ)/CollapsedOrthologyMatrix.o \
		$(OBJ)/ResidueOrthologyGraph.o \
		$(OBJ)/SparseMatrix3D.o \
		$(OBJ)/SparseMatrix2D.o \
		$(OBJ)/PosteriorMatrix.o \
		$(OBJ)/BandedFB_Base.o \
		$(OBJ)/BandedForward.o \
		$(OBJ)/BandedBackward.o \
		$(OBJ)/WigBinary.o \
		$(OBJ)/HirschPosteriors.o \
		$(OBJ)/FunctionalDollo.o \
		$(OBJ)/PosetBuilder.o \
		$(OBJ)/HirschThread.o \
		$(OBJ)/FunctionalElement.o \
		$(OBJ)/LinkParsimony.o \
		$(OBJ)/GainLossType.o \
		$(OBJ)/Posteriors.o \
		$(OBJ)/LinkSampler.o \
		$(OBJ)/LinkViterbi.o \
		$(OBJ)/LinkForward.o \
		$(OBJ)/SiblingViterbi.o \
		$(OBJ)/ResidueAddress.o \
		$(OBJ)/AlignmentBuilder.o \
		$(OBJ)/LinkBackward.o \
		$(OBJ)/HspMatrix.o \
		$(OBJ)/Bander.o \
		$(OBJ)/FunctionalParse.o \
		$(OBJ)/FunctionalClass.o \
		$(OBJ)/BranchHMM.o \
		$(OBJ)/TemplateInstantiator.o \
		$(OBJ)/TransducerTemplate.o \
		$(OBJ)/ModelCompiler.o \
		$(OBJ)/BranchAttributes.o \
		$(OBJ)/GapPattern.o \
		$(OBJ)/GapPatternAlphabet.o \
		$(OBJ)/ProfileSampler.o \
		$(OBJ)/AlignmentView.o \
		$(OBJ)/Transducer.o \
		$(OBJ)/ProfileBackward.o \
		$(OBJ)/ProfileFelsenstein.o \
		$(OBJ)/LinkFelsenstein.o \
		$(OBJ)/LossyFelsenstein.o \
		$(OBJ)/GainLossFelsenstein.o \
		$(OBJ)/ProfileViterbi.o \
		$(OBJ)/Taxon.o \
		$(OBJ)/State.o \
		$(OBJ)/StatePath.o \
		$(OBJ)/PairHMM.o \
		$(OBJ)/Viterbi.o \
		$(OBJ)/HirschbergFrame.o \
		$(OBJ)/Hirschberg.o \
		$(OBJ)/SiblingHirschberg.o \
		$(OBJ)/HirschPass.o \
		$(OBJ)/HirschForwardSum.o \
		$(OBJ)/HirschBackwardSum.o \
		$(OBJ)/HirschBackwardSumPair.o \
		$(OBJ)/HirschForwardMax.o \
		$(OBJ)/HirschBackwardMax.o \
		$(OBJ)/BandingPattern.o
#---------------------------------------------------------


#########################################################
#                   APPLICATIONS
#########################################################
$(OBJ)/AVES.o:\
		AVES.C
	$(MPICC) $(CFLAGS) -o $(OBJ)/AVES.o -c \
		AVES.C
#--------------------------------------------------------
$(OBJ)/MAFIA.o:\
		MAFIA.C\
		AVES.H
	$(MPICC) $(CFLAGS) -o $(OBJ)/MAFIA.o -c \
		MAFIA.C
#--------------------------------------------------------
MAFIA: \
		$(OBJ) \
		$(OBJ)/AVES.o \
		$(OBJ)/MAFIA.o \
		EGGS \
		$(CLASSES)
	$(MPICC) $(LDFLAGS) -o MAFIA \
		$(OBJ)/AVES.o \
		$(OBJ)/MAFIA.o \
		$(CLASSES) \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/convert-to-aln.o:\
		convert-to-aln.C\
		AVES.H
	$(MPICC) $(CFLAGS) -o $(OBJ)/convert-to-aln.o -c \
		convert-to-aln.C
#--------------------------------------------------------
convert-to-aln: \
		$(OBJ)/AVES.o \
		$(OBJ)/convert-to-aln.o \
		$(CLASSES)
	$(MPICC) $(LDFLAGS) -o convert-to-aln \
		$(OBJ)/AVES.o \
		$(OBJ)/convert-to-aln.o \
		$(CLASSES) \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/debug-lengths.o:\
		debug-lengths.C
	$(CC) $(CFLAGS) -o $(OBJ)/debug-lengths.o -c \
		debug-lengths.C
#--------------------------------------------------------
debug-lengths: \
		$(OBJ)/debug-lengths.o \
		$(OBJ)/Sampler.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(CLASSES)
	$(CC) $(LDFLAGS) -o debug-lengths \
		$(OBJ)/debug-lengths.o \
		$(OBJ)/Sampler.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(CLASSES) \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/test-mpi.o:\
		test-mpi.C
	$(CC) $(CFLAGS) -o $(OBJ)/test-mpi.o -c \
		test-mpi.C
#--------------------------------------------------------
test-mpi: \
		$(OBJ)/test-mpi.o \
		$(CLASSES)
	mpiCC $(LDFLAGS) -o test-mpi \
		$(OBJ)/test-mpi.o \
		$(CLASSES) \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/train-background.o:\
		train-background.C
	$(CC) $(CFLAGS) -o $(OBJ)/train-background.o -c \
		train-background.C
#---------------------------------------------------------
train-background: \
		$(OBJ)/train-background.o \
		$(CLASSES)
	$(CC) $(LDFLAGS) -o train-background \
		$(OBJ)/train-background.o \
		$(CLASSES) \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/evos.o:\
		evos.C
	$(CC) $(CFLAGS) -o $(OBJ)/evos.o -c \
		evos.C
#---------------------------------------------------------
evos: \
		$(CLASSES) \
		$(OBJ)/evos.o
	$(CC) $(LDFLAGS) -o evos \
		$(CLASSES) \
		$(OBJ)/evos.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/evaluate.o:\
		evaluate.C
	$(CC) $(CFLAGS) -o $(OBJ)/evaluate.o -c \
		evaluate.C
#---------------------------------------------------------
evaluate: \
		$(OBJ)/evaluate.o
	$(CC) $(LDFLAGS) -o evaluate \
		$(OBJ)/evaluate.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/matrix-entropy.o:\
		matrix-entropy.C
	$(CC) $(CFLAGS) -o $(OBJ)/matrix-entropy.o -c \
		matrix-entropy.C
#---------------------------------------------------------
matrix-entropy: \
		$(OBJ)/matrix-entropy.o
	$(CC) $(LDFLAGS) -o matrix-entropy \
		$(OBJ)/matrix-entropy.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/dump-model.o:\
		dump-model.C
	$(CC) $(CFLAGS) -o $(OBJ)/dump-model.o -c \
		dump-model.C
#---------------------------------------------------------
dump-model: \
		$(CLASSES) \
		$(OBJ)/dump-model.o
	$(CC) $(LDFLAGS) -o dump-model \
		$(CLASSES) \
		$(OBJ)/dump-model.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/halpern-bruno.o:\
		halpern-bruno.C
	$(CC) $(CFLAGS) -o $(OBJ)/halpern-bruno.o -c \
		halpern-bruno.C
#---------------------------------------------------------
halpern-bruno: \
		$(OBJ)/halpern-bruno.o
	$(CC) $(LDFLAGS) -o halpern-bruno \
		$(OBJ)/halpern-bruno.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/generate-background-seq.o:\
		generate-background-seq.C
	$(CC) $(CFLAGS) -o $(OBJ)/generate-background-seq.o -c \
		generate-background-seq.C
#---------------------------------------------------------
generate-background-seq: \
		$(OBJ)/generate-background-seq.o
	$(CC) $(LDFLAGS) -o generate-background-seq \
		$(OBJ)/generate-background-seq.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/test-mixture-model.o:\
		test-mixture-model.C
	$(CC) $(CFLAGS) -o $(OBJ)/test-mixture-model.o -c \
		test-mixture-model.C
#---------------------------------------------------------
test-mixture-model: \
		$(OBJ)/test-mixture-model.o
	$(CC) $(LDFLAGS) -o test-mixture-model \
		$(OBJ)/test-mixture-model.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/alternate-maf-elements.o:\
		alternate-maf-elements.C
	$(CC) $(CFLAGS) -o $(OBJ)/alternate-maf-elements.o -c \
		alternate-maf-elements.C
#---------------------------------------------------------
alternate-maf-elements: \
		$(OBJ)/alternate-maf-elements.o
	$(CC) $(LDFLAGS) -o alternate-maf-elements \
		$(OBJ)/alternate-maf-elements.o \
		$(LIBS)
#--------------------------------------------------------
restore-nj-root: \
		$(OBJ)/restore-nj-root.o
	$(CC) $(LDFLAGS) -o restore-nj-root \
		$(OBJ)/restore-nj-root.o \
		$(LIBS)
#--------------------------------------------------------



#########################################################
#                     CLASSES
#########################################################
$(OBJ)/ConstrainedBackward.o:\
		ConstrainedBackward.C\
		ConstrainedBackward.H
	$(CC) $(CFLAGS) -o $(OBJ)/ConstrainedBackward.o -c \
		ConstrainedBackward.C
#--------------------------------------------------------
$(OBJ)/PostLeafBackward.o:\
		PostLeafBackward.C\
		PostLeafBackward.H
	$(CC) $(CFLAGS) -o $(OBJ)/PostLeafBackward.o -c \
		PostLeafBackward.C
#--------------------------------------------------------
$(OBJ)/LinkBackward.o:\
		LinkBackward.C\
		LinkBackward.H
	$(CC) $(CFLAGS) -o $(OBJ)/LinkBackward.o -c \
		LinkBackward.C
#--------------------------------------------------------
$(OBJ)/PairHMM.o:\
		$(PHMM)/PairHMM.C \
		$(PHMM)/PairHMM.H
	$(CC) $(CFLAGS) -o $(OBJ)/PairHMM.o -c \
		$(PHMM)/PairHMM.C
#--------------------------------------------------------
$(OBJ)/Transducer.o:\
		$(PHMM)/Transducer.C \
		$(PHMM)/Transducer.H
	$(CC) $(CFLAGS) -o $(OBJ)/Transducer.o -c \
		$(PHMM)/Transducer.C
#--------------------------------------------------------
$(OBJ)/PosetBuilder.o:\
		PosetBuilder.C \
		PosetBuilder.H
	$(CC) $(CFLAGS) -o $(OBJ)/PosetBuilder.o -c \
		PosetBuilder.C
#--------------------------------------------------------
$(OBJ)/PosteriorBackward.o:\
		PosteriorBackward.C \
		PosteriorBackward.H
	$(CC) $(CFLAGS) -o $(OBJ)/PosteriorBackward.o -c \
		PosteriorBackward.C
#--------------------------------------------------------
$(OBJ)/State.o:\
		$(PHMM)/State.C \
		$(PHMM)/State.H
	$(CC) $(CFLAGS) -o $(OBJ)/State.o -c \
		$(PHMM)/State.C
#--------------------------------------------------------
$(OBJ)/StatePath.o:\
		$(PHMM)/StatePath.C \
		$(PHMM)/StatePath.H
	$(CC) $(CFLAGS) -o $(OBJ)/StatePath.o -c \
		$(PHMM)/StatePath.C
#---------------------------------------------------------
$(OBJ)/Viterbi.o:\
		$(PHMM)/Viterbi.C \
		$(PHMM)/Viterbi.H
	$(CC) $(CFLAGS) -o $(OBJ)/Viterbi.o -c \
		$(PHMM)/Viterbi.C
#---------------------------------------------------------
$(OBJ)/Sampler.o:\
		$(PHMM)/Sampler.C \
		$(PHMM)/Sampler.H
	$(CC) $(CFLAGS) -o $(OBJ)/Sampler.o -c \
		$(PHMM)/Sampler.C
#---------------------------------------------------------
$(OBJ)/BackwardAlgorithm.o:\
		$(PHMM)/BackwardAlgorithm.C \
		$(PHMM)/BackwardAlgorithm.H
	$(CC) $(CFLAGS) -o $(OBJ)/BackwardAlgorithm.o -c \
		$(PHMM)/BackwardAlgorithm.C
#--------------------------------------------------------
$(OBJ)/Taxon.o:\
		Taxon.C \
		Taxon.H
	$(CC) $(CFLAGS) -o $(OBJ)/Taxon.o -c \
		Taxon.C
#--------------------------------------------------------
$(OBJ)/FunctionalDollo.o:\
		FunctionalDollo.C \
		FunctionalDollo.H
	$(CC) $(CFLAGS) -o $(OBJ)/FunctionalDollo.o -c \
		FunctionalDollo.C
#--------------------------------------------------------
$(OBJ)/FunctionalParse.o:\
		FunctionalParse.C \
		FunctionalParse.H
	$(CC) $(CFLAGS) -o $(OBJ)/FunctionalParse.o -c \
		FunctionalParse.C
#--------------------------------------------------------
$(OBJ)/FunctionalClass.o:\
		FunctionalClass.C \
		FunctionalClass.H
	$(CC) $(CFLAGS) -o $(OBJ)/FunctionalClass.o -c \
		FunctionalClass.C
#--------------------------------------------------------
$(OBJ)/ProfileFelsenstein.o:\
		ProfileFelsenstein.C \
		ProfileFelsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/ProfileFelsenstein.o -c \
		ProfileFelsenstein.C
#--------------------------------------------------------
$(OBJ)/LinkFelsenstein.o:\
		LinkFelsenstein.C \
		LinkFelsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/LinkFelsenstein.o -c \
		LinkFelsenstein.C
#--------------------------------------------------------
$(OBJ)/LossyFelsenstein.o:\
		LossyFelsenstein.C \
		LossyFelsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/LossyFelsenstein.o -c \
		LossyFelsenstein.C
#--------------------------------------------------------
$(OBJ)/GainLossFelsenstein.o:\
		GainLossFelsenstein.C \
		GainLossFelsenstein.H
	$(CC) $(CFLAGS) -o $(OBJ)/GainLossFelsenstein.o -c \
		GainLossFelsenstein.C
#--------------------------------------------------------
$(OBJ)/ProfileViterbi.o:\
		ProfileViterbi.C \
		ProfileViterbi.H
	$(CC) $(CFLAGS) -o $(OBJ)/ProfileViterbi.o -c \
		ProfileViterbi.C
#---------------------------------------------------------
$(OBJ)/ProfileBackward.o:\
		ProfileBackward.C\
		ProfileBackward.H
	$(CC) $(CFLAGS) -o $(OBJ)/ProfileBackward.o -c \
		ProfileBackward.C
#---------------------------------------------------------
$(OBJ)/AlignmentView.o:\
		AlignmentView.C\
		AlignmentView.H
	$(CC) $(CFLAGS) -o $(OBJ)/AlignmentView.o -c \
		AlignmentView.C
#---------------------------------------------------------
$(OBJ)/ProfileSampler.o:\
		ProfileSampler.C\
		ProfileSampler.H
	$(CC) $(CFLAGS) -o $(OBJ)/ProfileSampler.o -c \
		ProfileSampler.C
#---------------------------------------------------------
$(OBJ)/GapPatternAlphabet.o:\
		GapPatternAlphabet.C\
		GapPatternAlphabet.H
	$(CC) $(CFLAGS) -o $(OBJ)/GapPatternAlphabet.o -c \
		GapPatternAlphabet.C
#--------------------------------------------------------
$(OBJ)/GapPattern.o:\
		GapPattern.C\
		GapPattern.H
	$(CC) $(CFLAGS) -o $(OBJ)/GapPattern.o -c \
		GapPattern.C
#--------------------------------------------------------
$(OBJ)/BranchAttributes.o:\
		BranchAttributes.C\
		BranchAttributes.H
	$(CC) $(CFLAGS) -o $(OBJ)/BranchAttributes.o -c \
		BranchAttributes.C
#--------------------------------------------------------
$(OBJ)/ModelCompiler.o:\
		ModelCompiler.C\
		ModelCompiler.H
	$(CC) $(CFLAGS) -o $(OBJ)/ModelCompiler.o -c \
		ModelCompiler.C
#--------------------------------------------------------
$(OBJ)/TransducerTemplate.o:\
		TransducerTemplate.C\
		TransducerTemplate.H
	$(CC) $(CFLAGS) -o $(OBJ)/TransducerTemplate.o -c \
		TransducerTemplate.C
#--------------------------------------------------------
$(OBJ)/TemplateInstantiator.o:\
		TemplateInstantiator.C\
		TemplateInstantiator.H
	$(CC) $(CFLAGS) -o $(OBJ)/TemplateInstantiator.o -c \
		TemplateInstantiator.C
#--------------------------------------------------------
$(OBJ)/BranchHMM.o:\
		BranchHMM.C\
		BranchHMM.H
	$(CC) $(CFLAGS) -o $(OBJ)/BranchHMM.o -c \
		BranchHMM.C
#--------------------------------------------------------
$(OBJ)/Bander.o:\
		Bander.C\
		Bander.H
	$(CC) $(CFLAGS) -o $(OBJ)/Bander.o -c \
		Bander.C
#--------------------------------------------------------
$(OBJ)/HspMatrix.o:\
		HspMatrix.C\
		HspMatrix.H
	$(CC) $(CFLAGS) -o $(OBJ)/HspMatrix.o -c \
		HspMatrix.C
#--------------------------------------------------------
$(OBJ)/AlignIndexMap.o:\
		AlignIndexMap.C\
		AlignIndexMap.H
	$(CC) $(CFLAGS) -o $(OBJ)/AlignIndexMap.o -c \
		AlignIndexMap.C
#---------------------------------------------------------
$(OBJ)/restore-nj-root.o:\
		restore-nj-root.C
	$(CC) $(CFLAGS) -o $(OBJ)/restore-nj-root.o -c \
		restore-nj-root.C
#---------------------------------------------------------
$(OBJ)/LinkForward.o:\
		LinkForward.C\
		LinkForward.H
	$(CC) $(CFLAGS) -o $(OBJ)/LinkForward.o -c \
		LinkForward.C
#---------------------------------------------------------
$(OBJ)/BandedForward.o:\
		BandedForward.C\
		BandedForward.H
	$(CC) $(CFLAGS) -o $(OBJ)/BandedForward.o -c \
		BandedForward.C
#---------------------------------------------------------
$(OBJ)/BandedBackward.o:\
		BandedBackward.C\
		BandedBackward.H
	$(CC) $(CFLAGS) -o $(OBJ)/BandedBackward.o -c \
		BandedBackward.C
#--------------------------------------------------------
$(OBJ)/AlignmentBuilder.o:\
		AlignmentBuilder.C\
		AlignmentBuilder.H
	$(CC) $(CFLAGS) -o $(OBJ)/AlignmentBuilder.o -c \
		AlignmentBuilder.C
#--------------------------------------------------------
$(OBJ)/ResidueAddress.o:\
		ResidueAddress.C\
		ResidueAddress.H
	$(CC) $(CFLAGS) -o $(OBJ)/ResidueAddress.o -c \
		ResidueAddress.C
#--------------------------------------------------------
$(OBJ)/LinkSampler.o:\
		LinkSampler.C\
		LinkSampler.H
	$(CC) $(CFLAGS) -o $(OBJ)/LinkSampler.o -c \
		LinkSampler.C
#--------------------------------------------------------
$(OBJ)/LinkViterbi.o:\
		LinkViterbi.C\
		LinkViterbi.H
	$(CC) $(CFLAGS) -o $(OBJ)/LinkViterbi.o -c \
		LinkViterbi.C
#--------------------------------------------------------
$(OBJ)/SiblingViterbi.o:\
		SiblingViterbi.C\
		SiblingViterbi.H
	$(CC) $(CFLAGS) -o $(OBJ)/SiblingViterbi.o -c \
		SiblingViterbi.C
#---------------------------------------------------------
$(OBJ)/GainLossType.o:\
		GainLossType.C\
		GainLossType.H
	$(CC) $(CFLAGS) -o $(OBJ)/GainLossType.o -c \
		GainLossType.C
#---------------------------------------------------------
$(OBJ)/Posteriors.o:\
		Posteriors.C\
		Posteriors.H
	$(CC) $(CFLAGS) -o $(OBJ)/Posteriors.o -c \
		Posteriors.C
#---------------------------------------------------------
$(OBJ)/LinkParsimony.o:\
		LinkParsimony.C\
		LinkParsimony.H
	$(CC) $(CFLAGS) -o $(OBJ)/LinkParsimony.o -c \
		LinkParsimony.C
#--------------------------------------------------------
$(OBJ)/FunctionalElement.o:\
		FunctionalElement.C\
		FunctionalElement.H
	$(CC) $(CFLAGS) -o $(OBJ)/FunctionalElement.o -c \
		FunctionalElement.C
#---------------------------------------------------------
$(OBJ)/PrecomputedEmissions.o:\
		PrecomputedEmissions.C\
		PrecomputedEmissions.H
	$(CC) $(CFLAGS) -o $(OBJ)/PrecomputedEmissions.o -c \
		PrecomputedEmissions.C
#---------------------------------------------------------
$(OBJ)/BandingPattern.o:\
		BandingPattern.C\
		BandingPattern.H
	$(CC) $(CFLAGS) -o $(OBJ)/BandingPattern.o -c \
		BandingPattern.C
#---------------------------------------------------------
$(OBJ)/HirschbergFrame.o:\
		HirschbergFrame.C\
		HirschbergFrame.H
	$(CC) $(CFLAGS) -o $(OBJ)/HirschbergFrame.o -c \
		HirschbergFrame.C
#---------------------------------------------------------
$(OBJ)/HirschPosteriors.o:\
		HirschPosteriors.C\
		HirschPosteriors.H
	$(CC) $(CFLAGS) -o $(OBJ)/HirschPosteriors.o -c \
		HirschPosteriors.C
#---------------------------------------------------------
$(OBJ)/Hirschberg.o:\
		Hirschberg.C\
		Hirschberg.H
	$(CC) $(CFLAGS) -o $(OBJ)/Hirschberg.o -c \
		Hirschberg.C
#---------------------------------------------------------
$(OBJ)/HirschForwardMax.o:\
		HirschForwardMax.C\
		HirschForwardMax.H
	$(CC) $(CFLAGS) -o $(OBJ)/HirschForwardMax.o -c \
		HirschForwardMax.C
#---------------------------------------------------------
$(OBJ)/HirschBackwardMax.o:\
		HirschBackwardMax.C\
		HirschBackwardMax.H
	$(CC) $(CFLAGS) -o $(OBJ)/HirschBackwardMax.o -c \
		HirschBackwardMax.C
#---------------------------------------------------------
$(OBJ)/HirschForwardSum.o:\
		HirschForwardSum.C\
		HirschForwardSum.H
	$(CC) $(CFLAGS) -o $(OBJ)/HirschForwardSum.o -c \
		HirschForwardSum.C
#---------------------------------------------------------
$(OBJ)/HirschBackwardSum.o:\
		HirschBackwardSum.C\
		HirschBackwardSum.H
	$(CC) $(CFLAGS) -o $(OBJ)/HirschBackwardSum.o -c \
		HirschBackwardSum.C
#---------------------------------------------------------
$(OBJ)/HirschBackwardSumPair.o:\
		HirschBackwardSumPair.C\
		HirschBackwardSumPair.H
	$(CC) $(CFLAGS) -o $(OBJ)/HirschBackwardSumPair.o -c \
		HirschBackwardSumPair.C
#---------------------------------------------------------
$(OBJ)/HirschPass.o:\
		HirschPass.C\
		HirschPass.H
	$(CC) $(CFLAGS) -o $(OBJ)/HirschPass.o -c \
		HirschPass.C
#---------------------------------------------------------
$(OBJ)/SiblingHirschberg.o:\
		SiblingHirschberg.C\
		SiblingHirschberg.H
	$(CC) $(CFLAGS) -o $(OBJ)/SiblingHirschberg.o -c \
		SiblingHirschberg.C
#---------------------------------------------------------
$(OBJ)/ViterbiInterface.o:\
		ViterbiInterface.C\
		ViterbiInterface.H
	$(CC) $(CFLAGS) -o $(OBJ)/ViterbiInterface.o -c \
		ViterbiInterface.C
#---------------------------------------------------------
$(OBJ)/HirschThread.o:\
		HirschThread.C\
		HirschThread.H
	$(CC) $(CFLAGS) -o $(OBJ)/HirschThread.o -c \
		HirschThread.C
#---------------------------------------------------------
$(OBJ)/eval-alignment.o:\
		eval-alignment.C
	$(CC) $(CFLAGS) -o $(OBJ)/eval-alignment.o -c \
		eval-alignment.C
#---------------------------------------------------------
eval-alignment: \
		$(OBJ)/eval-alignment.o \
		$(CLASSES)
	$(CC) $(LDFLAGS) -o eval-alignment \
		$(OBJ)/eval-alignment.o \
		$(CLASSES) \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/discretize-phastcons.o:\
		discretize-phastcons.C
	$(CC) $(CFLAGS) -o $(OBJ)/discretize-phastcons.o -c \
		discretize-phastcons.C
#---------------------------------------------------------
discretize-phastcons: \
		$(OBJ)/discretize-phastcons.o
	$(CC) $(LDFLAGS) -o discretize-phastcons \
		$(OBJ)/discretize-phastcons.o \
		$(CLASSES) \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/fine-tune.o:\
		fine-tune.C
	$(MPICC) $(CFLAGS) -o $(OBJ)/fine-tune.o -c \
		fine-tune.C
#---------------------------------------------------------
fine-tune: \
		$(OBJ)/fine-tune.o \
		$(CLASSES) \
		$(OBJ)/AVES.o
	$(MPICC) $(LDFLAGS) -o fine-tune \
		$(OBJ)/fine-tune.o \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/build-library.o:\
		build-library.C
	$(MPICC) $(CFLAGS) -o $(OBJ)/build-library.o -c \
		build-library.C
#---------------------------------------------------------
build-library: \
		$(OBJ)/build-library.o \
		$(CLASSES) \
		$(OBJ)/AVES.o
	$(MPICC) $(LDFLAGS) -o build-library \
		$(OBJ)/build-library.o \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/library-info.o:\
		library-info.C
	$(MPICC) $(CFLAGS) -o $(OBJ)/library-info.o -c \
		library-info.C
#---------------------------------------------------------
library-info: \
		$(OBJ)/library-info.o \
		$(OBJ)/AVES.o \
		$(CLASSES)
	$(MPICC) $(LDFLAGS) -o library-info \
		$(OBJ)/library-info.o \
		$(OBJ)/AVES.o \
		$(CLASSES) \
		$(LIBS)
#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/make-phylogibbs-tree.o:\
		make-phylogibbs-tree.C
	$(CC) $(CFLAGS) -o $(OBJ)/make-phylogibbs-tree.o -c \
		make-phylogibbs-tree.C
#---------------------------------------------------------
make-phylogibbs-tree: \
		$(OBJ)/make-phylogibbs-tree.o
	$(CC) $(LDFLAGS) -o make-phylogibbs-tree \
		$(OBJ)/make-phylogibbs-tree.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/shuffle-wmm.o:\
		shuffle-wmm.C
	$(CC) $(CFLAGS) -o $(OBJ)/shuffle-wmm.o -c \
		shuffle-wmm.C
#---------------------------------------------------------
shuffle-wmm: \
		$(CLASSES) \
		$(OBJ)/shuffle-wmm.o
	$(CC) $(LDFLAGS) -o shuffle-wmm \
		$(OBJ)/shuffle-wmm.o \
		$(CLASSES) \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/wmm-entropy.o:\
		wmm-entropy.C
	$(CC) $(CFLAGS) -o $(OBJ)/wmm-entropy.o -c \
		wmm-entropy.C
#---------------------------------------------------------
wmm-entropy: \
		$(CLASSES) \
		$(OBJ)/wmm-entropy.o
	$(CC) $(LDFLAGS) -o wmm-entropy \
		$(OBJ)/wmm-entropy.o \
		$(CLASSES) \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/regress-lambda-mu.o:\
		regress-lambda-mu.C
	$(CC) $(CFLAGS) -o $(OBJ)/regress-lambda-mu.o -c \
		regress-lambda-mu.C
#---------------------------------------------------------
regress-lambda-mu: \
		$(CLASSES) \
		$(OBJ)/regress-lambda-mu.o
	$(CC) $(LDFLAGS) -o regress-lambda-mu \
		$(OBJ)/regress-lambda-mu.o \
		$(CLASSES) \
		$(LIBS)
#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/WigBinary.o:\
		WigBinary.C\
		WigBinary.H
	$(CC) $(CFLAGS) -o $(OBJ)/WigBinary.o -c \
		WigBinary.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/compile-wig.o:\
		compile-wig.C
	$(CC) $(CFLAGS) -o $(OBJ)/compile-wig.o -c \
		compile-wig.C
#---------------------------------------------------------
compile-wig: \
		$(OBJ)/compile-wig.o
	$(CC) $(LDFLAGS) -o compile-wig \
		$(OBJ)/compile-wig.o \
		$(LIBS)
#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/dump-bwig.o:\
		dump-bwig.C
	$(CC) $(CFLAGS) -o $(OBJ)/dump-bwig.o -c \
		dump-bwig.C
#---------------------------------------------------------
dump-bwig: \
		$(OBJ)/dump-bwig.o
	$(CC) $(LDFLAGS) -o dump-bwig \
		$(OBJ)/dump-bwig.o \
		$(LIBS)
#---------------------------------------------

#--------------------------------------------------------
$(OBJ)/roc-wig.o:\
		roc-wig.C
	$(CC) $(CFLAGS) -o $(OBJ)/roc-wig.o -c \
		roc-wig.C
#---------------------------------------------------------
roc-wig: \
		$(OBJ)/roc-wig.o \
		$(OBJ)/WigBinary.o
	$(CC) $(LDFLAGS) -o roc-wig \
		$(OBJ)/roc-wig.o \
		$(OBJ)/WigBinary.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/correlate-wig.o:\
		correlate-wig.C
	$(CC) $(CFLAGS) -o $(OBJ)/correlate-wig.o -c \
		correlate-wig.C
#---------------------------------------------------------
correlate-wig: \
		$(OBJ)/correlate-wig.o \
		$(OBJ)/WigBinary.o
	$(CC) $(LDFLAGS) -o correlate-wig \
		$(OBJ)/correlate-wig.o \
		$(OBJ)/WigBinary.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/normalize-wig.o:\
		normalize-wig.C
	$(CC) $(CFLAGS) -o $(OBJ)/normalize-wig.o -c \
		normalize-wig.C
#---------------------------------------------------------
normalize-wig: \
		$(OBJ)/normalize-wig.o \
		$(OBJ)/WigBinary.o
	$(CC) $(LDFLAGS) -o normalize-wig \
		$(OBJ)/normalize-wig.o \
		$(OBJ)/WigBinary.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/scan-wig.o:\
		scan-wig.C
	$(CC) $(CFLAGS) -o $(OBJ)/scan-wig.o -c \
		scan-wig.C
#---------------------------------------------------------
scan-wig: \
		$(OBJ)/scan-wig.o \
		$(OBJ)/WigBinary.o
	$(CC) $(LDFLAGS) -o scan-wig \
		$(OBJ)/scan-wig.o \
		$(OBJ)/WigBinary.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/find-pileups.o:\
		find-pileups.C
	$(CC) $(CFLAGS) -o $(OBJ)/find-pileups.o -c \
		find-pileups.C
#---------------------------------------------------------
find-pileups: \
		$(OBJ)/find-pileups.o
	$(CC) $(LDFLAGS) -o find-pileups \
		$(OBJ)/find-pileups.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/analyze-conservation.o:\
		analyze-conservation.C
	$(CC) $(CFLAGS) -o $(OBJ)/analyze-conservation.o -c \
		analyze-conservation.C
#---------------------------------------------------------
analyze-conservation: \
		$(OBJ)/analyze-conservation.o \
		$(CLASSES)
	$(CC) $(LDFLAGS) -o analyze-conservation \
		$(CLASSES) \
		$(OBJ)/analyze-conservation.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/BandedFB_Base.o:\
		BandedFB_Base.C\
		BandedFB_Base.H
	$(CC) $(CFLAGS) -o $(OBJ)/BandedFB_Base.o -c \
		BandedFB_Base.C
#---------------------------------------------------------
$(OBJ)/PosteriorMatrix.o:\
		PosteriorMatrix.C\
		PosteriorMatrix.H
	$(CC) $(CFLAGS) -o $(OBJ)/PosteriorMatrix.o -c \
		PosteriorMatrix.C
#---------------------------------------------------------
$(OBJ)/SparseMatrix3D.o:\
		SparseMatrix3D.C\
		SparseMatrix3D.H
	$(CC) $(CFLAGS) -o $(OBJ)/SparseMatrix3D.o -c \
		SparseMatrix3D.C
#---------------------------------------------------------
$(OBJ)/SparseMatrix2D.o:\
		SparseMatrix2D.C\
		SparseMatrix2D.H
	$(CC) $(CFLAGS) -o $(OBJ)/SparseMatrix2D.o -c \
		SparseMatrix2D.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/blimp.o:\
		blimp.C
	$(MPICC) $(CFLAGS) -o $(OBJ)/blimp.o -c \
		blimp.C
#---------------------------------------------------------
blimp: \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/blimp.o
	$(MPICC) $(LDFLAGS) -o blimp \
		$(OBJ)/blimp.o \
		$(OBJ)/AVES.o \
		$(CLASSES) \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/consistency.o:\
		consistency.C
	$(CC) $(CFLAGS) -o $(OBJ)/consistency.o -c \
		consistency.C
#---------------------------------------------------------
consistency: \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/consistency.o
	$(CC) $(LDFLAGS) -o consistency \
		$(OBJ)/consistency.o \
		$(OBJ)/AVES.o \
		$(CLASSES) \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/marginalize-states.o:\
		marginalize-states.C
	$(CC) $(CFLAGS) -o $(OBJ)/marginalize-states.o -c \
		marginalize-states.C
#---------------------------------------------------------
marginalize-states: \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/marginalize-states.o
	$(CC) $(LDFLAGS) -o marginalize-states \
		$(OBJ)/marginalize-states.o \
		$(OBJ)/AVES.o \
		$(CLASSES) \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/combine-factors.o:\
		combine-factors.C
	$(CC) $(CFLAGS) -o $(OBJ)/combine-factors.o -c \
		combine-factors.C
#---------------------------------------------------------
combine-factors: \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/combine-factors.o
	$(CC) $(LDFLAGS) -o combine-factors \
		$(OBJ)/combine-factors.o \
		$(OBJ)/AVES.o \
		$(CLASSES) \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/sample-alignments.o:\
		sample-alignments.C
	$(MPICC) $(CFLAGS) -o $(OBJ)/sample-alignments.o -c \
		sample-alignments.C
#---------------------------------------------------------
sample-alignments: \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/sample-alignments.o
	$(MPICC) $(LDFLAGS) -o sample-alignments \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/sample-alignments.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/marge.o:\
		marge.C
	$(MPICC) $(CFLAGS) -o $(OBJ)/marge.o -c \
		marge.C
#---------------------------------------------------------
marge: \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/marge.o
	$(MPICC) $(LDFLAGS) -o marge \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/marge.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/sample-annotations.o:\
		sample-annotations.C
	$(MPICC) $(CFLAGS) -o $(OBJ)/sample-annotations.o -c \
		sample-annotations.C
#---------------------------------------------------------
sample-annotations: \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/sample-annotations.o
	$(MPICC) $(LDFLAGS) -o sample-annotations \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/sample-annotations.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/merge-slave-alignments.o:\
		merge-slave-alignments.C
	$(MPICC) $(CFLAGS) -o $(OBJ)/merge-slave-alignments.o -c \
		merge-slave-alignments.C
#---------------------------------------------------------
merge-slave-alignments: \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/merge-slave-alignments.o
	$(MPICC) $(LDFLAGS) -o merge-slave-alignments \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/merge-slave-alignments.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/ResidueOrthologyGraph.o:\
		ResidueOrthologyGraph.C\
		ResidueOrthologyGraph.H
	$(CC) $(CFLAGS) -o $(OBJ)/ResidueOrthologyGraph.o -c \
		ResidueOrthologyGraph.C
#--------------------------------------------------------
$(OBJ)/CollapsedOrthologyMatrix.o:\
		CollapsedOrthologyMatrix.C\
		CollapsedOrthologyMatrix.H
	$(CC) $(CFLAGS) -o $(OBJ)/CollapsedOrthologyMatrix.o -c \
		CollapsedOrthologyMatrix.C
#--------------------------------------------------------
$(OBJ)/correlate-posteriors.o:\
		correlate-posteriors.C
	$(CC) $(CFLAGS) -o $(OBJ)/correlate-posteriors.o -c \
		correlate-posteriors.C
#---------------------------------------------------------
correlate-posteriors: \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/correlate-posteriors.o
	$(CC) $(LDFLAGS) -o correlate-posteriors \
		$(CLASSES) \
		$(OBJ)/AVES.o \
		$(OBJ)/correlate-posteriors.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/compound.o:\
		compound.C
	$(CC) $(CFLAGS) -o $(OBJ)/compound.o -c \
		compound.C
#---------------------------------------------------------
compound: \
		$(CLASSES) \
		$(OBJ)/compound.o
	$(CC) $(LDFLAGS) -o compound \
		$(CLASSES) \
		$(OBJ)/compound.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/BirthDeathMatrix.o:\
		BirthDeathMatrix.C\
		BirthDeathMatrix.H
	$(CC) $(CFLAGS) -o $(OBJ)/BirthDeathMatrix.o -c \
		BirthDeathMatrix.C
#---------------------------------------------------------
$(OBJ)/test-gainloss-matrix.o:\
		test-gainloss-matrix.C
	$(CC) $(CFLAGS) -o $(OBJ)/test-gainloss-matrix.o -c \
		test-gainloss-matrix.C
#---------------------------------------------------------
test-gainloss-matrix: \
		$(CLASSES) \
		$(OBJ)/test-gainloss-matrix.o
	$(CC) $(LDFLAGS) -o test-gainloss-matrix \
		$(CLASSES) \
		$(OBJ)/test-gainloss-matrix.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/GainLossEvent.o:\
		GainLossEvent.C\
		GainLossEvent.H
	$(CC) $(CFLAGS) -o $(OBJ)/GainLossEvent.o -c \
		GainLossEvent.C
#---------------------------------------------------------
$(OBJ)/sample-wmm.o:\
		sample-wmm.C
	$(CC) $(CFLAGS) -o $(OBJ)/sample-wmm.o -c \
		sample-wmm.C
#---------------------------------------------------------
sample-wmm: \
		$(CLASSES) \
		$(OBJ)/sample-wmm.o
	$(CC) $(LDFLAGS) -o sample-wmm \
		$(CLASSES) \
		$(OBJ)/sample-wmm.o \
		$(LIBS)
#---------------------------------------------
