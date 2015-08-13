CC		= g++
DEBUG		= -g -w -fpermissive
OPTIMIZE	= -O -w -fpermissive
CFLAGS		= $(OPTIMIZE)
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
LIBS		= -LBOOM -lBOOM
#---------------------------------------------------------
$(OBJ)/pair-align.o:\
		pair-align.C
	$(CC) $(CFLAGS) -o $(OBJ)/pair-align.o -c \
		pair-align.C
#---------------------------------------------------------
$(OBJ)/StatePath.o:\
		StatePath.C \
		StatePath.H
	$(CC) $(CFLAGS) -o $(OBJ)/StatePath.o -c \
		StatePath.C
#---------------------------------------------------------
$(OBJ)/PairHMM.o:\
		PairHMM.C \
		PairHMM.H
	$(CC) $(CFLAGS) -o $(OBJ)/PairHMM.o -c \
		PairHMM.C
#---------------------------------------------------------
$(OBJ)/Viterbi.o:\
		Viterbi.C \
		Viterbi.H
	$(CC) $(CFLAGS) -o $(OBJ)/Viterbi.o -c \
		Viterbi.C
#---------------------------------------------------------
$(OBJ)/Sampler.o:\
		Sampler.C \
		Sampler.H
	$(CC) $(CFLAGS) -o $(OBJ)/Sampler.o -c \
		Sampler.C
#---------------------------------------------------------
$(OBJ)/BackwardAlgorithm.o:\
		BackwardAlgorithm.C \
		BackwardAlgorithm.H
	$(CC) $(CFLAGS) -o $(OBJ)/BackwardAlgorithm.o -c \
		BackwardAlgorithm.C
#---------------------------------------------------------
pair-align: \
		$(OBJ)/StatePath.o \
		$(OBJ)/State.o \
		$(OBJ)/Sampler.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/Viterbi.o \
		$(OBJ)/PairHMM.o \
		$(OBJ)/pair-align.o
	$(CC) $(LDFLAGS) -o pair-align \
		$(OBJ)/StatePath.o \
		$(OBJ)/State.o \
		$(OBJ)/Sampler.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/Viterbi.o \
		$(OBJ)/PairHMM.o \
		$(OBJ)/pair-align.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/train-3state.o:\
		train-3state.C
	$(CC) $(CFLAGS) -o $(OBJ)/train-3state.o -c \
		train-3state.C
#---------------------------------------------------------
train-3state: \
		$(OBJ)/PairHMM.o \
		$(OBJ)/train-3state.o
	$(CC) $(LDFLAGS) -o train-3state \
		$(OBJ)/PairHMM.o \
		$(OBJ)/train-3state.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/State.o:\
		State.C\
		State.H
	$(CC) $(CFLAGS) -o $(OBJ)/State.o -c \
		State.C
#---------------------------------------------------------


