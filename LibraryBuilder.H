#ifndef INCL_LibraryBuilder_H
#define INCL_LibraryBuilder_H

#include "AVES.H"

enum MessageTag {
  COMPUTE_TABLE,
  DONE_COMPUTING
};

/****************************************************************
                         class LibraryBuilder
 ****************************************************************/
class LibraryBuilder : public AVES {
  String lambdaFile;
  Vector<STATE> QI, QD, QM;
  int numI, numD, numM;
  int numAlpha;
  String fastaFile, outDir;
  MpiSlaveDriver *slaveDriver;
  String factor;
  double sparsenessThreshold;
  String targetTaxon;
  int targetID;
  bool reestimateEq;
  Array1D<double> bgFreqs;
  void estimateBgEq();
  virtual void initBranches();
  void buildLibrary();
  int master(int argc,char *argv[]);
  int slave(int argc,char *argv[]);
  void slave_serveTheMaster();
  void slave_computeTable(MpiFixedMsg &);
  void processCmdLine(int argc,char *argv[]);
public:
  LibraryBuilder();
  virtual int main(int argc,char *argv[]);
};


#endif
