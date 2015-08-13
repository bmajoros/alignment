/****************************************************************
 shuffle-wmm.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Random.H"
#include "BOOM/Entropy.H"
#include "BOOM/PureDnaAlphabet.H"
#include "EGGS/WMM.H"
using namespace std;
using namespace BOOM;

const int NUM_ITERATIONS=50000;
const float MAX_INF_DIFF=0.4;
BOOM::Alphabet alphabet;

class Application
{
  void shuffle(WMM &);
  float infContent(int pos,BOOM::Array2D<float> &);
  void swap(BOOM::Array2D<float> &,int pos1,int pos2);
public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
  {
    try
      {
	Application app;
	return app.main(argc,argv);
      }
    catch(const char *p)
      {
	cerr << p << endl;
      }
    catch(const string &msg)
      {
	cerr << msg.c_str() << endl;
      }
    catch(const exception &e)
      {
	cerr << "STL exception caught in main:\n" << e.what() << endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



Application::Application()
  {
  }



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"s:");
  if(cmd.numArgs()!=2)
    throw String("shuffle-wmm [-s <RNSEED>] <infile> <outfile>");
  String infile=cmd.arg(0);
  String outfile=cmd.arg(1);
  GarbageIgnorer GC;

  if(cmd.option('s')) SeedRandomizer(cmd.optParm('s').asUnsigned());
  else randomize();

  // Load matrix
  WMM *wmm=dynamic_cast<WMM*>(SignalSensor::load(infile,GC));
  if(!wmm) throw "Error loading matrix file";
  
  // Perform shuffling
  shuffle(*wmm);

  // Save output
  wmm->save(outfile);

  return 0;
}



float Application::infContent(int pos,BOOM::Array2D<float> &M)
{
  Vector<float> row;
  M.getRow(pos,row);
  int n=row.size();
  Vector<double> dRow(n);
  for(int i=0 ; i<n ; ++i) dRow[i]=exp(row[i]);
  return Entropy<double>::informationContent(dRow);
}


/*
void Application::shuffle(WMM &wmm)
{
  BOOM::Array2D<float> &M=wmm.getMatrix(); // indexed as M[pos][symbol]
  int numPositions=M.getFirstDim();
  for(int i=0 ; i<NUM_ITERATIONS ; ++i) {
    int pos1=RandomNumber(numPositions);
    int pos2=(pos1+RandomNumber(numPositions-1))%numPositions;
    while(1) {
      float infDiff=fabs(infContent(pos2,M)-infContent(pos1,M));
      if(infDiff>MAX_INF_DIFF) pos2=(pos2+1)%numPositions;
      else break;
      if(pos2==pos1) {
	throw "Can't swap columns";
      }
    }
    swap(M,pos1,pos2);
  }
}
*/


void Application::shuffle(WMM &wmm)
{
  BOOM::Array2D<float> &M=wmm.getMatrix(); // indexed as M[pos][symbol]
  int numPositions=M.getFirstDim();
  int retries=0;
  Array1D<bool> colsShuffled(numPositions);
  colsShuffled.setAllTo(false);
  for(int i=0 ; i<NUM_ITERATIONS ; ++i) {
    int pos1=RandomNumber(numPositions);
    int pos2=(pos1+RandomNumber(numPositions-1))%numPositions;
    float infDiff=fabs(infContent(pos2,M)-infContent(pos1,M));
    if(infDiff>MAX_INF_DIFF) {
      --i;
      ++retries;
      //cout<<"retrying..."<<endl;
      if(retries>1000) throw "Can't swap columns";
      continue;
    }
    retries=0;
    swap(M,pos1,pos2);
    colsShuffled[pos1]=colsShuffled[pos2]=true;
  }
  int numColsShuffled=0;
  for(int i=0 ; i<numPositions ; ++i) if(colsShuffled[i]) ++numColsShuffled;
  cout<<numColsShuffled<<" columns were shuffled (out of "
      <<numPositions<<")"<<endl;
}



void Application::swap(BOOM::Array2D<float> &M,int pos1,int pos2)
{
  int n=M.getSecondDim();
  BOOM::Array2D<float>::RowIn2DArray<float> row1=M[pos1], row2=M[pos2];
  for(int i=0 ; i<n ; ++i) {
    float tmp=row1[i];
    row1[i]=row2[i];
    row2[i]=tmp;
  }
}


