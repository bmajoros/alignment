/****************************************************************
 eval-alignment.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/MultSeqAlignment.H"
#include "BOOM/DnaDashDotAlphabet.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;

extern Alphabet alphabet;
Alphabet alphabet;

/****************************************************************
                          class Edge
 ****************************************************************/
struct Edge {
  int fromIndex;
  int toIndex;
  int toSpecies;
  Edge(int from,int to,int toS) 
    : fromIndex(from), toIndex(to), toSpecies(toS) {}
  bool operator==(const Edge &other) {
    return other.toSpecies==toSpecies &&
      other.fromIndex==fromIndex &&
      other.toIndex==toIndex;
  }
  bool constrainedEquality(const Edge &other) {
    return other.toSpecies==toSpecies && other.toIndex==toIndex;
  }
};



/****************************************************************
                        class EdgeArray
 ****************************************************************/
class EdgeArray {
  Array1D< List<Edge> > edges;
  int genomicLength;
public:
  void resize(int maxLength) { edges.resize(maxLength); }
  List<Edge> &operator[](int i) {return edges[i];}
  void setLength(int L) {genomicLength=L;}
  int getLength() {return genomicLength;}
};



/****************************************************************
                        class MafGraph
 ****************************************************************/
class MafGraph {
  Array1D<EdgeArray> tracks;
public:
  MafGraph(int numTracks,int length) 
    : tracks(numTracks)
  { for(int i=0 ; i<numTracks ; ++i) tracks[i].resize(length); }
  EdgeArray &operator[](int i) {return tracks[i];}
  int getNumTracks() const {return tracks.size();}
};



/****************************************************************
                          Application
 ****************************************************************/
class Application {
public:
  Application();
  int main(int argc,char *argv[]);
protected:
  Symbol gapSymbol;
  MultSeqAlignment *loadAlignment(const String &filename);
  MafGraph *computeGraph(MultSeqAlignment &);
  void score(MafGraph &correctGraph,MafGraph &predictedGraph);
  bool findEdge(const Edge &,List<Edge> &);
};



/****************************************************************
                            main()
 ****************************************************************/
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



/****************************************************************
                       Application methods
 ****************************************************************/
Application::Application()
{
  ::alphabet=DnaDashDotAlphabet::global();
  gapSymbol=alphabet.lookup('-');
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("eval-alignment <correct.maf> <predicted.maf>");
  String infileCorrect=cmd.arg(0);
  String infilePredicted=cmd.arg(1);

  // Load alignments
  MultSeqAlignment *correctAln=loadAlignment(infileCorrect);
  MultSeqAlignment *predictedAln=loadAlignment(infilePredicted);

  // Collect edges
  MafGraph *correctGraph=computeGraph(*correctAln);
  MafGraph *predictedGraph=computeGraph(*predictedAln);

  // Compute scores
  score(*correctGraph,*predictedGraph);
  
  return 0;
}



MultSeqAlignment *Application::loadAlignment(const String &filename)
{
  ifstream is(filename.c_str());
  if(!is.good()) throw String("Error opening file: ")+filename;
  MultiAlignment *maf=MultiAlignment::nextAlignmentFromMAF(is);
  MultSeqAlignment *msa=new MultSeqAlignment(*maf,alphabet,gapSymbol);
  is.close();
  delete maf;
  return msa;
}



MafGraph *Application::computeGraph(MultSeqAlignment &maf)
{
  maf.sortTracksByName();
  int numTracks=maf.getNumTracks();
  int numTracksMinus1=numTracks-1;
  int L=maf.getLength();
  MafGraph &graph=*new MafGraph(numTracks,L);
  for(int from=0 ; from<numTracksMinus1 ; ++from) {
    AlignmentSeq &fromTrack=maf.getIthTrack(from);
    for(int to=from+1 ; to<numTracks ; ++to) {
      AlignmentSeq &toTrack=maf.getIthTrack(to);
      int fromGenPos=0, toGenPos=0;
      for(int col=0 ; col<L ; ++col) {
	Symbol a=fromTrack[col], b=toTrack[col];
	if(a!=gapSymbol && b!=gapSymbol) {
	  Edge edge(fromGenPos,toGenPos,to);
	  graph[from][fromGenPos].push_back(edge);
	}
	if(a!=gapSymbol) ++fromGenPos;
	if(b!=gapSymbol) ++toGenPos;
      }
      graph[from].setLength(fromGenPos);
      graph[to].setLength(toGenPos);
    }
  }
  return &graph;
}



void Application::score(MafGraph &correctGraph,MafGraph &predictedGraph)
{
  int numTracks=correctGraph.getNumTracks();
  int numTracksMinus1=numTracks-1;
  int TP=0, FP=0, FN=0;
  for(int from=0 ; from<numTracksMinus1 ; ++from) {
    EdgeArray &predFromArray=predictedGraph[from];
    EdgeArray &knownFromArray=correctGraph[from];
    int L=predFromArray.getLength();
    for(int pos=0 ; pos<L ; ++pos) {
      List<Edge> &predEdges=predFromArray[pos];
      List<Edge> &knownEdges=knownFromArray[pos];
      List<Edge>::iterator cur=predEdges.begin(), end=predEdges.end();
      int found=0;
      for(; cur!=end ; ++cur)
	if(findEdge(*cur,knownEdges)) ++found;
      TP+=found;
      FP+=predEdges.size()-found;
      FN+=knownEdges.size()-found;
    }
  }
  double Sn=TP+FN>0 ? int(1000*TP/double(TP+FN)+5/9.0)/10.0: 0.0;
  double Sp=TP+FP>0 ? int(1000*TP/double(TP+FP)+5/9.0)/10.0 : 0.0;
  double F=2*Sn*Sp/(Sn+Sp);
  if(isNaN(F) || isInfinity(F)) F=0.0;
  cout<<"F="<<F<<" Sn="<<Sn<<" Sp="<<Sp<<endl;
}



bool Application::findEdge(const Edge &edge,List<Edge> &edges)
{
  List<Edge>::iterator cur=edges.begin(), end=edges.end();
  for(; cur!=end ; ++cur)
    if(edge.constrainedEquality(*cur)) return true;
  return false;
}

