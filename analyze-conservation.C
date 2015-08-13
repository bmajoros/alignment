/****************************************************************
 analyze-conservation.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <math.h>
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Regex.H"
#include "BOOM/File.H"
#include "BOOM/MultiAlignment.H"
#include "BOOM/GffReader.H"
#include "BOOM/Exceptions.H"
#include "BOOM/SequenceEntropy.H"
using namespace std;
using namespace BOOM;

static double logOf2=log(2.0);
inline double lg(double x) {return log(x)/logOf2;}

class Application {
  Regex numeric;
  String targetSpecies;
  Vector<GffFeature*> *regions;
  void processFile(const String &mafFile,const String &gffFile,
		   const String &regionID);
  float computeScore(MultiAlignment &,GffFeature *);
  float scoreForeground(MultiAlignment &,int begin,int end);
  float scoreBackground(MultiAlignment &,int begin,int end);
  float occupancy(MultiAlignment &,int begin,int end);
  void toChrCoords(GffFeature &f,const String &regionID);
  GffFeature *findRegion(const String &regionID);
  float information(MultiAlignment &,int begin,int end);
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
  : numeric("\\d+")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=5)
    throw String("analyze-conservation <chunks-dir> <maf-filestem> <gff-filestem> <species> <regions.gff>");
  String chunksDir=cmd.arg(0);
  String mafFilestem=cmd.arg(1);
  String gffFilestem=cmd.arg(2);
  targetSpecies=cmd.arg(3);
  String regionsFilename=cmd.arg(4);

  GffReader reader(regionsFilename);
  regions=reader.loadFeatures();

  Vector<String> subdirs;
  File::getFileList(chunksDir,subdirs);
  int numDirs=subdirs.size();
  for(int i=0 ; i<numDirs ; ++i) {
    String subdir=subdirs[i];
    if(!numeric.match(subdir)) continue;
    String mafFile=chunksDir+"/"+subdir+"/"+mafFilestem+".maf";
    String gffFile=chunksDir+"/"+subdir+"/"+gffFilestem+".gff";
    processFile(mafFile,gffFile,subdir);
  }
  
  return 0;
}



void Application::processFile(const String &mafFile,const String &gffFile,
			      const String &regionID)
{
  MultiAlignment A;
  A.loadMAF(mafFile);
  GffReader reader(gffFile);
  Vector<GffFeature*> &features=*reader.loadFeatures();
  int numFeatures=features.size();
  for(int i=0 ; i<numFeatures ; ++i) {
    GffFeature *f=features[i];
    if(f->getSubstrate()!=targetSpecies) continue;
    float score=computeScore(A,f);
    Vector<String> &extra=f->getExtraFields();
    String field=String("/cons=")+score;
    extra.push_back(field);
    toChrCoords(*f,regionID);
    cout<<*f<<endl;
  }
}



GffFeature *Application::findRegion(const String &regionID)
{
  int n=regions->size();
  for(int i=0 ; i<n ; ++i) {
    GffFeature *region=(*regions)[i];
    if(region->getFeatureType()==regionID)
      return region;
  }
  INTERNAL_ERROR;
}



void Application::toChrCoords(GffFeature &f,const String &regionID)
{
  GffFeature *region=findRegion(regionID);
  int chunkBegin=region->getBegin();
  String chr=region->getSubstrate();
  //cout<<"CHANGING "<<f.getSubstrate()<<" -> "<<regionID<<"="<<chr<<endl;
  f.setSubstrate(chr);
  f.setBegin(f.getBegin()+chunkBegin+1);
  f.setEnd(f.getEnd()+chunkBegin);
}



float Application::computeScore(MultiAlignment &A,GffFeature *f)
{
  int begin=f->getBegin(), end=f->getEnd();
  String species=f->getSubstrate();
  AlignmentTrack &track=A.getTrackByName(species);
  int L=A.getLength();
  int b=track.mapUngappedCoordToGapped(begin);
  int e=b+f->length();
  float fg=scoreForeground(A,b,e);
  float bg=scoreBackground(A,b,e);
  //cout<<"fg="<<fg<<" bg="<<bg<<endl;
  float score=fg-bg; // log(fg)-0.1*log(bg);
  return score;
}



float Application::scoreForeground(MultiAlignment &A,int begin,int end) 
{
  return information(A,begin,end);
  //return occupancy(A,begin,end);
}



float Application::scoreBackground(MultiAlignment &A,int begin,int end) 
{
  int L=A.getLength();
  int margin=10;
  int left=begin-margin, right=end+margin;
  if(left<0) left=0;
  if(right>L) right=L;
  return (information(A,left,begin)+information(A,end,right))/2;
  //return occupancy(A,left,begin)+occupancy(A,end,right);
}



float Application::information(MultiAlignment &A,int begin,int end)
{
  int N=A.getNumTracks();
  float sum=0;
  int sampleSize=0;
  for(int pos=begin ; pos<end ; ++pos) {
    String S;
    for(int i=0 ; i<N ; ++i) {
      AlignmentTrack &track=A.getIthTrack(i);
      if(track.getName()[0]=='A') continue; // skip ancestors
      char c=track[pos];
      S.append(1,c);
    }
    double Hmax;
    double H=SequenceEntropy::entropy(S,Hmax);
    Hmax=-lg(1/5.0);
    //cout<<"\""<<S<<"\" -> "<<H<<" "<<Hmax<<endl;
    float inf=(Hmax-H)/Hmax;
    sum+=inf;
    ++sampleSize;
  }
  return sum/sampleSize;
}



float Application::occupancy(MultiAlignment &A,int begin,int end) 
{
  int residues=0, positions=0;
  int N=A.getNumTracks();
  for(int pos=begin ; pos<end ; ++pos) {
    for(int i=0 ; i<N ; ++i) {
      AlignmentTrack &track=A.getIthTrack(i);
      if(track.getName()[0]=='A') continue; // skip ancestors
      char c=track[pos];
      if(c!='-') ++residues;
      ++positions;
    }
  }
  return residues/float(positions);
}








