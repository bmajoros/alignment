/****************************************************************
 roc-wig.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Regex.H"
#include "BOOM/Set.H"
#include "BOOM/Map.H"
#include "BOOM/Constants.H"
#include "BOOM/GffReader.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Array1D.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/SumLogProbs.H"
#include <gsl/gsl_cdf.h>
#include "WigBinary.H"
using namespace std;
using namespace BOOM;

const int MIN_SAMPLE_SIZE=5;
const int chunkSize=10000;

class Application {
  Regex numeric;
  Regex wigFilenameRegex; // chr-factor.bwig
  Set<String> factors, chromosomes;
  Map<String,WigBinary*> wigFiles; // chr x factor -> WigBinary*
  Map<String,float> minima, maxima; // factor -> extrema
  Vector<GffFeature*> regions; // substrate=chr, type=chunk# (=regionID)
  void allPairs();
  double processPair(const String &factor1,const String &factor2,double &P);
  void loadFactorList(const String &filename);
  void attachWigFiles(const String &wigDir);
  void loadRegions(const String &filename);
  double computeR(double sumX,double sumXX,double sumY,
		  double sumYY,double sumXY,int n);
  double subtractLogProbs(double logA,double logB);
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
  : numeric("\\d+"), wigFilenameRegex("(\\S+)-(\\S+)\\.bwig")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()<4)
    throw String("\n\
roc-wig <factors.txt> <target-species> <chunks-dir> <wig-dir> <regions>\n\
");
  String factorsFile=cmd.arg(0);
  String targetSpecies=cmd.arg(1);
  String chunksDir=cmd.arg(2);
  String wigDir=cmd.arg(3);
  loadRegions(cmd.arg(4));

  loadFactorList(factorsFile);
  attachWigFiles(wigDir);
  allPairs();
  
  return 0;
}



void Application::loadRegions(const String &filename)
{
  GffReader reader(filename);
  Vector<GffFeature*> *features=reader.loadFeatures();
  Vector<GffFeature*>::iterator cur=features->begin(), end=features->end();
  for(; cur!=end ; ++cur) {
    GffFeature *feature=*cur;
    regions.push_back(feature);
    chromosomes.insert(feature->getSubstrate());
  }
  delete features;
}



double Application::processPair(const String &factor1,const String &factor2,
				double &P)
{
  LogProbsSummer<float> sumX(chunkSize), sumXX(chunkSize), sumY(chunkSize), 
    sumYY(chunkSize), sumXY(chunkSize);
  int n=0;
  Vector<GffFeature*>::iterator cur=regions.begin(), end=regions.end();
  //int counter=0;
  for(; cur!=end ; ++cur) {
    GffFeature* feature=*cur;
    String chr=feature->getSubstrate();
    String key1=chr+" "+factor1, key2=chr+" "+factor2;
    WigBinary *wig1=wigFiles[key1], *wig2=wigFiles[key2];
    int begin=feature->getBegin(), end=feature->getEnd();
    for(int pos=begin ; pos<end ; ++pos) {
      float x=wig1->read(pos), y=wig2->read(pos);
      //if(factor1=="CAD" && factor2=="KR") cout<<x<<" "<<y<<endl;
      //cout<<x<<" "<<y<<endl;
      //++counter;
      //if(counter%100) continue;
      x=log(x); y=log(y);
      sumX.add(x);
      sumXX.add(x+x);
      sumY.add(y);
      sumYY.add(y+y);
      sumXY.add(x+y);
      ++n;
    }
  }
  double r=computeR(sumX.getSum(),sumXX.getSum(),sumY.getSum(),
		    sumYY.getSum(),sumXY.getSum(),n);
  double t=r*sqrt(double(n-2))/sqrt(double(1-r*r));
  P=gsl_cdf_tdist_Q(t,n-2);
  return r;
}



void Application::allPairs()
{
  Set<String>::iterator cur=factors.begin(), end=factors.end();
  for(; cur!=end ; ++cur) {
    String factor1=*cur;
    Set<String>::iterator other=cur;
    ++other;
    for(; other!=end ; ++other) {
      String factor2=*other;
      double P;
      double r=processPair(factor1,factor2,P);
      //cout<<factor1<<" vs "<<factor2<<" = "<<r<<" P="<<P<<endl;
    }
  }
}



void Application::attachWigFiles(const String &wigDir)
{
  Vector<String> wigFilenames;
  if(!File::getFileList(wigDir,wigFilenames)) 
    throw String("Can't get file listing for ")+wigDir;
  int nFiles=wigFilenames.size();
  float smallestMin=POSITIVE_INFINITY, largestMax=NEGATIVE_INFINITY;
  for(int i=0 ; i<nFiles ; ++i) {
    const String &filename=wigFilenames[i];
    if(!wigFilenameRegex.match(filename)) continue;
    String chr=wigFilenameRegex[1], factor=wigFilenameRegex[2];
    String key=chr+" "+factor;
    WigBinary *wig=new WigBinary(wigDir+"/"+filename);
    wigFiles[key]=wig;
    float thisMin, thisMax;

    /*
    wig->getExtrema(thisMin,thisMax);
    float &min=minima[factor], &max=maxima[factor];
    if(!isFinite(min)) {min=thisMin; max=thisMax;}
    else {
      if(thisMin<min) min=thisMin;
      if(thisMax>max) max=thisMax;
    }
    if(!isFinite(smallestMin)) {
      if(min<smallestMin) smallestMin=min;
      if(max>largestMax) largestMax=max;
    }
    cout<<"smallestMin="<<smallestMin<<" largestMax="<<largestMax<<endl;
    */
  }
}





void Application::loadFactorList(const String &filename)
{
  File f(filename);
  while(!f.eof()) {
    String line=f.getline();
    line.trimWhitespace();
    if(line.length()>0) {
      factors.insert(line);
      minima[line]=POSITIVE_INFINITY;
      maxima[line]=NEGATIVE_INFINITY;
    }
  }
}




double Application::computeR(double sumX,double sumXX,double sumY,
			     double sumYY,double sumXY,int n)
{
  double logN=log((double)n);
  const double yBar=sumY-logN;
  const double xBar=sumX-logN;

  const double Sxx=subtractLogProbs(sumXX,sumX+sumX-logN);
  if(Sxx==0.0) throw String("Error in BOOM::LinRegressor: Sxx is 0");

  const double Syy=subtractLogProbs(sumYY,sumY+sumY-logN);
  const double Sxy=subtractLogProbs(sumXY,sumX+sumY-logN);

  //cout<<Sxx<<" "<<Syy<<" "<<Sxy<<" "<<n<<endl;

  // return Sxy/sqrt(Syy*Sxx)  =  sqrt(Sxy*Sxy / (Syy*Sxx))
  return sqrt(exp(2*Sxy-(Syy+Sxx)));
}



double Application::subtractLogProbs(double logA,double logB)
{
  double smallerValue, largerValue;
  if(logA<logB) {smallerValue=logA; largerValue=logB;}
  else {smallerValue=logB; largerValue=logA;}
  if(!isFinite(logA) && !isFinite(logB)) return logA;
  if(smallerValue==NEGATIVE_INFINITY) return log(-exp(largerValue));

  return largerValue+log(1-exp(smallerValue-largerValue));
}


