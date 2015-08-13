/****************************************************************
 evaluate.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
#include "BOOM/Map.H"
#include "BOOM/Vector.H"
#include "BOOM/Array1D.H"
#include "BOOM/Set.H"
#include "BOOM/FastaReader.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


const int BACKGROUND=0;

enum UniquenessType {
  UT_NONE,
  UT_MIN,
  UT_MAX
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  UniquenessType uniquenessType;
  Map<String, Vector<GffFeature*> > known, predicted;
  Set<String> substrates;
  Map<String,int> substrateLengths;
  Map<String,int> types;
  int gTP, gTN, gFP, gFN;
  bool strandAgnostic;
  int loadGff(const String &filename,
	      Map<String, Vector<GffFeature*> > &,
	      bool allowUniquify,const String &targetSpecies);
  void uniquify(Vector<GffFeature*> &);
  void evalNuc(const String &substrate,Vector<GffFeature*> &known,
	       Vector<GffFeature*> &predicted);
  void evalSite(const String &substrate,Vector<GffFeature*> &known,
		Vector<GffFeature*> &predicted);
  void getSubstrateLengths(const String &filename);
  void render(Vector<GffFeature*> &features,Array1D<int> &anno,
	      int substrateLength);
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
  : gTP(0), gFP(0), gTN(0), gFN(0), uniquenessType(UT_NONE)
{
  types["<background>"]=BACKGROUND;
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"uUS");
  if(cmd.numArgs()!=4)
    throw String("\n\
evaluate [options] <in.fasta> <target-species> <known.gff> <predicted.gff>\n\
    where options=\n\
      -u : uniqify overlapping sites (keeping the minimum scoring site)\n\
      -U : uniqify overlapping sites (keeping the maximum scoring site)\n\
      -S : strand agnostic\n\
");
  String seqFile=cmd.arg(0);
  String targetSpecies=cmd.arg(1);
  String correctGffFile=cmd.arg(2);
  String predictedGffFile=cmd.arg(3);
  if(cmd.option('u')) uniquenessType=UT_MIN;
  if(cmd.option('U')) uniquenessType=UT_MAX;
  strandAgnostic=cmd.option('S');
  
  // Load GFF
  getSubstrateLengths(seqFile);
  substrateLengths.getKeys(substrates);
  int numSites=loadGff(correctGffFile,known,false,targetSpecies);
  int numPredictedSites=loadGff(predictedGffFile,predicted,true,targetSpecies);
  cout<<numPredictedSites<<" predicted sites, "<<numSites
      <<" actual sites"<<endl;
  //known.getKeys(substrates);
  //predicted.getKeys(substrates);
  
  // Process each substrate (species)
  Set<String>::iterator cur=substrates.begin(), end=substrates.end();
  for(; cur!=end ; ++cur) {
    String substrate=*cur;
    if(substrate!=targetSpecies) continue;
    {//if(known.isDefined(substrate) && predicted.isDefined(substrate)) {
      Vector<GffFeature*> knownFeatures=known[substrate];
      Vector<GffFeature*> predictedFeatures=predicted[substrate];
      evalNuc(substrate,knownFeatures,predictedFeatures);
      evalSite(substrate,knownFeatures,predictedFeatures);
    }
  }

  /*
  double Sn=int(1000*gTP/double(gTP+gFN)+5/9.0)/10.0;
  double Sp=int(1000*gTP/double(gTP+gFP)+5/9.0)/10.0;
  double F=2*Sn*Sp/(Sn+Sp);
  cout<<"OVERALL: F="<<F<<" Sn="<<Sn<<"% Sp="<<Sp<<"%"<<endl;
  */

  return 0;
}



int Application::loadGff(const String &filename,
			 Map<String,Vector< GffFeature*> > &M,
			 bool allowUniquify,const String &targetSpecies) 
{
  GffReader reader(filename);
  Vector<GffFeature*> *features=reader.loadFeatures();
  Vector<GffFeature*>::iterator cur=features->begin(), end=features->end();
  Vector<GffFeature*> temp;
  for(; cur!=end ; ++cur) {
    GffFeature *feature=*cur;
    feature->setBegin(feature->getBegin()+1); // ###
    if(feature->getSubstrate()==targetSpecies) temp.push_back(feature);
    else delete feature;
  }
  *features=temp;
  if(uniquenessType!=UT_NONE && allowUniquify) uniquify(*features);
  cur=features->begin();
  end=features->end();
  int n=0;
  for(; cur!=end ; ++cur) {
    GffFeature *feature=*cur;
    const String &substrate=feature->getSubstrate();
    if(!substrates.isMember(substrate)) continue;
    if(!M.isDefined(substrate)) M[substrate]=Vector<GffFeature*>();
    String ftype=feature->getFeatureType();
    if(ftype[0]=='-') {
      ftype=ftype.substring(1,ftype.length()-1);
      feature->setFeatureType(ftype);
      feature->setStrand('-');
    }
    M[substrate].push_back(feature);
    ++n;
  }
  delete features;
  return n;
}



void Application::uniquify(Vector<GffFeature*> &features)
{
  // Precondition: Application::uniquenessType!=UT_NONE

  Set<GffFeature*> blacklist;
  int n=features.size();
  for(int i=0 ; i<n ; ++i) {
    GffFeature *a=features[i];
    double aScore=a->getScore();
    for(int j=0 ; j<n ; ++j) {
      if(j==i) continue;
      GffFeature *b=features[j];
      double bScore=b->getScore();
      if(a->overlapBases(*b)>0) {
	//cout<<"overlap: "<<a->overlapBases(*b)<<" bp"<<endl;//###
	switch(uniquenessType) 
	  {
	  case UT_MIN: 
	    if(bScore<aScore) blacklist.insert(a);
	    //j=n;
	    break;
	  case UT_MAX:
	    if(bScore<aScore) blacklist.insert(b);
	    //j=n;
	    break;
	  }	
      }
    }
  }
  Vector<GffFeature*> temp;
  for(int i=0 ; i<n ; ++i) {
    GffFeature *feature=features[i];
    if(blacklist.isMember(feature))
      delete feature;
    else
      temp.push_back(feature);
  }
  features=temp;
}



void Application::getSubstrateLengths(const String &filename)
{
  FastaReader reader(filename);
  String defline, sequence, substrate, junk;
  while(reader.nextSequence(defline,sequence)) {
    reader.parseDefline(defline,substrate,junk);
    substrateLengths[substrate]=sequence.length();
  }
}



void Application::render(Vector<GffFeature*> &features,Array1D<int> &anno,
			 int substrateLen)
{
  anno.setAllTo(BACKGROUND);
  Vector<GffFeature*>::iterator cur=features.begin(), end=features.end();
  for(; cur!=end ; ++cur) {
    GffFeature *feature=*cur;
    String type=feature->getFeatureType();
    type+=feature->getStrand();
    if(!types.isDefined(type)) types[type]=types.size();
    int typeID=types[type];
    int from=feature->getBegin()/*+1*/, to=feature->getEnd()+1;//###
    if(from<0 || to<0 || from>=substrateLen || to>substrateLen)
      throw String("Invalid coordinates: (")+from+","+to+"), L="+substrateLen;
    for(int i=from ; i<to ; ++i) anno[i]=typeID;
    //cout<<type<<" ["<<from<<","<<to<<")"<<endl;
  }
}



void Application::evalNuc(const String &substrate,Vector<GffFeature*>
			  &known,Vector<GffFeature*> &predicted)
{
  int L=substrateLengths[substrate];
  if(L==0) return;
  Array1D<int> knownAnno(L), predictedAnno(L);
  render(known,knownAnno,L);
  render(predicted,predictedAnno,L);
  int FP=0, FN=0, TP=0, TN=0, wrongFactor=0;
  for(int i=0 ; i<L ; ++i) {
    int o=predictedAnno[i], e=knownAnno[i];
    if(e==BACKGROUND)
      if(o==BACKGROUND) ++TN;
      else ++FP;
    else if(o==BACKGROUND) ++FN;
    else {
      ++TP;
      //if(o==e) ++TP; else ++FP;
      if(o!=e) ++wrongFactor;
    }
  }
  gTP+=TP;
  gFP+=FP;
  gTN+=TN;
  gFN+=FN;
  double Sn=TP+FN>0 ? int(1000*TP/double(TP+FN)+5/9.0)/10.0 : 0.0;
  double Sp=TP+FP>0 ? int(1000*TP/double(TP+FP)+5/9.0)/10.0 : 0.0;
  double F=2*Sn*Sp/(Sn+Sp);
  if(isNaN(F) || isInfinity(F)) F=0.0;
  double foregroundConfusion=
    TP ? int(1000*wrongFactor/double(TP)+5/9.0)/10.0 : 0;
  cout<<"nuc: TP="<<TP<<" TN="<<TN<<" FP="<<FP<<" FN="<<FN<<" F="<<F
      <<" Sn="<<Sn<<" Sp="<<Sp<<" fg_conf="<<foregroundConfusion<<endl;
}



void Application::evalSite(const String &substrate,Vector<GffFeature*>
			  &known,Vector<GffFeature*> &predicted)
{
  // Count FP, FN, TP
  if(substrateLengths[substrate]==0) return;
  int numKnown=known.size(), numPredicted=predicted.size();
  int FP=0, FN=0, TP=0, wrongFactor=0;
  Array1D<bool> seen(numKnown);
  seen.setAllTo(false);
  for(int i=0 ; i<numPredicted ; ++i) {
    GffFeature *pFeature=predicted[i];
    const String &pFeatureType=pFeature->getFeatureType();
    char pFeatureStrand=pFeature->getStrand();
    const float L=pFeature->length();
    bool found=false;
    for(int j=0 ; j<numKnown ; ++j) {
      GffFeature *kFeature=known[j];
      if(kFeature->getEnd()<pFeature->getBegin()) continue;
      if(kFeature->getBegin()>pFeature->getEnd()) break;
      if(kFeature->getFeatureType()!=pFeatureType) continue;
      if(kFeature->getStrand()!=pFeatureStrand && !strandAgnostic) continue;
      int overlap=pFeature->overlapBases(*kFeature);
      if(overlap/L>=0.5) {
	found=true;
	seen[j]=true;
	break;
      }
    }
    ++(found ? TP : FP);
  }

  // compute FN
  for(int i=0 ; i<numKnown ; ++i)
    if(!seen[i]) ++FN;

  // Compute Sn, Sp, F
  double Sn=TP+FN>0 ? int(1000*TP/double(TP+FN)+5/9.0)/10.0: 0.0;
  double Sp=TP+FP>0 ? int(1000*TP/double(TP+FP)+5/9.0)/10.0 : 0.0;
  double F=2*Sn*Sp/(Sn+Sp);
  if(isNaN(F) || isInfinity(F)) F=0.0;
  double foregroundConfusion=
    TP ? int(1000*wrongFactor/double(TP)+5/9.0)/10.0 : 0;
  cout<<"site: TP="<<TP<<" FP="<<FP<<" FN="<<FN<<" F="<<F
      <<" Sn="<<Sn<<" Sp="<<Sp<<endl;
}


