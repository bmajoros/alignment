/****************************************************************
 normalize-wig.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Regex.H"
#include "BOOM/Set.H"
#include "BOOM/Map.H"
#include "BOOM/Constants.H"
#include "BOOM/GffReader.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Array1D.H"
#include "WigBinary.H"
using namespace std;
using namespace BOOM;


struct RocPoint {
  int TP, FP, TN, FN;
  float Sn, Sp, FDR;
  RocPoint() : TP(0), FP(0), TN(0), FN(0), Sn(NEGATIVE_INFINITY), 
	       Sp(NEGATIVE_INFINITY) {}
};


class Application {
  int numPoints;
  Regex numeric;
  Regex wigFilenameRegex; // chr-factor.bwig
  Set<String> factors, chromosomes;
  Vector<GffFeature*> regions; // substrate=chr, type=chunk# (=regionID)
  Map<String,Vector<GffFeature*> > predictions; 
                               // factor x region# -> vect<feature>
  Map<String,WigBinary*> wigFiles; // chr x factor -> WigBinary*
  Map<String,float> minima, maxima; // factor -> extrema
  Array1D<RocPoint> ROC;
  int TP, TN, FP, FN;
  bool useFixedRanges;
  float fixedMax;
  void loadFactorList(const String &filename);
  int loadPredictions(const String &predictionsFilestem,
		      const String &filename,
		      const String &targetSpecies);
  int loadPredictions(const String &filename,int regionID,
		      const String &targetSpecies);
  void loadRegions(const String &filename);
  void attachWigFiles(const String &wigDir);
  void computeROC();
  void computeSnSp();
  void updateCounts(Vector<GffFeature*> &predictions,
		    Vector<WigInterval> &sites,RocPoint &);
  void updateCounts(Vector<GffFeature*> &predictions,
		    WigBinary *wig,const String &factor,
		    GffFeature *chunk);
  GffFeature *findRegion(int regionID);
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
  : numeric("\\d+"), wigFilenameRegex("(\\S+)-(\\S+)\\.bwig"),
    numPoints(100), useFixedRanges(false)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4)
    throw String(
          "normalize-wig <wig-dir> <output-dir> <bg-filestem> <counts.txt>");
  String wigDir=cmd.arg(0);
  String outputDir=cmd.arg(1);
  String bgFilestem=cmd.arg(2);
  String countsFilename=cmd.arg(3);

  Map<String,int> counts;
  File countsFile(countsFilename);
  while(!countsFile.eof()) {
    String line=countsFile.getline();
    Vector<String> *fields=line.getFields();
    if(fields->size()==2) {
      const String &factor=(*fields)[0], &count=(*fields)[1];
      counts[factor]=count.asInt();
    }
    delete fields;
  }

  Vector<String> wigFilenames;
  if(!File::getFileList(wigDir,wigFilenames)) 
    throw String("Can't get file listing for ")+wigDir;
  int nFiles=wigFilenames.size();
  float smallestMin=POSITIVE_INFINITY, largestMax=NEGATIVE_INFINITY;
  for(int i=0 ; i<nFiles ; ++i) {
    const String &filename=wigFilenames[i];
    if(!wigFilenameRegex.match(filename)) continue;
    String chr=wigFilenameRegex[1], factor=wigFilenameRegex[2];
    float totalReads=counts[factor]/1000000;
    float totalBgReads=counts[bgFilestem]/1000000;
    float readsTerm=totalBgReads/totalReads;
    WigBinary wig(wigDir+"/"+filename);
    String bgFilename=wigDir+"/"+chr+"-"+bgFilestem+".bwig";
    WigBinary bgWig(bgFilename);
    int L=wig.getLength(), bgL=bgWig.getLength();
    String outName=outputDir+"/"+filename;
    File outFile(outName,"w");
    for(int i=0 ; i<L ; ++i) {
      float fg=wig.read(i), bg=i<bgL ? bgWig.read(i) : 0.0;
      float newValue=fg/totalReads-bg/totalBgReads;
      if(newValue<0) newValue=0;
      outFile.write(newValue);
    }
  }  
  
  return 0;
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
    //throw filename+String(": BWIG files must be named: chrXX-factor.bwig");
    String chr=wigFilenameRegex[1], factor=wigFilenameRegex[2];
    String key=chr+" "+factor;
    WigBinary *wig=new WigBinary(wigDir+"/"+filename);
    wigFiles[key]=wig;
    float thisMin, thisMax;

    if(useFixedRanges) {
      minima[factor]=0;
      maxima[factor]=fixedMax;
    }
    else {
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
    }
  }

  // ###
  //cout<<"min="<<smallestMin<<" max="<<largestMax<<endl;
  //exit(0);
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



int Application::loadPredictions(const String &predictionsFilestem,
				 const String &chunksDir,
				 const String &targetSpecies) 
{
  Vector<String> subdirs;
  if(!File::getFileList(chunksDir,subdirs)) 
    throw String("Can't get file listing for ")+chunksDir;
  int nDirs=subdirs.size();
  for(int i=0 ; i<nDirs ; ++i) {
    const String &subdir=subdirs[i];
    if(!numeric.match(subdir)) continue;
    int regionID=subdir.asInt();
    String filename=chunksDir+"/"+subdir+"/"+predictionsFilestem+".gff";
    loadPredictions(filename,regionID,targetSpecies);
  }
}


GffFeature *Application::findRegion(int regionID)
{
  Vector<GffFeature*>::iterator cur=regions.begin(), end=regions.end();
  for(; cur!=end ; ++cur) {
    GffFeature *region=*cur;
    if(region->getFeatureType().asInt()==regionID) return region;
  }	
  INTERNAL_ERROR;
}



int Application::loadPredictions(const String &filename,
				 int regionID,
				 const String &targetSpecies)
{
  GffFeature *region=findRegion(regionID);
  int regionOffset=region->getBegin();
  GffReader reader(filename);
  Vector<GffFeature*> *features=reader.loadFeatures();
  Vector<GffFeature*>::iterator cur=features->begin(), end=features->end();
  int n=0;
  for(; cur!=end ; ++cur) {
    GffFeature *feature=*cur;
    feature->setBegin(regionOffset+feature->getBegin()+1); // ###
    feature->setEnd(regionOffset+feature->getEnd());
    if(feature->getSubstrate()!=targetSpecies) {
      delete feature;
      continue;
    }
    const String &substrate=feature->getSubstrate();
    String ftype=feature->getFeatureType();
    ftype.toupper(); // ###
    if(ftype[0]=='-') {
      ftype=ftype.substring(1,ftype.length()-1);
      feature->setFeatureType(ftype);
      feature->setStrand('-');
    }
    String key=ftype+" "+regionID;
    predictions[key].push_back(feature);
    ++n;
  }
  delete features;
  return n;
}



/*
void Application::evalSite(const String &substrate,Vector<GffFeature*>
			  &known,Vector<GffFeature*> &predicted)
{
  // Count FP, FN, TP
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
*/



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



void Application::computeROC()
{
  int numRegions=regions.size();
  for(int i=0 ; i<numRegions ; ++i) {
    GffFeature *region=regions[i];
    String chr=region->getSubstrate();
    int regionID=region->getFeatureType().asInt();
    Set<String>::iterator cur=factors.begin(), end=factors.end();
    for(; cur!=end ; ++cur) {
      const String &factor=*cur;
      String key=factor+" "+regionID;
      //if(!predictions.isDefined(key)) INTERNAL_ERROR;
      Vector<GffFeature*> &preds=predictions[key];
      key=chr+" "+factor;

      if(factor=="HB") key=chr+" HB1"; // ### HACK
      if(factor=="KR") key=chr+" KR1"; // ### HACK

      WigBinary *wig=wigFiles.isDefined(key) ? wigFiles[key] : NULL;
      if(!wig) continue;
      updateCounts(preds,wig,factor,region);
    }
  }
  computeSnSp();
}



void Application::updateCounts(Vector<GffFeature*> &predictions,
			       WigBinary *wig,const String &factor,
			       GffFeature *chunk)
{
  int chunkBegin=chunk->getBegin(), chunkEnd=chunk->getEnd();
  float min=minima[factor], max=maxima[factor];
  float delta=(max-min)/(numPoints-1);
  Vector<WigInterval> boundingRegions;
  boundingRegions.push_back(WigInterval(chunkBegin,chunkEnd));
  for(int i=0 ; i<numPoints ; ++i) {
    float threshold=min+i*delta;
    Vector<WigInterval> sites;
    //cout<<"  asking for regions above "<<threshold<<" in ("<<chunkBegin<<","<<chunkEnd<<")"<<Endl;
    wig->regionsAbove(threshold,sites,boundingRegions);
    //int n=sites.size();  for(int i=0 ; i<n ; ++i) cout<<"      "<<sites[i].begin<<"-"<<sites[i].end<<endl;
    //cout<<"   "<<sites.size()<<" found"<<endl;
    RocPoint &rocPoint=ROC[i];
    updateCounts(predictions,sites,rocPoint);
    boundingRegions=sites;
  }
}



void Application::updateCounts(Vector<GffFeature*> &predictions,
			       Vector<WigInterval> &sites,
			       RocPoint &rocPoint)
{
  // Count TP & FP
  int numKnown=sites.size(), numPredicted=predictions.size();
  int &FP=rocPoint.FP, &TP=rocPoint.TP, &FN=rocPoint.FN;
  Array1D<bool> seen(numKnown);
  seen.setAllTo(false);
  for(int i=0 ; i<numPredicted ; ++i) {
    GffFeature *pFeature=predictions[i];
    const float L=pFeature->length();
    bool found=false;
    for(int j=0 ; j<numKnown ; ++j) {
      WigInterval &kFeature=sites[j];
      if(kFeature.end<pFeature->getBegin()) continue;
      if(kFeature.begin>pFeature->getEnd()) break;
      found=true;
      seen[j]=true;
      break;
    }
    ++(found ? TP : FP);
  }

  // compute FN
  for(int i=0 ; i<numKnown ; ++i)
    if(!seen[i]) ++FN;
}



void Application::computeSnSp()
{
  for(int i=0 ; i<numPoints ; ++i) {
    RocPoint &pt=ROC[i];
    pt.Sn=pt.TP+pt.FN>0 ? pt.TP/float(pt.TP+pt.FN) : 0.0;
    pt.Sp=pt.TP+pt.FP>0 ? pt.TP/float(pt.TP+pt.FP) : 0.0;
    pt.FDR=pt.FP+pt.TP>0 ? pt.FP/float(pt.TP+pt.FP) : 0.0;
    //cout<<pt.Sn<<"\t"<<pt.Sp<<endl;
    cout<<pt.Sn<<"\t"<<pt.FDR<<endl;
  }
}






