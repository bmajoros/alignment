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
#include "WigBinary.H"
using namespace std;
using namespace BOOM;

const int MIN_SAMPLE_SIZE=5;

struct SiteComparator : public Comparator<GffFeature*> {
  bool equal(GffFeature* &a,GffFeature* &b)
    {return a->getScore()==b->getScore();}
  bool greater(GffFeature* &a,GffFeature* &b)
    {return a->getScore()>b->getScore();}
  bool less(GffFeature* &a,GffFeature* &b)
    {return a->getScore()<b->getScore();}
};


enum UniquenessType {
  UT_NONE,
  UT_MIN,
  UT_MAX
};

struct RocPoint {
  int TP, FP, TN, FN, N;
  float Sn, Sp, FDR, FPR, threshold;
  RocPoint() : TP(0), FP(0), TN(0), FN(0), Sn(NEGATIVE_INFINITY), 
	       Sp(NEGATIVE_INFINITY) {}
};


class Application {
  SiteComparator cmp;
  UniquenessType uniquenessType;
  int numPoints;
  Regex numeric;
  Regex wigFilenameRegex; // chr-factor.bwig
  Set<String> factors, chromosomes;
  Vector<GffFeature*> regions; // substrate=chr, type=chunk# (=regionID)
  Map<String,Vector<GffFeature*> > predictions;
                               // factor x region# -> vect<feature>
  Map<String,Vector<GffFeature*> > goldStandard;
                               // factor x region# -> vect<feature>
  Map<String,Vector<GffFeature*> > blacklist;
                               // chr x factor -> vect<feature>
  Map<String,WigBinary*> wigFiles; // chr x factor -> WigBinary*
  Map<String,float> minima, maxima; // factor -> extrema
  Array1D<RocPoint> ROC;
  String dumpFile;
  int TP, TN, FP, FN;
  bool useFixedRanges;
  bool useGold;
  float fixedMax;
  int smoothingIterations, smoothingNeighborhood;
  int skipChunk;
  void sort(Vector<GffFeature*> &);
  void loadBlacklist(const String &filename);
  void maskBlacklistSites(Vector<WigInterval> &sites,GffFeature *region,
			  const String &factor);
  void loadFactorList(const String &filename);
  int loadPredictions(const String &predictionsFilestem,
		      const String &filename,
		      const String &targetSpecies,
		      Map<String,Vector<GffFeature*> > &);
  int loadPredictionsRegion(const String &filename,const String &regionID,
			    const String &targetSpecies,
			    Map<String,Vector<GffFeature*> > &);
  void loadRegions(const String &filename);
  void attachWigFiles(const String &wigDir);
  void computeROC();
  void smoothROC();
  void computeSnSp();
  void updateCounts(Vector<GffFeature*> &predictions,
		    Vector<WigInterval> &sites,RocPoint &,int chunkBegin,
		    int chunkEnd);
  void updateCounts(Vector<GffFeature*> &predictions,
		    WigBinary *wig,const String &factor,
		    GffFeature *chunk);
  GffFeature *findRegion(const String &regionID);
  void getNonSites(Vector<WigInterval> &sites,Vector<WigInterval> &nonSites,
		   int chunkBegin,int chunkEnd);
  void uniquify(Vector<GffFeature*> &features);
  void uniquifyByClustering(Vector<GffFeature*> &features);
  void contrast();
  float getScore(GffFeature &site,GffFeature &region);
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
    numPoints(100), useFixedRanges(false), uniquenessType(UT_NONE),
    smoothingIterations(0), smoothingNeighborhood(1),
    useGold(false)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"t:m:uUs:N:d:g:b:S:");
  if(cmd.numArgs()<3)
    throw String("\n\
roc-wig [options] <factors.txt> <predictions-filestem> <target-species> \n\
                  <regions.gff> <chunks-dir> <wig-dir>\n\
   options:\n\
      -t # : use this many threshold points (default=100)\n\
      -m # : use min=0 max=# for threshold ranges\n\
      -u : uniqify overlapping sites (keeping the minimum scoring site)\n\
      -U : uniqify overlapping sites (keeping the maximum scoring site)\n\
      -s # : apply # iterations of smoothing (default=0)\n\
      -N # : smoothing neighborhood (default=1)\n\
      -d file : dump extra curve data into file\n\
      -g filestem : use gold standard GFFs to contrast with predictions\n\
      -b file : mask blacklist intervals from GFF file\n\
      -S N : skip chunk #N\n\
");
  String factorsFile=cmd.arg(0);
  String predictionsFilestem=cmd.arg(1);
  String targetSpecies=cmd.arg(2);
  String regionsFile=cmd.arg(3);
  String chunksDir=cmd.arg(4);
  String wigDir=cmd.arg(5);
  if(cmd.option('u')) uniquenessType=UT_MIN;
  if(cmd.option('U')) uniquenessType=UT_MAX;
  if(cmd.option('t')) numPoints=cmd.optParm('t').asInt();
  ROC.resize(numPoints);
  if(cmd.option('m')) {
    useFixedRanges=true;
    fixedMax=cmd.optParm('m').asFloat();
  }
  if(cmd.option('N')) smoothingNeighborhood=cmd.optParm('N').asInt();
  if(cmd.option('s')) smoothingIterations=cmd.optParm('s').asInt();
  if(cmd.option('d')) dumpFile=cmd.optParm('d');
  if(cmd.option('b')) loadBlacklist(cmd.optParm('b'));
  skipChunk=cmd.option('S') ? cmd.optParm('S').asInt() : LARGEST_INTEGER;

  loadFactorList(factorsFile);
  loadRegions(regionsFile);
  loadPredictions(predictionsFilestem,chunksDir,targetSpecies,predictions);
  if(cmd.option('g')) {
    loadPredictions(cmd.optParm('g'),chunksDir,targetSpecies,goldStandard);
    useGold=true;
  }
  attachWigFiles(wigDir);
  computeROC();

  if(useGold) contrast();
  
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
      //cout<<"smallestMin="<<smallestMin<<" largestMax="<<largestMax<<endl;
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
				 const String &targetSpecies,
			 Map<String,Vector<GffFeature*> > &predictions) 
{
  Vector<String> subdirs;
  if(!File::getFileList(chunksDir,subdirs)) 
    throw String("Can't get file listing for ")+chunksDir;
  int nDirs=subdirs.size();
  for(int i=0 ; i<nDirs ; ++i) {
    const String &subdir=subdirs[i];
    //if(!numeric.match(subdir)) continue;
    String regionID=subdir; //subdir.asInt();
    String filename=chunksDir+"/"+subdir+"/"+predictionsFilestem+".gff";
    if(!File::exists(filename)) continue;
    int n=loadPredictionsRegion(filename,regionID,targetSpecies,predictions);
    //cout<<n<<" predictions loaded"<<endl;
  }
}


GffFeature *Application::findRegion(const String &regionID)
{
  Vector<GffFeature*>::iterator cur=regions.begin(), end=regions.end();
  for(; cur!=end ; ++cur) {
    GffFeature *region=*cur;
    if(region->getFeatureType()==regionID) return region;
  }
  cout<<"can't find region: "<<regionID<<endl;
  INTERNAL_ERROR;
}



int Application::loadPredictionsRegion(const String &filename,
				       const String &regionID,
				       const String &targetSpecies,
		       Map<String,Vector<GffFeature*> > &predictions)
{
  GffFeature *region=findRegion(regionID);
  int regionOffset=region->getBegin();
  GffReader reader(filename);
  Vector<GffFeature*> *features=reader.loadFeatures();
  if(uniquenessType!=UT_NONE) 
    uniquifyByClustering(*features);
    //uniquify(*features);
  Vector<GffFeature*>::iterator cur=features->begin(), end=features->end();
  int n=0;
  for(; cur!=end ; ++cur) {
    GffFeature *feature=*cur;
    feature->setBegin(regionOffset+feature->getBegin()+1);
    feature->setEnd(regionOffset+feature->getEnd());
    if(feature->getSubstrate()!=targetSpecies) {
      delete feature;
      continue;
    }
    const String &substrate=feature->getSubstrate();
    String ftype=feature->getFeatureType();
    ftype.toupper();
    if(ftype[0]=='-') {
      ftype=ftype.substring(1,ftype.length()-1);
      feature->setFeatureType(ftype);
      feature->setStrand('-');
    }
    else feature->setFeatureType(ftype);
    String key=ftype+" "+regionID;
    predictions[key].push_back(feature);
    ++n;
  }
  delete features;
  return n;
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



void Application::computeROC()
{
  int numRegions=regions.size();
  for(int i=0 ; i<numRegions ; ++i) {
    GffFeature *region=regions[i];
    String chr=region->getSubstrate();
    String regionID=region->getFeatureType();
    if(regionID.asInt()==skipChunk) continue;
    Set<String>::iterator cur=factors.begin(), end=factors.end();
    for(; cur!=end ; ++cur) {
      const String &factor=*cur;
      String key=factor+" "+regionID;
      //if(!predictions.isDefined(key)) INTERNAL_ERROR;
      Vector<GffFeature*> &preds=predictions[key];
      key=chr+" "+factor;
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
    wig->regionsAbove(threshold,sites,boundingRegions);
    maskBlacklistSites(sites,chunk,factor);
    RocPoint &rocPoint=ROC[i];
    rocPoint.threshold=threshold;
    //if(sites.size()>=MIN_SAMPLE_SIZE)
      updateCounts(predictions,sites,rocPoint,chunkBegin,chunkEnd);
    boundingRegions=sites;
  }
}



void Application::updateCounts(Vector<GffFeature*> &predictions,
			       Vector<WigInterval> &sites,
			       RocPoint &rocPoint,int chunkBegin,
			       int chunkEnd)
{
  // Count TP & FP
  int numKnown=sites.size(), numPredicted=predictions.size();
  int &FP=rocPoint.FP, &TP=rocPoint.TP, &FN=rocPoint.FN, &TN=rocPoint.TN;
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

  // compute TN
  Vector<WigInterval> nonSites;
  getNonSites(sites,nonSites,chunkBegin,chunkEnd);
  int numNon=nonSites.size();
  for(int i=0 ; i<numNon ; ++i) {
    WigInterval &nonSite=nonSites[i];
    bool found=false;
    for(int j=0 ; j<numPredicted ; ++j) {
      GffFeature *pred=predictions[j];
      if(pred->getEnd()<nonSite.begin) continue;
      if(pred->getBegin()>nonSite.end) break;
      found=true;
      break;
    }
    if(!found) ++TN;
  }

  // compute FN
  for(int i=0 ; i<numKnown ; ++i)
    if(!seen[i]) ++FN;
}



void Application::getNonSites(Vector<WigInterval> &sites,
			      Vector<WigInterval> &nonSites,
			      int chunkBegin,int chunkEnd)
{
  int n=sites.size();
  for(int i=0 ; i<n ; ++i) {
    WigInterval &site=sites[i];
    if(i==0 && site.begin>chunkBegin+5) 
      nonSites.push_back(WigInterval(chunkBegin,site.begin-1));
    if(i<n-1 && site.end<sites[i+1].begin-5)
      nonSites.push_back(WigInterval(site.end+1,sites[i+1].begin-1));
    if(i==n-1 && site.end<chunkEnd-5)
      nonSites.push_back(WigInterval(site.end+1,chunkEnd));
  }
}



void Application::smoothROC()
{
  for(int j=0 ; j<smoothingIterations ; ++j) {
    for(int i=0 ; i<numPoints ; ++i) {
      RocPoint &pt=ROC[i];
      float sumFPR=0.0, sumSn=0.0;
      int SS=0;
      for(int k=-smoothingNeighborhood ; k<=smoothingNeighborhood ; ++k) {
	if(i+k<0) {
	  ++SS;
	}
	else if(i+k>=numPoints) {
	  sumFPR+=1;
	  sumSn+=1;
	  ++SS;
	}
	else {
	//if(i+k>=0 && i+k<numPoints) {
	  RocPoint &other=ROC[i+k];
	  sumFPR+=other.FPR;
	  sumSn+=other.Sn;
	  ++SS;
	}
      }
      pt.FPR=sumFPR/SS;
      pt.Sn=sumSn/SS;
    }
  }
}



void Application::computeSnSp()
{
  for(int i=0 ; i<numPoints ; ++i) {
    RocPoint &pt=ROC[i];
    pt.Sn=pt.TP+pt.FN>0 ? pt.TP/float(pt.TP+pt.FN) : 0.0;
    pt.Sp=pt.TP+pt.FP>0 ? pt.TP/float(pt.TP+pt.FP) : 0.0;
    pt.FDR=pt.FP+pt.TP>0 ? pt.FP/float(pt.TP+pt.FP) : 0.0;
    pt.FPR=pt.FP+pt.TN>0 ? pt.FP/float(pt.FP+pt.TN) : 0.0;
    pt.N=pt.TP+pt.FN;
    if(pt.FP+pt.TN==0) continue;
    //cout<<pt.Sp<<"\t"<<pt.Sn<<endl;
    //cout<<pt.FDR<<"\t"<<pt.Sn<<endl;
    //cout<<pt.FPR<<"\t"<<pt.Sn<<endl;
    //cout<<pt.threshold<<"\t"<<pt.FPR<<"\t"<<pt.Sn<<endl;
  }
  smoothROC();
  ofstream f;
  if(!dumpFile.isEmpty()) {
    f.open(dumpFile.c_str());
    f<<"thres\tN\tSn\tSp\tFPR\tFDR\n";
  }
  cout<<"0\t0"<<endl;
  for(int i=0 ; i<numPoints ; ++i) {
    RocPoint &pt=ROC[i];
    if(pt.FP+pt.TN==0) continue;
    //cout<<pt.Sp<<"\t"<<pt.Sn<<endl;
    //cout<<pt.FDR<<"\t"<<pt.Sn<<endl;
    cout<<pt.FPR<<"\t"<<pt.Sn<<endl;
    //cout<<pt.threshold<<"\t"<<pt.FPR<<"\t"<<pt.Sn<<endl;
    if(!dumpFile.isEmpty())
      f<<pt.threshold<<"\t"<<pt.N<<"\t"<<pt.Sn<<"\t"<<pt.Sp
       <<"\t"<<pt.FPR<<"\t"<<pt.FDR<<"\n";
  }
  cout<<"1\t1"<<endl;
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



void Application::contrast()
{
  Set<String>::iterator cur=factors.begin(), end=factors.end();
  for(; cur!=end ; ++cur) {
    String factor=*cur;
    int numRegions=regions.size();
    for(int i=0 ; i<numRegions ; ++i) {
      GffFeature *region=regions[i];
      String key=factor+" "+region->getFeatureType();
      Vector<GffFeature*> &pred=predictions[key];
      Vector<GffFeature*> &gold=goldStandard[key];
      int numPred=pred.size(), numGold=gold.size();
      for(int i=0 ; i<numPred ; ++i) {
	GffFeature *p=pred[i];
	bool found=false;
	for(int i=0 ; i<numGold ; ++i) {
	  GffFeature *g=gold[i];
	  if(g->overlapBases(*p)>0) {found=true;break;}
	}
	if(!found) {
	  float score=getScore(*p,*region);
	  Vector<String> &extra=p->getExtraFields();
	  extra.push_back(String("/wigscore=")+score);
	  p->setSubstrate(region->getSubstrate());
	  cerr<<*p<<endl;
	}
      }
    }
  }
}



float Application::getScore(GffFeature &site,GffFeature &region)
{
  String chr=region.getSubstrate();
  String factor=site.getFeatureType();
  String key=chr+" "+factor;
  WigBinary *wig=wigFiles[key];
  if(!wig) return 0;
  int begin=site.getBegin(), end=site.getEnd();
  float score=0.0;
  int L=wig->getLength();
  if(L<end) return 0;
  for(int pos=begin ; pos<=end ; ++pos) score+=wig->read(pos);
  score/=(end-begin+1);
  return score;
}


void Application::loadBlacklist(const String &filename) {
  /*
  Map<String,Vector<GffFeature*> > blacklist;
                               // chr x factor -> vect<feature>
			       */
  GffReader reader(filename);
  Vector<GffFeature*> *features=reader.loadFeatures();
  Vector<GffFeature*>::iterator cur=features->begin(), end=features->end();
  for(; cur!=end ; ++cur) {
    GffFeature *feature=*cur;
    String key=feature->getSubstrate()+" "+feature->getFeatureType();
    blacklist[key].push_back(feature);
  }  
}



void Application::maskBlacklistSites(Vector<WigInterval> &sites,
				     GffFeature *region,
				     const String &factor) {
  /*
  Map<String,Vector<GffFeature*> > blacklist;
                               // chr x factor -> vect<feature>
   */
  const String &chr=region->getSubstrate();
  String key=chr+" "+factor;
  Vector<GffFeature*> &intervals=blacklist[key];
  int numSites=sites.size();
  for(int i=0 ; i<numSites ; ++i) {
    WigInterval &w=sites[i];
    Vector<GffFeature*>::iterator cur=intervals.begin(), end=intervals.end();
    for(; cur!=end ; ++cur) {
      GffFeature *mask=*cur;
      if(!w.overlaps(mask)) continue;
      int mb=mask->getBegin(), me=mask->getEnd();
      if(mb>w.begin && me<w.end) {
	int oldEnd=w.end;
	w.end=mb;
	sites.insertByIndex(WigInterval(me,oldEnd),i+1);
      }
      else if(mb>w.begin) {
	w.begin=mb;
      }
      else {
	w.end=me;
      }
    }
  }
}



void Application::sort(Vector<GffFeature*> &v)
{
  VectorSorter<GffFeature*> sorter(v,cmp);
  switch(uniquenessType)
    {
    case UT_MIN:
      sorter.sortAscendInPlace();
      break;
    case UT_MAX:
      sorter.sortDescendInPlace();
      break;
    }
}



void Application::uniquifyByClustering(Vector<GffFeature*> &features)
{
  sort(features);
  int n=features.size();
  for(int i=0 ; i<n ; ++i) {
    GffFeature *best=features[i];
    if(!best) continue;
    for(int j=i+1 ; j<n ; ++j) {
      GffFeature *&other=features[j];
      if(!other) continue;
      //if(best->overlapBases(*other)>0) {
      if(best->overlapBases(*other)>=best->length()/2) {
	delete other;
	other=NULL;
      }
    }
  }
  Vector<GffFeature*> temp;
  for(int i=0 ; i<n ; ++i) {
    GffFeature *feature=features[i];
    if(feature) temp.push_back(feature);
  }
  features=temp;
}


