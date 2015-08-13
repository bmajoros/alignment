#include <iostream>
#include "BOOM/Exceptions.H"
#include "BOOM/RouletteWheel.H"
#include "BOOM/Random.H"
using namespace std;
using namespace BOOM;

void go();
void dump(float *D);

int main(int argc,char *argv[])
  {
    try
      {
	go();
	return 0;
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
    catch(const RootException &e)
      {
	cerr << "Exception: "<<e.getMessage()<<endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }



void go()
{
  randomize();
  float D1[10], D2[10], sum=0;
  for(int i=0 ; i<10 ; ++i) {
    D1[i]=1.0/(i+6);
    sum+=D1[i];
  }
  for(int i=0 ; i<10 ; ++i) D1[i]/=sum;
  sum=0;
  for(int i=0 ; i<10 ; ++i) {
    D2[i]=1.0/(10-i);
    sum+=D2[i];
  }
  for(int i=0 ; i<10 ; ++i) D2[i]/=sum;
  RouletteWheel W1, W2;
  for(int i=0 ; i<10 ; ++i) {
    W1.addSector(D1[i]);
    W2.addSector(D2[i]);
  }	
  W1.doneAddingSectors();
  W2.doneAddingSectors();
  float D3[10];
  sum=0;
  for(int i=0 ; i<10 ; ++i) {
    D3[i]=D1[i]*D2[i];
    sum+=D3[i];
  }
  for(int i=0 ; i<10 ; ++i) D3[i]/=sum;
  float D4[10];
  sum=0;
  for(int i=0 ; i<10 ; ++i) {
    D4[i]=D1[i]+D2[i];
    sum+=D4[i];
  }
  for(int i=0 ; i<10 ; ++i) D4[i]/=sum;
  
  dump(D2);  cout<<endl;
  dump(D1);  cout<<endl;
  dump(D3);  cout<<endl;
  //dump(D4);  cout<<endl;
  
  
  float observed[10];
  for(int i=0 ; i<10 ; ++i) observed[i]=0;
  const int S2=100000;
  for(int i=0 ; i<S2 ; ++i) {
    const int S=100;
    int pool[S];
    for(int i=0 ; i<S ; ++i) pool[i]=W1.spin();
    float hist[S];
    sum=0;
    for(int i=0 ; i<S ; ++i) {
      hist[i]=D2[pool[i]];
      sum+=hist[i];
    }
    RouletteWheel W3;
    for(int i=0 ; i<S ; ++i) {
      hist[i]/=sum;
      W3.addSector(hist[i]);
    }
    W3.doneAddingSectors();
    
    int x=pool[W3.spin()];
    ++observed[x];
    ++sum;
  }
  for(int i=0 ; i<10 ; ++i) {
    observed[i]/=float(S2);
    cout<<i<<"\t"<<observed[i]<<endl;
  }
}


void dump(float *D)
{
  for(int i=0 ; i<10 ; ++i)
    cout<<i<<"\t"<<D[i]<<endl;
}
