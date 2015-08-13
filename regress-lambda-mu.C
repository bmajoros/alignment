/****************************************************************
 regress-lambda-mu.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GSL/Optimizer.H"
#include "BOOM/Array1D.H"
#include "ModelCompiler.H"
using namespace std;
using namespace BOOM;

const double EPSILON=0.1; //2.2e-5; //### BEST VALUE IS 0.1 ###
const double stepSize=0.001;//0.01;
const double tolerance=0.1;//0.1;
const double gradientThreshold=0.001;//0.01; //1.0;
const int maxIterations=1000;


/****************************************************************
                         struct ErrorFunc
****************************************************************/
struct ErrorFunc : public GSL::ObjectiveFunction
{
  Array1D<ModelParm> &parms;
  double S, L, Ng, Nl, t;
  ErrorFunc(Array1D<ModelParm> &,double S, double L, double Ng, 
	    double Nl,double t);
  virtual double f(const GSL::Vector &currentPoint);
  virtual void gradient(const GSL::Vector &currentPoint,
			GSL::Vector &gradient);
};



/****************************************************************
                         class Application
****************************************************************/
class Application
{
  Array1D<ModelParm> parms;
  double S, L, Ng, Nl, t;
  void regress();
public:
  Application();
  int main(int argc,char *argv[]);
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
                             methods
****************************************************************/

Application::Application()
  : parms(2)
{
  parms[0]=ModelParm("lambda",0.00001,100.0);
  parms[1]=ModelParm("mu",0.00001,100.0);
}



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=5)
      throw String("regress-lambda-mu <sites-in-anc> <seq-len> <#gains> <#losses> <branchlen>");
    S=(double) cmd.arg(0).asInt();
    L=(double) cmd.arg(1).asInt();
    Ng=(double) cmd.arg(2).asInt();
    Nl=(double) cmd.arg(3).asInt();
    t=cmd.arg(4).asFloat();

    regress();

    return 0;
  }



void Application::regress()
{
  ErrorFunc f(parms,S,L,Ng,Nl,t);
  GSL::Vector initialPoint(2);
  initialPoint[0]=parms[0].unconstrain(0.01); // lambda
  initialPoint[1]=parms[1].unconstrain(0.08);  // mu
  GSL::Optimizer optimizer(GSL::BFGS,f,initialPoint,stepSize,GSL::BY_EITHER,
			   tolerance,gradientThreshold,maxIterations);
  optimizer.run();
  const GSL::Vector &optimalPoint=optimizer.getOptimalPoint();
  cout<<"lambda="<<parms[0].constrain(optimalPoint[0])
      <<" mu="<<parms[1].constrain(optimalPoint[1])<<endl;
}



ErrorFunc::ErrorFunc(Array1D<ModelParm> &parms,double S, double L, 
		     double Ng, double Nl, double t)
  : S(S), L(L), Ng(Ng), Nl(Nl), t(t), parms(parms)
{
  // ctor
}



double ErrorFunc::f(const GSL::Vector &currentPoint)
{
  double lambda=parms[0].constrain(currentPoint[0]);
  double mu=parms[1].constrain(currentPoint[1]);
  double lm=lambda+mu;
  double exp_Ng=lambda/lm*L*(1-exp(-lm*t));
  double exp_Nl=mu/lm*S*(1-exp(-lm*t));
  double diff_Ng=Ng-exp_Ng, diff_Nl=Nl-exp_Nl;
  double y=diff_Ng*diff_Ng + diff_Nl*diff_Nl;
  //cout<<"currentPoint="<<currentPoint<<endl;
  //cout<<"error="<<y<<" Ng:"<<exp_Ng<<":"<<Ng<<" Nl:"<<exp_Nl<<":"<<Nl<<endl;
  return y;
}



void ErrorFunc::gradient(const GSL::Vector &currentPt,
			 GSL::Vector &gradient)
{
  const double epsilon=EPSILON;
  GSL::Vector perturbed=currentPt;
  int n=gradient.getDim();
  for(int i=0 ; i<n ; ++i)
    {
      double &x=perturbed[i];
      const double x0=x;
      double dx=epsilon*x;
      double temp=x+dx; 
      dx=temp-x; // make sure it's an "exact value" in the machine
      double twoDX=2*dx;
      x-=dx;
      double y1=f(perturbed);
      x+=twoDX;
      double y2=f(perturbed);
      gradient[i]=(y2-y1)/twoDX;
      x=x0;
    }
  //cout<<"gradient = "<<gradient<<endl;
}



