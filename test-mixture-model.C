/****************************************************************
 test-mixture-model.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Random.H"
#include "BOOM/Exceptions.H"
#include "PhyLib/JukesCantor.H"
#include "PhyLib/REV.H"
#include "PhyLib/SubstitutionMixtureModel.H"
using namespace std;
using namespace BOOM;




class Application
{
  void checkSumTo1(SubstitutionMatrix &M);
  SubstitutionMatrix *integrate(const RateMatrix &Q1,const RateMatrix &Q2,
				double t,int steps);
  SubstitutionMatrix *compose(const RateMatrix &Q1,const RateMatrix &Q2,
				double t,double s);
  void install(const GSL::Matrix &from,GSL::Matrix &to,int leftCol,
	       int topRow,float p);
  void pullOut(const GSL::Matrix &from,GSL::Matrix &to,int leftCol,
	       int topRow);
  SubstitutionMatrix *megaMatrixMethod(const RateMatrix &Q1,
				       const RateMatrix &Q2,double t);
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
    // ctor
  }



int Application::main(int argc,char *argv[])
  {
    // Process command line
    CommandLine cmd(argc,argv,"");
    if(cmd.numArgs()!=3)
      throw String("test-mixture-model <t> <steps> <seed>");
    double t=cmd.arg(0).asFloat();
    int steps=cmd.arg(1).asInt();
    int seed=cmd.arg(2).asInt();
    SeedRandomizer(seed);//randomize();

    RateMatrix &Q1=*RateMatrix::random(DNA,MT_REV);
    RateMatrix &Q2=*RateMatrix::random(DNA,MT_REV);
    SubstitutionMixtureModel M(Q1,Q2,t);
    /*
    SubstitutionMixtureModel Mz(Q2,Q1,t);
    SubstitutionMatrix *M1=Q1.instantiate(t), *M2=Q2.instantiate(t);

    cout<<"Pt(Q1)=\n"<<*M1<<endl;
    cout<<"Pt(Q2)=\n"<<*M2<<endl;
    cout<<"reverse mixture=\n"<<Mz<<endl;
    cout<<"mixture=\n"<<M<<endl;

    SubstitutionMatrix *C=compose(Q1,Q2,t,t/2);
    cout<<"half-and-half:\n"<<*C<<endl;

    M1->peek().add(M2->peek(),C->peek());
    C->peek().scale(0.5);
    cout<<"simple average:\n"<<*C<<endl;
    */
    SubstitutionMatrix *I=integrate(Q1,Q2,t,steps);
    //cout<<"numerically integrated:\n"<<*I<<endl;
    GSL::Matrix residuals;
    M.peek().subtract(I->peek(),residuals);
    cout<<"residuals:\n"<<residuals<<endl;
    //checkSumTo1(M); // looked OK when I ran it

    SubstitutionMatrix *mm=megaMatrixMethod(Q1,Q2,t);
    M.peek().subtract(mm->peek(),residuals);
    cout<<"MM residuals:\n"<<residuals<<endl;
    SubstitutionMatrix &P1=*Q1.instantiate(t), &P2=*Q2.instantiate(t);
    cout<<*mm<<"\n"<<P1<<"\n"<<P2<<"\n"<<M<<endl;

    return 0;
  }



void Application::checkSumTo1(SubstitutionMatrix &M)
{
  int n=M.getAlphabet().size();
  for(Symbol s=0 ; s<n ; ++s) {
    double sum=0.0;
    for(Symbol t=0 ; t<n ; ++t)
      sum+=M(s,t);
    cout<<s<<": "<<sum<<endl;
  }
}



SubstitutionMatrix *Application::integrate(const RateMatrix &Q1,
					   const RateMatrix &Q2,
					   double t,int steps)
{
  SubstitutionMatrix *r=NULL;
  if(steps<2) throw "steps must be >1";
  double inc=t/(steps-1);
  for(int i=0 ; i<steps ; ++i) {
    double s=inc*i;
    SubstitutionMatrix *Pt=compose(Q1,Q2,t,s);
    if(!r) { r=Pt; continue; }
    GSL::Matrix temp;
    r->peek().add(Pt->peek(),temp);
    r->peek()=temp;
    delete Pt;
  }
  r->peek().scale(1.0/steps);
  return r;
}



SubstitutionMatrix *Application::compose(const RateMatrix &Q1,
					 const RateMatrix &Q2,
					 double t,double s)
{
  SubstitutionMatrix *M1=Q1.instantiate(s), *M2=Q2.instantiate(t-s);
  GSL::Matrix M3;
  M1->peek().times(M2->peek(),M3);
  M1->peek()=M3;
  delete M2;
  return M1;
}



void Application::install(const GSL::Matrix &from,GSL::Matrix &to,
			  int leftCol,int topRow,float p)
{
  int width=from.getNumColumns(), height=from.getNumRows();
  for(int x=0 ; x<width ; ++x)
    for(int y=0 ; y<height ; ++y)
      to(leftCol+x,topRow+y)=p*from(x,y);
}



void Application::pullOut(const GSL::Matrix &from,GSL::Matrix &to,
			  int leftCol,int topRow)
{
  int width=to.getNumColumns(), height=to.getNumRows();
  for(int x=0 ; x<width ; ++x)
    for(int y=0 ; y<height ; ++y)
      to(x,y)=from(leftCol+x,topRow+y);
}



SubstitutionMatrix *Application::megaMatrixMethod(const RateMatrix &Q1,
						  const RateMatrix &Q2,
						  double t)
{
  const float EPSILON=1.0;

  // First, pack matrices and their products into a "mega" matrix
  const GSL::Matrix &M1=Q1.peek(), &M2=Q2.peek();
  const int m=M1.getNumColumns(), doubleM=2*m;
  SubstitutionMatrix &P1=*Q1.instantiate(EPSILON);
  SubstitutionMatrix &P2=*Q2.instantiate(EPSILON);
  GSL::Matrix PstarEps(doubleM,doubleM);
  GSL::Matrix M12(m,m), M21(m,m);
  P1.peek().times(P2.peek(),M12);
  P2.peek().times(P1.peek(),M21);
  install(P1.peek(),PstarEps,0,0,0.9);
  install(P2.peek(),PstarEps,m,0,0.1);
  install(P1.peek(),PstarEps,0,m,0.1);
  install(P2.peek(),PstarEps,m,m,0.9);
  for(int i=0 ; i<doubleM ; ++i) {
    float sum=0.0;
    for(int j=0 ; j<doubleM ; ++j) sum+=PstarEps(i,j);
    for(int j=0 ; j<doubleM ; ++j) PstarEps(i,j)/=sum;
  }
  cout<<"P*(eps)=\n"<<PstarEps<<endl;

  // Decompose into eigenvectors and eigenvalues
  GSL::Matrix G; // eigenvectors
  GSL::Vector eigenvalues;
  GSL::Vector imagParts; // should be all zeros for practical rate matrices
  //P1.peek().getEigenVectors(G,eigenvalues,imagParts);
  //cout<<"eigenvectors=\n"<<G<<"\neigenvalues="<<eigenvalues<<"\nimagParts="<<imagParts<<endl;
  PstarEps.getEigenVectors(G,eigenvalues,imagParts);
  int n=PstarEps.getNumRows();
  cout<<"eigenvectors=\n"<<G<<"\neigenvalues="<<eigenvalues<<"\nimagParts="<<imagParts<<endl;

  // Prepare Ginverse
  GSL::Matrix Ginverse;
  G.invert(Ginverse);

  // Prepare matrix U (which consists of exp(t*lambda) diagonal values)
  GSL::Matrix U(n,n);
  U.setAllTo(0.0);
  for(int i=0 ; i<n ; ++i)
    U(i,i)=exp(t*log(eigenvalues[i])/EPSILON);

  // Form the matrix product P(t) = G * U * Ginverse
  GSL::Matrix PstarT(n,n);
  GSL::Matrix GU;
  G.times(U,GU);
  GU.times(Ginverse,PstarT);
  cout<<"P*(t)=\n"<<PstarT<<endl;

  // Pull out the desired submatrix
  const Alphabet &alphabet=Q1.getAlphabet();
  AlphabetMap &alphabetMap=Q1.getAlphabetMap();
  SubstitutionMatrix *substMatrix=
    new SubstitutionMatrix(alphabet,alphabetMap);
  pullOut(PstarT,substMatrix->peek(),m,0);
  for(int i=0 ; i<m ; ++i) {
    float sum=0.0;
    for(int j=0 ; j<m ; ++j) sum+=substMatrix->peek()(i,j);
    for(int j=0 ; j<m ; ++j) substMatrix->peek()(i,j)/=sum;
  }
  return substMatrix;
}



