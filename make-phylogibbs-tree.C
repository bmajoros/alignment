/****************************************************************
 make-phylogibbs-tree.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "PhyLib/Phylogeny.H"
#include "PhyLib/RateMatrix.H"
using namespace std;
using namespace BOOM;


class Application {
  RateMatrix *Q;
  double recompute(double);
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
    throw String("make-phylogibbs-tree <in.tree> <in.matrix> <out.tree>");
  String treeFile=cmd.arg(0);
  String matrixFile=cmd.arg(1);
  String outfile=cmd.arg(2);
  
  Phylogeny tree(treeFile);
  Q=RateMatrix::load(matrixFile);

  Vector<PhylogenyNode*> &nodes=*tree.gatherNodes();
  int n=nodes.size();
  for(int i=0 ; i<n ; ++i) {
    PhylogenyNode *node=nodes[i];
    switch(node->getNodeType()) 
      {
      case ROOT_NODE:
	{
	  RootNode *root=static_cast<RootNode*>(node);
	  double t=root->getBranchLength();
	  t=recompute(t);
	  root->setBranchLength(t);
	  break;
	}
      case LEAF_NODE:
	break;
      case INTERNAL_NODE:
	{
	  InternalNode *inode=static_cast<InternalNode*>(node);
	  double t=inode->getLeftDistance();
	  t=recompute(t);
	  inode->setLeftDistance(t);
	  t=inode->getRightDistance();
	  t=recompute(t);
	  inode->setRightDistance(t);
	}
      }
  }
  delete &nodes;

  tree.save(outfile);
  return 0;
}



double Application::recompute(double t)
{
  SubstitutionMatrix &M=*Q->instantiate(t);
  Array1D<double> eq;
  M.getEqFreqs(eq);
  int n=eq.size();
  double r=0.0;
  for(BOOM::Symbol s=0 ; s<n ; ++s) {
    r+=eq[s]*M(s,s);
  }
  delete &M;
  return r;
}




