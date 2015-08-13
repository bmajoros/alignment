/****************************************************************
 restore-nj-root.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Map.H"
#include "BOOM/Constants.H"
#include "PhyLib/Phylogeny.H"
using namespace std;
using namespace BOOM;


class Application
{
  Map<String,double> branchLengths;
  void transfer(Phylogeny &rerootedTree,Phylogeny &njTree);
  void changeBranchLen(PhylogenyNode *parent,PhylogenyNode *child,
		       Phylogeny &,double len);
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
    throw String("restore-nj-root <nj.phy> <rerooted.phy <out.phy>");
  String njFilename=cmd.arg(0);
  String rerootedFilename=cmd.arg(1);
  String outFilename=cmd.arg(2);
  
  // Load phylogenies
  Phylogeny *njTree=new Phylogeny(njFilename);
  Phylogeny *rerootedTree=new Phylogeny(rerootedFilename);
  
  // Transfer branch lengths from the rerooted tree to the original nj tree
  transfer(*rerootedTree,*njTree);
  
  // Output the rescaled nj tree
  ofstream os(outFilename.c_str());
  os<<*njTree<<endl;
  
  return 0;
}



void Application::transfer(Phylogeny &rerootedTree,Phylogeny &njTree)
{
  /*
  struct Visitor : public Treevisitor {
    void processNode(InternalNode &) {}
    void processNode(LeafNode &) {}
    void processNode(RootNode &) {}
  } V();
  */

  // Get rerootedTree's branch lengrhs into a map for quick lookup
  Vector<PhylogenyBranch> &rerootedBranches=*rerootedTree.gatherBranches();
  int numRerootedBranches=rerootedBranches.size();
  for(int i=0 ; i<numRerootedBranches ; ++i) {
    PhylogenyBranch branch=rerootedBranches[i];
    String parentName=branch.getParent()->getName();
    String childName=branch.getChild()->getName();
    String key=parentName+" "+childName;
    branchLengths[key]=branch.getLength();
  }

  // Process the original njTree's branches and change their lengths
  Vector<PhylogenyBranch> &njBranches=*njTree.gatherBranches();
  int numNjBranches=njBranches.size();
  for(int i=0 ; i<numNjBranches ; ++i) {
    PhylogenyBranch branch=njBranches[i];
    String parentName=branch.getParent()->getName();
    String childName=branch.getChild()->getName();
    String key=parentName+" "+childName, revKey=childName+" "+parentName;
    double newLength=NEGATIVE_INFINITY;
    if(branchLengths.isDefined(key)) newLength=branchLengths[key];
    else if(branchLengths.isDefined(revKey)) newLength=branchLengths[revKey];
    else {
      PhylogenyNode *root=njTree.getRoot();
      InternalNode *iRoot=dynamic_cast<InternalNode*>(root);
      if(!iRoot) throw "Root is not an internal node";
      PhylogenyNode *left=iRoot->getLeft(), *right=iRoot->getRight();
      String leftName=left->getName(), rightName=right->getName();
      String key=leftName+" "+rightName, revKey=rightName+" "+leftName;
      double dLen=NEGATIVE_INFINITY;
      if(branchLengths.isDefined(key)) dLen=branchLengths[key];
      else if(branchLengths.isDefined(revKey)) dLen=branchLengths[revKey];
      else throw String("can't find ")+leftName+":"+rightName;
      newLength=dLen/2.0;
    }

    // Change the branch's length
    changeBranchLen(branch.getParent(),branch.getChild(),njTree,newLength);
  }
}



void Application::changeBranchLen(PhylogenyNode *parent,PhylogenyNode *child,
				  Phylogeny &tree,double len)
{
  switch(parent->getNodeType())
    {
    case ROOT_NODE:
      static_cast<RootNode*>(parent)->setBranchLength(len);
      break;
    case LEAF_NODE:
      throw "changeBranchLen: LEAF not possible!";
      break;
    case INTERNAL_NODE:
      {
	InternalNode *iNode=static_cast<InternalNode*>(parent);
	switch(iNode->whichChildIsThis(child)) 
	  {
	  case LEFT:   iNode->setLeftDistance(len); break;
	  case RIGHT:  iNode->setRightDistance(len); break;
	  }
      }
      break;
    }
}


