/****************************************************************
 SparseMatrix3D.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "SparseMatrix3D.H"
#include "BOOM/File.H"
#include "BOOM/Exceptions.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


SparseMatrix3D::SparseMatrix3D(short L1,short L2,short numStates)
  : M(L1,numStates), L1(L1), L2(L2), numStates(numStates), 
    iter(L1,numStates)
{
  // ctor
}



SparseMatrix3D::SparseMatrix3D(Array3D<float> &M,float threshold)
{
  initFrom(M,log(threshold));
}



SparseMatrix3D::SparseMatrix3D(SparseMatrix3D &other,float threshold)
  : M(other.L1,other.numStates), L1(other.L1),
    L2(other.L2), numStates(other.numStates), 
    iter(other.L1,other.numStates)
{
  threshold=log(threshold);
  for(int i=0 ; i<L1 ; ++i)
    for(int j=0 ; j<numStates ; ++j) {
      EntryList &row=M[i][j], &otherRow=other.M[i][j];
      EntryList::iterator cur=otherRow.begin(), end=otherRow.end();
      for(; cur!=end ; ++cur) {
	Entry &e=*cur;
	if(e.value>=threshold) 
	  row.append(e);
      }
    }
}



SparseMatrix3D::EntryList &SparseMatrix3D::operator()(short x,short q)
{
  return M[x][q];
}



void SparseMatrix3D::addEntry(short x,short y,short q,float value)
{
  M[x][q].push_back(Entry(y,value));
}



SparseMatrix3D *SparseMatrix3D::transpose()
{
  SparseMatrix3D &T=*new SparseMatrix3D(L2,L1,numStates);
  //Array2D<EntryList> &TM=T.M;
  for(short x=0 ; x<L1 ; ++x) {
    Array2D<EntryList>::RowIn2DArray<EntryList> row=M[x];
    for(short q=0 ; q<numStates ; ++q) {
      EntryList &entries=row[q];
      EntryList::iterator cur=entries.begin(), end=entries.end();
      for(; cur!=end ; ++cur) {
	Entry &entry=*cur;
	short y=entry.y;
	float value=entry.value;
	short newQ=q;
	if(numStates==3) { if(newQ>0) newQ=3-newQ; }
	else INTERNAL_ERROR; // must swap ins/del!
	T.addEntry(y,x,newQ,value);
      }	
    }
  }
  return &T;
}



SparseMatrix3D *SparseMatrix3D::loadBinary(const String &filename)
{
  //bool debug=true;//filename.contains("dm3") && filename.contains("dp4");
  File file(filename);
  int L1=file.readShort();
  int L2=file.readShort();
  int numStates=file.readShort();
  int numEntries=file.readInt();
  //cout<<L1<<" "<<L2<<" "<<numStates<<" "<<numEntries<<endl;
  SparseMatrix3D *PM=new SparseMatrix3D(L1,L2,numStates);
  PM->iter.resize(L1,numStates);
  Array2D<EntryList> &M=PM->M;
  M.resize(L1,numStates);
  for(int e=0 ; e<numEntries ; ++e) {
    int i=file.readShort();
    int j=file.readShort();
    int k=file.readShort();
    float value=file.readFloat();
    //if(debug /*&& value>-1*/) cout<<"xxx "<<i<<" "<<j<<" "<<k<<" "<<value<<" "<<filename<<endl;
    PM->addEntry(i,j,k,value);
  }
  return PM;
}



void SparseMatrix3D::saveBinary(const String &filename)
{
  File file(filename,"w");
  file.write(L1);
  file.write(L2);
  file.write(numStates);
  int numEntries=0;
  for(short i=0 ; i<L1 ; ++i) {
    Array2D<EntryList>::RowIn2DArray<EntryList> row=M[i];
    for(short q=0 ; q<numStates ; ++q)
      numEntries+=row[q].size();
  }
  file.write(numEntries);
  for(short i=0 ; i<L1 ; ++i) {
    Array2D<EntryList>::RowIn2DArray<EntryList> row=M[i];
    for(short q=0 ; q<numStates ; ++q) {
      EntryList &entries=row[q];
      EntryList::iterator cur=entries.begin(), end=entries.end();
      for(; cur!=end ; ++cur) {
	Entry &entry=*cur;
	file.write(i);
	file.write(entry.y);
	file.write(q);
	file.write(entry.value);
      }
    }
  }
}



void SparseMatrix3D::initFrom(Array3D<float> &A,float threshold)
{
  L1=A.getFirstDim();
  L2=A.getSecondDim();
  numStates=A.getThirdDim();
  iter.resize(L1,numStates);
  M.resize(0,0);
  M.resize(L1,numStates);
  for(int i=0 ; i<L1 ; ++i)
    for(int j=0 ; j<L2 ; ++j)
      for(int q=0 ; q<numStates ; ++q) {
	float value=A[i][j][q];
	if(value>=threshold) addEntry(i,j,q,value);
      }
}



float SparseMatrix3D::find(short x,short y,short state)
{
  EntryList &row=M(x,state);
  EntryList::iterator cur=row.begin(), end=row.end();
  for(; cur!=end ; ++cur) {
    Entry e=*cur;
    if(e.y<y) continue;
    if(e.y==y) return e.value;
    break;
  }
  return NEGATIVE_INFINITY;
}



List<SparseMatrix3D::Entry>::iterator &SparseMatrix3D::getIter(int x,STATE q)
{
  return iter[x][q];
}



List<SparseMatrix3D::Entry>::iterator SparseMatrix3D::begin(int x,STATE q)
{
  return M[x][q].begin();
}



List<SparseMatrix3D::Entry>::iterator SparseMatrix3D::end(int x,STATE q)
{
  return M[x][q].end();
}



