/****************************************************************
 SparseMatrix2D.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <iostream>
#include "SparseMatrix2D.H"
#include "BOOM/File.H"
#include "BOOM/Constants.H"
using namespace std;
using namespace BOOM;


SparseMatrix2D::SparseMatrix2D(short L1,short L2)
  : M(L1), L1(L1), L2(L2)
{
  // ctor
}



SparseMatrix2D::SparseMatrix2D(Array2D<float> &M,float threshold)
{
  initFrom(M,log(threshold));
}



SparseMatrix2D::SparseMatrix2D(short L1)
  : M(L1), L1(L1), L2(L2)
{
  becomeIdentity();
}



SparseMatrix2D::EntryList &SparseMatrix2D::operator()(short x)
{
  return M[x];
}



void SparseMatrix2D::addEntry(short x,short y,float value)
{
  M[x].push_back(Entry(y,value));
}



void SparseMatrix2D::becomeIdentity()
{
  for(int i=0 ; i<L1 ; ++i)
    addEntry(i,i,LOG_1);
}



SparseMatrix2D *SparseMatrix2D::transpose()
{
  SparseMatrix2D &T=*new SparseMatrix2D(L2,L1);
  for(short x=0 ; x<L1 ; ++x) {
    EntryList &row=M[x];
    EntryList::iterator cur=row.begin(), end=row.end();
    for(; cur!=end ; ++cur) {
      Entry &entry=*cur;
      short y=entry.y;
      float value=entry.value;
      T.addEntry(y,x,value);
    }	
  }
  return &T;
}



SparseMatrix2D *SparseMatrix2D::loadBinary(const String &filename)
{
  File file(filename);
  int L1=file.readShort();
  int L2=file.readShort();
  short numEntries=file.readShort();
  SparseMatrix2D *PM=new SparseMatrix2D(L1,L2);
  Array1D<EntryList> &M=PM->M;
  M.resize(L1);
  for(int e=0 ; e<numEntries ; ++e) {
    short i=file.readShort();
    short j=file.readShort();
    float value=file.readFloat();
    PM->addEntry(i,j,value);
  }
  return PM;
}



void SparseMatrix2D::saveBinary(const String &filename)
{
  File file(filename,"w");
  file.write(L1);
  file.write(L2);
  short numEntries=0;
  for(int i=0 ; i<L1 ; ++i) {
    EntryList &row=M[i];
    numEntries+=row.size();
  }
  file.write(numEntries);
  //cout<<"numEntries="<<numEntries<<endl;
  for(short i=0 ; i<L1 ; ++i) {
    EntryList &entries=M[i];
    EntryList::iterator cur=entries.begin(), end=entries.end();
    for(; cur!=end ; ++cur) {
      Entry &entry=*cur;
      file.write(i);
      file.write(entry.y);
      file.write(entry.value);
      //cout<<"\t"<<i<<" "<<entry.y<<" "<<entry.value<<endl;
    }
  }
}



void SparseMatrix2D::initFrom(Array2D<float> &A,float threshold)
{
  L1=A.getFirstDim();
  L2=A.getSecondDim();
  M.resize(0);
  M.resize(L1);
  for(int i=0 ; i<L1 ; ++i) {
    for(int j=0 ; j<L2 ; ++j) {
      float value=A[i][j];
      if(value>=threshold) addEntry(i,j,value);
    }
  }
}



float SparseMatrix2D::find(short x,short y)
{
  EntryList &entries=M[x];
  EntryList::iterator cur=entries.begin(), end=entries.end();
  for(; cur!=end ; ++cur) {
    Entry &entry=*cur;
    if(entry.y==y) return entry.value;
  }
  return LOG_0;
}




