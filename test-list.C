
#include <iostream>
#include <list>
using namespace std;

int main(int argc,char *argv[])
{
  list<int> x;
  x.push_back(2);
  x.push_back(4);
  x.push_back(6);
  list<int>::iterator i=x.end();
  while(1) {
    --i;
    cout<<*i<<endl;
    if(i==x.begin()) break;
  }
  return 0;
}

