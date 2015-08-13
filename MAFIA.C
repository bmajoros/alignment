/****************************************************************
 AVES-progressive.C
 william.majoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3 
 (see www.opensource.org).
 ****************************************************************/
#include <fstream>
#include "AVES.H"
using namespace std;
using namespace BOOM;
using BOOM::Symbol;


/****************************************************************
                           main()
 ****************************************************************/
int main(int argc,char *argv[])
  {
    try
      {
	AVES app;
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
    catch(const BOOM::RootException &e)
      {
	cerr<<e.getMessage()<<endl;
      }
    catch(...)
      {
	cerr << "Unknown exception caught in main" << endl;
      }
    return -1;
  }




