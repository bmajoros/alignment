/****************************************************************
 test-mpi.C
 bmajoros@duke.edu

 This is open-source software, governed by the Gnu General Public License (GPL) version 3.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Vector.H"
#include "BOOM/Random.H"
#include "BOOM/Constants.H"
#include "BOOM/Environment.H"
#include "BOOM/Time.H"
#include "BOOM/Exceptions.H"
#include "BOOM/MPI.H"
using namespace std;
using namespace BOOM;


/****************************************************************
                         enum MessageType
 ****************************************************************/
enum MessageTag {
  MSG_WORK,
  MSG_ANSWER
};


/****************************************************************
                       class Application
 ****************************************************************/
class Application
{
  MPI *mpi;
  MpiSlaveDriver *slaveDriver;
  int master(int argc,char *argv[]);
  int slave(int argc,char *argv[]);
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
  catch(const String &msg)
    {
      cerr << msg.c_str() << endl;
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
      cerr << "exception: "<< e.getMessage() << endl;
    }
  catch(...)
    {
      cerr << "Unknown exception caught in main" << endl;
    }
  return -1;
}




/****************************************************************
                      Application methods
 ****************************************************************/

Application::Application()
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Initialize MPI
  mpi=new MPI(argc,&argv);
  slaveDriver=new MpiSlaveDriver(*mpi);
  int exitCode;

  CommandLine cmd(argc,argv,"b:vp:t:e:G");
  cout<<"number of args: "<<cmd.numArgs()<<endl;

  // Determine whether this process is a MASTER or a SLAVE
  if(mpi->getProcessID()==0) exitCode=master(argc,argv);
  else exitCode=slave(argc,argv);

  //MPI_Abort(MPI_COMM_WORLD,exitCode);
  return exitCode;
}



/****************************************************************
                             MASTER
 ****************************************************************/
int Application::master(int argc,char *argv[])
{
  // Process command line
  //CommandLine cmd(argc,argv,"b:vp:t:e:G");
  //cout<<"master args: "<<cmd.numArgs()<<endl;
  //if(cmd.numArgs()!=0)
  //  throw String("test-mpi <no args>");

  // Initialize the slaves
  cout<<"MASTER init slaves"<<endl;
  for(int i=0 ; i<100 ; ++i) {
    MpiVariableMsg &message=*new MpiVariableMsg(MSG_WORK);
    message<<i;
    message<<endmsg;
    slaveDriver->addWork(&message);
  }
  Vector<MpiFixedMsg*> results;
  slaveDriver->waitForResults(results);
  int n=results.size();
  for(int i=0 ; i<n ; ++i) {
    MpiFixedMsg *msg=results[i];
    int x, x2;
    (*msg)>>x>>x2;
    cout<<"RESULT: "<<x<<" "<<x2<<endl;
    delete msg;
  }
  cout<<"terminating slaves..."<<endl;
  slaveDriver->terminateSlaves();
  cout<<"MASTER exiting..."<<endl;
  return 0;
}



/****************************************************************
                              SLAVE
 ****************************************************************/
int Application::slave(int argc,char *argv[])
{
  // Process command line
  //CommandLine cmd(argc,argv,"b:vp:t:e:G");
  //cout<<"slave args: "<<cmd.numArgs()<<endl;
  //if(cmd.numArgs()!=0)
  //  throw String("test-mpi <no args>");

  // Serve the master...
  bool done=false;
  while(!done) {
    cout<<mpi->getProcessID()<<" waiting for work"<<endl;
    MpiFixedMsg &work=*slaveDriver->acceptWork();
    cout<<mpi->getProcessID()<<" received work"<<endl;
    switch(work.getTag()) 
      {
      case MSG_WORK:
	{
	  int x;
	  work>>x;
	  int x2=x*x;
	  cout<<x<<"^2="<<x2<<" (slave #"<<mpi->getProcessID()<<")"<<endl;
	  MpiVariableMsg &reply=*new MpiVariableMsg(MSG_ANSWER);
	  reply<<x<<x2;
	  reply<<endmsg;
	  cout<<"replying to master "<<mpi->getProcessID()<<endl;
	  slaveDriver->replyToMaster(&reply);
	  cout<<"done replying "<<mpi->getProcessID()<<endl;
	}
	break;
      case TERMINATE_SLAVE:
	cout<<"slave "<<mpi->getProcessID()<<" quitting"<<endl;
	done=true;
	break;
      default:
	cout<<"SLAVE: unknown message!"<<endl;
	break;
      }
    delete &work;
  }
  return 0;
}



/*

int Application::master(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=0)
    throw String("test-mpi <no args>");

  // Initialize the slaves
  cout<<"MASTER init slaves"<<endl;
  //MpiFixedMsg message(MSG_EXAMPLE_1,sizeof(int)+sizeof(double));
  MpiVariableMsg message(MSG_EXAMPLE_1);
  int x=800;
  double y=5.6;
  message<<x<<y;
  message<<endmsg;
  int numSlaves=mpi->getNumSlaves();
  for(int i=0 ; i<numSlaves ; ++i) {
    int slaveID=mpi->getIthSlaveID(i);
    mpi->send(message,slaveID);
  }

  cout<<"MASTER terminating slaves"<<endl;
  // Terminate the slaves
  for(int i=0 ; i<numSlaves ; ++i) {
    int slaveID=mpi->getIthSlaveID(i);
    mpi->send(MSG_QUIT,slaveID);
  }

  cout<<"MASTER waiting for slaves to terminate"<<endl;
  for(int i=0 ; i<numSlaves ; ++i) {
    MpiFixedMsg *msg=mpi->waitForMessage();
    cout<<"MASTER received another slave termination signal"<<endl;
    cout<<"    from: "<<msg->getSource()<<" tag="<<msg->getTag()<<endl;
    delete msg;
  }
  cout<<"MASTER terminating"<<endl;

  return 0;
}



int Application::slave(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=0)
    throw String("test-mpi <no args>");

  // Serve the master...
  bool done=false;
  while(!done) {
    cout<<"slave #"<<mpi->getProcessID()<<" calling MPI_Probe"<<endl;
    MpiFixedMsg *message=mpi->waitForMessage();
    switch(message->getTag()) 
      {
      case MSG_EXAMPLE_1: 
	{
	  cout<<"slave "<<mpi->getProcessID()<<" received MSG_EXAMPLE_1 from process "<<message->getSource()<<endl;
	  int x;
	  double y;
	  (*message)>>x>>y;
	  cout<<"x="<<x<<", y="<<y<<endl;
	}
	break;
      case MSG_QUIT:
	cout<<"slave "<<mpi->getProcessID()<<" quitting"<<endl;
	done=true;
	break;
      default:
	cout<<"SLAVE: unknown message!"<<endl;
	break;
      }
    delete message;
  }

  cout<<"slave #"<<mpi->getProcessID()<<"returning QUIT to master"<<endl;
  mpi->send(MSG_QUIT,0);

  return 0;
}


*/
