
#include "GainLossType.H"
#include "BOOM/Exceptions.H"
using namespace std;

ostream &operator<<(ostream &os,GainLossType t)
{
  switch(t)
    {
    case GLT_VOID: os<<"GLT_VOID"; break;
    case GLT_GAIN: os<<"GLT_GAIN"; break;
    case GLT_LOSS: os<<"GLT_LOSS"; break;
    case GLT_RETENTION: os<<"GLT_RETENTION"; break;
    }
}



String toString(GainLossType t)
{
  switch(t)
    {
    case GLT_VOID: return "GLT_VOID";
    case GLT_GAIN: return "gain";
    case GLT_LOSS: return "loss";
    case GLT_RETENTION: return "retention";
    }
  INTERNAL_ERROR;
}



