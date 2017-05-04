#include "IonMovers/IonDriverBase.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
IonDriverBase::IonDriverBase()
: ions(NULL), bosurface(NULL), numSteps(0), tau(0), t(0)
{
  m_param.add(t,"temperature","double");
  m_param.add(tau,"timestep","double");
  m_param.add(numSteps,"steps","int");
}

IonDriverBase::IonDriverBase(IonSystem* i, BOSurfaceBase* bo, RandomGenerator_t& m)
: ions(i), bosurface(bo), numSteps(0), tau(0), t(0),Rng(m)
{
  m_param.add(t,"temperature","double");
  m_param.add(tau,"timestep","double");
  m_param.add(numSteps,"steps","int");
}

bool IonDriverBase::put(xmlNodePtr cur)
{
  m_param.put(cur); 
  return 1; 
}

}
