#ifndef QMCPLUSPLUS_ION_SYSTEM_H
#define QMCPLUSPLUS_ION_SYSTEM_H

#include "OhmmsData/libxmldefs.h"

namespace qmcplusplus
{
class IonSystem
{
  public:
    IonSystem();
    ~IonSystem();
    
    virtual bool parseCell(xmlNodePtr cur);
    virtual bool parseIons(xmlNodePtr cur);
    virtual bool put(xmlNodePtr cur);
    
    virtual void updateAll();
  private:
    ParticleSet P; //for the ions
    ParticleSet::ParticleLayout_t simulation_cell;

    
    
}
}
#endif

