#ifndef QMCPLUSPLUS_ION_MOVER_FACTORY_H
#define QMCPLUSPLUS_ION_MOVER_FACTORY_H

#include "OhmmsData/libxmldefs.h"
#include "IonMovers/IonSystem.h"

namespace qmcplusplus
{

class BOSurfaceBase;
class IonSystem;
class IonDriverBase;

class IonMoverFactory
{
  public:
    IonMoverFactory();
    ~IonMoverFactory(){};
  
    bool parse(xmlNodePtr cur);
    bool execute();

  private:

    BOSurfaceBase* bosurface;
    IonSystem* ions;
    IonDriverBase* iondriver;

    xmlNodePtr boxml;
    xmlNodePtr iondriverxml;
    xmlNodePtr ionsysxml;  
  
};
}

#endif
