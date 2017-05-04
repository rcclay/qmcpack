#ifndef QMCPLUSPLUS_ION_MOVER_FACTORY_H
#define QMCPLUSPLUS_ION_MOVER_FACTORY_H

#include "OhmmsData/libxmldefs.h"
#include "IonMovers/IonSystem.h"
#include "OhmmsApp/RandomNumberControl.h"

namespace qmcplusplus
{

class BOSurfaceBase;
class IonSystem;
class IonDriverBase;
class RandomNumberControl;

class IonMoverFactory
{
  public:
    IonMoverFactory(RandomNumberControl& rand);
    ~IonMoverFactory(){};
  
    bool parse(xmlNodePtr cur);
    bool execute();

  private:
    RandomNumberControl& myRandControl;
 
    BOSurfaceBase* bosurface;
    IonSystem* ions;
    IonDriverBase* iondriver;

    xmlNodePtr boxml;
    xmlNodePtr iondriverxml;
    xmlNodePtr ionsysxml;  
  
};
}

#endif
