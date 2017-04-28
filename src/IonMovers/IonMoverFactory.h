#ifndef QMCPLUSPLUS_ION_MOVER_FACTORY_H
#define QMCPLUSPLUS_ION_MOVER_FACTORY_H

#include "OhmmsData/libxmldefs.h"
#include "IonMovers/IonSystem.h"

namespace qmcplusplus
{
class IonMoverFactory
{
  public:
    IonMoverFactory();
    ~IonMoverFactory(){};
  
    bool parse(xmlNodePtr cur);
    bool execute(xmlNodePtr cur){return 0;};

  private:
    bool parseCell(xmlNodePtr cur){return 0;};
    bool parseIons(xmlNodePtr cur){return 0;};
    bool parseIonDriver(xmlNodePtr cur){return 0;};
    bool parseBOSurface(xmlNodePtr cur){return 0;};

   // BOSurfaceBase* bosurface;
   // ParticleSet* ion0;
   // IonDriver* 
  
};
}

#endif
