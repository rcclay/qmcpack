#ifndef ION_UPDATE_BASE_H
#define ION_UPDATE_BASE_H

#include "Configuration.h"
#include "OhmmsApp/RandomNumberControl.h"

namespace qmcplusplus
{

class IonSystem;
class BOSurfaceBase;

//We derive from QMCTraits to get access to all the nice data containers and types
//particular to QMCPACK.
class IonUpdateBase: public QMCTraits
{
  public:
    IonUpdateBase(IonSystem* i,BOSurfaceBase* bo,RandomNumberControl& m):
        ions(i),bosurface(bo),Rng(m)
    {};
    ~IonUpdateBase(){};

    virtual bool advanceIons()=0;
   
  protected:
    IonSystem* ions;
    BOSurfaceBase* bosurface;
    
    RandomNumberControl& Rng;

};

}

#endif
