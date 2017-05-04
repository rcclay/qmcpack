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
    IonUpdateBase(IonSystem* i,BOSurfaceBase* bo,RandomGenerator_t& m):
        ions(i),bosurface(bo),Rng(m),tau(0.0)
    {};
    ~IonUpdateBase(){};

    virtual bool advanceIons()=0;
    virtual bool resetRun()=0;

    inline RealType getTau(){return tau;};
    inline RealType setTau(RealType t){tau=t;};

  protected:
    IonSystem* ions;
    BOSurfaceBase* bosurface;
    
    RandomGenerator_t Rng;
   
    RealType tau; 
    std::vector<RealType> SqrtTauOverM;
};

}

#endif
