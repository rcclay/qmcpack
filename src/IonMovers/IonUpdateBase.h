#ifndef ION_UPDATE_BASE_H
#define ION_UPDATE_BASE_H

#include "Configuration.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "IonMovers/physicalconstants.h"

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
        ions(i),bosurface(bo),Rng(m),tau(0.0),beta(0.0)
    {};
    ~IonUpdateBase(){};

    virtual bool advanceIons()=0;
    virtual bool resetRun()=0;
    virtual bool finalizeRun()=0;

    inline RealType getTau(){return tau;};
    inline void setTau(RealType t){tau=t;};
    inline void setTemperature(RealType T){ t=T*KELVIN_TO_HARTREE;
                                            if(t!=0) beta=1.0/t;};

  protected:
    IonSystem* ions;
    BOSurfaceBase* bosurface;
    
    RandomGenerator_t Rng;
   
    RealType tau;
    RealType t;
    RealType beta; 

    std::vector<RealType> SqrtTauOverM;
};

}

#endif
