#ifndef QMCPLUSPLUC_CEIMC_UPDATE_ALL_H
#define QMCPLUSPLUC_CEIMC_UPDATE_ALL_H

#include "Particle/ParticleSet.h"

namespace qmcplusplus
{

class IonUpdateBase;
class IonSystem;
class BOSurfaceBase;
class RandomNumberControl;

class CEIMCUpdateAll: public IonUpdateBase
{
  typedef ParticleSet::ParticlePos_t ParticlePos_t;

  public:

    CEIMCUpdateAll(IonSystem* i, BOSurfaceBase* bo, RandomGenerator_t& m);
    ~CEIMCUpdateAll(){};
  
    bool advanceIons();

  protected:
    //Temporary data for ion moves;
    ParticlePos_t Rnew;
    ParticlePos_t deltaR;

  //Inherited from IonUpdateBase:
  
  //IonSystem* ions
  //BOSurfaceBase* bosurface.
  //RandomNumberControl Rng;
    
};

}

#endif
