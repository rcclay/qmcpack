#ifndef QMCPLUSPLUS_ION_SYSTEM_H
#define QMCPLUSPLUS_ION_SYSTEM_H

#include "OhmmsData/libxmldefs.h"
#include "Particle/ParticleSet.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "Configuration.h"

namespace qmcplusplus
{

class IonSystem: public QMCTraits
{
  public:
    typedef ParticleSet::ParticleLayout_t CellType;
    IonSystem();
    ~IonSystem(){};
    
    virtual bool parseCell(xmlNodePtr cur);
    virtual bool parseIons(xmlNodePtr cur);
    virtual bool put(xmlNodePtr cur);
    
    virtual void updateAll(){};

    //accessor function
    ParticleSet* getIonSet(){return &P;};

  private:
    ParticleSet P; //for the ions
    CellType simulation_cell;

    //We are going to assume that the tile matrix at the IonDriver is
    //the identity.  We will potentially allow a different tile matrix
    //for the QMC driver to handle electronic & other finite size effects.
    Tensor<int,OHMMS_DIM> TileMatrix;
    
    
};
}
#endif

