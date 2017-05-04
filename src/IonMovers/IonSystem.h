#ifndef QMCPLUSPLUS_ION_SYSTEM_H
#define QMCPLUSPLUS_ION_SYSTEM_H

#include "OhmmsData/libxmldefs.h"
#include "Particle/ParticleSet.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "Configuration.h"
#include "LongRange/StructFact.h"
namespace qmcplusplus
{

class IonSystem: public QMCTraits
{
  public:
    typedef ParticleSet::ParticleLayout_t CellType;
    typedef ParticleSet::ParticlePos_t ParticlePos_t;
    typedef ParticleSet::ParticleGradient_t ParticleGradient_t;
    typedef ParticleSet::SingleParticlePos_t SingleParticlePos_t;
    IonSystem();
    ~IonSystem(){};
    
    virtual bool parseCell(xmlNodePtr cur);
    virtual bool parseIons(xmlNodePtr cur);
    virtual bool put(xmlNodePtr cur);
    
    virtual void updateAll(){};

    virtual bool makeMove(ParticlePos_t& Rnew,
                          const ParticlePos_t& dR,
                          const std::vector<RealType>& taueff);
    virtual bool acceptMove();
    inline int getNumAtoms(){return Natom;};
    inline ParticleSet* getCurParticleSet(){if(P0.size()>0) return &P0[curID]; else return NULL;};
    inline ParticleSet* getTmpParticleSet(){if(P0.size()>0) return &P0[tmpID]; else return NULL;};
    //accessor function
    
    ParticlePos_t R; //The current position
    ParticlePos_t Rnew;  //Container for proposed position.
    ParticlePos_t deltaR; //Container for proposed dR
    ParticleGradient_t G; //Container for gradient

  private:
    std::vector<ParticleSet> P0; //for the ions
    
    int Natom; 
    int curID;
    int tmpID;

    CellType simulation_cell;
 //Might want to add another simulation cell...  
    //We are going to assume that the tile matrix at the IonDriver is
    //the identity.  We will potentially allow a different tile matrix
    //for the QMC driver to handle electronic & other finite size effects.
    Tensor<int,OHMMS_DIM> TileMatrix;
    
    
};
}
#endif

