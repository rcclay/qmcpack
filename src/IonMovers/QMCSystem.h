#ifndef QMC_SYSTEM_BASE_H
#define QMC_SYSTEM_BASE_H

namespace qmcplusplus
{

class QMCSystemBase
{ 
  public:
    QMCSystemBase();
    ~QMCSystemBase();
    QMCSystemBase(QMCSystemBase& q);
    QMCSystemBase operator=(QMCSystemBase& q);

    setPWH(ParticleSet& P, TrialWavefunction& Psi, Hamiltonian& H);
    
  private:
    //QMCSystemBase should also include MPI information 
    ParticleSet::ParticleLayout_t* SimulationCell;

    Hamiltonian* H;
    TrialWavefuntioni* Psi;
    ParticleSet* P; //for the ions
    ParticleSet::ParticleLayout_t* SimulationCell;
    
};


}
#endif
