#ifndef QMC_SYSTEM_BASE_H
#define QMC_SYSTEM_BASE_H

namespace qmcplusplus
{

class QMCSystem: public IonSystem
{ 
  public:
    QMCSystem();
    ~QMCSystem();
    QMCSystem(QMCSystem& q);
    QMCSystem operator=(QMCSystem& q);

    setPWH(ParticleSet& P, TrialWavefunction& Psi, Hamiltonian& H);
    
  private:
    //QMCSystem should also include MPI information 
    ParticleSet::ParticleLayout_t* SimulationCell;

    Hamiltonian* H;
    TrialWavefuntioni* Psi;
    ParticleSet* P; //for the ions
    ParticleSet::ParticleLayout_t* SimulationCell;
    
};


}
#endif
