#ifndef QMC_SYSTEM_BASE_H
#define QMC_SYSTEM_BASE_H

namespace qmcplusplus
{

class QMCSystemBase: public X
{ 
  public:
    QMCSystemBase();
    ~QMCSystemBase();
    QMCSystemBase(QMCSystemBase& q);
    QMCSystemBase operator=(QMCSystemBase& q);

    setPWH(ParticleSet& P, TrialWavefunction& Psi, Hamiltonian& H);
    
  private:
    //We need to think about where we want P,Psi,H data to live.  Objects here?  Pointers here?  What?
    Hamiltonian H;
    TrialWavefuntion Psi;
    ParticleSet P;

    
}


}
#endif
