#include "IonMovers/IonUpdateBase.h"
#include "IonMovers/CEIMC/CEIMCUpdateAll.h"
#include "IonMovers/BOSurfaceBase.h"
#include "IonMovers/IonSystem.h"
#include "Configuration.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include <typeinfo>

namespace qmcplusplus
{

CEIMCUpdateAll::CEIMCUpdateAll(IonSystem* i, BOSurfaceBase* bo, RandomGenerator_t& m)
:IonUpdateBase(i,bo,m)
{
  int Natom = ions->getNumAtoms();
  SqrtTauOverM.resize(Natom);

//  Rnew.resize(Natom);
//  deltaR.resize(Natom);
//  G.resize(Natom);
}

bool CEIMCUpdateAll::advanceIons()
{
//  ParticleSetPos_t& ions->R
  ParticleSet* Pcur = ions->getCurParticleSet();
  ParticleSet* Ptmp = ions->getTmpParticleSet();
 // Rnew=P->R;  //Backup R with the old coordinates.
  
  //Generate a random 3N-dim vector.
  ions->Rnew=ions->R;

  makeGaussRandomWithEngine(ions->deltaR,Rng);
  ions->makeMove(ions->Rnew,ions->deltaR,SqrtTauOverM);
  
  ParticlePos_t dR(ions->R);
  dR=ions->Rnew- ions->R;
  
  app_log()<<"deltar = "<<dR<<std::endl;
  app_log()<<"  Rold = "<<ions->R<<std::endl;
  app_log()<<"  Rnew = "<<ions->Rnew<<std::endl;
  app_log()<<" ceimcupdateall::accept move()\n";
  double dE;
  double dE_err;
  bosurface->dE(*Pcur, *Ptmp, dE, dE_err);
  
  app_log()<<"   dE = "<<dE<<" err = "<<dE_err<<std::endl;
  ions->acceptMove();
   
   
  return 0;
}

bool CEIMCUpdateAll::resetRun()
{
  int Natom = ions->getNumAtoms();
  ParticleSet* P = ions->getCurParticleSet();

  SqrtTauOverM.resize(Natom);
  app_log()<<"Reset CEIMCUpdateAll\n";
  for(int iat=0; iat<Natom; iat++)
  {
    SqrtTauOverM[iat]=std::sqrt(tau/RealType(P->Mass[iat]));
    app_log()<<"  "<<iat<<" "<<P->Mass[iat]<<" "<<SqrtTauOverM[iat]<<std::endl;
  }  

}

}
