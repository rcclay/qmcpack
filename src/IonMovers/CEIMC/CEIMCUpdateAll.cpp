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
  ParticleSet* P = ions->getIonSet();
  int Natom = P->getTotalNum();
  Rnew.resize(Natom);
  deltaR.resize(Natom);
}

bool CEIMCUpdateAll::advanceIons()
{
  ParticleSet* P = ions->getIonSet();
  Rnew=P->R;
  app_log()<< "type of deltaR :  "<<typeid(deltaR).name()<<std::endl;
  app_log()<< "type of Rng :  "<<typeid(Rng).name()<<std::endl;
//  makeGaussWithRandomEngine(deltaR,Rng);
//  makeGaussRandom(deltaR);
  makeGaussRandomWithEngine(deltaR,Rng);
  app_log()<<"deltar = "<<deltaR<<std::endl;
  std::cout<<" ceimcupdateall::advanceIons()\n";
   
  return 0;
}

}
