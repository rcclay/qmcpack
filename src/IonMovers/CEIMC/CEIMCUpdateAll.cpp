#include "IonMovers/IonUpdateBase.h"
#include "IonMovers/CEIMC/CEIMCUpdateAll.h"
#include "IonMovers/BOSurfaceBase.h"
#include "IonMovers/IonSystem.h"
#include "Configuration.h"
#include "OhmmsApp/RandomNumberControl.h"

namespace qmcplusplus
{
CEIMCUpdateAll::CEIMCUpdateAll(IonSystem* i, BOSurfaceBase* bo, RandomNumberControl& m)
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

//  app_log()<<"dr = "<<dr<<std::endl;
  std::cout<<" ceimcupdateall::advanceIons()\n";
   
  return 0;
}

}
