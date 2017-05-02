#include "IonMovers/IonUpdateBase.h"
#include "IonMovers/CEIMC/CEIMCUpdateAll.h"
#include "IonMovers/BOSurfaceBase.h"
#include "IonMovers/IonSystem.h"

namespace qmcplusplus
{
CEIMCUpdateAll::CEIMCUpdateAll(IonSystem* i, BOSurfaceBase* bo)
:IonUpdateBase(i,bo)
{
}

bool CEIMCUpdateAll::advanceIons()
{
  std::cout<<" ceimcupdateall::advanceIons()\n";
  return 0;
}

}
