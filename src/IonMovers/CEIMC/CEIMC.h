#ifndef QMCPLUSPLUS_CEIMC_H
#define QMCPLUSPLUS_CEIMC_H

#include "IonMovers/IonDriverBase.h"

namespace qmcplusplus
{

//forward declaration
class IonUpdateBase;

class CEIMC: public IonDriverBase
{
  public:
    CEIMC(IonSystem* i,BOSurfaceBase*,RandomNumberControl& m); //default constructor.
    ~CEIMC(){};

    void test();
    bool put(xmlNodePtr cur);
    bool run();
  private:
    
    IonUpdateBase* mover;

    bool usedrift;
    bool prereject;
  protected:
  //These derive from IonDriverBase.

    //BOSurfaceBase bosurface;
    //IonSystem ions;   

  
};

}
#endif
