#ifndef QMCPLUSPLUS_CEIMC_H
#define QMCPLUSPLUS_CEIMC_H

#include "IonMovers/IonDriverBase.h"

namespace qmcplusplus
{


class CEIMC: public IonDriverBase
{
  public:
    CEIMC(){}; //default constructor.
    ~CEIMC(){};

    void test();
    //CEIMC(const CEIMC&); //copy constructor
    //CEIMC operator=(const CEIMC&);
  
};

}
#endif
