#ifndef QMCPLUSPLUS_BO_SURFACE_FACTORY_H
#define QMCPLUSPLUS_BO_SURFACE_FACTORY_H

#include "OhmmsData/libxmldefs.h"

namespace qmcplusplus
{

class BOSurfaceFactory
{
  public:
    BOSurfaceFactory(){};
    ~BOSurfaceFactory(){};

    bool parse(xmlNodePtr cur);

  private:

};

}
#endif
