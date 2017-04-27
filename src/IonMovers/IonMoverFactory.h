#ifndef QMCPLUSPLUS_ION_MOVER_FACTORY_H
#define QMCPLUSPLUS_ION_MOVER_FACTORY_H

namespace qmcplusplus
{
class IonMoverFactory
{
  IonMoverFactory();
  ~IonMoverFactory();
  
  bool parse(xmlNodePtr cur);
  
}
}

#endif
