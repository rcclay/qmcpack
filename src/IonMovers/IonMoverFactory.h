#ifndef QMCPLUSPLUS_ION_MOVER_FACTORY_H
#define QMCPLUSPLUS_ION_MOVER_FACTORY_H

namespace qmcplusplus
{
class IonMoverFactory
{
  public:
    IonMoverFactory();
    ~IonMoverFactory();
  
    bool parse(xmlNodePtr cur);
    bool execute(xmlNodePtr cur);

  private:
    bool parse_cell(xmlNodePtr cur);
    bool parse_ionset(xmlNodePtr cur);
    bool parse_iondriver(xmlNodePtr cur);
    bool parse_bosurface(xmlNodePtr cur);

   // BOSurfaceBase* bosurface;
   // ParticleSet* ion0;
   // IonDriver* 
  
};
}

#endif
