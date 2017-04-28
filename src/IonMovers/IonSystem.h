#ifndef QMCPLUSPLUS_ION_SYSTEM_H
#define

namespace qmcplusplus
{
class IonSystem
{
  public:
    IonSystem();
    ~IonSystem();
    
    virtual bool parse_cell(xmlNodePtr cur);
    virtual bool parse_ions(xmlNodePtr cur);
    
    virtual void update_all();
  private:
    ParticleSet* P; //for the ions
    ParticleSet::ParticleLayout_t* SimulationCell;

    
    
}
}
#endif

