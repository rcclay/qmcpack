#ifndef QMCPLUSPLUS_ION_DRIVER_BASE_H
#define QMCPLUSPLUS_ION_DRIVER_BASE_H

namespace qmcplusplus
{
class IonDriverBase
{
  public:
    IonDriverBase(){};
    ~IonDriverBase(){};
    IonDriverBase(IonDriverBase&);
    IonDriverBase operator=(IonDriverBase& u);

    virtual bool put(xmlNodePtr cur){}; //to parse options from XML
    virtual bool run(){};
    virtual bool update_ions(){};  //to update the ions
    virtual bool reset(){};
    virtual bool set_bosurface(BOSurfaceBase* bo){bosurface=bo};; 
  private:
    int numSteps; //number of steps
    double tau; //timestep
    double t; //temperature

    BOSurfaceBase* bosurface;
    //probably going to need QMCsystem1, QMCSystem2.  
     
}

}
#endif
