#ifndef QMCPLUSPLUS_ION_DRIVER_BASE_H
#define QMCPLUSPLUS_ION_DRIVER_BASE_H

#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
namespace qmcplusplus
{
//Forward declaration
class BOSurfaceBase;
class IonSystem;

class IonDriverBase
{
  public:
    IonDriverBase();
    IonDriverBase(IonSystem* ions, BOSurfaceBase* bosurface);
    ~IonDriverBase(){};
    //IonDriverBase(IonDriverBase&);
    //IonDriverBase operator=(IonDriverBase& u);

    virtual bool put(xmlNodePtr cur); //to parse options from XML
//    virtual bool run(){};
//    virtual bool updateIons(){};  //to update the ions
//    virtual bool reset(){};
//    virtual bool setBOSurface(BOSurfaceBase* bo){bosurface=bo};; 
    virtual bool reset(){return 0;};
    virtual bool run()=0;

  protected:
    int numSteps; //number of steps
    double tau; //timestep
    double t; //temperature

    BOSurfaceBase* bosurface;
    IonSystem* ions;
    //probably going to need QMCsystem1, QMCSystem2.  
 
    OhmmsAttributeSet pAttrib;
    ParameterSet m_param;  

    bool moveall;  
};

}
#endif
