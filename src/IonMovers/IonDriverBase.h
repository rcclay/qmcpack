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

    virtual bool put()=0; to parse options from XML
    virtual bool update_ions()=0;  //to update the ions
    
  private:
    //probably going to need QMCsystem1, QMCSystem2.  
     
}

}
#endif
