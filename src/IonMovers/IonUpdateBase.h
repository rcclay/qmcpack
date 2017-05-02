#ifndef ION_UPDATE_BASE_H
#define ION_UPDATE_BASE_H

namespace qmcplusplus
{

class IonSystem;
class BOSurfaceBase;

class IonUpdateBase
{
  public:
    IonUpdateBase(IonSystem* i,BOSurfaceBase* bo):ions(i),bosurface(bo){};
    ~IonUpdateBase(){};

    virtual bool advanceIons()=0;
   
  protected:
    IonSystem* ions;
    BOSurfaceBase* bosurface;


};

}

#endif
