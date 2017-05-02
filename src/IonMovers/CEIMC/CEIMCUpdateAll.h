#ifndef QMCPLUSPLUC_CEIMC_UPDATE_ALL_H
#define QMCPLUSPLUC_CEIMC_UPDATE_ALL_H

namespace qmcplusplus
{

class IonUpdateBase;
class IonSystem;
class BOSurfaceBase;

class CEIMCUpdateAll: public IonUpdateBase
{
  public:
    CEIMCUpdateAll(IonSystem* i, BOSurfaceBase* bo);
    ~CEIMCUpdateAll(){};
  
    bool advanceIons();

  protected:
  //Inherited from IonUpdateBase:
  //IonSystem* ions
  //BOSurfaceBase* bosurface.
    
};

}

#endif
