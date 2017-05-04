#ifndef QMCPLUSPLUS_BOUPDATE_BASE_H
#define QMCPLUSPLUS_BOUPDATE_BASE_H

#include "OhmmsData/libxmldefs.h"
#include "Configuration.h"

namespace qmcplusplus
{

class BOSurfaceBase: public QMCTraits 
{
//  typedef std::vector<std::vector<double> > ParticlePos_t;
  public:
    BOSurfaceBase(){};
    ~BOSurfaceBase(){};  //destructor

//    BOSurfaceBase(const BOSurfaceBase& u){}; //copy constructor
//    BOSurfaceBase& operator=(const BOSurfaceBase& u){}; //copy assignment
    
    virtual bool put(xmlNodePtr){return 0;}
    virtual bool initialize(){return 0;}
    virtual bool reset(){return 0;}
    virtual bool compute(){return 0;};
    virtual bool update(){return 0;};

//    virtual bool E(const ParticlePos_t& R1, double& dE, double& err){return 0;};
//    virtual bool dE(const ParticlePos_t& R0, const ParticlePos_t& R1, double& dE, double& err){return 0;};
//    virtual bool dE(const ParticlePos_t& R0, const ParticlePos_t& R1, int iat, double& dE, double& err){return 0;};
//    virtual bool F(const ParticlePos_t& R1, ParticlePos_t& F){return 0;};
  private:
    //nothing for now
};

}

#endif
