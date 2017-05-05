#ifndef QMCPLUSPLUS_BOSURFACE_LJ_H
#define QMCPLUSPLUS_BOSURFACE_LJ_H

#include "IonMovers/BOSurfaceBase.h"
namespace qmcplusplus
{

class ParticleSet;

class LJPot: public BOSurfaceBase
{
  public: 
    //epsilon is prefactor, sig is the reduction in r, err is a fictious error
    LJPot(const double sig, const double eps, const double err );
    ~LJPot();
    bool E(const ParticleSet& P, double& e, double &err);
    bool dE(const ParticleSet& P0, const ParticleSet& P1, double& e, double& err);

  private:
    double epsilon;
    double sigma;
    double error;
};

}
#endif

