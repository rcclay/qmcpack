#include "IonMovers/BOSurfaces/LJPot.h"
#include "IonMovers/BOSurfaceBase.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Message/Communicate.h"
#include "Configuration.h"

namespace qmcplusplus
{

LJPot::LJPot(double sig, double eps, double err)
: sigma(sig),epsilon(eps), error(err)
{

}

bool LJPot::E(const ParticleSet& P, double & E, double & err)
{
  const DistanceTableData * d_aa(P.DistTables[0]);
  int Natom=P.getTotalNum();
  E=0;
  err=error;
  double d=0;

  for(int ipart=0; ipart<Natom; ipart++)
  {
    for(int nn=d_aa->M[ipart],jpart=ipart+1; nn<d_aa->M[ipart+1]; nn++,jpart++)
    {
      app_log()<<"ipart="<<ipart<<" nn="<<nn<<std::endl;
      app_log()<<" d="<<d_aa->r(nn)<<" eps="<<epsilon<<" sig="<<sigma<<std::endl;
      d = d_aa->r(nn);
      E += epsilon*(std::pow(sigma/d,12.0)-std::pow(sigma/d,6.0));
    }
  }
}

bool LJPot::dE(const ParticleSet& P1, const ParticleSet& P2, double &dE, double & derr)
{
  app_log()<<" LJPot::dE called\n";
  double E1(0);
  double E2(0);
  double err1(error);
  double err2(error);

  dE=0;
  derr=0;
  
  app_log()<<"Calling E1\n";
  E(P1,E1,err1);
  app_log()<<"Calling E2\n";
  E(P2,E2,err2);
  app_log()<<" E1 = "<<E1<<" E2 = "<<E2<<std::endl;
  app_log()<<"done.  combine\n";
  dE=E2-E1;
  derr=std::sqrt(err1*err1+err2*err2);
  
  return 0;  
  
}
}
