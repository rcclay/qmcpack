#ifndef QMCPLUSPLUS_CEIMC_HELPER_FUNCTIONS_H
#define QMCPLUSPLUS_CEIMC_HELPER_FUNCTIONS_H

namespace qmcplusplus
{
  //takes in an energy errorbar and inverse temperature (make sure units are consistent)
  //returns penalty in limit of infinite samples.

  inline double ceimc_penalty(const double err, const double beta)
  {
    double berr=beta*err;
    double bevar=berr*berr;
    return bevar/2.0;
  }

  //takes in energy error bar and inverse temperature.  Expands up to 2rd order in number of samples
  inline double ceimc_penalty(const double err, const double beta, const int nsamps)
  {
    double berr=beta*err;
    double bevar=berr*berr;
 
    return bevar/2.0 + bevar*bevar/(4.0*double(nsamps+1))+bevar*bevar*bevar/(3.0*double(nsamps+1)*double(nsamps+3));

  }
}
#endif
