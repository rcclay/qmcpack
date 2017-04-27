#include "IonMovers/BOSurfaceFactory.h"
#include <iostream>
#include <string>
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/OhmmsElementBase.h"

namespace qmcplusplus
{

bool BOSurfaceFactory::parse(xmlNodePtr cur)
{
  std::string curName((const char*)cur->name);
  std::string imtype("none");
  
  OhmmsAttributeSet aAttrib;
  aAttrib.add(imtype,"type");
  aAttrib.put(cur);

  if(imtype=="classicalpot")
  {
    std::cout<<" Found a classical potential\n";
  }
  else if(imtype=="qmc")
  {
    std::cout<<" Found a qmc specification\n";
    std::cout<<" Need to specify a Particle set, Hamiltonian, Wavefunction, and QMCDriver";
  }
  else if(imtype=="dft")
  {
    std::cout<<" Found a DFT specification\n";
  }
  else if(imtype=="none")
  {
    std::cout<<" Error:  No BOSurface specified, or type=\"none\"\n";
  }
  else
  {
    std::cout<<" Error:  Unrecognized option\n";
    return 0;
  }
  
  return 1;
}

}
