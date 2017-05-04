#include "Message/Communicate.h"
#include "IonMovers/CEIMC/CEIMC.h" 
#include "IonMovers/IonDriverBase.h"
#include "OhmmsData/AttributeSet.h"
#include "IonMovers/IonUpdateBase.h"
#include "IonMovers/CEIMC/CEIMCUpdateAll.h"

#include <iostream>
#include <algorithm>
#include <string>

namespace qmcplusplus
{

CEIMC::CEIMC(IonSystem* i, BOSurfaceBase* bo, RandomNumberControl & m)
:IonDriverBase(i,bo,m),usedrift(false),prereject(false)
{

}
void CEIMC::test()
{
  std::cout<<" CEIMC::test()\n";
}  

bool CEIMC::put(xmlNodePtr cur)
{
  std::string usedriftstr("no");
  std::string prerejectstr("no");

  m_param.add(usedriftstr,"useDrift","string");
  m_param.add(usedriftstr,"usedrift","string");

  m_param.add(prereject,"prereject","string");
  m_param.put(cur);

  tolower(usedriftstr);
  tolower(prerejectstr);

  usedrift=(usedriftstr=="yes");
  prereject=(prerejectstr=="yes");

  if (!usedrift && !prereject)
  {
    std::cout<<"Yeehaw.  Would have built a CEIMCUpdateAllDriver here\n";
    mover=new CEIMCUpdateAll(ions,bosurface,Rng);
  }
 
  return 0; 
}

bool CEIMC::run()
{
  std::cout<<"I would run a CEIMC run here.  Here's the parameters I got:\n";
  std::cout<<" nsteps = "<<numSteps<<std::endl;
  std::cout<<" temperature = "<<t<<" K"<<std::endl;
  std::cout<<" ionic timestep = "<<tau<<" a.u.\n";
  for(int i=0; i<numSteps; i++)
  {
    std::cout<<i<<" step.  Advancing ions\n";
    mover->advanceIons();
  }

  return 0;
}

}
