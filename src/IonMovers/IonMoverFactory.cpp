#include "IonMovers/IonMoverFactory.h"
#include "IonMovers/BOSurfaceFactory.h"
#include "Configuration.h"
#include <string>


namespace qmcplusplus
{
bool IonMoverFactory::parse(xmlNodePtr cur)
{
  cur=cur->xmlChildrenNode;
  while (cur!=NULL)
  {
    std::string nodename((const char*)cur->name);
    std::cout<<"  "<<nodename<<"\n";
    if(nodename=="simulationcell")
    {
      app_log()<<"I would build the simulation cell here.\n";
    }
    else if(nodename=="particleset")
    {
      app_log()<<"I would build the particle set here.\n"; 
    }
    if(nodename=="bosurface")
    {
      app_log()<<"Initializing BOSurfaceBase\n";
      BOSurfaceFactory x;
      x.parse(cur);
   //   BOSurfaceBase* pes= new BOSurfaceBase();
    }
    cur=cur->next;
  } 
  //all's good
  return 0; 
}

 
}
