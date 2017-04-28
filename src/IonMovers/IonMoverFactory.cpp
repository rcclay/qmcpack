#include "IonMovers/IonMoverFactory.h"
#include "IonMovers/BOSurfaceFactory.h"
#include "Configuration.h"
#include <string>


namespace qmcplusplus
{

IonMoverFactory::IonMoverFactory()
{
}

bool IonMoverFactory::parse(xmlNodePtr cur)
{
  IonSystem ionsystem;
  cur=cur->xmlChildrenNode;
  bool foundIons(false);
  bool foundCell(false);
  bool foundPES(false);

  while (cur!=NULL)
  {
    std::string nodename((const char*)cur->name);
    std::cout<<"  "<<nodename<<"\n";
    if(nodename=="simulationcell")
    {
      foundCell=ionsystem.parseCell(cur);
    }
    else if(nodename=="particleset")
    {
      foundIons=ionsystem.parseIons(cur);
    }
    else if(nodename=="bosurface")
    {
      //I need to work on the bosurfacefactory logic here.  
      app_log()<<"Initializing BOSurfaceBase\n";
      BOSurfaceFactory x;
      x.parse(cur);
      //for debugging, set to true.
      foundPES=true;
    }
    else;

    cur=cur->next;
  } 
  
  if(!foundCell)
  {
    app_log()<<"Error.  Unable to find a simulation cell in XML file.  Please check \
               that it exists and that the syntax is correct.\n";
  }

  if(!foundIons)
  {
    app_log()<<"Error.  Unable to find an ion ParticleSet in XML file.  Please check \
               that it exists and that the syntax is correct.\n";
  } 
  if(!foundPES)
  {
    app_log()<<"Error.  Unable to construct potential energy surface from XML file.  Please check \
               that it exists and that the syntax is correct.\n";
  } 
  //We need all the pieces for a run.
  if( foundCell && foundPES && foundIons )
  { 
    ionsystem.put(cur);
  }
  else
  {
    APP_ABORT("IonMoverFactory::parse(cur):  Unable to build IonMover from xml file.\n");
    return 1;
  }
  //all's good
  return 0; 
}

 
}
