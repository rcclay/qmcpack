#include "IonMovers/IonMoverFactory.h"
#include "IonMovers/BOSurfaceFactory.h"
#include "IonMovers/BOSurfaceBase.h"
#include "IonMovers/IonDriverBase.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "IonMovers/CEIMC/CEIMC.h"

#include "Configuration.h"
#include <string>
#include <iostream>

namespace qmcplusplus
{

IonMoverFactory::IonMoverFactory(RandomNumberControl& m)
: boxml(NULL), ionsysxml(NULL), iondriverxml(NULL), iondriver(NULL),bosurface(NULL),
ions(0),myRandControl(m)
{
  std::cout<<"IonMoverFactory constructor\n";
}

bool IonMoverFactory::parse(xmlNodePtr cur)
{
  ions = new IonSystem();
  bosurface = new BOSurfaceBase();

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
      foundCell=ions->parseCell(cur);
    }
    else if(nodename=="particleset")
    {
      foundIons=ions->parseIons(cur);
      ionsysxml=cur;
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
    else if(nodename=="ionmover")
    {
      //we're going to parse this after initial data is set up and consistent.
      //for now, just store the xml node.
      iondriverxml=cur;
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
    ions->put(cur);
    
  }
  else
  {
    APP_ABORT("IonMoverFactory::parse(cur):  Unable to build IonMover from xml file.\n");
    return 1;
  }

  ///////////////////////////////////////////////////////////////////////////////
  //Now that we have the system set up.  build the factory.
  ///////////////////////////////////////////////////////////////////////////////
  if (!iondriverxml)
  {
    APP_ABORT("IonMoverFactory::parse(cur):  Missing an ion driver.  Please check <ionmover type=...\n");
    return 1;
  }
  std::string ionmovertype("none");
  OhmmsAttributeSet imAttribute;
  imAttribute.add(ionmovertype,"type");
  imAttribute.put(iondriverxml);
  

  if (ionmovertype=="ceimc")
  {
    CEIMC* ceimc = new CEIMC(ions,bosurface,*(RandomNumberControl::Children[0]));
    ceimc->put(iondriverxml);
    iondriver = static_cast<IonDriverBase*>(ceimc);
    app_log()<<" iondriver="<<iondriver<<std::endl; 
  }
  else
  {
    APP_ABORT("IonMoverFactory::parse(cur):  Unrecognized ion driver.  Please check <ionmover type=...\n");
    return 1; 
  }  
   
  //all's good
  return 0; 
}

bool IonMoverFactory::execute()
{
  app_log()<<" In IonMoverFactory::execute()\n";
  app_log()<<" iondriver = "<<iondriver<<std::endl;
  iondriver->run();
  return 0;
}
 
}
