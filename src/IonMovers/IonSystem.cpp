#include "IonMovers/IonSystem.h"
#include "Utilities/ProgressReportEngine.h"
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/OhmmsParameter.h"
#include "ParticleIO/XMLParticleIO.h"
namespace qmcplusplus
{

IonSystem::IonSystem():
simulation_cell(), P(), TileMatrix()
{
  app_log()<<"--Calling IonSystem constructor\n";
  TileMatrix.diagonal(1);
}

bool IonSystem::parseCell(xmlNodePtr cur)
{
  ReportEngine PRE("ParticleSetPool","putLattice");
  
  LatticeParser a(simulation_cell);
  bool success=a.put(cur);
  simulation_cell.print(app_log());
  
  return success;

}

bool IonSystem::parseIons(xmlNodePtr cur)
{
  ReportEngine PRE("ParticleSetPool","put");
  //const ParticleSet::ParticleLayout_t* sc=DistanceTable::getSimulationCell();
  //ParticleSet::ParticleLayout_t* sc=0;
  std::string id("e");
  std::string role("none");
  std::string randomR("no");
  std::string randomsrc;
  OhmmsAttributeSet pAttrib;
  pAttrib.add(id,"id");
  pAttrib.add(id,"name");
  pAttrib.add(role,"role");
  pAttrib.add(randomR,"random");
  pAttrib.add(randomsrc,"randomsrc");
  pAttrib.add(randomsrc,"random_source");
  pAttrib.put(cur);
  //backward compatibility
  if(id == "e" && role=="none")
    role="MC";
  
  app_log() << "  Creating " << id << " particleset" << std::endl;

  //THIS WILL BE CALLED AFTER INITIALIZATION OF SIMULATION CELL.  
//    if(SimulationCell)
//    {
//      app_log() << "  Initializing the lattice of " << id << " by the global supercell" << std::endl;
//      pTemp->Lattice.copy(*SimulationCell);
//    }
//    myPool[id] = pTemp;
   XMLParticleParser pread(P,TileMatrix);
   bool success = pread.put(cur);
    //if random_source is given, create a node <init target="" soruce=""/>
//   if(randomR=="yes" && !randomsrc.empty())
//   {
//     xmlNodePtr anode = xmlNewNode(NULL,(const xmlChar*)"init");
//     xmlNewProp(anode,(const xmlChar*)"source",(const xmlChar*)randomsrc.c_str());
//     xmlNewProp(anode,(const xmlChar*)"target",(const xmlChar*)id.c_str());
//     randomize_nodes.push_back(anode);
//   }
    P.setName(id);
    app_log() << P.getName() << std::endl;
   
    return success;
}

bool IonSystem::put(xmlNodePtr cur)
{
  P.Lattice.copy(simulation_cell);
  P.get(app_log());

  return 0;
}

}
