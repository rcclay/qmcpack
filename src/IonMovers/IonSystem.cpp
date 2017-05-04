#include "IonMovers/IonSystem.h"
#include "Utilities/ProgressReportEngine.h"
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/OhmmsParameter.h"
#include "ParticleIO/XMLParticleIO.h"

namespace qmcplusplus
{

IonSystem::IonSystem():
simulation_cell(), TileMatrix(),curID(0),tmpID(1),Natom(0)
{
  P0.resize(0);
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
  P0.resize(1);
  XMLParticleParser pread(P0[0],TileMatrix);
  bool success = pread.put(cur);
   
  //override whatever the user calls the particle set.
  P0[0].setName("ionset0");
  app_log() << P0[0].getName() << std::endl;

    //Now to store the particle set.    
    //I'm assuming that adding a vector deep copies.  Hopefully I'm not mistaken.

     
    //Now allocate temporary containers for positions,gradients, etc.
  Natom=P0[0].getTotalNum();
  R.resize(Natom);
  Rnew.resize(Natom);
  deltaR.resize(Natom);
  G.resize(Natom);

    //Lastly, store the particle set.  
    return success;
}

bool IonSystem::put(xmlNodePtr cur)
{
  P0[0].Lattice.copy(simulation_cell);
  P0[0].get(app_log());

  //We wait this long to initialize the secondary particleset
  //because we want the lattice to get copied over first.  
  //Now we finish.
  ParticleSet P(P0[0]);
  P.setName("ionset1");
  
  P0.push_back(P);

  R=P0[0].R;
  return 0;
}

bool IonSystem::makeMove(ParticlePos_t& Rnew, const ParticlePos_t& dR,
                         const std::vector<RealType>& dt)
{
  //In this section, we do the following:
  //1.)  Take the random vector dR and the info on masses + time steps.
  //2.)  Using the current position defined in Rnew, create a new ion configuration.
  ///          --by adding sqrt(tau/m)*xi to R0.  While doing this, we check to 
  ///            make sure the move is not pathological
  //3.)  We update the temporary particleset position to Rnew, and update the
  //       distance tables and SK.
  //4.)  finish.

  //Grab the "current" and "proposed" particle sets.
  ParticleSet& P(P0[curID]);
  ParticleSet& Pnew(P0[tmpID]);

  P.activePtcl=-1;
  // Now propose the move based on the current particleset.  Check all moves.
  if (P.UseBoundBox)
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
    {
      SingleParticlePos_t displ(dt[iat]*deltaR[iat]);
      if (P.Lattice.outOfBound(P.Lattice.toUnit(displ)))
        return false;
      SingleParticlePos_t newpos(Rnew[iat]+displ);
      if (!P.Lattice.isValid(P.Lattice.toUnit(newpos)))
        return false;
      Rnew[iat]=newpos;
    }
  }
  else
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
      Rnew[iat]=Rnew[iat]+dt[iat]*deltaR[iat];
  }
  //Update R in the new particle set.
  Pnew.R=Rnew;

  //Update stuffs.
  for (int i=0; i< Pnew.DistTables.size(); i++)
    Pnew.DistTables[i]->evaluate(Pnew);
  if (Pnew.SK)
    Pnew.SK->UpdateAllPart(Pnew);
  //every move is valid
  return true; 
}

bool IonSystem::acceptMove()
{
  //Swap tmpID and curID.
  int tmp;
  tmp=tmpID;
  tmpID=curID;
  curID=tmp;

  //Update the R and Rnew vectors.  Sure we recopy, but we'll make this more efficient when
  //a profiler says its a bottleneck.

  ParticleSet* Pnew = getCurParticleSet();
  R=Rnew=Pnew->R;
  
}

}
