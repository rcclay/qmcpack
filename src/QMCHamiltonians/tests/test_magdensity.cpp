//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Numerics/Quadrature.h"
#include "QMCHamiltonians/MagDensityEstimator.h"
//for wavefunction
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/SpinorSet.h"
#include "Particle/Walker.h"
//for nonlocal moves
#include "QMCHamiltonians/NonLocalTOperator.h"

#include "QMCWaveFunctions/ConstantSPOSet.h"
//for Hamiltonian manipulations.
#include "Particle/ParticleSet.h"
#include "LongRange/EwaldHandler3D.h"

//This is for the spinor test.
#include "QMCWaveFunctions/ElectronGas/FreeOrbital.h"

using std::string;

namespace qmcplusplus
{
#ifdef QMC_COMPLEX
TEST_CASE("MagDensity", "[hamiltonian]")
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!!   Evaluate MagDensity   !!!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  using RealType          = QMCTraits::RealType;
  using ValueType         = QMCTraits::ValueType;
  using PosType           = QMCTraits::PosType;
  using GradType          = QMCTraits::GradType;
  using LogValueType      = WaveFunctionComponent::LogValueType;
  using ParticlePos       = ParticleSet::ParticlePos;
  using ParticleGradient  = ParticleSet::ParticleGradient;
  using ParticleLaplacian = ParticleSet::ParticleLaplacian;
  using PropertySetType   = OperatorBase::PropertySetType;
  using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;
  using WP        = WalkerProperties::Indexes;
// O2 test example from pwscf non-collinear calculation.
  ParticleSet::ParticleLayout lattice;
  lattice.R(0, 0) = 5.10509515;
  lattice.R(0, 1) = -3.23993545;
  lattice.R(0, 2) = 0.00000000;
  lattice.R(1, 0) = 5.10509515;
  lattice.R(1, 1) = 3.23993545;
  lattice.R(1, 2) = 0.00000000;
  lattice.R(2, 0) = -6.49690625;
  lattice.R(2, 1) = 0.00000000;
  lattice.R(2, 2) = 7.08268015;
  
  lattice.BoxBConds = true; //periodic
  lattice.LR_dim_cutoff = 15; 
  lattice.reset();
  //Shamelessly stealing this from test_einset.cpp.  3 particles though.
  const SimulationCell simulation_cell(lattice);
  ParticleSet ions_(simulation_cell);
  ParticleSet elec_(simulation_cell);
  ions_.setName("ion");
  ions_.create({2});

  //ConstantSPOSet dnchannel(
  ions_.R[0][0] = 0.00000000;
  ions_.R[0][1] = 0.00000000;
  ions_.R[0][2] = 1.08659253;
  ions_.R[1][0] = 0.00000000;
  ions_.R[1][1] = 0.00000000;
  ions_.R[1][2] = -1.08659253;

  elec_.setName("elec");
  elec_.create({2});

  elec_.R[0][0] = 0;
  elec_.R[0][1] = 0;
  elec_.R[0][2] = 0;
  elec_.R[1][0] = 5;
  elec_.R[1][1] = 0;
  elec_.R[1][2] = 0;

  elec_.spins[0] = 1.9;
  elec_.spins[1] = 2.5410;
  elec_.setSpinor(true);

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  elec_.addTable(ions_);
  elec_.resetGroups();
  elec_.update();
  // </steal>


  const auto nelec = elec_.R.size();
  //Our test case is going to be three electron gas orbitals distinguished by 3 different kpoints.
  //Independent SPO's for the up and down channels.
  //
 /*
  std::vector<PosType> kup, kdn;
  std::vector<RealType> k2up, k2dn;


  kup.resize(nelec);
  kup[0] = PosType(0.25, 0.25, 0.25);
  kup[1] = PosType(0.1, 0.2, 0.3);

  kdn.resize(nelec);
  kdn[0] = PosType(0.25, 0.25, 0.25);
  kdn[1] = PosType(0.1, 0.2, 0.3);

  auto spo_up = std::make_unique<FreeOrbital>(kup);
  auto spo_dn = std::make_unique<FreeOrbital>(kdn);

  */
  using ValueVector = OrbitalSetTraits<ValueType>::ValueVector;
  ValueVector uprow0{ValueType(0.92387953,0),ValueType(0.92387953,0.)};
  ValueVector dnrow0{ValueType(0.27059805,0.27059805),ValueType(0.27059805,0.27059805)};
  ValueVector uprow1{ValueType(0.29131988,0.87674747),ValueType(0.81078057,0.44293144)};
  ValueVector dnrow1{ValueType(-0.17146777,0.342119),ValueType(0.10774051,0.36720375)};

  using ValueMatrix = OrbitalSetTraits<ValueType>::ValueMatrix;
  ValueMatrix mup,mdn;
  mup.resize(nelec,nelec);
  mdn.resize(nelec,nelec);

  for(int iorb=0; iorb<nelec; iorb++)
  {
    mup(0,iorb)=uprow0[iorb];
    mdn(0,iorb)=dnrow0[iorb];
    mup(1,iorb)=uprow1[iorb];
    mdn(1,iorb)=dnrow1[iorb];
  }
  auto spo_up = std::make_unique<ConstantSPOSet>(mup);
  auto spo_dn = std::make_unique<ConstantSPOSet>(mdn);

  auto spinor_set = std::make_unique<SpinorSet>();
  spinor_set->set_spos(std::move(spo_up), std::move(spo_dn));

  app_log() << " nelec=" << nelec << std::endl;

  auto dd = std::make_unique<DiracDeterminant<>>(std::move(spinor_set), 0, nelec);

  TrialWaveFunction psi;
  psi.addComponent(std::move(dd));
  app_log()<<"psi_val = "<<psi.evaluateLog(elec_)<<std::endl;
  MagDensityEstimator magdensity(elec_, psi);

  const char* magtxt = "<tmp> \
  <estimator name=\"MagDensity\" type=\"magdensity\" delta=\"0.5 0.5 0.5\" nsamples=\"9\"/> \
  </tmp> \
  ";

  Libxml2Document doc;
  bool okay = doc.parseFromString(magtxt);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr magxml = xmlFirstElementChild(root);
  magdensity.put(magxml);
  PropertySetType Observables;
  MCPWalker my_walker(nelec);
  my_walker.Weight=1.0;
  elec_.saveWalker(my_walker);
  magdensity.setHistories(my_walker);
  magdensity.addObservables(Observables,elec_.Collectables);
  magdensity.evaluate(elec_);
  app_log()<<"Collectables = \n";
  for (int i=0; i<elec_.Collectables.size(); i++) 
    app_log()<<i<<" "<<elec_.Collectables[i]<<std::endl;

}
#endif
}
