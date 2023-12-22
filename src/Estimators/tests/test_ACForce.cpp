//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "Configuration.h"
#include "Estimators/ACForceInput.h"
#include "Estimators/ACForce.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/WaveFunctionTypes.hpp"
#include "ParticleIO/XMLParticleIO.h"

#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"

using std::string;


namespace qmcplusplus
{
namespace testing
{
class ACForceTests
{
  using WF       = WaveFunctionTypes<QMCTraits::ValueType, QMCTraits::FullPrecValueType>;
  using Position = QMCTraits::PosType;
  using Data     = ACForce::Data;
  using Real     = WF::Real;
  using Value    = WF::Value;

public:
  void testCopyConstructor(const ACForce& magdens)
  {
  }

  void testData(const ACForce& magdens, const Data& data)
  {
  }

};
void create_CN_particlesets(ParticleSet& elec, ParticleSet& ions)
{
  const char* particles = R"(<tmp>
  <particleset name="ion0" size="2">
    <group name="C">
      <parameter name="charge">4</parameter>
      <parameter name="valence">2</parameter>
      <parameter name="atomicnumber">6</parameter>
    </group>
    <group name="N">
      <parameter name="charge">5</parameter>
      <parameter name="valence">3</parameter>
      <parameter name="atomicnumber">7</parameter>
    </group>
    <attrib name="position" datatype="posArray">
  0.0000000000e+00  0.0000000000e+00  0.0000000000e+00
  0.0000000000e+00  0.0000000000e+00  2.0786985865e+00
</attrib>
    <attrib name="ionid" datatype="stringArray">
 C N
</attrib>
  </particleset>
  <particleset name="e">
    <group name="u" size="5">
      <parameter name="charge">-1</parameter>
    <attrib name="position" datatype="posArray">
-0.55936725 -0.26942464 0.14459603
0.19146719 1.40287983 0.63931251
1.14805915 -0.52057335 3.49621107
0.28293870 -0.10273952 0.01707021
0.60626935 -0.25538121 1.75750740
</attrib>
    </group>
    <group name="d" size="4">
      <parameter name="charge">-1</parameter>
    <attrib name="position" datatype="posArray">
-0.47405939 0.59523171 -0.59778601
0.03150661 -0.27343474 0.56279442
-1.32648025 0.00970226 2.26944242
2.42944286 0.64884151 1.87505288
</attrib>
    </group>
  </particleset>
  </tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root  = doc.getRoot();
  xmlNodePtr part1 = xmlFirstElementChild(root);
  xmlNodePtr part2 = xmlNextElementSibling(part1);

  XMLParticleParser parse_ions(ions);
  parse_ions.readXML(part1);

  XMLParticleParser parse_electrons(elec);
  parse_electrons.readXML(part2);
  ions.addTable(ions);
  elec.addTable(ions);
  elec.addTable(elec);
  elec.update();
  ions.update();
}
QMCHamiltonian& create_CN_Hamiltonian(HamiltonianFactory& hf)
{
  //Incantation to build hamiltonian
  const char* hamiltonian_xml = R"(<hamiltonian name="h0" type="generic" target="e">
         <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
         <pairpot type="coulomb" name="IonIon" source="ion0" target="ion0"/>
         <pairpot name="PseudoPot" type="pseudo" source="ion0" wavefunction="psi0" format="xml" algorithm="non-batched">
           <pseudo elementType="C" href="C.ccECP.xml"/>
           <pseudo elementType="N" href="N.ccECP.xml"/>
         </pairpot>
         </hamiltonian>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(hamiltonian_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  hf.put(root);

  return *hf.getH();
}
} //namespace testing

TEST_CASE("ACForce::IntegrationTest", "[estimators]")
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!!   Evaluate ACForce   !!!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

  //For now, do a small square case.
  const int nelec       = 2;
  const int norb        = 2;
  using WF              = WaveFunctionTypes<QMCTraits::ValueType, QMCTraits::FullPrecValueType>;
  using Real            = WF::Real;
  using Value           = WF::Value;
  using Grad            = WF::Grad;
  using ValueVector     = Vector<Value>;
  using GradVector      = Vector<Grad>;
  using ValueMatrix     = Matrix<Value>;
  using PropertySetType = OperatorBase::PropertySetType;
  using MCPWalker       = Walker<QMCTraits, PtclOnLatticeTraits>;
  using Data            = ACForce::Data;
  using GradMatrix      = Matrix<Grad>;
  using namespace testing;
  app_log() << "========================================================================================\n";
  app_log() << "========================================================================================\n";
  app_log() << "====================Ion Derivative Test:  Prototype Single Slater+Jastrow===========\n";
  app_log() << "========================================================================================\n";
  app_log() << "========================================================================================\n";
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;

  Communicate* c = OHMMS::Controller;

  const SimulationCell simulation_cell;
  auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto &ions(*ions_ptr), elec(*elec_ptr);

  //Build a CN test molecule.
  testing::create_CN_particlesets(elec, ions);
  //////////////////////////////////
  /////////////////////////////////
  //Incantation to read and build a TWF from cn.wfnoj//
  Libxml2Document doc2;
  bool okay = doc2.parse("cn.wfj.xml");
  REQUIRE(okay);
  xmlNodePtr root2 = doc2.getRoot();

  WaveFunctionComponentBuilder::PSetMap particle_set_map;
  particle_set_map.emplace("e", std::move(elec_ptr));
  particle_set_map.emplace("ion0", std::move(ions_ptr));

  RuntimeOptions runtime_options;
  WaveFunctionFactory wff(elec, particle_set_map, c);
  HamiltonianFactory::PsiPoolType psi_map;
  psi_map.emplace("psi0", wff.buildTWF(root2, runtime_options));

  TrialWaveFunction* psi = psi_map["psi0"].get();
  REQUIRE(psi != nullptr);
  //end incantation

  TWFFastDerivWrapper twf;

  psi->initializeTWFFastDerivWrapper(elec, twf);
  SPOSet::ValueVector values;
  SPOSet::GradVector dpsi;
  SPOSet::ValueVector d2psi;
  values.resize(9);
  dpsi.resize(9);
  d2psi.resize(9);

  HamiltonianFactory hf("h0", elec, particle_set_map, psi_map, c);

  QMCHamiltonian& ham = testing::create_CN_Hamiltonian(hf);

  //This is already defined in QMCHamiltonian, but keep it here for easy access.
  enum observ_id
  {
    KINETIC = 0,
    LOCALECP,
    NONLOCALECP,
    ELECELEC,
    IONION
  };

  using ValueMatrix = SPOSet::ValueMatrix;
  PropertySetType Observables;
  MCPWalker my_walker(nelec);
  my_walker.Weight = 1.0;

  elec.saveWalker(my_walker);

  //Now to create the crowd related quantities:
  int nwalkers = 2;
  std::vector<MCPWalker> walkers(nwalkers);
  std::vector<ParticleSet> psets{elec, elec};

 // psets[1].R[0][0] = 0;
 // psets[1].R[0][1] = 0.5;
 // psets[1].R[0][2] = 0;
 // psets[1].R[1][0] = 1;
 // psets[1].R[1][1] = 3;
 // psets[1].R[1][2] = 1;

  std::vector<UPtr<TrialWaveFunction>> twfcs(2);
  for (int iw = 0; iw < nwalkers; iw++)
    twfcs[iw] = psi->makeClone(psets[iw]);

  auto updateWalker = [](auto& walker, auto& pset_target, auto& trial_wavefunction) {
    pset_target.update(true);
    pset_target.donePbyP();
    trial_wavefunction.evaluateLog(pset_target);
    pset_target.saveWalker(walker);
  };
  for (int iw = 0; iw < nwalkers; iw++)
    updateWalker(walkers[iw], psets[iw], *(twfcs[iw]));

  std::vector<QMCHamiltonian> hams;
  auto ref_walkers(makeRefVector<MCPWalker>(walkers));
  auto ref_psets(makeRefVector<ParticleSet>(psets));
  auto ref_twfcs(convertUPtrToRefVector(twfcs));
  auto ref_hams(makeRefVector<QMCHamiltonian>(hams));

  FakeRandom rng;

  std::string_view input_xml = R"XML(<estimator name="doop" type="ACForce" epsilon="0.01" spacewarp="yes"/>)XML";

  Libxml2Document doc;
  xmlNodePtr node = doc.getRoot();
  ACForceInput acfin(node);
  ACForce acf(std::move(acfin),2);
 

  
}
} // namespace qmcplusplus
