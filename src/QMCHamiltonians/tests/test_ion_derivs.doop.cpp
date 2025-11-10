//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "type_traits/template_types.hpp"
#include "type_traits/ConvertToReal.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
#include "ParticleIO/XMLParticleIO.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/RuntimeOptions.h"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "LongRange/EwaldHandler3D.h"
#include "Utilities/ProjectData.h"

#include "Utilities/RuntimeOptions.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "BsplineFactory/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
namespace qmcplusplus
{
using DiracDet = DiracDeterminant<DelayedUpdate<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
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

void create_C_pbc_particlesets(ParticleSet& elec, ParticleSet& ions)
{
  ions.setName("ion0");
  ions.create({2});
  ions.R[0] = {0.1, 0.1, 0.1};
  ions.R[1] = {1.6865805750, 1.6865805750, 1.6865805750};
  ions.R[0][0] += -1e-5;
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("C");
  int pChargeIdx                = ion_species.addAttribute("charge");
  int iatnumber                 = ion_species.addAttribute("atomic_number");
  ion_species(pChargeIdx, pIdx) = 4;
  ion_species(iatnumber, pIdx)  = 6;

  elec.setName("e");
  elec.create({4, 4});
  elec.R[0] = {3.6006741306e+00, 1.0104445324e+00, 3.9141099719e+00};
  elec.R[1] = {2.6451694427e+00, 3.4448681473e+00, 5.8351296103e+00};
  elec.R[2] = {2.5458446692e+00, 4.5219372791e+00, 4.4785209995e+00};
  elec.R[3] = {2.8301650128e+00, 1.5351128324e+00, 1.5004137310e+00};
  elec.R[4] = {5.6422291182e+00, 2.9968904592e+00, 3.3039907052e+00};
  elec.R[5] = {2.6062992989e+00, 4.0493925313e-01, 2.5900053291e+00};
  elec.R[6] = {8.1001577415e-01, 9.7303865512e-01, 1.3901383112e+00};
  elec.R[7] = {1.6343332400e+00, 6.1895704609e-01, 1.2145253306e+00};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int dnIdx                  = tspecies.addSpecies("d");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(chargeIdx, dnIdx) = -1;
  tspecies(massIdx, dnIdx)   = 1.0;


  ions.resetGroups();
  elec.resetGroups();
  ions.createSK();
  elec.createSK();
  elec.addTable(elec);
  elec.addTable(ions);
  elec.update();
}

//Takes a HamiltonianFactory and handles the XML I/O to get a QMCHamiltonian pointer.  For CN molecule with pseudopotentials.
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

#include "OhmmsPETE/Tensor.h"
TEST_CASE("Eloc_Derivatives:slater_fastderiv_stress_complex_pbc_spline", "[hamiltonian]")
{
  Communicate* c = OHMMS::Controller;

  using RealType = QMCTraits::RealType;
  using ValueMatrix = SPOSet::ValueMatrix;

  app_log()<<"!!!!!!! Eloc_Derivatives:slater_fastderiv_stress_complex_pbc_spline !!!!!!!!!!!!\n";
  const DynamicCoordinateKind kind_selected = DynamicCoordinateKind::DC_POS;
  // diamondC_1x1x1
  Lattice lattice;
  lattice.R         = {3.37316115, 3.37316115, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};
  lattice.BoxBConds = {1, 1, 1};
  lattice.LR_dim_cutoff            = 40;
  lattice.LR_tol = 0.1;
  lattice.reset();

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell(), kind_selected);
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell(), kind_selected);
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({2});
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};
  SpeciesSet& ion_species       = ions_.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("C");
  int pChargeIdx                = ion_species.addAttribute("charge");
  int iatnumber                 = ion_species.addAttribute("atomic_number");
  ion_species(pChargeIdx, pIdx) = 4;
  ion_species(iatnumber, pIdx)  = 6;
 
  int Nions=ions_.getTotalNum();

  ions_.addTable(ions_);
  ions_.resetGroups();
  ions_.createSK(); 
  ions_.update();


  elec_.setName("e");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({4, 4});
  elec_.R[0] = {0.16, 0.51, 0.04};
  elec_.R[1] = {1.05, 1.1, 0.38};
  elec_.R[2] = {1.91, 1.33, 1.15};
  elec_.R[3] = {1.22, 1.73, 1.38};

  elec_.R[4] = {0.22,  0 , 0.51};
  elec_.R[5] = {0.53, 0.0, 0.25};
  elec_.R[6] = {1.71, 1.78, 1.72};
  elec_.R[7] = {1.01, 1.32, 1.25};

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;

  elec_.addTable(ions_);
  elec_.addTable(elec_);
  elec_.resetGroups();
  elec_.createSK(); 
  elec_.update();

  // make a ParticleSet Clone
  ParticleSet elec_clone(elec_);

  //diamondC_1x1x1
  const char* spo_xml = R"(
<sposet_collection type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" meshfactor="1.0" precision="float">
  <sposet name="updet" size="2"/>
</sposet_collection>
)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml);
  REQUIRE(okay);

  xmlNodePtr spo_root = doc.getRoot();
  xmlNodePtr ein1     = xmlFirstElementChild(spo_root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, spo_root);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo != nullptr);

  std::vector<std::unique_ptr<DiracDeterminantBase>> dets;
  dets.push_back(std::make_unique<DiracDet>(spo->makeClone(), 0, 2));
  dets.push_back(std::make_unique<DiracDet>(spo->makeClone(), 2, 4));

  auto slater_det = std::make_unique<SlaterDet>(elec_, std::move(dets));

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  psi.addComponent(std::move(slater_det));

  psi.evaluateLog(elec_);
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::EWALD;
  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions_);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions_);
  LRCoulombSingleton::CoulombDerivHandler = std::make_unique<EwaldHandler3D>(ions_);
  LRCoulombSingleton::CoulombDerivHandler->initBreakup(ions_);

 
  HamiltonianFactory::PSetMap particle_set_map;

  std::unique_ptr<ParticleSet> elec_ptr = std::make_unique<ParticleSet>(elec_);
  std::unique_ptr<ParticleSet> ions_ptr = std::make_unique<ParticleSet>(ions_);
  particle_set_map.emplace("e", std::move(elec_ptr));
  particle_set_map.emplace("ion0", std::move(ions_ptr));

  HamiltonianFactory::PsiPoolType psi_map;
  auto psi_ptr = psi.makeClone(elec_); 

  psi_map.emplace("psi0",std::move(psi_ptr));

  HamiltonianFactory hf("h0", elec_, particle_set_map, psi_map, c);

  const char* hamiltonian_xml = "<hamiltonian name=\"h0\" pbc=\"yes\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
         <pairpot type=\"coulomb\" name=\"IonIon\" source=\"ion0\" target=\"ion0\"/> \
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
           <pseudo elementType=\"C\" href=\"C.ccECP.xml\"/> \
         </pairpot> \
         </hamiltonian>";


  Libxml2Document hdoc;
  bool okay2 = hdoc.parseFromString(hamiltonian_xml);
  REQUIRE(okay2);

  xmlNodePtr hroot = hdoc.getRoot();
  hf.put(hroot);

  QMCHamiltonian& ham = *hf.getH();


  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  RealType logpsi = psi.evaluateLog(elec_);
  app_log() << "SRESS LOGPSI = " << std::setprecision(16)<< logpsi << std::endl;
  RealType eloc = ham.evaluateDeterministic(elec_);
  enum observ_id
  {
    KINETIC = 0,
    LOCALECP,
    NONLOCALECP,
    ELECELEC,
    IONION
  };

  app_log() << "STRESS LocalEnergy = " << std::setprecision(16) << eloc << std::endl;
  app_log() << "STRESS Kinetic = " << std::setprecision(16) << ham.getObservable(KINETIC) << std::endl;
  app_log() << "STRESS LocalECP = " << std::setprecision(16) << ham.getObservable(LOCALECP) << std::endl;
  app_log() << "STRESS NonLocalECP = " << std::setprecision(16) << ham.getObservable(NONLOCALECP) << std::endl;
  app_log() << "STRESS ELECELEC = " << std::setprecision(16) << ham.getObservable(ELECELEC) << std::endl;
  app_log() << "STRESS IonIon = " << std::setprecision(16) << ham.getObservable(IONION) << std::endl;

  //Now for the derivative tests
  ParticleSet::ParticleGradient wfgradraw;
  ParticleSet::ParticlePos hf_term;
  ParticleSet::ParticlePos pulay_term;
  ParticleSet::ParticlePos wf_grad;

  wfgradraw.resize(Nions);
  wf_grad.resize(Nions);
  hf_term.resize(Nions);
  pulay_term.resize(Nions);

  TWFFastDerivWrapper twf;

  psi.initializeTWFFastDerivWrapper(elec_, twf);


  //This builds and initializes all the auxiliary matrices needed to do fast derivative evaluation.
  //These matrices are not necessarily square to accomodate orb opt and multidets.

  ValueMatrix upmat; //Up slater matrix.
  ValueMatrix dnmat; //Down slater matrix.
  int Nup  = 4;      //These are hard coded until the interface calls get implemented/cleaned up.
  int Ndn  = 4;
  int Norb = 26;
  upmat.resize(Nup, Norb);
  dnmat.resize(Ndn, Norb);
  


}

} // namespace qmcplusplus
