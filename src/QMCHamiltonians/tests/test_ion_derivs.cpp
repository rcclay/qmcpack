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
  Communicate* comm = OHMMS::Controller;

  using RealType = QMCTraits::RealType;
  using ValueMatrix = SPOSet::ValueMatrix;

  RuntimeOptions runtime_options;
  
  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(runtime_options, comm, particle_pool);

  ParticleSet* elec_ = particle_pool.getParticleSet("e"); 
  ParticleSet* ions_ = particle_pool.getParticleSet("ion");

  int Nelec = elec_->getTotalNum();
  int Nions = ions_->getTotalNum();

  elec_->R[0] = {3.6006741306e+00, 1.0104445324e+00, 3.9141099719e+00};
  elec_->R[1] = {2.6451694427e+00, 3.4448681473e+00, 5.8351296103e+00};
  elec_->R[2] = {2.5458446692e+00, 4.5219372791e+00, 4.4785209995e+00};
  elec_->R[3] = {2.8301650128e+00, 1.5351128324e+00, 1.5004137310e+00};
  elec_->R[4] = {5.6422291182e+00, 2.9968904592e+00, 3.3039907052e+00};
  elec_->R[5] = {2.6062992989e+00, 4.0493925313e-01, 2.5900053291e+00};
  elec_->R[6] = {8.1001577415e-01, 9.7303865512e-01, 1.3901383112e+00};
  elec_->R[7] = {1.6343332400e+00, 6.1895704609e-01, 1.2145253306e+00};
 
  SimulationCell simulation_cell =particle_pool.getSimulationCell();

  RealType lat_const               = 3.37316115000000e+00;
  using Tensor_t = Tensor<RealType, 3>;
  Tensor_t strans = {1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0};
  Tensor_t refR = {lat_const, lat_const, 0, 0, lat_const, lat_const, lat_const, 0, lat_const};
//  Tensor_t refR = {lat_const, 0, 0, 0, lat_const, 0, 0, 0, lat_const};

  RealType sdelta=-1e-7;
  int si=2,sj=2;
  strans(si,sj)+=sdelta;
  
  app_log()<<"Transformation matrix\n";
  app_log()<<strans<<std::endl; 
  app_log()<<"Before L\n";
  app_log()<<refR<<std::endl;
  lattice.R                        = dot(strans,refR);
  app_log()<<"After L\n";
  app_log()<<lattice.R<<std::endl;

  app_log()<<"before elec=\n"<<elec_->R<<std::endl;
  app_log()<<"before ion=\n"<<ions_->R<<std::endl;
  for (int e=0; e<Nelec; e++)
    elec_->R[e] = dot(strans,elec_->R[e]);

  for (int d=0; d<Nions; d++)
    ions_->R[d] = dot(strans,ions_->R[d]);

  app_log()<<"after elec=\n"<<elec_->R<<std::endl;
  app_log()<<"after ion=\n"<<ions_->R<<std::endl;


  ions_->addTable(*ions_);
  elec_->addTable(*ions_);
  elec_->addTable(*elec_);
  ions_->update(); 
  elec_->update(); 

  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::EWALD;
  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(*ions_);
  LRCoulombSingleton::CoulombHandler->initBreakup(*ions_);
  LRCoulombSingleton::CoulombDerivHandler = std::make_unique<EwaldHandler3D>(*ions_);
  LRCoulombSingleton::CoulombDerivHandler->initBreakup(*ions_);

  const char* hamiltonian_xml = "<hamiltonian name=\"h0\" pbc=\"yes\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
         <pairpot type=\"coulomb\" name=\"IonIon\" source=\"ion\" target=\"ion\"/> \
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
           <pseudo elementType=\"C\" href=\"C.ccECP.xml\"/> \
         </pairpot> \
         </hamiltonian>";

  
  TrialWaveFunction* psi = wavefunction_pool.getPrimary();
  HamiltonianPool hamiltonian_pool(particle_pool, wavefunction_pool, comm);

  Libxml2Document doch;
  doch.parseFromString(hamiltonian_xml);

  xmlNodePtr rooth = doch.getRoot();
  hamiltonian_pool.put(rooth);

  QMCHamiltonian* ham = hamiltonian_pool.getPrimary();
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  RealType logpsi = psi->evaluateLog(*elec_);
  app_log() << "SRESS LOGPSI = " << std::setprecision(16)<< logpsi << std::endl;
  RealType eloc = ham->evaluateDeterministic(*elec_);
  enum observ_id
  {
    KINETIC = 0,
    LOCALECP,
    NONLOCALECP,
    ELECELEC,
    IONION
  };

  app_log() << "STRESS LocalEnergy = " << std::setprecision(16) << eloc << std::endl;
  app_log() << "STRESS Kinetic = " << std::setprecision(16) << ham->getObservable(KINETIC) << std::endl;
  app_log() << "STRESS LocalECP = " << std::setprecision(16) << ham->getObservable(LOCALECP) << std::endl;
  app_log() << "STRESS NonLocalECP = " << std::setprecision(16) << ham->getObservable(NONLOCALECP) << std::endl;
  app_log() << "STRESS ELECELEC = " << std::setprecision(16) << ham->getObservable(ELECELEC) << std::endl;
  app_log() << "STRESS IonIon = " << std::setprecision(16) << ham->getObservable(IONION) << std::endl;

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

  psi->initializeTWFFastDerivWrapper(*elec_, twf);


  //This builds and initializes all the auxiliary matrices needed to do fast derivative evaluation.
  //These matrices are not necessarily square to accomodate orb opt and multidets.

  ValueMatrix upmat; //Up slater matrix.
  ValueMatrix dnmat; //Down slater matrix.
  int Nup  = 4;      //These are hard coded until the interface calls get implemented/cleaned up.
  int Ndn  = 4;
  int Norb = 26;
  upmat.resize(Nup, Norb);
  dnmat.resize(Ndn, Norb);
  
  //The first two lines consist of vectors of matrices.  The vector index corresponds to the species ID.
  //For example, matlist[0] will be the slater matrix for up electrons, matlist[1] will be for down electrons.
  std::vector<ValueMatrix> matlist; //Vector of slater matrices.
  std::vector<ValueMatrix> B, X;    //Vector of B matrix, and auxiliary X matrix.

  std::vector<std::vector<ValueMatrix>> dM; //Derivative of slater matrix.
  std::vector<std::vector<ValueMatrix>> dB; //Derivative of B matrices.
  matlist.push_back(upmat);
  matlist.push_back(dnmat);

  dM.push_back(matlist);
  dM.push_back(matlist);
  dM.push_back(matlist);

  dB.push_back(matlist);
  dB.push_back(matlist);
  dB.push_back(matlist);

  B.push_back(upmat);
  B.push_back(dnmat);

  X.push_back(upmat);
  X.push_back(dnmat);

  twf.getM(*elec_, matlist);

  OperatorBase* kinop = ham->getHamiltonian(KINETIC);
  app_log() << kinop << std::endl;
  kinop->evaluateOneBodyOpMatrix(*elec_, twf, B);


  std::vector<ValueMatrix> minv;
  std::vector<ValueMatrix> B_gs, M_gs; //We are creating B and M matrices for assumed ground-state occupations.
                                       //These are N_s x N_s square matrices (N_s is number of particles for species s).
  B_gs.push_back(upmat);
  B_gs.push_back(dnmat);
  M_gs.push_back(upmat);
  M_gs.push_back(dnmat);
  minv.push_back(upmat);
  minv.push_back(dnmat);

  //twf.getStrainGradM(elec,0,0,dM);
  //  twf.getM(elec, matlist);
  std::vector<ValueMatrix> tmp_gs;
  twf.getGSMatrices(B, B_gs);
  twf.getGSMatrices(matlist, M_gs);
  twf.invertMatrices(M_gs, minv);
  twf.buildX(minv, B_gs, X);
  for (int id = 0; id < matlist.size(); id++)
  {
    //    int ptclnum = twf.numParticles(id);
    int ptclnum = (id == 0 ? Nup : Ndn); //hard coded until twf interface comes online.
    ValueMatrix gs_m;
    gs_m.resize(ptclnum, ptclnum);
    tmp_gs.push_back(gs_m);
    app_log()<<"matid="<<id<<std::endl;
    for (int p=0; p<ptclnum; p++)
      for(int q=0; q<Norb; q++)
      {
        app_log()<<"i="<<p<<" j="<<q<<" B= "<<B[id][p][q]<<" M= "<<matlist[id][p][q]<<" dM="<<dM[id][p][q]<<std::endl;
      }
  }


  //Finally, we have all the data structures with the right dimensions.  Continue.

  ValueType keval = 0.0;
  RealType keobs  = 0.0;
  keval           = twf.trAB(minv, B_gs);
  convertToReal(keval, keobs);
  CHECK(keobs == Approx(8.7183990166));

  app_log() << " KEVal = " << keval << std::endl;

  app_log() << " Now evaluating nonlocalecp\n";
  OperatorBase* nlppop = ham->getHamiltonian(NONLOCALECP);
  app_log() << "  Evaluated.  Calling evaluteOneBodyOpMatrix\n";


  twf.wipeMatrices(B);
  twf.wipeMatrices(B_gs);
  twf.wipeMatrices(X);
  nlppop->evaluateOneBodyOpMatrix(*elec_, twf, B);
  twf.getGSMatrices(B, B_gs);
  twf.buildX(minv, B_gs, X);

  ValueType nlpp    = 0.0;
  RealType nlpp_obs = 0.0;
  nlpp              = twf.trAB(minv, B_gs);
  convertToReal(nlpp, nlpp_obs);

  app_log() << "NLPP = " << nlpp << std::endl;

  CHECK(nlpp_obs == Approx(0.0021584195));
}

} // namespace qmcplusplus
