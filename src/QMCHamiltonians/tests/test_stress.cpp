//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//
// File created by: Yubo "Paul" Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "LongRange/EwaldHandler3D.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/StressPBC.h"
#include "Utilities/RuntimeOptions.h"
#include "OhmmsPETE/Tensor.h"

//This is for the Assaraf-Caffarel stress estimator
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "type_traits/ConvertToReal.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
// PBC case
TEST_CASE("Stress BCC H Ewald3D", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(3.24957306);
  lattice.LR_dim_cutoff = 40;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({2});
  ions.R[0]                     = {0.0, 0.0, 0.0};
  ions.R[1]                     = {1.62478653, 1.62478653, 1.62478653};
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  ion_species(pChargeIdx, pIdx) = 1;
  ions.resetGroups();
  ions.createSK();
  ions.update();

  elec.setName("elec");
  elec.create({2});
  elec.R[0]                  = {0.4, 0.4, 0.4};
  elec.R[1]                  = {2.02478653, 2.02478653, 2.02478653};
  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();
  elec.createSK();

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);

  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);
  LRCoulombSingleton::CoulombDerivHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombDerivHandler->initBreakup(ions);

  StressPBC est(ions, elec, psi);

  elec.update();
  est.evaluate(elec);

  // i-i = e-e stress is validated against Quantum Espresso's ewald method
  //  they are also double checked using finite-difference
  CHECK(est.getStressEE()(0, 0) == Approx(-0.01087883));
  CHECK(est.getStressEE()(0, 1) == Approx(0.0));
  CHECK(est.getStressEE()(0, 2) == Approx(0.0));
  CHECK(est.getStressEE()(1, 0) == Approx(0.0));
  CHECK(est.getStressEE()(1, 1) == Approx(-0.01087883));
  CHECK(est.getStressEE()(1, 2) == Approx(0.0));
  CHECK(est.getStressEE()(2, 0) == Approx(0.0));
  CHECK(est.getStressEE()(2, 1) == Approx(0.0));
  CHECK(est.getStressEE()(2, 2) == Approx(-0.01087883));

  //Electron-Ion stress diagonal is internally validated using fd.
  CHECK(est.getStressEI()(0, 0) == Approx(-0.00745376));
  //CHECK(est.getStressEI()(0, 1) == Approx(0.0));
  //CHECK(est.getStressEI()(0, 2) == Approx(0.0));
  //CHECK(est.getStressEI()(1, 0) == Approx(0.0));
  CHECK(est.getStressEI()(1, 1) == Approx(-0.00745376));
  //CHECK(est.getStressEI()(1, 2) == Approx(0.0));
  //CHECK(est.getStressEI()(2, 0) == Approx(0.0));
  //CHECK(est.getStressEI()(2, 1) == Approx(0.0));
  CHECK(est.getStressEI()(2, 2) == Approx(-0.00745376));

  LRCoulombSingleton::CoulombHandler.reset(nullptr);

  LRCoulombSingleton::CoulombDerivHandler.reset(nullptr);
}
TEST_CASE("Eloc_Derivatives:slater_fastderiv_stress_complex_pbc", "[hamiltonian]")
{
  using RealType = QMCTraits::RealType;
  using ValueMatrix = SPOSet::ValueMatrix;
  app_log() << "====Ion Derivative Test: Single Slater No Jastrow Stress Complex PBC====\n";
  std::ostringstream section_name;
  section_name << "Carbon diamond off gamma unit test: ";

  Communicate* c = OHMMS::Controller;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds[0]             = 1; // periodic
  lattice.BoxBConds[1]             = 1; // periodic
  lattice.BoxBConds[2]             = 1; // periodic
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
  lattice.LR_dim_cutoff            = 60;
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::EWALD;
// LRCoulombSingleton::this_lr_type = LRCoulombSingleton::ESLER;
  lattice.reset();
  const SimulationCell simulation_cell(lattice);
  auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto &ions(*ions_ptr), elec(*elec_ptr);

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
  int Nions = ions.getTotalNum();
  int Nelec = elec.getTotalNum();

  app_log()<<"before elec=\n"<<elec.R<<std::endl;
  app_log()<<"before ion=\n"<<ions.R<<std::endl;
  for (int e=0; e<Nelec; e++)
    elec.R[e] = dot(strans,elec.R[e]);

  for (int d=0; d<Nions; d++)
    ions.R[d] = dot(strans,ions.R[d]);

  app_log()<<"after elec=\n"<<elec.R<<std::endl;
  app_log()<<"after ion=\n"<<ions.R<<std::endl;

  ions.resetGroups();
  elec.resetGroups();
  ions.createSK();
  elec.createSK();
  elec.addTable(elec);
  elec.addTable(ions);
  elec.update();
  
//create_C_pbc_particlesets(elec, ions);
  
 //Oh my goodness this is important.  Otherwise all the numbers come out wrong :) (: :)
  LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombHandler->initBreakup(ions);
  LRCoulombSingleton::CoulombDerivHandler = std::make_unique<EwaldHandler3D>(ions);
  LRCoulombSingleton::CoulombDerivHandler->initBreakup(ions);

  HamiltonianFactory::PSetMap particle_set_map;
  particle_set_map.emplace("e", std::move(elec_ptr));
  particle_set_map.emplace("ion0", std::move(ions_ptr));

  WaveFunctionFactory wff(elec, particle_set_map, c);

  Libxml2Document wfdoc;
  bool wfokay = wfdoc.parse("C_diamond-twist-third-cart.wfj.xml");
  REQUIRE(wfokay);

  RuntimeOptions runtime_options;
  xmlNodePtr wfroot = wfdoc.getRoot();

  OhmmsXPathObject wfnode("//wavefunction[@name='psi0']", wfdoc.getXPathContext());
  HamiltonianFactory::PsiPoolType psi_map;
  psi_map.emplace("psi0", wff.buildTWF(wfnode[0], runtime_options));

  TrialWaveFunction* psi = psi_map["psi0"].get();
  REQUIRE(psi != nullptr);

  HamiltonianFactory hf("h0", elec, particle_set_map, psi_map, c);

  const char* hamiltonian_xml = "<hamiltonian name=\"h0\" pbc=\"yes\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
         <pairpot type=\"coulomb\" name=\"IonIon\" source=\"ion0\" target=\"ion0\"/> \
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
           <pseudo elementType=\"C\" href=\"C.BFD.xml\"/> \
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
  RealType logpsi = psi->evaluateLog(elec);
  RealType eloc = ham.evaluateDeterministic(elec);
  enum observ_id
  {
    KINETIC = 0,
    LOCALECP,
    NONLOCALECP,
    ELECELEC,
    IONION
  };


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

  psi->initializeTWFFastDerivWrapper(elec, twf);


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

  //The first index corresponds to the x,y,z force derivative.  Current interface assumes that the ion index is fixed,
  // so these vectors of vectors of matrices store the derivatives of the M and B matrices.
  // dB[0][0] is the x component of the iat force derivative of the up B matrix, dB[0][1] is for the down B matrix.

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

  twf.getM(elec, matlist);

  OperatorBase* kinop = ham.getHamiltonian(KINETIC);
  app_log() << kinop << std::endl;
  kinop->evaluateOneBodyOpMatrix(elec, twf, B);


  std::vector<ValueMatrix> minv;
  std::vector<ValueMatrix> B_gs, M_gs; //We are creating B and M matrices for assumed ground-state occupations.
                                       //These are N_s x N_s square matrices (N_s is number of particles for species s).
  B_gs.push_back(upmat);
  B_gs.push_back(dnmat);
  M_gs.push_back(upmat);
  M_gs.push_back(dnmat);
  minv.push_back(upmat);
  minv.push_back(dnmat);


  //  twf.getM(elec, matlist);
  std::vector<std::vector<ValueMatrix>> dB_gs;
  std::vector<std::vector<ValueMatrix>> dM_gs;
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
      for(int q=0; q<ptclnum; q++)
      {
        app_log()<<"i="<<p<<" j="<<q<<B[id][p][q]<<std::endl;
      }
  }


  dB_gs.push_back(tmp_gs);
  dB_gs.push_back(tmp_gs);
  dB_gs.push_back(tmp_gs);

  dM_gs.push_back(tmp_gs);
  dM_gs.push_back(tmp_gs);
  dM_gs.push_back(tmp_gs);

  //Finally, we have all the data structures with the right dimensions.  Continue.

  ValueType keval = 0.0;
  RealType keobs  = 0.0;
  keval           = twf.trAB(minv, B_gs);
  convertToReal(keval, keobs);


  CHECK(keobs == Approx(7.3599525590e+00));

  app_log() << " KEVal = " << keval << std::endl;

  app_log() << " Now evaluating nonlocalecp\n";
  OperatorBase* nlppop = ham.getHamiltonian(NONLOCALECP);
  app_log() << "  Evaluated.  Calling evaluteOneBodyOpMatrix\n";

  kinop->evaluateOneBodyOpMatrixStrainDeriv(elec,twf,0,0,B);

  twf.wipeMatrices(B);
  twf.wipeMatrices(B_gs);
  twf.wipeMatrices(X);
  nlppop->evaluateOneBodyOpMatrix(elec, twf, B);
  twf.getGSMatrices(B, B_gs);
  twf.buildX(minv, B_gs, X);

  ValueType nlpp    = 0.0;
  RealType nlpp_obs = 0.0;
  nlpp              = twf.trAB(minv, B_gs);
  convertToReal(nlpp, nlpp_obs);

  app_log() << "NLPP = " << nlpp << std::endl;

  CHECK(nlpp_obs == Approx(-2.4018757000e-02));

  //Reference values are from finite difference (step size 1e-5) from a different code path.

  LRCoulombSingleton::CoulombHandler.reset(nullptr);

  LRCoulombSingleton::CoulombDerivHandler.reset(nullptr);
} 
} // namespace qmcplusplus
