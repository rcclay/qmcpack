//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by:
//
//
// File created by:
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "type_traits/template_types.hpp"
#include "type_traits/ConvertToReal.h"
#include "QMCHamiltonians/HamiltonianFactory.h"
#include <MinimalParticlePool.h>
#include <MinimalWaveFunctionPool.h>
#include <MinimalHamiltonianPool.h>
#include "ParticleIO/XMLParticleIO.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/RuntimeOptions.h"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDetTableMethod.h"
#include "QMCHamiltonians/CoulombPBCAB.h"
#include "LongRange/EwaldHandler3D.h"
#include "OhmmsPETE/SymTensor.h"

namespace qmcplusplus
{

void create_C_pbc_strained_particlesets(ParticleSet& elec, ParticleSet& ions)
{
  ions.setName("ion0");
  ions.create({2});
  ions.R[0] = {0.0, 0.0, 0.0};
  ions.R[1] = {1.68320741, 1.68995374, 1.68320741}; 
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("C");
  int pChargeIdx                = ion_species.addAttribute("charge");
  int iatnumber                 = ion_species.addAttribute("atomic_number");
  ion_species(pChargeIdx, pIdx) = 4;
  ion_species(iatnumber, pIdx)  = 6;

  elec.setName("e");
  elec.create({4, 4});

  elec.R[0] = {2.72492921, 2.52144024, 0.44105349};
  elec.R[1] = {3.38401004, 5.58872726, 3.38995410};
  elec.R[2] = {2.75148528, 1.55528232, 2.12699312};
  elec.R[3] = {0.67498584, 1.54704628, 1.11808223};
  elec.R[4] = {6.50921692, 5.93675975, 5.93981184};
  elec.R[5] = {2.87389015, 3.48275016, 1.03841721};
  elec.R[6] = {2.23066251, 0.66463443, 1.75609752};
  elec.R[7] = {1.90372709, 2.81801419, 1.87069193};

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

TEST_CASE("ZVZB stress test PBC no-Jastrow", "[hamiltonian]")
{

  using RealType = QMCTraits::RealType;
  using ValueMatrix = SPOSet::ValueMatrix;
  app_log() << "====Ion Derivative Test: Single Slater No Jastrow Stress Complex PBC====\n";
  std::ostringstream section_name;
  section_name << "Carbon diamond off gamma unit test: ";

  Communicate* c = OHMMS::Controller;

  Lattice lattice;
  lattice.BoxBConds[0]             = 1; // periodic
  lattice.BoxBConds[1]             = 1; // periodic
  lattice.BoxBConds[2]             = 1; // periodic
  lattice.R                        = {3.37653431,3.37653431,-0.00674632,-0.00674632,3.37653431,3.37653431,3.36304166,0.00674632,3.36304166};
  lattice.LR_dim_cutoff            = 30;
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::EWALD;
  lattice.reset();
  const SimulationCell simulation_cell(lattice);
  auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto &ions(*ions_ptr), elec(*elec_ptr);

  create_C_pbc_strained_particlesets(elec, ions);


  int Nions = ions.getTotalNum();
  int Nelec = elec.getTotalNum();

  HamiltonianFactory::PSetMap particle_set_map;
  particle_set_map.emplace("e", std::move(elec_ptr));
  particle_set_map.emplace("ion0", std::move(ions_ptr));

  WaveFunctionFactory wff(elec, particle_set_map, c);

  Libxml2Document wfdoc;
  bool wfokay = wfdoc.parse("qmc_strained.noj.wfj.xml");
  REQUIRE(wfokay);

  RuntimeOptions runtime_options;
  xmlNodePtr wfroot = wfdoc.getRoot();

  OhmmsXPathObject wfnode("//wavefunction[@name='psi0']", wfdoc.getXPathContext());
  auto psi_ptr = wff.buildTWF(wfnode[0], runtime_options);
  auto& psi(*psi_ptr);

  HamiltonianFactory hf("h0", elec, particle_set_map, psi, c);

  const char* hamiltonian_xml = "<hamiltonian name=\"h0\" pbc=\"yes\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
         <pairpot type=\"coulomb\" name=\"IonIon\" source=\"ion0\" target=\"ion0\"/> \
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
           <pseudo elementType=\"C\" href=\"C.ccECP.xml\"/> \
         </pairpot> \
         </hamiltonian>";


  Libxml2Document hdoc;
  REQUIRE(hdoc.parseFromString(hamiltonian_xml));

  xmlNodePtr hroot = hdoc.getRoot();
  hf.put(hroot);
  auto ham_ptr        = hf.releaseHamiltonian();
  QMCHamiltonian& ham = *ham_ptr;

  app_log()<<" R ="<<std::endl;
  app_log()<<elec.R<<std::endl;
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  RealType logpsi = psi.evaluateLog(elec);
  app_log() << "SRESS LOGPSI = " << std::setprecision(16)<< logpsi << std::endl;
  RealType eloc = ham.evaluateDeterministic(psi, elec);
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

  psi.initializeTWFFastDerivWrapper(elec, twf);


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

  std::vector<ValueMatrix> dM; //Derivative of slater matrix.
  std::vector<ValueMatrix> dB; //Derivative of B matrices.
  matlist.push_back(upmat);
  matlist.push_back(dnmat);

  dM.push_back(upmat);
  dM.push_back(dnmat);

  dB.push_back(upmat);
  dB.push_back(dnmat);

  B.push_back(upmat);
  B.push_back(dnmat);

  X.push_back(upmat);
  X.push_back(dnmat);

  twf.getM(elec, matlist);

  OperatorBase* kinop = ham.getComponent(KINETIC);
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
  std::vector<ValueMatrix> dB_gs;
  std::vector<ValueMatrix> dM_gs;
  std::vector<ValueMatrix> tmp_gs;
  twf.getGSMatrices(B, B_gs);
  twf.getGSMatrices(matlist, M_gs);
  twf.invertMatrices(M_gs, minv);
  twf.buildX(minv, B_gs, X);

  //app_log()<<" B_MATRIX_UP = "<<std::endl;
  //app_log()<<B_gs[0]<<std::endl;
  //app_log()<<" B_MATRIX_DN = "<<std::endl;
  //app_log()<<B_gs[1]<<std::endl;
  for (int id = 0; id < matlist.size(); id++)
  {
    //    int ptclnum = twf.numParticles(id);
    int ptclnum = (id == 0 ? Nup : Ndn); //hard coded until twf interface comes online.
    ValueMatrix gs_m;
    gs_m.resize(ptclnum, ptclnum);
    tmp_gs.push_back(gs_m);
  }


  dB_gs=tmp_gs;
  dB_gs=tmp_gs;

  dM_gs=tmp_gs;
  dM_gs=tmp_gs;

  //Finally, we have all the data structures with the right dimensions.  Continue.

  ValueType keval = 0.0;
  RealType keobs  = 0.0;
  keval           = twf.trAB(minv, B_gs);
  convertToReal(keval, keobs);
  CHECK(keobs == Approx(7.3599525590e+00));

  app_log() << " KEVal = " << keval << std::endl;
 // app_log()<<"M_gs_up=\n";
 // app_log()<<M_gs[0]<<std::endl;
 // app_log()<<"minv_up=\n";
 // app_log()<<minv[0]<<std::endl;
 // app_log()<<"M_gs_dn=\n";
 // app_log()<<M_gs[1]<<std::endl;
 // app_log()<<"minv_dn=\n";
 // app_log()<<minv[1]<<std::endl;

  ValueMatrix sref_kin,s_kin;
  ValueMatrix sref_logpsi,s_logpsi;

  s_kin.resize(3,3);
  s_logpsi.resize(3,3);
  sref_kin.resize(3,3);
  sref_logpsi.resize(3,3);

  sref_kin[0][0]=-6.746469250000331;
  sref_kin[0][1]=-2.839323449999931;
  sref_kin[0][2]=-0.6752923999999716;
  sref_kin[1][1]=-10.569138249999721;
  sref_kin[1][2]=0.2669740999996506;
  sref_kin[2][2]=-2.027908250000099;
  
  sref_logpsi[0][0]=-0.1336388499999508;
  sref_logpsi[0][1]=-0.20182555000003433;
  sref_logpsi[0][2]=1.6868885000000944;
  sref_logpsi[1][1]=-2.051093900000023;
  sref_logpsi[1][2]=-0.326234050000096;
  sref_logpsi[2][2]=-1.4712421500000517;

  for (int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      twf.wipeMatrices(dM);
      twf.wipeMatrices(dB);
      twf.wipeMatrices(dM_gs);
      twf.wipeMatrices(dB_gs);

      twf.getStrainGradM(elec, i, j, dM);
      twf.getGSMatrices(dM, dM_gs);
     // app_log()<<"ORBITAL INFO   STRAIN"<<" "<<i<<" "<<j<<"\n";
     // app_log()<<"  M_up="<<std::endl;
     // app_log()<<M_gs[0]<<std::endl;
     // app_log()<<"  dM_up="<<std::endl;
     // app_log()<<dM_gs[0]<<std::endl;
     // app_log()<<"  M_dn="<<std::endl;
     // app_log()<<M_gs[1]<<std::endl;
     // app_log()<<"  dM_dn="<<std::endl;
     // app_log()<<dM_gs[1]<<std::endl;
      //app_log()<<std::endl;

      kinop->evaluateOneBodyOpMatrixStrainDeriv(elec, twf, i,j, dB);

      twf.getGSMatrices(dB, dB_gs);
      //app_log()<<"  dBkin_up=\n";
      //app_log()<<dB_gs[0]<<std::endl;
      //app_log()<<"  dBkin_dn=\n";
      //app_log()<<dB_gs[1]<<std::endl;
      s_kin[i][j] = twf.computeGSDerivative(minv, X, dM_gs, dB_gs);
      s_logpsi[i][j] = twf.trAB(minv,dM_gs);
      app_log()<<"S="<<i<<" "<<j<<" "<<s_logpsi[i][j]<<" "<<sref_kin[i][j]<<" "<<s_kin[i][j]<<std::endl;
    }
  
  for (int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      CHECK(s_kin[i][j] == ComplexApprox(sref_kin[i][j]));
      CHECK(s_logpsi[i][j] == ComplexApprox(sref_logpsi[i][j]));
    }
 
   

  app_log() << " Now evaluating nonlocalecp\n";
  OperatorBase* nlppop = ham.getComponent(NONLOCALECP);
  app_log() << "  Evaluated.  Calling evaluteOneBodyOpMatrix\n";


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
  ValueMatrix s_nlpp;
  s_nlpp.resize(3,3);

  for (int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
  
      twf.wipeMatrices(dM);
      twf.wipeMatrices(dB);
      twf.wipeMatrices(dM_gs);
      twf.wipeMatrices(dB_gs);

      twf.getStrainGradM(elec, i, j, dM);
      twf.getGSMatrices(dM, dM_gs);
      app_log()<<"ORBITAL INFO   STRAIN"<<" "<<i<<" "<<j<<"\n";
      app_log()<<"  M_up="<<std::endl;
      app_log()<<M_gs[0]<<std::endl;
      app_log()<<"  dM_up="<<std::endl;
      app_log()<<dM_gs[0]<<std::endl;
      app_log()<<"  M_dn="<<std::endl;
      app_log()<<M_gs[1]<<std::endl;
      app_log()<<"  dM_dn="<<std::endl;
      app_log()<<dM_gs[1]<<std::endl;
      app_log()<<std::endl;

      nlppop->evaluateOneBodyOpMatrixStrainDeriv(elec, twf, i,j, dB);

      twf.getGSMatrices(dB, dB_gs);
      app_log()<<"  dBnlpp_up=\n";
      app_log()<<dB_gs[0]<<std::endl;
      app_log()<<"  dBnlpp_dn=\n";
      app_log()<<dB_gs[1]<<std::endl;
      s_nlpp[i][j] = twf.computeGSDerivative(minv, X, dM_gs, dB_gs);
      app_log()<<"S="<<i<<" "<<j<<" "<<s_nlpp[i][j]<<std::endl;
   }
  CHECK(nlpp_obs == Approx(-2.4018757000e-02));
}

TEST_CASE("ZVZB stress test PBC Jastrow", "[hamiltonian]")
{

  using RealType = QMCTraits::RealType;
  using ValueMatrix = SPOSet::ValueMatrix;
  app_log() << "====Ion Derivative Test: Single Slater No Jastrow Stress Complex PBC====\n";
  std::ostringstream section_name;
  section_name << "Carbon diamond off gamma unit test: ";

  Communicate* c = OHMMS::Controller;

  Lattice lattice;
  lattice.BoxBConds[0]             = 1; // periodic
  lattice.BoxBConds[1]             = 1; // periodic
  lattice.BoxBConds[2]             = 1; // periodic
  lattice.R                        = {3.37653431,3.37653431,-0.00674632,-0.00674632,3.37653431,3.37653431,3.36304166,0.00674632,3.36304166};
  lattice.LR_dim_cutoff            = 30;
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::EWALD;
  lattice.reset();
  const SimulationCell simulation_cell(lattice);
  auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto &ions(*ions_ptr), elec(*elec_ptr);

  create_C_pbc_strained_particlesets(elec, ions);


  int Nions = ions.getTotalNum();
  int Nelec = elec.getTotalNum();

  HamiltonianFactory::PSetMap particle_set_map;
  particle_set_map.emplace("e", std::move(elec_ptr));
  particle_set_map.emplace("ion0", std::move(ions_ptr));

  WaveFunctionFactory wff(elec, particle_set_map, c);

  Libxml2Document wfdoc;
  bool wfokay = wfdoc.parse("qmc_strained.wfj.xml");
  REQUIRE(wfokay);

  RuntimeOptions runtime_options;
  xmlNodePtr wfroot = wfdoc.getRoot();

  OhmmsXPathObject wfnode("//wavefunction[@name='psi0']", wfdoc.getXPathContext());
  auto psi_ptr = wff.buildTWF(wfnode[0], runtime_options);
  auto& psi(*psi_ptr);

  HamiltonianFactory hf("h0", elec, particle_set_map, psi, c);

  const char* hamiltonian_xml = "<hamiltonian name=\"h0\" pbc=\"yes\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
         <pairpot type=\"coulomb\" name=\"IonIon\" source=\"ion0\" target=\"ion0\"/> \
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"non-batched\"> \
           <pseudo elementType=\"C\" href=\"C.ccECP.xml\"/> \
         </pairpot> \
         </hamiltonian>";


  Libxml2Document hdoc;
  REQUIRE(hdoc.parseFromString(hamiltonian_xml));

  xmlNodePtr hroot = hdoc.getRoot();
  hf.put(hroot);
  auto ham_ptr        = hf.releaseHamiltonian();
  QMCHamiltonian& ham = *ham_ptr;

  app_log()<<" R ="<<std::endl;
  app_log()<<elec.R<<std::endl;
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  RealType logpsi = psi.evaluateLog(elec);
  app_log() << "SRESS LOGPSI = " << std::setprecision(16)<< logpsi << std::endl;
  RealType eloc = ham.evaluateDeterministic(psi, elec);
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

  psi.initializeTWFFastDerivWrapper(elec, twf);


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

  std::vector<ValueMatrix> dM; //Derivative of slater matrix.
  std::vector<ValueMatrix> dB; //Derivative of B matrices.
  matlist.push_back(upmat);
  matlist.push_back(dnmat);

  dM.push_back(upmat);
  dM.push_back(dnmat);

  dB.push_back(upmat);
  dB.push_back(dnmat);

  B.push_back(upmat);
  B.push_back(dnmat);

  X.push_back(upmat);
  X.push_back(dnmat);

  twf.getM(elec, matlist);

  OperatorBase* kinop = ham.getComponent(KINETIC);
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
  std::vector<ValueMatrix> dB_gs;
  std::vector<ValueMatrix> dM_gs;
  std::vector<ValueMatrix> tmp_gs;
  twf.getGSMatrices(B, B_gs);
  twf.getGSMatrices(matlist, M_gs);
  twf.invertMatrices(M_gs, minv);
  twf.buildX(minv, B_gs, X);

  app_log()<<" B_MATRIX_UP = "<<std::endl;
  app_log()<<B_gs[0]<<std::endl;
  app_log()<<" B_MATRIX_DN = "<<std::endl;
  app_log()<<B_gs[1]<<std::endl;
  for (int id = 0; id < matlist.size(); id++)
  {
    //    int ptclnum = twf.numParticles(id);
    int ptclnum = (id == 0 ? Nup : Ndn); //hard coded until twf interface comes online.
    ValueMatrix gs_m;
    gs_m.resize(ptclnum, ptclnum);
    tmp_gs.push_back(gs_m);
  }


  ParticleSet::ParticleGradient dG;
  ParticleSet::ParticleLaplacian dL;
  dG.resize(8);
  dL.resize(8);
  ValueType mystrain=0.0;

  RealType jval = twf.evaluateJastrowVGL(elec,dG, dL);

  dG=0;dL=0;
  app_log()<<"dL before = "<<dL<<std::endl;
  twf.getStrainGradJ(elec,0,0,mystrain,dG,dL);
	
  app_log()<<" J = "<<jval<<" dJ/dstrain = "<<mystrain<<std::endl;
  app_log()<<" dG/dstrain = "<<dG<<std::endl;
  app_log()<<" dL/dstrain = "<<dL<<std::endl;
  dB_gs=tmp_gs;
  dB_gs=tmp_gs;

  dM_gs=tmp_gs;
  dM_gs=tmp_gs;

  //Finally, we have all the data structures with the right dimensions.  Continue.

  ValueType keval = 0.0;
  RealType keobs  = 0.0;
  keval           = twf.trAB(minv, B_gs);
  convertToReal(keval, keobs);
  CHECK(keobs == Approx(7.3599525590e+00));

  app_log() << " KEVal = " << keval << std::endl;
  app_log()<<"M_gs_up=\n";
  app_log()<<M_gs[0]<<std::endl;
  app_log()<<"minv_up=\n";
  app_log()<<minv[0]<<std::endl;
  app_log()<<"M_gs_dn=\n";
  app_log()<<M_gs[1]<<std::endl;
  app_log()<<"minv_dn=\n";
  app_log()<<minv[1]<<std::endl;

  ValueMatrix sref_kin,s_kin;
  ValueMatrix sref_logpsi,s_logpsi;

  s_kin.resize(3,3);
  s_logpsi.resize(3,3);
  sref_kin.resize(3,3);
  sref_logpsi.resize(3,3);

  sref_kin[0][0]=-6.746469250000331;
  sref_kin[0][1]=-2.839323449999931;
  sref_kin[0][2]=-0.6752923999999716;
  sref_kin[1][1]=-10.569138249999721;
  sref_kin[1][2]=0.2669740999996506;
  sref_kin[2][2]=-2.027908250000099;
  
  sref_logpsi[0][0]=-0.1336388499999508;
  sref_logpsi[0][1]=-0.20182555000003433;
  sref_logpsi[0][2]=1.6868885000000944;
  sref_logpsi[1][1]=-2.051093900000023;
  sref_logpsi[1][2]=-0.326234050000096;
  sref_logpsi[2][2]=-1.4712421500000517;

  for (int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      twf.wipeMatrices(dM);
      twf.wipeMatrices(dB);
      twf.wipeMatrices(dM_gs);
      twf.wipeMatrices(dB_gs);

      twf.getStrainGradM(elec, i, j, dM);
      twf.getGSMatrices(dM, dM_gs);
      app_log()<<"ORBITAL INFO   STRAIN"<<" "<<i<<" "<<j<<"\n";
      app_log()<<"  M_up="<<std::endl;
      app_log()<<M_gs[0]<<std::endl;
      app_log()<<"  dM_up="<<std::endl;
      app_log()<<dM_gs[0]<<std::endl;
      app_log()<<"  M_dn="<<std::endl;
      app_log()<<M_gs[1]<<std::endl;
      app_log()<<"  dM_dn="<<std::endl;
      app_log()<<dM_gs[1]<<std::endl;
      app_log()<<std::endl;

      kinop->evaluateOneBodyOpMatrixStrainDeriv(elec, twf, i,j, dB);

      twf.getGSMatrices(dB, dB_gs);
      app_log()<<"  dBkin_up=\n";
      app_log()<<dB_gs[0]<<std::endl;
      app_log()<<"  dBkin_dn=\n";
      app_log()<<dB_gs[1]<<std::endl;
      s_kin[i][j] = twf.computeGSDerivative(minv, X, dM_gs, dB_gs);
      s_logpsi[i][j] = twf.trAB(minv,dM_gs);
      app_log()<<"S="<<i<<" "<<j<<" "<<s_logpsi[i][j]<<" "<<sref_kin[i][j]<<" "<<s_kin[i][j]<<std::endl;
    }
  
  for (int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      CHECK(s_kin[i][j] == ComplexApprox(sref_kin[i][j]));
      CHECK(s_logpsi[i][j] == ComplexApprox(sref_logpsi[i][j]));
    }
 
   

  app_log() << " Now evaluating nonlocalecp\n";
  OperatorBase* nlppop = ham.getComponent(NONLOCALECP);
  app_log() << "  Evaluated.  Calling evaluteOneBodyOpMatrix\n";


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

  ValueMatrix s_nlpp;
  s_nlpp.resize(3,3);

  for (int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
  
      twf.wipeMatrices(dM);
      twf.wipeMatrices(dB);
      twf.wipeMatrices(dM_gs);
      twf.wipeMatrices(dB_gs);

      twf.getStrainGradM(elec, i, j, dM);
      twf.getGSMatrices(dM, dM_gs);
      app_log()<<"ORBITAL INFO   STRAIN"<<" "<<i<<" "<<j<<"\n";
      app_log()<<"  M_up="<<std::endl;
      app_log()<<M_gs[0]<<std::endl;
      app_log()<<"  dM_up="<<std::endl;
      app_log()<<dM_gs[0]<<std::endl;
      app_log()<<"  M_dn="<<std::endl;
      app_log()<<M_gs[1]<<std::endl;
      app_log()<<"  dM_dn="<<std::endl;
      app_log()<<dM_gs[1]<<std::endl;
      app_log()<<std::endl;

      nlppop->evaluateOneBodyOpMatrixStrainDeriv(elec, twf, i,j, dB);

      twf.getGSMatrices(dB, dB_gs);
      app_log()<<"  dBnlpp_up=\n";
      app_log()<<dB_gs[0]<<std::endl;
      app_log()<<"  dBnlpp_dn=\n";
      app_log()<<dB_gs[1]<<std::endl;
      s_nlpp[i][j] = twf.computeGSDerivative(minv, X, dM_gs, dB_gs);
      app_log()<<"S="<<i<<" "<<j<<" "<<s_nlpp[i][j]<<std::endl;
   }
  CHECK(nlpp_obs == Approx(-2.4018757000e-02));
}

} // namespace qmcplusplus
