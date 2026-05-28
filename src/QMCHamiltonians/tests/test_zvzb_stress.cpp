//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "type_traits/template_types.hpp"
#include "type_traits/ConvertToReal.h"

#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDetTableMethod.h"
#include <MinimalWaveFunctionPool.h>
#include "LongRange/EwaldHandler3D.h"
#include "QMCHamiltonians/CoulombPBCAB.h"

#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"

#include "OhmmsData/Libxml2Doc.h"
#include "Utilities/RuntimeOptions.h"
#include "Utilities/RandomGenerator.h"

#include <memory>
#include <string>

namespace qmcplusplus
{

/** Build the strained carbon test particlesets used by the ZVZB stress tests. */
void create_C_pbc_strained_particlesets(ParticleSet& elec, ParticleSet& ions)
{
  ions.setName("ion0");
  ions.create({2});
  ions.R[0] = {0.0, 0.0, 0.0};
  ions.R[1] = {1.68320741, 1.68995374, 1.68320741};

  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("C");
  int pChargeIdx                = ion_species.addAttribute("charge");
  int atomic_number_idx         = ion_species.addAttribute("atomic_number");
  ion_species(pChargeIdx, pIdx) = 4;
  ion_species(atomic_number_idx, pIdx) = 6;

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

/** Small ownership-safe bundle for particles, wavefunction, and Hamiltonian. */
struct StressTestSystem
{
  using PSetMap = HamiltonianFactory::PSetMap;

  PSetMap particle_set_map;

  ParticleSet* ions_ptr = nullptr;
  ParticleSet* elec_ptr = nullptr;

  std::unique_ptr<TrialWaveFunction> psi_ptr;
  std::unique_ptr<QMCHamiltonian> ham_ptr;

  ParticleSet& ions() { return *ions_ptr; }
  ParticleSet& elec() { return *elec_ptr; }
  TrialWaveFunction& psi() { return *psi_ptr; }
  QMCHamiltonian& ham() { return *ham_ptr; }
};

/** Collection of the matrices needed for fast-derivative stress checks. */
struct StressMatrixBundle
{
  using ValueMatrix = SPOSet::ValueMatrix;

  std::vector<ValueMatrix> M;
  std::vector<ValueMatrix> M_gs;
  std::vector<ValueMatrix> Minv;

  std::vector<ValueMatrix> X_kin;
  std::vector<ValueMatrix> X_nlpp;

  std::vector<ValueMatrix> B_kin;
  std::vector<ValueMatrix> B_kin_gs;

  std::vector<ValueMatrix> B_nlpp;
  std::vector<ValueMatrix> B_nlpp_gs;

  std::vector<ValueMatrix> dM;
  std::vector<ValueMatrix> dM_gs;

  std::vector<ValueMatrix> dB_kin;
  std::vector<ValueMatrix> dB_kin_gs;

  std::vector<ValueMatrix> dB_nlpp;
  std::vector<ValueMatrix> dB_nlpp_gs;
};

/** Scalar contractions used to validate the matrix kernels. */
struct StressContractions
{
  QMCTraits::ValueType tr_Minv_Bkin;
  QMCTraits::ValueType tr_Minv_dM;
  QMCTraits::ValueType dBkin_minus_XdM;
  QMCTraits::ValueType tr_Minv_Bnlpp;
  QMCTraits::ValueType dBnlpp_minus_XdM;
};

/** Allocate per-group rectangular matrices [nptcl x norb]. */
inline std::vector<SPOSet::ValueMatrix> make_group_value_matrices(const TWFFastDerivWrapper& twf,
                                                                  const ParticleSet& P)
{
  std::vector<SPOSet::ValueMatrix> mats(twf.numGroups());
  for (int gid = 0; gid < twf.numGroups(); ++gid)
  {
    const int sid    = twf.getTWFGroupIndex(gid);
    const int first  = P.first(gid);
    const int last   = P.last(gid);
    const int nptcls = last - first;
    const int norbs  = twf.numOrbitals(sid);
    mats[sid].resize(nptcls, norbs);
    mats[sid] = 0.0;
  }
  return mats;
}

/** Allocate per-group square matrices [nptcl x nptcl]. */
inline std::vector<SPOSet::ValueMatrix> make_group_square_matrices(const TWFFastDerivWrapper& twf,
                                                                   const ParticleSet& P)
{
  std::vector<SPOSet::ValueMatrix> mats(twf.numGroups());
  for (int gid = 0; gid < twf.numGroups(); ++gid)
  {
    const int sid    = twf.getTWFGroupIndex(gid);
    const int first  = P.first(gid);
    const int last   = P.last(gid);
    const int nptcls = last - first;
    mats[sid].resize(nptcls, nptcls);
    mats[sid] = 0.0;
  }
  return mats;
}

/** Build a wavefunction + Hamiltonian test system from XML inputs. */
inline StressTestSystem build_stress_test_system(std::unique_ptr<ParticleSet>&& ions_uptr,
                                                 std::unique_ptr<ParticleSet>&& elec_uptr,
                                                 const std::string& wfn_xml_file,
                                                 const std::string& hamiltonian_xml,
                                                 Communicate* comm)
{
  StressTestSystem sys;

  // Move ownership into the particle pool map first
  sys.particle_set_map.emplace("e", std::move(elec_uptr));
  sys.particle_set_map.emplace("ion0", std::move(ions_uptr));

  // Recover stable pointers after transfer of ownership
  sys.elec_ptr = sys.particle_set_map.find("e")->second.get();
  sys.ions_ptr = sys.particle_set_map.find("ion0")->second.get();

  WaveFunctionFactory wff(sys.elec(), sys.particle_set_map, comm);

  Libxml2Document wfdoc;
  bool wfokay = wfdoc.parse(wfn_xml_file);
  REQUIRE(wfokay);

  RuntimeOptions runtime_options;
  OhmmsXPathObject wfnode("//wavefunction[@name='psi0']", wfdoc.getXPathContext());
  sys.psi_ptr = wff.buildTWF(wfnode[0], runtime_options);

  HamiltonianFactory hf("h0", sys.elec(), sys.particle_set_map, sys.psi(), comm);

  Libxml2Document hdoc;
  REQUIRE(hdoc.parseFromString(hamiltonian_xml));
  xmlNodePtr hroot = hdoc.getRoot();
  hf.put(hroot);
  sys.ham_ptr = hf.releaseHamiltonian();

  return sys;
}

/** Extract M, dM, B, dB, and GS slices for a chosen strain component. */
inline StressMatrixBundle extract_stress_matrices(ParticleSet& elec,
                                                  TrialWaveFunction& psi,
                                                  QMCHamiltonian& ham,
                                                  TWFFastDerivWrapper& twf,
                                                  int mu,
                                                  int nu,
                                                  int kinetic_index,
                                                  int nlpp_index)
{
  using ValueMatrix = TWFFastDerivWrapper::ValueMatrix;
  using GradMatrix  = TWFFastDerivWrapper::GradMatrix;

  StressMatrixBundle out;

  out.M         = make_group_value_matrices(twf, elec);
  out.M_gs      = make_group_square_matrices(twf, elec);
  out.Minv      = make_group_square_matrices(twf, elec);

  out.X_kin     = make_group_square_matrices(twf, elec);
  out.X_nlpp    = make_group_square_matrices(twf, elec);

  out.B_kin     = make_group_value_matrices(twf, elec);
  out.B_kin_gs  = make_group_square_matrices(twf, elec);

  out.B_nlpp    = make_group_value_matrices(twf, elec);
  out.B_nlpp_gs = make_group_square_matrices(twf, elec);

  out.dM        = make_group_value_matrices(twf, elec);
  out.dM_gs     = make_group_square_matrices(twf, elec);

  out.dB_kin    = make_group_value_matrices(twf, elec);
  out.dB_kin_gs = make_group_square_matrices(twf, elec);

  out.dB_nlpp    = make_group_value_matrices(twf, elec);
  out.dB_nlpp_gs = make_group_square_matrices(twf, elec);

  // M and inverse
  twf.getM(elec, out.M);
  twf.getGSMatrices(out.M, out.M_gs);
  twf.invertMatrices(out.M_gs, out.Minv);

  // B matrices
  ham.getComponent(kinetic_index)->evaluateOneBodyOpMatrix(elec, twf, out.B_kin);
  ham.getComponent(nlpp_index)->evaluateOneBodyOpMatrix(elec, twf, out.B_nlpp);

  twf.getGSMatrices(out.B_kin, out.B_kin_gs);
  twf.getGSMatrices(out.B_nlpp, out.B_nlpp_gs);

  // X matrices
  twf.buildX(out.Minv, out.B_kin_gs, out.X_kin);
  twf.buildX(out.Minv, out.B_nlpp_gs, out.X_nlpp);

  // dM
  std::vector<GradMatrix> dgrad_dummy(twf.numGroups());
  std::vector<ValueMatrix> dlapl_dummy = make_group_value_matrices(twf, elec);

  for (int gid = 0; gid < twf.numGroups(); ++gid)
  {
    const int sid    = twf.getTWFGroupIndex(gid);
    const int first  = elec.first(gid);
    const int last   = elec.last(gid);
    const int nptcls = last - first;
    const int norbs  = twf.numOrbitals(sid);

    dgrad_dummy[sid].resize(nptcls, norbs);
    dgrad_dummy[sid] = 0.0;
    dlapl_dummy[sid] = 0.0;
  }

  twf.getStrainGradM(elec, mu, nu, out.dM, dgrad_dummy, dlapl_dummy);
  twf.getGSMatrices(out.dM, out.dM_gs);

  // dB matrices
  ham.getComponent(kinetic_index)->evaluateOneBodyOpMatrixStrainDeriv(elec, twf, mu, nu, out.dB_kin);
  ham.getComponent(nlpp_index)->evaluateOneBodyOpMatrixStrainDeriv(elec, twf, mu, nu, out.dB_nlpp);

  twf.getGSMatrices(out.dB_kin, out.dB_kin_gs);
  twf.getGSMatrices(out.dB_nlpp, out.dB_nlpp_gs);

  return out;
}

/** Compute the scalar contractions used in the fast derivative estimator. */
inline StressContractions compute_stress_contractions(TWFFastDerivWrapper& twf,
                                                      const StressMatrixBundle& mats)
{
  StressContractions out;
  out.tr_Minv_Bkin      = twf.trAB(mats.Minv, mats.B_kin_gs);
  out.tr_Minv_dM        = twf.trAB(mats.Minv, mats.dM_gs);
  out.dBkin_minus_XdM   = twf.computeGSDerivative(mats.Minv, mats.X_kin, mats.dM_gs, mats.dB_kin_gs);
  out.tr_Minv_Bnlpp     = twf.trAB(mats.Minv, mats.B_nlpp_gs);
  out.dBnlpp_minus_XdM  = twf.computeGSDerivative(mats.Minv, mats.X_nlpp, mats.dM_gs, mats.dB_nlpp_gs);
  return out;
}
inline void print_group_value_matrices(const std::string& label,
                                       const std::vector<SPOSet::ValueMatrix>& mats)
{
  app_log() << label << "\n";
  for (size_t g = 0; g < mats.size(); ++g)
  {
    app_log() << "group " << g << ":\n";
    app_log() << mats[g] << "\n";
  }
  app_log() << std::endl;
}
inline void check_matrix_element(const std::vector<SPOSet::ValueMatrix>& mats,
                                 int group,
                                 int row,
                                 int col,
                                 const QMCTraits::ValueType& ref,
                                 double tol = 1e-7)
{
  REQUIRE(group < static_cast<int>(mats.size()));
  REQUIRE(row < mats[group].rows());
  REQUIRE(col < mats[group].cols());
  CHECK(mats[group](row, col) == ComplexApprox(ref).epsilon(tol));
}

/** Print contractions in a copy/paste-friendly form while gathering reference values. */
inline void print_stress_contractions(const std::string& label,
                                      int mu,
                                      int nu,
                                      const StressContractions& c)
{
  app_log() << label << "  mu=" << mu << "  nu=" << nu << "\n";
  app_log() << "  tr(Minv B_kin)         = " << c.tr_Minv_Bkin << "\n";
  app_log() << "  tr(Minv dM)            = " << c.tr_Minv_dM << "\n";
  app_log() << "  dBkin - XdM            = " << c.dBkin_minus_XdM << "\n";
  app_log() << "  tr(Minv B_nlpp)        = " << c.tr_Minv_Bnlpp << "\n";
  app_log() << "  dBnlpp - XdM           = " << c.dBnlpp_minus_XdM << "\n";
  app_log() << std::endl;
}

//Ray:  So there's a lot of data to check.  Every strain has a dB_kin and dB_nlpp matrix for
//each spin channel and for jastrows & no-jastrow systems.  I didn't feel like checking them all, 
//so here's the testing philosophy, which appeals to conspiracy.  
//
//1.) The Tr(Minv*dB - XdM) quantities should equal the finite difference strain derivatives of 
//    whatever the matrix refers to.  I checked local kinetic and NLPP derivatives with nexus
//    and confirmed that all the diagonals matched.  
//2.) I verified that all off-diagonal strain derivatives matched (actually, with the NLPP, 
//    0.5*(dO/de_i_j+dO/de_j_i) matches to finite-difference precision.  
//3.) I checked e_0_1, e_1_0 to provide a concrete finite difference check within this test.  
//4.) On the basis of this, I assume that the current state of the code is correct, 
//    and so the remainder of the checks are to catch changes from state of the current implementation.
//    I spot check the e_1_2 element and intermittent
//    entries from all the dB matrices at e_0_0, e_1_1, e_2_2, e_0_1, e_1_0, e_1_2.  If 
//    there's a breaking change that is not caught by this unit test, then they deserve a medal for
//    breaking the code in such a heroically subtle way.  
//
TEST_CASE("ZVZB stress matrix", "[hamiltonian]")
{
  Communicate* c = OHMMS::Controller;

  Lattice lattice;
  lattice.BoxBConds[0] = 1;
  lattice.BoxBConds[1] = 1;
  lattice.BoxBConds[2] = 1;
  lattice.R = {3.37653431, 3.37653431, -0.00674632,
               -0.00674632, 3.37653431, 3.37653431,
               3.36304166, 0.00674632, 3.36304166};
  lattice.LR_dim_cutoff            = 30;
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::EWALD;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);

  const char* hamiltonian_xml =
      "<hamiltonian name=\"h0\" pbc=\"yes\" type=\"generic\" target=\"e\"> \
         <pairpot type=\"coulomb\" name=\"ElecElec\" source=\"e\" target=\"e\"/> \
         <pairpot type=\"coulomb\" name=\"IonIon\" source=\"ion0\" target=\"ion0\"/> \
         <pairpot name=\"PseudoPot\" type=\"pseudo\" source=\"ion0\" wavefunction=\"psi0\" format=\"xml\" algorithm=\"batched\"> \
           <pseudo elementType=\"C\" href=\"C.ccECP.xml\"/> \
         </pairpot> \
       </hamiltonian>";

  constexpr int KINETIC     = 0;
  constexpr int NONLOCALECP = 2;

  SECTION("no jastrow")
  {
    auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
    create_C_pbc_strained_particlesets(*elec_ptr, *ions_ptr);

    auto sys = build_stress_test_system(std::move(ions_ptr),
                                        std::move(elec_ptr),
                                        "qmc_strained.noj.wfj.xml",
                                        hamiltonian_xml,
                                        c);

    TWFFastDerivWrapper twf;
    sys.psi().initializeTWFFastDerivWrapper(sys.elec(), twf);

    SECTION("mu=0 nu=0")
    {
      const int mu = 0;
      const int nu = 0;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("no-jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(7.858499059));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-6.321344104e-07));

      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);

      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(-8.78582608)); //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(-5.938670978e-07));

      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001490944851));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-7.949896997e-10));

      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(-0.01264947503)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(1.369739239e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-0.6825570039,  0.9940778847));
      check_matrix_element(mats.dB_kin, 0, 0, 1, QMCTraits::ValueType(-0.02674419388, -0.4466225522));
      check_matrix_element(mats.dB_kin, 0, 3, 3, QMCTraits::ValueType(-1.041067764,  -2.403566657));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.4551743677,  0.6629172209));
      check_matrix_element(mats.dB_kin, 1, 1, 3, QMCTraits::ValueType( 1.302426378,   3.00697872));
      check_matrix_element(mats.dB_kin, 1, 3, 1, QMCTraits::ValueType( 0.06385873173, 1.066432947));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.006776717005,  0.009869628217));
      check_matrix_element(mats.dB_nlpp, 0, 0, 1, QMCTraits::ValueType( 0.0003691487078, 0.006164753842));
      check_matrix_element(mats.dB_nlpp, 0, 2, 1, QMCTraits::ValueType(-0.001013086731, -0.01691845143));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.0006181845969,  0.0009003256482));
      check_matrix_element(mats.dB_nlpp, 1, 1, 1, QMCTraits::ValueType( 0.0002146511614, 0.003584657354));
      check_matrix_element(mats.dB_nlpp, 1, 2, 1, QMCTraits::ValueType(-0.0002569342475,-0.004290777602));
    }

    SECTION("mu=0 nu=1")
    {
      const int mu = 0;
      const int nu = 1;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("no-jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(7.858499059));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-6.321344104e-07));

      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);

      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(0.8103845785)); //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(2.07494008e-07));

      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001490944851));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-7.949896997e-10));

      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.01999140106)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(1.437201169e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-0.01569453196,  0.02285745942));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType(-0.7485492757,  -1.816014684));
      check_matrix_element(mats.dB_kin, 0, 3, 2, QMCTraits::ValueType(-1.014580338,   -2.461418074));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.1810250001,   0.2636452665));
      check_matrix_element(mats.dB_kin, 1, 1, 2, QMCTraits::ValueType( 0.2672963804,   0.6484731646));
      check_matrix_element(mats.dB_kin, 1, 3, 3, QMCTraits::ValueType( 0.7079872393,   1.634566521));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.009098676496,  0.01325133606));
      check_matrix_element(mats.dB_nlpp, 0, 0, 2, QMCTraits::ValueType(-0.001760505665, -0.004271066924));
      check_matrix_element(mats.dB_nlpp, 0, 3, 0, QMCTraits::ValueType(-0.002268208503,  0.003303424697));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.002993404912,  0.004359602709));
      check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType(-0.001773512318, -0.004302621674));
      check_matrix_element(mats.dB_nlpp, 1, 2, 0, QMCTraits::ValueType( 0.009329166367, -0.01358702184));
    }

    SECTION("mu=1 nu=0")
    {
      const int mu = 1;
      const int nu = 0;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("no-jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(7.858499059));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-6.321344104e-07));
      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(0.8103845785)); //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(2.07494008e-07));
      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001490944851));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-7.949896997e-10));
      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.001824137201)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(1.483571418e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-0.01569453196,  0.02285745942));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType(-0.7485492757,  -1.816014684));
      check_matrix_element(mats.dB_kin, 0, 3, 2, QMCTraits::ValueType(-1.014580338,   -2.461418074));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.1810250001,   0.2636452665));
      check_matrix_element(mats.dB_kin, 1, 1, 2, QMCTraits::ValueType( 0.2672963804,   0.6484731646));
      check_matrix_element(mats.dB_kin, 1, 3, 3, QMCTraits::ValueType( 0.7079872393,   1.634566521));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.009042032798,  0.01316884));
      check_matrix_element(mats.dB_nlpp, 0, 0, 2, QMCTraits::ValueType( 0.00188846374,    0.004581499089));
      check_matrix_element(mats.dB_nlpp, 0, 3, 0, QMCTraits::ValueType(-0.002289004268,  0.003333711683));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.002877869225,  0.004191336173));
      check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType( 0.001875456511,  0.00454994294));
      check_matrix_element(mats.dB_nlpp, 1, 2, 0, QMCTraits::ValueType( 0.009308370602, -0.01355673485));
    }

    SECTION("mu=1 nu=1")
    {
      const int mu = 1;
      const int nu = 1;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("no-jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(7.858499059));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-6.321344104e-07));
      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(-3.148339991)); //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(4.373013349e-07));
      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001490944851));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-7.949896997e-10));
      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.01314733922)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(8.027728474e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-0.6720004152,  0.9787031405));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType( 1.20640289,    2.926788426));
      check_matrix_element(mats.dB_kin, 0, 3, 1, QMCTraits::ValueType(-0.09344256437, -1.560477419));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.4217772863,  0.6142774822));
      check_matrix_element(mats.dB_kin, 1, 1, 3, QMCTraits::ValueType( 0.7723622621,  1.783192419));
      check_matrix_element(mats.dB_kin, 1, 3, 2, QMCTraits::ValueType( 0.4369013228,  1.059942496));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.01181444185,  0.01720658375));
      check_matrix_element(mats.dB_nlpp, 0, 1, 1, QMCTraits::ValueType(-0.001205988092,-0.02013990378));
      check_matrix_element(mats.dB_nlpp, 0, 3, 0, QMCTraits::ValueType(-0.0001288676075, 0.0001876831184));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.01093048106,  0.01591918095));
      //check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType(-0.00004590211803,-0.0007665468298));
      check_matrix_element(mats.dB_nlpp, 1, 3, 1, QMCTraits::ValueType(-0.001131254978,-0.01889184998));
    }

    SECTION("mu=2 nu=2")
    {
      const int mu = 2;
      const int nu = 2;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("no-jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(7.858499059));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-6.321344104e-07));
      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(-3.782832047)); //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(1.420834583e-06));
      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001490944851));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-7.949896997e-10));
      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.03088931028)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(7.338550546e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-0.260132384,  0.3788574735));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType( 0.7797416834, 1.891688907));
      check_matrix_element(mats.dB_kin, 0, 3, 3, QMCTraits::ValueType( 0.5585887391, 1.289642519));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.4499063373,  0.6552446291));
      check_matrix_element(mats.dB_kin, 1, 1, 2, QMCTraits::ValueType( 0.7835566901,  1.900944286));
      check_matrix_element(mats.dB_kin, 1, 3, 2, QMCTraits::ValueType( 0.4198713272,  1.018626959));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.003090290622,  0.004500707285));
      check_matrix_element(mats.dB_nlpp, 0, 1, 1, QMCTraits::ValueType( 0.00001267080187, 0.000211600992));
      check_matrix_element(mats.dB_nlpp, 0, 3, 0, QMCTraits::ValueType(-0.005044154367,  0.007346319362));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.01013278366,  0.01475741241));
      check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType( 0.00005500297419, 0.0001334397645));
      check_matrix_element(mats.dB_nlpp, 1, 3, 0, QMCTraits::ValueType(-0.0003713928989, 0.0005408975638));
    }

    SECTION("mu=1 nu=2")
    {
      const int mu = 1;
      const int nu = 2;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("no-jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(7.858499059));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-6.321344104e-07));
      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(0.4687736725));
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(7.056942694e-07));
      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001490944851));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-7.949896997e-10));
      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.0008575785036));
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(2.415482123e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-0.4851849574,  0.7066247985));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType(-0.2616531644, -0.6347825416));
      check_matrix_element(mats.dB_kin, 0, 3, 2, QMCTraits::ValueType(-0.6271180859, -1.52141705));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType( 0.1272257768,-0.1852921516));
      check_matrix_element(mats.dB_kin, 1, 1, 2, QMCTraits::ValueType(-0.8793126001,-2.133252535));
      check_matrix_element(mats.dB_kin, 1, 3, 3, QMCTraits::ValueType(-0.7396943358,-1.707770324));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType( 0.006952532632,-0.01012568647));
      check_matrix_element(mats.dB_nlpp, 0, 1, 1, QMCTraits::ValueType(-0.0001145367581,-0.001912752384));
      check_matrix_element(mats.dB_nlpp, 0, 3, 0, QMCTraits::ValueType(-0.001252325673, 0.00182389029));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.0105828161,  0.01541284069));
      check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType( 0.0006564847109,0.00159266181));
      check_matrix_element(mats.dB_nlpp, 1, 3, 1, QMCTraits::ValueType(-0.0001919644248,-0.003205787505));
    }
  }

  SECTION("jastrow")
  {
    auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
    create_C_pbc_strained_particlesets(*elec_ptr, *ions_ptr);

    auto sys = build_stress_test_system(std::move(ions_ptr),
                                        std::move(elec_ptr),
                                        "qmc_strained.wfj.xml",
                                        hamiltonian_xml,
                                        c);

    TWFFastDerivWrapper twf;
    sys.psi().initializeTWFFastDerivWrapper(sys.elec(), twf);

    SECTION("mu=0 nu=0")
    {
      const int mu = 0;
      const int nu = 0;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(8.719162572));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-5.564999444e-07));

      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);

      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(-7.450585969));  //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(-3.719495888e-07));

      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001809233814));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-8.42520724e-10));

      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(-0.01024410877)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(1.451906273e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-1.1131516,    1.621196991));
      check_matrix_element(mats.dB_kin, 0, 0, 1, QMCTraits::ValueType(-0.03536047762,-0.5905136148));
      check_matrix_element(mats.dB_kin, 0, 3, 3, QMCTraits::ValueType(-1.187512102,  -2.741670228));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.4168586782,  0.6071141512));
      check_matrix_element(mats.dB_kin, 1, 1, 3, QMCTraits::ValueType( 0.9450524569,  2.181891191));
      check_matrix_element(mats.dB_kin, 1, 3, 2, QMCTraits::ValueType( 0.5909010987,  1.433552921));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.007178744636,  0.01045514229));
      check_matrix_element(mats.dB_nlpp, 0, 0, 1, QMCTraits::ValueType( 0.0004379270061, 0.00731334549));
      check_matrix_element(mats.dB_nlpp, 0, 3, 3, QMCTraits::ValueType(-0.0006604937337,-0.001524915822));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.0005989787039, 0.0008723541389));
      check_matrix_element(mats.dB_nlpp, 1, 1, 1, QMCTraits::ValueType( 0.0002396718897, 0.004002501136));
      check_matrix_element(mats.dB_nlpp, 1, 2, 1, QMCTraits::ValueType(-0.0002372876735,-0.003962681659));
    }

    SECTION("mu=0 nu=1")
    {
      const int mu = 0;
      const int nu = 1;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(8.719162572));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-5.564999444e-07));

      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);

      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(0.3736770645)); //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(2.973051695e-07));

      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001809233814));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-8.42520724e-10));

      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.02208895135)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(1.645831167e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-0.7809267515,  1.137343666));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType(-0.3498672528, -0.8487939119));
      check_matrix_element(mats.dB_kin, 0, 3, 2, QMCTraits::ValueType(-1.114544483,  -2.703935646));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.5857263646,  0.85305339));
      check_matrix_element(mats.dB_kin, 1, 1, 3, QMCTraits::ValueType(-0.5987918831, -1.382461574));
      check_matrix_element(mats.dB_kin, 1, 3, 3, QMCTraits::ValueType( 0.6348175762,  1.465635959));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.009536936875,  0.01388961959));
      check_matrix_element(mats.dB_nlpp, 0, 0, 2, QMCTraits::ValueType(-0.002197484248, -0.005331196895));
      check_matrix_element(mats.dB_nlpp, 0, 3, 3, QMCTraits::ValueType(-0.0005109609361,-0.001179681768));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.002790984444,  0.004064796985));
      check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType(-0.001823902354, -0.004424870211));
      check_matrix_element(mats.dB_nlpp, 1, 3, 3, QMCTraits::ValueType(-0.0006946489118,-0.001603771634));
    }

    SECTION("mu=1 nu=0")
    {
      const int mu = 1;
      const int nu = 0;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(8.719162572));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-5.564999444e-07));
      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(0.3736770645)); //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(2.973051695e-07));
      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001809233814));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-8.42520724e-10));
      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.001691248473)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(1.331060512e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-0.7809267515,  1.137343666));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType(-0.3498672528, -0.8487939119));
      check_matrix_element(mats.dB_kin, 0, 3, 2, QMCTraits::ValueType(-1.114544483,  -2.703935646));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.5857263646,  0.85305339));
      check_matrix_element(mats.dB_kin, 1, 1, 3, QMCTraits::ValueType(-0.5987918831, -1.382461574));
      check_matrix_element(mats.dB_kin, 1, 3, 3, QMCTraits::ValueType( 0.6348175762,  1.465635959));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.009541067677,  0.01389563569));
      check_matrix_element(mats.dB_nlpp, 0, 0, 2, QMCTraits::ValueType( 0.001881721245,  0.004565141483));
      check_matrix_element(mats.dB_nlpp, 0, 3, 0, QMCTraits::ValueType(-0.001840217081,  0.002680096878));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.002843020935,  0.004140583059));
      check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType( 0.001935059921,  0.004694543544));
      check_matrix_element(mats.dB_nlpp, 1, 2, 0, QMCTraits::ValueType( 0.008243098272, -0.01200526951));
    }

    SECTION("mu=1 nu=1")
    {
      const int mu = 1;
      const int nu = 1;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(8.719162572));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-5.564999444e-07));
      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(-2.393652996)); //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(4.513888246e-07));
      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001809233814));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-8.42520724e-10));
      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.01019220967)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(7.881749281e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType(-1.00694976,  1.466524247));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType( 1.422618585, 3.451337571));
      check_matrix_element(mats.dB_kin, 0, 3, 1, QMCTraits::ValueType( 0.02927694199, 0.4889274797));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.2356651279,  0.3432232464));
      check_matrix_element(mats.dB_kin, 1, 1, 2, QMCTraits::ValueType( 0.154309589,  0.37436214));
      check_matrix_element(mats.dB_kin, 1, 3, 2, QMCTraits::ValueType( 0.8117272976, 1.969287385));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.01251318688, 0.01822423781));
      check_matrix_element(mats.dB_nlpp, 0, 1, 1, QMCTraits::ValueType(-0.001402115332, -0.02341521143));
      check_matrix_element(mats.dB_nlpp, 0, 3, 0, QMCTraits::ValueType(-0.00005833201142, 0.00008495489067));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.01074102498, 0.01564325662));
      check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType(-0.00004533200108, -0.0001099775102));
      check_matrix_element(mats.dB_nlpp, 1, 3, 1, QMCTraits::ValueType(-0.00100278741, -0.01674645409));
    }

    SECTION("mu=2 nu=2")
    {
      const int mu = 2;
      const int nu = 2;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(8.719162572));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-5.564999444e-07));
      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(-2.723835118)); //validated against finite difference
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(1.309284634e-06));
      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001809233814));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-8.42520724e-10));
      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.03844405173)); //validated against finite difference
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(8.460240796e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType( 0.2273138207,-0.3310604434));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType( 0.5738918478, 1.392287821));
      check_matrix_element(mats.dB_kin, 0, 3, 3, QMCTraits::ValueType( 0.2621335725, 0.6052012386));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.5611580446, 0.8172719429));
      check_matrix_element(mats.dB_kin, 1, 1, 2, QMCTraits::ValueType( 1.025830683,  2.488712048));
      check_matrix_element(mats.dB_kin, 1, 3, 2, QMCTraits::ValueType( 0.2538848112, 0.6159360987));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType(-0.003462536555, 0.005042847223));
      check_matrix_element(mats.dB_nlpp, 0, 1, 1, QMCTraits::ValueType( 0.000008297464266, 0.000138566606));
      check_matrix_element(mats.dB_nlpp, 0, 3, 0, QMCTraits::ValueType(-0.003958471461, 0.005765127996));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.01030887635, 0.01501387425));
      check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType(-0.0006555517895, -0.0015903985));
      check_matrix_element(mats.dB_nlpp, 1, 3, 0, QMCTraits::ValueType(-0.0003900410488, 0.0005680567787));
    }

    SECTION("mu=1 nu=2")
    {
      const int mu = 1;
      const int nu = 2;

      auto mats = extract_stress_matrices(sys.elec(), sys.psi(), sys.ham(), twf,
                                          mu, nu, KINETIC, NONLOCALECP);
      auto contractions = compute_stress_contractions(twf, mats);
      print_stress_contractions("jastrow", mu, nu, contractions);

      CHECK(std::real(contractions.tr_Minv_Bkin)    == Approx(8.719162572));
      CHECK(std::imag(contractions.tr_Minv_Bkin)    == Approx(-5.564999444e-07));
      CHECK(std::abs(std::real(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::abs(std::imag(contractions.tr_Minv_dM)) < 1e-12);
      CHECK(std::real(contractions.dBkin_minus_XdM) == Approx(1.091801783));
      CHECK(std::imag(contractions.dBkin_minus_XdM) == Approx(7.273143206e-07));
      CHECK(std::real(contractions.tr_Minv_Bnlpp)   == Approx(-0.001809233814));
      CHECK(std::imag(contractions.tr_Minv_Bnlpp)   == Approx(-8.42520724e-10));
      CHECK(std::real(contractions.dBnlpp_minus_XdM) == Approx(0.004683633961));
      CHECK(std::imag(contractions.dBnlpp_minus_XdM) == Approx(3.088295823e-09));

      check_matrix_element(mats.dB_kin, 0, 0, 0, QMCTraits::ValueType( 0.2751884822,-0.4007850418));
      check_matrix_element(mats.dB_kin, 0, 0, 2, QMCTraits::ValueType(-0.6858009681,-1.663784444));
      check_matrix_element(mats.dB_kin, 0, 3, 2, QMCTraits::ValueType(-0.5209571364,-1.26386575));

      check_matrix_element(mats.dB_kin, 1, 0, 0, QMCTraits::ValueType(-0.375043865, 0.5462146402));
      check_matrix_element(mats.dB_kin, 1, 1, 2, QMCTraits::ValueType(-0.8805968333,-2.13636815));
      check_matrix_element(mats.dB_kin, 1, 3, 3, QMCTraits::ValueType(-0.7771258513,-1.794190389));

      check_matrix_element(mats.dB_nlpp, 0, 0, 0, QMCTraits::ValueType( 0.007322662166,-0.01066474407));
      check_matrix_element(mats.dB_nlpp, 0, 1, 1, QMCTraits::ValueType(-0.0001353221575,-0.002259866654));
      check_matrix_element(mats.dB_nlpp, 0, 3, 0, QMCTraits::ValueType(-0.001011151348, 0.001472643398));

      check_matrix_element(mats.dB_nlpp, 1, 0, 0, QMCTraits::ValueType(-0.01066586367, 0.01553379138));
      check_matrix_element(mats.dB_nlpp, 1, 1, 2, QMCTraits::ValueType( 0.0005832389392, 0.001414964231));
      check_matrix_element(mats.dB_nlpp, 1, 3, 1, QMCTraits::ValueType(-0.0001693070321,-0.002827411262));
    }
  }
}

} // namespace qmcplusplus
