//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Raymond Clay, Sandia National Laboratories
//
// File created by: Raymond Clay, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "ACStress.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#include "CPU/VectorOps.h"

namespace qmcplusplus
{
constexpr std::array<std::pair<int, int>, 6> symm_pairs = {{
    {0, 0},
    {1, 1},
    {2, 2},
    {0, 1},
    {0, 2},
    {1, 2}
}};

ACStress::ACStress(ParticleSet& target, TrialWaveFunction& psi_in, QMCHamiltonian& H)
    : ham_(H),
      first_stress_index_(-1),
      reg_epsilon_(0.0),
      f_epsilon_(1.0)
{
  setName("ACStress");

  hf_stress_      = 0.0;
  pulay_stress_   = 0.0;
  wf_strain_grad_ = 0.0;

  psi_in.getOrCreateTWFFastDerivWrapper(target);
}

std::unique_ptr<OperatorBase> ACStress::makeClone(ParticleSet& qp, TrialWaveFunction& psi) const
{
  APP_ABORT("ACStress::makeClone(ParticleSet&,TrialWaveFunction&) shouldn't be called");
  return nullptr;
}

std::unique_ptr<OperatorBase> ACStress::makeClone(ParticleSet& qp, TrialWaveFunction& psi_in, QMCHamiltonian& ham_in) const
{
  std::unique_ptr<ACStress> myclone = std::make_unique<ACStress>(qp, psi_in, ham_in);
  myclone->first_stress_index_      = first_stress_index_;
  return myclone;
}
bool ACStress::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attr;
  attr.add(reg_epsilon_, "epsilon");
  attr.put(cur);

  if (reg_epsilon_ < 0)
    throw std::runtime_error("ACStress::put(): epsilon<0 not allowed.");

  return true;
}

bool ACStress::get(std::ostream& os) const
{
  os << "ACStress (fast-derivative pathway only)";
  return true;
}

void ACStress::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH) const
{
  std::unique_ptr<OperatorBase> myclone = makeClone(qp, psi, targetH);
  if (myclone)
    targetH.addOperator(std::move(myclone), name_, update_mode_[PHYSICAL]);
}

ACStress::Return_t ACStress::evaluate(TrialWaveFunction& psi, ParticleSet& P)
{
  hf_stress_      = 0.0;
  pulay_stress_   = 0.0;
  wf_strain_grad_ = 0.0;

  TWFFastDerivWrapper& psiwrapper_ = psi.getOrCreateTWFFastDerivWrapper(P);

  // For now, interpret the fast Hamiltonian kernel outputs as:
  //   hf_stress_      <- operator derivative contribution d/dε (O Psi / Psi)
  //   wf_strain_grad_ <- d/dε log(Psi)
  //
  // pulay_stress_ can then be formed later at the observable/reporting level
  // if desired as:
  //   pulay = hf - E * wf_grad
  //
  // but in this first implementation we store the raw outputs separately.
  ham_.evaluateStrainDerivsFast(P, psi, psiwrapper_, hf_stress_, wf_strain_grad_);

  // In this first fast-only implementation, pulay_stress_ is not separately
  // decomposed by the Hamiltonian kernel.  Leave it zero for now and report the
  // remaining AC-style pieces explicitly in setObservables().
  pulay_stress_ = 0.0;

  f_epsilon_ = compute_regularizer_f(psi.G, reg_epsilon_);
  return 0.0;
}


void ACStress::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if (first_stress_index_ < 0)
    first_stress_index_ = plist.size();

  for (const auto& [mu, nu] : symm_pairs)
  {
    const std::string muStr(std::to_string(mu));
    const std::string nuStr(std::to_string(nu));

    const std::string hfname("ACStress_hf_" + muStr + "_" + nuStr);
    const std::string pulayname("ACStress_pulay_" + muStr + "_" + nuStr);
    const std::string wfgradname1("ACStress_Ewfgrad_" + muStr + "_" + nuStr);
    const std::string wfgradname2("ACStress_wfgrad_" + muStr + "_" + nuStr);

    plist.add(hfname);
    plist.add(pulayname);
    plist.add(wfgradname1);
    plist.add(wfgradname2);
  }
}

ACStress::RealType ACStress::compute_regularizer_f(const ParticleSet::ParticleGradient& G, const RealType epsilon)
{
  if (std::abs(epsilon) < 1e-6)
    return 1.0;

  RealType gdotg = 0.0;
#if defined(QMC_COMPLEX)
  gdotg = Dot_CC(G, G);
#else
  gdotg = Dot(G, G);
#endif

  RealType gmag = std::sqrt(gdotg);
  double xovereps = 1.0 / (epsilon * gmag);

  RealType regvalue = 0.0;
  if (xovereps >= 1.0)
    regvalue = 1.0;
  else
    regvalue = 7.0 * std::pow(xovereps, 6.0)
             - 15.0 * std::pow(xovereps, 4.0)
             +  9.0 * std::pow(xovereps, 2.0);

  return regvalue;
}

void ACStress::setObservables(PropertySetType& plist)
{
  int myindex = first_stress_index_;

  for (const auto& [mu, nu] : symm_pairs)
  {
    RealType hf_val     = 0.0;
    RealType pulay_val  = 0.0;
    RealType wfgrad_val = 0.0;

    if (mu == nu)
    {
      hf_val     = hf_stress_(mu, nu);
      pulay_val  = pulay_stress_(mu, nu);
      wfgrad_val = wf_strain_grad_(mu, nu);
    }
    else
    {
      hf_val =
          0.5 * (hf_stress_(mu, nu) + hf_stress_(nu, mu));
      pulay_val =
          0.5 * (pulay_stress_(mu, nu) + pulay_stress_(nu, mu));
      wfgrad_val =
          0.5 * (wf_strain_grad_(mu, nu) + wf_strain_grad_(nu, mu));
    }

    plist[myindex++] = hf_val * f_epsilon_;
    plist[myindex++] = pulay_val * f_epsilon_;
    plist[myindex++] = ham_.getLocalEnergy() * wfgrad_val * f_epsilon_;
    plist[myindex++] = wfgrad_val * f_epsilon_;
  }
}

void ACStress::setParticlePropertyList(PropertySetType& plist, int offset)
{
  int myindex = first_stress_index_ + offset;

  for (const auto& [mu, nu] : symm_pairs)
  {
    RealType hf_val     = 0.0;
    RealType pulay_val  = 0.0;
    RealType wfgrad_val = 0.0;

    if (mu == nu)
    {
      hf_val     = hf_stress_(mu, nu);
      pulay_val  = pulay_stress_(mu, nu);
      wfgrad_val = wf_strain_grad_(mu, nu);
    }
    else
    {
      hf_val =
          0.5 * (hf_stress_(mu, nu) + hf_stress_(nu, mu));
      pulay_val =
          0.5 * (pulay_stress_(mu, nu) + pulay_stress_(nu, mu));
      wfgrad_val =
          0.5 * (wf_strain_grad_(mu, nu) + wf_strain_grad_(nu, mu));
    }

    plist[myindex++] = hf_val * f_epsilon_;
    plist[myindex++] = pulay_val * f_epsilon_;
    plist[myindex++] = ham_.getLocalEnergy() * wfgrad_val * f_epsilon_;
    plist[myindex++] = wfgrad_val * f_epsilon_;
  }
}

} // namespace qmcplusplus
