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

namespace qmcplusplus
{

ACStress::ACStress(ParticleSet& target, TrialWaveFunction& psi_in, QMCHamiltonian& H)
    : ham_(H), first_stress_index_(-1)
{
  setName("ACStress");

  hf_stress_       = 0.0;
  pulay_stress_    = 0.0;
  wf_strain_grad_  = 0.0;

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
  // No tunable options in the first implementation
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

  return 0.0;
}

void ACStress::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if (first_stress_index_ < 0)
    first_stress_index_ = plist.size();

  for (int mu = 0; mu < OHMMS_DIM; ++mu)
    for (int nu = 0; nu < OHMMS_DIM; ++nu)
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

void ACStress::setObservables(PropertySetType& plist)
{
  int myindex = first_stress_index_;

  for (int mu = 0; mu < OHMMS_DIM; ++mu)
    for (int nu = 0; nu < OHMMS_DIM; ++nu)
    {
      // Sign convention consistent with stress = -(1/V) dE/dε is left to the
      // downstream interpretation / normalization layer.
      // For now, follow the ACForce style and store the raw pieces.
      plist[myindex++] = hf_stress_(mu, nu);
      plist[myindex++] = pulay_stress_(mu, nu);
      plist[myindex++] = ham_.getLocalEnergy() * wf_strain_grad_(mu, nu);
      plist[myindex++] = wf_strain_grad_(mu, nu);
    }
}

void ACStress::setParticlePropertyList(PropertySetType& plist, int offset)
{
  int myindex = first_stress_index_ + offset;

  for (int mu = 0; mu < OHMMS_DIM; ++mu)
    for (int nu = 0; nu < OHMMS_DIM; ++nu)
    {
      plist[myindex++] = hf_stress_(mu, nu);
      plist[myindex++] = pulay_stress_(mu, nu);
      plist[myindex++] = ham_.getLocalEnergy() * wf_strain_grad_(mu, nu);
      plist[myindex++] = wf_strain_grad_(mu, nu);
    }
}

} // namespace qmcplusplus
