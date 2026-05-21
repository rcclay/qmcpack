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

#ifndef QMCPLUSPLUS_ACSTRESS_H
#define QMCPLUSPLUS_ACSTRESS_H

#include "QMCHamiltonians/OperatorBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "OhmmsPETE/Tensor.h"

namespace qmcplusplus
{
/** Assaraf-Caffarel style fast-derivative stress estimator
 *
 *  This first implementation is restricted to the fast-derivative pathway and
 *  assumes a Slater-Jastrow wavefunction. It reports the full unsymmetrized
 *  stress tensor.
 *
 *  Unlike most OperatorBase derived classes, it requires call Hamiltonian functions.
 *  To workaround the limitation of OperatorBase, it captures a QMCHamiltonian object
 *  reference via the constructor or add2Hamiltonian which wraps the usual makeClone
 *  function with a QMCHamiltonian object reference added.
 */
class ACStress : public OperatorBase
{
public:
  using StressTensor = Tensor<RealType, OHMMS_DIM>;

  /** Constructor **/
  ACStress(ParticleSet& target, TrialWaveFunction& psi, QMCHamiltonian& H);

  /** Destructor **/
  ~ACStress() override = default;

  std::string getClassName() const override { return "ACStress"; }

  /** I/O */
  bool put(xmlNodePtr cur) final;
  bool get(std::ostream& os) const final;

  /** Cloning **/
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) const final;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& H) const;

  /** observables */
  void addObservables(PropertySetType& plist, BufferType& collectables) final;
  void setObservables(PropertySetType& plist) final;
  void setParticlePropertyList(PropertySetType& plist, int offset) final;

  /** Since we store a reference to QMCHamiltonian, the baseclass method add2Hamiltonian
   *  isn't sufficient.  We override it here.
   */
  void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH) const final;

  /** Evaluate **/
  Return_t evaluate(TrialWaveFunction& psi, ParticleSet& P) final;

  /// accessors for debugging / testing
  const StressTensor& getHFStress() const { return hf_stress_; }
  const StressTensor& getPulayStress() const { return pulay_stress_; }
  const StressTensor& getWFStressGrad() const { return wf_strain_grad_; }

private:
  QMCHamiltonian& ham_;

  /// for indexing observables
  IndexType first_stress_index_;

  /// Hellmann-Feynman / operator-derivative contribution
  StressTensor hf_stress_;
  /// Pulay / mixed determinant derivative contribution
  StressTensor pulay_stress_;
  /// d/dε log(Psi)
  StressTensor wf_strain_grad_;
};

} // namespace qmcplusplus
#endif
