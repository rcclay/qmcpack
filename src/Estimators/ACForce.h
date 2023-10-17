//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ACFORCE_H
#define QMCPLUSPLUS_ACFORCE_H

#include "Estimators/OperatorEstBase.h"
#include "Estimators/ACForceInput.h"
#include "type_traits/complex_help.hpp"

namespace qmcplusplus
{

class ACForce : public OperatorEstBase
{
public:
  using Value              = QMCTraits::ValueType;
  using FullPrecValue      = QMCTraits::FullPrecValueType;
  using Real               = RealAlias<Value>;
  using FullPrecReal       = RealAlias<FullPrecValue>;
  using Grad               = TinyVector<Value, OHMMS_DIM>;
  using Lattice            = PtclOnLatticeTraits::ParticleLayout;
  using Position           = QMCTraits::PosType;
  static constexpr int DIM = QMCTraits::DIM;

  ACForce(ACForceInput&& acin);
  ACForce(const ACForce& magdens, DataLocality dl);

  void startBlock(int steps) override;

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomBase<FullPrecReal>& rng) override;

  void collect(const RefVector<OperatorEstBase>& operator_estimators) override;
  /**
  * This returns the total size of data object required for this estimator.
  * Right now, it's 3*number_of_realspace_gridpoints
  *
  * @return Size of data.
  */
  size_t getFullDataSize();
  std::unique_ptr<OperatorEstBase> spawnCrowdClone() const override;
  void registerOperatorEstimator(hdf_archive& file) override;


private:
  ACForce(const ACForce& magdens) = default;

  ACForceInput input_;
};
} // namespace qmcplusplus
#endif

