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
#include "Estimators/ACForce.h"
#include "Estimators/ACForceInput.h"

namespace qmcplusplus
{
ACForce::ACForce(ACForceInput&& acfinput)
    : OperatorEstBase(DataLocality::crowd), input_(acfinput)
{
  my_name_ = "ACForce";
  data_.resize(getFullDataSize(), 0);
}

ACForce::ACForce(const ACForce& magdens, DataLocality dl)
    : ACForce(magdens)
{
  my_name_       = "ACForce";
  data_locality_ = dl;
}
void ACForce::startBlock(int steps){};

size_t ACForce::getFullDataSize() { return 0; }

void ACForce::accumulate(const RefVector<MCPWalker>& walkers,
                                      const RefVector<ParticleSet>& psets,
                                      const RefVector<TrialWaveFunction>& wfns,
                                      const RefVector<QMCHamiltonian>& hamss,
                                      RandomBase<FullPrecReal>& rng)
{
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    MCPWalker& walker      = walkers[iw];
    ParticleSet& pset      = psets[iw];
    TrialWaveFunction& wfn = wfns[iw];

    QMCT::RealType weight = walker.Weight;
    assert(weight >= 0);

    walkers_weight_ += weight;
    //Ray:  Fill in
  }
};

void ACForce::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  if (data_locality_ == DataLocality::crowd)
  {
    OperatorEstBase::collect(type_erased_operator_estimators);
  }
  else
  {
    throw std::runtime_error("You cannot call collect on ACForce with this DataLocality");
  }
}


std::unique_ptr<OperatorEstBase> ACForce::spawnCrowdClone() const
{
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;

  //Everyone else has this attempt to set up a non-implemented memory saving optimization.
  //We won't rock the boat.
  if (data_locality_ == DataLocality::rank)
  {
    // This is just a stub until a memory saving optimization is deemed necessary
    spawn_data_locality = DataLocality::queue;
    data_size           = 0;
    throw std::runtime_error("There is no memory savings implementation for ACForce");
  }

  UPtr<ACForce> spawn(std::make_unique<ACForce>(*this, spawn_data_locality));
  spawn->get_data().resize(data_size);
  return spawn;
};


void ACForce::registerOperatorEstimator(hdf_archive& file)
{
  std::vector<size_t> my_indexes;

  std::vector<int> ng(1, getFullDataSize());

  hdf_path hdf_name{my_name_};
  h5desc_.emplace_back(hdf_name);
  auto& oh = h5desc_.back();
  oh.set_dimensions(ng, 0);
}



} //namespace qmcplusplus
