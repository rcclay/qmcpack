//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_FATWALKERTYPES_H
#define QMCPLUSPLUS_FATWALKERTYPES_H

#include <ResourceCollection.h>

namespace qmcplusplus
{
struct FatWalkerResourceCollection
{
  ResourceCollection pset_res;
  ResourceCollection twf_res;

  FatWalkerResourceCollection() : pset_res("ParticleSet"), twf_res("TrialWaveFunction") {}
};

class FatWalkerResourceCollectionLock
{
public:
  FatWalkerResourceCollectionLock(FatWalkerResourceCollection& fatwalker_res, ParticleSet& pset, TrialWaveFunction& twf)
      : pset_res_lock_(fatwalker_res.pset_res, pset), twf_res_lock_(fatwalker_res.twf_res, twf)
  {}

private:
  ResourceCollectionLock<ParticleSet> pset_res_lock_;
  ResourceCollectionLock<TrialWaveFunction> twf_res_lock_;
};
} // namespace qmcplusplus
#endif
