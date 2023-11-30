//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/BsplineFactory/createBsplineReader.h"
#include <PlatformSelector.hpp>
#include "SplineSetReader.h"
#include "HybridRepSetReader.h"

namespace qmcplusplus
{
std::unique_ptr<BsplineReader> createBsplineRealDouble(EinsplineSetBuilder* e,
                                                           bool hybrid_rep,
                                                           const std::string& useGPU)
{
  app_summary() << "    Using real valued spline SPOs with real double precision storage (R2R)." << std::endl;
  if (CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET)
    app_summary() << "OpenMP offload has not been implemented to support real valued spline SPOs with real storage!"
                  << std::endl;
  app_summary() << "    Running on CPU." << std::endl;

  std::unique_ptr<BsplineReader> aReader;
  if (hybrid_rep)
  {
    app_summary() << "    Using hybrid orbital representation." << std::endl;
    aReader = std::make_unique<HybridRepSetReader<HybridRepReal<SplineR2R<double>>>>(e);
  }
  else
    aReader = std::make_unique<SplineSetReader<SplineR2R<double>>>(e);
  return aReader;
}
} // namespace qmcplusplus
