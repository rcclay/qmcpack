//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "Configuration.h"
#include "Estimators/ACForceInput.h"
#include "Estimators/ACForce.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/WaveFunctionTypes.hpp"


using std::string;


namespace qmcplusplus
{
namespace testing
{
class ACForceTests
{
  using WF       = WaveFunctionTypes<QMCTraits::ValueType, QMCTraits::FullPrecValueType>;
  using Position = QMCTraits::PosType;
  using Data     = ACForce::Data;
  using Real     = WF::Real;
  using Value    = WF::Value;

public:
  void testCopyConstructor(const ACForce& magdens)
  {
  }

  void testData(const ACForce& magdens, const Data& data)
  {
  }

};
} //namespace testing

TEST_CASE("ACForce::IntegrationTest", "[estimators]")
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!!   Evaluate ACForce   !!!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

  //For now, do a small square case.
  const int nelec       = 2;
  const int norb        = 2;
  using WF              = WaveFunctionTypes<QMCTraits::ValueType, QMCTraits::FullPrecValueType>;
  using Real            = WF::Real;
  using Value           = WF::Value;
  using Grad            = WF::Grad;
  using ValueVector     = Vector<Value>;
  using GradVector      = Vector<Grad>;
  using ValueMatrix     = Matrix<Value>;
  using PropertySetType = OperatorBase::PropertySetType;
  using MCPWalker       = Walker<QMCTraits, PtclOnLatticeTraits>;
  using Data            = ACForce::Data;
  using GradMatrix      = Matrix<Grad>;
  using namespace testing;
}
} // namespace qmcplusplus
