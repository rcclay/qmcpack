//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by:  Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/UniformCommunicateError.h"

#include <iostream>
#include "ACForceInput.h"
namespace qmcplusplus
{
TEST_CASE("ACForceInput::from_xml", "[estimators]")
{
  using Real = QMCTraits::RealType;
  std::string_view input_xml = R"XML(<estimator name="doop" type="ACForce" epsilon="0.01" spacewarp="yes"/>)XML";

  Libxml2Document doc;
  bool okay      = doc.parseFromString(input_xml);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  ACForceInput acfin(node);
 
  Real eps=0.01;
  bool use_spacewarp=true;

  CHECK( eps == Approx(acfin.get_epsilon()));
  CHECK( use_spacewarp == acfin.get_use_spacewarp()); 
  
}

} // namespace qmcplusplus
