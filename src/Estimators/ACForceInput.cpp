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
#include "ACForceInput.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/UniformCommunicateError.h"
#include "EstimatorInput.h"
namespace qmcplusplus
{
ACForceInput::ACForceInput(xmlNodePtr cur)
{
  input_section_.readXML(cur);
  auto setIfInInput = LAMBDA_setIfInInput;
  setIfInInput(epsilon_, "epsilon");
  setIfInInput(use_spacewarp_,"spacewarp");

}

void ACForceInput::ACForceInputSection::checkParticularValidity()
{
}


};
