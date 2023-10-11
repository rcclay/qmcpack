
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

#ifndef QMCPLUSPLUS_ACFORCEINPUT_H
#define QMCPLUSPLUS_ACFORCEINPUT_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Estimators/InputSection.h"
namespace qmcplusplus
{
class ACForceInput
{
public:
  using Real               = QMCTraits::RealType;
  static constexpr int DIM = QMCTraits::DIM;

  class ACForceInputSection : public InputSection
  {
  public:
    /** parse time definition of input parameters */
    ACForceInputSection()
    {
      //Note.  There is a change from legacy input declaration.  In particular, instantiation
      // was required by specifying "type="Force" mode="ACForce"".  This is because there was a 
      // ForceBase instantiation required.  Instead, we deprecate the "mode" attribute, and require
      // "type="ACForce"".  This should make the code simpler, but will break backwards compatibility
      // in the xml.

      // clang-format off
      section_name  = "ACForce";
      attributes    = {"name","epsilon","spacewarp"};
      parameters    = {};
      bools         = {"spacewarp"};
      enums         = {};
      strings       = {"name"};
      multi_strings = {};
      integers      = {};
      reals         = {"epsilon"};
      positions     = {};
      required      = {"name"};
      // I'd much rather see the default defined in simple native c++ as below
      // clang-format on
    }
    void checkParticularValidity() override;
  };

public:
  ACForceInput(xmlNodePtr node);
  /** default copy constructor
   *  This is required due to SDI being part of a variant used as a vector element.
   */
  ACForceInput(const ACForceInput&) = default;
  Real get_epsilon(){ return epsilon_;};
  bool get_use_spacewarp(){ return use_spacewarp_;};
  

private:
  ACForceInputSection input_section_;
  Real epsilon_ = 0.0;
  bool use_spacewarp_ = false;
  
};
} //qmcplusplus
#endif
