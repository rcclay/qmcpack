//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef COULOMB_POT_H
#define COULOMB_POT_H

#include "PotentialBase.h"

class CoulombPot : public Potential
{
public:
  double Z1Z2;

  double V(double r) override;
  double dVdr(double r) override;
  double d2Vdr2(double r) override;
  double X_k(double rcut, double k) override;
  bool NeedsRel() override;
  void Write(IOSectionClass& out) override;
  void Read(IOSectionClass& in) override;
};

#endif
