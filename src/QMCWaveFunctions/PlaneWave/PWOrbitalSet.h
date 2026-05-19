//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file PWOrbitalSet.h
 * @brief Definition of member functions of Plane-wave basis set
 */
#ifndef QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H
#define QMCPLUSPLUS_PLANEWAVE_ORBITALSET_BLAS_H

#include "QMCWaveFunctions/PlaneWave/PWBasis.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "CPU/BLAS.hpp"

namespace qmcplusplus
{
class PWOrbitalSet : public SPOSet
{
public:
  using BasisSet_t = PWBasis;
  using PWBasisPtr = PWBasis*;

  /** inherit the enum of BasisSet_t */
  enum
  {
    PW_VALUE    = BasisSet_t::PW_VALUE,
    PW_LAP      = BasisSet_t::PW_LAP,
    PW_GRADX    = BasisSet_t::PW_GRADX,
    PW_GRADY    = BasisSet_t::PW_GRADY,
    PW_GRADZ    = BasisSet_t::PW_GRADZ,

    PW_HESS00   = BasisSet_t::PW_HESS00,
    PW_HESS01   = BasisSet_t::PW_HESS01,
    PW_HESS02   = BasisSet_t::PW_HESS02,
    PW_HESS10   = BasisSet_t::PW_HESS10,
    PW_HESS11   = BasisSet_t::PW_HESS11,
    PW_HESS12   = BasisSet_t::PW_HESS12,
    PW_HESS20   = BasisSet_t::PW_HESS20,
    PW_HESS21   = BasisSet_t::PW_HESS21,
    PW_HESS22   = BasisSet_t::PW_HESS22,

    PW_GHESS000 = BasisSet_t::PW_GHESS000,
    PW_GHESS001 = BasisSet_t::PW_GHESS001,
    PW_GHESS002 = BasisSet_t::PW_GHESS002,
    PW_GHESS010 = BasisSet_t::PW_GHESS010,
    PW_GHESS011 = BasisSet_t::PW_GHESS011,
    PW_GHESS012 = BasisSet_t::PW_GHESS012,
    PW_GHESS020 = BasisSet_t::PW_GHESS020,
    PW_GHESS021 = BasisSet_t::PW_GHESS021,
    PW_GHESS022 = BasisSet_t::PW_GHESS022,

    PW_GHESS100 = BasisSet_t::PW_GHESS100,
    PW_GHESS101 = BasisSet_t::PW_GHESS101,
    PW_GHESS102 = BasisSet_t::PW_GHESS102,
    PW_GHESS110 = BasisSet_t::PW_GHESS110,
    PW_GHESS111 = BasisSet_t::PW_GHESS111,
    PW_GHESS112 = BasisSet_t::PW_GHESS112,
    PW_GHESS120 = BasisSet_t::PW_GHESS120,
    PW_GHESS121 = BasisSet_t::PW_GHESS121,
    PW_GHESS122 = BasisSet_t::PW_GHESS122,

    PW_GHESS200 = BasisSet_t::PW_GHESS200,
    PW_GHESS201 = BasisSet_t::PW_GHESS201,
    PW_GHESS202 = BasisSet_t::PW_GHESS202,
    PW_GHESS210 = BasisSet_t::PW_GHESS210,
    PW_GHESS211 = BasisSet_t::PW_GHESS211,
    PW_GHESS212 = BasisSet_t::PW_GHESS212,
    PW_GHESS220 = BasisSet_t::PW_GHESS220,
    PW_GHESS221 = BasisSet_t::PW_GHESS221,
    PW_GHESS222 = BasisSet_t::PW_GHESS222,

    PW_MAXINDEX = BasisSet_t::PW_MAXINDEX
  };
  /** default constructor
  */
  PWOrbitalSet(const std::string& my_name, size_t size)
      : SPOSet(my_name, size), OwnBasisSet(false), myBasisSet(nullptr), BasisSetSize(0), C(nullptr), IsCloned(false)
  {}

  std::string getClassName() const override { return "PWOrbitalSet"; }


  /** delete BasisSet only it owns this
   *
   * Builder takes care of who owns what
   */
  ~PWOrbitalSet() override;

  std::unique_ptr<SPOSet> makeClone() const override;
  /** resize  the orbital base
   * @param bset PWBasis
   * @param cleaup if true, owns PWBasis. Will clean up.
   */
  void resize(PWBasisPtr bset, bool cleanup = false);

  /** Builder class takes care of the assertion
  */
  void addVector(const std::vector<ComplexType>& coefs, int jorb);
  void addVector(const std::vector<RealType>& coefs, int jorb);


  inline ValueType evaluate(int ib, const PosType& pos)
  {
    myBasisSet->evaluate(pos);
    return BLAS::dot(BasisSetSize, (*C)[ib], myBasisSet->Zv.data());
  }

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;
  void evaluateVGH(const ParticleSet& P,
                   int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   HessVector& grad_grad_psi) override;

  void evaluateVGHGH(const ParticleSet& P,
                     int iat,
                     ValueVector& psi,
                     GradVector& dpsi,
                     HessVector& grad_grad_psi,
                     GGGVector& grad_grad_grad_psi) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet,
                            GGGMatrix& grad_grad_grad_logdet) override;
  /** boolean
   *
   * If true, this has to delete the BasisSet
   */
  bool OwnBasisSet;
  ///TwistAngle of this PWOrbitalSet
  PosType TwistAngle;
  ///My basis set
  PWBasisPtr myBasisSet;
  ///number of basis
  IndexType BasisSetSize;
  /** pointer to matrix containing the coefficients
   *
   * makeClone makes a shallow copy and flag IsCloned
   */
  ValueMatrix* C;
  ///if true, do not clean up
  bool IsCloned;
  /////Plane-wave coefficients: (iband,g-vector)
  //Matrix<ValueType> Coefs;
  /** temporary array to perform gemm operation */
  Matrix<ValueType> Temp;
};
} // namespace qmcplusplus
#endif
