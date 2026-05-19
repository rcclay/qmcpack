//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "Message/Communicate.h"
#include "PWOrbitalSet.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{
PWOrbitalSet::~PWOrbitalSet()
{
  if (OwnBasisSet && myBasisSet)
    delete myBasisSet;
  if (!IsCloned && C != nullptr)
    delete C;
}

std::unique_ptr<SPOSet> PWOrbitalSet::makeClone() const
{
  auto myclone        = std::make_unique<PWOrbitalSet>(*this);
  myclone->myBasisSet = new PWBasis(*myBasisSet);
  myclone->IsCloned   = true;
  return myclone;
}

void PWOrbitalSet::resize(PWBasisPtr bset, bool cleanup)
{
  myBasisSet   = bset;
  OwnBasisSet  = cleanup;
  BasisSetSize = myBasisSet->NumPlaneWaves;
  C            = new ValueMatrix(OrbitalSetSize, BasisSetSize);
  Temp.resize(OrbitalSetSize, PW_MAXINDEX);
  app_log() << "  PWOrbitalSet::resize OrbitalSetSize =" << OrbitalSetSize << " BasisSetSize = " << BasisSetSize
            << std::endl;
}

void PWOrbitalSet::addVector(const std::vector<ComplexType>& coefs, int jorb)
{
  int ng = myBasisSet->inputmap.size();
  if (ng != coefs.size())
  {
    app_error() << "  Input G map does not match the basis size of wave functions " << std::endl;
    OHMMS::Controller->abort();
  }
  //drop G points for the given TwistAngle
  const std::vector<int>& inputmap(myBasisSet->inputmap);
  for (int ig = 0; ig < ng; ig++)
  {
    if (inputmap[ig] > -1)
      (*C)[jorb][inputmap[ig]] = coefs[ig];
  }
}

void PWOrbitalSet::addVector(const std::vector<RealType>& coefs, int jorb)
{
  int ng = myBasisSet->inputmap.size();
  if (ng != coefs.size())
  {
    app_error() << "  Input G map does not match the basis size of wave functions " << std::endl;
    OHMMS::Controller->abort();
  }
  //drop G points for the given TwistAngle
  const std::vector<int>& inputmap(myBasisSet->inputmap);
  for (int ig = 0; ig < ng; ig++)
  {
    if (inputmap[ig] > -1)
      (*C)[jorb][inputmap[ig]] = coefs[ig];
  }
}

void PWOrbitalSet::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  //Evaluate every orbital for particle iat.
  //Evaluate the basis-set at these coordinates:
  //myBasisSet->evaluate(P,iat);
  myBasisSet->evaluate(P.activeR(iat));
  MatrixOperators::product(*C, myBasisSet->Zv, psi);
}

void PWOrbitalSet::evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  //Evaluate the orbitals and derivatives for particle iat only.
  myBasisSet->evaluateAll(P, iat);
  MatrixOperators::product(*C, myBasisSet->Z, Temp);
  const ValueType* restrict tptr = Temp.data();
  for (int j = 0; j < OrbitalSetSize; j++, tptr += PW_MAXINDEX)
  {
    psi[j]   = tptr[PW_VALUE];
    d2psi[j] = tptr[PW_LAP];
    dpsi[j]  = GradType(tptr[PW_GRADX], tptr[PW_GRADY], tptr[PW_GRADZ]);
  }
}

void PWOrbitalSet::evaluate_notranspose(const ParticleSet& P,
                                        int first,
                                        int last,
                                        ValueMatrix& logdet,
                                        GradMatrix& dlogdet,
                                        ValueMatrix& d2logdet)
{
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    myBasisSet->evaluateAll(P, iat);
    MatrixOperators::product(*C, myBasisSet->Z, Temp);
    const ValueType* restrict tptr = Temp.data();
    for (int j = 0; j < OrbitalSetSize; j++, tptr += PW_MAXINDEX)
    {
      logdet(i, j)   = tptr[PW_VALUE];
      d2logdet(i, j) = tptr[PW_LAP];
      dlogdet(i, j)  = GradType(tptr[PW_GRADX], tptr[PW_GRADY], tptr[PW_GRADZ]);
    }
  }
}
void PWOrbitalSet::evaluateVGH(const ParticleSet& P,
                               int iat,
                               ValueVector& psi,
                               GradVector& dpsi,
                               HessVector& grad_grad_psi)
{
  myBasisSet->evaluateAll(P, iat);
  MatrixOperators::product(*C, myBasisSet->Z, Temp);

  const ValueType* restrict tptr = Temp.data();
  for (int j = 0; j < OrbitalSetSize; j++, tptr += PW_MAXINDEX)
  {
    psi[j] = tptr[PW_VALUE];

    dpsi[j] = GradType(tptr[PW_GRADX], tptr[PW_GRADY], tptr[PW_GRADZ]);

    grad_grad_psi[j](0, 0) = tptr[PW_HESS00];
    grad_grad_psi[j](0, 1) = tptr[PW_HESS01];
    grad_grad_psi[j](0, 2) = tptr[PW_HESS02];
    grad_grad_psi[j](1, 0) = tptr[PW_HESS10];
    grad_grad_psi[j](1, 1) = tptr[PW_HESS11];
    grad_grad_psi[j](1, 2) = tptr[PW_HESS12];
    grad_grad_psi[j](2, 0) = tptr[PW_HESS20];
    grad_grad_psi[j](2, 1) = tptr[PW_HESS21];
    grad_grad_psi[j](2, 2) = tptr[PW_HESS22];
  }
}

void PWOrbitalSet::evaluateVGHGH(const ParticleSet& P,
                                 int iat,
                                 ValueVector& psi,
                                 GradVector& dpsi,
                                 HessVector& grad_grad_psi,
                                 GGGVector& grad_grad_grad_psi)
{
  myBasisSet->evaluateAll(P, iat);
  MatrixOperators::product(*C, myBasisSet->Z, Temp);

  const ValueType* restrict tptr = Temp.data();
  for (int j = 0; j < OrbitalSetSize; j++, tptr += PW_MAXINDEX)
  {
    psi[j] = tptr[PW_VALUE];

    dpsi[j] = GradType(tptr[PW_GRADX], tptr[PW_GRADY], tptr[PW_GRADZ]);

    grad_grad_psi[j](0, 0) = tptr[PW_HESS00];
    grad_grad_psi[j](0, 1) = tptr[PW_HESS01];
    grad_grad_psi[j](0, 2) = tptr[PW_HESS02];
    grad_grad_psi[j](1, 0) = tptr[PW_HESS10];
    grad_grad_psi[j](1, 1) = tptr[PW_HESS11];
    grad_grad_psi[j](1, 2) = tptr[PW_HESS12];
    grad_grad_psi[j](2, 0) = tptr[PW_HESS20];
    grad_grad_psi[j](2, 1) = tptr[PW_HESS21];
    grad_grad_psi[j](2, 2) = tptr[PW_HESS22];

    grad_grad_grad_psi[j][0](0, 0) = tptr[PW_GHESS000];
    grad_grad_grad_psi[j][0](0, 1) = tptr[PW_GHESS001];
    grad_grad_grad_psi[j][0](0, 2) = tptr[PW_GHESS002];
    grad_grad_grad_psi[j][0](1, 0) = tptr[PW_GHESS010];
    grad_grad_grad_psi[j][0](1, 1) = tptr[PW_GHESS011];
    grad_grad_grad_psi[j][0](1, 2) = tptr[PW_GHESS012];
    grad_grad_grad_psi[j][0](2, 0) = tptr[PW_GHESS020];
    grad_grad_grad_psi[j][0](2, 1) = tptr[PW_GHESS021];
    grad_grad_grad_psi[j][0](2, 2) = tptr[PW_GHESS022];

    grad_grad_grad_psi[j][1](0, 0) = tptr[PW_GHESS100];
    grad_grad_grad_psi[j][1](0, 1) = tptr[PW_GHESS101];
    grad_grad_grad_psi[j][1](0, 2) = tptr[PW_GHESS102];
    grad_grad_grad_psi[j][1](1, 0) = tptr[PW_GHESS110];
    grad_grad_grad_psi[j][1](1, 1) = tptr[PW_GHESS111];
    grad_grad_grad_psi[j][1](1, 2) = tptr[PW_GHESS112];
    grad_grad_grad_psi[j][1](2, 0) = tptr[PW_GHESS120];
    grad_grad_grad_psi[j][1](2, 1) = tptr[PW_GHESS121];
    grad_grad_grad_psi[j][1](2, 2) = tptr[PW_GHESS122];

    grad_grad_grad_psi[j][2](0, 0) = tptr[PW_GHESS200];
    grad_grad_grad_psi[j][2](0, 1) = tptr[PW_GHESS201];
    grad_grad_grad_psi[j][2](0, 2) = tptr[PW_GHESS202];
    grad_grad_grad_psi[j][2](1, 0) = tptr[PW_GHESS210];
    grad_grad_grad_psi[j][2](1, 1) = tptr[PW_GHESS211];
    grad_grad_grad_psi[j][2](1, 2) = tptr[PW_GHESS212];
    grad_grad_grad_psi[j][2](2, 0) = tptr[PW_GHESS220];
    grad_grad_grad_psi[j][2](2, 1) = tptr[PW_GHESS221];
    grad_grad_grad_psi[j][2](2, 2) = tptr[PW_GHESS222];
  }
}

void PWOrbitalSet::evaluate_notranspose(const ParticleSet& P,
                                        int first,
                                        int last,
                                        ValueMatrix& logdet,
                                        GradMatrix& dlogdet,
                                        HessMatrix& grad_grad_logdet)
{
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector v(logdet[i], logdet.cols());
    GradVector g(dlogdet[i], dlogdet.cols());
    HessVector h(grad_grad_logdet[i], grad_grad_logdet.cols());
    evaluateVGH(P, iat, v, g, h);
  }
}

void PWOrbitalSet::evaluate_notranspose(const ParticleSet& P,
                                        int first,
                                        int last,
                                        ValueMatrix& logdet,
                                        GradMatrix& dlogdet,
                                        HessMatrix& grad_grad_logdet,
                                        GGGMatrix& grad_grad_grad_logdet)
{
  for (int iat = first, i = 0; iat < last; iat++, i++)
  {
    ValueVector v(logdet[i], logdet.cols());
    GradVector g(dlogdet[i], dlogdet.cols());
    HessVector h(grad_grad_logdet[i], grad_grad_logdet.cols());
    GGGVector gh(grad_grad_grad_logdet[i], grad_grad_grad_logdet.cols());
    evaluateVGHGH(P, iat, v, g, h, gh);
  }
}
} // namespace qmcplusplus
