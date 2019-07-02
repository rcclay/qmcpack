//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <OhmmsSoA/Container.h>
#include "spline2/MultiBspline.hpp"
#include "spline2/MultiBsplineEval.hpp"
#include "QMCWaveFunctions/BsplineFactory/contraction_helper.hpp"

namespace qmcplusplus
{
template<typename T>
void test_spline_bounds()
{
  T x = 2.2;
  T dx;
  int ind;
  int ng = 10;
  spline2::getSplineBound(x, dx, ind, ng);
  REQUIRE(dx == Approx(0.2));
  REQUIRE(ind == 2);

  // check clamping to a maximum index value
  x = 10.5;
  spline2::getSplineBound(x, dx, ind, ng);
  REQUIRE(dx == Approx(0.5));
  REQUIRE(ind == 10);

  x = 11.5;
  spline2::getSplineBound(x, dx, ind, ng);
  REQUIRE(dx == Approx(1.0));
  REQUIRE(ind == 10);

  // check clamping to a zero index value
  x = -1.3;
  spline2::getSplineBound(x, dx, ind, ng);
  REQUIRE(dx == Approx(0.0));
  REQUIRE(ind == 0);
}

TEST_CASE("getSplineBound double", "[spline2]") { test_spline_bounds<double>(); }


TEST_CASE("getSplineBound float", "[spline2]") { test_spline_bounds<float>(); }


TEST_CASE("SymTrace", "[spline2]")
{
  // Code from here to the end of the function generate by gen_trace.py

  double h00 = 1;
  double h01 = 2;
  double h02 = 3;
  double h11 = 4.4;
  double h12 = 1.1;
  double h22 = 0.9;

  double gg[6] = {0.1, 1.6, 1.2, 2.3, 9.4, 2.3};

  double tr = SymTrace(h00, h01, h02, h11, h12, h22, gg);
  REQUIRE(tr == Approx(29.43));
}


template<typename T>
void test_prefactors()
{
  T a[4];
  T da[4];
  T d2a[4];
  T tx;

  // Code from here to the end the of function generated by gen_prefactor.py

  tx = 0.1;
  spline2::MultiBsplineData<T>::compute_prefactors(a, tx);
  REQUIRE(a[0] == Approx(0.1215));
  REQUIRE(a[1] == Approx(0.657167));
  REQUIRE(a[2] == Approx(0.221167));
  REQUIRE(a[3] == Approx(0.000166667));

  tx = 0.8;
  spline2::MultiBsplineData<T>::compute_prefactors(a, da, d2a, tx);
  REQUIRE(a[0] == Approx(0.00133333));
  REQUIRE(da[0] == Approx(-0.02));
  REQUIRE(d2a[0] == Approx(0.2));
  REQUIRE(a[1] == Approx(0.282667));
  REQUIRE(da[1] == Approx(-0.64));
  REQUIRE(d2a[1] == Approx(0.4));
  REQUIRE(a[2] == Approx(0.630667));
  REQUIRE(da[2] == Approx(0.34));
  REQUIRE(d2a[2] == Approx(-1.4));
  REQUIRE(a[3] == Approx(0.0853333));
  REQUIRE(da[3] == Approx(0.32));
  REQUIRE(d2a[3] == Approx(0.8));
}

TEST_CASE("double prefactors", "[spline2]") { test_prefactors<double>(); }

TEST_CASE("float prefactors", "[spline2]") { test_prefactors<float>(); }

/** Supports testing many sizes of splines for benchmarking
 *  modified from einspline/tests/test_3d.cpp
 */

template<typename T, int GRID_SIZE, int NUM_SPLINES>
struct test_splines_base
{
  const int num_splines = NUM_SPLINES;
  int npad;
  BCtype_d bc[3];
  Ugrid grid[3];

  int N = GRID_SIZE;
  double delta;
  std::vector<double> data;

  test_splines_base()
  {
    data.resize(N * N * N);
    npad = getAlignedSize<T>(num_splines);

    grid[0].start = 0.0;
    grid[0].end   = 1.0;
    grid[0].num   = N;

    grid[1].start = 0.0;
    grid[1].end   = 1.0;
    grid[1].num   = N;

    grid[2].start = 0.0;
    grid[2].end   = 1.0;
    grid[2].num   = N;

    delta = (grid[0].end - grid[0].start) / grid[0].num;

    bc[0].lCode = PERIODIC;
    bc[0].rCode = PERIODIC;
    bc[0].lVal  = 0.0;
    bc[0].rVal  = 0.0;
    bc[1].lCode = PERIODIC;
    bc[1].rCode = PERIODIC;
    bc[1].lVal  = 0.0;
    bc[1].rVal  = 0.0;
    bc[2].lCode = PERIODIC;
    bc[2].rCode = PERIODIC;
    bc[2].lVal  = 0.0;
    bc[2].rVal  = 0.0;

    double tpi = 2 * M_PI;
    // Generate the data in double precision regardless of the target spline precision
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++)
        {
          double x                    = delta * i;
          double y                    = delta * j;
          double z                    = delta * k;
          data[i * N * N + j * N + k] = std::sin(tpi * x) + std::sin(3 * tpi * y) + std::sin(4 * tpi * z);
        }
  }
};

/** Unspecialized test_splines doesn't test eval values
 */
template<typename T, int GRID_SIZE = 5, int NUM_SPLINES = 1>
struct test_splines : public test_splines_base<T, GRID_SIZE, NUM_SPLINES>
{
  using base = test_splines_base<T, GRID_SIZE, NUM_SPLINES>;
  using base::bc;
  using base::data;
  using base::delta;
  using base::grid;
  using base::N;
  using base::npad;
  using base::num_splines;

  test_splines() : test_splines_base<T, GRID_SIZE, NUM_SPLINES>(){};

  void test()
  {
    MultiBspline<T> bs;

    bs.create(grid, bc, npad);

    REQUIRE(bs.num_splines() == npad);

    double tpi = 2 * M_PI;

    BsplineAllocator<double> mAllocator;
    UBspline_3d_d* aspline = mAllocator.allocateUBspline(grid[0], grid[1], grid[2], bc[0], bc[1], bc[2], data.data());

    for (int i = 0; i < num_splines; i++)
      bs.copy_spline(aspline, i);

    mAllocator.destroy(aspline);

    //  The values for N=5 are not good enough for finer grids so by default we don't do those checks

    TinyVector<T, 3> pos = {0, 0, 0};

    aligned_vector<T> v(npad);
    spline2::evaluate3d(bs.getSplinePtr(), pos, v);

    VectorSoaContainer<T, 3> dv(npad);
    VectorSoaContainer<T, 6> hess(npad);
    spline2::evaluate3d_vgh(bs.getSplinePtr(), pos, v, dv, hess);

    pos = {0.1, 0.2, 0.3};
    spline2::evaluate3d(bs.getSplinePtr(), pos, v);

    spline2::evaluate3d_vgh(bs.getSplinePtr(), pos, v, dv, hess);

    VectorSoaContainer<T, 3> lap(npad);
    spline2::evaluate3d_vgl(bs.getSplinePtr(), pos, v, dv, lap);

    VectorSoaContainer<T, 10> ghess(npad);
    spline2::evaluate3d_vghgh(bs.getSplinePtr(), pos, v, dv, hess, ghess);
  }
};

/** Partially specialized test_splines does test eval values
 * See gen_bspline_values.py
 */
template<typename T>
struct test_splines<T, 5, 1> : public test_splines_base<T, 5, 1>
{
  using base = test_splines_base<T, 5, 1>;
  using base::bc;
  using base::data;
  using base::delta;
  using base::grid;
  using base::N;
  using base::npad;
  using base::num_splines;

  test_splines() : test_splines_base<T, 5, 1>(){};

  void test()
  {
    MultiBspline<T> bs;

    bs.create(grid, bc, npad);

    REQUIRE(bs.num_splines() == npad);

    BsplineAllocator<double> mAllocator;
    UBspline_3d_d* aspline = mAllocator.allocateUBspline(grid[0], grid[1], grid[2], bc[0], bc[1], bc[2], data.data());

    for (int i = 0; i < num_splines; i++)
      bs.copy_spline(aspline, i);

    mAllocator.destroy(aspline);

    //  Code from here to the end of the function is generated by gen_bspline_values.py

    TinyVector<T, 3> pos = {0, 0, 0};

    // symbolic value at pos =  (cx[0]/6 + 2*cx[1]/3 + cx[2]/6)*(cy[0]/6 + 2*cy[1]/3 + cy[2]/6)*(cz[0]/6 + 2*cz[1]/3 + cz[2]/6)

    aligned_vector<T> v(npad);
    spline2::evaluate3d(bs.getSplinePtr(), pos, v);
    REQUIRE(v[0] == Approx(-3.529930688e-12));

    VectorSoaContainer<T, 3> dv(npad);
    VectorSoaContainer<T, 6> hess(npad);
    spline2::evaluate3d_vgh(bs.getSplinePtr(), pos, v, dv, hess);
    // Gradient
    REQUIRE(dv[0][0] == Approx(6.178320809));
    REQUIRE(dv[0][1] == Approx(-7.402942564));
    REQUIRE(dv[0][2] == Approx(-6.178320809));

    // Hessian
    for (int i = 0; i < 6; i++)
    {
      REQUIRE(hess[0][i] == Approx(0.0));
    }

    pos = {0.1, 0.2, 0.3};
    spline2::evaluate3d(bs.getSplinePtr(), pos, v);

    // Value
    REQUIRE(v[0] == Approx(-0.9476393279));

    spline2::evaluate3d_vgh(bs.getSplinePtr(), pos, v, dv, hess);
    // Value
    REQUIRE(v[0] == Approx(-0.9476393279));
    // Gradient
    REQUIRE(dv[0][0] == Approx(5.111042137));
    REQUIRE(dv[0][1] == Approx(5.989106342));
    REQUIRE(dv[0][2] == Approx(1.952244379));
    // Hessian
    REQUIRE(hess[0][0] == Approx(-21.34557341));
    REQUIRE(hess[0][1] == Approx(1.174505743e-09));
    REQUIRE(hess[0][2] == Approx(-1.1483271e-09));
    REQUIRE(hess[0][3] == Approx(133.9204891));
    REQUIRE(hess[0][4] == Approx(-2.15319293e-09));
    REQUIRE(hess[0][5] == Approx(34.53786329));


    VectorSoaContainer<T, 3> lap(npad);
    spline2::evaluate3d_vgl(bs.getSplinePtr(), pos, v, dv, lap);
    // Value
    REQUIRE(v[0] == Approx(-0.9476393279));
    // Gradient
    REQUIRE(dv[0][0] == Approx(5.111042137));
    REQUIRE(dv[0][1] == Approx(5.989106342));
    REQUIRE(dv[0][2] == Approx(1.952244379));
    // Laplacian
    REQUIRE(lap[0][0] == Approx(147.1127789));

    VectorSoaContainer<T, 10> ghess(npad);
    spline2::evaluate3d_vghgh(bs.getSplinePtr(), pos, v, dv, hess, ghess);
    // Value
    REQUIRE(v[0] == Approx(-0.9476393279));
    // Gradient
    REQUIRE(dv[0][0] == Approx(5.111042137));
    REQUIRE(dv[0][1] == Approx(5.989106342));
    REQUIRE(dv[0][2] == Approx(1.952244379));
    // Hessian
    REQUIRE(hess[0][0] == Approx(-21.34557341));
    REQUIRE(hess[0][1] == Approx(1.174505743e-09));
    REQUIRE(hess[0][2] == Approx(-1.1483271e-09));
    REQUIRE(hess[0][3] == Approx(133.9204891));

    REQUIRE(hess[0][4] == Approx(-2.15319293e-09));
    REQUIRE(hess[0][5] == Approx(34.53786329));


    // Catch default is 100*(float epsilson)
    double eps = 2000 * std::numeric_limits<float>::epsilon();

    // Gradient of Hessian
    REQUIRE(ghess[0][0] == Approx(-213.455734));
    REQUIRE(ghess[0][1] == Approx(2.311193459e-09).epsilon(eps));
    REQUIRE(ghess[0][2] == Approx(3.468205279e-09).epsilon(eps));
    REQUIRE(ghess[0][3] == Approx(1.58092329e-07).epsilon(eps));
    REQUIRE(ghess[0][4] == Approx(1.255694171e-08).epsilon(eps));
    REQUIRE(ghess[0][5] == Approx(4.78981157e-08).epsilon(eps));
    REQUIRE(ghess[0][6] == Approx(-1753.041961));
    REQUIRE(ghess[0][7] == Approx(-2.575826885e-09).epsilon(eps));
    REQUIRE(ghess[0][8] == Approx(-4.683496702e-09).epsilon(eps));
    REQUIRE(ghess[0][9] == Approx(-81.53283531));
  }
};

TEST_CASE("MultiBspline periodic double", "[spline2]") { test_splines<double>().test(); }

/** Purpose is to evaluate performance
 *  this seems to be the largest spline mesh in the examples
 *  the '#pragma omp parallel for's in BslpineAllocator are worth about 
 *  a 16% speedup here on 32 threads vs serial
 *  for a 500x500x500 grid the speed up is a factor of 2
 *
 *  This test is only run if you explicitly call the .spline2 tag.
 */
TEST_CASE("MultiBspline periodic double large", "[.spline2]")
{
  test_splines<double, 128, 10> test;
  for (int i = 0; i < 20; i++)
    test.test();
}

TEST_CASE("MultiBspline periodic float", "[spline2]") { test_splines<float>().test(); }

} // namespace qmcplusplus
