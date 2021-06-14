// Modified from Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2021, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

//----------------------------------------------------------------------------//
// SurrealD4_btest
// testing of SurrealD class w/ 4 derivitives

#include "Surreal/SurrealD.h"

#include <ostream>
using namespace std;

#define Real                        double
#define BOOST_CHECK                 assert
#define BOOST_CHECK_EQUAL(A,B)      assert((A) == (B))
#define BOOST_AUTO_TEST_CASE( fun ) void fun()


//############################################################################//
//BOOST_AUTO_TEST_SUITE( SurrealD4_test_suite )


//----------------------------------------------------------------------------//
bool
chkSurrealD4( const SurrealD& z, Real v, Real d0, Real d1, Real d2, Real d3 )
{
  bool isEqual = true;
  if ((z.value() != v) || (z.deriv(0) != d0) || (z.deriv(1) != d1) || (z.deriv(2) != d2) || (z.deriv(3) != d3))
  {
    isEqual = false;
    cout << "actual (" << z << ")  expected "
         << "((" << v << ", " << d0 << " " << d1 << " " << d2 << " " << d3 << ")) "
         << "diff (" << v - z.value()
         << " " << d0 - z.deriv(0) << " " << d1 - z.deriv(1)
         << " " << d2 - z.deriv(2) << " " << d3 - z.deriv(3) << ")" << endl;
  }
#if 0
  else
  {
    cout << "(" << z << ")" << endl;
  }
#endif
  return isEqual;
}


bool
chkSurrealD4( const SurrealD& z, Real v, Real d0, Real d1, Real d2, Real d3, Real tol )
{
  bool isEqual = true;
  if ((abs(z.value() - v) > tol) ||
      (abs(z.deriv(0) - d0) > tol) ||
      (abs(z.deriv(1) - d1) > tol) ||
      (abs(z.deriv(2) - d2) > tol) ||
      (abs(z.deriv(3) - d3) > tol))
  {
    isEqual = false;
    cout << "actual (" << z << ")  expected "
         << "((" << v << ", " << d0 << " " << d1 << " " << d2 << " " << d3 << "))"
         << "diff (" << v - z.value()
         << " " << d0 - z.deriv(0) << " " << d1 - z.deriv(1)
         << " " << d2 - z.deriv(2) << " " << d3 - z.deriv(3) << endl;
  }
#if 0
  else
  {
    cout << "(" << z << ")" << endl;
  }
#endif
  return isEqual;
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( ctors_access )
{
  Real d[4] = {1, 2, 3, 4};
  SurrealD v0(1, d, 4);
  SurrealD v1(1, 2, 4);
  SurrealD v2(v1);
  SurrealD v3 = v1;
  SurrealD v4(0,0.,4);
  SurrealD v5 = 1;
  SurrealD v6;

  // size
  BOOST_CHECK_EQUAL( 4, v0.size() );
  BOOST_CHECK_EQUAL( 4, v1.size() );
  BOOST_CHECK_EQUAL( 4, v2.size() );
  BOOST_CHECK_EQUAL( 4, v3.size() );
  BOOST_CHECK_EQUAL( 4, v4.size() );
  BOOST_CHECK_EQUAL( 0, v5.size() );
  BOOST_CHECK_EQUAL( 0, v6.size() );

  // accessors
  BOOST_CHECK_EQUAL( 1, v0.value() );
  BOOST_CHECK_EQUAL( 1, v0.deriv() );
  BOOST_CHECK_EQUAL( 1, v0.deriv(0) );
  BOOST_CHECK_EQUAL( 2, v0.deriv(1) );
  BOOST_CHECK_EQUAL( 3, v0.deriv(2) );
  BOOST_CHECK_EQUAL( 4, v0.deriv(3) );

  BOOST_CHECK_EQUAL( 1, v1.value() );
  BOOST_CHECK_EQUAL( 2, v1.deriv(0) );
  BOOST_CHECK_EQUAL( 2, v1.deriv(1) );
  BOOST_CHECK_EQUAL( 2, v1.deriv(2) );
  BOOST_CHECK_EQUAL( 2, v1.deriv(3) );

  BOOST_CHECK_EQUAL( 1, v2.value() );
  BOOST_CHECK_EQUAL( 2, v2.deriv(0) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(1) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(2) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(3) );

  BOOST_CHECK_EQUAL( 1, v3.value() );
  BOOST_CHECK_EQUAL( 2, v3.deriv(0) );
  BOOST_CHECK_EQUAL( 2, v3.deriv(1) );
  BOOST_CHECK_EQUAL( 2, v3.deriv(2) );
  BOOST_CHECK_EQUAL( 2, v3.deriv(3) );

  BOOST_CHECK_EQUAL( 0, v4.value() );
  BOOST_CHECK_EQUAL( 0, v4.deriv(0) );
  BOOST_CHECK_EQUAL( 0, v4.deriv(1) );
  BOOST_CHECK_EQUAL( 0, v4.deriv(2) );
  BOOST_CHECK_EQUAL( 0, v4.deriv(3) );

  BOOST_CHECK_EQUAL( 1, v5.value() );
  BOOST_CHECK_EQUAL( 0, v6.value() );

  v2.value() = 3;
  BOOST_CHECK_EQUAL( 3, v2.value() );
  BOOST_CHECK_EQUAL( 2, v2.deriv(0) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(1) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(2) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(3) );
  v2.deriv() = 5;
  BOOST_CHECK_EQUAL( 3, v2.value() );
  BOOST_CHECK_EQUAL( 5, v2.deriv(0) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(1) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(2) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(3) );
  v2.deriv(0) = 1;
  BOOST_CHECK_EQUAL( 3, v2.value() );
  BOOST_CHECK_EQUAL( 1, v2.deriv(0) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(1) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(2) );
  BOOST_CHECK_EQUAL( 2, v2.deriv(3) );

  SurrealD v7(Real(4), Real(2), 4);
  BOOST_CHECK_EQUAL( 4, v7.value() );
  BOOST_CHECK_EQUAL( 2, v7.deriv(0) );
  BOOST_CHECK_EQUAL( 2, v7.deriv(1) );
  BOOST_CHECK_EQUAL( 2, v7.deriv(2) );
  BOOST_CHECK_EQUAL( 2, v7.deriv(3) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( assign_unary_accum )
{
  Real d[4] = {1, 2, 3, 4};
  SurrealD v1(1, d, 4);
  SurrealD v2(1, 0., 4);
  SurrealD v3(v1);
  SurrealD v4, v5;

  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 0, 0, 0, 0 ) );

  // assignment
  v3 = 2;
  BOOST_CHECK( chkSurrealD4( v3,  2, 0, 0, 0, 0 ) );

  v3 = v2;
  BOOST_CHECK( chkSurrealD4( v2,  1, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v3,  1, 0, 0, 0, 0 ) );

  v3 = v2 = 2;
  BOOST_CHECK( chkSurrealD4( v2,  2, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v3,  2, 0, 0, 0, 0 ) );

  v3 = v2 = v1;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  1, 1, 2, 3, 4 ) );

  // unary
  v3 = +v2;
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  1, 1, 2, 3, 4 ) );

  v3 = -v2;
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, -1, -1, -2, -3, -4 ) );

  // unary for variables without derivatives
  v4 = 2;
  v4 = +v4;
  BOOST_CHECK_EQUAL(  0, v4.size() );
  BOOST_CHECK_EQUAL(  2, v4.value() );

  v5 = 2;
  v5 = -v5;
  BOOST_CHECK_EQUAL(  0, v5.size() );
  BOOST_CHECK_EQUAL( -2, v5.value() );

  // binary accumulation
  v4 = v3 = v2 = v1;
  v5 = 4*v1;
  v2 += 3;
  v3 -= 3;
  v4 *= 3;
  v5 /= 2;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  4, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, -2, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v4,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v5,  2, 2, 4, 6, 8 ) );

  v4 = v3 = v2 = v1;
  v5 = 4*v1;
  v2 += Real(3);
  v3 -= Real(3);
  v4 *= Real(3);
  v5 /= Real(2);
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  4, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, -2, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v4,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v5,  2, 2, 4, 6, 8 ) );

  v4 = v3 = v2 = v1;
  v5 = 4*v1;
  v2 += v1;
  v3 -= v1;
  v4 *= v1;
  v5 /= v1;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  0, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v4,  1, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v5,  4, 0, 0, 0, 0 ) );

  SurrealD v10(2.);
  v4 = v3 = v2 = v1;
  v5 = 4*v10;
  v2 += v10;
  v3 -= v10;
  v4 *= v10;
  v5 /= v10;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  3, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, -1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v4,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v5,  4, 0, 0, 0, 0 ) );

  // binary accumulation w/ default constructor
  SurrealD v6, v7, v8, v9;
  v6 += v1;
  v7 -= v1;
  v8 *= v1;
  v9 /= v1;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v6,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v7, -1, -1, -2, -3, -4 ) );
  BOOST_CHECK( chkSurrealD4( v8,  0, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v9,  0, 0, 0, 0, 0 ) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( addsubtract )
{
  Real d[4] = {1, 2, 3, 4};
  SurrealD v1(1, d, 4);
  SurrealD v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v4,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v5,  1, 1, 2, 3, 4 ) );

  // binary +/- operators

  SurrealD v8 = v1 + v2;
  BOOST_CHECK_EQUAL( 4, v8.size() );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v8,  2, 2, 4, 6, 8 ) );

  SurrealD v9;
  v9 = v1 + v2;
  BOOST_CHECK_EQUAL( 4, v9.size() );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v9,  2, 2, 4, 6, 8 ) );

  v9 = 0;
  BOOST_CHECK_EQUAL( 0, v9.size() );

  v9 += v1 + v2;
  BOOST_CHECK_EQUAL( 4, v9.size() );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v9,  2, 2, 4, 6, 8 ) );

  v9 = 0;
  BOOST_CHECK_EQUAL( 0, v9.size() );

  v9 -= v1 + v2;
  BOOST_CHECK_EQUAL( 4, v9.size() );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v9, -2,-2,-4,-6,-8 ) );

  SurrealD v10, v11(1), v12(2);
  v10 = v11 + v12;
  BOOST_CHECK_EQUAL( 0, v10.size() );
  BOOST_CHECK( chkSurrealD4( v11,  1, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v12,  2, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v10,  3, 0, 0, 0, 0 ) );

  v10 += v11 + v12;
  BOOST_CHECK_EQUAL( 0, v10.size() );
  BOOST_CHECK( chkSurrealD4( v11,  1, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v12,  2, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v10,  6, 0, 0, 0, 0 ) );

  v10 -= v11 + v12;
  BOOST_CHECK_EQUAL( 0, v10.size() );
  BOOST_CHECK( chkSurrealD4( v11,  1, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v12,  2, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v10,  3, 0, 0, 0, 0 ) );

  v2 = v1;
  v3 = v1 + v2;
  v4 = v1 + v2 + v3;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );

  v3 = v1 - v2;
  v4 = v1 - v2 - v3;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  0, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v4,  0, 0, 0, 0, 0 ) );

  v3 = 3;
  v4 = v1 + v2 - v3;
  v5 = v1 - v2 + v3;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v4, -1, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v5,  3, 0, 0, 0, 0 ) );

  v4 = v3 = v2 = v1;
  v2 += v1;
  v3 += v1 + v2;
  v4 += v1 + v2 + v3;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  4, 4, 8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v4,  8, 8, 16, 24, 32 ) );

  v4 = v3 = v2 = v1;
  v2 -= v1;
  v3 -= v1 - v2;
  v4 -= v1 - v2 - v3;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  0, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v3,  0, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v4,  0, 0, 0, 0, 0 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = (v1 + v2) + (v3 + v4);
  BOOST_CHECK( chkSurrealD4( v5,  10, 10, 20, 30, 40 ) );
  v5 = (v1 + v2) + (v3 - v4);
  BOOST_CHECK( chkSurrealD4( v5,  2, 2, 4, 6, 8 ) );
  v5 = (v1 + v2) - (v3 + v4);
  BOOST_CHECK( chkSurrealD4( v5, -4, -4, -8, -12, -16 ) );
  v5 = (v1 + v2) - (v3 - v4);
  BOOST_CHECK( chkSurrealD4( v5,  4, 4, 8, 12, 16 ) );
  v5 = (v1 - v2) + (v3 + v4);
  BOOST_CHECK( chkSurrealD4( v5,  6, 6, 12, 18, 24 ) );
  v5 = (v1 - v2) + (v3 - v4);
  BOOST_CHECK( chkSurrealD4( v5, -2, -2, -4, -6, -8 ) );
  v5 = (v1 - v2) - (v3 + v4);
  BOOST_CHECK( chkSurrealD4( v5, -8, -8, -16, -24, -32 ) );
  v5 = (v1 - v2) - (v3 - v4);
  BOOST_CHECK( chkSurrealD4( v5,  0, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );

  v5 += (v1 + v2) + (v3 + v4);
  v5 += (v1 + v2) + (v3 - v4);
  v5 += (v1 + v2) - (v3 + v4);
  v5 += (v1 + v2) - (v3 - v4);
  v5 += (v1 - v2) + (v3 + v4);
  v5 += (v1 - v2) + (v3 - v4);
  v5 += (v1 - v2) - (v3 + v4);
  v5 += (v1 - v2) - (v3 - v4);
  BOOST_CHECK( chkSurrealD4( v5,  8, 8, 16, 24, 32 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );

  v5 -= (v1 + v2) + (v3 + v4);
  v5 -= (v1 + v2) + (v3 - v4);
  v5 -= (v1 + v2) - (v3 + v4);
  v5 -= (v1 + v2) - (v3 - v4);
  v5 -= (v1 - v2) + (v3 + v4);
  v5 -= (v1 - v2) + (v3 - v4);
  v5 -= (v1 - v2) - (v3 + v4);
  v5 -= (v1 - v2) - (v3 - v4);
  BOOST_CHECK( chkSurrealD4( v5,  0, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );

  v2 = +v1;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1,  1,  2,  3,  4 ) );
  v2 = -v1;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v2, -1, -1, -2, -3, -4 ) );
  v3 = +(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v2, -1, -1, -2, -3, -4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  0,  0,  0,  0,  0 ) );
  v3 = +(v1 - v2);
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v2, -1, -1, -2, -3, -4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  2,  2,  4,  6,  8 ) );
  v3 = -(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v2, -1, -1, -2, -3, -4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  0,  0,  0,  0,  0 ) );
  v3 = -(v1 - v2);
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v2, -1, -1, -2, -3, -4 ) );
  BOOST_CHECK( chkSurrealD4( v3, -2, -2, -4, -6, -8 ) );

  // addition/subtraction with scalar quantities

  v3 = v1 + 3;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  4,  1,  2,  3,  4 ) );
  v3 = v1 + Real(3);
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  4,  1,  2,  3,  4 ) );
  v3 = 3 + v1;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  4,  1,  2,  3,  4 ) );
  v3 = Real(3) + v1;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  4,  1,  2,  3,  4 ) );
  v3 += v1 + 1;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  6,  2,  4,  6,  8 ) );
  v3 -= v1 + 1;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  4,  1,  2,  3,  4 ) );
  v3 = v1 - 3;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3, -2,  1,  2,  3,  4 ) );
  v3 = v1 - Real(3);
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3, -2,  1,  2,  3,  4 ) );
  v3 = 3 - v1;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  2, -1, -2, -3, -4 ) );
  v3 = Real(3) - v1;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  2, -1, -2, -3, -4 ) );

  v3 = +(v1 + v2) + 3;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v2, -1, -1, -2, -3, -4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3,  0,  0,  0,  0 ) );
  v3 = -(v1 + v2) + 3;
  BOOST_CHECK( chkSurrealD4( v1,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v2, -1, -1, -2, -3, -4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3,  0,  0,  0,  0 ) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( multiply )
{
  Real d[4] = {1, 2, 3, 4};
  SurrealD v1(1, d, 4);
  SurrealD v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v4,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v5,  1, 1, 2, 3, 4 ) );

  // binary * operators

  v2 = 3*v1;
  v3 = v2*2;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v3,  6, 6, 12, 18, 24 ) );

  v2 += 3*v1;
  v3 += v2*2;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  6, 6, 12, 18, 24 ) );
  BOOST_CHECK( chkSurrealD4( v3, 18, 18, 36, 54, 72 ) );

  v2 -= 3*v1;
  v3 -= v2*2;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v3, 12, 12, 24, 36, 48 ) );

  v2 = 3/v1;
  v3 = v2/1;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  3, -3, -6, -9, -12 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, -3, -6, -9, -12 ) );

  v2 = 2*v1;
  v3 = v1*v2;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  2, 4, 8, 12, 16 ) );
  v3 += v1*v2;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  4, 8, 16, 24, 32 ) );

  v3 = 2*v1;
  v3 = v1*v3; //Test when v3 is both on left and right
  BOOST_CHECK( chkSurrealD4( v3, 2, 4, 8, 12, 16 ) );
  v3 += v1*v3;
  BOOST_CHECK( chkSurrealD4( v3, 4, 10, 20, 30, 40 ) );
  v3 -= v1*v3;
  BOOST_CHECK( chkSurrealD4( v3, 0, -4, -8, -12, -16 ) );

  v2 = 2*(v1*2);
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  4, 4, 8, 12, 16 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = 3*(v1 + v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v5, 23, 23, 46, 69, 92 ) );
  v5 = 3*(v1 + v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5,  7, 7, 14, 21, 28 ) );
  v5 = 3*(v1 + v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -5, -5, -10, -15, -20 ) );
  v5 = 3*(v1 + v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, 11, 11, 22, 33, 44 ) );
  v5 = 3*(v1 - v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, 11, 11, 22, 33, 44 ) );
  v5 = 3*(v1 - v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -5, -5, -10, -15, -20 ) );
  v5 = 3*(v1 - v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -17, -17, -34, -51, -68 ) );
  v5 = 3*(v1 - v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -1, -1, -2, -3, -4 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = 5;
  v5 += 3*(v1 + v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v5, 28, 23, 46, 69, 92 ) );
  v5 += 3*(v1 + v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, 35, 30, 60, 90, 120 ) );
  v5 += 3*(v1 + v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, 30, 25, 50, 75, 100 ) );
  v5 += 3*(v1 + v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, 41, 36, 72, 108, 144 ) );
  v5 += 3*(v1 - v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, 52, 47, 94, 141, 188 ) );
  v5 += 3*(v1 - v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, 47, 42, 84, 126, 168 ) );
  v5 += 3*(v1 - v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, 30, 25, 50, 75, 100 ) );
  v5 += 3*(v1 - v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, 29, 24, 48, 72, 96 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = 5;
  v5 -= 3*(v1 + v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v5, -18, -23, -46, -69, -92 ) );
  v5 -= 3*(v1 + v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -25, -30, -60, -90, -120 ) );
  v5 -= 3*(v1 + v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -20, -25, -50, -75, -100 ) );
  v5 -= 3*(v1 + v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -31, -36, -72, -108, -144 ) );
  v5 -= 3*(v1 - v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -42, -47, -94, -141, -188 ) );
  v5 -= 3*(v1 - v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -37, -42, -84, -126, -168 ) );
  v5 -= 3*(v1 - v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -20, -25, -50, -75, -100 ) );
  v5 -= 3*(v1 - v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealD4( v5, -19, -24, -48, -72, -96 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  4, 4, 8, 12, 16 ) );

  v3 = 3*(v1 + v2)*2;
  BOOST_CHECK( chkSurrealD4( v3, 18, 18, 36, 54, 72 ) );
  v3 = 3*2*(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v3, 18, 18, 36, 54, 72 ) );
  v3 = (v1 + v2)*3*2;
  BOOST_CHECK( chkSurrealD4( v3, 18, 18, 36, 54, 72 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );

  v2 = +3*v1;
  BOOST_CHECK( chkSurrealD4( v2,  3,  3,  6,  9,  12 ) );
  v2 = -3*v1;
  BOOST_CHECK( chkSurrealD4( v2, -3, -3, -6, -9, -12 ) );
  v2 = +v1*3;
  BOOST_CHECK( chkSurrealD4( v2,  3,  3,  6,  9,  12 ) );
  v2 = -v1*3;
  BOOST_CHECK( chkSurrealD4( v2, -3, -3, -6, -9, -12 ) );
  v2 = +(3*v1);
  BOOST_CHECK( chkSurrealD4( v2,  3,  3,  6,  9,  12 ) );
  v2 = -(3*v1);
  BOOST_CHECK( chkSurrealD4( v2, -3, -3, -6, -9, -12 ) );
  v2 = +(v1*3);
  BOOST_CHECK( chkSurrealD4( v2,  3,  3,  6,  9,  12 ) );
  v2 = -(v1*3);
  BOOST_CHECK( chkSurrealD4( v2, -3, -3, -6, -9, -12 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = v1*v2*v3;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4,  6, 18, 36, 54, 72 ) );
  v4 = 2*v1*v2*v3;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4, 12, 36, 72, 108, 144 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = (v1 + v2)*v3;
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  4, 4, 8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v4, 12, 24, 48, 72, 96 ) );
  v4 += (v1 + v2)*v3;
  BOOST_CHECK( chkSurrealD4( v4, 24, 48, 96, 144, 192 ) );
  v4 = v3*(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v4, 12, 24, 48, 72, 96 ) );
  v4 += v3*(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v4, 24, 48, 96, 144, 192 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  4, 4, 8, 12, 16 ) );

  v2 = 2*v1;
  v3 = (v1 + v2)*(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  9, 18, 36, 54, 72 ) );
  v2 = 2*v1;
  v3 = 3*v1;
  v4 = (v1 + v2)*(v2 + v3);
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
  BOOST_CHECK( chkSurrealD4( v4, 15, 30, 60, 90, 120 ) );
  v4 += (v1 + v2)*(v2 + v3);
  BOOST_CHECK( chkSurrealD4( v4, 30, 60, 120, 180, 240 ) );
  v4 = 2*(v1 + v2)*(v2 + v3);
  BOOST_CHECK( chkSurrealD4( v4, 30, 60, 120, 180, 240 ) );
  BOOST_CHECK( chkSurrealD4( v1,  1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  3, 3, 6, 9, 12 ) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( divide )
{
  Real d[4] = {2, 4, 6, 8};
  SurrealD v1(2, d, 4);
  SurrealD v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealD4( v1,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v4,  2, 2, 4, 6, 8 ) );
  BOOST_CHECK( chkSurrealD4( v5,  2, 2, 4, 6, 8 ) );

  // binary / operators

  v2 = 4/v1;
  v3 = v2/2;
  BOOST_CHECK( chkSurrealD4( v1,  2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, -2, -4, -6, -8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  1, -1, -2, -3, -4 ) );

  v2 = Real(4)/v1;
  v3 = v2/Real(2);
  BOOST_CHECK( chkSurrealD4( v1,  2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,  2, -2, -4, -6, -8 ) );
  BOOST_CHECK( chkSurrealD4( v3,  1, -1, -2, -3, -4 ) );

  v2 += 4/v1;
  v3 += v2/2;
  BOOST_CHECK( chkSurrealD4( v1, 2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2, 4, -4, -8, -12, -16 ) );
  BOOST_CHECK( chkSurrealD4( v3, 3, -3, -6, -9, -12 ) );

  v2 -= 4/v1;
  v3 -= v2/2;
  BOOST_CHECK( chkSurrealD4( v1, 2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2, -2, -4, -6, -8 ) );
  BOOST_CHECK( chkSurrealD4( v3, 2, -2, -4, -6, -8 ) );

  v2 = 2*v1;
  v3 = v2/v1;
  BOOST_CHECK( chkSurrealD4( v1, 2, 2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2, 4, 4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3, 2, 0,  0,  0,  0 ) );

  v3 += v2/v1;
  BOOST_CHECK( chkSurrealD4( v1, 2, 2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2, 4, 4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3, 4, 0,  0,  0,  0 ) );

  v3 -= v2/v1;
  BOOST_CHECK( chkSurrealD4( v1, 2, 2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2, 4, 4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3, 2, 0,  0,  0,  0 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = (v2+v1)/v3;
  BOOST_CHECK( chkSurrealD4( v1, 2, 2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2, 4, 4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3, 6, 6, 12, 18, 24 ) );
  BOOST_CHECK( chkSurrealD4( v4, 1, 0,  0,  0,  0 ) );
  v4 += (v2+v1)/v3;
  BOOST_CHECK( chkSurrealD4( v4, 2, 0,  0,  0,  0 ) );
  v4 = v3/(v2+v1);
  BOOST_CHECK( chkSurrealD4( v4, 1, 0,  0,  0,  0 ) );
  v4 += v3/(v2+v1);
  BOOST_CHECK( chkSurrealD4( v4, 2, 0,  0,  0,  0 ) );
  BOOST_CHECK( chkSurrealD4( v1, 2, 2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2, 4, 4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3, 6, 6, 12, 18, 24 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = 6*v1;
  v5 = 12/(v1 + v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v1,   2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,   4,  4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3,   8,  8, 16, 24, 32 ) );
  BOOST_CHECK( chkSurrealD4( v4,  12, 12, 24, 36, 48 ) );
  BOOST_CHECK( chkSurrealD4( v5,  12,  8, 16, 24, 32 ) );
  v5 = 12/(v1 + v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,   0,  -4,  -8, -12, -16 ) );
  v5 = 12/(v1 + v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  -8, -12, -24, -36, -48 ) );
  v5 = 12/(v1 + v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,   4,   0,   0,   0,   0) );
  v5 = 12/(v1 - v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,   4,  16,  32,  48,  64 ) );
  v5 = 12/(v1 - v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  -8,   4,   8,  12,  16 ) );
  v5 = 12/(v1 - v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v5, -16,  -4,  -8, -12, -16 ) );
  v5 = 12/(v1 - v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  -4,   8,  16,  24,  32 ) );
  BOOST_CHECK( chkSurrealD4( v1,   2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,   4,  4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3,   8,  8, 16, 24, 32 ) );
  BOOST_CHECK( chkSurrealD4( v4,  12, 12, 24, 36, 48 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = 6*v1;
  v5 = 1;
  v5 += 12/(v1 + v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v1,   2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,   4,  4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3,   8,  8, 16, 24, 32 ) );
  BOOST_CHECK( chkSurrealD4( v4,  12, 12, 24, 36, 48 ) );
  BOOST_CHECK( chkSurrealD4( v5,  13,  8, 16, 24, 32 ) );
  v5 += 12/(v1 + v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  13,  4,  8, 12, 16 ) );
  v5 += 12/(v1 + v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,   5, -8, -16, -24, -32 ) );
  v5 += 12/(v1 + v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,   9, -8, -16, -24, -32 ) );
  v5 += 12/(v1 - v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  13,  8,  16,  24,  32 ) );
  v5 += 12/(v1 - v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,   5, 12,  24,  36,  48 ) );
  v5 += 12/(v1 - v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v5, -11,  8,  16,  24,  32 ) );
  v5 += 12/(v1 - v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5, -15, 16,  32,  48,  64 ) );
  BOOST_CHECK( chkSurrealD4( v1,   2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,   4,  4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3,   8,  8, 16, 24, 32 ) );
  BOOST_CHECK( chkSurrealD4( v4,  12, 12, 24, 36, 48 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = 6*v1;
  v5 = 1;
  v5 -= 12/(v1 + v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v1,   2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,   4,  4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3,   8,  8, 16, 24, 32 ) );
  BOOST_CHECK( chkSurrealD4( v4,  12, 12, 24, 36, 48 ) );
  BOOST_CHECK( chkSurrealD4( v5, -11, -8, -16, -24, -32 ) );
  v5 -= 12/(v1 + v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5, -11, -4, -8, -12, -16 ) );
  v5 -= 12/(v1 + v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  -3,  8, 16, 24, 32 ) );
  v5 -= 12/(v1 + v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  -7,  8, 16, 24, 32 ) );
  v5 -= 12/(v1 - v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v5, -11, -8, -16, -24, -32 ) );
  v5 -= 12/(v1 - v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  -3, -12, -24, -36, -48 ) );
  v5 -= 12/(v1 - v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  13,  -8, -16, -24, -32 ) );
  v5 -= 12/(v1 - v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealD4( v5,  17, -16, -32, -48, -64 ) );
  BOOST_CHECK( chkSurrealD4( v1,   2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,   4,  4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3,   8,  8, 16, 24, 32 ) );
  BOOST_CHECK( chkSurrealD4( v4,  12, 12, 24, 36, 48 ) );

  v5 = 12/(v1 + v2)/2;
  BOOST_CHECK( chkSurrealD4( v5,  1, -1, -2, -3, -4 ) );
  v5 = 12/2/(v1 - v2);
  BOOST_CHECK( chkSurrealD4( v5, -3,  3,  6,  9, 12 ) );
  v5 = (v1 + v2)/3/2;
  BOOST_CHECK( chkSurrealD4( v5,  1,  1,  2,  3,  4 ) );
  BOOST_CHECK( chkSurrealD4( v1,  2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,  4,  4,  8, 12, 16 ) );

  // cppcheck-suppress duplicateExpression
  v5 = (v1 + v2)/(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v5,  1, 0, 0, 0, 0 ) );
  v5 = 2*(v1 + v2)/(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v5,  2, 0, 0, 0, 0 ) );
  v5 += 2*(v1 + v2)/(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v5,  4, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v1,  2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,  4,  4,  8, 12, 16 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v5 = (v2 + v3)/(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v5,  2, 0, 0, 0, 0 ) );
  v5 = 2*(v2 + v3)/(v1 + v2);
  BOOST_CHECK( chkSurrealD4( v5,  4, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v1,  2,  2,  4,  6,  8 ) );
  BOOST_CHECK( chkSurrealD4( v2,  4,  4,  8, 12, 16 ) );
  BOOST_CHECK( chkSurrealD4( v3,  8,  8, 16, 24, 32 ) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( relational )
{
  SurrealD v1(1, 3, 4);
  SurrealD v2(1, 3, 4);
  SurrealD v3(2, 3, 4);

  BOOST_CHECK(   v2 == v1  );
  BOOST_CHECK( !(v3 == v1) );
  BOOST_CHECK(   v2 == 1   );
  BOOST_CHECK( !(v3 == 1)  );
  BOOST_CHECK(   1 == v2   );
  BOOST_CHECK( !(1 == v3)  );

  BOOST_CHECK( v2 == Real(1) );
  BOOST_CHECK( Real(1) == v2 );

  BOOST_CHECK( !(v2 != v1) );
  BOOST_CHECK(   v3 != v1  );
  BOOST_CHECK( !(v2 != 1)  );
  BOOST_CHECK(   v3 != 1   );
  BOOST_CHECK( !(1 != v2)  );
  BOOST_CHECK(   1 != v3   );

  BOOST_CHECK( !(v2 > v1) );
  BOOST_CHECK(   v3 > v1  );
  BOOST_CHECK( !(v2 > 1)  );
  BOOST_CHECK(   v3 > 1   );
  BOOST_CHECK( !(1 > v2)  );
  BOOST_CHECK( !(1 > v3)  );

  BOOST_CHECK( !(v2 < v1) );
  BOOST_CHECK( !(v3 < v1) );
  BOOST_CHECK( !(v2 < 1)  );
  BOOST_CHECK( !(v3 < 1)  );
  BOOST_CHECK( !( 1 < v2) );
  BOOST_CHECK(    1 < v3  );

  BOOST_CHECK(  v2 >= v1  );
  BOOST_CHECK(  v3 >= v1  );
  BOOST_CHECK(  v2 >= 1   );
  BOOST_CHECK(  v3 >= 1   );
  BOOST_CHECK(   1 >= v2  );
  BOOST_CHECK( !(1 >= v3) );

  BOOST_CHECK(  v2 <= v1 );
  BOOST_CHECK( !(v3 <= v1) );
  BOOST_CHECK(  v2 <= 1 );
  BOOST_CHECK( !(v3 <= 1) );
  BOOST_CHECK(  1 <= v2 );
  BOOST_CHECK(  1 <= v3 );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( cmath )
{
  const double tol = 1.e-13;
  Real d[4] = {1, 2, 3, 4};
  SurrealD v1(1, d, 4);
  SurrealD v2(v1), v3(v1), v4(v1), v5(v1);
  Real d0;

  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v4, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v5, 1, 1, 2, 3, 4 ) );

  // trig functions <cmath>

  v2 = cos(v1);  d0 = -sin(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, cos(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 += cos(v1);  d0 = -2*sin(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2*cos(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 = sin(v1);  d0 = cos(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, sin(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 += sin(v1);  d0 = 2*cos(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2*sin(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 = tan(v1);  d0 = 1/(cos(1.)*cos(1.));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, tan(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 += tan(v1);  d0 = 2/(cos(1.)*cos(1.));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2*tan(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );

  v2 = acos(cos(v1));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4, tol ) );
  v2 += acos(cos(v1));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2, 2, 4, 6, 8, tol ) );
  v2 = asin(sin(v1));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4, tol ) );
  v2 += asin(sin(v1));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2, 2, 4, 6, 8, tol ) );
  v2 = atan(tan(v1));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4, tol ) );
  v2 += atan(tan(v1));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2, 2, 4, 6, 8, tol ) );

  v2 = v1;
  v3 = v1;
  v4 = atan2(v2, v3);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v4, atan2(1., 1.), 0, 0, 0, 0, tol ) );
  v4 += atan2(v2, v3);
  BOOST_CHECK( chkSurrealD4( v4, 2*atan2(1., 1.), 0, 0, 0, 0, tol ) );

  // hyperbolic functions <cmath>

  v2 = cosh(v1);  d0 = sinh(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, cosh(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 += cosh(v1);  d0 = 2*sinh(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2*cosh(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 = sinh(v1);  d0 = cosh(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, sinh(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 += sinh(v1);  d0 = 2*cosh(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2*sinh(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 = tanh(v1);  d0 = 1/(cosh(1.)*cosh(1.));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, tanh(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 += tanh(v1);  d0 = 2/(cosh(1.)*cosh(1.));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2*tanh(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );

  // exp and log functions <cmath>

  v2 = exp(v1);  d0 = exp(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, exp(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 += exp(v1);  d0 = 2*exp(1.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2*exp(1.), d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 = log(v1);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 0, 1, 2, 3, 4, tol ) );
  v2 += log(v1);  d0 = 2;
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 0, d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 = log10(v1);  d0 = 1/log(10.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 0, d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 += log10(v1);  d0 = 2/log(10.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 0, d0, 2*d0, 3*d0, 4*d0, tol ) );

  // power functions <cmath>

  v2 = v1;
  v3 = pow(v2, v1);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 1, 1, 2, 3, 4, tol ) );
  v3 += pow(v2, v1);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 2, 2, 4, 6, 8, tol ) );
  v2 = pow(v1, 2);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 2, 4, 6, 8, tol ) );
  v2 = pow(v1, Real(2));
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 2, 4, 6, 8, tol ) );
  v2 += pow(v1, 2);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2, 4, 8, 12, 16, tol ) );
  v2 = pow(2, v1);  d0 = 2*log(2.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2, d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 = pow(Real(2), v1);  d0 = 2*log(2.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2, d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 += pow(2, v1);  d0 = 4*log(2.);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 4, d0, 2*d0, 3*d0, 4*d0, tol ) );
  v2 = v1;
  v3 = pow(v1+v2, v1+v2);  d0 = 13.54517744447956;
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 4, d0, 2*d0, 3*d0, 4*d0, tol ) );
  v3 += pow(v1+v2, v1+v2);  d0 = 2*13.54517744447956;
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 8, d0, 2*d0, 3*d0, 4*d0, tol ) );

  v2 = v1;
  v3 = pow(4*v1/v2, 0.5 );
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 2, 0, 0, 0, 0, tol ) );
  v3 += pow(4*v1/v2, 0.5 );
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 4, 0, 0, 0, 0, tol ) );

  v2 = pow(v1, 0);
  BOOST_CHECK( chkSurrealD4( v2, 1, 0, 0, 0, 0 ) );
  v2 = pow(0, v1);
  BOOST_CHECK( chkSurrealD4( v2, 0, 0, 0, 0, 0 ) );
  v2 = SurrealD(0, d, 4);
  v3 = pow(v2, 0);
  BOOST_CHECK( chkSurrealD4( v2, 0, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 1, 0, 0, 0, 0 ) );
  v3 = pow(v2, 1);
  BOOST_CHECK( chkSurrealD4( v2, 0, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 0, 1, 2, 3, 4 ) );
  Real d3[4] = {4, 3, 2, 1};
  v3 = SurrealD(0, d3, 4);
  v4 = pow(v2, v3);
  BOOST_CHECK( chkSurrealD4( v2, 0, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 0, 4, 3, 2, 1 ) );
  BOOST_CHECK( chkSurrealD4( v4, 1, 0, 0, 0, 0 ) );


  v2 = sqrt(v1);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 0.5, 1, 1.5, 2, tol ) );
  v2 += sqrt(v1);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 2, 1, 2, 3, 4, tol ) );
  v2 = 0;
  v3 = sqrt(v2);
  BOOST_CHECK( chkSurrealD4( v2, 0, 0, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealD4( v3, 0, 0, 0, 0, 0, tol ) );
  v3 = sqrt(4*v1);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 2, 1, 2, 3, 4, tol ) );
  v3 += sqrt(4*v1);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v3, 4, 2, 4, 6, 8, tol ) );

  // rounding functions <cmath>

  v2 = 1.5*v1;
  v3 = ceil(v2);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1.5, 1.5, 3, 4.5, 6 ) );
  BOOST_CHECK( chkSurrealD4( v3, 2, 0, 0, 0, 0, tol ) );
  v3 = floor(v2);
  BOOST_CHECK( chkSurrealD4( v2, 1.5, 1.5, 3, 4.5, 6 ) );
  BOOST_CHECK( chkSurrealD4( v3, 1, 0, 0, 0, 0, tol ) );

  // misc functions <cmath>

  v2 = abs(v1);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
  v2 = fabs(v1);
  BOOST_CHECK( chkSurrealD4( v1, 1, 1, 2, 3, 4 ) );
  BOOST_CHECK( chkSurrealD4( v2, 1, 1, 2, 3, 4 ) );
}

//----------------------------------------------------------------------------//
/*
BOOST_AUTO_TEST_CASE( IO )
{
  //Set the 2nd argument to false to regenerate the pattern file
  output_test_stream output( "IO/Surreal/SurrealD4_pattern.txt", true );

  Real d[4] = {1, 2, 3, 4};
  SurrealD v1(1, d, 4);
  SurrealD v2;

  output << v1 << std::endl;
  BOOST_CHECK( output.match_pattern() );
  v1.dump( 2, output );
  BOOST_CHECK( output.match_pattern() );

  output << v2 << std::endl;
  BOOST_CHECK( output.match_pattern() );
  v2.dump( 2, output );
  BOOST_CHECK( output.match_pattern() );

  output << v1 + v2 << std::endl;
  BOOST_CHECK( output.match_pattern() );
  v2.dump( 2, output );
  BOOST_CHECK( output.match_pattern() );

}
 */
//############################################################################//
// BOOST_AUTO_TEST_SUITE_END()
int main(int argc, char *argv[])
{
  ctors_access();
  assign_unary_accum();
  addsubtract();
  multiply();
  divide();
  relational();
  cmath();
  //  IO();
  std::cout << std::endl;
  std::cout << "SurrealD4_test_suite Complete!" << std::endl;
  std::cout << std::endl;
  return 0;
}
