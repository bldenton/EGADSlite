// Modified from Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2021, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

//----------------------------------------------------------------------------//
// SurrealS1_btest
// testing of SurrealS<1> class

#include "Surreal/SurrealS.h"

#include <ostream>
using namespace std;

#define Real                        double
#define BOOST_CHECK                 assert
#define BOOST_CHECK_EQUAL(A,B)      assert((A) == (B))
#define BOOST_AUTO_TEST_CASE( fun ) void fun()

// Explicitly instantiate the class to generate all the functions so that coverage
// information is correct
template class SurrealS<1>;


//############################################################################//
//BOOST_AUTO_TEST_SUITE( SurrealS1_test_suite )


typedef SurrealS<1> SurrealS1;


//----------------------------------------------------------------------------//
bool
chkSurrealS1( const SurrealS1& z, Real v, Real d )
{
  bool isEqual = true;
  if ((z.value() != v) || (z.deriv() != d))
  {
    isEqual = false;
    cout << "actual (" << z << ")  "
         << "expected ((" << v << ";" << d << "))  "
         << "diff ((" << v - z.value() << ";" << d - z.deriv() << "))" << endl;
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
bool
chkSurrealS1( const SurrealS<1,SurrealS1>& z, Real v, Real d, Real dd )
{
  bool isEqual = true;
  if ((z.value() != v) || (z.deriv().value() != d) || (z.deriv().deriv() != dd))
  {
    isEqual = false;
    cout << "actual (" << z << ")  "
         << "expected ((" << v << ";" << d << ";" << dd << "))  "
         << "diff ((" << v - z.value() << ";(" << d - z.deriv().value() << ";" << dd - z.deriv().deriv() << ")))" << endl;
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
chkSurrealS1( const SurrealS1& z, Real v, Real d, Real tol )
{
  bool isEqual = true;
  if ((abs(z.value() - v) > tol) || (abs(z.deriv() - d) > tol))
  {
    isEqual = false;
    cout << "actual (" << z << ")  "
         << "expected ((" << v << ";" << d << "))  "
         << "diff ((" << v - z.value() << ";" << d - z.deriv() << "))" << endl;
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
chkSurrealS1( const SurrealS<1,SurrealS1>& z, Real v, Real d, Real dd, Real tol )
{
  bool isEqual = true;
  if ((abs(z.value() - v) > tol) || (abs(z.deriv().value() - d) > tol) || (abs(z.deriv().deriv() - dd) > tol) )
  {
    isEqual = false;
    cout << "actual (" << z << ")  "
         << "expected ((" << v << ";" << d << ";" << dd << "))  "
         << "diff ((" << v - z.value() << ";" << d - z.deriv().value() << ";" << dd - z.deriv().deriv() << "))" << endl;
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
  Real d[1] = {1};
  SurrealS1 v0(4, d, 1);
  SurrealS1 v1(1, 2);
  SurrealS1 v2(v1);
  SurrealS1 v3 = v1;
  SurrealS1 v4 = 1;
  SurrealS1 v5;

  // size
  BOOST_CHECK_EQUAL( 1, v0.size() );
  BOOST_CHECK_EQUAL( 1, v1.size() );
  BOOST_CHECK_EQUAL( 1, v2.size() );
  BOOST_CHECK_EQUAL( 1, v3.size() );
  BOOST_CHECK_EQUAL( 1, v4.size() );
  BOOST_CHECK_EQUAL( 1, v5.size() );

  // accessors
  BOOST_CHECK_EQUAL( 4, v0.value() );
  BOOST_CHECK_EQUAL( 1, v0.deriv() );
  BOOST_CHECK_EQUAL( 1, v0.deriv(0) );

  BOOST_CHECK_EQUAL( 1, v1.value() );
  BOOST_CHECK_EQUAL( 2, v1.deriv() );

  BOOST_CHECK_EQUAL( 1, v2.value() );
  BOOST_CHECK_EQUAL( 2, v2.deriv() );

  BOOST_CHECK_EQUAL( 1, v3.value() );
  BOOST_CHECK_EQUAL( 2, v3.deriv() );

  BOOST_CHECK_EQUAL( 1, v4.value() );
  BOOST_CHECK_EQUAL( 0, v4.deriv() );

//  BOOST_CHECK_EQUAL( 0, v5.value() );
//  BOOST_CHECK_EQUAL( 0, v5.deriv() );

  v2.value() = 3;
  BOOST_CHECK_EQUAL( 3, v2.value() );
  BOOST_CHECK_EQUAL( 2, v2.deriv() );
  v2.deriv() = 5;
  BOOST_CHECK_EQUAL( 3, v2.value() );
  BOOST_CHECK_EQUAL( 5, v2.deriv() );
  v2.deriv(0) = 2;
  BOOST_CHECK_EQUAL( 3, v2.value() );
  BOOST_CHECK_EQUAL( 2, v2.deriv() );

  SurrealS1 v7(Real(4.), Real(2.));
  BOOST_CHECK_EQUAL( 4, v7.value() );
  BOOST_CHECK_EQUAL( 2, v7.deriv() );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( assign_unary_accum )
{
  Real d[1] = {2};
  SurrealS1 v1(1, d, 1);
  SurrealS1 v2(1, 0);
  SurrealS1 v3(v1);
  SurrealS1 v4, v5;

  BOOST_CHECK( chkSurrealS1( v1,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  1, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3,  1, 2 ) );
//  BOOST_CHECK( chkSurrealS1( v4,  0, 0 ) );
//  BOOST_CHECK( chkSurrealS1( v5,  0, 0 ) );

  // assignment
  v3 = 2;
  BOOST_CHECK_EQUAL( 2, v3.value() );
  BOOST_CHECK_EQUAL( 0, v3.deriv() );

  v3 = v2;
  BOOST_CHECK_EQUAL( 1, v2.value() );
  BOOST_CHECK_EQUAL( 0, v2.deriv() );
  BOOST_CHECK_EQUAL( 1, v3.value() );
  BOOST_CHECK_EQUAL( 0, v3.deriv() );

  v3 = v2 = 2;
  BOOST_CHECK_EQUAL( 2, v2.value() );
  BOOST_CHECK_EQUAL( 0, v2.deriv() );
  BOOST_CHECK_EQUAL( 2, v3.value() );
  BOOST_CHECK_EQUAL( 0, v3.deriv() );

  v3 = v2 = v1;
  BOOST_CHECK_EQUAL( 1, v1.value() );
  BOOST_CHECK_EQUAL( 2, v1.deriv() );
  BOOST_CHECK_EQUAL( 1, v2.value() );
  BOOST_CHECK_EQUAL( 2, v2.deriv() );
  BOOST_CHECK_EQUAL( 1, v3.value() );
  BOOST_CHECK_EQUAL( 2, v3.deriv() );

  // unary
  v3 = +v2;
  BOOST_CHECK_EQUAL( 1, v2.value() );
  BOOST_CHECK_EQUAL( 2, v2.deriv() );
  BOOST_CHECK_EQUAL( 1, v3.value() );
  BOOST_CHECK_EQUAL( 2, v3.deriv() );

  v3 = -v2;
  BOOST_CHECK_EQUAL( 1, v2.value() );
  BOOST_CHECK_EQUAL( 2, v2.deriv() );
  BOOST_CHECK_EQUAL( -1, v3.value() );
  BOOST_CHECK_EQUAL( -2, v3.deriv() );

  // binary accumulation
  v4 = v3 = v2 = v1;
  v5 = 4*v1;
  v2 += 3;
  v3 -= 3;
  v4 *= 3;
  v5 /= 2;
  BOOST_CHECK( chkSurrealS1( v1,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  4, 2 ) );
  BOOST_CHECK( chkSurrealS1( v3, -2, 2 ) );
  BOOST_CHECK( chkSurrealS1( v4,  3, 6 ) );
  BOOST_CHECK( chkSurrealS1( v5,  2, 4 ) );

  v4 = v3 = v2 = v1;
  v5 = 4*v1;
  v2 += Real(3);
  v3 -= Real(3);
  v4 *= Real(3);
  v5 /= Real(2.);
  BOOST_CHECK( chkSurrealS1( v1,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  4, 2 ) );
  BOOST_CHECK( chkSurrealS1( v3, -2, 2 ) );
  BOOST_CHECK( chkSurrealS1( v4,  3, 6 ) );
  BOOST_CHECK( chkSurrealS1( v5,  2, 4 ) );

  v4 = v3 = v2 = v1;
  v5 = 4*v1;
  v2 += v1;
  v3 -= v1;
  v4 *= v1;
  v5 /= v1;
  BOOST_CHECK( chkSurrealS1( v1,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3,  0, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4,  1, 4 ) );
  BOOST_CHECK( chkSurrealS1( v5,  4, 0 ) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( addsubtract )
{
  SurrealS1 v1(1, 2);
  SurrealS1 v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealS1( v1,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v3,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v4,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v5,  1, 2 ) );

  // binary +/- operators

  SurrealS1 v8 = v1 + v2;
  BOOST_CHECK_EQUAL( 2, v8.value() );
  BOOST_CHECK_EQUAL( 4, v8.deriv() );
  BOOST_CHECK_EQUAL( 1, v8.size() );

  v2 = v1;
  v3 = v1 + v2;
  v4 = v1 + v2 + v3;
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v4, 4, 8 ) );

  v3 = v1 - v2;
  v4 = v1 - v2 - v3;
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v3, 0, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, 0, 0 ) );

  v3 = 3;
  v4 = v1 + v2 - v3;
  v5 = v1 - v2 + v3;
  BOOST_CHECK( chkSurrealS1( v1,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v3,  3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, -1, 4 ) );
  BOOST_CHECK( chkSurrealS1( v5,  3, 0 ) );

  v4 = v3 = v2 = v1;
  v2 += v1;
  v3 += v1 + v2;
  v4 += v1 + v2 + v3;
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 4, 8 ) );
  BOOST_CHECK( chkSurrealS1( v4, 8, 16 ) );

  v4 = v3 = v2 = v1;
  v2 -= v1;
  v3 -= v1 - v2;
  v4 -= v1 - v2 - v3;
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 0, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 0, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, 0, 0 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = (v1 + v2) + (v3 + v4);
  BOOST_CHECK( chkSurrealS1( v5, 10,20 ) );
  v5 = (v1 + v2) + (v3 - v4);
  BOOST_CHECK( chkSurrealS1( v5, 2,4 ) );
  v5 = (v1 + v2) - (v3 + v4);
  BOOST_CHECK( chkSurrealS1( v5, -4,-8 ) );
  v5 = (v1 + v2) - (v3 - v4);
  BOOST_CHECK( chkSurrealS1( v5, 4,8 ) );
  v5 = (v1 - v2) + (v3 + v4);
  BOOST_CHECK( chkSurrealS1( v5, 6,12 ) );
  v5 = (v1 - v2) + (v3 - v4);
  BOOST_CHECK( chkSurrealS1( v5, -2,-4 ) );
  v5 = (v1 - v2) - (v3 + v4);
  BOOST_CHECK( chkSurrealS1( v5, -8,-16 ) );
  v5 = (v1 - v2) - (v3 - v4);
  BOOST_CHECK( chkSurrealS1( v5, 0, 0 ) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, 6 ) );
  BOOST_CHECK( chkSurrealS1( v4, 4, 8 ) );

  v5 += (v1 + v2) + (v3 + v4);
  v5 += (v1 + v2) + (v3 - v4);
  v5 += (v1 + v2) - (v3 + v4);
  v5 += (v1 + v2) - (v3 - v4);
  v5 += (v1 - v2) + (v3 + v4);
  v5 += (v1 - v2) + (v3 - v4);
  v5 += (v1 - v2) - (v3 + v4);
  v5 += (v1 - v2) - (v3 - v4);
  BOOST_CHECK( chkSurrealS1( v5, 8, 16 ) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, 6 ) );
  BOOST_CHECK( chkSurrealS1( v4, 4, 8 ) );

  v5 -= (v1 + v2) + (v3 + v4);
  v5 -= (v1 + v2) + (v3 - v4);
  v5 -= (v1 + v2) - (v3 + v4);
  v5 -= (v1 + v2) - (v3 - v4);
  v5 -= (v1 - v2) + (v3 + v4);
  v5 -= (v1 - v2) + (v3 - v4);
  v5 -= (v1 - v2) - (v3 + v4);
  v5 -= (v1 - v2) - (v3 - v4);
  BOOST_CHECK( chkSurrealS1( v5, 0, 0 ) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, 6 ) );
  BOOST_CHECK( chkSurrealS1( v4, 4, 8 ) );

  v2 = +v1;
  BOOST_CHECK( chkSurrealS1( v1,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  1, 2 ) );
  v2 = -v1;
  BOOST_CHECK( chkSurrealS1( v1,  1,  2 ) );
  BOOST_CHECK( chkSurrealS1( v2, -1, -2 ) );
  v3 = +(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v3,  0, 0 ) );
  v3 = +(v1 - v2);
  BOOST_CHECK( chkSurrealS1( v3,  2, 4 ) );
  v3 = -(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v3,  0, 0 ) );
  v3 = -(v1 - v2);
  BOOST_CHECK( chkSurrealS1( v3, -2, -4 ) );

  // addition/subtraction with scalar quantities

  v3 = v1 + 3;
  BOOST_CHECK( chkSurrealS1( v3,  4, 2 ) );
  v3 = v1 + Real(3);
  BOOST_CHECK( chkSurrealS1( v3,  4, 2 ) );
  v3 = 3 + v1;
  BOOST_CHECK( chkSurrealS1( v3,  4, 2 ) );
  v3 = Real(3) + v1;
  BOOST_CHECK( chkSurrealS1( v3,  4, 2 ) );
  v3 += v1 + 1;
  BOOST_CHECK( chkSurrealS1( v3,  6, 4 ) );
  v3 -= v1 + 1;
  BOOST_CHECK( chkSurrealS1( v3,  4, 2 ) );
  v3 = v1 - 3;
  BOOST_CHECK( chkSurrealS1( v3, -2, 2 ) );
  v3 = v1 - Real(3);
  BOOST_CHECK( chkSurrealS1( v3, -2, 2 ) );
  v3 = 3 - v1;
  BOOST_CHECK( chkSurrealS1( v3,  2, -2 ) );
  v3 = Real(3) - v1;
  BOOST_CHECK( chkSurrealS1( v3,  2, -2 ) );

  v3 = +(v1 + v2) + 3;
  BOOST_CHECK( chkSurrealS1( v3,  3, 0 ) );
  v3 = -(v1 + v2) + 3;
  BOOST_CHECK( chkSurrealS1( v3,  3,  0 ) );
  BOOST_CHECK( chkSurrealS1( v1,  1,  2 ) );
  BOOST_CHECK( chkSurrealS1( v2, -1, -2 ) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( multiply )
{
  SurrealS1 v1(1, 2);
  SurrealS1 v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealS1( v1,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v3,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v4,  1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v5,  1, 2 ) );

  // binary * operators

  v2 = 3*v1;
  v3 = v2*2;
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 3, 6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 6, 12 ) );

  v2 += 3*v1;
  v3 += v2*2;
  BOOST_CHECK( chkSurrealS1( v1,  1,  2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  6, 12 ) );
  BOOST_CHECK( chkSurrealS1( v3, 18, 36 ) );

  v2 -= 3*v1;
  v3 -= v2*2;
  BOOST_CHECK( chkSurrealS1( v1,  1,  2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  3,  6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 12, 24 ) );

  v2 = 3/v1;
  v3 = v2/1;
  BOOST_CHECK( chkSurrealS1( v1, 1,  2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 3, -6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, -6 ) );

  v2 = 2*v1;
  v3 = v1*v2;
  BOOST_CHECK( chkSurrealS1( v3, 2, 8 ) );
  v3 += v1*v2;
  BOOST_CHECK( chkSurrealS1( v3, 4, 16 ) );

  v3 = 2*v1;
  v3 = v1*v3; //Test when v3 is both on left and right
  BOOST_CHECK( chkSurrealS1( v3, 2, 8 ) );
  v3 += v1*v3;
  BOOST_CHECK( chkSurrealS1( v3, 4, 20 ) );
  v3 -= v1*v3;
  BOOST_CHECK( chkSurrealS1( v3, 0, -8 ) );

  v2 = 2*(v1*2);
  BOOST_CHECK( chkSurrealS1( v2, 4, 8 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = 3*(v1 + v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 23, 46 ) );
  v5 = 3*(v1 + v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 7, 14 ) );
  v5 = 3*(v1 + v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -5, -10 ) );
  v5 = 3*(v1 + v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 11, 22 ) );
  v5 = 3*(v1 - v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 11, 22 ) );
  v5 = 3*(v1 - v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -5, -10 ) );
  v5 = 3*(v1 - v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -17, -34 ) );
  v5 = 3*(v1 - v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -1, -2 ) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, 6 ) );
  BOOST_CHECK( chkSurrealS1( v4, 4, 8 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = 5;
  v5 += 3*(v1 + v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 28, 46 ) );
  v5 += 3*(v1 + v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 35, 60 ) );
  v5 += 3*(v1 + v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 30, 50 ) );
  v5 += 3*(v1 + v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 41, 72 ) );
  v5 += 3*(v1 - v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 52, 94 ) );
  v5 += 3*(v1 - v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 47, 84 ) );
  v5 += 3*(v1 - v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 30, 50 ) );
  v5 += 3*(v1 - v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 29, 48 ) );
  BOOST_CHECK( chkSurrealS1( v1,  1,  2 ) );
  BOOST_CHECK( chkSurrealS1( v2,  2,  4 ) );
  BOOST_CHECK( chkSurrealS1( v3,  3,  6 ) );
  BOOST_CHECK( chkSurrealS1( v4,  4,  8 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = 5;
  v5 -= 3*(v1 + v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -18, -46 ) );
  v5 -= 3*(v1 + v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -25, -60 ) );
  v5 -= 3*(v1 + v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -20, -50 ) );
  v5 -= 3*(v1 + v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -31, -72 ) );
  v5 -= 3*(v1 - v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -42, -94 ) );
  v5 -= 3*(v1 - v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -37, -84 ) );
  v5 -= 3*(v1 - v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -20, -50 ) );
  v5 -= 3*(v1 - v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -19, -48 ) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, 6 ) );
  BOOST_CHECK( chkSurrealS1( v4, 4, 8 ) );

  v3 = 3*(v1 + v2)*2;
  BOOST_CHECK( chkSurrealS1( v3, 18, 36 ) );
  v3 = 3*2*(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v3, 18, 36 ) );
  v3 = (v1 + v2)*3*2;
  BOOST_CHECK( chkSurrealS1( v3, 18, 36 ) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4 ) );

  v2 = +3*v1;
  BOOST_CHECK( chkSurrealS1( v2,  3, 6 ) );
  v2 = -3*v1;
  BOOST_CHECK( chkSurrealS1( v2, -3, -6 ) );
  v2 = +v1*3;
  BOOST_CHECK( chkSurrealS1( v2,  3, 6 ) );
  v2 = -v1*3;
  BOOST_CHECK( chkSurrealS1( v2, -3, -6 ) );
  v2 = +(3*v1);
  BOOST_CHECK( chkSurrealS1( v2,  3, 6 ) );
  v2 = -(3*v1);
  BOOST_CHECK( chkSurrealS1( v2, -3, -6 ) );
  v2 = +(v1*3);
  BOOST_CHECK( chkSurrealS1( v2,  3, 6 ) );
  v2 = -(v1*3);
  BOOST_CHECK( chkSurrealS1( v2, -3, -6 ) );
  BOOST_CHECK( chkSurrealS1( v1,  1, 2 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = v1*v2*v3;
  BOOST_CHECK( chkSurrealS1( v4, 6, 36 ) );
  v4 = 2*v1*v2*v3;
  BOOST_CHECK( chkSurrealS1( v4, 12, 72 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = (v1 + v2)*v3;
  BOOST_CHECK( chkSurrealS1( v4, 12, 48 ) );
  v4 += (v1 + v2)*v3;
  BOOST_CHECK( chkSurrealS1( v4, 24, 96 ) );
  v4 = v3*(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v4, 12, 48 ) );
  v4 += v3*(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v4, 24, 96 ) );

  v2 = 2*v1;
  v3 = (v1 + v2)*(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v3, 9, 36 ) );
  v2 = 2*v1;
  v3 = 3*v1;
  v4 = (v1 + v2)*(v2 + v3);
  BOOST_CHECK( chkSurrealS1( v4, 15, 60 ) );
  v4 += (v1 + v2)*(v2 + v3);
  BOOST_CHECK( chkSurrealS1( v4, 30, 120 ) );
  v4 = 2*(v1 + v2)*(v2 + v3);
  BOOST_CHECK( chkSurrealS1( v4, 30, 120 ) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( divide )
{
  SurrealS1 v1(2, 4);
  SurrealS1 v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealS1( v1,  2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v2,  2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3,  2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v4,  2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v5,  2, 4 ) );

  // binary / operators

  v2 = 4/v1;
  v3 = v2/2;
  BOOST_CHECK( chkSurrealS1( v1, 2,  4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, -4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, -2 ) );

  v2 = Real(4)/v1;
  v3 = v2/Real(2.);
  BOOST_CHECK( chkSurrealS1( v1, 2,  4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, -4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, -2 ) );

  v2 += 4/v1;
  v3 += v2/2;
  BOOST_CHECK( chkSurrealS1( v1, 2,  4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, -8 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, -6 ) );

  v2 -= 4/v1;
  v3 -= v2/2;
  BOOST_CHECK( chkSurrealS1( v1, 2,  4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, -4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, -4 ) );

  v2 = 2*v1;
  v3 = v2/v1;
  BOOST_CHECK( chkSurrealS1( v1, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0 ) );

  v3 += v2/v1;
  BOOST_CHECK( chkSurrealS1( v1, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8 ) );
  BOOST_CHECK( chkSurrealS1( v3, 4, 0 ) );

  v3 -= v2/v1;
  BOOST_CHECK( chkSurrealS1( v1, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = (v2+v1)/v3;
  BOOST_CHECK( chkSurrealS1( v1, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8 ) );
  BOOST_CHECK( chkSurrealS1( v3, 6, 12 ) );
  BOOST_CHECK( chkSurrealS1( v4, 1, 0 ) );
  v4 += (v2+v1)/v3;
  BOOST_CHECK( chkSurrealS1( v4, 2, 0 ) );
  v4 = v3/(v2+v1);
  BOOST_CHECK( chkSurrealS1( v4, 1, 0 ) );
  v4 += v3/(v2+v1);
  BOOST_CHECK( chkSurrealS1( v4, 2, 0 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = 6*v1;
  v5 = 12/(v1 + v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 12, 16 ) );
  v5 = 12/(v1 + v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 0, -8 ) );
  v5 = 12/(v1 + v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -8, -24 ) );
  v5 = 12/(v1 + v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 4, 0 ) );
  v5 = 12/(v1 - v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 4, 32 ) );
  v5 = 12/(v1 - v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -8, 8 ) );
  v5 = 12/(v1 - v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -16, -8 ) );
  v5 = 12/(v1 - v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -4, 16 ) );
  BOOST_CHECK( chkSurrealS1( v1, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8 ) );
  BOOST_CHECK( chkSurrealS1( v3, 8, 16 ) );
  BOOST_CHECK( chkSurrealS1( v4, 12, 24 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = 6*v1;
  v5 = 1;
  v5 += 12/(v1 + v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 13, 16 ) );
  v5 += 12/(v1 + v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 13, 8 ) );
  v5 += 12/(v1 + v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 5, -16 ) );
  v5 += 12/(v1 + v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 9, -16 ) );
  v5 += 12/(v1 - v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 13, 16 ) );
  v5 += 12/(v1 - v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 5, 24 ) );
  v5 += 12/(v1 - v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -11, 16 ) );
  v5 += 12/(v1 - v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -15, 32 ) );
  BOOST_CHECK( chkSurrealS1( v1, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8 ) );
  BOOST_CHECK( chkSurrealS1( v3, 8, 16 ) );
  BOOST_CHECK( chkSurrealS1( v4, 12, 24 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = 6*v1;
  v5 = 1;
  v5 -= 12/(v1 + v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -11, -16 ) );
  v5 -= 12/(v1 + v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -11, -8 ) );
  v5 -= 12/(v1 + v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -3, 16 ) );
  v5 -= 12/(v1 + v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -7, 16 ) );
  v5 -= 12/(v1 - v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -11, -16 ) );
  v5 -= 12/(v1 - v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -3, -24 ) );
  v5 -= 12/(v1 - v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 13, -16 ) );
  v5 -= 12/(v1 - v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 17, -32 ) );
  BOOST_CHECK( chkSurrealS1( v1, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8 ) );
  BOOST_CHECK( chkSurrealS1( v3, 8, 16 ) );
  BOOST_CHECK( chkSurrealS1( v4, 12, 24 ) );

  v5 = 12/(v1 + v2)/2;
  BOOST_CHECK( chkSurrealS1( v5, 1, -2 ) );
  v5 = 12/2/(v1 - v2);
  BOOST_CHECK( chkSurrealS1( v5, -3, 6 ) );
  v5 = (v1 + v2)/3/2;
  BOOST_CHECK( chkSurrealS1( v5, 1, 2 ) );
  BOOST_CHECK( chkSurrealS1( v1, 2, 4 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8 ) );

  // cppcheck-suppress duplicateExpression
  v5 = (v1 + v2)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 1, 0 ) );
  v5 = 2*(v1 + v2)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 2, 0 ) );
  v5 += 2*(v1 + v2)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 4, 0 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v5 = (v2 + v3)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 2, 0 ) );
  v5 = 2*(v2 + v3)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 4, 0 ) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( relational )
{
  SurrealS1 v1(1, 3);
  SurrealS1 v2(1, 3);
  SurrealS1 v3(2, 3);

  BOOST_CHECK( chkSurrealS1( v1,  1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2,  1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v3,  2, 3 ) );

  BOOST_CHECK(   v2 == v1  );
  BOOST_CHECK( !(v3 == v1) );
  BOOST_CHECK(   v2 == 1   );
  BOOST_CHECK( !(v3 == 1)  );
  BOOST_CHECK(   1 == v2   );
  BOOST_CHECK( !(1 == v3)  );

  BOOST_CHECK( v2 == Real(1.) );
  BOOST_CHECK( Real(1.) == v2 );

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
  const double tol = 1.e-14;
  SurrealS1 v1(1, 3);
  SurrealS1 v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealS1( v1,  1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2,  1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v3,  1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v4,  1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v5,  1, 3 ) );

  // trig functions <cmath>

  v2 = cos(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, cos(1.), -3*sin(1.), tol ) );
  v2 += cos(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*cos(1.), -2*3*sin(1.), tol ) );
  v2 = sin(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, sin(1.), 3*cos(1.), tol ) );
  v2 += sin(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*sin(1.), 2*3*cos(1.), tol ) );
  v2 = tan(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, tan(1.), 3/(cos(1.)*cos(1.)), tol ) );
  v2 += tan(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*tan(1.), 2*3/(cos(1.)*cos(1.)), tol ) );

  v2 = acos(cos(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, tol ) );
  v2 += acos(cos(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6, tol ) );
  v2 = asin(sin(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, tol ) );
  v2 += asin(sin(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6, tol ) );
  v2 = atan(tan(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, tol ) );
  v2 += atan(tan(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6, tol ) );

  v2 = v1;
  v3 = v1;
  v4 = atan2(v2, v3);
  BOOST_CHECK( chkSurrealS1( v2, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v4, atan2(1., 1.), 0, tol ) );
  v4 += atan2(v2, v3);
  BOOST_CHECK( chkSurrealS1( v4, 2*atan2(1., 1.), 0, tol ) );

  // hyperbolic functions <cmath>

  v2 = cosh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3 ) );
  BOOST_CHECK( chkSurrealS1( v2, cosh(1.), 3*sinh(1.), tol ) );
  v2 += cosh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*cosh(1.), 2*3*sinh(1.), tol ) );
  v2 = sinh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3 ) );
  BOOST_CHECK( chkSurrealS1( v2, sinh(1.), 3*cosh(1.), tol ) );
  v2 += sinh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*sinh(1.), 2*3*cosh(1.), tol ) );
  v2 = tanh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3 ) );
  BOOST_CHECK( chkSurrealS1( v2, tanh(1.), 3/(cosh(1.)*cosh(1.)), tol ) );
  v2 += tanh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*tanh(1.), 2*3/(cosh(1.)*cosh(1.)), tol ) );

  // exp and log functions <cmath>

  v2 = exp(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, exp(1.), 3*exp(1.), tol ) );
  v2 += exp(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*exp(1.), 2*3*exp(1.), tol ) );

  v2 = expm1(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, expm1(1.), 3*exp(1.), tol ) );
  v2 += expm1(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*expm1(1.), 2*3*exp(1.), tol ) );

  v2 = log(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 0, 3, tol ) );
  v2 += log(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 0, 2*3, tol ) );

  v2 = log10(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 0, 3/log(10.), tol ) );
  v2 += log10(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 0, 2*3/log(10.), tol ) );

  v2 = log1p(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, log(2.), 3./2., tol ) );
  v2 += log1p(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*log(2.), 2*3/2., tol ) );

  // power functions <cmath>

  v2 = v1;
  v3 = pow(v2, v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, 3, tol ) );
  v3 += pow(v2, v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 6, tol ) );
  v2 = pow(v1, 2);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 6, tol ) );
  v2 = pow(v1, Real(2.));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 6, tol ) );
  v2 += pow(v1, 2);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 12, tol ) );
  v2 = pow(2, v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 2*log(8.), tol ) );
  v2 = pow(Real(2.), v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 2*log(8.), tol ) );
  v2 += pow(2, v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 4*log(8.), tol ) );
  v2 = v1;
  v3 = pow(v1+v2, v1+v2);
  BOOST_CHECK( chkSurrealS1( v3, 4, 40.63553233343869, tol ) );
  v3 += pow(v1+v2, v1+v2);
  BOOST_CHECK( chkSurrealS1( v3, 8, 2*40.63553233343869, tol ) );

  v2 = v1;
  v3 = pow(4*v1/v2, Real(0.5) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0, tol ) );
  v3 += pow(4*v1/v2, Real(0.5) );
  BOOST_CHECK( chkSurrealS1( v3, 4, 0, tol ) );

  v2 = pow(v1, 0);
  BOOST_CHECK( chkSurrealS1( v2, 1, 0 ) );
  v2 = pow(0, v1);
  BOOST_CHECK( chkSurrealS1( v2, 0, 0 ) );
  v2 = SurrealS1(0, 1);
  v3 = pow(v2, 0);
  BOOST_CHECK( chkSurrealS1( v2, 0, 1 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, 0 ) );
  v3 = pow(v2, 1);
  BOOST_CHECK( chkSurrealS1( v2, 0, 1 ) );
  BOOST_CHECK( chkSurrealS1( v3, 0, 1 ) );
  v3 = SurrealS1(0, 3);
  v4 = pow(v2, v3);
  BOOST_CHECK( chkSurrealS1( v2, 0, 1 ) );
  BOOST_CHECK( chkSurrealS1( v3, 0, 3 ) );
  BOOST_CHECK( chkSurrealS1( v4, 1, 0 ) );

  v2 = SurrealS1(1e-16, 1);
  v3 = pow(v2, 1);
  BOOST_CHECK( chkSurrealS1( v2, 1e-16, 1 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1e-16, 1 ) );

  #if 0
  v2 = SurrealS1(0, 1);
  std::cout << "pow: a = " << v2 << "  pow(a,1) = " << pow(v2, 1) << std::endl;
  v2 = SurrealS1(1e-16, 1);
  std::cout << "pow: a = " << v2 << "  pow(a,1) = " << pow(v2, 1) << std::endl;
  #endif

  v2 = sqrt(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 1.5, tol ) );
  v2 += sqrt(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 3, tol ) );
  v2 = 0;
  v3 = sqrt(v2);
  BOOST_CHECK( chkSurrealS1( v2, 0, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 0, 0, tol ) );
  v3 = sqrt(4*v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 3, tol ) );
  v3 += sqrt(4*v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v3, 4, 6, tol ) );

  // rounding functions <cmath>

  v2 = 1.5*v1;
  v3 = ceil(v2);
  BOOST_CHECK( chkSurrealS1( v2, 1.5, 4.5 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0, tol ) );
  v3 = floor(v2);
  BOOST_CHECK( chkSurrealS1( v2, 1.5, 4.5 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1,0, tol ) );

  // misc functions <cmath>

  v2 = abs(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3 ) );
  v2 = fabs(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3 ) );

  v2 = 2*v1;
  v3 = max(v1, v2);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 6 ) );
  v3 = min(v1, v2);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, 3 ) );

  v3 = max(v1, 2.);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0 ) );
  v3 = min(v1, 2.);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, 3 ) );

  v3 = max(2., v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0 ) );
  v3 = min(2., v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, 3 ) );

  v2 = -2*v1;
  v3 = max( fabs(v1), fabs(v2) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, -2, -6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 6 ) );
  v3 = min( fabs(v1), fabs(v2) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 3 ) );
  BOOST_CHECK( chkSurrealS1( v2, -2, -6 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, 3 ) );
}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( multiply_deriv2 )
{
  SurrealS<1, SurrealS1> v1(1, 2);
  SurrealS<1, SurrealS1> v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealS1( v1,  1, 2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2,  1, 2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3,  1, 2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4,  1, 2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v5,  1, 2, 0 ) );

  // binary * operators

  v1.value().deriv() = 1;

  v2 = 3*v1;
  v3 = v2*2;
  BOOST_CHECK( chkSurrealS1( v1, 1, 2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 3, 6, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 6, 12, 0 ) );

  v2 += 3*v1;
  v3 += v2*2;
  BOOST_CHECK( chkSurrealS1( v1,  1,  2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2,  6, 12, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 18, 36, 0 ) );

  v2 -= 3*v1;
  v3 -= v2*2;
  BOOST_CHECK( chkSurrealS1( v1,  1,  2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2,  3,  6, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 12, 24, 0 ) );

  v2 = 3/v1;
  v3 = v2/1;
  BOOST_CHECK( chkSurrealS1( v1, 1,  2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 3, -6, 12 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, -6, 12 ) );

  v2 = 2*v1;
  v3 = v1*v2;
  BOOST_CHECK( chkSurrealS1( v3, 2, 8, 8 ) );
  v3 += v1*v2;
  BOOST_CHECK( chkSurrealS1( v3, 4, 16, 16 ) );

  v3 = 2*v1;
  v3 = v1*v3; //Test when v3 is both on left and right
  BOOST_CHECK( chkSurrealS1( v3, 2, 8, 8 ) );
  v3 += v1*v3;
  BOOST_CHECK( chkSurrealS1( v3, 4, 20, 32 ) );
  v3 -= v1*v3;
  BOOST_CHECK( chkSurrealS1( v3, 0, -8, -40 ) );

  v2 = 2*(v1*2);
  BOOST_CHECK( chkSurrealS1( v2, 4, 8, 0 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = 3*(v1 + v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 23, 46, 0 ) );
  v5 = 3*(v1 + v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 7, 14, 0 ) );
  v5 = 3*(v1 + v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -5, -10, 0 ) );
  v5 = 3*(v1 + v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 11, 22, 0 ) );
  v5 = 3*(v1 - v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 11, 22, 0 ) );
  v5 = 3*(v1 - v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -5, -10, 0 ) );
  v5 = 3*(v1 - v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -17, -34, 0 ) );
  v5 = 3*(v1 - v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -1, -2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, 6, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, 4, 8, 0 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = 5;
  v5 += 3*(v1 + v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 28, 46, 0 ) );
  v5 += 3*(v1 + v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 35, 60, 0 ) );
  v5 += 3*(v1 + v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 30, 50, 0 ) );
  v5 += 3*(v1 + v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 41, 72, 0 ) );
  v5 += 3*(v1 - v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 52, 94, 0 ) );
  v5 += 3*(v1 - v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 47, 84, 0 ) );
  v5 += 3*(v1 - v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 30, 50, 0 ) );
  v5 += 3*(v1 - v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, 29, 48, 0 ) );
  BOOST_CHECK( chkSurrealS1( v1,  1,  2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2,  2,  4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3,  3,  6, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4,  4,  8, 0 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = 4*v1;
  v5 = 5;
  v5 -= 3*(v1 + v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -18, -46, 0 ) );
  v5 -= 3*(v1 + v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -25, -60, 0 ) );
  v5 -= 3*(v1 + v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -20, -50, 0 ) );
  v5 -= 3*(v1 + v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -31, -72, 0 ) );
  v5 -= 3*(v1 - v2) + (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -42, -94, 0 ) );
  v5 -= 3*(v1 - v2) + (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -37, -84, 0 ) );
  v5 -= 3*(v1 - v2) - (v3 + v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -20, -50, 0 ) );
  v5 -= 3*(v1 - v2) - (v3 - v4)*2;
  BOOST_CHECK( chkSurrealS1( v5, -19, -48, 0 ) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, 6, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, 4, 8, 0 ) );

  v3 = 3*(v1 + v2)*2;
  BOOST_CHECK( chkSurrealS1( v3, 18, 36, 0 ) );
  v3 = 3*2*(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v3, 18, 36, 0 ) );
  v3 = (v1 + v2)*3*2;
  BOOST_CHECK( chkSurrealS1( v3, 18, 36, 0 ) );
  BOOST_CHECK( chkSurrealS1( v1, 1, 2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 4, 0 ) );

  v2 = +3*v1;
  BOOST_CHECK( chkSurrealS1( v2,  3, 6, 0 ) );
  v2 = -3*v1;
  BOOST_CHECK( chkSurrealS1( v2, -3, -6, 0 ) );
  v2 = +v1*3;
  BOOST_CHECK( chkSurrealS1( v2,  3, 6, 0 ) );
  v2 = -v1*3;
  BOOST_CHECK( chkSurrealS1( v2, -3, -6, 0 ) );
  v2 = +(3*v1);
  BOOST_CHECK( chkSurrealS1( v2,  3, 6, 0 ) );
  v2 = -(3*v1);
  BOOST_CHECK( chkSurrealS1( v2, -3, -6, 0 ) );
  v2 = +(v1*3);
  BOOST_CHECK( chkSurrealS1( v2,  3, 6, 0 ) );
  v2 = -(v1*3);
  BOOST_CHECK( chkSurrealS1( v2, -3, -6, 0 ) );
  BOOST_CHECK( chkSurrealS1( v1,  1, 2, 0 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = v1*v2*v3;
  BOOST_CHECK( chkSurrealS1( v4, 6, 36, 72 ) );
  v4 = 2*v1*v2*v3;
  BOOST_CHECK( chkSurrealS1( v4, 12, 72, 144 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = (v1 + v2)*v3;
  BOOST_CHECK( chkSurrealS1( v4, 12, 48, 48 ) );
  v4 += (v1 + v2)*v3;
  BOOST_CHECK( chkSurrealS1( v4, 24, 96, 96 ) );
  v4 = v3*(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v4, 12, 48, 48 ) );
  v4 += v3*(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v4, 24, 96, 96 ) );

  v2 = 2*v1;
  v3 = (v1 + v2)*(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v3, 9, 36, 36 ) );
  v2 = 2*v1;
  v3 = 3*v1;
  v4 = (v1 + v2)*(v2 + v3);
  BOOST_CHECK( chkSurrealS1( v4, 15, 60, 60 ) );
  v4 += (v1 + v2)*(v2 + v3);
  BOOST_CHECK( chkSurrealS1( v4, 30, 120, 120 ) );
  v4 = 2*(v1 + v2)*(v2 + v3);
  BOOST_CHECK( chkSurrealS1( v4, 30, 120, 120 ) );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( divide_deriv2 )
{
  SurrealS<1, SurrealS1> v1(2, 4);
  SurrealS<1, SurrealS1> v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealS1( v1,  2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2,  2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3,  2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4,  2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v5,  2, 4, 0 ) );

  v1.value().deriv() = 1;

  // binary / operators

  v2 = 4/v1;
  v3 = v2/2;
  BOOST_CHECK( chkSurrealS1( v1, 2,  4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, -4, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, -2, 2 ) );

  v2 = Real(4)/v1;
  v3 = v2/Real(2.);
  BOOST_CHECK( chkSurrealS1( v1, 2,  4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, -4, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, -2, 2 ) );

  v2 += 4/v1;
  v3 += v2/2;
  BOOST_CHECK( chkSurrealS1( v1, 2,  4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, -8, 8 ) );
  BOOST_CHECK( chkSurrealS1( v3, 3, -6, 6 ) );

  v2 -= 4/v1;
  v3 -= v2/2;
  BOOST_CHECK( chkSurrealS1( v1, 2,  4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, -4, 4 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, -4, 4 ) );

  v2 = 2*v1;
  v3 = v2/v1;
  BOOST_CHECK( chkSurrealS1( v1, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0, 0 ) );

  v3 += v2/v1;
  BOOST_CHECK( chkSurrealS1( v1, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 4, 0, 0 ) );

  v3 -= v2/v1;
  BOOST_CHECK( chkSurrealS1( v1, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0, 0 ) );

  v2 = 2*v1;
  v3 = 3*v1;
  v4 = (v2+v1)/v3;
  BOOST_CHECK( chkSurrealS1( v1, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 6, 12, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, 1, 0, 0 ) );
  v4 += (v2+v1)/v3;
  BOOST_CHECK( chkSurrealS1( v4, 2, 0, 0 ) );
  v4 = v3/(v2+v1);
  BOOST_CHECK( chkSurrealS1( v4, 1, 0, 0 ) );
  v4 += v3/(v2+v1);
  BOOST_CHECK( chkSurrealS1( v4, 2, 0, 0 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = 6*v1;
  v5 = 12/(v1 + v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 12, 16, 4 ) );
  v5 = 12/(v1 + v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 0, -8, 4 ) );
  v5 = 12/(v1 + v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -8, -24, 4 ) );
  v5 = 12/(v1 + v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 4, 0, 4 ) );
  v5 = 12/(v1 - v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 4, 32, -12 ) );
  v5 = 12/(v1 - v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -8, 8, -12 ) );
  v5 = 12/(v1 - v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -16, -8, -12 ) );
  v5 = 12/(v1 - v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -4, 16, -12 ) );
  BOOST_CHECK( chkSurrealS1( v1, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 8, 16, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, 12, 24, 0 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = 6*v1;
  v5 = 1;
  v5 += 12/(v1 + v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 13, 16, 4 ) );
  v5 += 12/(v1 + v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 13, 8, 8 ) );
  v5 += 12/(v1 + v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 5, -16, 12 ) );
  v5 += 12/(v1 + v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 9, -16, 16 ) );
  v5 += 12/(v1 - v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 13, 16, 4 ) );
  v5 += 12/(v1 - v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 5, 24, -8 ) );
  v5 += 12/(v1 - v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -11, 16, -20 ) );
  v5 += 12/(v1 - v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -15, 32, -32 ) );
  BOOST_CHECK( chkSurrealS1( v1, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 8, 16, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, 12, 24, 0 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v4 = 6*v1;
  v5 = 1;
  v5 -= 12/(v1 + v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -11, -16, -4 ) );
  v5 -= 12/(v1 + v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -11, -8, -8 ) );
  v5 -= 12/(v1 + v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -3, 16, -12 ) );
  v5 -= 12/(v1 + v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -7, 16, -16 ) );
  v5 -= 12/(v1 - v2) + (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -11, -16, -4 ) );
  v5 -= 12/(v1 - v2) + (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, -3, -24, 8 ) );
  v5 -= 12/(v1 - v2) - (v3 + v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 13, -16, 20 ) );
  v5 -= 12/(v1 - v2) - (v3 - v4)/2;
  BOOST_CHECK( chkSurrealS1( v5, 17, -32, 32 ) );
  BOOST_CHECK( chkSurrealS1( v1, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 8, 16, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, 12, 24, 0 ) );

  v5 = 12/(v1 + v2)/2;
  BOOST_CHECK( chkSurrealS1( v5, 1, -2, 2 ) );
  v5 = 12/2/(v1 - v2);
  BOOST_CHECK( chkSurrealS1( v5, -3, 6, -6 ) );
  v5 = (v1 + v2)/3/2;
  BOOST_CHECK( chkSurrealS1( v5, 1, 2, 0 ) );
  BOOST_CHECK( chkSurrealS1( v1, 2, 4, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 8, 0 ) );

  // cppcheck-suppress duplicateExpression
  v5 = (v1 + v2)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 1, 0, 0 ) );
  v5 = 2*(v1 + v2)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 2, 0, 0 ) );
  v5 += 2*(v1 + v2)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 4, 0, 0 ) );

  v2 = 2*v1;
  v3 = 4*v1;
  v5 = (v2 + v3)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 2, 0, 0 ) );
  v5 = 2*(v2 + v3)/(v1 + v2);
  BOOST_CHECK( chkSurrealS1( v5, 4, 0, 0 ) );
}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( cmath_deriv2 )
{
  const double tol = 1.e-13;
  SurrealS<1, SurrealS1> v1(1, 3);
  SurrealS<1, SurrealS1> v2(v1), v3(v1), v4(v1), v5(v1);

  BOOST_CHECK( chkSurrealS1( v1,  1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2,  1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3,  1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4,  1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v5,  1, 3, 0 ) );

  v1.value().deriv() = 1;

  // trig functions <cmath>

  v2 = cos(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, cos(1.), -3*sin(1.), -3*cos(1.), tol ) );
  v2 += cos(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*cos(1.), -2*3*sin(1.), -2*3*cos(1.), tol ) );
  v2 = sin(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, sin(1.), 3*cos(1.), -3*sin(1.), tol ) );
  v2 += sin(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*sin(1.), 2*3*cos(1.), -2*3*sin(1.), tol ) );
  v2 = tan(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, tan(1.), 3/(cos(1.)*cos(1.)), 6/(cos(1.)*cos(1.))*tan(1.), tol ) );
  v2 += tan(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*tan(1.), 2*3/(cos(1.)*cos(1.)), 2*6/(cos(1.)*cos(1.))*tan(1.), tol ) );

  v2 = acos(cos(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, 0, tol ) );
  v2 += acos(cos(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6, 0, tol ) );
  v2 = asin(sin(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, 0, tol ) );
  v2 += asin(sin(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6, 0, tol ) );
  v2 = atan(tan(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, 0, tol ) );
  v2 += atan(tan(v1));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6, 0, tol ) );

  v2 = v1;
  v3 = v1;
  v4 = atan2(v2, v3);
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v4, atan2(1., 1.), 0, 0, tol ) );
  v4 += atan2(v2, v3);
  BOOST_CHECK( chkSurrealS1( v4, 2*atan2(1., 1.), 0, 0, tol ) );

  // hyperbolic functions <cmath>

  v2 = cosh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, cosh(1.), 3*sinh(1.), 3*cosh(1.), tol ) );
  v2 += cosh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*cosh(1.), 2*3*sinh(1.), 2*3*cosh(1.), tol ) );
  v2 = sinh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, sinh(1.), 3*cosh(1.), 3*sinh(1.), tol ) );
  v2 += sinh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*sinh(1.), 2*3*cosh(1.), 2*3*sinh(1.), tol ) );
  v2 = tanh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, tanh(1.), 3/(cosh(1.)*cosh(1.)), -6/(cosh(1.)*cosh(1.))*tanh(1.), tol ) );
  v2 += tanh(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1,3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*tanh(1.), 2*3/(cosh(1.)*cosh(1.)), -2*6/(cosh(1.)*cosh(1.))*tanh(1.), tol ) );

  // exp and log functions <cmath>

  v2 = exp(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, exp(1.), 3*exp(1.), 3*exp(1.), tol ) );
  v2 += exp(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*exp(1.), 2*3*exp(1.), 2*3*exp(1.), tol ) );

  v2 = log(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 0, 3, -3, tol ) );
  v2 += log(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 0, 2*3, -2*3, tol ) );

  v2 = log10(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 0, 3/log(10.), -3/log(10.), tol ) );
  v2 += log10(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 0, 2*3/log(10.), -2*3/log(10.), tol ) );

  v2 = log1p(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, log(2.), 3./2., -3./4., tol ) );
  v2 += log1p(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2*log(2.), 2*3/2., -2*3./4., tol ) );

  // power functions <cmath>

  v2 = v1;
  v3 = pow(v2, v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1, 3, 6, tol ) );
  v3 += pow(v2, v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 6, 12, tol ) );
  v2 = pow(v1, 2);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 6, 6, tol ) );
  v2 = pow(v1, Real(2.));
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 6, 6, tol ) );
  v2 += pow(v1, 2);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 12, 12, tol ) );
  v2 = pow(2, v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6*log(2.), 6*log(2.)*log(2.), tol ) );
  v2 = pow(Real(2.), v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 6*log(2.), 6*log(2.)*log(2.), tol ) );
  v2 += pow(2, v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 4, 2*6*log(2.), 2*6*log(2.)*log(2.), tol ) );
  v2 = v1;
  v3 = pow(v1+v2, v1+v2);
  BOOST_CHECK( chkSurrealS1( v3, 4, 3*(8 + 8*log(2.)), 3*(24 + 32*log(2.) + 16*log(2.)*log(2.)), tol ) );
  v3 += pow(v1+v2, v1+v2);
  BOOST_CHECK( chkSurrealS1( v3, 8, 2*3*(8 + 8*log(2.)), 2*3*(24 + 32*log(2.) + 16*log(2.)*log(2.)), tol ) );

  v2 = v1;
  v3 = pow(4*v1/v2, Real(0.5) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0, 0, tol ) );
  v3 += pow(4*v1/v2, Real(0.5) );
  BOOST_CHECK( chkSurrealS1( v3, 4, 0, 0, tol ) );


  v2 = sqrt(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3./2., -3./4., tol ) );
  v2 += sqrt(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 2, 3, -3./2., tol ) );
  v2 = 0;
  v3 = sqrt(v2);
  BOOST_CHECK( chkSurrealS1( v2, 0, 0, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 0, 0, 0, tol ) );
  v3 = sqrt(4*v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 3, -3./2., tol ) );
  v3 += sqrt(4*v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 4, 6, -3, tol ) );

  // rounding functions <cmath>

  v2 = 1.5*v1;
  v3 = ceil(v2);
  BOOST_CHECK( chkSurrealS1( v2, 1.5, 4.5, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 2, 0, 0, tol ) );
  v3 = floor(v2);
  BOOST_CHECK( chkSurrealS1( v2, 1.5, 4.5, 0 ) );
  BOOST_CHECK( chkSurrealS1( v3, 1,0, 0, tol ) );

  // misc functions <cmath>

  v2 = abs(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, 0 ) );
  v2 = fabs(v1);
  BOOST_CHECK( chkSurrealS1( v1, 1, 3, 0 ) );
  BOOST_CHECK( chkSurrealS1( v2, 1, 3, 0 ) );
}

//----------------------------------------------------------------------------//
/*
BOOST_AUTO_TEST_CASE( IO )
{
  //Set the 2nd argument to false to regenerate the pattern file
  output_test_stream output( "IO/Surreal/SurrealS1_pattern.txt", true );

  SurrealS1 v1(1, 3);
  SurrealS1 v2(v1);

  output << v1 << std::endl;
  BOOST_CHECK( output.match_pattern() );
  //v1.dump( 2, output );
  //BOOST_CHECK( output.match_pattern() );

  output << (v1 + v2) << std::endl;
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
  multiply_deriv2();
  divide_deriv2();
  cmath_deriv2();
//  IO();
  std::cout << std::endl;
  std::cout << "SurrealS1_test_suite Complete!" << std::endl;
  std::cout << std::endl;
  return 0;
}
