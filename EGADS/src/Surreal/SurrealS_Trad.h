// Modified from Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2020, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALS_TRAD_H
#define SURREALS_TRAD_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <iostream>
#include <cassert>

#include <string>

#include "always_inline.h"

//----------------------------------------------------------------------------//
// SurrealS<N>:  value, N derivatives
//
// Traditional implementation of Operators
//
// statically defined derivative array
//----------------------------------------------------------------------------//

template<int N_, class T = double>
class SurrealS
{
public:
  static const int N = N_;

  //The default constructor is intentionally empty here. This means Surreal is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  ALWAYS_INLINE SurrealS() {}
  ALWAYS_INLINE SurrealS( const SurrealS& z );
  ALWAYS_INLINE SurrealS( const int v0 );
  ALWAYS_INLINE SurrealS( const double v0 );
  ALWAYS_INLINE SurrealS( const double& v0, const double d0[], const int n );
  ALWAYS_INLINE SurrealS( const double& v0, const double& d0 );

  ALWAYS_INLINE ~SurrealS() {}

  ALWAYS_INLINE int size() const { return N; }

  // value accessor operators
  ALWAYS_INLINE T& value()       { return v_; }
  ALWAYS_INLINE T  value() const { return v_; }

  // derivative accessor operators
  ALWAYS_INLINE T& deriv( int i=0 )       { return d_[i]; }
  ALWAYS_INLINE T  deriv( int i=0 ) const { return d_[i]; }

  // assignment
  SurrealS& operator=( const SurrealS& );
  SurrealS& operator=( const double& );

  // unary operators; no side effects
  const SurrealS& operator+() const;
  const SurrealS  operator-() const;

  // binary accumulation operators
  SurrealS& operator+=( const SurrealS& );
  SurrealS& operator+=( const double& );
  SurrealS& operator-=( const SurrealS& );
  SurrealS& operator-=( const double& );
  SurrealS& operator*=( const SurrealS& );
  SurrealS& operator*=( const double& );
  SurrealS& operator/=( const SurrealS& );
  SurrealS& operator/=( const double& );

  // binary operators
  template<int M, class R> friend SurrealS<M,R> operator+( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> operator+( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend SurrealS<M,R> operator+( const double&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> operator-( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> operator-( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend SurrealS<M,R> operator-( const double&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> operator*( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> operator*( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend SurrealS<M,R> operator*( const double&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> operator/( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> operator/( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend SurrealS<M,R> operator/( const double&, const SurrealS<M,R>& );

  // relational operators
  template<int M, class R> friend bool operator==( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator==( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend bool operator==( const double&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator!=( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator!=( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend bool operator!=( const double&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator>( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator>( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend bool operator>( const double&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator<( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator<( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend bool operator<( const double&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator>=( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator>=( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend bool operator>=( const double&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator<=( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend bool operator<=( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend bool operator<=( const double&, const SurrealS<M,R>& );

  // trig functions <cmath>
  template<int M, class R> friend SurrealS<M,R> cos( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> sin( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> tan( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> acos( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> asin( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> atan( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> atan2( const SurrealS<M,R>&, const SurrealS<M,R>& );

  // hyperbolic functions <cmath>
  template<int M, class R> friend SurrealS<M,R> cosh( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> sinh( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> tanh( const SurrealS<M,R>& );

  // exp and log functions <cmath>
  template<int M, class R> friend SurrealS<M,R> exp( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> expm1( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> log( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> log10( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> log1p( const SurrealS<M,R>& );

  // power functions <cmath>
  template<int M, class R> friend SurrealS<M,R> pow( const SurrealS<M,R>&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> pow( const SurrealS<M,R>&, const double& );
  template<int M, class R> friend SurrealS<M,R> pow( const SurrealS<M,R>&, const int& );
  template<int M, class R> friend SurrealS<M,R> pow( const double&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> pow( const int&, const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> sqrt( const SurrealS<M,R>& );

  // rounding functions <cmath>
  template<int M, class R> friend SurrealS<M,R> ceil( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> floor( const SurrealS<M,R>& );

  // misc functions <cmath>
  template<int M, class R> friend SurrealS<M,R> abs( const SurrealS<M,R>& );
  template<int M, class R> friend SurrealS<M,R> fabs( const SurrealS<M,R>& );

  // classification functions <cmath>
  template<int M, class R> friend bool isfinite( const SurrealS<M,R>& );
  template<int M, class R> friend bool isinf( const SurrealS<M,R>& );
  template<int M, class R> friend bool isnan( const SurrealS<M,R>& );

  // input/output
  template<int M, class R> friend std::istream& operator>>( std::istream&, SurrealS<M>& );
  template<int M, class R> friend std::ostream& operator<<( std::ostream&, const SurrealS<M>& );

  void dump( int indentSize=0 ) const;

private:
  T v_;      // value
  T d_[N];   // derivative array
};


// constructors

template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const SurrealS& z ) : v_(z.v_)
{
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
}
template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const int v0 )
{
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = 0;
}
template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const double v0 )
{
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = 0;
}
template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const double& v0, const double d0[], const int n ) : v_(v0)
{
  assert( n == N );
  for (int i = 0; i < N; i++)
    d_[i] = d0[i];
}

template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const double& v0, const double& d0 ) : v_(v0)
{
  for (int i = 0; i < N; i++)
    d_[i] = d0;
}


// assignment

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator=( const SurrealS& z )
{
  //Do nothing if assigning self to self
  if ( &z == this ) return *this;

  v_ = z.v_;
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator=( const double& r )
{
  v_ = r;
  for (int i = 0; i < N; i++)
    d_[i] = 0;

  return *this;
}


// unary operators; no side effects

template<int N, class T>
ALWAYS_INLINE const SurrealS<N,T>&
SurrealS<N,T>::operator+() const
{
  return *this;
}

template<int N, class T>
ALWAYS_INLINE const SurrealS<N,T>
SurrealS<N,T>::operator-() const
{
  SurrealS<N,T> c;
  c.v_ = -v_;
  for (int i = 0; i < N; i++)
    c.d_[i] = -d_[i];

  return c;
}


// binary accumulation operators

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator+=( const SurrealS& z )
{
  v_ += z.v_;
  for (int i = 0; i < N; i++)
    d_[i] += z.d_[i];
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator+=( const double& r )
{
  v_ += r;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator-=( const SurrealS& z )
{
  v_ -= z.v_;
  for (int i = 0; i < N; i++)
    d_[i] -= z.d_[i];
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator-=( const double& r )
{
  v_ -= r;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator*=( const SurrealS& z )
{
  for (int i = 0; i < N; i++)
    d_[i] = v_*z.d_[i] + d_[i]*z.v_;
  v_ *= z.v_;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator*=( const double& r )
{
  for (int i = 0; i < N; i++)
    d_[i] *= r;
  v_ *= r;
  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator/=( const SurrealS& z)
{
  double tmp = 1./(z.v_*z.v_);
  for (int i = 0; i < N; i++)
    d_[i] = (z.v_*d_[i] - v_*z.d_[i])*tmp;
  v_ /= z.v_;

  return *this;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator/=( const double& r )
{
  double tmp = 1./r;
  for (int i = 0; i < N; i++)
    d_[i] *= tmp;
  v_ *= tmp;
  return *this;
}


// debug dump of private data
template<int N, class T>
void
SurrealS<N,T>::dump( int indentSize ) const
{
  std::string indent(indentSize, ' ');
  std::cout << indent << "SurrealS<" << N << ">: v_ = " << v_;
  std::cout << "  d_[" << N << "] = (";
  for (int n = 0; n < N-1; n++)
    std::cout << d_[n] << ",";
  std::cout << d_[N-1] << ")" << std::endl;
}


// binary operators

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator+( const SurrealS<N,T>& a, const SurrealS<N,T>& b )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i] + b.d_[i];
  c.v_ = a.v_ + b.v_;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator+( const SurrealS<N,T>& a, const double& b )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i];
  c.v_ = a.v_ + b;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator+( const double& a, const SurrealS<N,T>& b )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = b.d_[i];
  c.v_ = a + b.v_;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator-( const SurrealS<N,T>& a, const SurrealS<N,T>& b )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i] - b.d_[i];
  c.v_ = a.v_ - b.v_;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator-( const SurrealS<N,T>& a, const double& b )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i];
  c.v_ = a.v_ - b;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator-( const double& a, const SurrealS<N,T>& b )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = -b.d_[i];
  c.v_ = a - b.v_;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator*( const SurrealS<N,T>& a, const SurrealS<N,T>& b )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.v_*b.d_[i] + a.d_[i]*b.v_;
  c.v_ = a.v_*b.v_;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator*( const SurrealS<N,T>& a, const double& b )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*b;
  c.v_ = a.v_*b;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator*( const double& a, const SurrealS<N,T>& b )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a*b.d_[i];
  c.v_ = a*b.v_;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator/( const SurrealS<N,T>& a, const SurrealS<N,T>& b )
{
  SurrealS<N,T> c;
  T tmp = 1./(b.v_*b.v_);
  for (int i = 0; i < N; i++)
    c.d_[i] = (b.v_*a.d_[i] - a.v_*b.d_[i])*tmp;
  c.v_ = a.v_/b.v_;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T> operator/( const SurrealS<N,T>& a, const double& b )
{
  SurrealS<N,T> c;
  double tmp = 1./b;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*tmp;
  c.v_ = a.v_/b;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
operator/( const double& a, const SurrealS<N,T>& b )
{
  T tmpv = a/(b.v_);
  T tmpd = -1./(b.v_);
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmpv*tmpd*b.d_[i];
  c.v_ = tmpv;
  return c;
}


// relational operators

template<int N, class T>
ALWAYS_INLINE bool
operator==( const SurrealS<N,T>& lhs, const SurrealS<N,T>& rhs )
{
  return lhs.v_ == rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator==( const SurrealS<N,T>& lhs, const double& rhs )
{
  return lhs.v_ == rhs;
}

template<int N, class T>
ALWAYS_INLINE bool
operator==( const double& lhs, const SurrealS<N,T>& rhs )
{
  return lhs == rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator!=( const SurrealS<N,T>& lhs, const SurrealS<N,T>& rhs )
{
  return lhs.v_ != rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator!=( const SurrealS<N,T>& lhs, const double& rhs )
{
  return lhs.v_ != rhs;
}

template<int N, class T>
ALWAYS_INLINE bool
operator!=( const double& lhs, const SurrealS<N,T>& rhs )
{
  return lhs != rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator>( const SurrealS<N,T>& lhs, const SurrealS<N,T>& rhs )
{
  return lhs.v_ > rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator>( const SurrealS<N,T>& lhs, const double& rhs )
{
  return lhs.v_ > rhs;
}

template<int N, class T>
ALWAYS_INLINE bool
operator>( const double& lhs, const SurrealS<N,T>& rhs )
{
  return lhs > rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator<( const SurrealS<N,T>& lhs, const SurrealS<N,T>& rhs )
{
  return lhs.v_ < rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator<( const SurrealS<N,T>& lhs, const double& rhs )
{
  return lhs.v_ < rhs;
}

template<int N, class T>
ALWAYS_INLINE bool
operator<( const double& lhs, const SurrealS<N,T>& rhs )
{
  return lhs < rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator>=( const SurrealS<N,T>& lhs, const SurrealS<N,T>& rhs )
{
  return lhs.v_ >= rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator>=( const SurrealS<N,T>& lhs, const double& rhs )
{
  return lhs.v_ >= rhs;
}

template<int N, class T>
ALWAYS_INLINE bool
operator>=( const double& lhs, const SurrealS<N,T>& rhs )
{
  return lhs >= rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator<=( const SurrealS<N,T>& lhs, const SurrealS<N,T>& rhs )
{
  return lhs.v_ <= rhs.v_;
}

template<int N, class T>
ALWAYS_INLINE bool
operator<=( const SurrealS<N,T>& lhs, const double& rhs )
{
  return lhs.v_ <= rhs;
}

template<int N, class T>
ALWAYS_INLINE bool
operator<=( const double& lhs, const SurrealS<N,T>& rhs )
{
  return lhs <= rhs.v_;
}

//Macros for functions

#define SURREALS_FUNC1( NAME, FUNC, DERIV ) \
template<int N, class T> \
ALWAYS_INLINE SurrealS<N,T> \
NAME( const SurrealS<N,T>& z ) \
{ \
  T tmp = DERIV; \
  SurrealS<N,T> c; \
  for (int i = 0; i < N; i++) \
    c.d_[i] = tmp*z.d_[i]; \
  c.v_ = FUNC; \
  return c; \
}

#define SURREALS_FUNC2( NAME, FUNC, DERIV ) \
template<int N, class T> \
ALWAYS_INLINE SurrealS<N,T> \
NAME( const SurrealS<N,T>& z1, const SurrealS<N,T>& z2) \
{ \
  T tmp = DERIV; \
  SurrealS<N,T> c; \
  for (int i = 0; i < N; i++) \
    c.d_[i] = tmp*(z2.v_*z1.d_[i] - z1.v_*z2.d_[i]); \
  c.v_ = FUNC; \
  return c; \
}

// trig functions <cmath>

SURREALS_FUNC1( cos, cos(z.v_), -sin(z.v_) )
SURREALS_FUNC1( sin, sin(z.v_),  cos(z.v_) )
SURREALS_FUNC1( tan, tan(z.v_),  double(1)/(cos(z.v_)*cos(z.v_)) )
SURREALS_FUNC1( acos, acos(z.v_), -double(1)/sqrt(1 - z.v_*z.v_) )
SURREALS_FUNC1( asin, asin(z.v_),  double(1)/sqrt(1 - z.v_*z.v_) )
SURREALS_FUNC1( atan, atan(z.v_),  double(1)/(1 + z.v_*z.v_) )

SURREALS_FUNC2( atan2, atan2(z1.v_, z2.v_),  double(1)/(z1.v_*z1.v_ + z2.v_*z2.v_) )

// hyperbolic functions <cmath>

SURREALS_FUNC1( cosh, cosh(z.v_), sinh(z.v_) )
SURREALS_FUNC1( sinh, sinh(z.v_), cosh(z.v_) )
SURREALS_FUNC1( tanh, tanh(z.v_), double(1)/(cosh(z.v_)*cosh(z.v_)) )

// exp and log functions <cmath>

SURREALS_FUNC1( exp, exp(z.v_), exp(z.v_) )
SURREALS_FUNC1( expm1, expm1(z.v_), exp(z.v_) )
SURREALS_FUNC1( log, log(z.v_), double(1)/z.v_ )
SURREALS_FUNC1( log10, log10(z.v_), double(1)/(z.v_*log(10.)) )
SURREALS_FUNC1( log1p, log1p(z.v_), double(1)/( 1 + z.v_ ) )

// power functions <cmath>

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
pow( const SurrealS<N,T>& a, const SurrealS<N,T>& b)
{
  // many sticky points were derivative is undefined or infinite
  // badness if 0 <= b < 1 and a == 0
  T powab = pow(a.v_,b.v_);
  T tmp1 = (a.v_ == 0) ? ((b.v_ == 1) ? T(1) : T(0)) : b.v_*pow(a.v_, b.v_ - 1);
  T tmp2 = (a.v_ == 0) ? T(0) : powab*log(a.v_);
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp1*a.d_[i] + tmp2*b.d_[i];
  c.v_ = powab;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
pow( const SurrealS<N,T>& a, const double& b)
{
  T tmp = (a.v_ == 0) ? ((b == 1) ? T(1) : T(0)) : b*pow(a.v_, b - 1);
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp*a.d_[i];
  c.v_ = pow(a.v_,b);
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
pow( const SurrealS<N,T>& a, const int& b)
{
  T tmp = (a.v_ == 0) ? ((b == 1) ? T(1) : T(0)) : b*pow(a.v_, static_cast<double>(b - 1));
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp*a.d_[i];
  c.v_ = pow(a.v_, static_cast<double>(b));
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
pow( const double& a, const SurrealS<N,T>& b)
{
  T powab = pow(a, b.v_);
  T tmp = (a == 0) ? T(0) : powab*std::log(a);
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp*b.d_[i];
  c.v_ = powab;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
pow( const int& a, const SurrealS<N,T>& b)
{
  T powab = pow(static_cast<double>(a), b.v_);
  T tmp = (a == 0) ? T(0) : powab*std::log(static_cast<double>(a));
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp*b.d_[i];
  c.v_ = powab;
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
sqrt( const SurrealS<N,T>& z )
{
  T sqrtv=sqrt(z.v_);
  if (sqrtv == 0)
  {
    return SurrealS<N,T>(0, 0);
  }
  else
  {
    T tmp = 0.5/sqrtv;
    SurrealS<N,T> c;
    for (int i = 0; i < N; i++)
      c.d_[i] = tmp*z.d_[i];
    c.v_ = sqrtv;
    return c;
  }
}


// rounding functions <cmath>

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
ceil( const SurrealS<N,T>& z )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = 0;
  c.v_ = ceil(z.v_);
  return c;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
floor( const SurrealS<N,T>& z )
{
  SurrealS<N,T> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = 0;
  c.v_ = floor(z.v_);
  return c;
}


// misc functions <cmath>

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
abs( const SurrealS<N,T>& z )
{
  return (z.v_ < 0) ? -z : SurrealS<N,T>(z);
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
fabs( const SurrealS<N,T>& z )
{
  return (z.v_ < 0) ? -z : SurrealS<N,T>(z);
}


template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
max( const SurrealS<N,T>& a, const SurrealS<N,T>& b )
{
  return a.value() > b.value() ? a : b;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
max( const double& a, const SurrealS<N,T>& b )
{
  if ( a > b.value() )
    return a;
  else
    return b;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
max( const SurrealS<N,T>& a, const double& b )
{
  if ( a.value() > b )
    return a;
  else
    return b;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
min( const SurrealS<N,T>& a, const SurrealS<N,T>& b )
{
  return a.value() < b.value() ? a : b;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
min( const double& a, const SurrealS<N,T>& b )
{
  if ( a < b.value() )
    return a;
  else
    return b;
}

template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>
min( const SurrealS<N,T>& a, const double& b )
{
  if ( a.value() < b )
    return a;
  else
    return b;
}

// I/O

template<int N, class T>
std::istream&
operator>>( std::istream& is, SurrealS<N,T>& z )
{
  double v = 0;
  double d[10] = {0};
  char c = 0;
  int n = 0;

  is >> c;
  if (c == '(')
  {
    is >> v;

    is >> c;
    bool done = false;
    while (! done)
    {
      if (c != ')') is.clear(std::ios::badbit);
      if (c == ',')
      {
        is >> d[n]; n++;
      }
      else if (c == ')')
      {
        done = true;
      }
    }
  }
  else
  {
    is.putback(c);
    is >> v;
  }

  if (is) z = SurrealS<N,T>(v, d, n);
  return is;
}

template<int N, class T>
std::ostream&
operator<<( std::ostream& os, const SurrealS<N,T>& z )
{
  os << '(' << z.value() << ';';
  for (int i = 0; i < N - 1; i++)
    os << z.deriv(i) << ',';
  os << z.deriv(N - 1) << ')';
  return os;
}


//Clean up macro definitions
#undef SURREALS_FUNC1
#undef SURREALS_FUNC2


#endif // SURREALS_TRAD_H
