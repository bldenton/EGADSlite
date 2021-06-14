// Modified from Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2021, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALS_LAZY_H
#define SURREALS_LAZY_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <iostream>
#include <cassert>

#include <string>

#include "always_inline.h"

class SurrealSTypeBase {};

template< class Derived, class T >
struct SurrealSType : SurrealSTypeBase
{
  //A convenient method for casting to the derived type
  ALWAYS_INLINE const Derived& cast() const { return static_cast<const Derived&>(*this); }

  //A simple way to call value without having to case first
  ALWAYS_INLINE T value() const { return cast().value(); }
};

namespace SurrealSExpr
{
template<class L, class R, class T > class OpMul;
class OpMul_impl;
}


//----------------------------------------------------------------------------//
// SurrealS:  value, N derivatives
//
// Operators with Lazy Expressions
//
// statically defined derivative array
//----------------------------------------------------------------------------//
template<int N_, class T = double>
class SurrealS : public SurrealSType< SurrealS<N_, T>, T >
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
  ALWAYS_INLINE SurrealS( const double v0, const double d0[], int n );
  ALWAYS_INLINE SurrealS( const double v0, const double& d0 );
  template<class Expr>
  ALWAYS_INLINE SurrealS( const SurrealSType<Expr, T>& r ) : v_(0) { operator=(r); }
  ALWAYS_INLINE ~SurrealS() {}

  ALWAYS_INLINE int size() const { return N; }

  // value accessor operators
  ALWAYS_INLINE       T& value()       { return v_; }
  ALWAYS_INLINE const T& value() const { return v_; }

  // derivative accessor operators
  ALWAYS_INLINE       T& deriv( int i=0 )       { return d_[i]; }
  ALWAYS_INLINE const T& deriv( int i=0 ) const { return d_[i]; }

  // assignment
  ALWAYS_INLINE SurrealS& operator=( const SurrealS& );
  ALWAYS_INLINE SurrealS& operator=( const int& );
  ALWAYS_INLINE SurrealS& operator=( const double& );

  template<class Expr> ALWAYS_INLINE SurrealS& operator= ( const SurrealSType<Expr, T>& );
  template<class Expr> ALWAYS_INLINE SurrealS& operator+=( const SurrealSType<Expr, T>& );
  template<class Expr> ALWAYS_INLINE SurrealS& operator-=( const SurrealSType<Expr, T>& );

  // unary operators; no side effects
  ALWAYS_INLINE const SurrealS& operator+() const;

  // binary accumulation operators
  ALWAYS_INLINE SurrealS& operator+=( const double& );
  ALWAYS_INLINE SurrealS& operator-=( const double& );
  ALWAYS_INLINE SurrealS& operator*=( const SurrealS& );
  ALWAYS_INLINE SurrealS& operator*=( const double& );
  ALWAYS_INLINE SurrealS& operator/=( const SurrealS& );
  ALWAYS_INLINE SurrealS& operator/=( const double& );

#if 0
  // classification functions <cmath>
  friend bool isfinite( const SurrealS& );
  friend bool isinf( const SurrealS& );
  friend bool isnan( const SurrealS& );
#endif

  // input/output
  template<int M> friend std::istream& operator>>( std::istream&, SurrealS<M, T>& );

protected:
  T d_[N];   // derivative array
  T v_;      // value
};

//Constructors

template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const SurrealS& z )
{
  v_ = z.v_;
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
SurrealS<N,T>::SurrealS( const double v0, const double d0[], int n )
{
  assert( n == N );
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = d0[i];
}
template<int N, class T>
ALWAYS_INLINE
SurrealS<N,T>::SurrealS( const double v0, const double& d0 )
{
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = d0;
}


namespace SurrealSExpr
{

// Lazy expressions


// Addition and Subtraction

template<class L, class R, class T>
class OpAdd : public SurrealSType< OpAdd<L,R,T>, T >
{
public:
  static const int N = L::N;
  static_assert( L::N == R::N, "Sizes must match" );

  ALWAYS_INLINE
  OpAdd(const L& Ll, const R& Rr) : Ll(Ll), Rr(Rr) {}

  ALWAYS_INLINE T value() const { return Ll.value() + Rr.value(); }
  ALWAYS_INLINE T deriv(const int& i) const { return Ll.deriv(i) + Rr.deriv(i); }

  ALWAYS_INLINE const OpAdd&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return Ll.size(); }
private:
  const L& Ll;
  const R& Rr;
};

}

template<class L, class R, class T>
ALWAYS_INLINE SurrealSExpr::OpAdd<L,R,T>
operator+(const SurrealSType<L,T>& Ll, const SurrealSType<R,T>& Rr)
{
  return SurrealSExpr::OpAdd<L,R,T>( Ll.cast(), Rr.cast() );
}

namespace SurrealSExpr
{

template<class L, class R, class T>
class OpSub : public SurrealSType< OpSub<L,R,T>, T >
{
public:
  static const int N = L::N;
  static_assert( L::N == R::N, "Sizes must match" );

  ALWAYS_INLINE
  OpSub(const L& Ll, const R& Rr) : Ll(Ll), Rr(Rr) {}

  ALWAYS_INLINE T value() const { return Ll.value() - Rr.value(); }
  ALWAYS_INLINE T deriv(const int& i) const { return Ll.deriv(i) - Rr.deriv(i); }

  ALWAYS_INLINE const OpSub&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return Ll.size(); }
private:
  const L& Ll;
  const R& Rr;
};

}

template<class L, class R, class T>
ALWAYS_INLINE SurrealSExpr::OpSub<L,R,T>
operator-(const SurrealSType<L,T>& Ll, const SurrealSType<R,T>& Rr)
{
  return SurrealSExpr::OpSub<L,R,T>( Ll.cast(), Rr.cast() );
}

//Addition and Subtraction with scalar quantities

namespace SurrealSExpr
{

template<class Expr, class T>
class OpScalar : public SurrealSType< OpScalar<Expr, T>, T >
{
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  OpScalar(const Expr& e, const double esgn, const double s) : e(e), esgn(esgn), s(s) {}

  ALWAYS_INLINE T value() const { return esgn*e.value() + s; }
  ALWAYS_INLINE T deriv(const int& i) const { return esgn*e.deriv(i); }

  ALWAYS_INLINE const OpScalar&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const double esgn;
  const double s;
};

}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr,T>
operator+(const SurrealSType<Expr,T>& e, const double& s)
{
  return SurrealSExpr::OpScalar<Expr,T>( e.cast(), 1, s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr,T>
operator+(const double& s, const SurrealSType<Expr,T>& e)
{
  return SurrealSExpr::OpScalar<Expr,T>( e.cast(), 1, s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr,T>
operator-(const SurrealSType<Expr,T>& e, const double& s)
{
  return SurrealSExpr::OpScalar<Expr,T>( e.cast(), 1, -s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr,T>
operator-(const double& s, const SurrealSType<Expr,T>& e)
{
  return SurrealSExpr::OpScalar<Expr,T>( e.cast(), -1, s );
}


//Multiplication with SurrealSs

namespace SurrealSExpr
{

template<class ExprL, class ExprR, class T>
class OpMul : public SurrealSType< OpMul<ExprL, ExprR, T>, T >
{
public:
  static const int N = ExprL::N;
  static_assert( ExprL::N == ExprR::N, "Sizes must match" );

  ALWAYS_INLINE
  OpMul(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), eL_val(eL.value()), eR_val(eR.value()) {}

  ALWAYS_INLINE T value() const { return eL_val*eR_val; }
  ALWAYS_INLINE T deriv(const int& i) const { return eL_val*eR.deriv(i) + eL.deriv(i)*eR_val; }

  ALWAYS_INLINE const OpMul&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const T eL_val, eR_val;
};

template<class ExprL, class T>
class OpMul<ExprL, double, T> : public SurrealSType< OpMul<ExprL, double, T>, T >
{
public:
  static const int N = ExprL::N;

  ALWAYS_INLINE
  OpMul(const ExprL& e, const double s) : e(e), s(s) {}

  ALWAYS_INLINE T value() const { return e.value()*s; }
  ALWAYS_INLINE T deriv(const int& i) const { return e.deriv(i)*s; }

  ALWAYS_INLINE const OpMul&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }

  const ExprL& e;
  const double s;
};
}


//=============================================================================
template<class ExprL, class ExprR, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<ExprL, ExprR, T>
operator*(const SurrealSType<ExprL, T>& z1, const SurrealSType<ExprR, T>& z2)
{
  return SurrealSExpr::OpMul<ExprL, ExprR, T>( z1.cast(), z2.cast() );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double, T>
operator*(const SurrealSType<Expr, T>& e, const double& s)
{
  return SurrealSExpr::OpMul<Expr, double, T>( e.cast(), s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double, T>
operator/(const SurrealSType<Expr, T>& e, const double& s)
{
  return SurrealSExpr::OpMul<Expr, double, T>( e.cast(), double(1)/s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double, T>
operator*(const double& s, const SurrealSType<Expr, T>& e)
{
  return SurrealSExpr::OpMul<Expr, double, T>( e.cast(), s );
}

//=============================================================================
// This is a special case when multiplies scalars are multiplying from two sides, i.e. B = 2*A*2;
// This reduces the complexity of the expression tree and hence reduces code bloat
template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double, T>
operator*(const SurrealSExpr::OpMul<Expr, double, T>& MulScal, const double& s)
{
  return SurrealSExpr::OpMul<Expr, double, T>( MulScal.e, MulScal.s*s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double, T>
operator/(const SurrealSExpr::OpMul<Expr, double, T>& MulScal, const double& s)
{
  return SurrealSExpr::OpMul<Expr, double, T>( MulScal.e, MulScal.s/s );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double, T>
operator*(const double& s, const SurrealSExpr::OpMul<Expr, double, T>& MulScal)
{
  return SurrealSExpr::OpMul<Expr, double, T>( MulScal.e, MulScal.s*s );
}


//=============================================================================
//Division with SurrealSs

namespace SurrealSExpr
{

template<class ExprL, class ExprR, class T>
class OpDiv : public SurrealSType< OpDiv<ExprL, ExprR, T>, T >
{
public:
  static const int N = ExprL::N;
  static_assert( ExprL::N == ExprR::N, "Sizes must match" );

  ALWAYS_INLINE
  OpDiv(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), eL_val(eL.value()), eR_val(eR.value())
                                          , vali(1/(eR_val*eR_val)) {}

  ALWAYS_INLINE T value() const { return eL_val/eR_val; }
  ALWAYS_INLINE T deriv(const int& i) const { return (eR_val*eL.deriv(i) - eR.deriv(i)*eL_val)*vali; }

  ALWAYS_INLINE const OpDiv&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const T eL_val, eR_val, vali;
};

}

template< class ExprL, class ExprR, class T >
ALWAYS_INLINE SurrealSExpr::OpDiv<ExprL, ExprR, T>
operator/(const SurrealSType<ExprL, T>& eL, const SurrealSType<ExprR, T>& eR)
{
  return SurrealSExpr::OpDiv<ExprL, ExprR, T>( eL.cast(), eR.cast() );
}


namespace SurrealSExpr
{

template<class Expr, class T>
class OpDivScalarNumerator : public SurrealSType< OpDivScalarNumerator<Expr, T>, T >
{
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  OpDivScalarNumerator(const Expr& e, const double& s) : e(e), s(s), e_val(e.value()), se_val2i(s/(e_val*e_val)) {}

  ALWAYS_INLINE T value() const { return s/e_val; }
  ALWAYS_INLINE T deriv(const int& i) const { return -se_val2i*e.deriv(i); }

  ALWAYS_INLINE const OpDivScalarNumerator&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const double s;
  const T e_val, se_val2i;
};

}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpDivScalarNumerator<Expr,T>
operator/(const double& s, const SurrealSType<Expr,T>& e)
{
  return SurrealSExpr::OpDivScalarNumerator<Expr,T>( e.cast(), s );
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
SurrealS<N,T>::operator=( const int& r )
{
  v_ = r;
  for (int i = 0; i < N; i++)
    d_[i] = 0;
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

template<int N, class T>
template< class Expr >
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator=( const SurrealSType<Expr,T>& r )
{
  const Expr& Tree = r.cast();

  static_assert( N == Expr::N, "Sizes must match" );

  for (int i = 0; i < N; ++i)
    d_[i] = Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ = Tree.value();

  return *this;
}

template<int N, class T>
template< class Expr >
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator+=( const SurrealSType<Expr,T>& r )
{
  const Expr& Tree = r.cast();

  static_assert( N == Expr::N, "Sizes must match" );

  for (int i = 0; i < N; ++i)
    d_[i] += Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ += Tree.value();

  return *this;
}

template<int N, class T>
template< class Expr >
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator-=( const SurrealSType<Expr,T>& r )
{
  const Expr& Tree = r.cast();

  static_assert( N == Expr::N, "Sizes must match" );

  for (int i = 0; i < N; ++i)
    d_[i] -= Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ -= Tree.value();

  return *this;
}


// unary operators; no side effects

template<int N, class T>
ALWAYS_INLINE const SurrealS<N,T>&
SurrealS<N,T>::operator+() const
{
  return *this;
}

template< class Expr, class T >
ALWAYS_INLINE const SurrealSExpr::OpMul<Expr, double, T>
operator-(SurrealSType<Expr,T> const& e)
{
  return SurrealSExpr::OpMul<Expr, double, T>( e.cast(), -1 );
}

template< class Expr, class T >
ALWAYS_INLINE const SurrealSExpr::OpMul<Expr, double, T>
operator-(SurrealSExpr::OpMul<Expr, double, T> const& Mul)
{
  return SurrealSExpr::OpMul<Expr, double, T>( Mul.e, -1*Mul.s );
}

// binary accumulation operators


template<int N, class T>
ALWAYS_INLINE SurrealS<N,T>&
SurrealS<N,T>::operator+=( const double& r )
{
  v_ += r;
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

// relational operators

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator==( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() == rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator==( const SurrealSType<Expr, T>& lhs, const double& rhs )
{
  return lhs.value() == rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator==( const double& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs == rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator!=( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() != rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator!=( const SurrealSType<Expr, T>& lhs, const double& rhs )
{
  return lhs.value() != rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator!=( const double& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs != rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator>( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() > rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator>( const SurrealSType<Expr, T>& lhs, const double& rhs )
{
  return lhs.value() > rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator>( const double& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs > rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator<( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() < rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator<( const SurrealSType<Expr, T>& lhs, const double& rhs )
{
  return lhs.value() < rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator<( const double& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs < rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator>=( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() >= rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator>=( const SurrealSType<Expr, T>& lhs, const double& rhs )
{
  return lhs.value() >= rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator>=( const double& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs >= rhs.value();
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE bool
operator<=( const SurrealSType<ExprL, T>& lhs, const SurrealSType<ExprR, T>& rhs )
{
  return lhs.value() <= rhs.value();
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator<=( const SurrealSType<Expr, T>& lhs, const double& rhs )
{
  return lhs.value() <= rhs;
}

template<class Expr, class T>
ALWAYS_INLINE bool
operator<=( const double& lhs, const SurrealSType<Expr, T>& rhs )
{
  return lhs <= rhs.value();
}


//Functions for SurrealSs
#define SURREALS_FUNC1( NAME, FUNC, DERIV ) \
namespace SurrealSExpr \
{  \
template<class Expr, class T> \
class SurrealS_ ## NAME : public SurrealSType< SurrealS_ ## NAME<Expr, T>, T > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  static const int N = Expr::N; \
  \
  ALWAYS_INLINE \
  SurrealS_ ## NAME(const Expr& e) : e(e), z(e.value()), der(DERIV) {} \
  \
  ALWAYS_INLINE T value() const { return FUNC; } \
  ALWAYS_INLINE T deriv(const int& i) const { return der*e.deriv(i); } \
  \
  ALWAYS_INLINE const SurrealS_ ## NAME& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return e.size(); } \
private: \
  const Expr& e; \
  const T z, der; \
}; \
} \
\
template<class Expr, class T> \
ALWAYS_INLINE SurrealSExpr::SurrealS_ ## NAME<Expr, T> \
NAME(const SurrealSType<Expr, T>& z) { return SurrealSExpr::SurrealS_ ## NAME<Expr, T>( z.cast() ); }


#define SURREALS_FUNC2( NAME, FUNC, DERIV ) \
namespace SurrealSExpr \
{  \
template<class ExprL, class ExprR, class T> \
class SurrealS_ ## NAME : public SurrealSType< SurrealS_ ## NAME<ExprL, ExprR, T>, T > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  static const int N = ExprL::N; \
  static_assert( ExprL::N == ExprR::N, "Sizes must match" ); \
  \
  ALWAYS_INLINE \
  SurrealS_ ## NAME(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), z1(eL.value()), z2(eR.value()), \
                                                                 der(DERIV) {} \
  \
  ALWAYS_INLINE T value() const { return FUNC; } \
  ALWAYS_INLINE T deriv(const int& i) const { return der*(z2*eL.deriv(i) - z1*eR.deriv(i)); } \
  \
  ALWAYS_INLINE const SurrealS_ ## NAME& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return eL.size(); } \
private: \
  const ExprL& eL; \
  const ExprR& eR; \
  const T z1, z2, der; \
}; \
  \
} \
\
template<class ExprL, class ExprR, class T> \
ALWAYS_INLINE SurrealSExpr::SurrealS_ ## NAME<ExprL, ExprR, T> \
NAME(const SurrealSType<ExprL, T>& z1, const SurrealSType<ExprR, T>& z2) \
{ return SurrealSExpr::SurrealS_ ## NAME<ExprL, ExprR, T>( z1.cast(), z2.cast() ); }

// trig functions <cmath>

SURREALS_FUNC1( cos, cos(z), -sin(z) )
SURREALS_FUNC1( sin, sin(z),  cos(z) )
SURREALS_FUNC1( tan, tan(z),  double(1)/(cos(z)*cos(z)) )
SURREALS_FUNC1( acos, acos(z), -double(1)/sqrt(1 - z*z) )
SURREALS_FUNC1( asin, asin(z),  double(1)/sqrt(1 - z*z) )
SURREALS_FUNC1( atan, atan(z),  double(1)/(1 + z*z) )

SURREALS_FUNC2( atan2, atan2(z1, z2),  double(1)/(z1*z1 + z2*z2) )

// hyperbolic functions <cmath>

SURREALS_FUNC1( cosh, cosh(z), sinh(z) )
SURREALS_FUNC1( sinh, sinh(z), cosh(z) )
SURREALS_FUNC1( tanh, tanh(z), double(1)/(cosh(z)*cosh(z)) )

// exp and log functions <cmath>

SURREALS_FUNC1( exp, exp(z), exp(z) )
SURREALS_FUNC1( expm1, expm1(z), exp(z) )
SURREALS_FUNC1( log, log(z), double(1)/z )
SURREALS_FUNC1( log10, log10(z), double(1)/(z*log(10.)) )
SURREALS_FUNC1( log1p, log1p(z), double(1)/( 1 + z ) )

// power functions <cmath>

namespace SurrealSExpr
{

template<class ExprL, class ExprR, class T>
class SurrealS_pow : public SurrealSType< SurrealS_pow<ExprL, ExprR, T>, T >
{ /*This is for functions when the argument is an expression*/
public:
  static const int N = ExprL::N;
  static_assert( ExprL::N == ExprR::N, "Sizes must match" );

  ALWAYS_INLINE
  SurrealS_pow(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), a(eL.value()), b(eR.value()),
                                                   powab(pow(a,b)),
                                                   tmp1( (a == T(0)) ? ((b == T(1)) ? T(1) : T(0)) : b*pow(a, b - 1) ),
                                                   tmp2( (a == T(0)) ? T(0) : powab*log(a) ) {}

  ALWAYS_INLINE T value() const { return powab; }
  ALWAYS_INLINE T deriv(const int& i) const { return tmp1*eL.deriv(i) + tmp2*eR.deriv(i); }

  ALWAYS_INLINE const SurrealS_pow&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const T a, b, powab, tmp1, tmp2;
};

template<class ExprL, class T>
class SurrealS_pow<ExprL, double, T> : public SurrealSType< SurrealS_pow<ExprL, double, T>, T >
{ /*This is optimized when the argument is SurrealS and double*/
public:
  static const int N = ExprL::N;

  ALWAYS_INLINE
  SurrealS_pow(const ExprL& eL, const double& b) : eL(eL), a(eL.value()),
                                                   powab(pow(a,b)),
                                                   tmp1( (a == T(0)) ? ((b == 1) ? T(1) : T(0)) : b*pow(a, b - 1) ) {}

  ALWAYS_INLINE T value() const { return powab; }
  ALWAYS_INLINE T deriv(const int& i) const { return tmp1*eL.deriv(i); }

  ALWAYS_INLINE const SurrealS_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const T a, powab, tmp1;
};


template<class ExprR, class T>
class SurrealS_pow<double, ExprR, T> : public SurrealSType< SurrealS_pow<double, ExprR, T>, T >
{ /*This is optimized when the argument is a double and SurrealS*/
public:
  static const int N = ExprR::N;

  ALWAYS_INLINE
  SurrealS_pow(const double& a, const ExprR& eR) : eR(eR), b(eR.value()),
                                                   powab( (b == T(0)) ? T(1) : pow(a,b) ),
                                                   tmp2( (a == 0) ? T(0) : powab*log(a) ) {}

  ALWAYS_INLINE T value() const { return powab; }
  ALWAYS_INLINE T deriv(const int& i) const { return tmp2*eR.deriv(i); }

  template<int I>
  ALWAYS_INLINE T    get_deriv() const { return tmp2*eR.template get_deriv<I>(); }

  ALWAYS_INLINE const SurrealS_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eR.size(); }
private:
  const ExprR& eR;
  const T b, powab, tmp2;

};

}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE SurrealSExpr::SurrealS_pow<ExprL, ExprR, T>
pow(const SurrealSType<ExprL, T>& a, const SurrealSType<ExprR, T>& b)
{
  return SurrealSExpr::SurrealS_pow<ExprL, ExprR, T>( a.cast(), b.cast() );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::SurrealS_pow<Expr, double, T>
pow(const SurrealSType<Expr,T>& a, const double& b )
{
  return SurrealSExpr::SurrealS_pow<Expr, double, T>( a.cast(), b );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::SurrealS_pow<double, Expr, T>
pow(const double& a, const SurrealSType<Expr, T>& b)
{
  return SurrealSExpr::SurrealS_pow<double, Expr, T>( a, b.cast() );
}


namespace SurrealSExpr
{

template<class Expr, class T>
class SurrealS_sqrt : public SurrealSType< SurrealS_sqrt<Expr, T>, T >
{ /*This is optimized when the argument is an Expression*/
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  SurrealS_sqrt(const Expr& e) : e(e), sqrtv( sqrt(e.value()) ), tmp( sqrtv == 0 ? sqrtv : 0.5/sqrtv ) {}

  ALWAYS_INLINE T value() const { return sqrtv; }
  ALWAYS_INLINE T deriv(const int& i) const { return tmp*e.deriv(i); }

  template<int I>
  ALWAYS_INLINE T get_deriv() const { return tmp*e.template get_deriv<I>(); }

  ALWAYS_INLINE const SurrealS_sqrt
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const T sqrtv, tmp;
};

}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::SurrealS_sqrt<Expr, T>
sqrt(const SurrealSType<Expr, T>& z)
{
  return SurrealSExpr::SurrealS_sqrt<Expr, T>( z.cast() );
}


// rounding functions <cmath>

SURREALS_FUNC1( ceil, ceil(z), 0 )
SURREALS_FUNC1( floor, floor(z), 0 )

// misc functions <cmath>

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double, T>
abs( const SurrealSType<Expr, T>& z )
{
  return (z.value() < 0) ?
         SurrealSExpr::OpMul<Expr, double, T>( z.cast(), -1 ) :
         SurrealSExpr::OpMul<Expr, double, T>( z.cast(),  1 );
}

template<class Expr, class T>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double, T>
fabs( const SurrealSType<Expr, T>& z )
{
  return (z.value() < 0) ?
         SurrealSExpr::OpMul<Expr, double, T>( z.cast(), -1 ) :
         SurrealSExpr::OpMul<Expr, double, T>( z.cast(),  1 );
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE SurrealS<ExprL::N,T>
max( const SurrealSType<ExprL, T>& a, const SurrealSType<ExprR, T>& b )
{
  static_assert( ExprL::N == ExprR::N, "Sizes must match" );
  return a.cast().value() > b.cast().value() ? a : b;
}

template<class ExprR, class T>
ALWAYS_INLINE SurrealS<ExprR::N,T>
max( const double& a, const SurrealSType<ExprR, T>& b )
{
  if ( a > b.cast().value() )
    return a;
  else
    return b;
}

template<class ExprL, class T>
ALWAYS_INLINE SurrealS<ExprL::N,T>
max( const SurrealSType<ExprL, T>& a, const double& b )
{
  if ( a.cast().value() > b )
    return a;
  else
    return b;
}

template<class ExprL, class ExprR, class T>
ALWAYS_INLINE SurrealS<ExprL::N,T>
min( const SurrealSType<ExprL, T>& a, const SurrealSType<ExprR, T>& b )
{
  static_assert( ExprL::N == ExprR::N, "Sizes must match" );
  return a.cast().value() < b.cast().value() ? a : b;
}

template<class ExprR, class T>
ALWAYS_INLINE SurrealS<ExprR::N,T>
min( const double& a, const SurrealSType<ExprR, T>& b )
{
  if ( a < b.cast().value() )
    return a;
  else
    return b;
}

template<class ExprL, class T>
ALWAYS_INLINE SurrealS<ExprL::N,T>
min( const SurrealSType<ExprL, T>& a, const double& b )
{
  if ( a.cast().value() < b )
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


template<class Expr, class T>
std::ostream&
operator<<( std::ostream& os, const SurrealSType<Expr, T>& ztype )
{
  const Expr& z = ztype.cast();
  os << '(' << z.value() << ';';
  for (int i = 0; i < Expr::N - 1; i++)
    os << z.deriv(i) << ',';
  os << z.deriv(Expr::N - 1) << ')';
  return os;
}

//Clean up macro definitions
#undef SURREALS_FUNC1
#undef SURREALS_FUNC2

#endif // SURREALS_LAZY_H
