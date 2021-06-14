// Modified from Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2021, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALS_H
#define SURREALS_H

// Old versions of visual studio are missing these functions
#if defined(_MSC_VER)
#if _MSC_VER < 1800
#include <cmath>
inline double expm1( double x ) { return exp(x) - 1; }
inline double log1p( double x ) { return log(x + 1); }
#endif
#endif

#define SURREAL_TRAD

#if defined(SURREAL_TRAD)
#include "SurrealS_Trad.h"
#elif defined(SURREAL_LAZY)
#include "SurrealS_Lazy.h"
#else
#error "Please define SURREAL_TRAD or SURREAL_LAZY"
#endif

#undef SURREAL_TRAD

#endif // SURREALS_H
