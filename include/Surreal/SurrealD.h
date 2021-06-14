// Modified from Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2021, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALD_H
#define SURREALD_H

#define SURREAL_LAZY

#if defined(SURREAL_LAZY)
#include "SurrealD_Lazy.h"
#elif defined(SURREAL_TRAD)
#include "SurrealD_Trad.h"
#elif defined(SURREAL_REVERSE)
#include "SurrealD_Lazy.h"
#else
#error "Please define SURREAL_TRAD, SURREAL_LAZY or SURREAL_REVERSE"
#endif

#undef SURREAL_LAZY

#endif // SURREALD_H
