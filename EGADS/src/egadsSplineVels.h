#ifndef EGADSSPLINEVELS_H
#define EGADSSPLINEVELS_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Blend/Rule Derivative Functions Header
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"

typedef struct {
  void *usrData;
  int (*velocityOfRange)( void* usrData, const ego *secs, int isec, ego edge,
                          double *trange, double *trange_dot );
  int (*velocityOfNode)( void* usrData, const ego *secs, int isec, ego node, ego edge,
                         double *xyz, double *xyz_dot );
  int (*velocityOfEdge)( void* usrData, const ego *secs, int isec, ego edge,
                         const int jmax, const double *ts, const double *ts_dot,
                         double *xyzs, double *xyzs_dot,
                         double *tbeg, double *tbeg_dot,
                         double *tend, double *tend_dot );
  int (*velocityOfBspline)( void* usrData, const ego *secs, int isec, ego edge, ego geom,
                            int **ivec, double **rvec, double **rvec_dot );
} egadsSplineVels;

#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__
#else
#define __ProtoExt__ extern
#endif

__ProtoExt__ int  EG_blend_vels( int nsec, const ego *secs,
                                 /*@null@*/ const double *rc1,
                                 /*@null@*/ const double *rc1_dot,
                                 /*@null@*/ const double *rcN,
                                 /*@null@*/ const double *rcN_dot,
                                 egadsSplineVels *vels, ego result );

__ProtoExt__ int  EG_ruled_vels( int nsec, const ego *secs,
                                 egadsSplineVels *vels, ego result );

#ifdef __cplusplus
}
#endif

#endif
