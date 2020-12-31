#ifndef EGADS_DOT_H
#define EGADS_DOT_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Function Prototypes for Sensitivities
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "egadsTypes.h"

#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__
#else
#define __ProtoExt__ extern
#endif

/* geometry functions */

__ProtoExt__ int  EG_setGeometry_dot( ego geom, int oclass, int mtype,
                                      /*@null@*/ const int *ivec,
                                      /*@null@*/ const double *rvec,
                                      /*@null@*/ const double *rvec_dot );
__ProtoExt__ int  EG_getGeometry_dot( const ego geom, double **rvec,
                                      double **rvec_dot );
__ProtoExt__ int  EG_hasGeometry_dot( const ego geom );
__ProtoExt__ int  EG_copyGeometry_dot( const ego obj,
                                       /*@null@*/ const double *xform,
                                       /*@null@*/ const double *xform_dot,
                                       ego copy );
__ProtoExt__ int  EG_evaluate_dot( const ego geom,
                                   /*@null@*/ const double *param,
                                   /*@null@*/ const double *param_dot,
                                   double *results, double *results_dot );
__ProtoExt__ int  EG_approximate_dot( ego bspline, int maxdeg, double tol,
                                      const int *sizes,
                                      const double *data, const double *data_dot );

/* topology functions */

__ProtoExt__ int  EG_makeSolidBody_dot( ego body, int stype, const double *rvec,
                                        const double *rvec_dot );

/* tessellation functions */

__ProtoExt__ int  EG_tessMassProps_dot( const ego tess, double *xyz_dot,
                                        double *props, double *props_dot );

/* high level functions */

__ProtoExt__ int  EG_extrude_dot( ego body, const ego src,
                                  double dist, double dist_dot,
                                  const double *dir, const double *dir_dot );
__ProtoExt__ int  EG_ruled_dot( ego body, int nsec, const ego *secs );
__ProtoExt__ int  EG_blend_dot( ego body, int nsec, const ego *secs,
                                /*@null@*/ double *rc1,
                                /*@null@*/ double *rc1_dot,
                                /*@null@*/ double *rcN,
                                /*@null@*/ double *rcN_dot );

__ProtoExt__ int  EG_getRange_dot( const egObject *geom,
                                   double *range, double *range_dot, int *periodic );

#ifdef __cplusplus
}

/* include Surreal for automatic differentiation */
#include "Surreal/SurrealS.h"

/* Surreal geometry functions */

int  EG_getGeometry( const ego geom, int *oclass, int *mtype,
                     ego *refGeom, int **ivec, SurrealS<1> **rvec );
int  EG_setGeometry_dot( ego geom, int oclass, int mtype,
                         /*@null@*/ const int *ivec, SurrealS<1> *rvec_dot );
int  EG_getGeometry_dot( const egObject *obj, SurrealS<1> **data_dot );
int  EG_copyGeometry_dot( const egObject *obj,
                          /*@null@*/ const SurrealS<1> *xform,
                          ego copy );
int  EG_evaluate( const egObject *geom, /*@null@*/ const SurrealS<1> *param,
                  SurrealS<1> *result );

int  EG_getRange( const egObject *geom, SurrealS<1> *range, int *periodic );
#endif

#endif /* EGADS_DOT_H */
