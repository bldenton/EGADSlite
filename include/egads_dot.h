#ifndef EGADS_DOT_H
#define EGADS_DOT_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Function Prototypes for Sensitivities
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
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

__ProtoExt__ int  EG_makeGeometry_dot(ego context, int oclass, int mtype,
                                      /*@null@*/ ego refGeom, /*@null@*/ const int *ints,
                                      const double *data, const double *data_dot, egObject **geom);
__ProtoExt__ int  EG_setGeometry_dot( ego geom, int oclass, int mtype,
                                      /*@null@*/ const int *ints,
                                      /*@null@*/ const double *reals,
                                      /*@null@*/ const double *reals_dot );
__ProtoExt__ int  EG_getGeometry_dot( const ego geom, double **reals,
                                      double **reals_dot );
__ProtoExt__ int  EG_hasGeometry_dot( const ego geom );
__ProtoExt__ int  EG_copyGeometry_dot( const ego obj,
                                       /*@null@*/ const double *mat,
                                       /*@null@*/ const double *mat_dot,
                                       ego copy );
__ProtoExt__ int  EG_evaluate_dot( const ego geom,
                                   /*@null@*/ const double *params,
                                   /*@null@*/ const double *params_dot,
                                   double *results, double *results_dot );
__ProtoExt__ int  EG_approximate_dot( ego bspline, int mDeg, double tol,
                                      const int *sizes,
                                      const double *xyzs, const double *xyzs_dot );
__ProtoExt__ int EG_skinning_dot(ego surface, int nCurves, ego *curves);

/* topology functions */

__ProtoExt__ int  EG_makeTopology_dot(egObject *context, /*@null@*/ egObject *geom,
                                      int oclass, int mtype,
                                      /*@null@*/ double *limits, /*@null@*/ double *limits_dot,
                                      int nChildren, /*@null@*/ egObject **children,
                                      /*@null@*/ int *senses, egObject **topo);

__ProtoExt__ int  EG_makeSolidBody_dot( ego body, int stype, const double *data,
                                        const double *data_dot );

__ProtoExt__ int  EG_setRange_dot( ego object, int oclass,
                                   const double *range, const double *range_dot );
__ProtoExt__ int  EG_getRange_dot( const ego geom,
                                   double *range, double *range_dot, int *periodic );

__ProtoExt__ int EG_makeFace_dot(ego face, ego object,
                                 const double *limits,
                                 const double *limits_dot);

/* tessellation functions */

__ProtoExt__ int  EG_tessMassProps_dot( const ego tess, double *xyz_dot,
                                        double *props, double *props_dot );

/* high level functions */

__ProtoExt__ int  EG_extrude_dot( ego body, const ego src,
                                  double dist, double dist_dot,
                                  const double *dir, const double *dir_dot );
__ProtoExt__ int  EG_ruled_dot( ego body, int nSection, const ego *sections );
__ProtoExt__ int  EG_blend_dot( ego body, int nSection, const ego *sections,
                                /*@null@*/ double *rc1,
                                /*@null@*/ double *rc1_dot,
                                /*@null@*/ double *rcN,
                                /*@null@*/ double *rcN_dot );

#ifdef __cplusplus
}

/* include Surreal for automatic differentiation */
#include "Surreal/SurrealS.h"

/* Surreal geometry functions */

int  EG_makeGeometry( ego context, int oclass, int mtype,
                      ego refGeo, const int *ivec,
                      const SurrealS<1> *rvec, egObject **geom );
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
int  EG_approximate_dot( ego bspline, int maxdeg, double tol,
                         const int *sizes,
                         const SurrealS<1> *data );

/* Surreal topology functions */

int  EG_setRange_dot( ego geom, int oclass, const SurrealS<1> *range);
int  EG_getRange( const egObject *geom, SurrealS<1> *range, int *periodic );

#endif

#endif /* EGADS_DOT_H */
