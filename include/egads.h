#ifndef EGADS_H
#define EGADS_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Function Prototypes
 *
 *      Copyright 2011-2018, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <petscsys.h>
#include "egadsTypes.h"

/* memory functions */

PETSC_EXTERN /*@null@*/ /*@out@*/ /*@only@*/
             void *EG_alloc( int nbytes );
PETSC_EXTERN /*@null@*/ /*@only@*/
             void *EG_calloc( int nele, int size );
PETSC_EXTERN /*@null@*/ /*@only@*/
             void *EG_reall( /*@null@*/ /*@only@*/ /*@returned@*/ void *ptr,
                             int nbytes );
PETSC_EXTERN /*@null@*/ /*@only@*/
             char *EG_strdup( /*@null@*/ const char *str );
PETSC_EXTERN void  EG_free( /*@null@*/ /*@only@*/ void *pointer );

/* base-level object functions */

PETSC_EXTERN void EG_revision( int *major, int *minor, const char **OCCrev );
PETSC_EXTERN int  EG_open( ego *context );
PETSC_EXTERN int  EG_loadModel( ego context, int bflg, const char *name,
                                ego *model );
PETSC_EXTERN int  EG_saveModel( const ego model, const char *name );
PETSC_EXTERN int  EG_deleteObject( ego object );
PETSC_EXTERN int  EG_makeTransform( ego context, const double *xform,
                                    ego *oform );
PETSC_EXTERN int  EG_getTransformation( const ego oform, double *xform );
PETSC_EXTERN int  EG_getContext( ego object, ego *context );
PETSC_EXTERN int  EG_setOutLevel( ego context, int outLevel );
PETSC_EXTERN int  EG_getInfo( const ego object, int *oclass, int *mtype,
                              ego *topObj, ego *prev, ego *next );
PETSC_EXTERN int  EG_copyObject( const ego object, /*@null@*/ void *oform,
                                 ego *copy );
PETSC_EXTERN int  EG_flipObject( const ego object, ego *flippedCopy );
PETSC_EXTERN int  EG_close( ego context );
PETSC_EXTERN int  EG_setUserPointer( ego context, void *ptr );
PETSC_EXTERN int  EG_getUserPointer( const ego context, void **ptr );

/* attribute functions */

PETSC_EXTERN int  EG_attributeAdd( ego obj, const char *name, int type, int len,
                                  /*@null@*/ const int    *ints,
                                  /*@null@*/ const double *reals,
                                  /*@null@*/ const char   *str );
PETSC_EXTERN int  EG_attributeDel( ego object, /*@null@*/ const char *name );
PETSC_EXTERN int  EG_attributeNum( const ego obj, int *num );
PETSC_EXTERN int  EG_attributeGet( const ego obj, int index, const char **name,
                                   int *atype, int *len,
                                   /*@null@*/ const int    **ints,
                                   /*@null@*/ const double **reals,
                                   /*@null@*/ const char   **str );
PETSC_EXTERN int  EG_attributeRet( const ego obj, const char *name, int *atype,
                                   int *len, /*@null@*/ const int    **ints,
                                             /*@null@*/ const double **reals,
                                             /*@null@*/ const char   **str );
PETSC_EXTERN int  EG_attributeDup( const ego src, ego dst );

/* geometry functions */

PETSC_EXTERN int  EG_getGeometry( const ego geom, int *oclass, int *mtype,
                                  ego *refGeom, int **ivec, double **rvec );
PETSC_EXTERN int  EG_makeGeometry( ego context, int oclass, int mtype,
                                  /*@null@*/ ego refGeom,
                                  /*@null@*/ const int *ivec,
                                  const double *rvec, ego *geom );
PETSC_EXTERN int  EG_getRange( const ego geom, double *range, int *periodic );
PETSC_EXTERN int  EG_evaluate( const ego geom, const double *param,
                               double *results );
PETSC_EXTERN int  EG_invEvaluate( const ego geom, const double *xyz, double *param,
                                  double *results );
PETSC_EXTERN int  EG_invEvaluateGuess( const ego geom, double *xyz,
                                       double *param, double *results );
PETSC_EXTERN int  EG_arcLength( const ego geom, double t1, double t2,
                                double *alen );
PETSC_EXTERN int  EG_curvature( const ego geom, const double *param,
                                double *results );
PETSC_EXTERN int  EG_approximate( ego context, int maxdeg, double tol,
                                  const int *sizes, const double *xyzs,
                                  ego *bspline );
PETSC_EXTERN int  EG_fitTriangles( ego context, int npts, double *xyzs,
                                   int ntris, const int *tris,
                                   /*@null@*/ const int *tric, double tol,
                                   ego *bspline );
PETSC_EXTERN int  EG_otherCurve( const ego surface, const ego curve,
                                 double tol, ego *newcurve );
PETSC_EXTERN int  EG_isSame( const ego geom1, const ego geom2 );
PETSC_EXTERN int  EG_isoCline( const ego surface, int UV, double val,
                               ego *newcurve );
PETSC_EXTERN int  EG_convertToBSpline( ego geom, ego *bspline );
PETSC_EXTERN int  EG_convertToBSplineRange( ego geom, const double *range,
                                            ego *bspline );
PETSC_EXTERN int  EG_skinning ( ego context, int nCurves, ego *curves,
                                int skinning_degree, ego *surface );

/* topology functions */

PETSC_EXTERN int  EG_tolerance( const ego topo, double *tol );
PETSC_EXTERN int  EG_getTolerance( const ego topo, double *tol );
PETSC_EXTERN int  EG_setTolerance(       ego topo, double  tol );
PETSC_EXTERN int  EG_getTopology( const ego topo, ego *geom, int *oclass,
                                  int *type, /*@null@*/ double *limits,
                                  int *nChildren, ego **children, int **sense );
PETSC_EXTERN int  EG_makeTopology( ego context, /*@null@*/ ego geom, int oclass,
                                   int mtype, /*@null@*/ double *limits,
                                   int nChildren, /*@null@*/ ego *children,
                                   /*@null@*/ int *senses, ego *topo );
PETSC_EXTERN int  EG_makeLoop( int nedge, ego *edges, /*@null@*/ ego geom,
                               double toler, ego *result );
PETSC_EXTERN int  EG_getArea( ego object, /*@null@*/ const double *limits,
                              double *area );
PETSC_EXTERN int  EG_makeFace( ego object, int mtype,
                               /*@null@*/ const double *limits, ego *face );
PETSC_EXTERN int  EG_getBodyTopos( const ego body, /*@null@*/ ego src,
                                   int oclass, int *ntopo,
                                   /*@null@*/ ego **topos );
PETSC_EXTERN int  EG_indexBodyTopo( const ego body, const ego src );
PETSC_EXTERN int  EG_objectBodyTopo( const ego body, int oclass, int index,
                                     ego *obj );
PETSC_EXTERN int  EG_inTopology( const ego topo, const double *xyz );
PETSC_EXTERN int  EG_inFace( const ego face, const double *uv );
PETSC_EXTERN int  EG_getEdgeUV( const ego face, const ego edge, int sense,
                                double t, double *UV );
PETSC_EXTERN int  EG_getEdgeUVs( const ego face, const ego edge, int sense,
                                 int nts, const double *ts, double *UV );
PETSC_EXTERN int  EG_getPCurve( const ego face, const ego edge, int sense,
                                int *mtype, int **ivec, double **rvec );
PETSC_EXTERN int  EG_getWindingAngle( ego edge, double t, double *angle );
PETSC_EXTERN int  EG_getBody( const ego object, ego *body );
PETSC_EXTERN int  EG_makeSolidBody( ego context, int stype, const double *rvec,
                                    ego *body );
PETSC_EXTERN int  EG_getBoundingBox( const ego topo, double *bbox );
PETSC_EXTERN int  EG_getMassProperties( const ego topo, double *result );
PETSC_EXTERN int  EG_isEquivalent( const ego topo1, const ego topo2 );
PETSC_EXTERN int  EG_sewFaces( int nobj, const ego *objs, double toler,
                               int flag, ego *result );
PETSC_EXTERN int  EG_replaceFaces( const ego body, int nobj, ego *objs,
                                   ego *result );
PETSC_EXTERN int  EG_mapBody( const ego sBody,   const ego dBody,
                              const char *fAttr, ego *mapBody );
PETSC_EXTERN int  EG_matchBodyFaces( const ego body1, const ego body2,
                                     double toler, int *nmatch, int **match );

/* tessellation functions */

PETSC_EXTERN int  EG_setTessParam( ego context, int iparam, double value,
                                   double *oldvalue );
PETSC_EXTERN int  EG_makeTessGeom( ego obj, double *params, int *sizes,
                                   ego *tess );
PETSC_EXTERN int  EG_getTessGeom( const ego tess, int *sizes, double **xyz );

PETSC_EXTERN int  EG_makeTessBody( ego object, double *params, ego *tess );
PETSC_EXTERN int  EG_remakeTess( ego tess, int nobj, ego *objs,
                                 double *params );
PETSC_EXTERN int  EG_mapTessBody( ego tess, ego body, ego *mapTess );
PETSC_EXTERN int  EG_locateTessBody( const ego tess, int npt, const int *ifaces,
                                     const double *uv, /*@null@*/ int *itri,
                                     double *results );

PETSC_EXTERN int  EG_getTessEdge( const ego tess, int eIndex, int *len,
                                  const double **xyz, const double **t );
PETSC_EXTERN int  EG_getTessFace( const ego tess, int fIndex, int *len,
                                  const double **xyz, const double **uv,
                                  const int **ptype, const int **pindex,
                                  int *ntri, const int **tris,
                                  const int **tric );
PETSC_EXTERN int  EG_getTessLoops( const ego tess, int fIndex, int *nloop,
                                   const int **lIndices );
PETSC_EXTERN int  EG_getTessQuads( const ego tess, int *nquad,
                                   int **fIndices );
PETSC_EXTERN int  EG_makeQuads( ego tess, double *params, int fIndex );
PETSC_EXTERN int  EG_getQuads( const ego tess, int fIndex, int *len,
                                  const double **xyz, const double **uv,
                                  const int **ptype, const int **pindex,
                                  int *npatch );
PETSC_EXTERN int  EG_getPatch( const ego tess, int fIndex, int patch,
                               int *nu, int *nv, const int **ipts,
                               const int **bounds );
PETSC_EXTERN int  EG_quadTess( const ego tess, ego *quadTess );

PETSC_EXTERN int  EG_insertEdgeVerts( ego tess, int eIndex, int vIndex,
                                      int npts, double *t );
PETSC_EXTERN int  EG_deleteEdgeVert( ego tess, int eIndex, int vIndex,
                                     int dir );
PETSC_EXTERN int  EG_moveEdgeVert( ego tess, int eIndex, int vIndex,
                                   double t );

PETSC_EXTERN int  EG_openTessBody( ego tess );
PETSC_EXTERN int  EG_initTessBody( ego object, ego *tess );
PETSC_EXTERN int  EG_statusTessBody( ego tess, ego *body, int *state, int *np );
PETSC_EXTERN int  EG_setTessEdge( ego tess, int eIndex, int len,
                                  const double *xyz, const double *t );
PETSC_EXTERN int  EG_setTessFace( ego tess, int fIndex, int len,
                                  const double *xyz, const double *uv,
                                  int ntri, const int *tris );
PETSC_EXTERN int  EG_localToGlobal( const ego tess, int index, int local,
                                    int *global );
PETSC_EXTERN int  EG_getGlobal( const ego tess, int global, int *ptype,
                                int *pindex, /*@null@*/ double *xyz );

/* high level functions */

PETSC_EXTERN int  EG_fuseSheets( const ego src, const ego tool, ego *sheet );
PETSC_EXTERN int  EG_solidBoolean( const ego src, const ego tool, int oper,
                                   ego *model );
PETSC_EXTERN int  EG_intersection( const ego src, const ego tool, int *nedge,
                                   /*@null@*/ ego **facEdg, ego *model );
PETSC_EXTERN int  EG_imprintBody( const ego src, int nedge, const ego *facEdg,
                                  ego *result );
PETSC_EXTERN int  EG_filletBody( const ego src, int nedge, const ego *edges,
                                 double radius,
                                 ego *result, /*@null@*/ int **facemap );
PETSC_EXTERN int  EG_chamferBody( const ego src, int nedge, const ego *edges,
                                  const ego *faces, double dis1, double dis2,
                                  ego *result, /*@null@*/ int **facemap );
PETSC_EXTERN int  EG_hollowBody( const ego src, int nface,
                                 /*@null@*/ const ego *faces, double offset,
                                 int join, ego *result,
                                 /*@null@*/ int **facemap );
PETSC_EXTERN int  EG_extrude( const ego src, double dist, const double *dir,
                                    ego *result );
PETSC_EXTERN int  EG_rotate( const ego src, double angle, const double *axis,
                                   ego *result );
PETSC_EXTERN int  EG_sweep( const ego src, const ego spine, int mode,
                                  ego *result );
PETSC_EXTERN int  EG_loft( int nsec, const ego *secs, int opt, ego *result );
PETSC_EXTERN int  EG_blend( int nsec, const ego *secs, /*@null@*/ double *rc1,
                            /*@null@*/ double *rcN, ego *result );
PETSC_EXTERN int  EG_ruled( int nsec, const ego *secs, ego *result );

/* other functions */

PETSC_EXTERN int  EG_inTriExact( double *t1, double *t2, double *t3, double *p,
                                 double *w );

#endif
