#ifndef EGADS_H
#define EGADS_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Function Prototypes
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "egadsTypes.h"

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif
#ifdef __DEVICE__
#undef __DEVICE__
#endif

#ifdef __CUDACC__
# define __HOST_AND_DEVICE__ __host__ __device__
# define __DEVICE__ __device__
#else
# define __HOST_AND_DEVICE__
# define __DEVICE__
#endif

#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__ __HOST_AND_DEVICE__
#else
#define __ProtoExt__ __HOST_AND_DEVICE__ extern
#endif

/* memory functions */

__ProtoExt__ /*@null@*/ /*@out@*/ /*@only@*/
             void *EG_alloc( size_t nbytes );
__ProtoExt__ /*@null@*/ /*@only@*/
             void *EG_calloc( size_t nele, size_t size );
__ProtoExt__ /*@null@*/ /*@only@*/
             void *EG_reall( /*@null@*/ /*@only@*/ /*@returned@*/ void *ptr,
                             size_t nbytes );
__ProtoExt__ /*@null@*/ /*@only@*/
             char *EG_strdup( /*@null@*/ const char *str );
__ProtoExt__ void  EG_free( /*@null@*/ /*@only@*/ void *pointer );

/* base-level object functions */

__ProtoExt__ void EG_revision( int *major, int *minor, const char **OCCrev );
__ProtoExt__ int  EG_open( ego *context );
__ProtoExt__ int  EG_loadModel( ego context, int bflg, const char *name, 
                                ego *model );
__ProtoExt__ int  EG_saveModel( const ego model, const char *name );
__ProtoExt__ int  EG_exportModel( ego model, size_t *nbytes, char **stream );
__ProtoExt__ int  EG_importModel( ego context, const size_t nbytes,
                                  const char *stream, ego *model );
__ProtoExt__ int  EG_deleteObject( ego object );
__ProtoExt__ int  EG_makeTransform( ego context, const double *xform, 
                                    ego *oform );
__ProtoExt__ int  EG_getTransformation( const ego oform, double *xform );
__ProtoExt__ int  EG_getContext( ego object, ego *context );
__ProtoExt__ int  EG_setOutLevel( ego context, int outLevel );
__ProtoExt__ int  EG_updateThread( ego context );
__ProtoExt__ int  EG_getInfo( const ego object, int *oclass, int *mtype, 
                              ego *topObj, ego *prev, ego *next );
__ProtoExt__ int  EG_copyObject( const ego object, /*@null@*/ void *oform,
                                 ego *copy );
__ProtoExt__ int  EG_flipObject( const ego object, ego *flippedCopy );
__ProtoExt__ int  EG_close( ego context );
__ProtoExt__ int  EG_setUserPointer( ego context, void *ptr );
__ProtoExt__ int  EG_getUserPointer( const ego context, void **ptr );

/* attribute functions */

__ProtoExt__ int  EG_setFullAttrs( ego context, int flag );
__ProtoExt__ int  EG_attributeAdd( ego obj, const char *name, int type, int len,
                                  /*@null@*/ const int    *ints, 
                                  /*@null@*/ const double *reals,
                                  /*@null@*/ const char   *str );
__ProtoExt__ int  EG_attributeDel( ego object, /*@null@*/ const char *name );
__ProtoExt__ int  EG_attributeNum( const ego obj, int *num );
__ProtoExt__ int  EG_attributeGet( const ego obj, int index, const char **name,
                                   int *atype, int *len, 
                                   /*@null@*/ const int    **ints,
                                   /*@null@*/ const double **reals, 
                                   /*@null@*/ const char   **str );
__ProtoExt__ int  EG_attributeRet( const ego obj, const char *name, int *atype, 
                                   int *len, /*@null@*/ const int    **ints,
                                             /*@null@*/ const double **reals, 
                                             /*@null@*/ const char   **str );
__ProtoExt__ int  EG_attributeDup( const ego src, ego dst );
__ProtoExt__ int  EG_attributeAddSeq( ego obj, const char *name, int type,
                                      int len, /*@null@*/ const int    *ints,
                                               /*@null@*/ const double *reals,
                                               /*@null@*/ const char   *str );
__ProtoExt__ int  EG_attributeNumSeq( const ego obj, const char *name,
                                      int *num );
__ProtoExt__ int  EG_attributeRetSeq( const ego obj, const char *name,
                                      int index, int *atype, int *len,
                                      /*@null@*/ const int    **ints,
                                      /*@null@*/ const double **reals,
                                      /*@null@*/ const char   **str );

/* geometry functions */

__ProtoExt__ int  EG_getGeometry( const ego geom, int *oclass, int *mtype,
                                  ego *refGeom, /*@null@*/ int **ivec,
                                  /*@null@*/ double **rvec );
__ProtoExt__ int  EG_makeGeometry( ego context, int oclass, int mtype,
                                  /*@null@*/ ego refGeom,
                                  /*@null@*/ const int *ivec,
                                  const double *rvec, ego *geom );
__ProtoExt__ int  EG_getRange( const ego geom, double *range, int *periodic );
__ProtoExt__ int  EG_evaluate( const ego geom, /*@null@*/ const double *param, 
                               double *results );
__ProtoExt__ int  EG_invEvaluate( const ego geom, double *xyz, double *param,
                                  double *results );
__ProtoExt__ int  EG_invEvaluateGuess( const ego geom, double *xyz, 
                                       double *param, double *results );
__ProtoExt__ int  EG_arcLength( const ego geom, double t1, double t2,
                                double *alen );
__ProtoExt__ int  EG_curvature( const ego geom, const double *param, 
                                double *results );
__ProtoExt__ int  EG_approximate( ego context, int maxdeg, double tol,
                                  const int *sizes, const double *xyzs,
                                  ego *bspline );
__ProtoExt__ int  EG_fitTriangles( ego context, int npts, double *xyzs,
                                   int ntris, const int *tris,
                                   /*@null@*/ const int *tric, double tol,
                                   ego *bspline );
__ProtoExt__ int  EG_otherCurve( const ego surface, const ego curve,
                                 double tol, ego *newcurve );
__ProtoExt__ int  EG_isSame( const ego geom1, const ego geom2 );
__ProtoExt__ int  EG_isoCline( const ego surface, int UV, double val,
                               ego *newcurve );
__ProtoExt__ int  EG_convertToBSpline( ego geom, ego *bspline );
__ProtoExt__ int  EG_convertToBSplineRange( ego geom, const double *range,
                                            ego *bspline );
__ProtoExt__ int  EG_skinning( int nCurves, ego *curves, int skinning_degree,
                               ego *surface );

/* topology functions */

__ProtoExt__ int  EG_tolerance( const ego topo, double *tol );
__ProtoExt__ int  EG_getTolerance( const ego topo, double *tol );
__ProtoExt__ int  EG_setTolerance(       ego topo, double  tol );
__ProtoExt__ int  EG_getTopology( const ego topo, ego *geom, int *oclass, 
                                  int *type, /*@null@*/ double *limits, 
                                  int *nChildren, ego **children, int **sense );
__ProtoExt__ int  EG_makeTopology( ego context, /*@null@*/ ego geom, int oclass,
                                   int mtype, /*@null@*/ double *limits,
                                   int nChildren, /*@null@*/ ego *children,
                                   /*@null@*/ int *senses, ego *topo );
__ProtoExt__ int  EG_makeLoop( int nedge, ego *edges, /*@null@*/ ego geom,
                               double toler, ego *result );
__ProtoExt__ int  EG_getArea( ego object, /*@null@*/ const double *limits,
                              double *area );
__ProtoExt__ int  EG_makeFace( ego object, int mtype,
                               /*@null@*/ const double *limits, ego *face );
__ProtoExt__ int  EG_getBodyTopos( const ego body, /*@null@*/ ego src,
                                   int oclass, int *ntopo,
                                   /*@null@*/ ego **topos );
__ProtoExt__ int  EG_indexBodyTopo( const ego body, const ego src );
__ProtoExt__ int  EG_objectBodyTopo( const ego body, int oclass, int index,
                                     ego *obj );
__ProtoExt__ int  EG_inTopology( const ego topo, const double *xyz );
__ProtoExt__ int  EG_inFace( const ego face, const double *uv );
__ProtoExt__ int  EG_getEdgeUV( const ego face, const ego edge, int sense,
                                double t, double *UV );
__ProtoExt__ int  EG_getEdgeUVs( const ego face, const ego edge, int sense,
                                 int nts, const double *ts, double *UV );
__ProtoExt__ int  EG_getPCurve( const ego face, const ego edge, int sense,
                                int *mtype, int **ivec, double **rvec );
__ProtoExt__ int  EG_getWindingAngle( ego edge, double t, double *angle );
__ProtoExt__ int  EG_getBody( const ego object, ego *body );
__ProtoExt__ int  EG_makeSolidBody( ego context, int stype, const double *rvec,
                                    ego *body );
__ProtoExt__ int  EG_getBoundingBox( const ego topo, double *bbox );
__ProtoExt__ int  EG_getMassProperties( const ego topo, 
                                        /*@null@*/ double *result );
__ProtoExt__ int  EG_isEquivalent( const ego topo1, const ego topo2 );
__ProtoExt__ int  EG_sewFaces( int nobj, const ego *objs, double toler,
                               int flag, ego *result );
__ProtoExt__ int  EG_replaceFaces( const ego body, int nobj, ego *objs,
                                   ego *result );
__ProtoExt__ int  EG_mapBody( const ego sBody,   const ego dBody,
                              const char *fAttr, ego *mapBody );
__ProtoExt__ int  EG_matchBodyEdges( const ego body1, const ego body2,
                                     double toler, int *nmatch, int **match );
__ProtoExt__ int  EG_matchBodyFaces( const ego body1, const ego body2,
                                     double toler, int *nmatch, int **match );

/* tessellation functions */

__ProtoExt__ int  EG_setTessParam( ego context, int iparam, double value,
                                   double *oldvalue );
__ProtoExt__ int  EG_makeTessGeom( ego obj, double *params, int *sizes, 
                                   ego *tess );
__ProtoExt__ int  EG_getTessGeom( const ego tess, int *sizes, double **xyz );

__ProtoExt__ int  EG_makeTessBody( ego object, double *params, ego *tess );
__ProtoExt__ int  EG_remakeTess( ego tess, int nobj, ego *objs, 
                                 double *params );
__ProtoExt__ int  EG_finishTess( ego tess, double *params );
__ProtoExt__ int  EG_mapTessBody( ego tess, ego body, ego *mapTess );
__ProtoExt__ int  EG_locateTessBody( const ego tess, int npt, const int *ifaces,
                                     const double *uv, /*@null@*/ int *itri, 
                                     double *results );

__ProtoExt__ int  EG_getTessEdge( const ego tess, int eIndex, int *len, 
                                  const double **xyz, const double **t );
__ProtoExt__ int  EG_getTessFace( const ego tess, int fIndex, int *len, 
                                  const double **xyz, const double **uv, 
                                  const int **ptype, const int **pindex, 
                                  int *ntri, const int **tris, 
                                  const int **tric );
__ProtoExt__ int  EG_getTessLoops( const ego tess, int fIndex, int *nloop,
                                   const int **lIndices );
__ProtoExt__ int  EG_getTessQuads( const ego tess, int *nquad,
                                   int **fIndices );
__ProtoExt__ int  EG_makeQuads( ego tess, double *params, int fIndex );
__ProtoExt__ int  EG_getQuads( const ego tess, int fIndex, int *len, 
                                  const double **xyz, const double **uv, 
                                  const int **ptype, const int **pindex, 
                                  int *npatch );
__ProtoExt__ int  EG_getPatch( const ego tess, int fIndex, int patch, 
                               int *nu, int *nv, const int **ipts, 
                               const int **bounds );
__ProtoExt__ int  EG_quadTess( const ego tess, ego *quadTess );
                               
__ProtoExt__ int  EG_insertEdgeVerts( ego tess, int eIndex, int vIndex, 
                                      int npts, double *t );
__ProtoExt__ int  EG_deleteEdgeVert( ego tess, int eIndex, int vIndex, 
                                     int dir );
__ProtoExt__ int  EG_moveEdgeVert( ego tess, int eIndex, int vIndex, 
                                   double t );
 
__ProtoExt__ int  EG_openTessBody( ego tess );
__ProtoExt__ int  EG_initTessBody( ego object, ego *tess );
__ProtoExt__ int  EG_statusTessBody( ego tess, ego *body, int *state, int *np );
__ProtoExt__ int  EG_setTessEdge( ego tess, int eIndex, int len,
                                  const double *xyz, const double *t );
__ProtoExt__ int  EG_setTessFace( ego tess, int fIndex, int len,
                                  const double *xyz, const double *uv,
                                  int ntri, const int *tris );
__ProtoExt__ int  EG_localToGlobal( const ego tess, int index, int local,
                                    int *global );
__ProtoExt__ int  EG_getGlobal( const ego tess, int global, int *ptype,
                                int *pindex, /*@null@*/ double *xyz );
__ProtoExt__ int  EG_saveTess( ego tess, const char *name );
__ProtoExt__ int  EG_loadTess( ego body, const char *name, ego *tess );
  
__ProtoExt__ int  EG_tessMassProps( const ego tess, double *props );

/* top down build functions */

__ProtoExt__ int  EG_fuseSheets( const ego src, const ego tool, ego *sheet );
__ProtoExt__ int  EG_generalBoolean( ego src, ego tool, int oper, double tol,
                                     ego *model );
__ProtoExt__ int  EG_solidBoolean( const ego src, const ego tool, int oper, 
                                   ego *model );
__ProtoExt__ int  EG_intersection( const ego src, const ego tool, int *nedge, 
                                   /*@null@*/ ego **facEdg, ego *model );
__ProtoExt__ int  EG_imprintBody( const ego src, int nedge, const ego *facEdg, 
                                  ego *result );
__ProtoExt__ int  EG_filletBody( const ego src, int nedge, const ego *edges, 
                                 double radius,
                                 ego *result, /*@null@*/ int **facemap );
__ProtoExt__ int  EG_chamferBody( const ego src, int nedge, const ego *edges, 
                                  const ego *faces, double dis1, double dis2, 
                                  ego *result, /*@null@*/ int **facemap );
__ProtoExt__ int  EG_hollowBody( const ego src, int nface,
                                 /*@null@*/ const ego *faces, double offset,
                                 int join, ego *result,
                                 /*@null@*/ int **facemap );
__ProtoExt__ int  EG_extrude( const ego src, double dist, const double *dir, 
                                    ego *result );
__ProtoExt__ int  EG_rotate( const ego src, double angle, const double *axis, 
                                   ego *result );
__ProtoExt__ int  EG_sweep( const ego src, const ego spine, int mode,
                                  ego *result );
__ProtoExt__ int  EG_loft( int nsec, const ego *secs, int opt, ego *result );
__ProtoExt__ int  EG_blend( int nsec, const ego *secs, /*@null@*/ double *rc1,
                            /*@null@*/ double *rcN, ego *result );
__ProtoExt__ int  EG_ruled( int nsec, const ego *secs, ego *result );

/* effective topology functions */

__ProtoExt__ int  EG_initEBody( ego tess, double angle, ego *ebody );
__ProtoExt__ int  EG_finishEBody( ego EBody );
__ProtoExt__ int  EG_makeEFace( ego EBody, int nFace, ego *Faces, ego *EFace );
__ProtoExt__ int  EG_makeAttrEFaces( ego EBody, const char *attrName, int *nEFace,
                                     /*@null@*/ ego **EFaces );
__ProtoExt__ int  EG_effectiveMap( ego EObject, double *eparam, ego *Object,
                                   double *param );
__ProtoExt__ int  EG_effectiveEdgeList( ego EEdge, int *nedge, ego **edges,
                                        int **senses, double **tstart );
__ProtoExt__ int  EG_effectiveTri( ego EObj, double *uv, int *fIndex, int *itri,
                                   double *w );

/* other functions */

__ProtoExt__ int  EG_inTriExact( double *t1, double *t2, double *t3, double *p,
                                 double *w );
#ifdef __cplusplus
}
#endif

#endif
