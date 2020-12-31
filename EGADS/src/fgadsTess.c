/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             FORTRAN Bindings for Tessellation Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>

#define INT8 unsigned long long

#include "egadsTypes.h"
#include "egadsInternals.h"


  extern int EG_makeTessGeom(egObject *obj, double *params, int *sizes,
                             egObject **tess);
  extern int EG_getTessGeom(const egObject *tess, int *sizes, double **xyz);

  extern int EG_makeTessBody(egObject *object, double *params, egObject **tess);
  extern int EG_remakeTess(egObject *tess, int nobj, /*@null@*/ egObject **objs,
                           double *params);
  extern int EG_finishTess(egObject *tess, double *params);
  extern int EG_mapTessBody(egObject *tess, egObject *body, egObject **mapTess);
  extern int EG_locateTessBody(const egObject *tess, int npt, const int *ifaces,
                               const double *uvs, /*@null@*/ int *itris, 
                               double *weights);

  extern int EG_getTessEdge(const egObject *tess, int index, int *len,
                            const double **xyz, const double **t);
  extern int EG_getTessFace(const egObject *tess, int index, int *len, 
                            const double **xyz, const double **uv, 
                            const int **ptype, const int **pindex, 
                            int *ntri, const int **tris, const int **tric);
  extern int EG_getTessLoops(const egObject *tess, int index, int *nloop,
                             const int **lIndices);

  extern int EG_getTessQuads(const egObject *tess, int *nquad, int **fIndices);
  extern int EG_makeQuads(egObject *tess, double *params, int fIndex);
  extern int EG_getQuads(const egObject *tess, int fIndex, int *len, 
                         const double **xyz, const double **uv, 
                         const int **ptype, const int **pindex, int *npatch);
  extern int EG_getPatch(const egObject *tess, int fIndex, int patch, int *nu, 
                         int *nv, const int **ipts, const int **bounds);
  extern int EG_quadTess(const egObject *tess, egObject **quadTess);

  extern int EG_insertEdgeVerts(egObject *tess, int eIndex, int vIndex, 
                                int npts, double *t);
  extern int EG_deleteEdgeVert(egObject *tess, int eIndex, int vIndex, int dir);
  extern int EG_moveEdgeVert(egObject *tess, int eIndex, int vIndex, double t);

  extern int EG_openTessBody(egObject *tess);
  extern int EG_initTessBody(egObject *object, egObject **tess);
  extern int EG_statusTessBody(egObject *tess, egObject **body, int *state,
                               int *npts);
  extern int EG_setTessEdge(const egObject *tess, int index, int len,
                            const double *xyz, const double *t);
  extern int EG_setTessFace(const egObject *tess, int index, int len,
                            const double *xyz, const double *uv,
                            int ntri, const int *tris);
  extern int EG_localToGlobal(const egObject *tess, int index, int local,
                              int *global);
  extern int EG_getGlobal(const egObject *tess, int global, int *ptype,
                          int *pindex, /*@null@*/ double *xyz);

  extern int  EG_tessMassProps(const ego tess, double *props);
  extern int  EG_tessMassProps_dot(const ego tess, double *xyz_dot,
                                   double *props, double *props_dot);


int
#ifdef WIN32
IG_MAKETESSGEOM (INT8 *obj, double *params, int *sizes, INT8 *itess)  
#else
ig_maketessgeom_(INT8 *obj, double *params, int *sizes, INT8 *itess)
#endif
{
  int      stat;
  egObject *object, *tess;

  *itess = 0;
  object = (egObject *) *obj;
  stat   = EG_makeTessGeom(object, params, sizes, &tess); 
  if (stat == EGADS_SUCCESS) *itess = (INT8) tess;
  return stat;
}


int
#ifdef WIN32
IG_GETTESSGEOM (INT8 *obj, int *sizes, double **xyz)
#else
ig_gettessgeom_(INT8 *obj, int *sizes, double **xyz)
#endif
{
  egObject *object;

  sizes[0] = 0;
  *xyz     = NULL;
  object   = (egObject *) *obj;
  return EG_getTessGeom(object, sizes, xyz);
}


int
#ifdef WIN32
IG_MAKETESSBODY (INT8 *obj, double *params, INT8 *itess)  
#else
ig_maketessbody_(INT8 *obj, double *params, INT8 *itess)
#endif
{
  int      stat;
  egObject *object, *tess;

  *itess = 0;
  object = (egObject *) *obj;
  stat   = EG_makeTessBody(object, params, &tess); 
  if (stat == EGADS_SUCCESS) *itess = (INT8) tess;
  return stat;
}


int
#ifdef WIN32
IG_REMAKETESS (INT8 *itess, int *nobj, INT8 *objs, double *params)  
#else
ig_remaketess_(INT8 *itess, int *nobj, INT8 *objs, double *params)
#endif
{
  int      i, stat;
  egObject *tess, **objects = NULL;

  tess = (egObject *) *itess;
  if (*nobj > 0) {
    objects = (egObject **) EG_alloc(*nobj*sizeof(egObject *));
    if (objects == NULL) return EGADS_MALLOC;
    for (i = 0; i < *nobj; i++)
      objects[i] = (egObject *) objs[i];
  }
  stat = EG_remakeTess(tess, *nobj, objects, params); 
  if (objects != NULL) EG_free(objects);
  return stat;
}


int
#ifdef WIN32
IG_FINISHTESS (INT8 *itess, double *params)
#else
ig_finishtess_(INT8 *itess, double *params)
#endif
{
  egObject *tess;
  
  tess = (egObject *) *itess;
  return EG_finishTess(tess, params);
}


int
#ifdef WIN32
IG_MAPTESSBODY (INT8 *tobj, INT8 *obj, INT8 *itess)
#else
ig_maptessbody_(INT8 *tobj, INT8 *obj, INT8 *itess)
#endif
{
  int      stat;
  egObject *object, *tobjct, *tess;
  
  *itess = 0;
  tobjct = (egObject *) *tobj;
  object = (egObject *) *obj;
  stat   = EG_mapTessBody(tobjct, object, &tess);
  if (stat == EGADS_SUCCESS) *itess = (INT8) tess;
  return stat;
}


int
#ifdef WIN32
IG_LOCATETESSBODY (INT8 *obj, int *npt, const int *iface, const double *uvs,
                   int *itri, double *results)
#else
ig_locatetessbody_(INT8 *obj, int *npt, const int *iface, const double *uvs,
                   int *itri, double *results)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_locateTessBody(object, *npt, iface, uvs, itri, results);
}


int
#ifdef WIN32
IG_GETTESSEDGE (INT8 *obj, int *index, int *len, const double **xyz,
                                                 const double **t)
#else
ig_gettessedge_(INT8 *obj, int *index, int *len, const double **xyz, 
                                                 const double **t)
#endif
{
  egObject *object;

  *len   = 0;
  *xyz   = NULL;
  *t     = NULL;
  object = (egObject *) *obj;
  return EG_getTessEdge(object, *index, len, xyz, t);
}


int
#ifdef WIN32
IG_GETTESSFACE (INT8 *obj, int *index, int *len, const double **xyz, 
                const double **uv, const int **ptype, const int **pindex,
                int *ntri, const int **tris, const int **tric)
#else
ig_gettessface_(INT8 *obj, int *index, int *len, const double **xyz, 
                const double **uv, const int **ptype, const int **pindex,
                int *ntri, const int **tris, const int **tric)
#endif
{
  egObject *object;

  *len    = *ntri = 0;
  *xyz    = NULL;
  *uv     = NULL;
  *ptype  = NULL;
  *pindex = NULL;
  *tris   = NULL;
  *tric   = NULL;
  object  = (egObject *) *obj;
  return EG_getTessFace(object, *index, len, xyz, uv, ptype, pindex, ntri,
                        tris, tric);
}


int
#ifdef WIN32
IG_GETTESSLOOPS (INT8 *obj, int *index, int *nloop, const int **lIndices)
#else
ig_gettessloops_(INT8 *obj, int *index, int *nloop, const int **lIndices)
#endif
{
  egObject *object;
  
  *nloop    = 0;
  *lIndices = NULL;
  object    = (egObject *) *obj;
  return EG_getTessLoops(object, *index, nloop, lIndices);
}


int
#ifdef WIN32
IG_GETTESSQUADS (INT8 *obj, int *nquad, int **fIndices)
#else
ig_gettessquads_(INT8 *obj, int *nquad, int **fIndices)
#endif
{
  egObject *object;

  *nquad    = 0;
  *fIndices = NULL;
  object    = (egObject *) *obj;
  return EG_getTessQuads(object, nquad, fIndices);
}


int
#ifdef WIN32
IG_MAKEQUADS (INT8 *obj, double *parms, int *fIndex)
#else
ig_makequads_(INT8 *obj, double *parms, int *fIndex)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_makeQuads(object, parms, *fIndex);
}


int
#ifdef WIN32
IG_GETQUADS (INT8 *obj, int *index, int *len, const double **xyz, 
            const double **uv, const int **ptype, const int **pindex,
            int *npatch)
#else
ig_getquads_(INT8 *obj, int *index, int *len, const double **xyz, 
            const double **uv, const int **ptype, const int **pindex,
            int *npatch)
#endif
{
  egObject *object;

  *len    = *npatch = 0;
  *xyz    = NULL;
  *uv     = NULL;
  *ptype  = NULL;
  *pindex = NULL;
  object  = (egObject *) *obj;
  return EG_getQuads(object, *index, len, xyz, uv, ptype, pindex, npatch);
}


int
#ifdef WIN32
IG_GETPATCH (INT8 *obj, int *index, int *patch, int *nu, int *nv, 
             const int **ipts, const int **bounds)
#else
ig_getpatch_(INT8 *obj, int *index, int *patch, int *nu, int *nv, 
             const int **ipts, const int **bounds)
#endif
{
  egObject *object;

  *nu     = *nv = 0;
  *ipts   = NULL;
  *bounds = NULL;
  object  = (egObject *) *obj;
  return EG_getPatch(object, *index, *patch, nu, nv, ipts, bounds);
}


int
#ifdef WIN32
IG_QUADTESS (INT8 *tobj, INT8 *itess)
#else
ig_quadtess_(INT8 *tobj, INT8 *itess)
#endif
{
  int      stat;
  egObject *tobjct, *tess;
  
  *itess = 0;
  tobjct = (egObject *) *tobj;
  stat   = EG_quadTess(tobjct, &tess);
  if (stat == EGADS_SUCCESS) *itess = (INT8) tess;
  return stat;
}


int
#ifdef WIN32
IG_INSERTEDGEVERTS (INT8 *obj, int *index, int *vert, int *npts, double *ts)
#else
ig_insertedgeverts_(INT8 *obj, int *index, int *vert, int *npts, double *ts)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_insertEdgeVerts(object, *index, *vert, *npts, ts);
}


int
#ifdef WIN32
IG_DELETEEDGEVERT (INT8 *obj, int *index, int *vert, int *dir)
#else
ig_deleteedgevert_(INT8 *obj, int *index, int *vert, int *dir)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_deleteEdgeVert(object, *index, *vert, *dir);
}


int
#ifdef WIN32
IG_MOVEEDGEVERT (INT8 *obj, int *index, int *vert, double *t)
#else
ig_moveedgevert_(INT8 *obj, int *index, int *vert, double *t)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_moveEdgeVert(object, *index, *vert, *t);
}


int
#ifdef WIN32
IG_OPENTESSBODY (INT8 *obj)
#else
ig_opentessbody_(INT8 *obj)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_openTessBody(object);
}


int
#ifdef WIN32
IG_INITTESSBODY (INT8 *obj, INT8 *itess)
#else
ig_inittessbody_(INT8 *obj, INT8 *itess)
#endif
{
  int      stat;
  egObject *object, *tess;
  
  *itess = 0;
  object = (egObject *) *obj;
  stat   = EG_initTessBody(object, &tess);
  if (stat == EGADS_SUCCESS) *itess = (INT8) tess;
  return stat;
}


int
#ifdef WIN32
IG_STATUSTESSBODY (INT8 *obj, INT8 *ibody, int *state, int *npts)
#else
ig_statustessbody_(INT8 *obj, INT8 *ibody, int *state, int *npts)
#endif
{
  int      stat;
  egObject *object, *body;
  
  *ibody = 0;
  object = (egObject *) *obj;
  stat   = EG_statusTessBody(object, &body, state, npts);
  if (stat >= EGADS_SUCCESS) *ibody = (INT8) body;
  return stat;
}


int
#ifdef WIN32
IG_SETTESSEDGE (INT8 *obj, int *index, int *len, const double *xyz,
                const double *t)
#else
ig_settessedge_(INT8 *obj, int *index, int *len, const double *xyz,
                const double *t)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_setTessEdge(object, *index, *len, xyz, t);
}


int
#ifdef WIN32
IG_SETTESSFACE (INT8 *obj, int *index, int *len, const double *xyz,
                const double *uv, int *ntri, const int *tris)
#else
ig_settessface_(INT8 *obj, int *index, int *len, const double *xyz,
                const double *uv, int *ntri, const int *tris)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_setTessFace(object, *index, *len, xyz, uv, *ntri, tris);
}


int
#ifdef WIN32
IG_LOCALTOGLOBAL (INT8 *obj, int *index, int *local, int *global)
#else
ig_localtoglobal_(INT8 *obj, int *index, int *local, int *global)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_localToGlobal(object, *index, *local, global);
}


int
#ifdef WIN32
IG_GETGLOBAL (INT8 *obj, int *index, int *ptype, int *pindex, double *xyz)
#else
ig_getglobal_(INT8 *obj, int *index, int *ptype, int *pindex, double *xyz)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_getGlobal(object, *index, ptype, pindex, xyz);
}


int
#ifdef WIN32
IG_TESSMASSPROPS (INT8 *obj, double *props)
#else
ig_tessmassprops_(INT8 *obj, double *props)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_tessMassProps(object, props);
}


int
#ifdef WIN32
IG_TESSMASSPROPS_DOT (INT8 *obj, double *xyz_dot, double *props, double *p_dot)
#else
ig_tessmassprops_dot_(INT8 *obj, double *xyz_dot, double *props, double *p_dot)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_tessMassProps_dot(object, xyz_dot, props, p_dot);
}
