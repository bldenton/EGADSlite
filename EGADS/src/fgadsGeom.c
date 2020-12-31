/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             FORTRAN Bindings for Geometrical Functions
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


  extern int EG_getGeometry(const egObject *geom, int *oclass, int *type, 
                                  egObject **ref, int **ivec, double **rvec);
  extern int EG_makeGeometry(egObject *context, int oclass, int mtype,
                             /*@null@*/ ego refGeom,
                             /*@null@*/ const int *ivec, 
                             const double *rvec, egObject **geom);
  extern int EG_getRange(const egObject *geom, double *range, int *pflag);
  extern int EG_evaluate(const egObject *geom, const double *param, 
                         double *results);
  extern int EG_invEvaluate(const egObject *geom, double *xyz, double *param, 
                            double *results);
  extern int EG_invEvaluateGuess(const egObject *geom, double *xyz, 
                                 double *param, double *results);
  extern int EG_arcLength(const egObject *geom, double t1, double t2,
                          double *alen);
  extern int EG_curvature(const egObject *geom, const double *param,
                          double *crva);
  extern int EG_approximate(egObject *context, int maxdeg, double tol,
                            const int *sizes, const double *xyzs, 
                            egObject **bspline);
  extern int EG_otherCurve(const egObject *surface, const egObject *curve, 
                           double tol, egObject **newcrv);
  extern int EG_isoCline(const egObject *surface, int iUV, double value, 
                               egObject **newcrv);
  extern int EG_isSame(const egObject *geom1, const egObject *geom2);
  extern int EG_convertToBSpline(egObject *geom, egObject **bspline); 
  extern int EG_addKnots(const egObject *object, int nU, /*@null@*/ double *Us,
                                                 int nV, /*@null@*/ double *Vs,
                         egObject **result);
  extern int EG_adjustCPs(const egObject *body, const egObject *face,
                          double *CPs, egObject **newBody, egObject **newFace);

int
#ifdef WIN32
IG_GETGEOMETRY (INT8 *obj, int *oclass, int *mtype, INT8 *igeom,  
                int **ivec, double **rvec)
#else
ig_getgeometry_(INT8 *obj, int *oclass, int *mtype, INT8 *igeom,  
                int **ivec, double **rvec)
#endif
{
  int      stat;
  egObject *object, *geom;

  *ivec   = NULL;
  *rvec   = NULL;
  *oclass = *mtype = 0;
  *igeom  = 0;
  object  = (egObject *) *obj;
  stat    = EG_getGeometry(object, oclass, mtype, &geom, ivec, rvec); 
  if (stat == EGADS_SUCCESS) *igeom = (INT8) geom;
  return stat;
}


int
#ifdef WIN32
IG_MAKEGEOMETRY (INT8 *cntx, int *oclass, int *mtype, INT8 *rgeom,  
                const int *ivec, const double *rvec,  INT8 *igeom)
#else
ig_makegeometry_(INT8 *cntx, int *oclass, int *mtype, INT8 *rgeom,  
                const int *ivec, const double *rvec,  INT8 *igeom)
#endif
{
  int      stat;
  egObject *context, *geom, *ref;

  *igeom  = 0;
  context = (egObject *) *cntx;
  ref     = (egObject *) *rgeom;
  stat    = EG_makeGeometry(context, *oclass, *mtype, ref, ivec, rvec,
                            &geom); 
  if (stat == EGADS_SUCCESS) *igeom = (INT8) geom;
  return stat;
}


int
#ifdef WIN32
IG_GETRANGE (INT8 *obj, double *range, int *pflag)
#else
ig_getrange_(INT8 *obj, double *range, int *pflag)
#endif
{
  egObject *object;

  *pflag = 0;
  object = (egObject *) *obj;
  return EG_getRange(object, range, pflag);
}


int
#ifdef WIN32
IG_EVALUATE (INT8 *obj, double *param, double *results)
#else
ig_evaluate_(INT8 *obj, double *param, double *results)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_evaluate(object, param, results);
}


int
#ifdef WIN32
IG_INVEVALUATE (INT8 *obj, double *xyz, double *param, double *results)
#else
ig_invevaluate_(INT8 *obj, double *xyz, double *param, double *results)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_invEvaluate(object, xyz, param, results);
}


int
#ifdef WIN32
IG_INVEVALUATEGUESS (INT8 *obj, double *xyz, double *param, double *results)
#else
ig_invevaluateguess_(INT8 *obj, double *xyz, double *param, double *results)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_invEvaluateGuess(object, xyz, param, results);
}


int
#ifdef WIN32
IG_ARCLENGTH (INT8 *obj, double *t1, double *t2, double *alen)
#else
ig_arclength_(INT8 *obj, double *t1, double *t2, double *alen)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_arcLength(object, *t1, *t2, alen);
}


int
#ifdef WIN32
IG_CURVATURE (INT8 *obj, double *param, double *crva)
#else
ig_curvature_(INT8 *obj, double *param, double *crva)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_curvature(object, param, crva);
}


int
#ifdef WIN32
IG_APPROXIMATE (INT8 *cntx, int *maxdeg, double *tol, const int *size, 
                const double *xyzs, INT8 *igeom)
#else
ig_approximate_(INT8 *cntx, int *maxdeg, double *tol, const int *size, 
                const double *xyzs, INT8 *igeom)
#endif
{
  int      stat;
  egObject *context, *geom;

  *igeom  = 0;
  context = (egObject *) *cntx;
  stat    = EG_approximate(context, *maxdeg, *tol, size, xyzs, &geom); 
  if (stat == EGADS_SUCCESS) *igeom = (INT8) geom;
  return stat;
}


int
#ifdef WIN32
IG_OTHERCURVE (INT8 *isurf, INT8 *icrv, double *tol, INT8 *igeom)
#else
ig_othercurve_(INT8 *isurf, INT8 *icrv, double *tol, INT8 *igeom)
#endif
{
  int      stat;
  egObject *surf, *curv, *geom;

  *igeom = 0;
  surf   = (egObject *) *isurf;
  curv   = (egObject *) *icrv;
  stat   = EG_otherCurve(surf, curv, *tol, &geom); 
  if (stat == EGADS_SUCCESS) *igeom = (INT8) geom;
  return stat;
}


int
#ifdef WIN32
IG_ISOCLINE (INT8 *isurf, int *iUV, double *value, INT8 *igeom)
#else
ig_isocline_(INT8 *isurf, int *iUV, double *value, INT8 *igeom)
#endif
{
  int      stat;
  egObject *surf, *geom;

  *igeom = 0;
  surf   = (egObject *) *isurf;
  stat   = EG_isoCline(surf, *iUV, *value, &geom); 
  if (stat == EGADS_SUCCESS) *igeom = (INT8) geom;
  return stat;
}


int
#ifdef WIN32
IG_ISSAME (INT8 *igeom1, INT8 *igeom2)
#else
ig_issame_(INT8 *igeom1, INT8 *igeom2)
#endif
{
  egObject *geom1, *geom2;
  
  geom1 = (egObject *) *igeom1;
  geom2 = (egObject *) *igeom2;
  return EG_isSame(geom1, geom2);
}


int
#ifdef WIN32
IG_CONVERTTOBSPLINE (INT8 *iobj, INT8 *igeom)
#else
ig_converttobspline_(INT8 *iobj, INT8 *igeom)
#endif
{
  int      stat;
  egObject *obj, *geom;

  *igeom = 0;
  obj    = (egObject *) *iobj;
  stat   = EG_convertToBSpline(obj, &geom); 
  if (stat == EGADS_SUCCESS) *igeom = (INT8) geom;
  return stat;
}


int
#ifdef WIN32
IG_ADDKNOTS (INT8 *iobj, int *nu, double *Us, int *nv, double *Vs, INT8 *igeom)
#else
ig_addknots_(INT8 *iobj, int *nu, double *Us, int *nv, double *Vs, INT8 *igeom)
#endif
{
  int      stat;
  egObject *obj, *geom;
  
  *igeom = 0;
  obj    = (egObject *) *iobj;
  stat   = EG_addKnots(obj, *nu, Us, *nv, Vs, &geom);
  if (stat == EGADS_SUCCESS) *igeom = (INT8) geom;
  return stat;
}


int
#ifdef WIN32
IG_ADJUSTCPS (INT8 *ibody, INT8 *iface, double *CPs, INT8 *nbody, INT8 *nface)
#else
ig_adjustcps_(INT8 *ibody, INT8 *iface, double *CPs, INT8 *nbody, INT8 *nface)
#endif
{
  int      stat;
  egObject *body, *face, *newBody, *newFace;
  
  *nbody = *nface = 0;
  body   = (egObject *) *ibody;
  face   = (egObject *) *iface;
  stat   = EG_adjustCPs(body, face, CPs, &newBody, &newFace);
  if (stat == EGADS_SUCCESS) {
    *nbody = (INT8) newBody;
    *nface = (INT8) newFace;
  }
  return stat;
}
