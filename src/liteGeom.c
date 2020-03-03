/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Lite Geometry Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "liteClasses.h"



#define PARAMACC         1.0e-4         /* parameter accuracy */
#define PI               3.1415926535897931159979635
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))
#define ABS(x)           ((x) < 0 ? -(x) : (x))
#define CROSS(a,b,c)       a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                           a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                           a[2] = (b[0]*c[1]) - (b[1]*c[0])


  extern int  EG_evaluateGeom( const egObject *geom, const double *param,
                               double *result );
  extern int  EG_invEvaGeomLimits( const egObject *geom,
                                   /*@null@*/ const double *limits,
                                   const double *xyz, double *param,
                                   double toler, double *result );
  extern int  EG_invEvaluateGeomGuess( const egObject *geom,
                                       /*@null@*/ const double *limits,
                                       double *xyz, double *param,
                                       double *result );

  extern int  EG_inFaceX( const egObject *face, const double *uv,
                          /*@null@*/ double *pt, /*@null@*/ double *uvx );
  extern int  EG_getEdgeUV ( const egObject *face, const egObject *edge,
                             int sense, double t, double *result );



int
EG_getGeometry(const egObject     *geom, int *oclass, int *type,
                     egObject **refGeom, /*@null@*/ int    **ivec,
                                         /*@null@*/ double **rvec)
{
  int          i, nreal, nint, *ints;
  double       *reals;
  liteGeometry *lgeom;
  
  if (ivec != NULL) *ivec = NULL;
  if (rvec != NULL) *rvec = NULL;
  *refGeom = NULL;
  *oclass  = *type = 0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass < PCURVE) || (geom->oclass > SURFACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  
  lgeom    = (liteGeometry *) geom->blind;
  *oclass  = geom->oclass;
  *type    = geom->mtype;
  *refGeom = lgeom->ref;
  if ((rvec == NULL) || (ivec == NULL)) return EGADS_SUCCESS;
  
  nreal    = nint = 0;
  if (geom->oclass == PCURVE) {
    switch (geom->mtype) {
        
      case LINE:
        nreal = 4;
        break;
        
      case CIRCLE:
        nreal = 7;
        break;
        
      case ELLIPSE:
        nreal = 8;
        break;
        
      case PARABOLA:
        nreal = 7;
        break;
        
      case HYPERBOLA:
        nreal = 8;
        break;
        
      case TRIMMED:
        nreal = 2;
        break;
        
      case BEZIER:
        nint  = 3;
        nreal = 2*lgeom->header[2];
        if ((lgeom->header[0]&2) != 0) nreal += lgeom->header[2];
        break;
        
      case BSPLINE:
        nint  = 4;
        nreal = lgeom->header[3] + 2*lgeom->header[2];
        if ((lgeom->header[0]&2) != 0) nreal += lgeom->header[2];
        break;
        
      case OFFSET:
        nreal = 1;
        break;
    }
    
  } else if (geom->oclass == CURVE) {
    switch (geom->mtype) {
        
      case LINE:
        nreal = 6;
        break;
        
      case CIRCLE:
        nreal = 10;
        break;
        
      case ELLIPSE:
        nreal = 11;
        break;
        
      case PARABOLA:
        nreal = 10;
        break;
        
      case HYPERBOLA:
        nreal = 11;
        break;
        
      case TRIMMED:
        nreal = 2;
        break;
        
      case BEZIER:
        nint  = 3;
        nreal = 3*lgeom->header[2];
        if ((lgeom->header[0]&2) != 0) nreal += lgeom->header[2];
        break;
        
      case BSPLINE:
        nint  = 4;
        nreal = lgeom->header[3] + 3*lgeom->header[2];
        if ((lgeom->header[0]&2) != 0) nreal += lgeom->header[2];
        break;
        
      case OFFSET:
        nreal = 4;
        break;
    }
    
  } else {
    /* surface */
    switch (geom->mtype) {
        
      case PLANE:
        nreal = 9;
        break;
        
      case SPHERICAL:
        nreal = 10;
        break;
        
      case CONICAL:
        nreal = 14;
        break;
        
      case CYLINDRICAL:
        nreal = 13;
        break;
        
      case TOROIDAL:
        nreal = 14;
        break;
        
      case REVOLUTION:
        nreal = 6;
        break;
        
      case EXTRUSION:
        nreal = 3;
        break;
        
      case TRIMMED:
        nreal = 4;
        break;
        
      case BEZIER:
        nint  = 5;
        nreal = 3*lgeom->header[2]*lgeom->header[4];
        if ((lgeom->header[0]&2) != 0)
          nreal += lgeom->header[2]*lgeom->header[4];
        break;
        
      case BSPLINE:
        nint  = 7;
        nreal = lgeom->header[3] + lgeom->header[6] +
                                 3*lgeom->header[2]*lgeom->header[5];
        if ((lgeom->header[0]&2) != 0)
          nreal += lgeom->header[2]*lgeom->header[5];
        break;
        
      case OFFSET:
        nreal = 1;
        break;
    }
  }
  if (nreal == 0) {
    printf(" EG_getGeometry: OCLASS = %d, MTYPE = %d not found!\n",
           geom->oclass, geom->mtype);
    return EGADS_GEOMERR;
  }

  if (nint != 0) {
    ints = (int *) EG_alloc(nint*sizeof(int));
    if (ints == NULL) return EGADS_MALLOC;
    for (i = 0; i < nint; i++) ints[i] = lgeom->header[i];
    *ivec = ints;
  }
  reals = (double *) EG_alloc(nreal*sizeof(double));
  if (reals == NULL) {
    if (*ivec != NULL) {
      EG_free(*ivec);
      *ivec = NULL;
    }
    return EGADS_MALLOC;
  }
  for (i = 0; i < nreal; i++) reals[i] = lgeom->data[i];
  *rvec = reals;
  
  return EGADS_SUCCESS;
}


int
EG_getRange(const egObject *geom, double *range, int *periodic)
{
  int          stat, mtype;
  egObject     *ref;
  liteEdge     *ledge;
  liteFace     *lface;
  liteGeometry *lgeom;
  
  *periodic = 0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  
  stat = EGADS_NOTFOUND;
  if (geom->oclass == PCURVE) {
    lgeom = (liteGeometry *) geom->blind;
    switch (geom->mtype) {
      case LINE:
        range[0]  = -2.e100;
        range[1]  =  2.e100;
        *periodic =  0;
        stat      = EGADS_SUCCESS;
        break;
        
      case CIRCLE:
        range[0]  = 0.0;
        range[1]  = 2.0*PI;
        *periodic = 1;
        stat      = EGADS_SUCCESS;
        break;
        
      case ELLIPSE:
        range[0]  = 0.0;
        range[1]  = 2.0*PI;
        *periodic = 1;
        stat      = EGADS_SUCCESS;
        break;
        
      case PARABOLA:
        range[0]  = -2.e100;
        range[1]  =  2.e100;
        *periodic =  0;
        stat      = EGADS_SUCCESS;
        break;
        
      case HYPERBOLA:
        range[0]  = -2.e100;
        range[1]  =  2.e100;
        *periodic =  0;
        stat      = EGADS_SUCCESS;
        break;
        
      case TRIMMED:
        stat = EG_getRange(lgeom->ref, range, periodic);
        if (stat == EGADS_SUCCESS) {
          range[0] = lgeom->data[0];
          range[1] = lgeom->data[1];
        }
        break;
        
      case BEZIER:
        range[0]  = 0.0;
        range[1]  = 1.0;
        *periodic = (lgeom->header[0]&4)/4;
        stat      = EGADS_SUCCESS;
        break;
        
      case BSPLINE:
        range[0]  =  lgeom->data[0];
        range[1]  =  lgeom->data[lgeom->header[3]-1];
        *periodic = (lgeom->header[0]&4)/4;
        stat      =  EGADS_SUCCESS;
        break;
        
      case OFFSET:
        stat = EG_getRange(lgeom->ref, range, periodic);
        break;
    }
  } else if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) {
    lgeom = (liteGeometry *) geom->blind;
    mtype = geom->mtype;
    if (geom->oclass == EDGE) {
      if (geom->mtype == DEGENERATE) {
        ledge    = (liteEdge *) geom->blind;
        range[0] = ledge->trange[0];
        range[1] = ledge->trange[1];
        return EGADS_SUCCESS;
      }
      ledge = (liteEdge *) geom->blind;
      ref   = ledge->curve;
      if (ref == NULL)        return EGADS_NULLOBJ;
      if (ref->blind == NULL) return EGADS_NODATA;
      lgeom = (liteGeometry *) ref->blind;
      mtype = ref->mtype;
    }
    switch (mtype) {
      case LINE:
        range[0]  = -2.e100;
        range[1]  =  2.e100;
        *periodic =  0;
        stat      = EGADS_SUCCESS;
        break;
        
      case CIRCLE:
        range[0]  = 0.0;
        range[1]  = 2.0*PI;
        *periodic = 1;
        stat      = EGADS_SUCCESS;
        break;
        
      case ELLIPSE:
        range[0]  = 0.0;
        range[1]  = 2.0*PI;
        *periodic = 1;
        stat      = EGADS_SUCCESS;
        break;
        
      case PARABOLA:
        range[0]  = -2.e100;
        range[1]  =  2.e100;
        *periodic =  0;
        stat      = EGADS_SUCCESS;
        break;
        
      case HYPERBOLA:
        range[0]  = -2.e100;
        range[1]  =  2.e100;
        *periodic =  0;
        stat      = EGADS_SUCCESS;
        break;
        
      case TRIMMED:
        stat = EG_getRange(lgeom->ref, range, periodic);
        if (stat == EGADS_SUCCESS) {
          range[0] = lgeom->data[0];
          range[1] = lgeom->data[1];
        }
        break;
        
      case BEZIER:
        range[0]  = 0.0;
        range[1]  = 1.0;
        *periodic = (lgeom->header[0]&4)/4;
        stat      = EGADS_SUCCESS;
        break;
        
      case BSPLINE:
        range[0]  =  lgeom->data[0];
        range[1]  =  lgeom->data[lgeom->header[3]-1];
        *periodic = (lgeom->header[0]&4)/4;
        stat      =  EGADS_SUCCESS;
        break;
        
      case OFFSET:
        stat = EG_getRange(lgeom->ref, range, periodic);
        break;
    }
    if (geom->oclass == EDGE) {
      ledge    = (liteEdge *) geom->blind;
      range[0] = ledge->trange[0];
      range[1] = ledge->trange[1];
    }
  } else {
    lgeom = (liteGeometry *) geom->blind;
    mtype = geom->mtype;
    if (geom->oclass == FACE) {
      lface = (liteFace *) geom->blind;
      ref   = lface->surface;
      if (ref == NULL)        return EGADS_NULLOBJ;
      if (ref->blind == NULL) return EGADS_NODATA;
      lgeom = (liteGeometry *) ref->blind;
      mtype = ref->mtype;
    }
    switch (mtype) {
      case PLANE:
        range[0]  = -2.e100;
        range[1]  =  2.e100;
        range[2]  = -2.e100;
        range[3]  =  2.e100;
        *periodic =  0;
        stat      = EGADS_SUCCESS;
        break;
        
      case SPHERICAL:
        range[0]  =  0.0;
        range[1]  =  2.0*PI;
        range[2]  = -0.5*PI;
        range[3]  =  0.5*PI;
        *periodic = 1;
        stat      = EGADS_SUCCESS;
        break;
        
      case CONICAL:
        range[0]  =  0.0;
        range[1]  =  2.0*PI;
        range[2]  = -2.e100;
        range[3]  =  2.e100;
        *periodic = 1;
        stat      = EGADS_SUCCESS;
        break;
        
      case CYLINDRICAL:
        range[0]  =  0.0;
        range[1]  =  2.0*PI;
        range[2]  = -2.e100;
        range[3]  =  2.e100;
        *periodic = 1;
        stat      = EGADS_SUCCESS;
        break;
        
      case TOROIDAL:
        range[0]  = 0.0;
        range[1]  = 2.0*PI;
        range[2]  = 0.0;
        range[3]  = 2.0*PI;
        *periodic = 3;
        stat      = EGADS_SUCCESS;
        break;
        
      case REVOLUTION:
        stat = EG_getRange(lgeom->ref, &range[2], &mtype);
        if (stat == EGADS_SUCCESS) {
          range[0]  = 0.0;
          range[1]  = 2.0*PI;
          *periodic = 1 + 2*mtype;
        }
        break;
        
      case EXTRUSION:
        stat = EG_getRange(lgeom->ref, range, periodic);
        if (stat == EGADS_SUCCESS) {
          range[2] = -2.e100;
          range[3] =  2.e100;
        }
        break;
        
      case TRIMMED:
        stat = EG_getRange(lgeom->ref, range, periodic);
        if (stat == EGADS_SUCCESS) {
          range[0] = lgeom->data[0];
          range[1] = lgeom->data[1];
          range[2] = lgeom->data[2];
          range[3] = lgeom->data[3];
        }
        break;
        
      case BEZIER:
        range[0]  = 0.0;
        range[1]  = 1.0;
        range[2]  = 0.0;
        range[3]  = 1.0;
        *periodic = (lgeom->header[0]&12)/4;
        stat      = EGADS_SUCCESS;

        break;
        
      case BSPLINE:
        range[0]  =  lgeom->data[0];
        range[1]  =  lgeom->data[lgeom->header[3]-1];
        range[2]  =  lgeom->data[lgeom->header[3]];
        range[3]  =  lgeom->data[lgeom->header[3]+lgeom->header[6]-1];
        *periodic = (lgeom->header[0]&12)/4;
        stat      =  EGADS_SUCCESS;
        break;
        
      case OFFSET:
        stat = EG_getRange(lgeom->ref, range, periodic);
        break;
    }
    if (geom->oclass == FACE) {
      lface = (liteFace *) geom->blind;
      range[0] = lface->urange[0];
      range[1] = lface->urange[1];
      range[2] = lface->vrange[0];
      range[3] = lface->vrange[1];
    }
  }

  return stat;
}


int
EG_evaluate(const egObject *geom, /*@null@*/ const double *param, double *result)
{
  const egObject *ref;
  liteNode       *lnode;
  liteEdge       *ledge;
  liteFace       *lface;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != NODE)  && (geom->oclass != PCURVE)  &&
      (geom->oclass != CURVE) && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)  && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;
  
  /* special Node section */
  if (geom->oclass == NODE) {
    lnode     = (liteNode *) geom->blind;
    result[0] = lnode->xyz[0];
    result[1] = lnode->xyz[1];
    result[2] = lnode->xyz[2];
    return EGADS_SUCCESS;
  }
  if (param == NULL)               return EGADS_NODATA;
  
  /* geometry */
  if ((geom->oclass == PCURVE) || (geom->oclass == CURVE) ||
      (geom->oclass == SURFACE)) return EG_evaluateGeom(geom, param, result);

  /* topology */
  if (geom->oclass == EDGE) {
    ledge = (liteEdge *) geom->blind;
    ref   = ledge->curve;
  } else {
    lface = (liteFace *) geom->blind;
    ref   = lface->surface;
  }
  
  if (ref == NULL)        return EGADS_NULLOBJ;
  if (ref->blind == NULL) return EGADS_NODATA;
  return EG_evaluateGeom(ref, param, result);
}


int
EG_invEvaLimits(const egObject *geom, /*@null@*/ const double *limits,
                const double *xyz, double *param, double *result)
{
  int            stat, per;
  double         range[4], srange[4], pt[3], uvs[2], period;
  liteEdge       *pedge;
  liteFace       *pface;
  liteGeometry   *lgeom;
  const egObject *ref;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  
  if ((geom->oclass == PCURVE) || (geom->oclass == CURVE) ||
      (geom->oclass == SURFACE))
    return EG_invEvaGeomLimits(geom, limits, xyz, param, 0.0, result);
  
  stat = EG_getRange(geom, range, &per);
  if (stat != EGADS_SUCCESS) return stat;
  
  /* do we re-limit? */
  if (limits != NULL) {
    range[0] = limits[0];
    range[1] = limits[1];
    if (geom->oclass == FACE) {
      range[2] = limits[2];
      range[3] = limits[3];
    }
  }
  
  if (geom->oclass == EDGE) {
    pedge = (liteEdge *) geom->blind;
    ref   = pedge->curve;
    if (ref == NULL)        return EGADS_NULLOBJ;
    if (ref->blind == NULL) return EGADS_NODATA;
    return EG_invEvaGeomLimits(ref, range, xyz, param, 0.0, result);
  }

  /* do the Face */
  pface = (liteFace *) geom->blind;
  ref   = pface->surface;
  if (ref == NULL)        return EGADS_NULLOBJ;
  if (ref->blind == NULL) return EGADS_NODATA;
  if (ref->mtype == TRIMMED) {
    lgeom = (liteGeometry *) ref->blind;
    ref   = lgeom->ref;
    if (ref->blind == NULL) return EGADS_NODATA;
    if (ref->mtype == TRIMMED)
      printf(" EGADSlite Internal: TRIMMED TRIMMED Surface!\n");
  }
  stat = EG_invEvaGeomLimits(ref, range, xyz, param, pface->tol, result);
  if (stat != EGADS_SUCCESS) return stat;
  
  stat = EG_getRange(ref, srange, &per);
  if (stat != EGADS_SUCCESS) return stat;

  stat = EG_inFaceX(geom, param, pt, uvs);
  if (stat < EGADS_SUCCESS) return stat;
  if (stat == EGADS_OUTSIDE) {
/*  printf(" Info: Point labelled outside!\n");  */
    param[0]  = uvs[0];
    param[1]  = uvs[1];
    result[0] = pt[0];
    result[1] = pt[1];
    result[2] = pt[2];
    if ((per&1) != 0) {
      period = srange[1] - srange[0];
      if ((param[0]+PARAMACC < srange[0]) || (param[0]-PARAMACC > srange[1]))
        if (param[0]+PARAMACC < srange[0]) {
          if (param[0]+period-PARAMACC < srange[1]) param[0] += period;
        } else {
          if (param[0]-period+PARAMACC > srange[0]) param[0] -= period;
        }
    }
    if ((per&2) != 0) {
      period = srange[3] - srange[2];
      if ((param[1]+PARAMACC < srange[2]) || (param[1]-PARAMACC > srange[3]))
        if (param[1]+PARAMACC < srange[2]) {
          if (param[1]+period-PARAMACC < srange[3]) param[1] += period;
        } else {
          if (param[1]-period+PARAMACC > srange[2]) param[1] -= period;
        }
    }
  }

  return EGADS_SUCCESS;
}


int
EG_invEvaluate(const egObject *geom, const double *xyz, double *param,
               double *result)
{
  return EG_invEvaLimits(geom, NULL, xyz, param, result);
}


int
EG_invEvaluateGuess(const egObject *geom, double *xyz,
                    double *param, double *result)
{
  int            stat, per;
  double         range[2];
  const egObject *ref;
  liteEdge       *ledge;
  liteFace       *lface;
  
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  
  if ((geom->oclass == PCURVE) || (geom->oclass == CURVE) ||
      (geom->oclass == SURFACE))
    return EG_invEvaluateGeomGuess(geom, NULL, xyz, param, result);
  
  if (geom->oclass == EDGE) {
    
    stat = EG_getRange(geom, range, &per);
    if (stat != EGADS_SUCCESS) return stat;
    ledge = (liteEdge *) geom->blind;
    ref   = ledge->curve;
    if (ref == NULL)        return EGADS_NULLOBJ;
    if (ref->blind == NULL) return EGADS_NODATA;
    stat = EG_invEvaluateGeomGuess(ref, range, xyz, param, result);
    
  } else {
    
    lface = (liteFace *) geom->blind;
    ref   = lface->surface;
    if (ref == NULL)        return EGADS_NULLOBJ;
    if (ref->blind == NULL) return EGADS_NODATA;
    stat = EG_invEvaluateGeomGuess(ref, NULL, xyz, param, result);
    
  }

  return stat;
}


int
EG_arcLength(const egObject *geom, double t1, double t2, double *alen)
{
  int            i, stat;
  double         t, d, ur, mid, result[9];
  liteEdge       *ledge;
  const egObject *ref;
/*
  static int     ngauss   = 5;
  static double  wg[2*5]  = { 0.5688888888888889,  0.0000000000000000,
                              0.4786286704993665, -0.5384693101056831,
                              0.4786286704993665,  0.5384693101056831,
                              0.2369268850561891, -0.9061798459386640,
                              0.2369268850561891,  0.9061798459386640 };
 */
  static int     ngauss   = 15;
  static double  wg[2*15] = { 0.2025782419255613,  0.0000000000000000,
                              0.1984314853271116, -0.2011940939974345,
                              0.1984314853271116,  0.2011940939974345,
                              0.1861610000155622, -0.3941513470775634,
                              0.1861610000155622,  0.3941513470775634,
                              0.1662692058169939, -0.5709721726085388,
                              0.1662692058169939,  0.5709721726085388,
                              0.1395706779261543, -0.7244177313601701,
                              0.1395706779261543,  0.7244177313601701,
                              0.1071592204671719, -0.8482065834104272,
                              0.1071592204671719,  0.8482065834104272,
                              0.0703660474881081, -0.9372733924007060,
                              0.0703660474881081,  0.9372733924007060,
                              0.0307532419961173, -0.9879925180204854,
                              0.0307532419961173,  0.9879925180204854 };
  
  *alen = 0.0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) && (geom->oclass != CURVE) &&
      (geom->oclass != EDGE))      return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;

  if (geom->oclass == PCURVE) {
    
    ur  =      t2 - t1;
    mid = 0.5*(t2 + t1);
    for (i = 0; i < ngauss; i++) {
      t    = t1 + 0.5*wg[2*i+1]*ur + mid;
      stat = EG_evaluate(geom, &t, result);
      if (stat != EGADS_SUCCESS) return stat;
      d = sqrt(result[2]*result[2] + result[3]*result[3]);
      *alen += d*wg[2*i];
    }
    *alen *= 0.5*ur;
    
  } else {
    
    ref = geom;
    if (geom->oclass == EDGE) {
      ledge = (liteEdge *) geom->blind;
      ref   = ledge->curve;
      if (ref == NULL)        return EGADS_NULLOBJ;
      if (ref->blind == NULL) return EGADS_NODATA;
    }
    ur  =      t2 - t1;
    mid = 0.5*(t2 + t1);
    for (i = 0; i < ngauss; i++) {
      t    = t1 + 0.5*wg[2*i+1]*ur + mid;
      stat = EG_evaluate(ref, &t, result);
      if (stat != EGADS_SUCCESS) return stat;
      d = sqrt(result[3]*result[3] + result[4]*result[4] + result[5]*result[5]);
      *alen += d*wg[2*i];
    }
    *alen *= 0.5*ur;
  }

  return EGADS_SUCCESS;
}


int
EG_curvature(const egObject *geom, const double *param, double *result)
{
  int    i, stat;
  double data[18], dir[3], d, s, *d1, *d2;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;
  
  stat = EG_evaluate(geom, param, data);
  if (stat != EGADS_SUCCESS) return stat;
  
  if (geom->oclass == PCURVE) {
    
    for (i = 0; i < 3; i++) result[i] = 0.0;
    
    s         = sqrt(data[2]*data[2] + data[3]*data[3]);
    if (s == 0.0) return EGADS_DEGEN;
    d         = ABS(data[2]*data[5] - data[3]*data[4]);
    result[0] = d/(s*s*s);
    result[1] = data[2]/s;
    result[2] = data[3]/s;
  
  } else if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) {
    
    for (i = 0; i < 4; i++) result[i] = 0.0;

    s         = sqrt(data[3]*data[3] + data[4]*data[4] + data[5]*data[5]);
    if (s == 0.0) return EGADS_DEGEN;
    d1        = &data[3];
    d2        = &data[6];
    CROSS(dir, d1, d2);
    d         = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    result[0] = d/(s*s*s);
    result[1] = data[3]/s;
    result[2] = data[4]/s;
    result[3] = data[5]/s;
    
  } else {
    
    double norm[3], *der1[2], *der2[3];
    double a, b, c, d11, d12, d21, d22, g11, g12, g21, g22, ud, vd, len;
    
    for (i = 0; i < 8; i++) result[i] = 0.0;
    der1[0]  = &data[ 3];
    der1[1]  = &data[ 6];
    der2[0]  = &data[ 9];
    der2[1]  = &data[15];
    der2[2]  = &data[12];
    norm[0]  = der1[0][1]*der1[1][2] - der1[0][2]*der1[1][1];
    norm[1]  = der1[0][2]*der1[1][0] - der1[0][0]*der1[1][2];
    norm[2]  = der1[0][0]*der1[1][1] - der1[0][1]*der1[1][0];
    len = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    if (len == 0.0) return EGADS_DEGEN;
    norm[0] /= len;
    norm[1] /= len;
    norm[2] /= len;
    if ((geom->oclass == FACE) && (geom->mtype == SREVERSE)) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
      norm[2] = -norm[2];
    }
    
    g11 = g12 = g21 = g22 = 0.0;
    d11 = d12 = d21 = d22 = 0.0;
    for (i = 0; i < 3; i++) {
      g11 += der1[0][i]*der1[0][i];
      g12 += der1[0][i]*der1[1][i];
      g21 += der1[1][i]*der1[0][i];
      g22 += der1[1][i]*der1[1][i];
      
      d11 += norm[i]*der2[0][i];
      d12 += norm[i]*der2[2][i];
      d21 += norm[i]*der2[2][i];
      d22 += norm[i]*der2[1][i];
    }
    a   =   g11*g22 - g21*g12;
    b   = -(g11*d22 + d11*g22 - 2.0*g12*d21);
    c   =   d11*d22 - d21*d12;
    len = b*b - 4.0*a*c;
    if (len < 0.0) len = 0.0;
    
    result[0] = (-b + sqrt(len))/(2.0*a);
    result[4] = (-b - sqrt(len))/(2.0*a);
    
    if (ABS(result[0]-result[4]) > 1.e-12) {
      
      /* find principal direction 1 */
      ud =  (d12 - result[0]*g12);
      vd = -(d11 - result[0]*g11);
      if ((ABS(ud) < 1.e-12) && (ABS(vd) < 1.e-12)) {
        ud =  (d22 - result[0]*g22);
        vd = -(d21 - result[0]*g21);
      }
      for (i = 0; i < 3; i++) result[i+1] = der1[0][i]*ud + der1[1][i]*vd;
      len = result[1]*result[1] + result[2]*result[2] + result[3]*result[3];
      if (len != 0.0) {
        len = sqrt(1.0/len);
        for (i = 1; i < 4; i++) result[i] *= len;
      }
 
    } else {

      /*  Align principal direction 1 with isocurves */
      len = 0;
      for (i = 0; i < 3; i++) {
        result[i+1] = der1[0][i];
        len        += der1[0][i]*der1[0][i];
      }
      if (len == 0.0) {
        for (i = 0; i < 3; i++) {
          result[i+1] = der1[1][i];
          len        += der1[1][i]*der1[1][i];
        }
      }
      if (len == 0.0) return EGADS_DEGEN;
      len = 1.0/sqrt(len);
      for (i = 1; i < 4; i++) result[i] *= len;
      
    }
    
    /* find principal direction 2 -- make orthogonal */
    d1 = &result[1];
    d2 = &result[5];
    CROSS(d2, d1, norm);
  
  }
  
  return EGADS_SUCCESS;
}


int
EG_relPosTs(egObject *geom, int n, /*@null@*/ const double *rel,
            double *ts, double *xyzs)
{
  int      i, j, stat;
  double   alen, len, t, tbeg, tend, frac, result[9];
  egObject *ref;
  liteEdge *ledge;
  
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != CURVE)  &&
      (geom->oclass != EDGE))      return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  ref = geom;
  if (geom->oclass == EDGE) {
    ledge = (liteEdge *) geom->blind;
    ref   = ledge->curve;
    if (ref == NULL)               return EGADS_NULLOBJ;
    if (ref->blind == NULL)        return EGADS_NODATA;
  }

  /* is t arc-length based? */
  if ((ref->mtype != TRIMMED) && (ref->mtype != OFFSET) &&
      (ref->mtype != BEZIER)  && (ref->mtype != BSPLINE)) {
    for (i = 1; i < n-1; i++) {
      if (rel == NULL) {
        t  = i;
        t  = ts[0] + t*(ts[n-1]-ts[0])/(n-1);
      } else {
        t  = ts[0] + rel[i-1]*(ts[n-1]-ts[0]);
      }
      stat = EG_evaluate(geom, &t, result);
      if (stat != EGADS_SUCCESS) return stat;
      ts[i]       = t;
      xyzs[3*i  ] = result[0];
      xyzs[3*i+1] = result[1];
      xyzs[3*i+2] = result[2];
    }
    return EGADS_SUCCESS;
  }
  
  /* not sure if t is arc-length based */
  stat = EG_arcLength(geom, ts[0], ts[n-1], &alen);
  if (stat != EGADS_SUCCESS) return stat;
  if (alen == 0.0) return EGADS_DEGEN;

  for (i = 1; i < n-1; i++) {
    if (rel == NULL) {
      frac  = i;
      frac /= n-1;
    } else {
      frac  = rel[i-1];
    }
    /* this needs to be improved -- only interval halving! */
    tbeg = ts[0];
    tend = ts[n-1];
    t    = 0.5*(tbeg + tend);
    for (j = 0; j < 20; j++) {
      stat = EG_arcLength(geom, ts[0], t, &len);
      if (stat != EGADS_SUCCESS) return stat;
/*    printf("  %d/%d:  %lf  %lf (%lf)\n", i, j, t, len/alen, rel[i-1]);  */
      if (len/alen == frac) break;
      if (len/alen >  frac) {
        tend = t;
      } else {
        tbeg = t;
      }
      t = 0.5*(tbeg + tend);
    }
    stat = EG_evaluate(geom, &t, result);
    if (stat != EGADS_SUCCESS) return stat;
    ts[i]       = t;
    xyzs[3*i  ] = result[0];
    xyzs[3*i+1] = result[1];
    xyzs[3*i+2] = result[2];
  }
  
  return EGADS_SUCCESS;
}
