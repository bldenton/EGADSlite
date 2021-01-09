/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Geometry Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsClasses.h"
#define TEMPLATE template<class TT>
#define DOUBLE TT
#define CROSS(a,b,c)       a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                           a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                           a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)          (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


#define PARAMACC 1.0e-4         // parameter accuracy
#define KNACC    1.0e-12	// knot accuracy

/* OCC can change p-curves through transformations...
 * so sensitivities for p-curves is optional for now.
 */
//#define REQUIRE_PCURVE_SENSITIVITIES


  extern "C" int  EG_fixedKnots( const ego object );
  extern "C" int  EG_destroyGeometry( egObject *geom );
  extern "C" int  EG_copyGeometry( /*@null@*/ egObject *cxt, const egObject *geo,
                                   /*@null@*/ double *xform, egObject **copy );
  extern "C" int  EG_flipGeometry( const egObject *geom, egObject **copy );
  extern "C" int  EG_spline1d( egObject *context, int endc, int imax,
                               const double *xyz, double tol, egObject **ecrv );
  extern "C" int  EG_spline2d( egObject *context, int endc,
                               /*@null@*/ const double **dr, int imax, int jmax,
                               const double *xyz, double tol, egObject **esrf );
  extern "C" int  EG_spline2dFit( egObject *context, const double *crosT,
                                  int imax, const double *uknot,
                                  const double *souT, const double *norT,
                                  int jmax, const double *vknot,
                                  const double *wesT, const double *easT,
                                  const double *xyz, double tol,
                                  egObject **esurf );
  extern "C" int EG_spline1dFit_dot( int endx, int imaxx,
                                     const double *xyz, const double *xyz_dot,
                                     const double *kn, const double *kn_dot,
                                     double tol, int *ivec,
                                     double **rdata, double **rdata_dot );
  template<class T>
  int EG_spline2dAppr( int endc, int imax, int jmax, const T *xyz,
                       /*@null@*/ const T   *uknot, /*@null@*/ const T *vknot,
                       /*@null@*/ const int *vdata,
                       /*@null@*/ const T   *wesT,  /*@null@*/ const T *easT,
                       /*@null@*/ const T   *south, /*@null@*/       T *snor,
                       /*@null@*/ const T   *north, /*@null@*/       T *nnor,
                       double tol, int *header, T **rdata );
  extern     int  EG_copyAttrTopo( egadsBody *pbody, /*@null@*/ double *xform,
                                   gp_Trsf form, const egObject *src,
                                   egObject *dst, egObject *topObj );
  extern     void EG_checkStatus( const Handle_BRepCheck_Result tResult );
  TEMPLATE   int  EG_evaluateGeom( const egObject *geom, const DOUBLE *param,
                                   DOUBLE *result );
  extern     int  EG_invEvaGeomLimits( const egObject *geom,
                                       /*@null@*/ const double *limits,
                                       const double *xyz, double *param,
                                       double toler, double *result );
  extern     int  EG_invEvaluateGeomGuess( const egObject *geom,
                                           /*@null@*/ const double *limits,
                                           double *xyz, double *param,
                                           double *result );
  extern     int  EG_inFaceX( const egObject *face, const double *uva,
                              /*@null@*/ double *pt, /*@null@*/ double *uvx );

  extern "C" int  EG_getGeometry( const egObject *geom, int *oclass, int *type,
                                  egObject **rGeom, /*@null@*/ int **ivec,
                                  /*@null@*/ double **rvec );
  extern "C" int  EG_setGeometry_dot( egObject *geom, int oclass, int mtype,
                                      /*@null@*/ const int *ivec,
                                      const double *rvec,
                                      const double *rvec_dot );
  extern "C" int  EG_getGeometry_dot( const egObject *geom, double **rvec,
                                      double **rvec_dot );
             int  EG_getGeometry_dot( const egObject *obj,
                                      SurrealS<1> **data_dot );
  extern "C" int  EG_makeGeometry( egObject *context, int oclass, int mtype,
                                   /*@null@*/ egObject *refGeo, const int *ivec,
                                   const double *rvec, egObject **geom );
  extern "C"  int EG_copyGeometry_dot(const egObject *obj,
                                      /*@null@*/ const double *xform,
                                      /*@null@*/ const double *xform_dot,
                                      egObject *copy);
  extern "C" int  EG_hasGeometry_dot(const egObject *obj);
  extern "C" int  EG_getRange( const egObject *geom, double *range, int *pflg );
  extern "C" int  EG_getRange_dot( const egObject *geom, double *range,
                                   double *range_dot, int *pflg );
  extern "C" int  EG_curvature( const egObject *geom, const double *param,
                                double *result );
  extern "C" int  EG_evaluate( const egObject *geom, const double *param,
                               double *result );
  extern "C" int  EG_evaluate_dot( const egObject *geom,
                                   const double *param, const double *param_dot,
                                   double *result, double *result_dot );
  extern "C" int  EG_invEvaluate( const egObject *geom, double *xyz,
                                  double *param, double *result );
  extern "C" int  EG_invEvaluateGuess( const egObject *geom, double *xyz,
                                       double *param, double *result );
  extern "C" int  EG_arcLength( const egObject *geom, double t1, double t2,
                                double *alen );
  extern "C" int  EG_approximate( egObject *context, int maxdeg, double tol,
                                  const int *sizes, const double *xyzs,
                                  egObject **bspline );
  extern "C" int  EG_approximate_dot( egObject *bspline, int maxdeg, double tol,
                                      const int *sizes,
                                      const double *data, const double *data_dot );
  extern "C" int  EG_otherCurve( const egObject *surface, const egObject *curve,
                                 double tol, egObject **newcurve );
  extern "C" int  EG_isoCline( const egObject *surface, int UV, double value,
                               egObject **newcurve );
  extern "C" int  EG_convertToBSplineRange( egObject *geom, const double *range,
                                            egObject **bspline );
  extern "C" int  EG_convertToBSpline( egObject *geom, egObject **bspline );
  extern "C" int  EG_flattenBSpline( egObject *geom, egObject **bspline );
  extern "C" int  EG_addKnots( const egObject *object, int nU,
                               /*@null@*/ double *Us, int nV,
                               /*@null@*/ double *Vs, egObject **result );
  extern "C" int  EG_mapKnots( egObject *src, egObject *dst, egObject **rslt );
  extern "C" int  EG_mapSequen( egObject *src, egObject *dst, egObject **rslt );
  extern "C" void EG_mapTessTs( egTess1D src, egTess1D dst );
  extern "C" int  EG_relPosTs( egObject *geom, int n, const double *rel,
                               double *ts, double *xyzs );


int
EG_addStrAttr(egObject *obj, const char *name, const char *str)
{
  int     i, length, find = -1;
  egAttr  *attr;
  egAttrs *attrs;

  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;
  if (EG_sameThread(obj))        return EGADS_CNTXTHRD;

  if ((name == NULL) || (str == NULL)) {
    printf(" EGADS Internal: NULL Name/Value (EG_addStrAttr)!\n");
    return EGADS_NONAME;
  }
  length = strlen(name);
  for (i = 0; i < length; i++)
    if (name[i] <= ' ') {
      length = 0;
      break;
    }
  if (length == 0) {
    printf(" EGADS Internal: BAD Name (EG_addStrAttr)!\n");
    return EGADS_INDEXERR;
  }
  attrs = (egAttrs *) obj->attrs;

  if (attrs != NULL)
    for (i = 0; i < attrs->nattrs; i++)
      if (strcmp(attrs->attrs[i].name,name) == 0) {
        find = i;
        break;
      }

  if ((find != -1) && (attrs != NULL)) {

    /* an existing attribute -- reset the values */

    if (attrs->attrs[find].type == ATTRINT) {
      if (attrs->attrs[find].length != 1)
        EG_free(attrs->attrs[find].vals.integers);
    } else if (attrs->attrs[find].type == ATTRREAL) {
      if (attrs->attrs[find].length != 1)
        EG_free(attrs->attrs[find].vals.reals);
    } else if (attrs->attrs[find].type == ATTRCSYS) {
      EG_free(attrs->attrs[find].vals.reals);
    } else if (attrs->attrs[find].type == ATTRSTRING) {
      EG_free(attrs->attrs[find].vals.string);
    }

  } else {

    if (attrs == NULL) {
      attrs = (egAttrs *) EG_alloc(sizeof(egAttrs));
      if (attrs == NULL) {
        printf(" EGADS Internal: Attrs MALLOC for %s (EG_addStrAttr)!\n",
               name);
        return EGADS_MALLOC;
      }
      attrs->nattrs = 0;
      attrs->attrs  = NULL;
      attrs->nseqs  = 0;
      attrs->seqs   = NULL;
      obj->attrs    = attrs;
    }
    if (attrs->attrs == NULL) {
      attr = (egAttr *) EG_alloc((attrs->nattrs+1)*sizeof(egAttr));
    } else {
      attr = (egAttr *) EG_reall(attrs->attrs,
                                 (attrs->nattrs+1)*sizeof(egAttr));
    }
    if (attr == NULL) {
      printf(" EGADS Internal: Attr MALLOC for %s (EG_addStrAttr)!\n",
             name);
      return EGADS_MALLOC;
    }
    attrs->attrs = attr;
    find = attrs->nattrs;
    attrs->attrs[find].vals.string = NULL;
    attrs->attrs[find].name        = EG_strdup(name);
    if (attrs->attrs[find].name == NULL) return EGADS_MALLOC;
    attrs->nattrs += 1;
  }

  attrs->attrs[find].type        = ATTRSTRING;
  attrs->attrs[find].length      = 0;
  attrs->attrs[find].vals.string = EG_strdup(str);
  if (attrs->attrs[find].vals.string != NULL)
    attrs->attrs[find].length = strlen(attrs->attrs[find].vals.string);

  return EGADS_SUCCESS;
}


int
EG_destroyGeometry(egObject *geom)
{
  egObject *obj = NULL;

  if (geom->oclass == PCURVE) {

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    if (ppcurv != NULL) {
      if (ppcurv->header   != NULL) EG_free(ppcurv->header);
      if (ppcurv->data     != NULL) EG_free(ppcurv->data);
      if (ppcurv->data_dot != NULL) EG_free(ppcurv->data_dot);
      obj = ppcurv->ref;
    }
    if (obj    != NULL)
      if (ppcurv->topFlg == 0) {
        EG_dereferenceObject(obj, geom);
      } else {
        EG_dereferenceTopObj(obj, geom);
      }
    if (ppcurv != NULL) delete ppcurv;

  } else if (geom->oclass == CURVE) {

    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    if (pcurve != NULL) {
      if (pcurve->header   != NULL) EG_free(pcurve->header);
      if (pcurve->data     != NULL) EG_free(pcurve->data);
      if (pcurve->data_dot != NULL) EG_free(pcurve->data_dot);
      obj = pcurve->ref;
    }
    if (obj    != NULL)
      if (pcurve->topFlg == 0) {
        EG_dereferenceObject(obj, geom);
      } else {
        EG_dereferenceTopObj(obj, geom);
      }
    if (pcurve != NULL) delete pcurve;

  } else {

    egadsSurface *psurf = (egadsSurface *) geom->blind;
    if (psurf != NULL) {
      if (psurf->header   != NULL) EG_free(psurf->header);
      if (psurf->data     != NULL) EG_free(psurf->data);
      if (psurf->data_dot != NULL) EG_free(psurf->data_dot);
      obj = psurf->ref;
    }
    if (obj   != NULL)
      if (psurf->topFlg == 0) {
        EG_dereferenceObject(obj, geom);
      } else {
        EG_dereferenceTopObj(obj, geom);
      }
    if (psurf != NULL) delete psurf;

  }
  return EGADS_SUCCESS;
}


static int
EG_getPCurveType(Handle(Geom2d_Curve) &hCurve)
{

  Handle(Geom2d_Line) hLine = Handle(Geom2d_Line)::DownCast(hCurve);
  if (!hLine.IsNull()) return LINE;

  Handle(Geom2d_Circle) hCirc = Handle(Geom2d_Circle)::DownCast(hCurve);
  if (!hCirc.IsNull()) return CIRCLE;

  Handle(Geom2d_Ellipse) hEllip = Handle(Geom2d_Ellipse)::DownCast(hCurve);
  if (!hEllip.IsNull()) return ELLIPSE;

  Handle(Geom2d_Parabola) hParab = Handle(Geom2d_Parabola)::DownCast(hCurve);
  if (!hParab.IsNull()) return PARABOLA;

  Handle(Geom2d_Hyperbola) hHypr = Handle(Geom2d_Hyperbola)::DownCast(hCurve);
  if (!hHypr.IsNull()) return HYPERBOLA;

  Handle(Geom2d_BezierCurve) hBezier = Handle(Geom2d_BezierCurve)::DownCast(hCurve);
  if (!hBezier.IsNull()) return BEZIER;

  Handle(Geom2d_BSplineCurve) hBSpline = Handle(Geom2d_BSplineCurve)::DownCast(hCurve);
  if (!hBSpline.IsNull()) return BSPLINE;

  Handle(Geom2d_TrimmedCurve) hTrim = Handle(Geom2d_TrimmedCurve)::DownCast(hCurve);
  if (!hTrim.IsNull()) return TRIMMED;

  Handle(Geom2d_OffsetCurve) hOffst = Handle(Geom2d_OffsetCurve)::DownCast(hCurve);
  if (!hOffst.IsNull()) return OFFSET;

  return 0;
}


void
EG_getGeometryLen(const egObject *geom, int *nivec, int *nrvec)
{
  int nint, nreal;

  *nivec = *nrvec = nint = nreal = 0;
  if  (geom == NULL)               return;
  if  (geom->magicnumber != MAGIC) return;
  if ((geom->oclass < PCURVE) || (geom->oclass > SURFACE))
                                   return;
  if (geom->blind == NULL)         return;

  if (geom->oclass == PCURVE) {

    egadsPCurve *lgeom = (egadsPCurve *) geom->blind;
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

    egadsCurve *lgeom = (egadsCurve *) geom->blind;
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
    egadsSurface *lgeom = (egadsSurface *) geom->blind;
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

  *nivec = nint;
  *nrvec = nreal;
}


int
EG_getGeometry(const egObject     *geom, int *oclass, int *type,
                     egObject **refGeom, /*@null@*/ int    **ivec,
                                         /*@null@*/ double **rvec)
{
  int    *ints = NULL, i, j, len, outLevel;
  double *data = NULL;

  if (ivec != NULL) *ivec = NULL;
  if (rvec != NULL) *rvec = NULL;
  *refGeom = NULL;
  *oclass  = *type = 0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass < PCURVE) || (geom->oclass > SURFACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(geom);
  *oclass = geom->oclass;
  *type   = geom->mtype;

  if (geom->oclass == PCURVE) {

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    *refGeom = ppcurv->ref;
    if ((rvec == NULL) || (ivec == NULL)) return EGADS_SUCCESS;

    switch (geom->mtype) {

      case LINE:
        data = (double *) EG_alloc(4*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PLine (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Line) hLine = Handle(Geom2d_Line)::DownCast(hCurve);
          gp_Dir2d direct = hLine->Direction();
          gp_Pnt2d locat  = hLine->Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = direct.X();
          data[3] = direct.Y();
          *rvec   = data;
        }
        break;

      case CIRCLE:
        data = (double *) EG_alloc(7*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PCircle (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Circle) hCirc = Handle(Geom2d_Circle)::DownCast(hCurve);
          gp_Circ2d circ  = hCirc->Circ2d();
          gp_Ax2d   xaxis = circ.XAxis();
          gp_Ax2d   yaxis = circ.YAxis();
          gp_Pnt2d  locat = circ.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = xaxis.Direction().X();
          data[3] = xaxis.Direction().Y();
          data[4] = yaxis.Direction().X();
          data[5] = yaxis.Direction().Y();
          data[6] = circ.Radius();
          *rvec   = data;
        }
        break;

      case ELLIPSE:
        data = (double *) EG_alloc(8*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PEllipse (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Ellipse) hEllip = Handle(Geom2d_Ellipse)::DownCast(hCurve);
          gp_Elips2d elips = hEllip->Elips2d();
          gp_Ax2d    xaxis = elips.XAxis();
          gp_Ax2d    yaxis = elips.YAxis();
          gp_Pnt2d   locat = elips.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = xaxis.Direction().X();
          data[3] = xaxis.Direction().Y();
          data[4] = yaxis.Direction().X();
          data[5] = yaxis.Direction().Y();
          data[6] = elips.MajorRadius();
          data[7] = elips.MinorRadius();
          *rvec   = data;
        }
        break;

      case PARABOLA:
        data = (double *) EG_alloc(7*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PParabola (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Parabola) hParab=Handle(Geom2d_Parabola)::DownCast(hCurve);
          gp_Parab2d parab = hParab->Parab2d();
          gp_Ax22d   axes  = parab.Axis();
          gp_Ax2d    xaxis = axes.XAxis();
          gp_Ax2d    yaxis = axes.YAxis();
          gp_Pnt2d   locat = parab.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = xaxis.Direction().X();
          data[3] = xaxis.Direction().Y();
          data[4] = yaxis.Direction().X();
          data[5] = yaxis.Direction().Y();
          data[6] = parab.Focal();
          *rvec   = data;
        }
        break;

      case HYPERBOLA:
        data = (double *) EG_alloc(8*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PHyperbola (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Hyperbola) hHypr =
            Handle(Geom2d_Hyperbola)::DownCast(hCurve);
          gp_Hypr2d hypr  = hHypr->Hypr2d();
          gp_Ax2d   xaxis = hypr.XAxis();
          gp_Ax2d   yaxis = hypr.YAxis();
          gp_Pnt2d  locat = hypr.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = xaxis.Direction().X();
          data[3] = xaxis.Direction().Y();
          data[4] = yaxis.Direction().X();
          data[5] = yaxis.Direction().Y();
          data[6] = hypr.MajorRadius();
          data[7] = hypr.MinorRadius();
          *rvec   = data;
        }
        break;

      case TRIMMED:
        data = (double *) EG_alloc(2*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PTrimmed Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_TrimmedCurve) hTrim =
            Handle(Geom2d_TrimmedCurve)::DownCast(hCurve);
          data[0] = hTrim->FirstParameter();
          data[1] = hTrim->LastParameter();
          *rvec   = data;
        }
        break;

      case BEZIER:
        ints = (int *) EG_alloc(3*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PBezier Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_BezierCurve) hBezier =
            Handle(Geom2d_BezierCurve)::DownCast(hCurve);
          int rational = 0;
          if (hBezier->IsRational()) rational = 1;
          ints[0] = rational*2;
          if (hBezier->IsPeriodic()) ints[0] |= 4;
          ints[1] = hBezier->Degree();
          ints[2] = hBezier->NbPoles();
          len = ints[2]*2;
          if (rational == 1) len += ints[2];
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc on PBezier Data (EG_getGeometry)!\n");
            EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= ints[2]; i++, len+=2) {
            gp_Pnt2d P  = hBezier->Pole(i);
            data[len  ] = P.X();
            data[len+1] = P.Y();
          }
          if (rational == 1)
            for (i = 1; i <= ints[2]; i++,len++) data[len] = hBezier->Weight(i);
          *ivec = ints;
          *rvec = data;
        }
        break;

      case BSPLINE:
        ints = (int *) EG_alloc(4*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PBSpline Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_BSplineCurve) hBSpline =
            Handle(Geom2d_BSplineCurve)::DownCast(hCurve);
          int rational = 0;
          if (hBSpline->IsRational()) rational = 1;
          ints[0] = rational*2;
          if (hBSpline->IsPeriodic()) ints[0] |= 4;
          ints[1] = hBSpline->Degree();
          ints[2] = hBSpline->NbPoles();
          ints[3] = 0;
          for (i = 1; i <= hBSpline->NbKnots(); i++)
            ints[3] += hBSpline->Multiplicity(i);
          len = ints[3] + ints[2]*2;
          if (rational == 1) len += ints[2];
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc PBSpline Data (EG_getGeometry)!\n");
            EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= hBSpline->NbKnots(); i++) {
            int km = hBSpline->Multiplicity(i);
            for (j = 1; j <= km; j++, len++) data[len] = hBSpline->Knot(i);
          }
          for (i = 1; i <= ints[2]; i++, len+=2) {
            gp_Pnt2d P  = hBSpline->Pole(i);
            data[len  ] = P.X();
            data[len+1] = P.Y();
          }
          if (rational == 1)
            for (i = 1; i <= ints[2]; i++,len++)
              data[len] = hBSpline->Weight(i);
          *ivec = ints;
          *rvec = data;
        }
        break;

      case OFFSET:
        data = (double *) EG_alloc(1*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on POffset Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_OffsetCurve) hOffst =
            Handle(Geom2d_OffsetCurve)::DownCast(hCurve);
          data[0] = hOffst->Offset();
          *rvec   = data;
        }
        break;
    }

  } else if (geom->oclass == CURVE) {

    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    *refGeom = pcurve->ref;
    if ((rvec == NULL) || (ivec == NULL)) return EGADS_SUCCESS;

    switch (geom->mtype) {

      case LINE:
        data = (double *) EG_alloc(6*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Line (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Line) hLine = Handle(Geom_Line)::DownCast(hCurve);
          gp_Lin line   = hLine->Lin();
          gp_Dir direct = line.Direction();
          gp_Pnt locat  = line.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = direct.X();
          data[4] = direct.Y();
          data[5] = direct.Z();
          *rvec   = data;
        }
        break;

      case CIRCLE:
        data = (double *) EG_alloc(10*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Circle (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Circle) hCirc = Handle(Geom_Circle)::DownCast(hCurve);
          gp_Circ circ  = hCirc->Circ();
          gp_Ax1  xaxis = circ.XAxis();
          gp_Ax1  yaxis = circ.YAxis();
          gp_Pnt  locat = circ.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = xaxis.Direction().X();
          data[4] = xaxis.Direction().Y();
          data[5] = xaxis.Direction().Z();
          data[6] = yaxis.Direction().X();
          data[7] = yaxis.Direction().Y();
          data[8] = yaxis.Direction().Z();
          data[9] = circ.Radius();
          *rvec   = data;
        }
        break;

      case ELLIPSE:
        data = (double *) EG_alloc(11*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Ellipse (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Ellipse) hEllip = Handle(Geom_Ellipse)::DownCast(hCurve);
          gp_Elips elips = hEllip->Elips();
          gp_Ax1   xaxis = elips.XAxis();
          gp_Ax1   yaxis = elips.YAxis();
          gp_Pnt   locat = elips.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = xaxis.Direction().X();
          data[ 4] = xaxis.Direction().Y();
          data[ 5] = xaxis.Direction().Z();
          data[ 6] = yaxis.Direction().X();
          data[ 7] = yaxis.Direction().Y();
          data[ 8] = yaxis.Direction().Z();
          data[ 9] = elips.MajorRadius();
          data[10] = elips.MinorRadius();
          *rvec    = data;
        }
        break;

      case PARABOLA:
        data = (double *) EG_alloc(10*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Parabola (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Parabola) hParab=Handle(Geom_Parabola)::DownCast(hCurve);
          gp_Parab parab = hParab->Parab();
          gp_Ax1   xaxis = parab.XAxis();
          gp_Ax1   yaxis = parab.YAxis();
          gp_Pnt   locat = parab.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = xaxis.Direction().X();
          data[4] = xaxis.Direction().Y();
          data[5] = xaxis.Direction().Z();
          data[6] = yaxis.Direction().X();
          data[7] = yaxis.Direction().Y();
          data[8] = yaxis.Direction().Z();
          data[9] = parab.Focal();
          *rvec   = data;
        }
        break;

      case HYPERBOLA:
        data = (double *) EG_alloc(11*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Hyperbola (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Hyperbola) hHypr =
            Handle(Geom_Hyperbola)::DownCast(hCurve);
          gp_Hypr hypr  = hHypr->Hypr();
          gp_Ax1  xaxis = hypr.XAxis();
          gp_Ax1  yaxis = hypr.YAxis();
          gp_Pnt  locat = hypr.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = xaxis.Direction().X();
          data[ 4] = xaxis.Direction().Y();
          data[ 5] = xaxis.Direction().Z();
          data[ 6] = yaxis.Direction().X();
          data[ 7] = yaxis.Direction().Y();
          data[ 8] = yaxis.Direction().Z();
          data[ 9] = hypr.MajorRadius();
          data[10] = hypr.MinorRadius();
          *rvec    = data;
        }
        break;

      case TRIMMED:
        data = (double *) EG_alloc(2*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Trimmed Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_TrimmedCurve) hTrim =
            Handle(Geom_TrimmedCurve)::DownCast(hCurve);
          data[0] = hTrim->FirstParameter();
          data[1] = hTrim->LastParameter();
          *rvec   = data;
        }
        break;

      case BEZIER:
        ints = (int *) EG_alloc(3*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Bezier Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_BezierCurve) hBezier =
            Handle(Geom_BezierCurve)::DownCast(hCurve);
          int rational = 0;
          if (hBezier->IsRational()) rational = 1;
          ints[0] = rational*2;
          if (hBezier->IsPeriodic()) ints[0] |= 4;
          ints[1] = hBezier->Degree();
          ints[2] = hBezier->NbPoles();
          len = ints[2]*3;
          if (rational == 1) len += ints[2];
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc on Bezier Data (EG_getGeometry)!\n");
            EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= ints[2]; i++, len+=3) {
            gp_Pnt P    = hBezier->Pole(i);
            data[len  ] = P.X();
            data[len+1] = P.Y();
            data[len+2] = P.Z();
          }
          if (rational == 1)
            for (i = 1; i <= ints[2]; i++,len++) data[len] = hBezier->Weight(i);
          *ivec = ints;
          *rvec = data;
        }
        break;

      case BSPLINE:
        ints = (int *) EG_alloc(4*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on BSpline Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_BSplineCurve) hBSpline =
            Handle(Geom_BSplineCurve)::DownCast(hCurve);
          int rational = 0;
          if (hBSpline->IsRational()) rational = 1;
          ints[0] = rational*2;
          if (hBSpline->IsPeriodic()) ints[0] |= 4;
          ints[1] = hBSpline->Degree();
          ints[2] = hBSpline->NbPoles();
          ints[3] = 0;
          for (i = 1; i <= hBSpline->NbKnots(); i++)
            ints[3] += hBSpline->Multiplicity(i);
          len = ints[3] + ints[2]*3;
          if (rational == 1) len += ints[2];
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc BSpline Data (EG_getGeometry)!\n");
            EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= hBSpline->NbKnots(); i++) {
            int km = hBSpline->Multiplicity(i);
            for (j = 1; j <= km; j++, len++) data[len] = hBSpline->Knot(i);
          }
          for (i = 1; i <= ints[2]; i++, len+=3) {
            gp_Pnt P    = hBSpline->Pole(i);
            data[len  ] = P.X();
            data[len+1] = P.Y();
            data[len+2] = P.Z();
          }
          if (rational == 1)
            for (i = 1; i <= ints[2]; i++,len++)
              data[len] = hBSpline->Weight(i);
          *ivec = ints;
          *rvec = data;
        }
        break;

      case OFFSET:
        data = (double *) EG_alloc(4*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Offset Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_OffsetCurve) hOffst =
            Handle(Geom_OffsetCurve)::DownCast(hCurve);
          gp_Dir direct = hOffst->Direction();
          data[0] = direct.X();
          data[1] = direct.Y();
          data[2] = direct.Z();
          data[3] = hOffst->Offset();
          *rvec   = data;
        }
        break;
    }

  } else {

    egadsSurface *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurf = psurf->handle;
    *refGeom = psurf->ref;
    if ((rvec == NULL) || (ivec == NULL)) return EGADS_SUCCESS;

    switch (geom->mtype) {

      case PLANE:
        data = (double *) EG_alloc(9*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Plane (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Plane) hPlane = Handle(Geom_Plane)::DownCast(hSurf);
          gp_Pln plane = hPlane->Pln();
          gp_Pnt locat = plane.Location();
          gp_Ax1 xaxis = plane.XAxis();
          gp_Ax1 yaxis = plane.YAxis();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = xaxis.Direction().X();
          data[4] = xaxis.Direction().Y();
          data[5] = xaxis.Direction().Z();
          data[6] = yaxis.Direction().X();
          data[7] = yaxis.Direction().Y();
          data[8] = yaxis.Direction().Z();
          *rvec   = data;
        }
        break;

      case SPHERICAL:
        data = (double *) EG_alloc(10*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Sphere (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_SphericalSurface) hSphere =
            Handle(Geom_SphericalSurface)::DownCast(hSurf);
          gp_Sphere sphere = hSphere->Sphere();
          gp_Pnt locat     = sphere.Location();
          gp_Ax1 xaxis     = sphere.XAxis();
          gp_Ax1 yaxis     = sphere.YAxis();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = xaxis.Direction().X();
          data[4] = xaxis.Direction().Y();
          data[5] = xaxis.Direction().Z();
          data[6] = yaxis.Direction().X();
          data[7] = yaxis.Direction().Y();
          data[8] = yaxis.Direction().Z();
          data[9] = sphere.Radius();
          if (!sphere.Direct()) data[9] = -data[9];
          *rvec   = data;
        }
        break;

      case CONICAL:
        data = (double *) EG_alloc(14*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Conical (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_ConicalSurface) hCone =
            Handle(Geom_ConicalSurface)::DownCast(hSurf);
          gp_Cone cone  = hCone->Cone();
          gp_Ax3  axes  = cone.Position();
          gp_Pnt  locat = cone.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = axes.XDirection().X();
          data[ 4] = axes.XDirection().Y();
          data[ 5] = axes.XDirection().Z();
          data[ 6] = axes.YDirection().X();
          data[ 7] = axes.YDirection().Y();
          data[ 8] = axes.YDirection().Z();
          data[ 9] = axes.Direction().X();
          data[10] = axes.Direction().Y();
          data[11] = axes.Direction().Z();
          data[12] = cone.SemiAngle();
          data[13] = cone.RefRadius();
          *rvec    = data;
        }
        break;

      case CYLINDRICAL:
        data = (double *) EG_alloc(13*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Cylinder (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_CylindricalSurface) hCyl =
            Handle(Geom_CylindricalSurface)::DownCast(hSurf);
          gp_Cylinder cyl   = hCyl->Cylinder();
          gp_Ax3      axes  = cyl.Position();
          gp_Pnt      locat = cyl.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = axes.XDirection().X();
          data[ 4] = axes.XDirection().Y();
          data[ 5] = axes.XDirection().Z();
          data[ 6] = axes.YDirection().X();
          data[ 7] = axes.YDirection().Y();
          data[ 8] = axes.YDirection().Z();
          data[ 9] = axes.Direction().X();
          data[10] = axes.Direction().Y();
          data[11] = axes.Direction().Z();
          data[12] = cyl.Radius();
          *rvec    = data;
        }
        break;

      case TOROIDAL:
        data = (double *) EG_alloc(14*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Cylinder (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_ToroidalSurface) hTorus =
            Handle(Geom_ToroidalSurface)::DownCast(hSurf);
          gp_Torus torus = hTorus->Torus();
          gp_Ax3   axes  = torus.Position();
          gp_Pnt   locat = torus.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = axes.XDirection().X();
          data[ 4] = axes.XDirection().Y();
          data[ 5] = axes.XDirection().Z();
          data[ 6] = axes.YDirection().X();
          data[ 7] = axes.YDirection().Y();
          data[ 8] = axes.YDirection().Z();
          data[ 9] = axes.Direction().X();
          data[10] = axes.Direction().Y();
          data[11] = axes.Direction().Z();
          data[12] = torus.MajorRadius();
          data[13] = torus.MinorRadius();
          *rvec    = data;
        }
        break;

      case BEZIER:
        ints = (int *) EG_alloc(5*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Bezier header (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_BezierSurface) hBezier =
            Handle(Geom_BezierSurface)::DownCast(hSurf);
          int rational = 0;
          if (hBezier->IsURational()) rational = 1;
          if (hBezier->IsVRational()) rational = 1;
          ints[0] = rational*2;
          if (hBezier->IsUPeriodic()) ints[0] |=  4;
          if (hBezier->IsVPeriodic()) ints[0] |=  8;
          ints[1] = hBezier->UDegree();
          ints[3] = hBezier->VDegree();
          ints[2] = hBezier->NbUPoles();
          ints[4] = hBezier->NbVPoles();
          int nCP = ints[2];
          nCP    *= ints[4];
          len     = nCP*3;
          if (rational == 1) len += nCP;
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc on Bezier Surf (EG_getGeometry)!\n");
            EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (j = 1; j <= ints[4]; j++)
            for (i = 1; i <= ints[2]; i++, len+=3) {
              gp_Pnt P    = hBezier->Pole(i, j);
              data[len  ] = P.X();
              data[len+1] = P.Y();
              data[len+2] = P.Z();
            }
          if (rational == 1)
            for (j = 1; j <= ints[4]; j++)
              for (i = 1; i <= ints[2]; i++,len++)
                data[len] = hBezier->Weight(i, j);
          *ivec = ints;
          *rvec = data;
        }
        break;

      case BSPLINE:
        ints = (int *) EG_alloc(7*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc BSpline header (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_BSplineSurface) hBSpline =
            Handle(Geom_BSplineSurface)::DownCast(hSurf);
          int rational = 0;
          if (hBSpline->IsURational()) rational = 1;
          if (hBSpline->IsVRational()) rational = 1;
          ints[0] = rational*2;
          if (hBSpline->IsUPeriodic()) ints[0] |=  4;
          if (hBSpline->IsVPeriodic()) ints[0] |=  8;
          ints[1] = hBSpline->UDegree();
          ints[4] = hBSpline->VDegree();
          ints[2] = hBSpline->NbUPoles();
          ints[5] = hBSpline->NbVPoles();
          ints[3] = ints[6] = 0;
          for (i = 1; i <= hBSpline->NbUKnots(); i++)
            ints[3] += hBSpline->UMultiplicity(i);
          for (i = 1; i <= hBSpline->NbVKnots(); i++)
            ints[6] += hBSpline->VMultiplicity(i);
          int nCP = ints[2];
          nCP    *= ints[5];
          len     = ints[3] + ints[6] + nCP*3;
          if (rational == 1) len += nCP;
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc BSpline Surf (EG_getGeometry)!\n");
            EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= hBSpline->NbUKnots(); i++) {
            int km = hBSpline->UMultiplicity(i);
            for (j = 1; j <= km; j++, len++) data[len] = hBSpline->UKnot(i);
          }
          for (i = 1; i <= hBSpline->NbVKnots(); i++) {
            int km = hBSpline->VMultiplicity(i);
            for (j = 1; j <= km; j++, len++) data[len] = hBSpline->VKnot(i);
          }
          for (j = 1; j <= ints[5]; j++)
            for (i = 1; i <= ints[2]; i++, len+=3) {
              gp_Pnt P    = hBSpline->Pole(i, j);
              data[len  ] = P.X();
              data[len+1] = P.Y();
              data[len+2] = P.Z();
            }
          if (rational == 1)
            for (j = 1; j <= ints[5]; j++)
              for (i = 1; i <= ints[2]; i++,len++)
                data[len] = hBSpline->Weight(i, j);
          *ivec = ints;
          *rvec = data;
        }
        break;

      case OFFSET:
        data = (double *) EG_alloc(sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc Offset Surface (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_OffsetSurface) hOffst =
            Handle(Geom_OffsetSurface)::DownCast(hSurf);
          data[0] = hOffst->Offset();
          *rvec   = data;
        }
        break;

      case TRIMMED:
        data = (double *) EG_alloc(4*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc Trimmed Surface (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_RectangularTrimmedSurface) hTrim =
            Handle(Geom_RectangularTrimmedSurface)::DownCast(hSurf);
          hTrim->Bounds(data[0], data[1], data[2], data[3]);
          *rvec = data;
        }
        break;

      case EXTRUSION:
        data = (double *) EG_alloc(3*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc Linear Extrusion (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_SurfaceOfLinearExtrusion) hSLExtr =
            Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(hSurf);
          gp_Dir direct = hSLExtr->Direction();
          data[0] = direct.X();
          data[1] = direct.Y();
          data[2] = direct.Z();
          *rvec   = data;
        }
        break;

      case REVOLUTION:
        data = (double *) EG_alloc(6*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc Revolved Surface (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_SurfaceOfRevolution) hSORev =
            Handle(Geom_SurfaceOfRevolution)::DownCast(hSurf);
          gp_Pnt locat     = hSORev->Location();
          gp_Ax1 axis      = hSORev->Axis();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = axis.Direction().X();
          data[4] = axis.Direction().Y();
          data[5] = axis.Direction().Z();
          *rvec   = data;
        }
        break;

    }

  }

  return EGADS_SUCCESS;
}


// Surreal version of getGeometry
int
EG_getGeometry(const ego geom, int *oclass, int *mtype, ego *refGeom,
               int **ivec, SurrealS<1> **rvec)
{
  int    stat;
  double *rinfo;

  stat = EG_getGeometry(geom, oclass, mtype, refGeom, ivec, &rinfo);
  EG_free(rinfo);
  if (stat != EGADS_SUCCESS) return stat;

  stat = EG_getGeometry_dot(geom, rvec);
  if (stat != EGADS_SUCCESS) return stat;

  return EGADS_SUCCESS;
}


namespace {
void EG_normalizeDir_dot(int dim, SurrealS<1>* data_dot)
{
  SurrealS<1> mag = 0;

  for (int i = 0; i < dim; i++)
    mag += data_dot[i]*data_dot[i];

  mag = sqrt(mag);

  if (mag == 0.0) {
    printf(" EGADS Warning: Zero magnitude (EG_setGeometry_dot)!\n");
    return;
  }

  for (int i = 0; i < dim; i++)
    data_dot[i] /= mag;
}


void EG_ortho_dot(SurrealS<1>* dirx, SurrealS<1>* diry)
{
  SurrealS<1> dirz[3];

  CROSS(dirz, dirx, diry);
  CROSS(diry, dirz, dirx);

  EG_normalizeDir_dot(3, diry);
}

  
/* taken from gp_Ax3::gp_Ax3(const gp_Pnt& P, const gp_Dir& N, const gp_Dir& Vx)
 * and gp_XYZ::CrossCross (const gp_XYZ& Coord1, const gp_XYZ& Coord2)
 */
void EG_CrossCross_dot(const SurrealS<1>* dirz, SurrealS<1>* dirx)
{
  SurrealS<1>& x1 = dirx[0];
  SurrealS<1>& y1 = dirx[1];
  SurrealS<1>& z1 = dirx[2];

  const SurrealS<1>& x2 = dirz[0];
  const SurrealS<1>& y2 = dirz[1];
  const SurrealS<1>& z2 = dirz[2];

  SurrealS<1> X =  y2 * (x1 * y2 - y1 * x2) -
                   z2 * (z1 * x2 - x1 * z2);
  SurrealS<1> Y =  z2 * (y1 * z2 - z1 * y2) -
                   x2 * (x1 * y2 - y1 * x2);
              z1 = x2 * (z1 * x2 - x1 * z2) -
                   y2 * (y1 * z2 - z1 * y2);
  y1 = Y;
  x1 = X;

  EG_normalizeDir_dot(3, dirx);
}
} //namespace


int
EG_setGeometry_dot(egObject *obj, int oclass, int mtype,
                   /*@null@*/ const int *ivec,
                   /*@null@*/ const double *rvec,
                   /*@null@*/ const double *rvec_dot)
{
  static
  const char *classType[27] = {"CONTEXT", "TRANSFORM", "TESSELLATION",
                               "NIL", "EMPTY", "REFERENCE", "", "",
                               "", "", "PCURVE", "CURVE", "SURFACE", "",
                               "", "", "", "", "", "", "NODE",
                               "EGDE", "LOOP", "FACE", "SHELL",
                               "BODY", "MODEL"};
  static
  const char *nodeType[1] = {""};
  static
  const char *edgeType[6] = {"", "ONENODE", "TWONODE",
                             "", "", "DEGENERATE"};
  static
  const char *curvType[9] = {"LINE", "CIRCLE", "ELLIPSE", "PARABOLA",
                             "HYPERBOLA", "TRIMMED", "BEZIER", "BSPLINE",
                             "OFFSET"};
  static
  const char *surfType[11] = {"PLANE", "SPHERICAL", "CYLINDER", "REVOLUTION",
                              "TOROIDAL", "TRIMMED" , "BEZIER", "BSPLINE",
                              "OFFSET", "CONICAL", "EXTRUSION"};
  const char  **geomType, **eType;
  int         i, j, len, ilen, outLevel, stat;
  double      *data, scale;
  SurrealS<1> **data_dot;
  egObject    *geom, *object;

  if  (obj == NULL)               return EGADS_NULLOBJ;
  if  (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((obj->oclass != PCURVE)  && (obj->oclass != CURVE) &&
      (obj->oclass != SURFACE) && (obj->oclass != EDGE)  &&
      (obj->oclass != LOOP)    && (obj->oclass != FACE)  &&
      (obj->oclass != NODE)    && (obj->oclass != SHELL) &&
      (obj->oclass != BODY))      return EGADS_NOTGEOM;
  if  (obj->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(obj))        return EGADS_CNTXTHRD;

  outLevel = EG_outLevel(obj);

  /* Node section */

  if (obj->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) obj->blind;
    if (rvec_dot == NULL) {
      if ((oclass != NODE) && (oclass != 0)) {
        if (outLevel > 0) {
          printf(" EGADS Error: Object Node is not expected %s (EG_setGeometry_dot)!\n",
                 classType[oclass]);
        }
        return EGADS_GEOMERR;
      }
      pnode->filled = 0;
      for (i = 0; i < 3; i++) {
        pnode->xyz_dot[i] = 0;
      }
    } else {
      if (oclass != NODE) {
        if (outLevel > 0) {
          printf(" EGADS Error: Object Node is not expected %s (EG_setGeometry_dot)!\n",
                 classType[oclass]);
        }
        return EGADS_GEOMERR;
      }
      pnode->filled = 1;
      for (i = 0; i < 3; i++) {
        pnode->xyz_dot[i].value() = pnode->xyz[i];
        pnode->xyz_dot[i].deriv() = rvec_dot[i];
      }

      /* check consistency in the data */
      stat = EGADS_SUCCESS;
      scale = 0.0;
      for (i = 0; i < 3; i++) {
        if (fabs(pnode->xyz[i]) > scale) scale = fabs(pnode->xyz[i]);
      }
      if (scale == 0.0) scale = 1.0;
      for (i = 0; i < 3; i++)
        if (fabs(pnode->xyz[i] - rvec[i]) > 1.e-14*scale) {
          stat++;
          break;
        }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0) {
          printf(" EGADS Error: Inconsistent Node geometry data! (EG_setGeometry_dot)\n");
          for (i = 0; i < 3; i++)
            printf("     data[%d] %lf : %lf\n", i, pnode->xyz[i], rvec[i]);
        }
        return EGADS_GEOMERR;
      }
    }
    return EGADS_SUCCESS;
  }

  /* Edge section */

  if (obj->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) obj->blind;
    if (rvec_dot == NULL) {
      if ((oclass != EDGE) && (oclass != 0)) {
        if (outLevel > 0) {
          printf(" EGADS Error: Object Edge is not expected %s (EG_setGeometry_dot)!\n",
                 classType[oclass]);
        }
        return EGADS_GEOMERR;
      }
      pedge->filled = 0;
      for (i = 0; i < 2; i++) {
        pedge->trange_dot[i] = 0;
      }
    } else {
      if (oclass != EDGE) {
        if (outLevel > 0) {
          printf(" EGADS Error: Object Edge is not expected %s (EG_setGeometry_dot)!\n",
                 classType[oclass]);
        }
        return EGADS_GEOMERR;
      }
      if (mtype != obj->mtype) {
        if (outLevel > 0) {
          printf(" EGADS Error: Object Edge %s is not expected Edge %s (EG_setGeometry_dot)!\n",
               edgeType[obj->mtype], edgeType[mtype]);
        }
        return EGADS_GEOMERR;
      }

      pedge->filled = 1;
      for (i = 0; i < 2; i++) {
        pedge->trange_dot[i].value() = pedge->trange[i];
        pedge->trange_dot[i].deriv() = rvec_dot[i];
      }

      /* check consistency in the data */
      stat = EGADS_SUCCESS;
      scale = 0.0;
      for (i = 0; i < 2; i++) {
        if (fabs(pedge->trange[i]) > scale) scale = fabs(pedge->trange[i]);
      }
      if (scale == 0.0) scale = 1.0;
      for (i = 0; i < 2; i++)
        if (fabs(pedge->trange[i] - rvec[i]) > 1.e-14*scale) {
          stat++;
          break;
        }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0) {
          printf(" EGADS Error: Inconsistent Edge t-range data! (EG_setGeometry_dot)\n");
          for (i = 0; i < 2; i++)
            printf("     data[%d] %lf : %lf\n", i, pedge->trange[i], rvec[i]);
        }
        return EGADS_GEOMERR;
      }
    }
    return EGADS_SUCCESS;
  }

  /* Loop section -- clear all dots */

  if (obj->oclass == LOOP) {
    if (oclass != LOOP && oclass != 0) {
      if (outLevel > 0) {
        printf(" EGADS Error: Object Loop is not expected %s (EG_setGeometry_dot)!\n",
               classType[oclass]);
      }
      return EGADS_GEOMERR;
    }
    if (rvec_dot != NULL) {
      if (outLevel > 0) {
        printf(" EGADS Error: Cannot set non-NULL sensitivity on Loop (EG_setGeometry_dot)!\n");
      }
      return EGADS_NODATA;
    }
    egadsLoop *ploop = (egadsLoop *) obj->blind;

    for (i = 0; i < ploop->nedges; i++) {
      stat = EG_setGeometry_dot(ploop->edges[i], EDGE, 0, NULL, NULL, NULL);
      if (stat != EGADS_SUCCESS) return stat;

      egadsEdge *pedge = (egadsEdge *) ploop->edges[i]->blind;

      stat = EG_setGeometry_dot(pedge->curve, CURVE, pedge->curve->mtype, NULL, NULL, NULL);
      if (stat != EGADS_SUCCESS) return stat;

      stat = EG_setGeometry_dot(pedge->nodes[0], NODE, 0, NULL, NULL, NULL);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_setGeometry_dot(pedge->nodes[1], NODE, 0, NULL, NULL, NULL);
      if (stat != EGADS_SUCCESS) return stat;
    }

    return EGADS_SUCCESS;
  }

  /* Face section -- clear all dots */

  if (obj->oclass == FACE) {
    if (oclass != FACE && oclass != 0) {
      if (outLevel > 0) {
        printf(" EGADS Error: Object Face is not expected %s (EG_setGeometry_dot)!\n",
               classType[oclass]);
      }
      return EGADS_GEOMERR;
    }
    if (rvec_dot != NULL) {
      if (outLevel > 0) {
        printf(" EGADS Error: Cannot set non-NULL sensitivity on Face (EG_setGeometry_dot)!\n");
      }
      return EGADS_NODATA;
    }
    egadsFace *pface = (egadsFace *) obj->blind;

    for (i = 0; i < pface->nloops; i++) {
      stat = EG_setGeometry_dot(pface->loops[i], LOOP, 0, NULL, NULL, NULL);
      if (stat != EGADS_SUCCESS) return stat;
    }

    return EGADS_SUCCESS;
  }

  /* Shell section -- clear all dots */

  if (obj->oclass == SHELL) {
    if (oclass != SHELL && oclass != 0) {
      if (outLevel > 0) {
        printf(" EGADS Error: Object Shell is not expected %s (EG_setGeometry_dot)!\n",
               classType[oclass]);
      }
      return EGADS_GEOMERR;
    }
    if (rvec_dot != NULL) {
      if (outLevel > 0) {
        printf(" EGADS Error: Cannot set non-NULL sensitivity on Shell (EG_setGeometry_dot)!\n");
      }
      return EGADS_NODATA;
    }
    egadsShell *pshell = (egadsShell *) obj->blind;

    for (i = 0; i < pshell->nfaces; i++) {
      stat = EG_setGeometry_dot(pshell->faces[i], FACE, 0, NULL, NULL, NULL);
      if (stat != EGADS_SUCCESS) return stat;
    }

    return EGADS_SUCCESS;
  }

  /* Body section -- clear all dots */

  if (obj->oclass == BODY) {
    if (oclass != BODY && oclass != 0) {
      if (outLevel > 0) {
        printf(" EGADS Error: Object Body is not expected %s (EG_setGeometry_dot)!\n",
               classType[oclass]);
      }
      return EGADS_GEOMERR;
    }
    if (rvec_dot != NULL) {
      if (outLevel > 0) {
        printf(" EGADS Error: Cannot set non-NULL sensitivity on Body (EG_setGeometry_dot)!\n");
      }
      return EGADS_NODATA;
    }
    egadsBody *pbody = (egadsBody *) obj->blind;

    /* Nodes */
    for (i = 0; i < pbody->nodes.map.Extent(); i++) {
      object = pbody->nodes.objs[i];
      egadsNode *pnode = (egadsNode *) object->blind;
      if (pnode == NULL) continue;
      pnode->filled = 0;
      for (j = 0; j < 3; j++) {
        pnode->xyz_dot[j] = 0;
      }
    }

    /* Curves through Edges */
    for (i = 0; i < pbody->edges.map.Extent(); i++) {
      object = pbody->edges.objs[i];
      egadsEdge *pedge = (egadsEdge *) object->blind;
      if (pedge == NULL) continue;
      object = pedge->curve;
      egadsCurve *lgeom = (egadsCurve *) object->blind;
      if (lgeom == NULL) continue;
      if (lgeom->data_dot != NULL) {
        EG_free(lgeom->data_dot);
        lgeom->data_dot = NULL;
      }
    }

    /* PCurves through Loops */
    for (i = 0; i < pbody->loops.map.Extent(); i++) {
      object = pbody->loops.objs[i];
      egadsLoop *ploop = (egadsLoop *) object->blind;
      if (ploop == NULL) continue;
      object = ploop->surface;
      if (object == NULL) continue;
      if (object->mtype == PLANE) continue;
      for (j = 0; j < ploop->nedges; j++) {
        object = ploop->edges[ploop->nedges+j];
        if (object == NULL) continue;
        egadsPCurve *lgeom = (egadsPCurve *) object->blind;
        if (lgeom == NULL) continue;
        if (lgeom->data_dot != NULL) {
          EG_free(lgeom->data_dot);
          lgeom->data_dot = NULL;
        }
      }
    }

    /* Surfaces through Faces */
    for (i = 0; i < pbody->faces.map.Extent(); i++) {
      object = pbody->faces.objs[i];
      egadsFace *pface = (egadsFace *) object->blind;
      if (pface == NULL) continue;
      object = pface->surface;
      egadsSurface *lgeom = (egadsSurface *) object->blind;
      if (lgeom == NULL) continue;
      if (lgeom->data_dot != NULL) {
        EG_free(lgeom->data_dot);
        lgeom->data_dot = NULL;
      }
    }

    return EGADS_SUCCESS;
  }

  /* Geometry section */
  geom = obj;

  if ((rvec_dot != NULL) || (oclass != 0)) {
    if ((oclass != geom->oclass) || (mtype != geom->mtype)) {
      if (outLevel > 0) {
        int emtype = mtype-1;
        if (oclass == SURFACE) {
          eType = surfType;
        } else if ((oclass == CURVE) || (oclass == PCURVE)) {
          eType = curvType;
        } else if (oclass == NODE) {
          eType = nodeType;
          emtype = 0;
        } else {
          printf(" EGADS Error: Unexpected oclass %s (EG_setGeometry_dot)!\n",
                 classType[oclass]);
          return EGADS_GEOMERR;
        }
        int gmtype = geom->mtype;
        if (geom->oclass == SURFACE) {
          geomType = surfType;
        } else if ((geom->oclass == CURVE) || (geom->oclass == PCURVE)) {
          geomType = curvType;
        } else if (oclass == NODE) {
          geomType = nodeType;
          gmtype = 0;
        } else {
          printf(" EGADS Error: Unexpected geom oclass %s (EG_setGeometry_dot)!\n",
                 classType[geom->oclass]);
          return EGADS_GEOMERR;
        }
        printf(" EGADS Error: Object %s %s is not expected %s %s (EG_setGeometry_dot)!\n",
               classType[geom->oclass], geomType[gmtype], classType[oclass],
               eType[emtype]);
      }
      return EGADS_GEOMERR;
    }
  }

  if (geom->oclass == PCURVE) {

#ifndef REQUIRE_PCURVE_SENSITIVITIES
    printf(" EGADS Error: PCURVE sensitivities are not yet supported (EG_setGeometry_dot)!\n");
    return EGADS_GEOMERR;
#endif

    egadsPCurve *lgeom = (egadsPCurve *) geom->blind;
    len                =  lgeom->dataLen;
    data               =  lgeom->data;
    data_dot           = &lgeom->data_dot;

    if (rvec_dot != NULL) {
      if ((geom->mtype == BEZIER) || (geom->mtype == BSPLINE)) {
        if (ivec == NULL) return EGADS_NODATA;
        ilen = geom->mtype == BEZIER ? 3 : 4;
        for (i = 0; i < ilen; i++) {
          if (lgeom->header[i] != ivec[i]) {
            if (outLevel > 0) {
              printf(" EGADS Error: Inconsistent %s %s ivec[%d](%d) != %d (EG_setGeometry_dot)!\n",
                     classType[geom->oclass], curvType[geom->mtype-1], i, ivec[i],
                     lgeom->header[i]);
            }
            return EGADS_GEOMERR;
          }
        }
      }
    }

  } else if (geom->oclass == CURVE) {

    egadsCurve *lgeom = (egadsCurve *) geom->blind;
    len               =  lgeom->dataLen;
    data              =  lgeom->data;
    data_dot          = &lgeom->data_dot;

    if (rvec_dot != NULL) {
      if ((geom->mtype == BEZIER) || (geom->mtype == BSPLINE)) {
        if (ivec == NULL) return EGADS_NODATA;
        ilen = geom->mtype == BEZIER ? 3 : 4;
        for (i = 0; i < ilen; i++) {
          if (lgeom->header[i] != ivec[i]) {
            if (outLevel > 0) {
              printf(" EGADS Error: Inconsistent %s %s ivec[%d](%d) != %d (EG_setGeometry_dot)!\n",
                     classType[geom->oclass], curvType[geom->mtype-1], i, ivec[i],
                     lgeom->header[i]);
            }
            return EGADS_GEOMERR;
          }
        }
      }
    }

  } else {

    egadsSurface *lgeom = (egadsSurface *) geom->blind;
    len                 =  lgeom->dataLen;
    data                =  lgeom->data;
    data_dot            = &lgeom->data_dot;

    if (rvec_dot != NULL) {
      if ((geom->mtype == BEZIER) || (geom->mtype == BSPLINE)) {
        if (ivec == NULL) return EGADS_NODATA;
        ilen = geom->mtype == BEZIER ? 5 : 7;
        for (i = 0; i < ilen; i++) {
          if (lgeom->header[i] != ivec[i]) {
            if (outLevel > 0) {
              printf(" EGADS Error: Inconsistent %s %s ivec[%d](%d) != %d (EG_setGeometry_dot)!\n",
                     classType[geom->oclass], curvType[geom->mtype-1], i, ivec[i],
                     lgeom->header[i]);
            }
            return EGADS_GEOMERR;
          }
        }
      }
    }
  }

  if (rvec_dot == NULL) {
    if (*data_dot != NULL) EG_free(*data_dot);
    *data_dot = NULL;
    return EGADS_SUCCESS;
  } else {
    if (*data_dot == NULL) {
      if (len == 0) return EGADS_GEOMERR;
      *data_dot = (SurrealS<1> *) EG_alloc(len*sizeof(SurrealS<1>));
      if (*data_dot == NULL) return EGADS_MALLOC;
    }
    for (i = 0; i < len; i++) {
      (*data_dot)[i].value() = rvec[i];
      (*data_dot)[i].deriv() = rvec_dot[i];
    }
  }

  if (geom->oclass == PCURVE) {

    switch (geom->mtype) {
    case LINE:
      {
        //gp_Dir2d dirl(data[2], data[3]);

        (*data_dot)[2].value() = rvec[2];
        (*data_dot)[3].value() = rvec[3];

        EG_normalizeDir_dot(2, &(*data_dot)[2]);
      }
      break;

    case CIRCLE:
    case ELLIPSE:
    case PARABOLA:
    case HYPERBOLA:
      {
        //gp_Dir2d dirx(data[2], data[3]);
        //gp_Dir2d diry(data[4], data[5]);

        (*data_dot)[2].value() = rvec[2];
        (*data_dot)[3].value() = rvec[3];

        (*data_dot)[4].value() = rvec[4];
        (*data_dot)[5].value() = rvec[5];

        EG_normalizeDir_dot(2, &(*data_dot)[2]);
        EG_normalizeDir_dot(2, &(*data_dot)[4]);
      }
      break;
    }

  } else if (geom->oclass == CURVE) {

    switch (geom->mtype) {
    case LINE:
      {
        //gp_Dir dirl(data[3], data[4], data[5]);

        (*data_dot)[3].value() = rvec[3];
        (*data_dot)[4].value() = rvec[4];
        (*data_dot)[5].value() = rvec[5];

        EG_normalizeDir_dot(3, &(*data_dot)[3]);
      }
      break;

    case CIRCLE:
    case ELLIPSE:
    case PARABOLA:
    case HYPERBOLA:
      {
        //gp_Dir dirx(data[3], data[4], data[5]);
        //gp_Dir diry(data[6], data[7], data[8]);
        //gp_Dir dirz = dirx.Crossed(diry);
        //gp_Ax2 axi2(pntc, dirz, dirx);

        (*data_dot)[3].value() = rvec[3];
        (*data_dot)[4].value() = rvec[4];
        (*data_dot)[5].value() = rvec[5];

        (*data_dot)[6].value() = rvec[6];
        (*data_dot)[7].value() = rvec[7];
        (*data_dot)[8].value() = rvec[8];

        EG_normalizeDir_dot(3, &(*data_dot)[3]);
        EG_normalizeDir_dot(3, &(*data_dot)[6]);

        EG_ortho_dot(&(*data_dot)[3], &(*data_dot)[6]);
      }
      break;

    case OFFSET:
      {
        //gp_Dir dir(data[0], data[1], data[2]);

        (*data_dot)[0].value() = rvec[0];
        (*data_dot)[1].value() = rvec[1];
        (*data_dot)[2].value() = rvec[2];

        EG_normalizeDir_dot(3, &(*data_dot)[0]);
      }
      break;
    }

  } else { // geom->oclass == SURFACE

    switch (geom->mtype) {

    case PLANE:
      {
        //gp_Dir dirx(data[3], data[4], data[5]);
        //gp_Dir diry(data[6], data[7], data[8]);
        //gp_Dir dirz = dirx.Crossed(diry);
        //gp_Ax3 axi3(pntp, dirz, dirx);
        SurrealS<1> dirz[3];

        (*data_dot)[3].value() = rvec[3];
        (*data_dot)[4].value() = rvec[4];
        (*data_dot)[5].value() = rvec[5];

        (*data_dot)[6].value() = rvec[6];
        (*data_dot)[7].value() = rvec[7];
        (*data_dot)[8].value() = rvec[8];

        EG_normalizeDir_dot(3, &(*data_dot)[3]);
        EG_normalizeDir_dot(3, &(*data_dot)[6]);

        CROSS(dirz, (&(*data_dot)[3]), (&(*data_dot)[6]));
        EG_normalizeDir_dot(3, dirz);

        EG_CrossCross_dot(dirz, (&(*data_dot)[3]));
        CROSS((&(*data_dot)[6]), dirz, (&(*data_dot)[3]));
        EG_normalizeDir_dot(3, (&(*data_dot)[6]));
      }
      break;

    case SPHERICAL:
      {
        //gp_Pnt pnts(data[0], data[1], data[2]);
        //gp_Dir dirx(data[3], data[4], data[5]);
        //gp_Dir diry(data[6], data[7], data[8]);
        //gp_Dir dirz = dirx.Crossed(diry);
        //if (data[9] < 0.0) dirz.SetCoord(-dirz.X(), -dirz.Y(), -dirz.Z());
        //gp_Ax3 axi3(pnts, dirz);

        (*data_dot)[3].value() = rvec[3];
        (*data_dot)[4].value() = rvec[4];
        (*data_dot)[5].value() = rvec[5];

        (*data_dot)[6].value() = rvec[6];
        (*data_dot)[7].value() = rvec[7];
        (*data_dot)[8].value() = rvec[8];

        (*data_dot)[9].value() = rvec[9];

        SurrealS<1> *vxdir = &(*data_dot)[3];
        SurrealS<1> *vydir = &(*data_dot)[6];
        SurrealS<1> axis[3];

        EG_normalizeDir_dot(3, vxdir);
        EG_normalizeDir_dot(3, vydir);

        CROSS(axis, vxdir, vydir);
        EG_normalizeDir_dot(3, axis);

        gp_Pnt pnts(rvec[0], rvec[1], rvec[2]);
        gp_Dir dirx(rvec[3], rvec[4], rvec[5]);
        gp_Dir diry(rvec[6], rvec[7], rvec[8]);
        gp_Dir dirz = dirx.Crossed(diry);
        if (rvec[9] < 0.0) {
          axis[0] = -axis[0];
          axis[1] = -axis[1];
          axis[2] = -axis[2];
          dirz.SetCoord(-dirz.X(), -dirz.Y(), -dirz.Z());
        }
        gp_Ax3 axi3(pnts, dirz);
        // axi3.SetXDirection(dirx);
        // axi3.SetYDirection(diry);

        // axi3.SetXDirection(dirx);
        //
        // Standard_Boolean direct = Direct();
        // vxdir = axis.Direction().CrossCrossed (Vx, axis.Direction());
        // if (direct) { vydir = axis.Direction().Crossed(vxdir); }
        // else        { vydir = vxdir.Crossed(axis.Direction()); }

        // Calling SetXDirection can only modify 'direct'. Otherwise it is irrelevant.
        axi3.SetXDirection(dirx);

        // axi3.SetYDirection(diry);
        //
        // Standard_Boolean direct = Direct();
        // vxdir = Vy.Crossed (axis.Direction());
        // vydir = (axis.Direction()).Crossed (vxdir);
        // if (!direct) { vxdir.Reverse(); }

        Standard_Boolean direct = axi3.Direct();
        //axi3.SetYDirection(diry);
        CROSS(vxdir, vydir, axis);
        EG_normalizeDir_dot(3, vxdir);
        CROSS(vydir, axis, vxdir);
        EG_normalizeDir_dot(3, vydir);

        if (!direct) {
          vxdir[3] = -vxdir[3];
          vxdir[4] = -vxdir[4];
          vxdir[5] = -vxdir[5];
        }

        // abs of radius
        (*data_dot)[9] = fabs((*data_dot)[9]);
      }
      break;

    case CONICAL:
    case CYLINDRICAL:
    case TOROIDAL:
      {
        //gp_Dir dirx(data[3], data[4],  data[5]);
        //gp_Dir diry(data[6], data[7],  data[8]);
        //gp_Dir dirz(data[9], data[10], data[11]);

        (*data_dot)[3].value() = rvec[3];
        (*data_dot)[4].value() = rvec[4];
        (*data_dot)[5].value() = rvec[5];

        (*data_dot)[6].value() = rvec[6];
        (*data_dot)[7].value() = rvec[7];
        (*data_dot)[8].value() = rvec[8];

        (*data_dot)[ 9].value() = rvec[ 9];
        (*data_dot)[10].value() = rvec[10];
        (*data_dot)[11].value() = rvec[11];

        SurrealS<1> *vxdir = &(*data_dot)[3];
        SurrealS<1> *vydir = &(*data_dot)[6];
        SurrealS<1> *axis  = &(*data_dot)[9];

        EG_normalizeDir_dot(3, vxdir);
        EG_normalizeDir_dot(3, vydir);
        EG_normalizeDir_dot(3, axis );

        gp_Pnt pnts(rvec[0], rvec[1], rvec[2]);
        gp_Dir dirx(rvec[3], rvec[4], rvec[5]);
        gp_Dir diry(rvec[6], rvec[7], rvec[8]);
        gp_Dir dirz(rvec[9], rvec[10], rvec[11]);
        gp_Ax3 axi3(pnts, dirz);
        // axi3.SetXDirection(dirx);
        // axi3.SetYDirection(diry);

        // body of axi3.SetXDirection(dirx):
        //
        // Standard_Boolean direct = Direct();
        // vxdir = axis.Direction().CrossCrossed (Vx, axis.Direction());
        // if (direct) { vydir = axis.Direction().Crossed(vxdir); }
        // else        { vydir = vxdir.Crossed(axis.Direction()); }

        // Calling SetXDirection can only modify 'direct'. Otherwise it is irrelevant.
        axi3.SetXDirection(dirx);

        // body of axi3.SetYDirection(diry):
        //
        // Standard_Boolean direct = Direct();
        // vxdir = Vy.Crossed (axis.Direction());
        // vydir = (axis.Direction()).Crossed (vxdir);
        // if (!direct) { vxdir.Reverse(); }

        Standard_Boolean direct = axi3.Direct();
        //axi3.SetYDirection(diry);
        CROSS(vxdir, vydir, axis);
        EG_normalizeDir_dot(3, vxdir);
        CROSS(vydir, axis, vxdir);
        EG_normalizeDir_dot(3, vydir);

        if (!direct) {
          vxdir[3] = -vxdir[3];
          vxdir[4] = -vxdir[4];
          vxdir[5] = -vxdir[5];
        }
      }
      break;

    case EXTRUSION:
      {
        // gp_Dir dir(data[0], data[1],  data[2]);

        (*data_dot)[0].value() = rvec[0];
        (*data_dot)[1].value() = rvec[1];
        (*data_dot)[2].value() = rvec[2];

        EG_normalizeDir_dot(3, &(*data_dot)[0]);
      }
      break;

    case REVOLUTION:
      {
        // gp_Dir dir(data[3], data[4], data[5]);

        (*data_dot)[3].value() = rvec[3];
        (*data_dot)[4].value() = rvec[4];
        (*data_dot)[5].value() = rvec[5];

        EG_normalizeDir_dot(3, &(*data_dot)[3]);
      }
      break;
    }
  }

  /* check consistency in the data */
  stat = EGADS_SUCCESS;
  scale = 0.0;
  for (i = 0; i < len; i++) {
    if (fabs(data[i]) > scale) scale = fabs(data[i]);
  }
  if (scale == 0.0) scale = 1.0;
  for (i = 0; i < len; i++)
    if (fabs(data[i] - (*data_dot)[i].value()) > 1.e-14*scale) {
      stat++;
      break;
    }
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0) {
      printf(" EGADS Error: Inconsistent geometry data! (EG_setGeometry_dot)\n");
      if (outLevel > 1) {
        for (i = 0; i < len; i++)
          printf("     data[%d] %lf : %lf\n", i, data[i], (*data_dot)[i].value());
      }
    }
    return EGADS_GEOMERR;
  }

  return EGADS_SUCCESS;
}


int
EG_setGeometry_dot(egObject *geom, int oclass, int mtype,
                   /*@null@*/ const int *ivec, SurrealS<1> *data_dot)
{
  int    stat, i, len;
  double *rvec=NULL, *rvec_dot=NULL;

  if (geom == NULL)               return EGADS_NULLOBJ;
  if (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (geom->blind == NULL)        return EGADS_NODATA;

  if (data_dot != NULL) {
    if (geom->oclass == NODE) {
      len = 3;
    } else if (geom->oclass == EDGE) {
      len = 2;
    } else {
      EG_getGeometryLen(geom, &i, &len);
      if (len == 0) return EGADS_GEOMERR;
    }

    rvec     = (double*)EG_alloc(len*sizeof(double));
    rvec_dot = (double*)EG_alloc(len*sizeof(double));

    for (i = 0; i < len; i++) {
      rvec[i]     = data_dot[i].value();
      rvec_dot[i] = data_dot[i].deriv();
    }
  }

  stat = EG_setGeometry_dot(geom, oclass, mtype, ivec, rvec, rvec_dot);
  EG_free(rvec);
  EG_free(rvec_dot);

  return stat;
}


int
EG_hasGeometry_dot(const egObject *obj)
{
  int            i, found;
  const egObject *geom, *object;

  if  (obj == NULL)               return EGADS_NULLOBJ;
  if  (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((obj->oclass != PCURVE)  && (obj->oclass != CURVE) &&
      (obj->oclass != SURFACE) && (obj->oclass != EDGE)  &&
      (obj->oclass != FACE)    && (obj->oclass != NODE)  &&
      (obj->oclass != LOOP)    && (obj->oclass != SHELL) &&
      (obj->oclass != BODY))
                                  return EGADS_NOTGEOM;
  if  (obj->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(obj))        return EGADS_CNTXTHRD;

  /* Body section -- check that all geometry has dot information */

  if (obj->oclass == BODY) {
    egadsBody *pbody = (egadsBody *) obj->blind;

    /* Nodes */
    for (i = 0; i < pbody->nodes.map.Extent(); i++) {
      object = pbody->nodes.objs[i];
      egadsNode *pnode = (egadsNode *) object->blind;
      if (pnode->filled == 0) return EGADS_NOTFOUND;
    }

    /* Curves through Edges */
    for (i = 0; i < pbody->edges.map.Extent(); i++) {
      object = pbody->edges.objs[i];
      egadsEdge *pedge = (egadsEdge *) object->blind;
      if (pedge == NULL) continue;
      if (pedge->filled == 0) return EGADS_NOTFOUND;
      object = pedge->curve;
      egadsCurve *lgeom = (egadsCurve *) object->blind;
      if (lgeom == NULL) continue;
      if (lgeom->data_dot == NULL) return EGADS_NOTFOUND;
    }

#ifdef REQUIRE_PCURVE_SENSITIVITIES
    /* PCurves through Loops */
    for (i = 0; i < pbody->loops.map.Extent(); i++) {
      object = pbody->loops.objs[i];
      egadsLoop *ploop = (egadsLoop *) object->blind;
      if (ploop == NULL) continue;
      object = ploop->surface;
      if (object == NULL) continue;
      if (object->mtype == PLANE) continue;
      for (int j = 0; j < ploop->nedges; j++) {
        object = ploop->edges[ploop->nedges+j];
        if (object == NULL) continue;
        egadsPCurve *lgeom = (egadsPCurve *) object->blind;
        if (lgeom == NULL) continue;
        if (lgeom->data_dot == NULL) return EGADS_NOTFOUND;
      }
    }
#endif

    /* Surfaces through Faces */
    for (i = 0; i < pbody->faces.map.Extent(); i++) {
      object = pbody->faces.objs[i];
      egadsFace *pface = (egadsFace *) object->blind;
      if (pface == NULL) continue;
      object = pface->surface;
      egadsSurface *lgeom = (egadsSurface *) object->blind;
      if (lgeom == NULL) continue;
      if (lgeom->data_dot == NULL) return EGADS_NOTFOUND;
    }

    return EGADS_SUCCESS;
  }

  /* Node section */

  if (obj->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) obj->blind;
    if (pnode->filled == 0) return EGADS_NOTFOUND;
    return EGADS_SUCCESS;
  }

  /* Edge section -- check curve and nodes */

  if (obj->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) obj->blind;
    if (pedge->filled == 0) return EGADS_NOTFOUND;
    if (obj->mtype == DEGENERATE) {
      return EG_hasGeometry_dot(pedge->nodes[0]);
    }
    geom = pedge->curve;
    if (geom->blind == NULL) return EGADS_NODATA;

    if ( (EG_hasGeometry_dot(pedge->curve)    == EGADS_SUCCESS) &&
         (EG_hasGeometry_dot(pedge->nodes[0]) == EGADS_SUCCESS) &&
         (EG_hasGeometry_dot(pedge->nodes[1]) == EGADS_SUCCESS) )
      return EGADS_SUCCESS;

    return EGADS_NOTFOUND;
  }

  /* Loop section -- check edges and pcurves*/

  if (obj->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) obj->blind;

    found = 1;
    for (i = 0; i < ploop->nedges; i++)
      if (EG_hasGeometry_dot(ploop->edges[i]) != EGADS_SUCCESS) found = 0;

#ifdef REQUIRE_PCURVE_SENSITIVITIES
    /* PCurves through Loops */
    if ((ploop->surface != NULL) && (ploop->surface->mtype != PLANE)) {
      for (int j = 0; j < ploop->nedges; j++) {
        object = ploop->edges[ploop->nedges+j];
        if (object == NULL) continue;
        egadsPCurve *lgeom = (egadsPCurve *) object->blind;
        if (lgeom == NULL) continue;
        if (lgeom->data_dot == NULL) found = 0;
      }
    }
#endif

    if (found == 0) return EGADS_NOTFOUND;
    return EGADS_SUCCESS;
  }

  /* Face section -- check loops and surface */

  if (obj->oclass == FACE) {
    egadsFace *pface = (egadsFace *) obj->blind;

    found = 1;
    for (i = 0; i < pface->nloops; i++)
      if (EG_hasGeometry_dot(pface->loops[i]) != EGADS_SUCCESS) found = 0;
    if (found == 0) return EGADS_NOTFOUND;

    if (EG_hasGeometry_dot(pface->surface) != EGADS_SUCCESS)
      return EGADS_NOTFOUND;

    return EGADS_SUCCESS;
  }

  /* Shell section -- check faces */

  if (obj->oclass == SHELL) {
    egadsShell *pshell = (egadsShell *) obj->blind;

    found = 1;
    for (i = 0; i < pshell->nfaces; i++)
      if (EG_hasGeometry_dot(pshell->faces[i]) != EGADS_SUCCESS) found = 0;
    if (found == 0) return EGADS_NOTFOUND;

    return EGADS_SUCCESS;
  }

  /* Geometry section */

  geom = obj;

  if (geom->oclass == PCURVE) {

    egadsPCurve *lgeom = (egadsPCurve *) geom->blind;
    if (lgeom->data_dot == NULL) return EGADS_NOTFOUND;

  } else if (geom->oclass == CURVE) {

    egadsCurve *lgeom = (egadsCurve *) geom->blind;
    if (lgeom->data_dot == NULL) return EGADS_NOTFOUND;

  } else {

    egadsSurface *lgeom = (egadsSurface *) geom->blind;
    if (lgeom->data_dot == NULL) return EGADS_NOTFOUND;

  }

  return EGADS_SUCCESS;
}


int
EG_getGeometry_dot(const egObject *obj, double **rvec, double **rvec_dot)
{
  int            i, len;
  SurrealS<1>    *rdata_dot;
  const egObject *geom;

  if  (obj == NULL)               return EGADS_NULLOBJ;
  if  (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((obj->oclass != PCURVE)  && (obj->oclass != CURVE) &&
      (obj->oclass != SURFACE) && (obj->oclass != EDGE)  &&
      (obj->oclass != FACE)    && (obj->oclass != NODE))
                                  return EGADS_NOTGEOM;
  if  (obj->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(obj))        return EGADS_CNTXTHRD;

  if  (rvec     == NULL)          return EGADS_NODATA;
  if  (rvec_dot == NULL)          return EGADS_NODATA;
  *rvec     = NULL;
  *rvec_dot = NULL;

  /* Node section */

  if (obj->oclass == NODE) {

    egadsNode *pnode = (egadsNode *) obj->blind;
    if (pnode->filled == 0)       return EGADS_NODATA;
    rdata_dot = pnode->xyz_dot;
    len = 3;

  } else {

    /* Geometry section */

    geom = obj;
    if (obj->oclass == EDGE) {
      egadsEdge *pedge = (egadsEdge *) obj->blind;
      geom = pedge->curve;
      if (geom->blind == NULL) return EGADS_NODATA;
    } else if (obj->oclass == FACE) {
      egadsFace *pface = (egadsFace *) obj->blind;
      geom = pface->surface;
      if (geom->blind == NULL) return EGADS_NODATA;
    }

    if (geom->oclass == PCURVE) {

      egadsPCurve *lgeom = (egadsPCurve *) geom->blind;
      if (lgeom->data_dot == NULL) return EGADS_NODATA;
      len                = lgeom->dataLen;
      rdata_dot          = lgeom->data_dot;

    } else if (geom->oclass == CURVE) {

      egadsCurve *lgeom = (egadsCurve *) geom->blind;
      if (lgeom->data_dot == NULL) return EGADS_NODATA;
      len               = lgeom->dataLen;
      rdata_dot         = lgeom->data_dot;

    } else {

      egadsSurface *lgeom = (egadsSurface *) geom->blind;
      if (lgeom->data_dot == NULL) return EGADS_NODATA;
      len                 = lgeom->dataLen;
      rdata_dot           = lgeom->data_dot;

    }
  }

  if (len == 0) return EGADS_GEOMERR;
  *rvec     = (double *) EG_alloc(len*sizeof(double));
  if (*rvec == NULL) return EGADS_MALLOC;
  *rvec_dot = (double *) EG_alloc(len*sizeof(double));
  if (*rvec_dot == NULL) {
    EG_free(*rvec);
    *rvec = NULL;
    return EGADS_MALLOC;
  }

  for (i = 0; i < len; i++) {
    (*rvec)[i]     = rdata_dot[i].value();
    (*rvec_dot)[i] = rdata_dot[i].deriv();
  }

  return EGADS_SUCCESS;
}


int
EG_getGeometry_dot(const egObject *obj, SurrealS<1> **data_dot)
{
  int    stat;
  int    i, len;
  double *rvec, *rvec_dot;

  if (obj == NULL)      return EGADS_NULLOBJ;
  if (data_dot == NULL) return EGADS_NODATA;

  if (obj->oclass == NODE) {
    len = 3;
  } else {
    EG_getGeometryLen(obj, &i, &len);
    if (len == 0) return EGADS_GEOMERR;
  }

  (*data_dot) = (SurrealS<1> *) EG_alloc(len*sizeof(SurrealS<1>));

  stat = EG_getGeometry_dot(obj, &rvec, &rvec_dot);

  for (i = 0; i < len; i++) {
    (*data_dot)[i].value() = rvec[i];
    (*data_dot)[i].deriv() = rvec_dot[i];
  }

  EG_free(rvec);
  EG_free(rvec_dot);

  return stat;
}


namespace { /* private to this file */

template<class T>
class egMat;

template<class T>
class egTrsf;

//! This is a strippped down version of eg_XYZ
template<class T>
class egXYZ
{
public:

  //! creates an XYZ uninitialized
  egXYZ() { }

  //! creates an XYZ with given coordinates
  egXYZ(const T& X, const T& Y, const T& Z) : x(X), y(Y), z(Z) { }

  //! Returns the X coordinate
  const T& X() const { return x; }

  //! Returns the Y coordinate
  const T& Y() const { return y; }

  //! Returns the Z coordinate
  const T& Z() const { return z; }

  //! computes sqrt(X*X + Y*Y + Z*Z) where X, Y and Z are the three coordinates of this XYZ object.
  T Modulus() const
  {
    return sqrt(x * x + y * y + z * z);
  }

  egXYZ& operator= (const egXYZ& Other)
  {
    x = Other.x;
    y = Other.y;
    z = Other.z;
    return *this;
  }

  //! computes the scalar product between <me> and Other
  T Dot (const egXYZ& Other) const
  {
    return x * Other.x + y * Other.y + z * Other.z;
  }

  //! <me>.X() = <me>.X() + Other.X()
  //! <me>.Y() = <me>.Y() + Other.Y()
  //! <me>.Z() = <me>.Z() + Other.Z()
  void Add (const egXYZ& Other)
  {
    x += Other.x;
    y += Other.y;
    z += Other.z;
  }

  //! <me>.X() = <me>.X() - Other.X()
  //! <me>.Y() = <me>.Y() - Other.Y()
  //! <me>.Z() = <me>.Z() - Other.Z()
  void Subtract (const egXYZ& Right)
  {
    x -= Right.x;
    y -= Right.y;
    z -= Right.z;
  }
  void operator -= (const egXYZ& Right)
  {
    Subtract(Right);
  }
  //! new.X() = <me>.X() + Other.X()
  //! new.Y() = <me>.Y() + Other.Y()
  //! new.Z() = <me>.Z() + Other.Z()
  egXYZ Added (const egXYZ& Other) const
  {
    return egXYZ(x + Other.x,y + Other.y,z + Other.z);
  }
  egXYZ operator + (const egXYZ& Other) const
  {
    return Added(Other);
  }

  //! <me>.X() = <me>.X() * Scalar;
  //! <me>.Y() = <me>.Y() * Scalar;
  //! <me>.Z() = <me>.Z() * Scalar;
  void Multiply (const T& Scalar)
  {
    x *= Scalar;
    y *= Scalar;
    z *= Scalar;
  }

  egXYZ Multiplied (const T& Scalar) const
  {
    return egXYZ(x * Scalar,y * Scalar,z * Scalar);
  }
  egXYZ operator * (const T& Scalar) const
  {
    return Multiplied(Scalar);
  }

  //! <me> = Matrix * <me>
  void Multiply (const egMat<T>& Matrix)
  {
    T Xresult = Matrix.matrix[0][0] * x + Matrix.matrix[0][1] * y +
                Matrix.matrix[0][2] * z;
    T Yresult = Matrix.matrix[1][0] * x + Matrix.matrix[1][1] * y +
                Matrix.matrix[1][2] * z;
    z         = Matrix.matrix[2][0] * x + Matrix.matrix[2][1] * y +
                Matrix.matrix[2][2] * z;
    x         = Xresult;
    y         = Yresult;
  }

  //! <me>.X() = <me>.X()/ <me>.Modulus()
  //! <me>.Y() = <me>.Y()/ <me>.Modulus()
  //! <me>.Z() = <me>.Z()/ <me>.Modulus()
  //! Raised if <me>.Modulus() <= Resolution from gp
  void Normalize()
  {
    T D = Modulus();
    Standard_ConstructionError_Raise_if (D <= gp::Resolution(),
                                         "egXYZ::Normalize() - vector has zero norm");
    x = x / D;  y = y / D;  z = z / D;
  }

private:
  T x;
  T y;
  T z;
};


//! Stripped down version of gp_Pnt
//! Stores reference of values for in-place update
template<class T>
class egPnt
{
public:

  //! Creates a direction with its 3 cartesian coordinates.
  //! Raises ConstructionError if Sqrt(Xv*Xv + Yv*Yv + Zv*Zv) <= Resolution
  //! Modification of the direction's coordinates
  //! If Sqrt (X*X + Y*Y + Z*Z) <= Resolution from gp where
  //! X, Y ,Z are the new coordinates it is not possible to
  //! construct the direction and the method raises the
  //! exception ConstructionError.
  egPnt(T& X, T& Y, T& Z) : x(X), y(Y), z(Z) {}

  //! computes Sqrt (X*X + Y*Y + Z*Z) where X, Y and Z are the three coordinates of this XYZ object.
  T Modulus() const
  {
    return sqrt (x * x + y * y + z * z);
  }

  //! <me>.X() = -<me>.X()
  //! <me>.Y() = -<me>.Y()
  //! <me>.Z() = -<me>.Z()
  void Reverse()
  {
    x = -x;
    y = -y;
    z = -z;
  }

  //! <me> = Matrix * <me>
  void Multiply (const egMat<T>& Matrix)
  {
    T Xresult = Matrix.matrix[0][0] * x + Matrix.matrix[0][1] * y +
                Matrix.matrix[0][2] * z;
    T Yresult = Matrix.matrix[1][0] * x + Matrix.matrix[1][1] * y +
                Matrix.matrix[1][2] * z;
    z         = Matrix.matrix[2][0] * x + Matrix.matrix[2][1] * y +
                Matrix.matrix[2][2] * z;
    x         = Xresult;
    y         = Yresult;
  }

  //! <me>.X() = <me>.X() / Scalar;
  //! <me>.Y() = <me>.Y() / Scalar;
  //! <me>.Z() = <me>.Z() / Scalar;
  void Divide (const T& Scalar)
  {
    x /= Scalar;
    y /= Scalar;
    z /= Scalar;
  }

  void Transform (const egTrsf<T>& TT)
  {
    TT.Transforms(x, y, z);
  }

private:
  T& x;
  T& y;
  T& z;
};


//! Stripped down version of gp_Dir
template<class T>
class egDir
{
public:

  //! Creates a direction with its 3 cartesian coordinates.
  //! Raises ConstructionError if Sqrt(Xv*Xv + Yv*Yv + Zv*Zv) <= Resolution
  //! Modification of the direction's coordinates
  //! If Sqrt (X*X + Y*Y + Z*Z) <= Resolution from gp where
  //! X, Y ,Z are the new coordinates it is not possible to
  //! construct the direction and the method raises the
  //! exception ConstructionError.
  egDir(T& X, T& Y, T& Z) : coord(X,Y,Z) {}

  void Transform (const egTrsf<T>& TT)
  {
    coord.Multiply (TT.HVectorialPart());
    T D = coord.Modulus();
    coord.Divide(D);
    if (TT.ScaleFactor() < 0.0) { coord.Reverse(); }
  }

private:
  egPnt<T> coord;
};


//! Describes a three column, three row matrix. This sort of
//! object is used in various vectorial or matrix computations.
// Stripped down implementation of gp_Mat
template<class T>
class egMat
{
public:

#define Mat00 matrix[0][0]
#define Mat01 matrix[0][1]
#define Mat02 matrix[0][2]
#define Mat10 matrix[1][0]
#define Mat11 matrix[1][1]
#define Mat12 matrix[1][2]
#define Mat20 matrix[2][0]
#define Mat21 matrix[2][1]
#define Mat22 matrix[2][2]

  //! creates  a matrix with null coefficients.
  egMat() {}

  egMat(const T& a11, const T& a12, const T& a13,
        const T& a21, const T& a22, const T& a23,
        const T& a31, const T& a32, const T& a33)
  {
    Mat00 = a11;
    Mat01 = a12;
    Mat02 = a13;
    Mat10 = a21;
    Mat11 = a22;
    Mat12 = a23;
    Mat20 = a31;
    Mat21 = a32;
    Mat22 = a33;
  }

  //! Creates a matrix.
  //! Col1, Col2, Col3 are the 3 columns of the matrix.
  Standard_EXPORT egMat(const egXYZ<T>& Col1, const egXYZ<T>& Col2,
                        const egXYZ<T>& Col3)
  {
    SetCols(Col1, Col2, Col3);
  }

  //! Assigns the number triples Col1, Col2, Col3 to the three
  //! columns of this matrix.
  Standard_EXPORT void SetCols(const egXYZ<T>& Col1, const egXYZ<T>& Col2,
                               const egXYZ<T>& Col3)
  {
    Mat00 = Col1.X(); Mat10 = Col1.Y(); Mat20 = Col1.Z();
    Mat01 = Col2.X(); Mat11 = Col2.Y(); Mat21 = Col2.Z();
    Mat02 = Col3.X(); Mat12 = Col3.Y(); Mat22 = Col3.Z();
  }

  //! Assigns the number triples Row1, Row2, Row3 to the three
  //! rows of this matrix.
  Standard_EXPORT void SetRows(const egXYZ<T>& Row1, const egXYZ<T>& Row2,
                               const egXYZ<T>& Row3)
  {
    Mat00 = Row1.X(); Mat01 = Row1.Y(); Mat02 = Row1.Z();
    Mat10 = Row2.X(); Mat11 = Row2.Y(); Mat12 = Row2.Z();
    Mat20 = Row3.X(); Mat21 = Row3.Y(); Mat22 = Row3.Z();
  }

  //! Returns the column of Col index.
  egXYZ<T> Column (const Standard_Integer Col) const
  {
    if (Col == 1) return egXYZ<T> (Mat00,Mat10,Mat20);
    if (Col == 2) return egXYZ<T> (Mat01,Mat11,Mat21);
    else          return egXYZ<T> (Mat02,Mat12,Mat22);
  }

  //! returns the row of Row index.
  egXYZ<T> Row (const Standard_Integer Row) const
  {
    if (Row == 1) return egXYZ<T> (Mat00,Mat01,Mat02);
    if (Row == 2) return egXYZ<T> (Mat10,Mat11,Mat12);
    else          return egXYZ<T> (Mat20,Mat21,Mat22);
  }

  //! Computes the determinant of the matrix.
  T Determinant() const
  {
    return Mat00 * (Mat11 * Mat22 - Mat21 * Mat12) -
           Mat01 * (Mat10 * Mat22 - Mat20 * Mat12) +
           Mat02 * (Mat10 * Mat21 - Mat20 * Mat11);
  }

  void Divide (const T& Scalar)
  {
    T val = Scalar;
    if (val < 0) val = - val;
    Standard_ConstructionError_Raise_if
        (val <= gp::Resolution(),"egMat : Divide by 0");
    T UnSurScalar = 1.0 / Scalar;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        matrix[i][j] *= UnSurScalar;
  }

#undef Mat00
#undef Mat01
#undef Mat02
#undef Mat10
#undef Mat11
#undef Mat12
#undef Mat20
#undef Mat21
#undef Mat22

friend class egXYZ<T>;
friend class egPnt<T>;

private:
  T matrix[3][3];
};


// Stripped down implementation of gp_Trsf
template<class T>
class egTrsf
{
public:

  //! Returns the identity transformation.
  egTrsf() : scale(1.0) {}

  //! Sets the coefficients  of the transformation.  The
  //! transformation  of the  point  x,y,z is  the point
  //! x',y',z' with :
  //!
  //! x' = a11 x + a12 y + a13 z + a14
  //! y' = a21 x + a22 y + a23 z + a24
  //! z' = a31 x + a32 y + a33 z + a34
  //!
  //! The method Value(i,j) will return aij.
  //! Raises ConstructionError if the determinant of  the aij is null.
  //! The matrix is orthogonalized before future using.
  Standard_EXPORT void SetValues (const T& a11, const T& a12, const T& a13,
                                  const T& a14, const T& a21, const T& a22,
                                  const T& a23, const T& a24, const T& a31,
                                  const T& a32, const T& a33, const T& a34)
  {
    egXYZ<T> col1(a11,a21,a31);
    egXYZ<T> col2(a12,a22,a32);
    egXYZ<T> col3(a13,a23,a33);
    egXYZ<T> col4(a14,a24,a34);
    // compute the determinant
    egMat<T> M(col1,col2,col3);
    T s = M.Determinant();
    T As = s;
    if (As < 0) As = - As;
    Standard_ConstructionError_Raise_if
      (As < gp::Resolution(),"gp_Trsf::SetValues, null determinant");
    if (s > 0)
      s = pow(s,1./3.);
    else
      s = -pow(-s,1./3.);
    M.Divide(s);

    scale = s;

    matrix = M;
    Orthogonalize();

    loc = col4;
  }

  //! Computes the homogeneous vectorial part of the transformation.
  //! It is a 3*3 matrix which doesn't include the scale factor.
  //! In other words, the vectorial part of this transformation is equal
  //! to its homogeneous vectorial part, multiplied by the scale factor.
  //! The coefficients of this matrix must be multiplied by the
  //! scale factor to obtain the coefficients of the transformation.
  const egMat<T>& HVectorialPart() const { return matrix; }

  //! Returns the scale factor.
  const T& ScaleFactor() const { return scale; }

  //! Transform a coordinate
  void Transforms (T& X, T& Y, T& Z) const
  {
    egXYZ<T> Triplet (X, Y, Z);
    Triplet.Multiply (matrix);
    Triplet.Multiply (scale);
    Triplet.Add(loc);
    X = Triplet.X();
    Y = Triplet.Y();
    Z = Triplet.Z();
  }

  //! Transformation of a triplet XYZ with a Trsf
  void Transforms (egXYZ<T>& Coord) const
  {
    Coord.Multiply (matrix);
    Coord.Multiply (scale);
    Coord.Add(loc);
  }

protected:

  //! Makes orthogonalization of "matrix"
  void Orthogonalize()
  {
    egMat<T> aTM(matrix);

    egXYZ<T> aV1 = aTM.Column(1);
    egXYZ<T> aV2 = aTM.Column(2);
    egXYZ<T> aV3 = aTM.Column(3);

    aV1.Normalize();

    aV2 -= aV1*(aV2.Dot(aV1));
    aV2.Normalize();

    aV3 -= aV1*(aV3.Dot(aV1)) + aV2*(aV3.Dot(aV2));
    aV3.Normalize();

    aTM.SetCols(aV1, aV2, aV3);

    aV1 = aTM.Row(1);
    aV2 = aTM.Row(2);
    aV3 = aTM.Row(3);

    aV1.Normalize();

    aV2 -= aV1*(aV2.Dot(aV1));
    aV2.Normalize();

    aV3 -= aV1*(aV3.Dot(aV1)) + aV2*(aV3.Dot(aV2));
    aV3.Normalize();

    aTM.SetRows(aV1, aV2, aV3);

    matrix = aTM;
  }

private:

  T scale;
  egMat<T> matrix;
  egXYZ<T> loc;

};
} // namespace


int
EG_copyGeometry_dot(const egObject *obj, /*@null@*/ const double *xform,
                    const egTrsf< SurrealS<1> >& form, egObject *copy)
{
  typedef SurrealS<1> T;
  int            stat, *ints = NULL, i, j, len, outLevel;
  double         *cdata, scale;
  const egObject *geom1;
  ego            geom2;
#ifdef REQUIRE_PCURVE_SENSITIVITIES
  int            j;
  const egObject *object1;
  ego            object2;
#endif
  SurrealS<1>    *rdata_dot, *cdata_dot;

  if (obj == NULL)                 return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if ((obj->oclass != PCURVE)  && (obj->oclass != CURVE) &&
      (obj->oclass != SURFACE) && (obj->oclass != EDGE)  &&
      (obj->oclass != FACE)    && (obj->oclass != NODE)  &&
      (obj->oclass != LOOP)    && (obj->oclass != BODY))
                                   return EGADS_NOTGEOM;
  if (obj->blind == NULL)          return EGADS_NODATA;
  if (EG_sameThread(obj))          return EGADS_CNTXTHRD;

  if (copy == NULL)                return EGADS_NULLOBJ;
  if (copy->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if (copy->blind == NULL)         return EGADS_NODATA;
  if (copy->oclass != obj->oclass) return EGADS_GEOMERR;
  if (EG_sameThread(copy))         return EGADS_CNTXTHRD;

  if (EG_context(obj) != EG_context(copy)) return EGADS_MIXCNTX;

  outLevel = EG_outLevel(obj);

  /* Body section -- recursively copy over sensitivities of all geometry */

  if (obj->oclass == BODY) {
    egadsBody *pbody1 = (egadsBody *) obj->blind;
    egadsBody *pbody2 = (egadsBody *) copy->blind;

    if (pbody1->nodes.map.Extent() !=
        pbody2->nodes.map.Extent() ) {
      if (outLevel > 0)
        printf(" EGADS Error: Node count mismatch (EG_copyGeometry_dot)!\n");
      return EGADS_TOPOERR;
    }

    if (pbody1->edges.map.Extent() !=
        pbody2->edges.map.Extent() ) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge count mismatch (EG_copyGeometry_dot)!\n");
      return EGADS_TOPOERR;
    }

    if (pbody1->loops.map.Extent() !=
        pbody2->loops.map.Extent() ) {
      if (outLevel > 0)
        printf(" EGADS Error: Loop count mismatch (EG_copyGeometry_dot)!\n");
      return EGADS_TOPOERR;
    }

    if (pbody1->faces.map.Extent() !=
        pbody2->faces.map.Extent() ) {
      if (outLevel > 0)
        printf(" EGADS Error: Face count mismatch (EG_copyGeometry_dot)!\n");
      return EGADS_TOPOERR;
    }

    /* Nodes */
    for (i = 0; i < pbody1->nodes.map.Extent(); i++) {
      stat = EG_copyGeometry_dot( pbody1->nodes.objs[i], xform, form,
                                  pbody2->nodes.objs[i] );
      if (stat != EGADS_SUCCESS) return stat;
    }

    /* Curves through Edges */
    for (i = 0; i < pbody1->edges.map.Extent(); i++) {
      /* make sure both edges are the same mtype */
      if (pbody1->edges.objs[i]->mtype != pbody2->edges.objs[i]->mtype) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d mtype mismatch (EG_copyGeometry_dot)!\n", i+1);
        return EGADS_TOPOERR;
      }

      egadsEdge *pedge1 = (egadsEdge *) pbody1->edges.objs[i]->blind;
      egadsEdge *pedge2 = (egadsEdge *) pbody2->edges.objs[i]->blind;

      /* copy curve if the the edge is not degenerate */
      if (pbody1->edges.objs[i]->mtype != DEGENERATE) {
        stat = EG_copyGeometry_dot( pedge1->curve, xform, form, pedge2->curve );
        if (stat != EGADS_SUCCESS) return stat;
      }

      if (pedge1->filled == 0) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge without data_dot (EG_copyGeometry_dot)!\n");
        return EGADS_NODATA;
      }
      pedge2->filled = 1;
      for (j = 0; j < 2; j++)
        pedge2->trange_dot[j] = pedge1->trange_dot[j];
      if (pedge2->curve->mtype == LINE ||
          pedge2->curve->mtype == PARABOLA)
        for (j = 0; j < 2; j++)
          pedge2->trange_dot[j] *= fabs(form.ScaleFactor());

      if (pedge2->curve->mtype == OFFSET) {
        egadsCurve *pcurve = (egadsCurve *) pedge2->curve->blind;
        if (pcurve->ref->mtype == LINE ||
            pcurve->ref->mtype == PARABOLA)
          for (j = 0; j < 2; j++)
            pedge2->trange_dot[j] *= fabs(form.ScaleFactor());
      }


      /* check consistency in the data */
      stat = EGADS_SUCCESS;
      scale = 0.0;
      for (j = 0; j < 2; j++)
        if (fabs(pedge2->trange[j]) > scale) scale = fabs(pedge2->trange[j]);
      if (scale == 0.0) scale = 1.0;
      for (j = 0; j < 2; j++)
        if (fabs(pedge2->trange[j] - pedge2->trange_dot[j].value()) > 1.e-14*scale) {
          stat++;
          break;
        }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0) {
          printf(" EGADS Error: Inconsistent t-range data! (EG_copyGeometry_dot)\n");
          if (outLevel > 1) {
            for (j = 0; j < 2; j++)
              printf("     data[%d] %lf : %lf\n",
                     j, pedge2->trange[j], pedge2->trange_dot[j].value());
          }
        }
        return EGADS_GEOMERR;
      }
    }

#ifdef REQUIRE_PCURVE_SENSITIVITIES
    /* PCurves through Loops */
    for (i = 0; i < pbody1->loops.map.Extent(); i++) {
      object1 = pbody1->loops.objs[i];
      object2 = pbody2->loops.objs[i];
      egadsLoop *ploop1 = (egadsLoop *) object1->blind;
      egadsLoop *ploop2 = (egadsLoop *) object2->blind;
      if (ploop1 == NULL) {
        if (ploop2 != NULL)                 return EGADS_TOPOERR;
        continue;
      }
      if (ploop2 == NULL)                   return EGADS_TOPOERR;
      object1 = ploop1->surface;
      object2 = ploop2->surface;
      if (object1 == NULL) {
        if (object2 != NULL)                return EGADS_TOPOERR;
        continue;
      }
      if (object2 == NULL)                  return EGADS_TOPOERR;
      if (object1->mtype != object2->mtype) return EGADS_TOPOERR;
      if (object1->mtype == PLANE) continue;
      for (j = 0; j < ploop1->nedges; j++) {
        object1 = ploop1->edges[ploop1->nedges+j];
        object2 = ploop2->edges[ploop2->nedges+j];
        if (object1 == NULL) {
          if (object2 != NULL)              return EGADS_TOPOERR;
          continue;
        }
        if (object2 == NULL)                return EGADS_TOPOERR;
        //if (EG_hasGeometry_dot(object1) != EGADS_SUCCESS) continue;
        stat = EG_copyGeometry_dot( object1, NULL, form, object2 );
        if (stat != EGADS_SUCCESS)          return stat;
      }
    }
#endif

    /* Surfaces through Faces */
    for (i = 0; i < pbody1->faces.map.Extent(); i++) {

      egadsFace *pface1 = (egadsFace *) pbody1->faces.objs[i]->blind;
      egadsFace *pface2 = (egadsFace *) pbody2->faces.objs[i]->blind;

      stat = EG_copyGeometry_dot(pface1->surface, xform, form, pface2->surface);
      if (stat != EGADS_SUCCESS) return stat;
    }

    return EGADS_SUCCESS;
  }

  /* Node section */

  if (obj->oclass == NODE) {

    egadsNode *pnode1 = (egadsNode *) obj->blind;
    egadsNode *pnode2 = (egadsNode *) copy->blind;
    if (pnode1->filled == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Node without data_dot (EG_copyGeometry_dot)!\n");
      return EGADS_NODATA;
    }
    for (i = 0; i < 3; i++)
      pnode2->xyz_dot[i] = pnode1->xyz_dot[i];

    pnode2->filled = 1;
    if (xform != NULL) {
      egPnt<T> pos(pnode2->xyz_dot[0], pnode2->xyz_dot[1], pnode2->xyz_dot[2]);
      pos.Transform(form);
    }

    /* check consistency in the data */
    stat = EGADS_SUCCESS;
    scale = 0.0;
    for (i = 0; i < 3; i++) {
      if (fabs(pnode2->xyz[i]) > scale) scale = fabs(pnode2->xyz[i]);
    }
    if (scale == 0.0) scale = 1.0;
    for (i = 0; i < 3; i++)
      if (fabs(pnode2->xyz[i] - pnode2->xyz_dot[i].value()) > 1.e-14*scale) {
        stat++;
        break;
      }
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0) {
        printf(" EGADS Error: Inconsistent geometry data! (EG_copyGeometry_dot)\n");
        if (outLevel > 1) {
          for (i = 0; i < 3; i++)
            printf("     data[%d] %lf : %lf\n",
                   i, pnode2->xyz[i], pnode2->xyz_dot[i].value());
        }
      }
      return EGADS_GEOMERR;
    }

    return EGADS_SUCCESS;
  }

  /* Edge section -- recursively copy over curve and nodes */

  if (obj->oclass == EDGE) {
    egadsEdge *pedge1 = (egadsEdge *) obj->blind;
    egadsEdge *pedge2 = (egadsEdge *) copy->blind;

    /* make sure both edges are the same mtype */
    if (obj->mtype != copy->mtype) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge mtype mismatch (EG_copyGeometry_dot)!\n");
      return EGADS_TOPOERR;
    }

    if (pedge1->filled == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge without data_dot (EG_copyGeometry_dot)!\n");
      return EGADS_NODATA;
    }
    pedge2->filled = 1;
    for (i = 0; i < 2; i++)
      pedge2->trange_dot[i] = pedge1->trange_dot[i];
    if (pedge2->curve->mtype == LINE ||
        pedge2->curve->mtype == PARABOLA)
      for (i = 0; i < 2; i++)
        pedge2->trange_dot[i] *= fabs(form.ScaleFactor());

    if (pedge2->curve->mtype == OFFSET) {
      egadsCurve *pcurve = (egadsCurve *) pedge2->curve->blind;
      if (pcurve->ref->mtype == LINE ||
          pcurve->ref->mtype == PARABOLA)
        for (i = 0; i < 2; i++)
          pedge2->trange_dot[i] *= fabs(form.ScaleFactor());
    }

    /* check consistency in the data */
    stat = EGADS_SUCCESS;
    scale = 0.0;
    for (i = 0; i < 2; i++) {
      if (fabs(pedge2->trange[i]) > scale) scale = fabs(pedge2->trange[i]);
    }
    if (scale == 0.0) scale = 1.0;
    for (i = 0; i < 2; i++)
      if (fabs(pedge2->trange[i] - pedge2->trange_dot[i].value()) > 1.e-14*scale) {
        stat++;
        break;
      }
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0) {
        printf(" EGADS Error: Inconsistent t-range data! (EG_copyGeometry_dot)\n");
        if (outLevel > 1) {
          for (i = 0; i < 2; i++)
            printf("     data[%d] %lf : %lf\n",
                   i, pedge2->trange[i], pedge2->trange_dot[i].value());
        }
      }
      return EGADS_GEOMERR;
    }


    if (obj->mtype == DEGENERATE) {

      /* Node */
      stat = EG_copyGeometry_dot(pedge1->nodes[0], xform, form, pedge2->nodes[0]);
      if (stat != EGADS_SUCCESS) return stat;

    } else {

      /* Curve */
      stat = EG_copyGeometry_dot(pedge1->curve, xform, form, pedge2->curve);
      if (stat != EGADS_SUCCESS) return stat;

      /* Nodes */
      stat = EG_copyGeometry_dot(pedge1->nodes[0], xform, form, pedge2->nodes[0]);
      if (stat != EGADS_SUCCESS) return stat;

      if (obj->mtype == TWONODE) {

        stat = EG_copyGeometry_dot(pedge1->nodes[1], xform, form, pedge2->nodes[1]);
        if (stat != EGADS_SUCCESS) return stat;
      }
    }

    return EGADS_SUCCESS;
  }

  /* Loop section -- recursively copy over edges (not the surface) */

  if (obj->oclass == LOOP) {

    egadsLoop *ploop1 = (egadsLoop *) obj->blind;
    egadsLoop *ploop2 = (egadsLoop *) copy->blind;

    if (ploop1->nedges != ploop2->nedges) {
      if (outLevel > 0)
        printf(" EGADS Error: obj and copy Loop edge count mismatch (EG_copyGeometry_dot)!\n");
      return EGADS_TOPOERR;
    }

    for (i = 0; i < ploop1->nedges; i++) {
      /* Edges */
      stat = EG_copyGeometry_dot(ploop1->edges[i], xform, form, ploop2->edges[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }

    /* surface if it exists */
//    if (ploop1->surface != NULL) {
//      if (ploop2->surface == NULL) return EGADS_GEOMERR;
//
//      stat = EG_copyGeometry_dot( ploop1->surface, xform, form, ploop2->surface );
//      if (stat != EGADS_SUCCESS) return stat;
//    }

    return EGADS_SUCCESS;
  }

  /* Face section -- recursively copy over loops and surface */

  if (obj->oclass == FACE) {
    egadsFace *pface1 = (egadsFace *) obj->blind;
    egadsFace *pface2 = (egadsFace *) copy->blind;

    /* Loops */
    for (i = 0; i < pface1->nloops; i++) {
      stat = EG_copyGeometry_dot(pface1->loops[i], xform, form, pface2->loops[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }

    if (pface1->surface != NULL) {
      if (pface2->surface == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: copy Surface is NULL (EG_copyGeometry_dot)!\n");
        return EGADS_GEOMERR;
      }

      stat = EG_copyGeometry_dot(pface1->surface, xform, form, pface2->surface);
      if (stat != EGADS_SUCCESS) return stat;
    }

    return EGADS_SUCCESS;
  }

  geom1 = obj;
  geom2 = copy;

  if (geom1->oclass == PCURVE) {

    if (xform != NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot transform PCurve (EG_copyGeometry_dot)!\n");
      return EGADS_CONSTERR;
    }

#ifndef REQUIRE_PCURVE_SENSITIVITIES
    printf(" EGADS Error: PCurve semsitivties are not yet supported (EG_copyGeometry_dot)!\n");
    return EGADS_CONSTERR;
#endif

    /* Need code from BRepTools_TrsfModification::NewCurve2d */

    if (geom1->mtype == geom2->mtype) {
      if (outLevel > 0) {
        printf(" EGADS Error: Inconsistent PCurve types %d != %d (EG_copyGeometry_dot)!\n",
               geom1->mtype, geom2->mtype);
      }
      return EGADS_GEOMERR;
    }

    egadsPCurve *ppcurv1 = (egadsPCurve *) geom1->blind;
    if (ppcurv1->data_dot == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: P-Curve without data_dot (EG_copyGeometry_dot)!\n");
      return EGADS_NODATA;
    }
    len       = ppcurv1->dataLen;
    rdata_dot = ppcurv1->data_dot;

    egadsPCurve *ppcurv2 = (egadsPCurve *) geom2->blind;
    if (ppcurv2->data_dot == NULL) {
      ppcurv2->data_dot = (SurrealS<1> *) EG_alloc(len*sizeof(SurrealS<1>));
      if (ppcurv2->data_dot == NULL) return EGADS_MALLOC;
    }
    cdata     = ppcurv2->data;
    cdata_dot = ppcurv2->data_dot;

    for (i = 0; i < len; i++) {
      cdata_dot[i] = rdata_dot[i];
    }

  } else if (geom1->oclass == CURVE) {

    if (geom1->mtype != geom2->mtype) {
      if (outLevel > 0) {
        printf(" EGADS Error: Inconsistent Curve types %d != %d (EG_copyGeometry_dot)!\n",
               geom1->mtype, geom2->mtype);
      }
      return EGADS_GEOMERR;
    }

    egadsCurve *pcurve1 = (egadsCurve *) geom1->blind;
    if (pcurve1->data_dot == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Curve without data_dot (EG_copyGeometry_dot)!\n");
      return EGADS_NODATA;
    }
    len       = pcurve1->dataLen;
    rdata_dot = pcurve1->data_dot;

    egadsCurve *pcurve2 = (egadsCurve *) geom2->blind;
    if (pcurve2->data_dot == NULL) {
      pcurve2->data_dot = (SurrealS<1> *) EG_alloc(len*sizeof(SurrealS<1>));
      if (pcurve2->data_dot == NULL) return EGADS_MALLOC;
    }
    cdata     = pcurve2->data;
    cdata_dot = pcurve2->data_dot;

    for (i = 0; i < len; i++) cdata_dot[i] = rdata_dot[i];

    if (xform != NULL) {
      switch (geom1->mtype) {

        case LINE:
          {
            // Geom_Line::Transform
            egPnt<T> pos(cdata_dot[0], cdata_dot[1], cdata_dot[2]);
            egDir<T> dir(cdata_dot[3], cdata_dot[4], cdata_dot[5]);
            pos.Transform(form);
            dir.Transform(form);
          }
          break;

        case CIRCLE:
          {
            // Geom_Circle::Transform
            egPnt<T> pos (cdata_dot[0], cdata_dot[1], cdata_dot[2]);
            egDir<T> dir1(cdata_dot[3], cdata_dot[4], cdata_dot[5]);
            egDir<T> dir2(cdata_dot[6], cdata_dot[7], cdata_dot[8]);
            pos.Transform(form);
            dir1.Transform(form);
            dir2.Transform(form);
            cdata_dot[ 9] *= fabs(form.ScaleFactor()); // radius
           }
          break;

        case ELLIPSE:
        case HYPERBOLA:
          {
            // Geom_Ellipse::Transform
            // Geom_Hyperbola::Transform
            egPnt<T> pos (cdata_dot[0], cdata_dot[1], cdata_dot[2]);
            egDir<T> dir1(cdata_dot[3], cdata_dot[4], cdata_dot[5]);
            egDir<T> dir2(cdata_dot[6], cdata_dot[7], cdata_dot[8]);
            pos.Transform(form);
            dir1.Transform(form);
            dir2.Transform(form);
            cdata_dot[ 9] *= fabs(form.ScaleFactor()); // minor radius
            cdata_dot[10] *= fabs(form.ScaleFactor()); // major radius
          }
          break;

        case PARABOLA:
          {
            // Geom_Parabola::Transform
            egPnt<T> pos (cdata_dot[0], cdata_dot[1], cdata_dot[2]);
            egDir<T> dir1(cdata_dot[3], cdata_dot[4], cdata_dot[5]);
            egDir<T> dir2(cdata_dot[6], cdata_dot[7], cdata_dot[8]);
            pos.Transform(form);
            dir1.Transform(form);
            dir2.Transform(form);
            cdata_dot[ 9] *= fabs(form.ScaleFactor()); // focal
          }
          break;

        case TRIMMED:
          {
            // Geom_TrimmedCurve::Transform

            // Geom_Line::TransformedParameter
            Handle(Geom_Line) hLine = Handle(Geom_Line)::DownCast(pcurve1->handle);
            if (!hLine.IsNull()) {
              cdata_dot[0] *= fabs(form.ScaleFactor());
              cdata_dot[1] *= fabs(form.ScaleFactor());
            }

            // Geom_Parabola::TransformedParameter
            Handle(Geom_Parabola) hParab = Handle(Geom_Parabola)::DownCast(pcurve1->handle);
            if (!hParab.IsNull()) {
              cdata_dot[0] *= fabs(form.ScaleFactor());
              cdata_dot[1] *= fabs(form.ScaleFactor());
            }

            // Nothing to do for any other curve
          }
          break;

        case BEZIER:
          {
            // Geom_BezierCurve::Transform
            ints = pcurve2->header;
            len = 0;
            for (i = 1; i <= ints[2]; i++, len+=3) {
              egPnt<T> pos(cdata_dot[len  ], cdata_dot[len+1], cdata_dot[len+2]);
              pos.Transform(form);
            }
          }
          break;

        case BSPLINE:
          {
            // Geom_BSplineCurve::Transform
            ints = pcurve2->header;
            len = ints[3];
            for (i = 0; i < ints[2]; i++, len+=3) {
              egPnt<T> pos(cdata_dot[len  ], cdata_dot[len+1], cdata_dot[len+2]);
              pos.Transform(form);
            }
          }
          break;

        case OFFSET:
          {
            // Geom_OffsetCurve::Transform
            egDir<T> dir(cdata_dot[0], cdata_dot[1], cdata_dot[2]);
            dir.Transform(form);
            cdata_dot[3] *= form.ScaleFactor(); // offset
          }
          break;
      }
    } // xform != NULL

    if (pcurve1->ref != NULL) {
      if (pcurve2->ref == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Inconsistent CURVE ref geometry! (EG_copyGeometry_dot)\n");
        return EGADS_GEOMERR;
      }

      /* transfer the data_dot information on the reference geometry */
      stat = EG_copyGeometry_dot(pcurve1->ref, xform, form, pcurve2->ref);
      if (stat != EGADS_SUCCESS) return stat;
    }

  } else {

    if (geom1->mtype != geom2->mtype) {
      if (outLevel > 0) {
        printf(" EGADS Error: Inconsistent Surface types %d != %d (EG_copyGeometry_dot)!\n",
               geom1->mtype, geom2->mtype);
      }
      return EGADS_GEOMERR;
    }

    egadsSurface *psurf1 = (egadsSurface *) geom1->blind;
    if (psurf1->data_dot == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface without data_dot (EG_copyGeometry_dot)!\n");
      return EGADS_NODATA;
    }
    len       = psurf1->dataLen;
    rdata_dot = psurf1->data_dot;

    egadsSurface *psurf2 = (egadsSurface *) geom2->blind;
    if (psurf2->data_dot == NULL) {
      psurf2->data_dot = (SurrealS<1> *) EG_alloc(len*sizeof(SurrealS<1>));
      if (psurf2->data_dot == NULL) return EGADS_MALLOC;
    }
    cdata     = psurf2->data;
    cdata_dot = psurf2->data_dot;

    for (i = 0; i < len; i++) {
      cdata_dot[i] = rdata_dot[i];
    }

    if (xform != NULL) {
      switch (geom1->mtype) {

        case PLANE:
          {
            // Geom_Plane::Transform
            egPnt<T> pos (cdata_dot[0], cdata_dot[1], cdata_dot[2]);
            egDir<T> dir1(cdata_dot[3], cdata_dot[4], cdata_dot[5]);
            egDir<T> dir2(cdata_dot[6], cdata_dot[7], cdata_dot[8]);
            pos.Transform(form);
            dir1.Transform(form);
            dir2.Transform(form);
          }
          break;

        case SPHERICAL:
          {
            // Geom_SphericalSurface::Transform
            egPnt<T> pos (cdata_dot[0], cdata_dot[1], cdata_dot[2]);
            egDir<T> dir1(cdata_dot[3], cdata_dot[4], cdata_dot[5]);
            egDir<T> dir2(cdata_dot[6], cdata_dot[7], cdata_dot[8]);
            pos.Transform(form);
            dir1.Transform(form);
            dir2.Transform(form);
            cdata_dot[ 9] *= fabs(form.ScaleFactor()); // radius
          }
          break;

        case CONICAL:
          {
            // Geom_ConicalSurface::Transform
            egPnt<T> pos (cdata_dot[ 0], cdata_dot[ 1], cdata_dot[ 2]);
            egDir<T> dir1(cdata_dot[ 3], cdata_dot[ 4], cdata_dot[ 5]);
            egDir<T> dir2(cdata_dot[ 6], cdata_dot[ 7], cdata_dot[ 8]);
            egDir<T> dir3(cdata_dot[ 9], cdata_dot[10], cdata_dot[11]);
            pos.Transform(form);
            dir1.Transform(form);
            dir2.Transform(form);
            dir3.Transform(form);
            cdata_dot[13] *= fabs(form.ScaleFactor()); // radius
          }
          break;

        case CYLINDRICAL:
          {
            // Geom_CylindricalSurface::Transform
            egPnt<T> pos (cdata_dot[ 0], cdata_dot[ 1], cdata_dot[ 2]);
            egDir<T> dir1(cdata_dot[ 3], cdata_dot[ 4], cdata_dot[ 5]);
            egDir<T> dir2(cdata_dot[ 6], cdata_dot[ 7], cdata_dot[ 8]);
            egDir<T> dir3(cdata_dot[ 9], cdata_dot[10], cdata_dot[11]);
            pos.Transform(form);
            dir1.Transform(form);
            dir2.Transform(form);
            dir3.Transform(form);
            cdata_dot[12] *= fabs(form.ScaleFactor()); // radius
          }
          break;

        case TOROIDAL:
          {
            // Geom_ToroidalSurface::Transform
            egPnt<T> pos (cdata_dot[ 0], cdata_dot[ 1], cdata_dot[ 2]);
            egDir<T> dir1(cdata_dot[ 3], cdata_dot[ 4], cdata_dot[ 5]);
            egDir<T> dir2(cdata_dot[ 6], cdata_dot[ 7], cdata_dot[ 8]);
            egDir<T> dir3(cdata_dot[ 9], cdata_dot[10], cdata_dot[11]);
            pos.Transform(form);
            dir1.Transform(form);
            dir2.Transform(form);
            dir3.Transform(form);
            cdata_dot[12] *= fabs(form.ScaleFactor()); // minor radius
            cdata_dot[13] *= fabs(form.ScaleFactor()); // major radius
          }
          break;

        case BEZIER:
          {
            // Geom_BezierSurface::Transform
            ints = psurf2->header;
            len = 0;
            for (i = 0; i < ints[2]*ints[4]; i++, len+=3) {
              egPnt<T> pos(cdata_dot[len  ], cdata_dot[len+1], cdata_dot[len+2]);
              pos.Transform(form);
            }
          }
          break;

        case BSPLINE:
          {
            // Geom_BSplineSurface::Transform
            ints = psurf2->header;
            int nCP = ints[2];
            nCP    *= ints[5];
            len = ints[3] + ints[6];
            for (i = 0; i < nCP; i++, len+=3) {
              egPnt<T> pos(cdata_dot[len  ], cdata_dot[len+1], cdata_dot[len+2]);
              pos.Transform(form);
            }
          }
          break;

        case OFFSET:
          {
            // Geom_OffsetSurface::Transform
            cdata_dot[0] *= fabs(form.ScaleFactor());
          }
          break;

        case TRIMMED:
          {
            // Geom_RectangularTrimmedSurface::Transform

            // Geom_Plane::TransformedParameter
            Handle(Geom_Plane) hPlane = Handle(Geom_Plane)::DownCast(psurf1->handle);
            if (!hPlane.IsNull()) {
              cdata_dot[0] *= fabs(form.ScaleFactor());
              cdata_dot[1] *= fabs(form.ScaleFactor());
              cdata_dot[2] *= fabs(form.ScaleFactor());
              cdata_dot[3] *= fabs(form.ScaleFactor());
            }

            // Geom_ConicalSurface::TransformedParameter
            Handle(Geom_ConicalSurface) hConical =
              Handle(Geom_ConicalSurface)::DownCast(psurf1->handle);
            if (!hConical.IsNull()) {
              cdata_dot[2] *= fabs(form.ScaleFactor());
              cdata_dot[3] *= fabs(form.ScaleFactor());
            }

            // Geom_ConicalSurface::TransformedParameter
            Handle(Geom_CylindricalSurface) hCyl = Handle(Geom_CylindricalSurface)::DownCast(psurf1->handle);
            if (!hConical.IsNull()) {
              cdata_dot[2] *= fabs(form.ScaleFactor());
              cdata_dot[3] *= fabs(form.ScaleFactor());
            }
          }
          break;

        case EXTRUSION:
          {
            // Geom_SurfaceOfLinearExtrusion::Transform
            egDir<T> dir(cdata_dot[0], cdata_dot[1], cdata_dot[2]);
            dir.Transform(form);
          }
          break;

        case REVOLUTION:
          {
            // Geom_SurfaceOfRevolution::Transform
            egPnt<T> pos(cdata_dot[0], cdata_dot[1], cdata_dot[2]);
            egDir<T> dir(cdata_dot[3], cdata_dot[4], cdata_dot[5]);
            pos.Transform(form);
            dir.Transform(form);
          }
          break;

      }
    } // xform != NULL

    if (psurf1->ref != NULL) {
      if (psurf2->ref == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Inconsistent SURFACE ref geometry! (EG_copyGeometry_dot)\n");
        return EGADS_GEOMERR;
      }

      /* transfer the data_dot information on the reference geometry */
      stat = EG_copyGeometry_dot(psurf1->ref, xform, form, psurf2->ref);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }

  /* check consistency in the data */
  stat = EGADS_SUCCESS;
  scale = 0.0;
  for (i = 0; i < len; i++) {
    if (fabs(cdata[i]) > scale) scale = fabs(cdata[i]);
  }
  if (scale == 0.0) scale = 1.0;
  for (i = 0; i < len; i++)
    if (fabs(cdata[i] - cdata_dot[i].value()) > 1.e-14*scale) {
      stat++;
      break;
    }
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0) {
      printf(" EGADS Error: Inconsistent geometry data! (EG_copyGeometry_dot)\n");
      if (outLevel > 1) {
        for (i = 0; i < len; i++)
          printf("     data[%d] %lf : %lf\n", i, cdata[i], cdata_dot[i].value());
      }
    }
    return EGADS_GEOMERR;
  }

  return EGADS_SUCCESS;
}


int
EG_copyGeometry_dot(const egObject *obj,
                    /*@null@*/ const double *xform,
                    /*@null@*/ const double *xform_dot,
                    egObject *copy)
{
  typedef SurrealS<1> T;
  int outLevel;

  if  (obj == NULL)                 return EGADS_NULLOBJ;
  if  (obj->magicnumber != MAGIC)   return EGADS_NOTOBJ;

  outLevel = EG_outLevel(obj);

  if (xform == NULL && xform_dot != NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: xform_dot also requires xform (EG_copyGeometry_dot)!\n");
    return EGADS_NODATA;
  }

  /* this avoids repeated calculations of the transformation when obj is a BODY */
  egTrsf<T> form;

  if (xform != NULL) {
    if (xform_dot == NULL) {
      form.SetValues(xform[ 0], xform[ 1], xform[ 2], xform[ 3],
                     xform[ 4], xform[ 5], xform[ 6], xform[ 7],
                     xform[ 8], xform[ 9], xform[10], xform[11]);

    } else {
      form.SetValues(T(xform[ 0], xform_dot[ 0]), T(xform[ 1], xform_dot[ 1]),
                     T(xform[ 2], xform_dot[ 2]), T(xform[ 3], xform_dot[ 3]),
                     T(xform[ 4], xform_dot[ 4]), T(xform[ 5], xform_dot[ 5]),
                     T(xform[ 6], xform_dot[ 6]), T(xform[ 7], xform_dot[ 7]),
                     T(xform[ 8], xform_dot[ 8]), T(xform[ 9], xform_dot[ 9]),
                     T(xform[10], xform_dot[10]), T(xform[11], xform_dot[11]));
    }
  }

  return EG_copyGeometry_dot(obj, xform, form, copy);
}


int
EG_copyGeometry_dot(const egObject *obj,
                    /*@null@*/ const SurrealS<1> *xform,
                    egObject *copy)
{
  typedef SurrealS<1> T;
  double *xform_flag=NULL, xform_val;

  if  (obj == NULL)                 return EGADS_NULLOBJ;
  if  (obj->magicnumber != MAGIC)   return EGADS_NOTOBJ;

  /* this avoids repeated calculations of the transformation when obj is a BODY */
  egTrsf<T> form;

  if (xform != NULL) {
    form.SetValues(xform[ 0], xform[ 1], xform[ 2], xform[ 3],
                   xform[ 4], xform[ 5], xform[ 6], xform[ 7],
                   xform[ 8], xform[ 9], xform[10], xform[11]);

    /* this is just a flag, so the value does not matter */
    xform_flag = &xform_val;
  }

  return EG_copyGeometry_dot(obj, xform_flag, form, copy);
}


void
EG_completePCurve(egObject *geom, Handle(Geom2d_Curve) &hCurve)
{
  int                  m, n, stat, oclass, mtype;
  double               d, x0[2], x1[2];
  egObject             *obj, *ref;
  Handle(Geom2d_Curve) base;

  geom->oclass        = PCURVE;
  egadsPCurve *ppcurv = new egadsPCurve;
  ppcurv->handle      = hCurve;
  ppcurv->ref         = NULL;
  ppcurv->topFlg      = 0;
  ppcurv->header      = NULL;
  ppcurv->data        = NULL;
  ppcurv->data_dot    = NULL;
  geom->blind         = ppcurv;

  // stand alone geometry
  Handle(Geom2d_Line) hLine = Handle(Geom2d_Line)::DownCast(hCurve);
  if (!hLine.IsNull()) {
    geom->mtype = LINE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &ppcurv->header, &ppcurv->data);
    EG_getGeometryLen(geom, &m, &ppcurv->dataLen);
    ppcurv->trange[0] = hCurve->FirstParameter();
    ppcurv->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom2d_Circle) hCirc = Handle(Geom2d_Circle)::DownCast(hCurve);
  if (!hCirc.IsNull()) {
    geom->mtype = CIRCLE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &ppcurv->header, &ppcurv->data);
    EG_getGeometryLen(geom, &m, &ppcurv->dataLen);
    ppcurv->trange[0] = hCurve->FirstParameter();
    ppcurv->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom2d_Ellipse) hEllip = Handle(Geom2d_Ellipse)::DownCast(hCurve);
  if (!hEllip.IsNull()) {
    geom->mtype = ELLIPSE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &ppcurv->header, &ppcurv->data);
    EG_getGeometryLen(geom, &m, &ppcurv->dataLen);
    ppcurv->trange[0] = hCurve->FirstParameter();
    ppcurv->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom2d_Parabola) hParab = Handle(Geom2d_Parabola)::DownCast(hCurve);
  if (!hParab.IsNull()) {
    geom->mtype = PARABOLA;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &ppcurv->header, &ppcurv->data);
    EG_getGeometryLen(geom, &m, &ppcurv->dataLen);
    ppcurv->trange[0] = hCurve->FirstParameter();
    ppcurv->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom2d_Hyperbola) hHypr = Handle(Geom2d_Hyperbola)::DownCast(hCurve);
  if (!hHypr.IsNull()) {
    geom->mtype = HYPERBOLA;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &ppcurv->header, &ppcurv->data);
    EG_getGeometryLen(geom, &m, &ppcurv->dataLen);
    ppcurv->trange[0] = hCurve->FirstParameter();
    ppcurv->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom2d_BezierCurve) hBezier = Handle(Geom2d_BezierCurve)::DownCast(hCurve);
  if (!hBezier.IsNull()) {
    geom->mtype = BEZIER;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &ppcurv->header, &ppcurv->data);
    EG_getGeometryLen(geom, &m, &ppcurv->dataLen);
    ppcurv->trange[0] = hCurve->FirstParameter();
    ppcurv->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom2d_BSplineCurve) hBSpline = Handle(Geom2d_BSplineCurve)::DownCast(hCurve);
  if (!hBSpline.IsNull()) {
    geom->mtype = BSPLINE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &ppcurv->header, &ppcurv->data);
    EG_getGeometryLen(geom, &m, &ppcurv->dataLen);
    ppcurv->trange[0] = hCurve->FirstParameter();
    ppcurv->trange[1] = hCurve->LastParameter();
    for (n = 2; n < ppcurv->header[2]; n++) {
      m     = ppcurv->header[3] + 2*n - 4;
      x0[0] = ppcurv->data[m+2] - ppcurv->data[m  ];
      x0[1] = ppcurv->data[m+3] - ppcurv->data[m+1];
      d     = sqrt(x0[0]*x0[0] + x0[1]*x0[1]);
      if (d != 0.0) {
        x0[0] /= d;
        x0[1] /= d;
      }
      x1[0] = ppcurv->data[m+4] - ppcurv->data[m+2];
      x1[1] = ppcurv->data[m+5] - ppcurv->data[m+3];
      d     = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
      if (d != 0.0) {
        x1[0] /= d;
        x1[1] /= d;
      }
      d = x0[0]*x1[0] + x0[1]*x1[1];
      if (d < -0.95) {
#ifdef DEBUG
        printf(" EGADS Info: PCurve dot flip at %d/%d (%lf) -- %lf %lf!\n",
               n-2, ppcurv->header[2]-2, d, ppcurv->data[0],
               ppcurv->data[ppcurv->header[3]-1]);
#endif
        stat = EG_addStrAttr(geom, ".Bad", "CPrev");
        if (stat != EGADS_SUCCESS)
          printf("             EG_addStrAttr CPrev= %d\n", stat);
      }
    }
    return;
  }

  // referencing geometry
  Handle(Geom2d_TrimmedCurve) hTrim = Handle(Geom2d_TrimmedCurve)::DownCast(hCurve);
  if (!hTrim.IsNull()) {
    geom->mtype = TRIMMED;
    base = hTrim->BasisCurve();
  }
  Handle(Geom2d_OffsetCurve) hOffst = Handle(Geom2d_OffsetCurve)::DownCast(hCurve);
  if (!hOffst.IsNull()) {
    geom->mtype = OFFSET;
    base = hOffst->BasisCurve();
  }
  if (geom->mtype == 0) {
    printf(" EGADS Error: Unknown PCurve Type!\n");
    return;
  }

  // make the reference curve
  stat = EG_makeObject(EG_context(geom), &obj);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: make Object = %d (EG_completePCurve)!\n", stat);
    return;
  }
  ppcurv->ref = obj;
  if (geom->topObj == EG_context(geom)) {
    obj->topObj = geom;
  } else {
    obj->topObj = geom->topObj;
  }
  EG_completePCurve(obj,  base);
  EG_referenceObject(obj, geom);
  EG_getGeometry(geom, &oclass, &mtype, &ref, &ppcurv->header, &ppcurv->data);
  EG_getGeometryLen(geom, &m, &ppcurv->dataLen);
  ppcurv->trange[0] = hCurve->FirstParameter();
  ppcurv->trange[1] = hCurve->LastParameter();
}


void
EG_completeCurve(egObject *geom, Handle(Geom_Curve) &hCurve)
{
  int                stat, m, oclass, mtype;
  egObject           *obj, *ref;
  Handle(Geom_Curve) base;

  geom->oclass       = CURVE;
  egadsCurve *pcurve = new egadsCurve;
  pcurve->handle     = hCurve;
  pcurve->ref        = NULL;
  pcurve->topFlg     = 0;
  pcurve->header     = NULL;
  pcurve->data       = NULL;
  pcurve->data_dot   = NULL;
  geom->blind        = pcurve;

  // stand alone geometry
  Handle(Geom_Line) hLine = Handle(Geom_Line)::DownCast(hCurve);
  if (!hLine.IsNull()) {
    geom->mtype = LINE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &pcurve->header, &pcurve->data);
    EG_getGeometryLen(geom, &m, &pcurve->dataLen);
    pcurve->trange[0] = hCurve->FirstParameter();
    pcurve->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom_Circle) hCirc = Handle(Geom_Circle)::DownCast(hCurve);
  if (!hCirc.IsNull()) {
    geom->mtype = CIRCLE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &pcurve->header, &pcurve->data);
    EG_getGeometryLen(geom, &m, &pcurve->dataLen);
    pcurve->trange[0] = hCurve->FirstParameter();
    pcurve->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom_Ellipse) hEllip = Handle(Geom_Ellipse)::DownCast(hCurve);
  if (!hEllip.IsNull()) {
    geom->mtype = ELLIPSE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &pcurve->header, &pcurve->data);
    EG_getGeometryLen(geom, &m, &pcurve->dataLen);
    pcurve->trange[0] = hCurve->FirstParameter();
    pcurve->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom_Parabola) hParab = Handle(Geom_Parabola)::DownCast(hCurve);
  if (!hParab.IsNull()) {
    geom->mtype = PARABOLA;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &pcurve->header, &pcurve->data);
    EG_getGeometryLen(geom, &m, &pcurve->dataLen);
    pcurve->trange[0] = hCurve->FirstParameter();
    pcurve->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom_Hyperbola) hHypr = Handle(Geom_Hyperbola)::DownCast(hCurve);
  if (!hHypr.IsNull()) {
    geom->mtype = HYPERBOLA;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &pcurve->header, &pcurve->data);
    EG_getGeometryLen(geom, &m, &pcurve->dataLen);
    pcurve->trange[0] = hCurve->FirstParameter();
    pcurve->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom_BezierCurve) hBezier = Handle(Geom_BezierCurve)::DownCast(hCurve);
  if (!hBezier.IsNull()) {
    geom->mtype = BEZIER;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &pcurve->header, &pcurve->data);
    EG_getGeometryLen(geom, &m, &pcurve->dataLen);
    pcurve->trange[0] = hCurve->FirstParameter();
    pcurve->trange[1] = hCurve->LastParameter();
    return;
  }
  Handle(Geom_BSplineCurve) hBSpline = Handle(Geom_BSplineCurve)::DownCast(hCurve);
  if (!hBSpline.IsNull()) {
    geom->mtype = BSPLINE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &pcurve->header, &pcurve->data);
    EG_getGeometryLen(geom, &m, &pcurve->dataLen);
    pcurve->trange[0] = hCurve->FirstParameter();
    pcurve->trange[1] = hCurve->LastParameter();
    return;
  }

  // referencing geometry
  Handle(Geom_TrimmedCurve) hTrim = Handle(Geom_TrimmedCurve)::DownCast(hCurve);
  if (!hTrim.IsNull()) {
    geom->mtype = TRIMMED;
    base = hTrim->BasisCurve();
  }
  Handle(Geom_OffsetCurve) hOffst = Handle(Geom_OffsetCurve)::DownCast(hCurve);
  if (!hOffst.IsNull()) {
    geom->mtype = OFFSET;
    base = hOffst->BasisCurve();
  }
  if (geom->mtype == 0) {
    printf(" EGADS Error: Unknown Curve Type!\n");
    return;
  }

  // make the reference curve
  stat = EG_makeObject(EG_context(geom), &obj);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: make Object = %d (EG_completeCurve)!\n", stat);
    return;
  }
  pcurve->ref = obj;
  if (geom->topObj == EG_context(geom)) {
    obj->topObj = geom;
  } else {
    obj->topObj = geom->topObj;
  }
  EG_completeCurve(obj, base);
  EG_referenceObject(obj, geom);
  EG_getGeometry(geom, &oclass, &mtype, &ref, &pcurve->header, &pcurve->data);
  EG_getGeometryLen(geom, &m, &pcurve->dataLen);
  pcurve->trange[0] = hCurve->FirstParameter();
  pcurve->trange[1] = hCurve->LastParameter();
}


void EG_completeSurf(egObject *geom, Handle(Geom_Surface) &hSurf)
{
  int                  stat, m, oclass, mtype;
  egObject             *obj, *ref;
  Handle(Geom_Surface) base;
  Handle(Geom_Curve)   curve;

  geom->oclass        = SURFACE;
  egadsSurface *psurf = new egadsSurface;
  psurf->handle       = hSurf;
  psurf->ref          = NULL;
  psurf->topFlg       = 0;
  psurf->header       = NULL;
  psurf->data         = NULL;
  psurf->data_dot     = NULL;
  geom->blind         = psurf;

  // stand alone geometry
  Handle(Geom_Plane) hPlane = Handle(Geom_Plane)::DownCast(hSurf);
  if (!hPlane.IsNull()) {
    geom->mtype = PLANE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(geom, &m, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    return;
  }
  Handle(Geom_SphericalSurface) hSphere =
      Handle(Geom_SphericalSurface)::DownCast(hSurf);
  if (!hSphere.IsNull()) {
    geom->mtype = SPHERICAL;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(geom, &m, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    return;
  }
  Handle(Geom_ConicalSurface) hCone =
    Handle(Geom_ConicalSurface)::DownCast(hSurf);
  if (!hCone.IsNull()) {
    geom->mtype = CONICAL;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(geom, &m, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    return;
  }
  Handle(Geom_CylindricalSurface) hCyl =
    Handle(Geom_CylindricalSurface)::DownCast(hSurf);
  if (!hCyl.IsNull()) {
    geom->mtype = CYLINDRICAL;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(geom, &m, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    return;
  }
  Handle(Geom_ToroidalSurface) hTorus =
    Handle(Geom_ToroidalSurface)::DownCast(hSurf);
  if (!hTorus.IsNull()) {
    geom->mtype = TOROIDAL;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(geom, &m, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    return;
  }
  Handle(Geom_BezierSurface) hBezier =
    Handle(Geom_BezierSurface)::DownCast(hSurf);
  if (!hBezier.IsNull()) {
    geom->mtype = BEZIER;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(geom, &m, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    return;
  }
  Handle(Geom_BSplineSurface) hBSpline =
    Handle(Geom_BSplineSurface)::DownCast(hSurf);
  if (!hBSpline.IsNull()) {
    geom->mtype = BSPLINE;
    EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(geom, &m, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    return;
  }

  // referencing geometry -- surface
  Handle(Geom_OffsetSurface) hOffst =
    Handle(Geom_OffsetSurface)::DownCast(hSurf);
  if (!hOffst.IsNull()) {
    geom->mtype = OFFSET;
    base = hOffst->BasisSurface();
    stat = EG_makeObject(EG_context(geom), &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Object = %d (EG_completeSurface)!\n", stat);
      return;
    }
    psurf->ref = obj;
    if (geom->topObj == EG_context(geom)) {
      obj->topObj = geom;
    } else {
      obj->topObj = geom->topObj;
    }
    EG_completeSurf(obj, base);
    EG_referenceObject(obj, geom);
    EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(geom, &m, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    return;
  }
  Handle(Geom_RectangularTrimmedSurface) hTrim =
    Handle(Geom_RectangularTrimmedSurface)::DownCast(hSurf);
  if (!hTrim.IsNull()) {
    geom->mtype = TRIMMED;
    base = hTrim->BasisSurface();
    stat = EG_makeObject(EG_context(geom), &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Object = %d (EG_completeSurface)!\n", stat);
      return;
    }
    psurf->ref = obj;
    if (geom->topObj == EG_context(geom)) {
      obj->topObj = geom;
    } else {
      obj->topObj = geom->topObj;
    }
    EG_completeSurf(obj, base);
    EG_referenceObject(obj, geom);
    EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(geom, &m, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    return;
  }

  // referencing geometry -- curve
  Handle(Geom_SurfaceOfLinearExtrusion) hSLExtr =
    Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(hSurf);
  if (!hSLExtr.IsNull()) {
    geom->mtype = EXTRUSION;
    curve       = hSLExtr->BasisCurve();
  }
  Handle(Geom_SurfaceOfRevolution) hSORev =
    Handle(Geom_SurfaceOfRevolution)::DownCast(hSurf);
  if (!hSORev.IsNull()) {
    geom->mtype = REVOLUTION;
    curve       = hSORev->BasisCurve();
  }
  if (geom->mtype == 0) {
    printf(" EGADS Error: Unknown Surface Type!\n");
    return;
  }

  // make the reference curve
  stat = EG_makeObject(EG_context(geom), &obj);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: make Curve = %d (EG_completeSurface)!\n", stat);
    return;
  }
  psurf->ref = obj;
  if (geom->topObj == EG_context(geom)) {
    obj->topObj = geom;
  } else {
    obj->topObj = geom->topObj;
  }
  EG_completeCurve(obj, curve);
  EG_referenceObject(obj, geom);
  EG_getGeometry(geom, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
  EG_getGeometryLen(geom, &m, &psurf->dataLen);
  hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                psurf->vrange[0], psurf->vrange[1]);
}


int
EG_copyGeometry(/*@null@*/ egObject *context, const egObject *geom,
                /*@null@*/ double *xform, egObject **copy)
{
  int      stat, outLevel;
  egObject *obj;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (context == NULL)            return EGADS_NOTCNTX;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != CURVE) && (geom->oclass != SURFACE))
                                   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(geom);

  gp_Trsf form = gp_Trsf();
  if (xform != NULL)
    form.SetValues(xform[ 0], xform[ 1], xform[ 2], xform[ 3],
                   xform[ 4], xform[ 5], xform[ 6], xform[ 7],
                   xform[ 8], xform[ 9], xform[10], xform[11]);

  if (geom->oclass == CURVE) {

    egadsCurve           *pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve)    hCurve = pcurve->handle;
    Handle(Geom_Geometry) nGeom  = hCurve->Transformed(form);
    Handle(Geom_Curve)    nCurve = Handle(Geom_Curve)::DownCast(nGeom);
    if (nCurve.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Error: XForm Curve Failed (EG_copyGeometry)!\n");
      return EGADS_CONSTERR;
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeObject = %d (EG_copyGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeCurve(obj, nCurve);

  } else {

    egadsSurface         *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface)  hSurf = psurf->handle;
    Handle(Geom_Geometry) nGeom = hSurf->Transformed(form);
    Handle(Geom_Surface)  nSurf = Handle(Geom_Surface)::DownCast(nGeom);
    if (nSurf.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Error: XForm Surface Failed (EG_copyGeometry)!\n");
      return EGADS_CONSTERR;
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeObject = %d (EG_copyGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeSurf(obj, nSurf);

  }

  EG_referenceObject(obj, context);
  *copy = obj;
  return EGADS_SUCCESS;
}


int
EG_flipGeometry(const egObject *geom, egObject **copy)
{
  int      stat, outLevel;
  egObject *obj, *context;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE)     && (geom->oclass != CURVE) &&
      (geom->oclass != SURFACE))
                                   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(geom))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(geom);
  context  = EG_context(geom);

  if (geom->oclass == PCURVE) {

    egadsPCurve          *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve)  hPCurv = ppcurv->handle;
    Handle(Geom2d_Curve)  nPCurv = hPCurv->Reversed();

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeObject = %d (EG_flipGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completePCurve(obj, nPCurv);

  } else if (geom->oclass == CURVE) {

    egadsCurve         *pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve)  hCurve = pcurve->handle;
    Handle(Geom_Curve)  nCurve = hCurve->Reversed();

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeObject = %d (EG_flipGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeCurve(obj, nCurve);

  } else {

    egadsSurface         *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface)  hSurf = psurf->handle;
    Handle(Geom_Surface)  nSurf = hSurf->UReversed();

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeObject = %d (EG_flipGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeSurf(obj, nSurf);

  }

  EG_referenceObject(obj, context);
  *copy = obj;
  return EGADS_SUCCESS;
}


int
EG_makeGeometry(egObject *context, int oclass, int mtype,
                /*@null@*/ egObject *refGeom, /*@null@*/ const int *ints,
                const double *data, egObject **geom)
{
  int      i, j, m, n, len, stat, outLevel, rational, nmult;
  double   d, x0[2], x1[2];
  egObject *ref, *obj = NULL, *basis = NULL;

  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (EG_sameThread(context))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(context);

  if ((oclass < PCURVE) || (oclass > SURFACE)) {
    if (outLevel > 0)
      printf(" EGADS Error: oclass = %d (EG_makeGeometry)!\n", oclass);
    return EGADS_NOTGEOM;
  }

  if (oclass == PCURVE) {

    if ((mtype < LINE) || (mtype > OFFSET)) {
      if (outLevel > 0)
        printf(" EGADS Error: PCurve mtype = %d (EG_makeGeometry)!\n",
               mtype);
      return EGADS_RANGERR;
    }
    if ((mtype == TRIMMED) || (mtype == OFFSET)) {
      if (refGeom == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: PCrv mtype = %d Ref is NULL (EG_makeGeometry)!\n",
                 mtype);
        return EGADS_NULLOBJ;
      }
      if (refGeom->oclass != PCURVE) {
        if (outLevel > 0)
          printf(" EGADS Error: PCrv mtype = %d Ref is %d (EG_makeGeometry)!\n",
                 mtype, refGeom->oclass);
        return EGADS_NOTGEOM;
      }
      basis = refGeom;
      if (basis->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: PCrv mtype = %d Ref has no data (EG_makeGeometry)!\n",
                 mtype);
        return EGADS_NODATA;
      }
    }

    Handle(Geom2d_Curve) hCurve;
    try {
      switch (mtype) {

        case LINE:
          {
            gp_Pnt2d pntl(data[0], data[1]);
            gp_Dir2d dirl(data[2], data[3]);
            hCurve = new Geom2d_Line(pntl, dirl);
          }
          break;

        case CIRCLE:
          {
            gp_Pnt2d pntc(data[0], data[1]);
            gp_Dir2d dirx(data[2], data[3]);
            gp_Dir2d diry(data[4], data[5]);
            gp_Ax22d axi2(pntc, dirx, diry);
            hCurve = new Geom2d_Circle(axi2, data[6]);
          }
          break;

        case ELLIPSE:
          {
            gp_Pnt2d pnte(data[0], data[1]);
            gp_Dir2d dirx(data[2], data[3]);
            gp_Dir2d diry(data[4], data[5]);
            gp_Ax22d axi2(pnte, dirx, diry);
            hCurve = new Geom2d_Ellipse(axi2, data[6], data[7]);
          }
          break;

        case PARABOLA:
          {
            gp_Pnt2d pntp(data[0], data[1]);
            gp_Dir2d dirx(data[2], data[3]);
            gp_Dir2d diry(data[4], data[5]);
            gp_Ax22d axi2(pntp, dirx, diry);
            hCurve = new Geom2d_Parabola(axi2, data[6]);
          }
          break;

        case HYPERBOLA:
          {
            gp_Pnt2d pnth(data[0], data[1]);
            gp_Dir2d dirx(data[2], data[3]);
            gp_Dir2d diry(data[4], data[5]);
            gp_Ax22d axi2(pnth, dirx, diry);
            hCurve = new Geom2d_Hyperbola(axi2, data[6], data[7]);
          }
          break;

        case TRIMMED:
          {
            egadsPCurve *ppcurv = (egadsPCurve *) basis->blind;
            hCurve = new Geom2d_TrimmedCurve(ppcurv->handle, data[0], data[1]);
          }
          break;

        case BEZIER:
          {
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            len = ints[2];
            TColgp_Array1OfPnt2d aPoles(1, len);
            for (i = 1; i <= ints[2]; i++)
              aPoles(i) = gp_Pnt2d(data[2*i-2], data[2*i-1]);
            if (rational == 0) {
              hCurve = new Geom2d_BezierCurve(aPoles);
            } else {
              TColStd_Array1OfReal aWeights(1, len);
              len = 2*ints[2];
              for (i = 1; i <= ints[2]; i++, len++)
                aWeights(i) = data[len];
              hCurve = new Geom2d_BezierCurve(aPoles, aWeights);
            }
          }
          break;

        case BSPLINE:
          {
            Standard_Boolean periodic = Standard_False;
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            if ((ints[0]&4) != 0) periodic = Standard_True;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) len++;
            TColStd_Array1OfReal    aKnots(1, len);
            TColStd_Array1OfInteger aMults(1, len);
            aKnots(1) = data[0];
            aMults(1) = nmult = 1;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) {
                len++;
                aKnots(len) = data[i];
                aMults(len) = nmult = 1;
              } else {
                nmult++;
                aMults(len) = nmult;
              }
            len = ints[3];
            TColgp_Array1OfPnt2d aPoles(1, ints[2]);
            for (i = 1; i <= ints[2]; i++, len+=2)
              aPoles(i) = gp_Pnt2d(data[len], data[len+1]);
            if (rational == 0) {
              hCurve = new Geom2d_BSplineCurve(aPoles, aKnots, aMults,
                                               ints[1], periodic);
            } else {
              TColStd_Array1OfReal aWeights(1, ints[2]);
              for (i = 1; i <= ints[2]; i++, len++)
                aWeights(i) = data[len];
              hCurve = new Geom2d_BSplineCurve(aPoles, aWeights, aKnots,
                                               aMults, ints[1], periodic);
            }
          }
          break;

        case OFFSET:
          {
            egadsPCurve *ppcurv = (egadsPCurve *) basis->blind;
            hCurve = new Geom2d_OffsetCurve(ppcurv->handle, data[0]);
          }
          break;
      }
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      printf("                %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make PCurve = %d (EG_makeGeometry)!\n", stat);
      return stat;
    }
    obj->oclass         = PCURVE;
    obj->mtype          = mtype;
    egadsPCurve *ppcurv = new egadsPCurve;
    ppcurv->handle      = hCurve;
    ppcurv->ref         = basis;
    ppcurv->topFlg      = 1;
    ppcurv->header      = NULL;
    ppcurv->data        = NULL;
    ppcurv->data_dot    = NULL;
    obj->blind          = ppcurv;
    EG_getGeometry(obj, &i, &j, &ref, &ppcurv->header, &ppcurv->data);
    EG_getGeometryLen(obj, &i, &ppcurv->dataLen);
    ppcurv->trange[0]   = hCurve->FirstParameter();
    ppcurv->trange[1]   = hCurve->LastParameter();
    EG_referenceObject(obj, context);
    if (basis != NULL) EG_referenceTopObj(basis, obj);
    if (mtype == BSPLINE) {
      for (n = 2; n < ppcurv->header[2]; n++) {
        m     = ppcurv->header[3] + 2*n - 4;
        x0[0] = ppcurv->data[m+2] - ppcurv->data[m  ];
        x0[1] = ppcurv->data[m+3] - ppcurv->data[m+1];
        d     = sqrt(x0[0]*x0[0] + x0[1]*x0[1]);
        if (d != 0.0) {
          x0[0] /= d;
          x0[1] /= d;
        }
        x1[0] = ppcurv->data[m+4] - ppcurv->data[m+2];
        x1[1] = ppcurv->data[m+5] - ppcurv->data[m+3];
        d     = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
        if (d != 0.0) {
          x1[0] /= d;
          x1[1] /= d;
        }
        d = x0[0]*x1[0] + x0[1]*x1[1];
        if (d < -0.95) {
#ifdef DEBUG
          printf(" EGADS Info: PCurve dot flip at %d/%d (%lf) -- %lf %lf!\n",
                 n-2, ppcurv->header[2]-2, d, ppcurv->data[0],
                 ppcurv->data[ppcurv->header[3]-1]);
#endif
          stat = EG_addStrAttr(obj, ".Bad", "CPrev");
          if (stat != EGADS_SUCCESS)
            printf("             EG_addStrAttr CPrev= %d\n", stat);
        }
      }
    }

  } else if (oclass == CURVE) {

    if ((mtype < LINE) || (mtype > OFFSET)) {
      if (outLevel > 0)
        printf(" EGADS Error: Curve mtype = %d (EG_makeGeometry)!\n",
               mtype);
      return EGADS_RANGERR;
    }
    if ((mtype == TRIMMED) || (mtype == OFFSET)) {
      if (refGeom == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Crv mtype = %d Ref is NULL (EG_makeGeometry)!\n",
                 mtype);
        return EGADS_NULLOBJ;
      }
      if (refGeom->oclass != CURVE) {
        if (outLevel > 0)
          printf(" EGADS Error: Crv mtype = %d Ref is %d (EG_makeGeometry)!\n",
                 mtype, refGeom->oclass);
        return EGADS_NOTGEOM;
      }
      basis = refGeom;
      if (basis->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Crv mtype = %d Ref has no data (EG_makeGeometry)!\n",
                 mtype);
        return EGADS_NODATA;
      }
    }

    Handle(Geom_Curve) hCurve;
    try {
      switch (mtype) {

        case LINE:
          {
            gp_Pnt pntl(data[0], data[1], data[2]);
            gp_Dir dirl(data[3], data[4], data[5]);
            hCurve = new Geom_Line(pntl, dirl);
          }
          break;

        case CIRCLE:
          {
            gp_Pnt pntc(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2(pntc, dirz, dirx);
            hCurve = new Geom_Circle(axi2, data[9]);
          }
          break;

        case ELLIPSE:
          {
            gp_Pnt pnte(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2(pnte, dirz, dirx);
            hCurve = new Geom_Ellipse(axi2, data[9], data[10]);
          }
          break;

        case PARABOLA:
          {
            gp_Pnt pntp(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2(pntp, dirz, dirx);
            hCurve = new Geom_Parabola(axi2, data[9]);
          }
          break;

        case HYPERBOLA:
          {
            gp_Pnt pnth(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2(pnth, dirz, dirx);
            hCurve = new Geom_Hyperbola(axi2, data[9], data[10]);
          }
          break;

        case TRIMMED:
          {
            egadsCurve *pcurve = (egadsCurve *) basis->blind;
            hCurve = new Geom_TrimmedCurve(pcurve->handle, data[0], data[1]);
          }
          break;

        case BEZIER:
          {
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            len = ints[2];
            TColgp_Array1OfPnt aPoles(1, len);
            for (i = 1; i <= ints[2]; i++)
              aPoles(i) = gp_Pnt(data[3*i-3], data[3*i-2], data[3*i-1]);
            if (rational == 0) {
              hCurve = new Geom_BezierCurve(aPoles);
            } else {
              TColStd_Array1OfReal aWeights(1, len);
              len = 3*ints[2];
              for (i = 1; i <= ints[2]; i++, len++)
                aWeights(i) = data[len];
              hCurve = new Geom_BezierCurve(aPoles, aWeights);
            }
          }
          break;

        case BSPLINE:
          {
            Standard_Boolean periodic = Standard_False;
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            if ((ints[0]&4) != 0) periodic = Standard_True;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) len++;
            TColStd_Array1OfReal    aKnots(1, len);
            TColStd_Array1OfInteger aMults(1, len);
            aKnots(1) = data[0];
            aMults(1) = nmult = 1;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) {
                len++;
                aKnots(len) = data[i];
                aMults(len) = nmult = 1;
              } else {
                nmult++;
                aMults(len) = nmult;
              }
            len = ints[3];
            TColgp_Array1OfPnt aPoles(1, ints[2]);
            for (i = 1; i <= ints[2]; i++, len+=3)
              aPoles(i) = gp_Pnt(data[len], data[len+1], data[len+2]);
            if (rational == 0) {
              hCurve = new Geom_BSplineCurve(aPoles, aKnots, aMults,
                                             ints[1], periodic);
            } else {
              TColStd_Array1OfReal aWeights(1, ints[2]);
              for (i = 1; i <= ints[2]; i++, len++)
                aWeights(i) = data[len];
              hCurve = new Geom_BSplineCurve(aPoles, aWeights, aKnots,
                                             aMults, ints[1], periodic);
            }
          }
          break;

        case OFFSET:
          {
            egadsCurve *pcurve = (egadsCurve *) basis->blind;
            gp_Dir dir(data[0], data[1], data[2]);
            hCurve = new Geom_OffsetCurve(pcurve->handle, data[3], dir);
          }
          break;
      }
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      printf("                %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_makeGeometry)!\n", stat);
      return stat;
    }
    obj->oclass        = CURVE;
    obj->mtype         = mtype;
    egadsCurve *pcurve = new egadsCurve;
    pcurve->handle     = hCurve;
    pcurve->ref        = basis;
    pcurve->topFlg     = 1;
    pcurve->header     = NULL;
    pcurve->data       = NULL;
    pcurve->data_dot   = NULL;
    obj->blind         = pcurve;
    EG_getGeometry(obj, &i, &j, &ref, &pcurve->header, &pcurve->data);
    EG_getGeometryLen(obj, &i, &pcurve->dataLen);
    pcurve->trange[0]  = hCurve->FirstParameter();
    pcurve->trange[1]  = hCurve->LastParameter();
    EG_referenceObject(obj, context);
    if (basis != NULL) EG_referenceTopObj(basis, obj);

  } else {

    if ((mtype < PLANE) || (mtype > EXTRUSION)) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface mtype = %d (EG_makeGeometry)!\n",
               mtype);
      return EGADS_RANGERR;
    }
    if ((mtype == EXTRUSION) || (mtype == REVOLUTION)) {
      if (refGeom == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref is NULL (EG_makeGeometry)!\n",
                 mtype);
        return EGADS_NULLOBJ;
      }
      if (refGeom->oclass != CURVE) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref is %d (EG_makeGeometry)!\n",
                 mtype, refGeom->oclass);
        return EGADS_NOTGEOM;
      }
      basis = refGeom;
      if (basis->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref has no data (EG_makeGeometry)!\n",
                 mtype);
        return EGADS_NODATA;
      }
    }
    if ((mtype == OFFSET) || (mtype == TRIMMED)) {
      if (refGeom == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref is NULL (EG_makeGeometry)!\n",
                 mtype);
        return EGADS_NULLOBJ;
      }
      if (refGeom->oclass != SURFACE) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref is %d (EG_makeGeometry)!\n",
                 mtype, refGeom->oclass);
        return EGADS_NOTGEOM;
      }
      basis = refGeom;
      if (basis->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref has no data (EG_makeGeometry)!\n",
                 mtype);
        return EGADS_NODATA;
      }
    }

    Handle(Geom_Surface) hSurf;
    try {
      switch (mtype) {

        case PLANE:
          {
            gp_Pnt pntp(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax3 axi3(pntp, dirz, dirx);
            hSurf = new Geom_Plane(axi3);
          }
          break;

        case SPHERICAL:
          {
            gp_Pnt pnts(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            if (data[9] < 0.0) dirz.SetCoord(-dirz.X(), -dirz.Y(), -dirz.Z());
            gp_Ax3 axi3(pnts, dirz);
            axi3.SetXDirection(dirx);
            axi3.SetYDirection(diry);
            hSurf = new Geom_SphericalSurface(axi3, fabs(data[9]));
          }
          break;

        case CONICAL:
          {
            gp_Pnt pntc(data[0], data[1],  data[2]);
            gp_Dir dirx(data[3], data[4],  data[5]);
            gp_Dir diry(data[6], data[7],  data[8]);
            gp_Dir dirz(data[9], data[10], data[11]);
            gp_Ax3 axi3(pntc, dirz);
            axi3.SetXDirection(dirx);
            axi3.SetYDirection(diry);
            hSurf = new Geom_ConicalSurface(axi3, data[12], data[13]);
          }
          break;

        case CYLINDRICAL:
          {
            gp_Pnt pntc(data[0], data[1],  data[2]);
            gp_Dir dirx(data[3], data[4],  data[5]);
            gp_Dir diry(data[6], data[7],  data[8]);
            gp_Dir dirz(data[9], data[10], data[11]);
            gp_Ax3 axi3(pntc, dirz);
            axi3.SetXDirection(dirx);
            axi3.SetYDirection(diry);
            hSurf = new Geom_CylindricalSurface(axi3, data[12]);
          }
          break;

        case TOROIDAL:
          {
            gp_Pnt pntt(data[0], data[1],  data[2]);
            gp_Dir dirx(data[3], data[4],  data[5]);
            gp_Dir diry(data[6], data[7],  data[8]);
            gp_Dir dirz(data[9], data[10], data[11]);
            gp_Ax3 axi3(pntt, dirz);
            axi3.SetXDirection(dirx);
            axi3.SetYDirection(diry);
            hSurf = new Geom_ToroidalSurface(axi3, data[12], data[13]);
          }
          break;

        case BEZIER:
          {
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            TColgp_Array2OfPnt aPoles(1, ints[2], 1, ints[4]);
            len = 0;
            for (j = 1; j <= ints[4]; j++)
              for (i = 1; i <= ints[2]; i++, len+=3)
                aPoles(i, j) = gp_Pnt(data[len], data[len+1], data[len+2]);
            if (rational == 0) {
              hSurf = new Geom_BezierSurface(aPoles);
            } else {
              TColStd_Array2OfReal aWeights(1, ints[2], 1, ints[4]);
              for (j = 1; j <= ints[4]; j++)
                for (i = 1; i <= ints[2]; i++, len++)
                  aWeights(i, j) = data[len];
              hSurf = new Geom_BezierSurface(aPoles, aWeights);
            }
          }
          break;

        case BSPLINE:
          {
            Standard_Boolean uPeriodic = Standard_False;
            Standard_Boolean vPeriodic = Standard_False;
            rational = 0;
            if ((ints[0]&2) != 0) rational  = 1;
            if ((ints[0]&4) != 0) uPeriodic = Standard_True;
            if ((ints[0]&8) != 0) vPeriodic = Standard_True;
//          uKnots
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) len++;
            TColStd_Array1OfReal    uKnots(1, len);
            TColStd_Array1OfInteger uMults(1, len);
            uKnots(1) = data[0];
            uMults(1) = nmult = 1;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) {
                len++;
                uKnots(len) = data[i];
                uMults(len) = nmult = 1;
              } else {
                nmult++;
                uMults(len) = nmult;
              }
//          vKnots
            for (len = i = 1; i < ints[6]; i++)
              if (fabs(data[ints[3]+i]-data[ints[3]+i-1]) > KNACC) len++;
            TColStd_Array1OfReal    vKnots(1, len);
            TColStd_Array1OfInteger vMults(1, len);
            vKnots(1) = data[ints[3]];
            vMults(1) = nmult = 1;
            for (len = i = 1; i < ints[6]; i++)
              if (fabs(data[ints[3]+i]-data[ints[3]+i-1]) > KNACC) {
                len++;
                vKnots(len) = data[ints[3]+i];
                vMults(len) = nmult = 1;
              } else {
                nmult++;
                vMults(len) = nmult;
              }
            len = ints[3]+ints[6];
            TColgp_Array2OfPnt aPoles(1, ints[2], 1, ints[5]);
            for (j = 1; j <= ints[5]; j++)
              for (i = 1; i <= ints[2]; i++, len+=3)
                aPoles(i, j) = gp_Pnt(data[len], data[len+1], data[len+2]);
            if (rational == 0) {
              hSurf = new Geom_BSplineSurface(aPoles, uKnots, vKnots,
                                              uMults, vMults, ints[1], ints[4],
                                              uPeriodic, vPeriodic);
            } else {
              TColStd_Array2OfReal aWeights(1, ints[2], 1, ints[5]);
              for (j = 1; j <= ints[5]; j++)
                for (i = 1; i <= ints[2]; i++, len++)
                  aWeights(i, j) = data[len];
              hSurf = new Geom_BSplineSurface(aPoles, aWeights, uKnots, vKnots,
                                              uMults, vMults, ints[1], ints[4],
                                              uPeriodic, vPeriodic);
            }
          }
          break;

        case OFFSET:
          {
            egadsSurface *psurf = (egadsSurface *) basis->blind;
            hSurf = new Geom_OffsetSurface(psurf->handle, data[0]);
          }
          break;

        case TRIMMED:
          {
            egadsSurface *psurf = (egadsSurface *) basis->blind;
            hSurf = new Geom_RectangularTrimmedSurface(psurf->handle,
                                        data[0], data[1], data[2], data[3]);
          }
          break;

      case EXTRUSION:
          {
            gp_Dir dir(data[0], data[1],  data[2]);
            egadsCurve *pcurve = (egadsCurve *) basis->blind;
            hSurf = new Geom_SurfaceOfLinearExtrusion(pcurve->handle, dir);
          }
          break;

        case REVOLUTION:
          {
            gp_Pnt pnt(data[0], data[1], data[2]);
            gp_Dir dir(data[3], data[4], data[5]);
            gp_Ax1 axi1(pnt, dir);
            egadsCurve *pcurve = (egadsCurve *) basis->blind;
            hSurf = new Geom_SurfaceOfRevolution(pcurve->handle, axi1);
          }
          break;

      }
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      printf("                %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Surface = %d (EG_makeGeometry)!\n", stat);
      return stat;
    }
    obj->oclass         = SURFACE;
    obj->mtype          = mtype;
    egadsSurface *psurf = new egadsSurface;
    psurf->handle       = hSurf;
    psurf->ref          = basis;
    psurf->topFlg       = 1;
    psurf->header       = NULL;
    psurf->data         = NULL;
    psurf->data_dot     = NULL;
    obj->blind          = psurf;
    EG_getGeometry(obj, &i, &j, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(obj, &i, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    EG_referenceObject(obj, context);
    if (basis != NULL) EG_referenceTopObj(basis, obj);

  }

  *geom = obj;
  return EGADS_SUCCESS;
}


int
EG_getRange(const egObject *geom, double *range, int *periodic)
{
  int per;

  *periodic = 0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;

  if (geom->oclass == PCURVE) {

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    if (hCurve->IsPeriodic()) *periodic = 1;
    range[0] = ppcurv->trange[0];
    range[1] = ppcurv->trange[1];

  } else if (geom->oclass == CURVE) {

    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    if (hCurve->IsPeriodic()) *periodic = 1;
    range[0] = pcurve->trange[0];
    range[1] = pcurve->trange[1];

  } else if (geom->oclass == SURFACE) {

    egadsSurface *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurf = psurf->handle;
    per = 0;
    if (hSurf->IsUPeriodic()) per  = 1;
    if (hSurf->IsVPeriodic()) per |= 2;
    *periodic = per;
    range[0] = psurf->urange[0];
    range[1] = psurf->urange[1];
    range[2] = psurf->vrange[0];
    range[3] = psurf->vrange[1];

  } else if (geom->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) geom->blind;
    range[0] = pedge->trange[0];
    range[1] = pedge->trange[1];
    if (geom->mtype != DEGENERATE) {
      const egObject *ref    = pedge->curve;
      egadsCurve     *pcurve = (egadsCurve *) ref->blind;
      Handle(Geom_Curve) hCurve = pcurve->handle;
      if (hCurve->IsPeriodic()) *periodic = 1;
    }

  } else {

    egadsFace *pface = (egadsFace *) geom->blind;
    range[0] = pface->urange[0];
    range[1] = pface->urange[1];
    range[2] = pface->vrange[0];
    range[3] = pface->vrange[1];
    const egObject *ref   = pface->surface;
    egadsSurface   *psurf = (egadsSurface *) ref->blind;
    Handle(Geom_Surface) hSurf = psurf->handle;
    per = 0;
    if (hSurf->IsUPeriodic()) per  = 1;
    if (hSurf->IsVPeriodic()) per |= 2;
    *periodic = per;

  }

  return EGADS_SUCCESS;
}


int
EG_getRange(const egObject *geom, SurrealS<1> *range, int *periodic)
{
  int per, outLevel;

  *periodic = 0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;

  outLevel = EG_outLevel(geom);

  if (geom->oclass == PCURVE) {

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    if (hCurve->IsPeriodic()) *periodic = 1;
    range[0] = ppcurv->trange[0];
    range[1] = ppcurv->trange[1];

  } else if (geom->oclass == CURVE) {

    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    if (hCurve->IsPeriodic()) *periodic = 1;
    range[0] = pcurve->trange[0];
    range[1] = pcurve->trange[1];

  } else if (geom->oclass == SURFACE) {

    egadsSurface *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurf = psurf->handle;
    per = 0;
    if (hSurf->IsUPeriodic()) per  = 1;
    if (hSurf->IsVPeriodic()) per |= 2;
    *periodic = per;
    range[0] = psurf->urange[0];
    range[1] = psurf->urange[1];
    range[2] = psurf->vrange[0];
    range[3] = psurf->vrange[1];

  } else if (geom->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) geom->blind;
    if (pedge->filled == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge without data_dot (EG_getRange_dot)!\n");
      return EGADS_NODATA;
    }
    range[0] = pedge->trange_dot[0];
    range[1] = pedge->trange_dot[1];
    if (geom->mtype != DEGENERATE) {
      const egObject *ref    = pedge->curve;
      egadsCurve     *pcurve = (egadsCurve *) ref->blind;
      Handle(Geom_Curve) hCurve = pcurve->handle;
      if (hCurve->IsPeriodic()) *periodic = 1;
    }

  } else {

    egadsFace *pface = (egadsFace *) geom->blind;
    range[0] = pface->urange[0];
    range[1] = pface->urange[1];
    range[2] = pface->vrange[0];
    range[3] = pface->vrange[1];
    const egObject *ref   = pface->surface;
    egadsSurface   *psurf = (egadsSurface *) ref->blind;
    Handle(Geom_Surface) hSurf = psurf->handle;
    per = 0;
    if (hSurf->IsUPeriodic()) per  = 1;
    if (hSurf->IsVPeriodic()) per |= 2;
    *periodic = per;

  }

  return EGADS_SUCCESS;
}


int
EG_getRange_dot(const egObject *geom, double *range, double *range_dot,
                int *periodic)
{
  int stat;
  SurrealS<1> rangeS[4];

  stat = EG_getRange(geom, rangeS, periodic);
  if (stat != EGADS_SUCCESS) return stat;

  if (geom->oclass == PCURVE ||
      geom->oclass == CURVE ||
      geom->oclass == EDGE) {

    range[0] = rangeS[0].value();
    range[1] = rangeS[1].value();

    range_dot[0] = rangeS[0].deriv();
    range_dot[1] = rangeS[1].deriv();

  } else if (geom->oclass == SURFACE || geom->oclass == FACE) {

    range[0] = rangeS[0].value();
    range[1] = rangeS[1].value();
    range[2] = rangeS[2].value();
    range[3] = rangeS[3].value();

    range_dot[0] = rangeS[0].deriv();
    range_dot[1] = rangeS[1].deriv();
    range_dot[2] = rangeS[2].deriv();
    range_dot[3] = rangeS[3].deriv();
  }

  return EGADS_SUCCESS;
}


int
EG_curvature(const egObject *geom, const double *param, double *result)
{
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  int outLevel = EG_outLevel(geom);

  if (geom->oclass == PCURVE) {

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    Geom2dLProp_CLProps2d aProp(hCurve, param[0], 2, Precision::Confusion());
    if (aProp.IsTangentDefined()) {
      gp_Dir2d tang;

      aProp.Tangent(tang);
      result[0] = aProp.Curvature();
      result[1] = tang.X();
      result[2] = tang.Y();
    } else {
      for (int i = 0; i < 3; i++) result[i] = 0.0;
      return EGADS_DEGEN;
    }

  } else if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) {

    // 1D -- curves & Edges
    Handle(Geom_Curve) hCurve;
    if (geom->oclass == CURVE) {
      egadsCurve *pcurve = (egadsCurve *) geom->blind;
      hCurve = pcurve->handle;
    } else {
      egadsEdge *pedge = (egadsEdge *) geom->blind;
      egObject   *curv = pedge->curve;
      if (curv == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Object for Edge (EG_curvature)!\n");
        return EGADS_NULLOBJ;
      }
      egadsCurve *pcurve = (egadsCurve *) curv->blind;
      if (pcurve == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Data for Edge (EG_curvature)!\n");
        return EGADS_NODATA;
      }
      hCurve = pcurve->handle;
    }
    GeomLProp_CLProps aProp(hCurve, param[0], 2, Precision::Confusion());
    if (aProp.IsTangentDefined()) {
      gp_Dir tang;

      aProp.Tangent(tang);
      result[0] = aProp.Curvature();
      result[1] = tang.X();
      result[2] = tang.Y();
      result[3] = tang.Z();
    } else {
      for (int i = 0; i < 4; i++) result[i] = 0.0;
      return EGADS_DEGEN;
    }

  } else {

    // 2D -- surfaces & Faces
    Handle(Geom_Surface) hSurface;
    if (geom->oclass == SURFACE) {
      egadsSurface *psurf = (egadsSurface *) geom->blind;
      hSurface = psurf->handle;
    } else {
      egadsFace *pface = (egadsFace *) geom->blind;
      egObject *surf = pface->surface;
      if (surf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Object for Face (EG_curvature)!\n");
        return EGADS_NULLOBJ;
      }
      egadsSurface *psurf = (egadsSurface *) surf->blind;
      if (psurf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Data for Face (EG_curvature)!\n");
        return EGADS_NODATA;
      }
      hSurface = psurf->handle;
    }
    GeomLProp_SLProps aProp(hSurface, param[0], param[1], 2,
                            Precision::Angular());
    if (aProp.IsCurvatureDefined()) {
      gp_Dir MaxD, MinD;

      aProp.CurvatureDirections(MaxD, MinD);
      result[0] = aProp.MaxCurvature();
      result[1] = MaxD.X();
      result[2] = MaxD.Y();
      result[3] = MaxD.Z();
      result[4] = aProp.MinCurvature();
      result[5] = MinD.X();
      result[6] = MinD.Y();
      result[7] = MinD.Z();
      if ((geom->oclass == FACE) && (geom->mtype == SREVERSE)) {
        result[0] = -result[0];
        result[4] = -result[4];
      }
    } else {
      for (int i = 0; i < 8; i++) result[i] = 0.0;
      return EGADS_DEGEN;
    }

  }

  return EGADS_SUCCESS;
}


int
EG_evaluate(const egObject *geom, /*@null@*/ const double *param, double *result)
{
  int            stat, outLevel, per, our = 1;
  double         range[4];
  const egObject *ref;
  gp_Pnt         P0;
  gp_Vec         V1, V2, U1, U2, UV;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != NODE)  && (geom->oclass != PCURVE)  &&
      (geom->oclass != CURVE) && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)  && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(geom);

  // special Node section
  if (geom->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) geom->blind;
    for (int i = 0; i < 3; i++) result[i] = pnode->xyz[i];
    return EGADS_SUCCESS;
  }
  if (param == NULL)               return EGADS_NODATA;

  // use our evaluators if the data exists and we are not a periodic BSpline
  ref = geom;
  if (geom->oclass == EDGE) {
    if (geom->mtype == DEGENERATE) {
      if (outLevel > 0)
        printf(" EGADS Warning: Degenerate Edge (EG_evaluate)!\n");
      return EGADS_DEGEN;
    }
    egadsEdge *pedge = (egadsEdge *) geom->blind;
    ref = pedge->curve;
    if (ref == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No curve Object for Edge (EG_evaluate)!\n");
      return EGADS_NULLOBJ;
    }
  }
  if (geom->oclass == FACE) {
    egadsFace *pface = (egadsFace *) geom->blind;
    ref = pface->surface;
    if (ref == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No surface Object for Face (EG_evaluate)!\n");
      return EGADS_NULLOBJ;
    }
  }
  if (ref->blind == NULL) return EGADS_NODATA;
  while (ref->mtype == TRIMMED) {
    if (ref->oclass == PCURVE) {
      egadsPCurve *ppcurv = (egadsPCurve *) ref->blind;
      ref = ppcurv->ref;
    } else if (ref->oclass == CURVE) {
      egadsCurve *pcurve = (egadsCurve *) ref->blind;
      ref = pcurve->ref;
    } else if (ref->oclass == SURFACE) {
      egadsSurface *psurf = (egadsSurface *) ref->blind;
      ref = psurf->ref;
    }
  }
  if (ref->oclass == PCURVE) {
    egadsPCurve *ppcurv = (egadsPCurve *) ref->blind;
    if (ppcurv->data == NULL) our = 0;
  } else if (ref->oclass == CURVE) {
    egadsCurve *pcurve = (egadsCurve *) ref->blind;
    if (pcurve->data == NULL) our = 0;
  } else if (ref->oclass == SURFACE) {
    egadsSurface *psurf = (egadsSurface *) ref->blind;
    if (psurf->data == NULL) our = 0;
  } else {
    if (outLevel > 0)
      printf(" EGADS Warning: Geom Object class = %d (EG_evaluate)!\n",
             ref->oclass);
    return EGADS_NOTGEOM;
  }
  if ((ref->mtype == BSPLINE) && (our == 1)) {
    stat = EG_getRange(ref, range, &per);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Warning: getRange = %d (EG_evaluate)!\n", stat);
      return stat;
    }
    if (per != 0) our = 0;
  }
  if (our == 1) return EG_evaluateGeom(ref, param, result);

  // use OpenCASCADE
  if (geom->oclass == PCURVE) {

    gp_Pnt2d P2d;
    gp_Vec2d V12d, V22d;

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    hCurve->D2(*param, P2d, V12d, V22d);
    result[0] = P2d.X();
    result[1] = P2d.Y();
    result[2] = V12d.X();
    result[3] = V12d.Y();
    result[4] = V22d.X();
    result[5] = V22d.Y();

  } else if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) {

    // 1D -- curves & Edges
    egadsCurve *pcurve = (egadsCurve *) ref->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    hCurve->D2(*param, P0, V1, V2);

    result[0] = P0.X();
    result[1] = P0.Y();
    result[2] = P0.Z();
    result[3] = V1.X();
    result[4] = V1.Y();
    result[5] = V1.Z();
    result[6] = V2.X();
    result[7] = V2.Y();
    result[8] = V2.Z();

  } else {

    // 2D -- surfaces & Faces
    egadsSurface *psurf = (egadsSurface *) ref->blind;
    Handle(Geom_Surface) hSurface = psurf->handle;
    hSurface->D2(param[0], param[1], P0, U1, V1, U2, V2, UV);

    result[ 0] = P0.X();
    result[ 1] = P0.Y();
    result[ 2] = P0.Z();
    result[ 3] = U1.X();
    result[ 4] = U1.Y();
    result[ 5] = U1.Z();
    result[ 6] = V1.X();
    result[ 7] = V1.Y();
    result[ 8] = V1.Z();
    result[ 9] = U2.X();
    result[10] = U2.Y();
    result[11] = U2.Z();
    result[12] = UV.X();
    result[13] = UV.Y();
    result[14] = UV.Z();
    result[15] = V2.X();
    result[16] = V2.Y();
    result[17] = V2.Z();

  }

  return EGADS_SUCCESS;
}


int
EG_evaluate(const egObject *geom, /*@null@*/ const SurrealS<1> *param,
            SurrealS<1> *result)
{
  int            i, stat, outLevel, len, per, our = 1;
  double         range[4];
  const egObject *ref;
  SurrealS<1>    data[18];

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != NODE)  && (geom->oclass != PCURVE)  &&
      (geom->oclass != CURVE) && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)  && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(geom);

  // special Node section
  if (geom->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) geom->blind;
    if (pnode->filled == 0) return EGADS_NODATA;
    for (i = 0; i < 3; i++) {
      result[i] = pnode->xyz_dot[i];
    }
    return EGADS_SUCCESS;
  }

  // get the return length
  len = 18;
  if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) len = 9;
  if  (geom->oclass == PCURVE) len = 6;

  // use our evaluators if the data exists and we are not a periodic BSpline
  ref = geom;
  if (geom->oclass == EDGE) {
    if (geom->mtype == DEGENERATE) {
      if (outLevel > 0)
        printf(" EGADS Warning: Degenerate Edge (EG_evaluate)!\n");
      return EGADS_DEGEN;
    }
    egadsEdge *pedge = (egadsEdge *) geom->blind;
    ref = pedge->curve;
    if (ref == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No curve Object for Edge (EG_evaluate)!\n");
      return EGADS_NULLOBJ;
    }
  }
  if (geom->oclass == FACE) {
    egadsFace *pface = (egadsFace *) geom->blind;
    ref = pface->surface;
    if (ref == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No surface Object for Face (EG_evaluate)!\n");
      return EGADS_NULLOBJ;
    }
  }
  if (ref->blind == NULL) return EGADS_NODATA;
  if (ref->oclass == PCURVE) {
    egadsPCurve *ppcurv = (egadsPCurve *) ref->blind;
    if (ppcurv->data == NULL) our = 0;
  } else if (ref->oclass == CURVE) {
    egadsCurve *pcurve = (egadsCurve *) ref->blind;
    if (pcurve->data == NULL) our = 0;
  } else if (ref->oclass == SURFACE) {
    egadsSurface *psurf = (egadsSurface *) ref->blind;
    if (psurf->data == NULL) our = 0;
  } else {
    if (outLevel > 0)
      printf(" EGADS Warning: Geom Object class = %d (EG_evaluate)!\n",
             ref->oclass);
    return EGADS_NOTGEOM;
  }
  if ((ref->mtype == BSPLINE) && (our == 1)) {
    stat = EG_getRange(ref, range, &per);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Warning: getRange = %d (EG_evaluate)!\n", stat);
      return stat;
    }
    if (per != 0) our = 0;
  }
  if (our == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Periodic BSpline (EG_evaluate)!\n");
    return EGADS_GEOMERR;
  }

  stat = EG_evaluateGeom(ref, param, data);
  if (stat != EGADS_SUCCESS) return stat;

  for (i = 0; i < len; i++) {
    result[i] = data[i];
  }
  return EGADS_SUCCESS;
}


int
EG_evaluate_dot(const egObject *geom, /*@null@*/ const double *param,
                /*@null@*/ const double *param_dot,
                double *result, double *result_dot)
{
  int            i, stat, outLevel, len, per, our = 1;
  double         range[4];
  const egObject *ref;
  SurrealS<1>    data[18], paramS[2];

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != NODE)  && (geom->oclass != PCURVE)  &&
      (geom->oclass != CURVE) && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)  && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(geom);

  // special Node section
  if (geom->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) geom->blind;
    if (pnode->filled == 0) return EGADS_NODATA;
    for (i = 0; i < 3; i++) {
      result[i]     = pnode->xyz_dot[i].value();
      result_dot[i] = pnode->xyz_dot[i].deriv();
    }
    return EGADS_SUCCESS;
  }
  if (param == NULL)               return EGADS_NODATA;

  // get the return length
  len = 18;
  if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) len = 9;
  if  (geom->oclass == PCURVE) len = 6;

  // use our evaluators if the data exists and we are not a periodic BSpline
  ref = geom;
  if (geom->oclass == EDGE) {
    if (geom->mtype == DEGENERATE) {
      if (outLevel > 0)
        printf(" EGADS Warning: Degenerate Edge (EG_evaluate_dot)!\n");
      return EGADS_DEGEN;
    }
    egadsEdge *pedge = (egadsEdge *) geom->blind;
    ref = pedge->curve;
    if (ref == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No curve Object for Edge (EG_evaluate_dot)!\n");
      return EGADS_NULLOBJ;
    }
  }
  if (geom->oclass == FACE) {
    egadsFace *pface = (egadsFace *) geom->blind;
    ref = pface->surface;
    if (ref == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No surface Object for Face (EG_evaluate_dot)!\n");
      return EGADS_NULLOBJ;
    }
  }
  if (ref->blind == NULL) return EGADS_NODATA;
  if (ref->oclass == PCURVE) {
    egadsPCurve *ppcurv = (egadsPCurve *) ref->blind;
    if (ppcurv->data == NULL) our = 0;
  } else if (ref->oclass == CURVE) {
    egadsCurve *pcurve = (egadsCurve *) ref->blind;
    if (pcurve->data == NULL) our = 0;
  } else if (ref->oclass == SURFACE) {
    egadsSurface *psurf = (egadsSurface *) ref->blind;
    if (psurf->data == NULL) our = 0;
  } else {
    if (outLevel > 0)
      printf(" EGADS Warning: Geom Object class = %d (EG_evaluate_dot)!\n",
             ref->oclass);
    return EGADS_NOTGEOM;
  }
  if ((ref->mtype == BSPLINE) && (our == 1)) {
    stat = EG_getRange(ref, range, &per);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Warning: getRange = %d (EG_evaluate_dot)!\n", stat);
      return stat;
    }
    if (per != 0) our = 0;
  }
  if (our == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Periodic BSpline (EG_evaluate_dot)!\n");
    return EGADS_GEOMERR;
  }

  if (ref->oclass == PCURVE || ref->oclass == CURVE) {
    paramS[0] = param[0];
    if (param_dot != NULL) {
      paramS[0].deriv() = param_dot[0];
    }
  } else {
    paramS[0] = param[0];
    paramS[1] = param[1];
    if (param_dot != NULL) {
      paramS[0].deriv() = param_dot[0];
      paramS[1].deriv() = param_dot[1];
    }
  }

  stat = EG_evaluateGeom(ref, paramS, data);
  if (stat != EGADS_SUCCESS) return stat;

  for (i = 0; i < len; i++) {
    result[i]     = data[i].value();
    result_dot[i] = data[i].deriv();
  }
  return EGADS_SUCCESS;
}


static void
EG_nearestPCurve(Handle_Geom2d_Curve hCurve, const double *coor,
                 double tmin, double tmax, int flag, double *t, double *uv)
{
  int      i;
  double   a, b, tx, pw[2];
  gp_Pnt2d pnt;
  gp_Vec2d t1, t2;
  static double ratios[5] = {0.02, 0.25, 0.5, 0.75, 0.98};

  // sample and pick closest
  if (flag == 0) {
    b = 0.0;
    for (i = 0; i < 5; i++) {
      tx = (1.0-ratios[i])*tmin + ratios[i]*tmax;
      hCurve->D0(tx, pnt);
      a = (pnt.X()-coor[0])*(pnt.X()-coor[0]) +
          (pnt.Y()-coor[1])*(pnt.Y()-coor[1]);
      if (i == 0) {
        *t = tx;
        b  = a;
      } else {
        if (a < b) {
          *t = tx;
          b  = a;
        }
      }
    }
  }

  // netwon-raphson from picked position
  for (i = 0; i < 20; i++) {
    if ((*t < tmin) || (*t > tmax)) break;
    hCurve->D2(*t, pnt, t1, t2);
    pw[0] = pnt.X() - coor[0];
    pw[1] = pnt.Y() - coor[1];
    b     = -( pw[0]*t1.X() +  pw[1]*t1.Y());
    a     =  (t1.X()*t1.X() + t1.Y()*t1.Y()) +
             ( pw[0]*t2.X() +  pw[1]*t2.Y());
    if (a == 0.0) break;
    b  /= a;
//  if (fabs(b) < 1.e-10*(tmax-tmin)) break;
    *t += b;
  }
  if (*t < tmin) *t = tmin;
  if (*t > tmax) *t = tmax;

  hCurve->D0(*t, pnt);
  uv[0] = pnt.X();
  uv[1] = pnt.Y();
}


void
EG_nearestCurve(Handle_Geom_Curve hCurve, const double *coor,
                double tmin, double tmax, int flag, double *t, double *xyz)
{
  int    i;
  double a, b, tx, pw[3];
  gp_Pnt pnt;
  gp_Vec t1, t2;
  static double ratios[5] = {0.02, 0.25, 0.5, 0.75, 0.98};

  // sample and pick closest
  if (flag == 0) {
    b = 0.0;
    for (i = 0; i < 5; i++) {
      tx = (1.0-ratios[i])*tmin + ratios[i]*tmax;
      hCurve->D0(tx, pnt);
      a = (pnt.X()-coor[0])*(pnt.X()-coor[0]) +
          (pnt.Y()-coor[1])*(pnt.Y()-coor[1]) +
          (pnt.Z()-coor[2])*(pnt.Z()-coor[2]);
      if (i == 0) {
        *t = tx;
        b  = a;
      } else {
        if (a < b) {
          *t = tx;
          b  = a;
        }
      }
    }
  }

  // netwon-raphson from picked position
  for (i = 0; i < 20; i++) {
    if ((*t < tmin) || (*t > tmax)) break;
    hCurve->D2(*t, pnt, t1, t2);
    pw[0] = pnt.X() - coor[0];
    pw[1] = pnt.Y() - coor[1];
    pw[2] = pnt.Z() - coor[2];
    b     = -( pw[0]*t1.X() +  pw[1]*t1.Y() +  pw[2]*t1.Z());
    a     =  (t1.X()*t1.X() + t1.Y()*t1.Y() + t1.Z()*t1.Z()) +
             ( pw[0]*t2.X() +  pw[1]*t2.Y() +  pw[2]*t2.Z());
    if (a == 0.0) break;
    b  /= a;
//  if (fabs(b) < 1.e-10*(tmax-tmin)) break;
    *t += b;
  }
  if (*t < tmin) *t = tmin;
  if (*t > tmax) *t = tmax;

  hCurve->D0(*t, pnt);
  xyz[0] = pnt.X();
  xyz[1] = pnt.Y();
  xyz[2] = pnt.Z();
}


static int
EG_nearestSurface(Handle_Geom_Surface hSurface, double *range,
                  const double *point, int flag, double *uv, double *coor)
{
  int    i, j, count;
  gp_Pnt P0;
  gp_Vec V1, V2, U1, U2, UV;
  double a00, a10, a11, b0, b1, det, dist, ldist, dx[3], uvs[2];
  static double ratios[5] = {0.02, 0.25, 0.5, 0.75, 0.98};

  if (flag == 0) {
    // find candidate starting point
    ldist = 1.e308;
    for (j = 0; j < 5; j++) {
      uvs[1] = (1.0-ratios[j])*range[2] + ratios[j]*range[3];
      for (i = 0; i < 5; i++) {
        uvs[0] = (1.0-ratios[i])*range[0] + ratios[i]*range[1];
        hSurface->D0(uvs[0], uvs[1], P0);
        dx[0] = P0.X() - point[0];
        dx[1] = P0.Y() - point[1];
        dx[2] = P0.Z() - point[2];
        dist  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        if (dist < ldist) {
          ldist = dist;
          uv[0] = uvs[0];
          uv[1] = uvs[1];
        }
      }
    }
  }

  // newton iteration
  ldist = uvs[0] = uvs[1] = 0.0;
  for (count = 0; count < 15; count++) {
    hSurface->D2(uv[0], uv[1], P0, U1, V1, U2, V2, UV);
    dx[0] = P0.X() - point[0];
    dx[1] = P0.Y() - point[1];
    dx[2] = P0.Z() - point[2];
    dist  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    if (dist < Precision::Confusion()) break;
    if (count != 0) {
      if (fabs(dist-ldist) < Precision::Confusion()) break;
      if (dist > ldist) {
        uv[0] = uvs[0];
        uv[1] = uvs[1];
        hSurface->D0(uv[0], uv[1], P0);
        coor[0] = P0.X();
        coor[1] = P0.Y();
        coor[2] = P0.Z();
        return EGADS_EMPTY;
      }
    }

    b0  = -dx[0]*U1.X() -  dx[1]*U1.Y() -  dx[2]*U1.Z();
    b1  = -dx[0]*V1.X() -  dx[1]*V1.Y() -  dx[2]*V1.Z();
    a00 = U1.X()*U1.X() + U1.Y()*U1.Y() + U1.Z()*U1.Z() +
           dx[0]*U2.X() +  dx[1]*U2.Y() +  dx[2]*U2.Z();
    a10 = U1.X()*V1.X() + U1.Y()*V1.Y() + U1.Z()*V1.Z() +
           dx[0]*UV.X() +  dx[1]*UV.Y() +  dx[2]*UV.Z();
    a11 = V1.X()*V1.X() + V1.Y()*V1.Y() + V1.Z()*V1.Z() +
           dx[0]*V2.X() +  dx[1]*V2.Y() +  dx[2]*V2.Z();

    det    = a00*a11 - a10*a10;
    if (det == 0.0) return EGADS_DEGEN;
    det    = 1.0/det;
    uvs[0] = uv[0];
    uvs[1] = uv[1];
    uv[0] += det*(b0*a11 - b1*a10);
    uv[1] += det*(b1*a00 - b0*a10);
    ldist  = dist;
//  printf("   %d: %lf %lf   %le\n", count, uv[0], uv[1], ldist);
  }

  hSurface->D0(uv[0], uv[1], P0);
  coor[0] = P0.X();
  coor[1] = P0.Y();
  coor[2] = P0.Z();
  if (count == 15) return EGADS_EMPTY;

  return EGADS_SUCCESS;
}


static int
EG_invEvalClip(const egObject *geom, double *xyz, double *param, double *result)
{
/*
  int            stat, outLevel;
  Standard_Real  tol, period, t, u, v, srange[4], coor[3], dist2;

  egadsFace      *pface = (egadsFace *)    geom->blind;
  const egObject *ref   = pface->surface;
  egadsSurface   *psurf = (egadsSurface *) ref->blind;

  tol = BRep_Tool::Tolerance(pface->face);
  gp_Pnt2d pnt2d(param[0], param[1]);
  BRepClass_FaceClassifier aClassifier(pface->face, pnt2d, tol);
  if (aClassifier.State() != TopAbs_OUT) return EGADS_SUCCESS;

  Handle(Geom_Surface) hSurface = psurf->handle;
  outLevel  = EG_outLevel(geom);
  srange[0] = psurf->urange[0];
  srange[1] = psurf->urange[1];
  srange[2] = psurf->vrange[0];
  srange[3] = psurf->vrange[1];

  dist2 = 1.e308;
  u     = param[0];
  v     = param[1];
  gp_Pnt pnt(xyz[0], xyz[1], xyz[2]);
  gp_Pnt pnts(xyz[0], xyz[1], xyz[2]);
  gp_Pnt pntt(xyz[0], xyz[1], xyz[2]);
  TopExp_Explorer ExpW;
  for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
    TopoDS_Shape shapw = ExpW.Current();
    TopoDS_Wire  wire  = TopoDS::Wire(shapw);
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shape = ExpWE.Current();
      TopoDS_Edge  wedge = TopoDS::Edge(shape);
      if (BRep_Tool::Degenerated(wedge)) continue;
      Standard_Real t1, t2;
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(wedge, t1, t2);
      GeomAPI_ProjectPointOnCurve projPnt(pnts, hCurve);
      if (projPnt.NbPoints() == 0) {
        EG_nearestCurve(hCurve, xyz, t1, t2, 0, &t, result);
        pnt.SetX(result[0]);
        pnt.SetY(result[1]);
        pnt.SetZ(result[2]);
      } else {
        pnt = projPnt.NearestPoint();
        t   = projPnt.LowerDistanceParameter();
      }
      if ((t < t1) || (t > t2)) {
        EG_nearestCurve(hCurve, xyz, t1, t2, 0, &t, result);
        pnt.SetX(result[0]);
        pnt.SetY(result[1]);
        pnt.SetZ(result[2]);
      }
      double d = (pnts.X()-pnt.X())*(pnts.X()-pnt.X()) +
                 (pnts.Y()-pnt.Y())*(pnts.Y()-pnt.Y()) +
                 (pnts.Z()-pnt.Z())*(pnts.Z()-pnt.Z());
      if (d < dist2) {
        pntt  = pnt;
        dist2 = d;
      }
    }
  }

  GeomAPI_ProjectPointOnSurf projPnt(pntt, hSurface);
  if (!projPnt.IsDone()) {
    coor[0] = pntt.X();
    coor[1] = pntt.Y();
    coor[2] = pntt.Z();
    stat = EG_nearestSurface(hSurface, srange, coor, 0, param, result);
    if (stat == EGADS_DEGEN) {
      if (outLevel > 0)
        printf(" EGADS Warning: Face Proj Incomplete - DEGEN (EG_invEvaluate)!\n");
      return stat;
    } else if (stat == EGADS_EMPTY) {
      if (outLevel > 1)
        printf(" EGADS Warning: Face Proj Incomplete (EG_invEvaluate)!\n");
    }
    u = param[0];
    v = param[1];
    pnt.SetX(result[0]);
    pnt.SetY(result[1]);
    pnt.SetZ(result[2]);
  } else {
    if (projPnt.NbPoints() == 0) {
      if (outLevel > 0)
        printf(" EGADS Warning: No projection on Face (EG_invEvaluate)!\n");
      return EGADS_NOTFOUND;
    }
    pnt = projPnt.NearestPoint();
    projPnt.LowerDistanceParameters(u, v);
    gp_Pnt pnt1;
    hSurface->D0(u,v, pnt1);
    if (pnt1.Distance(pnt) > 1.e-7) {
      stat = EG_nearestSurface(hSurface, srange, xyz, 0, param, result);
      if (stat == EGADS_DEGEN) {
        if (outLevel > 0)
          printf(" EGADS Warning: Surf Proj Incomplete - DEGEN (EG_invEvaluate)!\n");
        return stat;
      } else if (stat == EGADS_EMPTY) {
        if (outLevel > 1)
          printf(" EGADS Warning: Surf Proj Incomplete (EG_invEvaluate)!\n");
      }
      u = param[0];
      v = param[1];
      pnt.SetX(result[0]);
      pnt.SetY(result[1]);
      pnt.SetZ(result[2]);
    }
  }

  if (hSurface->IsUPeriodic()) {
    period = hSurface->UPeriod();
    if ((u+PARAMACC < srange[0]) || (u-PARAMACC > srange[1])) {
      if (period != 0.0)
        if (u+PARAMACC < srange[0]) {
          if (u+period-PARAMACC < srange[1]) u += period;
        } else {
          if (u-period+PARAMACC > srange[0]) u -= period;
        }
    }
  }
  if (hSurface->IsVPeriodic()) {
    period = hSurface->VPeriod();
    if ((v+PARAMACC < srange[2]) || (v-PARAMACC > srange[3])) {
      if (period != 0.0)
        if (v+PARAMACC < srange[2]) {
          if (v+period-PARAMACC < srange[3]) v += period;
        } else {
          if (v-period+PARAMACC > srange[2]) v -= period;
        }
    }
  }

  result[0] = pnt.X();
  result[1] = pnt.Y();
  result[2] = pnt.Z();
  param[0]  = u;
  param[1]  = v;
*/
  int    stat, per;
  double pt[3], uvs[2], srange[4], period;

  egadsFace      *pface = (egadsFace *)    geom->blind;
  const egObject *ref   = pface->surface;

  stat = EG_inFaceX(geom, param, pt, uvs);
  if (stat != EGADS_OUTSIDE) return stat;

  /*  printf(" Info: Point labelled outside!\n");  */
  param[0]  = uvs[0];
  param[1]  = uvs[1];
  result[0] = pt[0];
  result[1] = pt[1];
  result[2] = pt[2];

  stat = EG_getRange(ref, srange, &per);
  if (stat != EGADS_SUCCESS) return stat;

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

  return EGADS_SUCCESS;
}


int
EG_invEvaluate(const egObject *geom, double *xyz, double *param, double *result)
{
  int            outLevel, stat, per, our = 1;
  const egObject *ref;
  Standard_Real  period, t, u, v, range[4], srange[4], tol = 0.0;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(geom);

  // use our evaluators if the data exists and we are not a periodic BSpline
  ref = geom;
  if (geom->oclass == EDGE) {
    if (geom->mtype == DEGENERATE) {
      if (outLevel > 0)
        printf(" EGADS Warning: Degenerate Edge (EG_invEvaluate)!\n");
      return EGADS_DEGEN;
    }
    egadsEdge *pedge = (egadsEdge *) geom->blind;
    ref = pedge->curve;
    if (ref == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No curve Object for Edge (EG_invEvaluate)!\n");
      return EGADS_NULLOBJ;
    }
  }
  if (geom->oclass == FACE) {
    egadsFace *pface = (egadsFace *) geom->blind;
    tol = BRep_Tool::Tolerance(pface->face);
    ref = pface->surface;
    if (ref == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No surface Object for Face (EG_invEvaluate)!\n");
      return EGADS_NULLOBJ;
    }
  }
  if (ref->blind == NULL) return EGADS_NODATA;
  while (ref->mtype == TRIMMED) {
    if (ref->oclass == PCURVE) {
      egadsPCurve *ppcurv = (egadsPCurve *) ref->blind;
      ref = ppcurv->ref;
    } else if (ref->oclass == CURVE) {
      egadsCurve *pcurve = (egadsCurve *) ref->blind;
      ref = pcurve->ref;
    } else if (ref->oclass == SURFACE) {
      egadsSurface *psurf = (egadsSurface *) ref->blind;
      ref = psurf->ref;
    }
  }
  if (ref->oclass == PCURVE) {
    egadsPCurve *ppcurv = (egadsPCurve *) ref->blind;
    if (ppcurv->data == NULL) our = 0;
  } else if (ref->oclass == CURVE) {
    egadsCurve *pcurve = (egadsCurve *) ref->blind;
    if (pcurve->data == NULL) our = 0;
  } else if (ref->oclass == SURFACE) {
    egadsSurface *psurf = (egadsSurface *) ref->blind;
    if (psurf->data == NULL) our = 0;
  } else {
    if (outLevel > 0)
      printf(" EGADS Warning: Geom Object class = %d (EG_invEvaluate)!\n",
             ref->oclass);
    return EGADS_NOTGEOM;
  }
  if (ref->mtype == BSPLINE) {
    stat = EG_getRange(ref, srange, &per);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Warning: EG_getRange = %d (EG_invEvaluate)!\n", stat);
      return stat;
    }
    if (per != 0) our = 0;
  }
  stat = EG_getRange(geom, range, &per);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: getRange = %d (EG_invEvaluate)!\n", stat);
    return stat;
  }
  if (our == 1)
    if (ref->oclass == PCURVE) {
      return EG_invEvaGeomLimits(ref, NULL,  xyz, param, 0.0, result);
    } else if (geom->oclass == FACE) {
      stat = EG_invEvaGeomLimits(ref, range, xyz, param, tol, result);
      if (stat != EGADS_SUCCESS) return stat;
      return EG_invEvalClip(geom, xyz, param, result);
    } else {
      return EG_invEvaGeomLimits(ref, range, xyz, param, tol, result);
    }

  if (geom->oclass == PCURVE) {

    // 2D on PCurves
    gp_Pnt2d pnt(xyz[0], xyz[1]);

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    Geom2dAPI_ProjectPointOnCurve projPnt(pnt, hCurve);
    if (projPnt.NbPoints() == 0) {
      EG_nearestPCurve(hCurve, xyz, range[0], range[1], 0, &t, result);
      pnt.SetX(result[0]);
      pnt.SetY(result[1]);
    } else {
      pnt = projPnt.NearestPoint();
      t   = projPnt.LowerDistanceParameter();
    }
    if (hCurve->IsPeriodic()) {
      period = hCurve->Period();
      if ((t+PARAMACC < range[0]) || (t-PARAMACC > range[1]))
        if (period != 0.0)
          if (t+PARAMACC < range[0]) {
            if (t+period-PARAMACC < range[1]) t += period;
          } else {
            if (t-period+PARAMACC > range[0]) t -= period;
          }
    }
    result[0] = pnt.X();
    result[1] = pnt.Y();
    *param    = t;
    return EGADS_SUCCESS;
  }

  // make the point
  gp_Pnt pnt(xyz[0], xyz[1], xyz[2]);

  if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) {

    // 1D -- curves & Edges
    egadsCurve *pcurve = (egadsCurve *) ref->blind;
    if (pcurve == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No curve Data for Edge (EG_invEvaluate)!\n");
      return EGADS_NODATA;
    }
    Handle(Geom_Curve) hCurve = pcurve->handle;
    if (geom->oclass != CURVE) {
      srange[0] = hCurve->FirstParameter();
      srange[1] = hCurve->LastParameter();
    } else {
      srange[0] = range[0];
      srange[1] = range[1];
    }

    GeomAPI_ProjectPointOnCurve projPnt(pnt, hCurve);
    if (projPnt.NbPoints() == 0) {
      EG_nearestCurve(hCurve, xyz, range[0], range[1], 0, &t, result);
      pnt.SetX(result[0]);
      pnt.SetY(result[1]);
      pnt.SetZ(result[2]);
    } else {
      pnt = projPnt.NearestPoint();
      t   = projPnt.LowerDistanceParameter();
    }

    if (hCurve->IsPeriodic()) {
      period = hCurve->Period();
      if ((t+PARAMACC < srange[0]) || (t-PARAMACC > srange[1])) {
        if (period != 0.0)
          if (t+PARAMACC < srange[0]) {
            if (t+period-PARAMACC < srange[1]) t += period;
          } else {
            if (t-period+PARAMACC > srange[0]) t -= period;
          }
      }
    }

    /* clip it? */
    if (geom->oclass == EDGE)
      if ((t < range[0]) || (t > range[1])) {
/*      if (t < range[0]) t = range[0];
        if (t > range[1]) t = range[1];
        hCurve->D0(t, pnt);  */
        EG_nearestCurve(hCurve, xyz, range[0], range[1], 0, &t, result);
        pnt.SetX(result[0]);
        pnt.SetY(result[1]);
        pnt.SetZ(result[2]);
      }

    result[0] = pnt.X();
    result[1] = pnt.Y();
    result[2] = pnt.Z();
    *param    = t;

  } else {

    // 2D -- surfaces & Faces
    egadsSurface *psurf = (egadsSurface *) ref->blind;
    Handle(Geom_Surface) hSurface = psurf->handle;
    if (geom->oclass != SURFACE) {
      srange[0] = psurf->urange[0];
      srange[1] = psurf->urange[1];
      srange[2] = psurf->vrange[0];
      srange[3] = psurf->vrange[1];
    } else {
      srange[0] = range[0];
      srange[1] = range[1];
      srange[2] = range[2];
      srange[3] = range[3];
    }

    GeomAPI_ProjectPointOnSurf projPnt(pnt, hSurface);
    if (!projPnt.IsDone()) {
      stat = EG_nearestSurface(hSurface, srange, xyz, 0, param, result);
      if (stat == EGADS_DEGEN) {
        if (outLevel > 0)
          printf(" EGADS Warning: Surf Proj Incomplete - DEGEN (EG_invEvaluate)!\n");
        return stat;
      } else if (stat == EGADS_EMPTY) {
        if (outLevel > 1)
          printf(" EGADS Warning: Surf Proj Incomplete (EG_invEvaluate)!\n");
      }
      u = param[0];
      v = param[1];
      pnt.SetX(result[0]);
      pnt.SetY(result[1]);
      pnt.SetZ(result[2]);
    } else {
      if (projPnt.NbPoints() == 0) {
        if (outLevel > 0)
          printf(" EGADS Warning: No projection on Surf (EG_invEvaluate)!\n");
        return EGADS_NOTFOUND;
      }
      pnt = projPnt.NearestPoint();
      projPnt.LowerDistanceParameters(u, v);

      // not only is this expensive -- sometimes it returns inconsistent stuff!
      gp_Pnt pnt1;
      hSurface->D0(u,v, pnt1);
      if (pnt1.Distance(pnt) > 1.e-7) {
/*      printf("  *** Point distance = %le  %d %d ***!\n", pnt1.Distance(pnt),
               hSurface->IsUPeriodic(), hSurface->IsVPeriodic());  */
        stat = EG_nearestSurface(hSurface, srange, xyz, 0, param, result);
        if (stat == EGADS_DEGEN) {
          if (outLevel > 0)
            printf(" EGADS Warning: Surf Proj Incomplete - DEGEN (EG_invEvaluate)!\n");
          return stat;
        } else if (stat == EGADS_EMPTY) {
          if (outLevel > 1)
            printf(" EGADS Warning: Surf Proj Incomplete (EG_invEvaluate)!\n");
        }
        u = param[0];
        v = param[1];
        pnt.SetX(result[0]);
        pnt.SetY(result[1]);
        pnt.SetZ(result[2]);
      }
    }

    if (hSurface->IsUPeriodic()) {
      period = hSurface->UPeriod();
      if ((u+PARAMACC < srange[0]) || (u-PARAMACC > srange[1])) {
        if (period != 0.0)
          if (u+PARAMACC < srange[0]) {
            if (u+period-PARAMACC < srange[1]) u += period;
          } else {
            if (u-period+PARAMACC > srange[0]) u -= period;
          }
      }
    }
    if (hSurface->IsVPeriodic()) {
      period = hSurface->VPeriod();
      if ((v+PARAMACC < srange[2]) || (v-PARAMACC > srange[3])) {
        if (period != 0.0)
          if (v+PARAMACC < srange[2]) {
            if (v+period-PARAMACC < srange[3]) v += period;
          } else {
            if (v-period+PARAMACC > srange[2]) v -= period;
          }
      }
    }

    result[0] = pnt.X();
    result[1] = pnt.Y();
    result[2] = pnt.Z();
    param[0]  = u;
    param[1]  = v;

    // clip to Face bounds
    stat = EGADS_SUCCESS;
    if (geom->oclass == FACE) stat = EG_invEvalClip(geom, xyz, param, result);
    if (stat != EGADS_SUCCESS) return stat;

  }

  return EGADS_SUCCESS;
}


int
EG_invEvaluateGuess(const egObject *geom, double *xyz,
                    double *param, double *result)
{
  int            stat, per;
  double         range[4], grange[4];
  const egObject *ref;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
    return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;

  if ((geom->oclass == PCURVE) || (geom->oclass == CURVE) ||
      (geom->oclass == SURFACE)) {
    if (geom->mtype == BSPLINE) {
      stat = EG_getRange(geom, range, &per);
      if (stat != EGADS_SUCCESS) return stat;
      if (per  != 0) {
        return EG_invEvaluate(geom, xyz, param, result);
      }
    }
    return EG_invEvaluateGeomGuess(geom, NULL, xyz, param, result);
  }

  if (geom->oclass == EDGE) {

    if (geom->mtype == DEGENERATE) return EGADS_DEGEN;
    stat = EG_getRange(geom, range, &per);
    if (stat != EGADS_SUCCESS) return stat;
    egadsEdge *pedge = (egadsEdge *) geom->blind;
    ref   = pedge->curve;
    if (ref == NULL)        return EGADS_NULLOBJ;
    if (ref->blind == NULL) return EGADS_NODATA;
    if (ref->mtype == BSPLINE) {
      stat = EG_getRange(ref, grange, &per);
      if (stat != EGADS_SUCCESS) return stat;
      if (per  != 0) {
        return EG_invEvaluate(geom, xyz, param, result);
      }
    }
    stat = EG_invEvaluateGeomGuess(ref, range, xyz, param, result);

  } else {

    egadsFace *pface = (egadsFace *) geom->blind;
    ref              = pface->surface;
    if (ref == NULL)        return EGADS_NULLOBJ;
    if (ref->blind == NULL) return EGADS_NODATA;
    if (ref->mtype == BSPLINE) {
      stat = EG_getRange(ref, grange, &per);
      if (stat != EGADS_SUCCESS) return stat;
      if (per  != 0) {
        return EG_invEvaluate(geom, xyz, param, result);
      }
    }
    stat = EG_invEvaluateGeomGuess(ref, NULL, xyz, param, result);

  }

  return stat;
}


int
EG_arcLength(const egObject *geom, double t1, double t2, double *alen)
{

  *alen = 0.0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) && (geom->oclass != CURVE) &&
      (geom->oclass != EDGE))      return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;

  if (geom->oclass == PCURVE) {

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    if (ppcurv == NULL) return EGADS_NULLOBJ;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    Geom2dAdaptor_Curve AC(hCurve);
    *alen = GCPnts_AbscissaPoint::Length(AC, t1, t2);

  } else if (geom->oclass == CURVE) {

    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    if (pcurve == NULL) return EGADS_NULLOBJ;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    GeomAdaptor_Curve AC(hCurve);
    *alen = GCPnts_AbscissaPoint::Length(AC, t1, t2);

  } else {

    if (geom->mtype == DEGENERATE) return EGADS_SUCCESS;
    egadsEdge *pedge  = (egadsEdge *) geom->blind;
    if (pedge  == NULL) return EGADS_NULLOBJ;
    egObject  *curvo  = pedge->curve;
    if (curvo  == NULL) return EGADS_NULLOBJ;
    egadsCurve *pcurve = (egadsCurve *) curvo->blind;
    if (pcurve == NULL) return EGADS_NULLOBJ;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    GeomAdaptor_Curve AC(hCurve);
    *alen = GCPnts_AbscissaPoint::Length(AC, t1, t2);

  }

  return EGADS_SUCCESS;
}


int
EG_approximate(egObject *context, int maxdeg, double tol, const int *sizes,
               const double *data, egObject **bspline)
{
  int      i, j, outLevel, stat, fixed, imax, len = 0;
  egObject *obj, *ref;

  *bspline = NULL;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (EG_sameThread(context))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(context);
  fixed    = EG_fixedKnots(context);

  if ((maxdeg < 0) || (maxdeg > 8)) {
    if (outLevel > 0)
      printf(" EGADS Warning: maxDeg = %d (EG_approximate)!\n", maxdeg);
    return EGADS_RANGERR;
  }

  if (sizes[1] == -1)
    if ((maxdeg < 3) && (sizes[0] > 2))
      return EG_spline1d(context, maxdeg, -sizes[0], data, tol, bspline);

  if (sizes[1] == 0) {

    if ((maxdeg < 3) && (sizes[0] > 2)) {
      imax = sizes[0];
      if (fixed != 0) imax = -imax;
      return EG_spline1d(context, maxdeg, imax, data, tol, bspline);
    }

    // curve
    Handle(Geom_BSplineCurve) hCurve;
    if (sizes[0] == 2) {
      TColStd_Array1OfReal    aKnots(1, sizes[0]);
      TColStd_Array1OfInteger aMults(1, sizes[0]);
      TColgp_Array1OfPnt      aPoles(1, sizes[0]);
      aKnots(1) = 0.0;
      aMults(1) = 2;
      aKnots(2) = 1.0;
      aMults(2) = 2;
      for (i = 0; i < sizes[0]; i++)
        aPoles(i+1) = gp_Pnt(data[3*i], data[3*i+1], data[3*i+2]);
      hCurve = new Geom_BSplineCurve(aPoles, aKnots, aMults, 1);
      if (hCurve.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Error: Null Curve (EG_approximate)!\n");
        return EGADS_GEOMERR;
      }
    } else {
      try {
        TColgp_Array1OfPnt aPnts(1, sizes[0]);
        for (int i = 1; i <= sizes[0]; i++, len+=3)
          aPnts(i) = gp_Pnt(data[len], data[len+1], data[len+2]);
        hCurve = GeomAPI_PointsToBSpline(aPnts, 3, maxdeg, GeomAbs_C2,
                                         tol).Curve();
      }
      catch (const Standard_Failure& e) {
        if (outLevel > 0) {
          printf(" EGADS Warning: Internal Error (EG_approximate)!\n");
          printf("                %s\n", e.GetMessageString());
        }
        return EGADS_GEOMERR;
      }
      catch (...) {
        if (outLevel > 0)
          printf(" EGADS Warning: Internal Error (EG_approximate)!\n");
        return EGADS_GEOMERR;
      }
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: make Curve = %d (EG_approximate)!\n",
               stat);
      return stat;
    }
    obj->oclass        = CURVE;
    obj->mtype         = BSPLINE;
    egadsCurve *pcurve = new egadsCurve;
    pcurve->handle     = hCurve;
    pcurve->ref        = NULL;
    pcurve->topFlg     = 0;
    pcurve->header     = NULL;
    pcurve->data       = NULL;
    pcurve->data_dot   = NULL;
    obj->blind         = pcurve;
    EG_getGeometry(obj, &i, &j, &ref, &pcurve->header, &pcurve->data);
    EG_getGeometryLen(obj, &i, &pcurve->dataLen);
    pcurve->trange[0]  = hCurve->FirstParameter();
    pcurve->trange[1]  = hCurve->LastParameter();
    EG_referenceObject(obj, context);

  } else if ((sizes[0] <= 2) || (sizes[1] < 1)) {

    if (outLevel > 0)
      printf(" EGADS Error: Sizes = %d %d (EG_approximate)!\n",
             sizes[0], sizes[1]);
    return EGADS_RANGERR;

  } else {

    if (maxdeg < 3)
      return EG_spline2d(context, maxdeg, NULL, sizes[0], sizes[1], data, tol,
                         bspline);

    // surface
    Handle(Geom_BSplineSurface) hSurf;
    try {
      TColgp_Array2OfPnt aPnts(1, sizes[0], 1, sizes[1]);
      for (int j = 1; j <= sizes[1]; j++)
        for (int i = 1; i <= sizes[0]; i++, len+=3)
          aPnts(i, j) = gp_Pnt(data[len], data[len+1], data[len+2]);
      if (tol != 0.0) {
        hSurf = GeomAPI_PointsToBSplineSurface(aPnts, 3, maxdeg, GeomAbs_C2,
                                               tol).Surface();
      } else {
        GeomAPI_PointsToBSplineSurface P2BSpl;
        P2BSpl.Interpolate(aPnts);
        hSurf = P2BSpl.Surface();
      }
    }
    catch (const Standard_Failure& e) {
      if (outLevel > 0) {
        printf(" EGADS Warning: Internal Error (EG_approximate)!\n");
        printf("                %s\n", e.GetMessageString());
      }
      return EGADS_GEOMERR;
    }
    catch (...) {
      if (outLevel > 0)
        printf(" EGADS Warning: Internal Error (EG_approximate)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: make Surface = %d (EG_approximate)!\n",
               stat);
      return stat;
    }
    obj->oclass         = SURFACE;
    obj->mtype          = BSPLINE;
    egadsSurface *psurf = new egadsSurface;
    psurf->handle       = hSurf;
    psurf->ref          = NULL;
    psurf->topFlg       = 0;
    psurf->header       = NULL;
    psurf->data         = NULL;
    psurf->data_dot     = NULL;
    obj->blind          = psurf;
    EG_getGeometry(obj, &i, &j, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(obj, &i, &psurf->dataLen);
    hSurf->Bounds(psurf->urange[0], psurf->urange[1],
                  psurf->vrange[0], psurf->vrange[1]);
    EG_referenceObject(obj, context);

  }

  *bspline = obj;
  return EGADS_SUCCESS;
}


int
EG_approximate_dot(egObject *bspline, int maxdeg, double tol, const int *sizes,
                   const double *data, const double *data_dot)
{
  int      outLevel, stat, imax, fixed, header[7];
  double   *rdata=NULL, *rdata_dot=NULL;
  egObject *context;

  if (bspline == NULL)               return EGADS_NULLOBJ;
  if (bspline->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (bspline->mtype != BSPLINE)     return EGADS_NOTGEOM;
  context  = EG_context(bspline);
  if (context == NULL)               return EGADS_NOTCNTX;
  if (EG_sameThread(context))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(context);
  fixed    = EG_fixedKnots(context);

  if ((maxdeg < 0) || (maxdeg > 8)) {
    if (outLevel > 0)
      printf(" EGADS Warning: maxDeg = %d (EG_approximate)!\n", maxdeg);
    return EGADS_RANGERR;
  }

  if (sizes[1] == -1)
    if ((maxdeg < 3) && (sizes[0] > 2)) {

      stat = EG_spline1dFit_dot(maxdeg, -sizes[0], data, data_dot,
                                NULL, NULL, tol, header,
                                &rdata, &rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_setGeometry_dot(bspline, CURVE, BSPLINE, header, rdata, rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      goto cleanup;
    }

  if (sizes[1] == 0) {

    if ((maxdeg < 3) && (sizes[0] > 2)) {

      imax = sizes[0];
      if (fixed != 0) imax = -imax;
      stat = EG_spline1dFit_dot(maxdeg, imax, data, data_dot,
                                NULL, NULL, tol, header,
                                &rdata, &rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_setGeometry_dot(bspline, CURVE, BSPLINE, header, rdata, rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

    } else {
      if (outLevel > 0)
        printf(" EGADS Error: maxDeg = %d > 3 does not have sensitivities (EG_approximate)!\n",
               maxdeg);
      stat = EGADS_RANGERR;
      goto cleanup;
    }

  } else if ((sizes[0] <= 2) || (sizes[1] < 1)) {

    if (outLevel > 0)
      printf(" EGADS Error: Sizes = %d %d (EG_approximate)!\n",
             sizes[0], sizes[1]);
    stat = EGADS_RANGERR;
    goto cleanup;

  } else {

    if (maxdeg < 3) {

      SurrealS<1> *xyz = new SurrealS<1>[3*sizes[0]*sizes[1]];
      SurrealS<1> *sdata = NULL;

      for (int i = 0; i < 3*sizes[0]*sizes[1]; i++) {
        xyz[i].value() = data[i];
        xyz[i].deriv() = data_dot[i];
      }

      stat = EG_spline2dAppr< SurrealS<1> >(maxdeg, sizes[0], sizes[1], xyz,
                                            NULL, NULL, NULL, NULL, NULL,
                                            NULL, NULL, NULL, NULL,
                                            tol, header, &sdata);
      delete [] xyz;
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_setGeometry_dot(bspline, SURFACE, BSPLINE, header, sdata);
      EG_free(sdata);
      if (stat != EGADS_SUCCESS) goto cleanup;

    } else {
      if (outLevel > 0)
        printf(" EGADS Error: maxDeg = %d > 3 does not have sensitivities (EG_approximate)!\n",
               maxdeg);
      stat = EGADS_RANGERR;
      goto cleanup;
    }

  }

  stat = EGADS_SUCCESS;

cleanup:
  EG_free(rdata);
  EG_free(rdata_dot);

  return stat;
}


int
EG_otherCurve(const egObject *surface, const egObject *curve,
              double tol, egObject **newcurve)
{
  int      outLevel, stat;
  egObject *context, *obj;

  *newcurve = NULL;
  if  (surface == NULL)               return EGADS_NULLOBJ;
  if  (surface->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((surface->oclass != SURFACE) &&
      (surface->oclass != FACE))      return EGADS_NOTGEOM;
  if  (surface->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(surface))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(surface);
  context  = EG_context(surface);

  if (curve == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Input Curve (EG_otherCurve)!\n");
    return EGADS_NULLOBJ;
  }
  if ((curve->oclass != PCURVE) && (curve->oclass != CURVE) &&
      (curve->oclass != EDGE)) {
    if (outLevel > 0)
      printf(" EGADS Error: Not a PCurve/Curve or Edge (EG_otherCurve)!\n");
    return EGADS_NOTGEOM;
  }
  if (curve->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: PCurve has no data (EG_otherCurve)!\n");
    return EGADS_NODATA;
  }
  if (EG_context(curve) != context) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_otherCurve)!\n");
    return EGADS_MIXCNTX;
  }

  Handle(Geom_Surface) hSurface;
  Standard_Real        prec = tol;
  if (surface->oclass == SURFACE) {
    egadsSurface *psurf = (egadsSurface *) surface->blind;
    hSurface = psurf->handle;
  } else {
    egadsFace *pface = (egadsFace *) surface->blind;
    if (pface->surface == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face has no surface (EG_otherCurve)!\n");
      return EGADS_NODATA;
    }
    obj = pface->surface;
    if (obj->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface has no data (EG_otherCurve)!\n");
      return EGADS_NODATA;
    }
    egadsSurface *psurf = (egadsSurface *) obj->blind;
    hSurface = psurf->handle;
    double toler = BRep_Tool::Tolerance(pface->face);
    if (prec < toler) prec = toler;
  }
  if (prec < Precision::Confusion()) prec = Precision::Confusion();

  if (curve->oclass == PCURVE) {

    Standard_Real maxDev, aveDev;

    egadsPCurve *ppcurv         = (egadsPCurve *) curve->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    GeomAdaptor_Surface  aGAS   = hSurface;
    Handle(GeomAdaptor_HSurface) aHGAS = new GeomAdaptor_HSurface(aGAS);
    Handle(Geom2dAdaptor_HCurve) Crv   = new Geom2dAdaptor_HCurve(hCurve);
    Adaptor3d_CurveOnSurface ConS(Crv,aHGAS);

    Handle(Geom_Curve) newcrv;
    GeomLib::BuildCurve3d(prec, ConS, hCurve->FirstParameter(),
                          hCurve->LastParameter(), newcrv, maxDev, aveDev);

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_otherCurve)!\n", stat);
      return stat;
    }
    EG_completeCurve(obj, newcrv);

  } else {

    Handle(Geom2d_Curve) newcrv;

    if (curve->oclass == EDGE) {

      Standard_Real t1, t2;

      egadsEdge *pedge = (egadsEdge *) curve->blind;
      // is the PCurve already attached?
      if (surface->oclass == FACE) {
        egadsFace *pface = (egadsFace *) surface->blind;
        newcrv = BRep_Tool::CurveOnSurface(pedge->edge, pface->face, t1, t2);
      }
      if (newcrv.IsNull()) {
        double toler     = BRep_Tool::Tolerance(pedge->edge);
        if (prec < toler) prec = toler;
        egObject *geom   = pedge->curve;
        if (geom->blind == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: NULL Curve Data (EG_otherCurve)!\n");
          return EGADS_NODATA;
        }
        Handle(Geom_Curve) hCurve = BRep_Tool::Curve(pedge->edge, t1, t2);
        try {
          newcrv = GeomProjLib::Curve2d(hCurve, t1, t2, hSurface, prec);
        }
        catch (const Standard_Failure& e) {
          printf(" EGADS Warning: Geometry Creation Error (EG_otherCurve)!\n");
          printf("                %s\n", e.GetMessageString());
          return EGADS_GEOMERR;
        }
        catch (...) {
          printf(" EGADS Warning: Geometry Creation Error (EG_otherCurve)!\n");
          return EGADS_GEOMERR;
        }
      }

    } else {

      egadsCurve *pcurve        = (egadsCurve *) curve->blind;
      Handle(Geom_Curve) hCurve = pcurve->handle;
      try {
        newcrv = GeomProjLib::Curve2d(hCurve, hCurve->FirstParameter(),
                                      hCurve->LastParameter(), hSurface, prec);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Geometry Creation Error (EG_otherCurve)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Warning: Geometry Creation Error (EG_otherCurve)!\n");
        return EGADS_GEOMERR;
      }

    }

    if (EG_getPCurveType(newcrv) == 0) {
      if (outLevel > 0)
        printf(" EGADS Info: Cannot construct PCurve (EG_otherCurve)!\n");
      return EGADS_CONSTERR;
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make PCurve = %d (EG_otherCurve)!\n", stat);
      return stat;
    }
    EG_completePCurve(obj, newcrv);
  }

  EG_referenceObject(obj, context);
  *newcurve = obj;
  return EGADS_SUCCESS;
}


namespace // no-name namespace makes private to this file
{

/* Function from GeomFill_BSplineCurves.cxx used internally of GeomFill_BSplineCurves */
Standard_Integer EG_SetSameDistribution(Handle(Geom_BSplineCurve)& C1,
                                        Handle(Geom_BSplineCurve)& C2)
{
  Standard_Integer nbp1 = C1->NbPoles();
  Standard_Integer nbk1 = C1->NbKnots();
  TColgp_Array1OfPnt      P1(1,nbp1);
  TColStd_Array1OfReal    W1(1,nbp1);
  W1.Init(1.);
  TColStd_Array1OfReal    K1(1,nbk1);
  TColStd_Array1OfInteger M1(1,nbk1);

  C1->Poles(P1);
  if (C1->IsRational()) C1->Weights(W1);
  C1->Knots(K1);
  C1->Multiplicities(M1);

  Standard_Integer nbp2 = C2->NbPoles();
  Standard_Integer nbk2 = C2->NbKnots();
  TColgp_Array1OfPnt      P2(1,nbp2);
  TColStd_Array1OfReal    W2(1,nbp2);
  W2.Init(1.);
  TColStd_Array1OfReal    K2(1,nbk2);
  TColStd_Array1OfInteger M2(1,nbk2);

  C2->Poles(P2);
  if (C2->IsRational()) C2->Weights(W2);
  C2->Knots(K2);
  C2->Multiplicities(M2);

  Standard_Real K11 = K1( 1  );
  Standard_Real K12 = K1(nbk1);
  Standard_Real K21 = K2( 1  );
  Standard_Real K22 = K2(nbk2);

  if ((K12-K11) > (K22-K21)) {
    BSplCLib::Reparametrize( K11, K12, K2);
    C2->SetKnots(K2);
  } else if ((K12-K11) < (K22-K21)) {
    BSplCLib::Reparametrize( K21, K22, K1);
    C1->SetKnots(K1);
  } else if(Abs(K12-K11) > Precision::PConfusion()) {
    BSplCLib::Reparametrize( K11, K12, K2);
    C2->SetKnots(K2);
  }

  Standard_Integer NP,NK;
  if (BSplCLib::PrepareInsertKnots(C1->Degree(),Standard_False,
                                   K1,M1,K2,&M2,NP,NK,Precision::PConfusion(),
                                   Standard_False)) {
    TColgp_Array1OfPnt      NewP(1, NP);
    TColStd_Array1OfReal    NewW(1, NP);
    TColStd_Array1OfReal    NewK(1, NK);
    TColStd_Array1OfInteger NewM(1, NK);
    BSplCLib::InsertKnots(C1->Degree(),Standard_False,
                          P1,&W1,K1,M1,K2,&M2,
                          NewP,&NewW,NewK,NewM,Precision::PConfusion(),
                          Standard_False);
    if (C1->IsRational()) {
      C1 = new Geom_BSplineCurve(NewP,NewW,NewK,NewM,C1->Degree());
    } else {
      C1 = new Geom_BSplineCurve(NewP,NewK,NewM,C1->Degree());
    }
    BSplCLib::InsertKnots(C2->Degree(),Standard_False,
                          P2,&W2,K2,M2,K1,&M1,
                          NewP,&NewW,NewK,NewM,Precision::PConfusion(),
                          Standard_False);
    if (C2->IsRational()) {
      C2 = new Geom_BSplineCurve(NewP,NewW,NewK,NewM,C2->Degree());
    } else {
      C2 = new Geom_BSplineCurve(NewP,NewK,NewM,C2->Degree());
    }
  } else {
    throw Standard_ConstructionError(" ");
  }

  return C1->NbPoles();
}


/* Function taken from GeomFill_BSplineCurves.cxx
 *
 * The OCC function attempts to arrange the curves based on a tolerance at the end points
 * This version assumes the curves are properly arranged as
 *
 *              CC3
 *          ----->-----
 *         |           |
 *         |           |
 *         |           |
 *     CC4 ^           ^ CC2
 *         |           |
 *         |           |
 *          ----->-----
 *           CC1 = C1
 */
Handle(Geom_BSplineSurface)
EG_GeomFill_BSplineCurves(Handle(Geom_BSplineCurve)& CC1,
                          Handle(Geom_BSplineCurve)& CC2,
                          Handle(Geom_BSplineCurve)& CC3,
                          Handle(Geom_BSplineCurve)& CC4,
                          const GeomFill_FillingStyle Type)
{
  Standard_Integer Deg1 = CC1->Degree();
  Standard_Integer Deg2 = CC2->Degree();
  Standard_Integer Deg3 = CC3->Degree();
  Standard_Integer Deg4 = CC4->Degree();
  Standard_Integer DegU = Max( Deg1, Deg3);
  Standard_Integer DegV = Max( Deg2, Deg4);
  if (Deg1 < DegU) CC1->IncreaseDegree(DegU);
  if (Deg2 < DegV) CC2->IncreaseDegree(DegV);
  if (Deg3 < DegU) CC3->IncreaseDegree(DegU);
  if (Deg4 < DegV) CC4->IncreaseDegree(DegV);

  Standard_Integer NbUPoles = EG_SetSameDistribution(CC1,CC3);
  Standard_Integer NbVPoles = EG_SetSameDistribution(CC2,CC4);

  if (Type == GeomFill_CoonsStyle) {
    if (NbUPoles < 4 || NbVPoles < 4)
      throw Standard_ConstructionError("GeomFill_BSplineCurves: invalid filling style");
  }

  TColgp_Array1OfPnt P1(1,NbUPoles);
  TColgp_Array1OfPnt P2(1,NbVPoles);
  TColgp_Array1OfPnt P3(1,NbUPoles);
  TColgp_Array1OfPnt P4(1,NbVPoles);
  CC1->Poles(P1);
  CC2->Poles(P2);
  CC3->Poles(P3);
  CC4->Poles(P4);

  Standard_Boolean isRat = ( CC1->IsRational() || CC2->IsRational() ||
                             CC3->IsRational() || CC4->IsRational() );

  TColStd_Array1OfReal W1(1,NbUPoles);
  TColStd_Array1OfReal W3(1,NbUPoles);
  TColStd_Array1OfReal W2(1,NbVPoles);
  TColStd_Array1OfReal W4(1,NbVPoles);
  W1.Init(1.);
  W2.Init(1.);
  W3.Init(1.);
  W4.Init(1.);
  if (isRat) {
    if (CC1->IsRational()) {
      CC1->Weights(W1);
    }
    if (CC2->IsRational()) {
      CC2->Weights(W2);
    }
    if (CC3->IsRational()) {
      CC3->Weights(W3);
    }
    if (CC4->IsRational()) {
      CC4->Weights(W4);
    }
  }

  GeomFill_Filling Caro;
  if (isRat) {
    switch (Type)
    {
    case GeomFill_StretchStyle:
      Caro = GeomFill_Stretch(P1, P2, P3, P4, W1, W2, W3, W4);
      break;
    case GeomFill_CoonsStyle:
      Caro = GeomFill_Coons(P1, P4, P3, P2, W1, W4, W3, W2);
      break;
    case GeomFill_CurvedStyle:
      Caro = GeomFill_Curved(P1, P2, P3, P4, W1, W2, W3, W4);
      break;
    }
  }
  else {
    switch (Type)
    {
    case GeomFill_StretchStyle:
      Caro = GeomFill_Stretch(P1, P2, P3, P4);
      break;
    case GeomFill_CoonsStyle:
      Caro = GeomFill_Coons(P1, P4, P3, P2);
      break;
    case GeomFill_CurvedStyle:
      Caro = GeomFill_Curved(P1, P2, P3, P4);
      break;
    }
  }

  NbUPoles = Caro.NbUPoles();
  NbVPoles = Caro.NbVPoles();
  TColgp_Array2OfPnt Poles(1,NbUPoles,1,NbVPoles);


  Standard_Integer NbUKnot = CC1->NbKnots();
  TColStd_Array1OfReal    UKnots(1,NbUKnot);
  TColStd_Array1OfInteger UMults(1,NbUKnot);
  CC1->Knots(UKnots);
  CC1->Multiplicities(UMults);

  Standard_Integer NbVKnot = CC2->NbKnots();
  TColStd_Array1OfReal    VKnots(1,NbVKnot);
  TColStd_Array1OfInteger VMults(1,NbVKnot);
  CC2->Knots(VKnots);
  CC2->Multiplicities(VMults);

  Caro.Poles(Poles);

  if (Caro.isRational()) {
    TColStd_Array2OfReal Weights(1,NbUPoles, 1,NbVPoles);
    Caro.Weights(Weights);
    return new Geom_BSplineSurface(Poles        , Weights,
                                   UKnots       , VKnots,
                                   UMults       , VMults,
                                   CC1->Degree(), CC2->Degree());
  } else {
    return new Geom_BSplineSurface(Poles        ,
                                   UKnots       , VKnots,
                                   UMults       , VMults,
                                   CC1->Degree(), CC2->Degree());
  }

}


Handle(Geom_BSplineSurface)
EG_GeomFill_BSplineCurves(Handle(Geom_BSplineCurve)& C1,
                          Handle(Geom_BSplineCurve)& C2,
                          Handle(Geom_BSplineCurve)& C3,
                          const GeomFill_FillingStyle Type)
{
  Handle(Geom_BSplineCurve) C4;
  TColgp_Array1OfPnt      Poles(1,2);
  TColStd_Array1OfReal    Knots(1,2);
  TColStd_Array1OfInteger Mults(1,2);

  Poles( 1) = C1->StartPoint();
  Poles( 2) = C3->StartPoint();
  Knots( 1) = C2->Knot(C2->FirstUKnotIndex());
  Knots( 2) = C2->Knot(C2->LastUKnotIndex());
  Mults( 1) = Mults( 2) = 2;
  C4 = new Geom_BSplineCurve(Poles, Knots, Mults, 1);
  return EG_GeomFill_BSplineCurves(C1, C2, C3, C4, Type);
}

} // namespace


int
EG_isoCline(const egObject *surface, int UV, double value,
                  egObject **newcurve)
{
  int      i, j, stat, outLevel;
  egObject *context, *obj, *ref;

  *newcurve = NULL;
  if  (surface == NULL)               return EGADS_NULLOBJ;
  if  (surface->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((surface->oclass != SURFACE) &&
      (surface->oclass != LOOP))      return EGADS_NOTGEOM;
  if  (surface->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(surface))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(surface);
  context  = EG_context(surface);

  // special loop fitting code
  if (surface->oclass == LOOP) {
    if (surface->mtype != CLOSED) {
      if (outLevel > 1)
        printf(" EGADS Info: Input Loop is NOT closed (EG_isoCline)!\n");
      return EGADS_TOPOERR;
    }
    egadsLoop *ploop  = (egadsLoop *) surface->blind;
    if (ploop->surface != NULL) {
      if (outLevel > 1)
        printf(" EGADS Info: Input Loop has attached Surface (EG_isoCline)!\n");
      return EGADS_TOPOERR;
    }
    if ((ploop->nedges < 3) || (ploop->nedges > 4)) {
      if (outLevel > 0)
        printf(" EGADS Error: Loop has %d Edges (EG_isoCline)!\n",
               ploop->nedges);
      return EGADS_GEOMERR;
    }

    Handle(Geom_BSplineCurve) bsc[4];
    ShapeConstruct_Curve      ShapeCC;
    gp_Trsf                   form  = gp_Trsf();
    Standard_Real             prec  = Precision::Confusion();
    if (value > prec)         prec  = value;
    GeomFill_FillingStyle     style = GeomFill_StretchStyle;
    if (UV < 0)               style = GeomFill_CurvedStyle;
    if (UV > 0)               style = GeomFill_CoonsStyle;
    for (i = 0; i < ploop->nedges; i++) {
      obj = ploop->edges[i];
      egadsEdge *pedge = (egadsEdge *) obj->blind;
      if (pedge == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: %d/%d NULL pedge (EG_isoCline)!\n",
                 i+1, ploop->nedges);
        return EGADS_NODATA;
      }
      double range[2];
      BRep_Tool::Range(pedge->edge, range[0], range[1]);
      obj = pedge->curve;
      egadsCurve *pcurve = (egadsCurve *) obj->blind;
      if (pcurve == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: %d/%d NULL pcurve (EG_isoCline)!\n",
                 i+1, ploop->nedges);
        return EGADS_NODATA;
      }
      Handle(Geom_Curve)    hCurve = pcurve->handle;
      // copy curve so we don't change it!
      Handle(Geom_Geometry) nGeom  = hCurve->Transformed(form);
      Handle(Geom_Curve)    nCurve = Handle(Geom_Curve)::DownCast(nGeom);
      bsc[i] = ShapeCC.ConvertToBSpline(nCurve, range[0], range[1], prec);
      if (bsc[i].IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Warning: Failure to Convert %d/%d (EG_isoCline)!\n",
                 i+1, ploop->nedges);
        return EGADS_GEOMERR;
      }

      //    The Curves are arranged in this way:
      //
      //                      bsc[2]
      //                  ----->-----
      //                 |           |
      //                 |           |
      //                 |           |
      //          bsc[3] ^           ^ bsc[1]
      //                 |           |
      //                 |           |
      //                  ----->-----
      //                     bsc[0]

      if (i < 2) {
        if (ploop->senses[i] == SREVERSE)
          bsc[i] = Handle(Geom_BSplineCurve)::DownCast(bsc[i]->Reversed());
      }
      else {
        if (ploop->senses[i] == SFORWARD)
          bsc[i] = Handle(Geom_BSplineCurve)::DownCast(bsc[i]->Reversed());
      }
    }

    Handle(Geom_BSplineSurface) hSurface;
    try {
      if (ploop->nedges == 3) {
        hSurface = EG_GeomFill_BSplineCurves(bsc[0], bsc[1], bsc[2], style);
      } else {
        hSurface = EG_GeomFill_BSplineCurves(bsc[0], bsc[1], bsc[2], bsc[3],
                                             style);
      }
    }
    catch (const Standard_Failure& e) {
      if (outLevel > 0) {
        printf(" EGADS Warning: Geometry Creation Error (EG_isoCline)!\n");
        printf("                %s\n", e.GetMessageString());
      }
      return EGADS_GEOMERR;
    }
    catch (...) {
      if (outLevel > 0)
        printf(" EGADS Warning: Geometry Creation Error (EG_isoCline)!\n");
      return EGADS_GEOMERR;
    }
    if (hSurface.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Warning: Failure to Construct (EG_isoCline)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Surface = %d (EG_isoCline)!\n", stat);
      return stat;
    }
    obj->oclass         = SURFACE;
    obj->mtype          = BSPLINE;
    egadsSurface *psurf = new egadsSurface;
    psurf->handle       = hSurface;
    psurf->ref          = NULL;
    psurf->topFlg       = 0;
    psurf->header       = NULL;
    psurf->data         = NULL;
    psurf->data_dot     = NULL;
    obj->blind          = psurf;
    EG_getGeometry(obj, &i, &j, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(obj, &i, &psurf->dataLen);
    hSurface->Bounds(psurf->urange[0], psurf->urange[1],
                     psurf->vrange[0], psurf->vrange[1]);
    EG_referenceObject(obj, context);
    *newcurve = obj;

    return EGADS_SUCCESS;
  }

  // normal isocline code
  egadsSurface *psurf           = (egadsSurface *) surface->blind;
  Handle(Geom_Surface) hSurface = psurf->handle;
  Handle(Geom_Curve)   newcrv;
  if (UV == UISO) {
    newcrv = hSurface->UIso(value);
  } else {
    newcrv = hSurface->VIso(value);
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 1)
      printf(" EGADS Error: make Curve = %d (EG_isoCline)!\n", stat);
    return stat;
  }
  EG_completeCurve(obj, newcrv);
  EG_referenceObject(obj, context);
  *newcurve = obj;

  return EGADS_SUCCESS;
}


int
EG_convertToBSplineRange(egObject *object, const double *range,
                         egObject **bspline)
{
  int      i, j, n, m, outLevel, stat, header[4];
  double   data[8], d, x0[2], x1[2];
  gp_Pnt2d pnt;
  egObject *obj, *geom, *context, *ref;

  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != PCURVE) &&
      (object->oclass != CURVE)  && (object->oclass != SURFACE) &&
      (object->oclass != EDGE)   && (object->oclass != FACE))
    return EGADS_NOTGEOM;
  if  (object->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(object))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);
  geom     = object;

  if (object->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) object->blind;
    geom = pedge->curve;
    if (geom->blind == NULL) return EGADS_NODATA;
  } else if (object->oclass == FACE) {
    egadsFace *pface = (egadsFace *) object->blind;
    geom = pface->surface;
    if (geom->blind == NULL) return EGADS_NODATA;
  }

  if (geom->oclass == PCURVE) {

    egadsPCurve *ppcurv         = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;

    if (geom->mtype == LINE) {
      header[0] = 0;
      header[1] = 1;
      header[2] = 2;
      header[3] = 4;
      data[0]   = data[1] = range[0];
      data[2]   = data[3] = range[1];
      hCurve->D0(range[0], pnt);
      data[4]   = pnt.X();
      data[5]   = pnt.Y();
      hCurve->D0(range[1], pnt);
      data[6]   = pnt.X();
      data[7]   = pnt.Y();
      return EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                             bspline);
    }

    try {
      Handle(Geom2d_BSplineCurve) hBSpline;
      if (geom->mtype == BSPLINE) {
        ShapeConstruct_Curve ShapeCC;
        hBSpline = ShapeCC.ConvertToBSpline(hCurve, range[0], range[1],
                                            Precision::Confusion());
      } else {
        hBSpline = ShapeConstruct::ConvertCurveToBSpline(hCurve, range[0],
                                                         range[1],
                                                         Precision::Confusion(),
                                                         GeomAbs_C2, 100, 20);
      }
      if (hBSpline.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Warning: Failure to convert (EG_convertToBSplineRange)!\n");
        return EGADS_GEOMERR;
      }

      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: make PCurve = %d (EG_convertToBSplineRange)!\n",
               stat);
        return stat;
      }
      obj->oclass        = PCURVE;
      obj->mtype         = BSPLINE;
      egadsPCurve *ppcrv = new egadsPCurve;
      ppcrv->handle      = hBSpline;
      ppcrv->ref         = NULL;
      ppcrv->topFlg      = 0;
      ppcrv->header      = NULL;
      ppcrv->data        = NULL;
      ppcrv->data_dot    = NULL;
      obj->blind         = ppcrv;
      EG_getGeometry(obj, &i, &j, &ref, &ppcrv->header, &ppcrv->data);
      EG_getGeometryLen(obj, &i, &ppcrv->dataLen);
      ppcrv->trange[0]   = hBSpline->FirstParameter();
      ppcrv->trange[1]   = hBSpline->LastParameter();
      for (n = 2; n < ppcrv->header[2]; n++) {
        m     = ppcrv->header[3] + 2*n - 4;
        x0[0] = ppcrv->data[m+2] - ppcrv->data[m  ];
        x0[1] = ppcrv->data[m+3] - ppcrv->data[m+1];
        d     = sqrt(x0[0]*x0[0] + x0[1]*x0[1]);
        if (d != 0.0) {
          x0[0] /= d;
          x0[1] /= d;
        }
        x1[0] = ppcrv->data[m+4] - ppcrv->data[m+2];
        x1[1] = ppcrv->data[m+5] - ppcrv->data[m+3];
        d     = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
        if (d != 0.0) {
          x1[0] /= d;
          x1[1] /= d;
        }
        d = x0[0]*x1[0] + x0[1]*x1[1];
        if (d < -0.95) {
#ifdef DEBUG
          printf(" EGADS Info: PCurve dot flip at %d/%d (%lf) -- %lf %lf!\n",
                 n-2, ppcrv->header[2]-2, d, ppcrv->data[0],
                 ppcrv->data[ppcrv->header[3]-1]);
#endif
          stat = EG_addStrAttr(obj, ".Bad", "CPrev");
          if (stat != EGADS_SUCCESS)
            printf("             EG_addStrAttr CPrev= %d\n", stat);
        }
      }
    }
    catch (const Standard_Failure& e) {
      if (outLevel > 0) {
        printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSplineRange)!\n");
        printf("                %s\n", e.GetMessageString());
      }
      return EGADS_GEOMERR;
    }
    catch (...) {
      if (outLevel > 0)
        printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSplineRange)!\n");
      return EGADS_GEOMERR;
    }

  } else if (geom->oclass == CURVE) {

    egadsCurve *pcurve        = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    Handle(Geom_BSplineCurve) hBSpline;
    if (geom->mtype == BSPLINE) {
      ShapeConstruct_Curve ShapeCC;
      hBSpline = ShapeCC.ConvertToBSpline(hCurve, range[0], range[1],
                                          Precision::Confusion());
    } else {
      hBSpline = ShapeConstruct::ConvertCurveToBSpline(hCurve, range[0],
                                                       range[1],
                                                       Precision::Confusion(),
                                                       GeomAbs_C2, 100, 20);
    }
    if (hBSpline.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Warning: Failure to convert (EG_convertToBSplineRange)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_convertToBSplineRange)!\n",
             stat);
      return stat;
    }
    obj->oclass       = CURVE;
    obj->mtype        = BSPLINE;
    egadsCurve *pcurv = new egadsCurve;
    pcurv->handle     = hBSpline;
    pcurv->ref        = NULL;
    pcurv->topFlg     = 0;
    pcurv->header     = NULL;
    pcurv->data       = NULL;
    pcurv->data_dot   = NULL;
    obj->blind        = pcurv;
    EG_getGeometry(obj, &i, &j, &ref, &pcurv->header, &pcurv->data);
    EG_getGeometryLen(obj, &i, &pcurv->dataLen);
    pcurv->trange[0]  = hBSpline->FirstParameter();
    pcurv->trange[1]  = hBSpline->LastParameter();

  } else {

    egadsSurface *psurface        = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurface = psurface->handle;
    try {
      Handle(Geom_BSplineSurface)
        hBSpline = ShapeConstruct::ConvertSurfaceToBSpline(hSurface,
                                   range[0], range[1], range[2], range[3],
                                   Precision::Confusion(), GeomAbs_C2, 100, 20);
      if (hBSpline.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Warning: Failure to Convert (EG_convertToBSplineRange)!\n");
        return EGADS_GEOMERR;
      }

      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: make Surface = %d (EG_convertToBSplineRange)!\n",
               stat);
        return stat;
      }
      obj->oclass         = SURFACE;
      obj->mtype          = BSPLINE;
      egadsSurface *psurf = new egadsSurface;
      psurf->handle       = hBSpline;
      psurf->ref          = NULL;
      psurf->topFlg       = 0;
      psurf->header       = NULL;
      psurf->data         = NULL;
      psurf->data_dot     = NULL;
      obj->blind          = psurf;
      EG_getGeometry(obj, &i, &j, &ref, &psurf->header, &psurf->data);
      EG_getGeometryLen(obj, &i, &psurf->dataLen);
      hBSpline->Bounds(psurf->urange[0], psurf->urange[1],
                       psurf->vrange[0], psurf->vrange[1]);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSplineRange)!\n");
      printf("                %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSplineRange)!\n");
      return EGADS_GEOMERR;
    }

  }

  *bspline = obj;
  EG_referenceObject(obj, context);

  return EGADS_SUCCESS;
}


int
EG_convertToBSpline(egObject *object, egObject **bspline)
{
  int           i, j, n, m, outLevel, stat;
  double        d, x0[2], x1[2];
  egObject      *obj, *geom, *context, *ref;
  Standard_Real range[4];

  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != PCURVE) &&
      (object->oclass != CURVE)  && (object->oclass != SURFACE) &&
      (object->oclass != EDGE)   && (object->oclass != FACE))
                                     return EGADS_NOTGEOM;
  if  (object->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(object))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);
  geom     = object;

  if (object->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) object->blind;
    geom = pedge->curve;
    if (geom->blind == NULL) return EGADS_NODATA;
  } else if (object->oclass == FACE) {
    egadsFace *pface = (egadsFace *) object->blind;
    geom = pface->surface;
    if (geom->blind == NULL) return EGADS_NODATA;
    if (geom->mtype == BSPLINE) {
      egadsSurface *psurface = (egadsSurface *) geom->blind;
      if ((psurface->urange[0] == pface->urange[0]) &&
          (psurface->urange[1] == pface->urange[1]) &&
          (psurface->vrange[0] == pface->vrange[0]) &&
          (psurface->vrange[1] == pface->vrange[1])) {
        *bspline = geom;
        return EGADS_SUCCESS;
      }
    }
  } else {
    if (object->mtype == BSPLINE) {
      *bspline = geom;
      return EGADS_SUCCESS;
    }
  }

  if (geom->oclass == PCURVE) {

    egadsPCurve *ppcurv         = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    range[0] = hCurve->FirstParameter();
    range[1] = hCurve->LastParameter();
    Handle(Geom2d_BSplineCurve)
      hBSpline = ShapeConstruct::ConvertCurveToBSpline(hCurve, range[0],
                                                       range[1],
                                                       Precision::Confusion(),
                                                       GeomAbs_C2, 100, 20);
    if (hBSpline.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Warning: Failure to convert (EG_convertToBSpline)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make PCurve = %d (EG_convertToBSpline)!\n",
             stat);
      return stat;
    }
    obj->oclass        = PCURVE;
    obj->mtype         = BSPLINE;
    egadsPCurve *ppcrv = new egadsPCurve;
    ppcrv->handle      = hBSpline;
    ppcrv->ref         = NULL;
    ppcrv->topFlg      = 0;
    ppcrv->header      = NULL;
    ppcrv->data        = NULL;
    ppcrv->data_dot    = NULL;
    obj->blind         = ppcrv;
    EG_getGeometry(obj, &i, &j, &ref, &ppcrv->header, &ppcrv->data);
    EG_getGeometryLen(obj, &i, &ppcrv->dataLen);
    ppcrv->trange[0]   = hBSpline->FirstParameter();
    ppcrv->trange[1]   = hBSpline->LastParameter();
    for (n = 2; n < ppcrv->header[2]; n++) {
      m     = ppcrv->header[3] + 2*n - 4;
      x0[0] = ppcrv->data[m+2] - ppcrv->data[m  ];
      x0[1] = ppcrv->data[m+3] - ppcrv->data[m+1];
      d     = sqrt(x0[0]*x0[0] + x0[1]*x0[1]);
      if (d != 0.0) {
        x0[0] /= d;
        x0[1] /= d;
      }
      x1[0] = ppcrv->data[m+4] - ppcrv->data[m+2];
      x1[1] = ppcrv->data[m+5] - ppcrv->data[m+3];
      d     = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
      if (d != 0.0) {
        x1[0] /= d;
        x1[1] /= d;
      }
      d = x0[0]*x1[0] + x0[1]*x1[1];
      if (d < -0.95) {
#ifdef DEBUG
        printf(" EGADS Info: PCurve %d dot flip at %d/%d (%lf) -- %lf %lf!\n",
               i, n-2, ppcrv->header[2]-2, d,
               ppcrv->data[0], ppcrv->data[ppcrv->header[3]-1]);
#endif
        stat = EG_addStrAttr(obj, ".Bad", "CPrev");
        if (stat != EGADS_SUCCESS)
          printf("             EG_addStrAttr CPrev= %d\n", stat);
      }
    }

  } else if (geom->oclass == CURVE) {

    egadsCurve *pcurve        = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    range[0] = hCurve->FirstParameter();
    range[1] = hCurve->LastParameter();
    if (object->oclass == EDGE) {
      egadsEdge *pedge = (egadsEdge *) object->blind;
      BRep_Tool::Range(pedge->edge, range[0], range[1]);
    }

    Handle(Geom_BSplineCurve) hBSpline;
    if (geom->mtype == BSPLINE) {
      ShapeConstruct_Curve ShapeCC;
      hBSpline = ShapeCC.ConvertToBSpline(hCurve, range[0], range[1],
                                          Precision::Confusion());
    } else {
      hBSpline = ShapeConstruct::ConvertCurveToBSpline(hCurve, range[0],
                                                       range[1],
                                                       Precision::Confusion(),
                                                       GeomAbs_C2, 100, 20);
    }
    if (hBSpline.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Warning: Failure to convert (EG_convertToBSpline)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_convertToBSpline)!\n",
             stat);
      return stat;
    }
    obj->oclass       = CURVE;
    obj->mtype        = BSPLINE;
    egadsCurve *pcurv = new egadsCurve;
    pcurv->handle     = hBSpline;
    pcurv->ref        = NULL;
    pcurv->topFlg     = 0;
    pcurv->header     = NULL;
    pcurv->data       = NULL;
    pcurv->data_dot   = NULL;
    obj->blind        = pcurv;
    EG_getGeometry(obj, &i, &j, &ref, &pcurv->header, &pcurv->data);
    EG_getGeometryLen(obj, &i, &pcurv->dataLen);
    pcurv->trange[0]  = hBSpline->FirstParameter();
    pcurv->trange[1]  = hBSpline->LastParameter();

  } else {

    egadsSurface *psurface        = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurface = psurface->handle;
    hSurface->Bounds(range[0],range[1], range[2],range[3]);
    if (object->oclass == FACE) {
      egadsFace *pface = (egadsFace *) object->blind;
      BRepTools::UVBounds(pface->face, range[0],range[1], range[2],range[3]);
    }
    try {
      Handle(Geom_BSplineSurface)
        hBSpline = ShapeConstruct::ConvertSurfaceToBSpline(hSurface,
                                   range[0], range[1], range[2], range[3],
                                   Precision::Confusion(), GeomAbs_C2, 100, 20);
      if (hBSpline.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Warning: Failure to Convert (EG_convertToBSpline)!\n");
        return EGADS_GEOMERR;
      }

      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: make Surface = %d (EG_convertToBSpline)!\n",
               stat);
        return stat;
      }
      obj->oclass         = SURFACE;
      obj->mtype          = BSPLINE;
      egadsSurface *psurf = new egadsSurface;
      psurf->handle       = hBSpline;
      psurf->ref          = NULL;
      psurf->topFlg       = 0;
      psurf->header       = NULL;
      psurf->data         = NULL;
      psurf->data_dot     = NULL;
      obj->blind          = psurf;
      EG_getGeometry(obj, &i, &j, &ref, &psurf->header, &psurf->data);
      EG_getGeometryLen(obj, &i, &psurf->dataLen);
      hBSpline->Bounds(psurf->urange[0], psurf->urange[1],
                       psurf->vrange[0], psurf->vrange[1]);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSpline)!\n");
      printf("                %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSpline)!\n");
      return EGADS_GEOMERR;
    }

  }

  *bspline = obj;
  EG_referenceObject(obj, context);

  return EGADS_SUCCESS;
}


int
EG_flattenBSpline(egObject *object, egObject **result)
{
  int      stat, ot, mc, NewNbUKnots, NewNbUPoles, NewNbVKnots, NewNbVPoles;
  egObject *obj, *context, *ref;

  *result = NULL;
  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != PCURVE) && (object->oclass != CURVE) &&
      (object->oclass != SURFACE))   return EGADS_NOTGEOM;
  if  (object->mtype != BSPLINE)     return EGADS_NOTGEOM;
  if  (object->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(object))        return EGADS_CNTXTHRD;

  context = EG_context(object);

  if (object->oclass == PCURVE) {
    egadsPCurve *pcurve = (egadsPCurve *) object->blind;
    if ((pcurve->header[0]&4) == 0) {
      printf(" EGADS Warning: Already flat (EG_flattenBSpline)!\n");
      return EGADS_GEOMERR;
    }
    Handle(Geom2d_Curve)      hCurve   = pcurve->handle;
    Handle(Geom2d_BSplineCurve) hBSpline =
                                  Handle(Geom2d_BSplineCurve)::DownCast(hCurve);

    TColStd_Array1OfInteger Mults(1, hBSpline->NbKnots());
    hBSpline->Multiplicities(Mults);
    TColStd_Array1OfReal Knots(1, hBSpline->NbKnots());
    hBSpline->Knots(Knots);
    TColgp_Array1OfPnt2d Poles(1, hBSpline->NbPoles());
    hBSpline->Poles(Poles);
    TColStd_Array1OfReal Weights(1, hBSpline->NbPoles());
    hBSpline->Weights(Weights);

    int NewNbKnots, NewNbPoles;
    BSplCLib::PrepareUnperiodize(hBSpline->Degree(), Mults, NewNbKnots,
                                 NewNbPoles);

    TColStd_Array1OfInteger NewMults(1, NewNbKnots);
    TColStd_Array1OfReal NewKnots(1, NewNbKnots);
    TColgp_Array1OfPnt2d NewPoles(1, NewNbPoles);
    TColStd_Array1OfReal NewWeights(1, NewNbPoles);

    BSplCLib::Unperiodize(hBSpline->Degree(), Mults, Knots, Poles, &Weights,
                          NewMults, NewKnots, NewPoles, &NewWeights);

    Handle(Geom2d_Curve) hnCurve = new Geom2d_BSplineCurve(NewPoles, NewWeights,
                                                           NewKnots, NewMults,
                                                           hBSpline->Degree());
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_flattenBSpline)!\n", stat);
      return stat;
    }
    obj->oclass         = PCURVE;
    obj->mtype          = BSPLINE;
    egadsPCurve *pcurvn = new egadsPCurve;
    pcurvn->handle      = hnCurve;
    pcurvn->ref         = NULL;
    pcurvn->topFlg      = 0;
    pcurvn->header      = NULL;
    pcurvn->data        = NULL;
    pcurvn->data_dot    = NULL;
    obj->blind          = pcurvn;
    EG_getGeometry(obj, &ot, &mc, &ref, &pcurvn->header, &pcurvn->data);
    EG_getGeometryLen(obj, &mc, &pcurvn->dataLen);
    pcurvn->trange[0]   = hnCurve->FirstParameter();
    pcurvn->trange[1]   = hnCurve->LastParameter();
    EG_referenceObject(obj, context);
    *result = obj;

  } else if (object->oclass == CURVE) {

    egadsCurve *pcurve = (egadsCurve *) object->blind;
    if ((pcurve->header[0]&4) == 0) {
      printf(" EGADS Warning: Already flat (EG_flattenBSpline)!\n");
      return EGADS_GEOMERR;
    }
    Handle(Geom_Curve)        hCurve   = pcurve->handle;
    Handle(Geom_BSplineCurve) hBSpline =
                                    Handle(Geom_BSplineCurve)::DownCast(hCurve);

    TColStd_Array1OfInteger Mults(1, hBSpline->NbKnots());
    hBSpline->Multiplicities(Mults);
    TColStd_Array1OfReal Knots(1, hBSpline->NbKnots());
    hBSpline->Knots(Knots);
    TColgp_Array1OfPnt Poles(1, hBSpline->NbPoles());
    hBSpline->Poles(Poles);
    TColStd_Array1OfReal Weights(1, hBSpline->NbPoles());
    hBSpline->Weights(Weights);

    int NewNbKnots, NewNbPoles;
    BSplCLib::PrepareUnperiodize(hBSpline->Degree(), Mults, NewNbKnots,
                                 NewNbPoles);

    TColStd_Array1OfInteger NewMults(1, NewNbKnots);
    TColStd_Array1OfReal NewKnots(1, NewNbKnots);
    TColgp_Array1OfPnt NewPoles(1, NewNbPoles);
    TColStd_Array1OfReal NewWeights(1, NewNbPoles);

    BSplCLib::Unperiodize(hBSpline->Degree(), Mults, Knots, Poles, &Weights,
                          NewMults, NewKnots, NewPoles, &NewWeights);

    Handle(Geom_Curve) hnCurve = new Geom_BSplineCurve(NewPoles, NewWeights,
                                                       NewKnots, NewMults,
                                                       hBSpline->Degree());
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_flattenBSpline)!\n", stat);
      return stat;
    }
    obj->oclass        = CURVE;
    obj->mtype         = BSPLINE;
    egadsCurve *pcurvn = new egadsCurve;
    pcurvn->handle     = hnCurve;
    pcurvn->ref        = NULL;
    pcurvn->topFlg     = 0;
    pcurvn->header     = NULL;
    pcurvn->data       = NULL;
    pcurvn->data_dot   = NULL;
    obj->blind         = pcurvn;
    EG_getGeometry(obj, &ot, &mc, &ref, &pcurvn->header, &pcurvn->data);
    EG_getGeometryLen(obj, &mc, &pcurvn->dataLen);
    pcurvn->trange[0]  = hnCurve->FirstParameter();
    pcurvn->trange[1]  = hnCurve->LastParameter();
    EG_referenceObject(obj, context);
    *result = obj;

  } else {

    egadsSurface *psurface = (egadsSurface *) object->blind;
    if ((psurface->header[0]&12) == 0) {
      printf(" EGADS Warning: Already flat (EG_flattenBSpline)!\n");
      return EGADS_GEOMERR;
    }
    Handle(Geom_Surface) hSurf = psurface->handle;
    Handle(Geom_BSplineSurface)
                        hBSpline = Handle(Geom_BSplineSurface)::DownCast(hSurf);

    TColStd_Array1OfInteger UMults(1, hBSpline->NbUKnots());
    hBSpline->UMultiplicities(UMults);
    TColStd_Array1OfReal UKnots(1, hBSpline->NbUKnots());
    hBSpline->UKnots(UKnots);
    TColStd_Array1OfInteger VMults(1, hBSpline->NbVKnots());
    hBSpline->VMultiplicities(VMults);
    TColStd_Array1OfReal VKnots(1, hBSpline->NbVKnots());
    hBSpline->VKnots(VKnots);
    TColgp_Array2OfPnt Poles(1, hBSpline->NbUPoles(), 1, hBSpline->NbVPoles());
    hBSpline->Poles(Poles);
    TColStd_Array2OfReal Weights(1, hBSpline->NbUPoles(), 1, hBSpline->NbVPoles());
    hBSpline->Weights(Weights);

    NewNbUKnots = hBSpline->NbUKnots();
    NewNbUPoles = hBSpline->NbUPoles();
    NewNbVKnots = hBSpline->NbVKnots();
    NewNbVPoles = hBSpline->NbVPoles();
    if ((psurface->header[0]&4) != 0) {
      BSplCLib::PrepareUnperiodize(hBSpline->UDegree(), UMults, NewNbUKnots,
                                   NewNbUPoles);

      TColStd_Array1OfInteger newUMults(1, NewNbUKnots);
      TColStd_Array1OfReal newUKnots(1, NewNbUKnots);
      TColgp_Array2OfPnt newPoles(1, NewNbUPoles, 1, hBSpline->NbVPoles());
      TColStd_Array2OfReal newWeights(1, NewNbUPoles, 1, hBSpline->NbVPoles());

      BSplSLib::Unperiodize(Standard_True, hBSpline->UDegree(), UMults,
                            UKnots, Poles, &Weights, newUMults, newUKnots,
                            newPoles, &newWeights);
      UMults  = newUMults;
      UKnots  = newUKnots;
      Poles   = newPoles;
      Weights = newWeights;
    }

    if ((psurface->header[0]&8) != 0) {
      BSplCLib::PrepareUnperiodize(hBSpline->VDegree(), VMults, NewNbVKnots,
                                   NewNbVPoles);

      TColStd_Array1OfInteger newVMults(1, NewNbVKnots);
      TColStd_Array1OfReal newVKnots(1, NewNbVKnots);
      TColgp_Array2OfPnt newPoles(1, NewNbUPoles, 1, NewNbVPoles);
      TColStd_Array2OfReal newWeights(1, NewNbUPoles, 1, NewNbVPoles);

      BSplSLib::Unperiodize(Standard_False, hBSpline->VDegree(), VMults,
                            VKnots, Poles, &Weights, newVMults, newVKnots,
                            newPoles, &newWeights);
      VMults  = newVMults;
      VKnots  = newVKnots;
      Poles   = newPoles;
      Weights = newWeights;
    }

    Handle(Geom_Surface) hSurfn;
    try {
      hSurfn = new Geom_BSplineSurface(Poles, Weights, UKnots, VKnots,
                                      UMults, VMults, hBSpline->UDegree(),
                                                      hBSpline->VDegree());
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Warning: Geometry Creation Error (EG_flattenBSpline)!\n");
      printf("                %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Warning: Geometry Creation Error (EG_flattenBSpline)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Surface = %d (EG_flattenBSpline)!\n", stat);
      return stat;
    }
    obj->oclass         = SURFACE;
    obj->mtype          = BSPLINE;
    egadsSurface *psurf = new egadsSurface;
    psurf->handle       = hSurfn;
    psurf->ref          = NULL;
    psurf->topFlg       = 0;
    psurf->header       = NULL;
    psurf->data         = NULL;
    psurf->data_dot     = NULL;
    obj->blind          = psurf;
    EG_getGeometry(obj, &ot, &mc, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(obj, &mc, &psurf->dataLen);
    hSurfn->Bounds(psurf->urange[0], psurf->urange[1],
                   psurf->vrange[0], psurf->vrange[1]);
    EG_referenceObject(obj, context);
    *result = obj;
  }

  return EGADS_SUCCESS;
}


int
EG_addKnots(const egObject *object, int nU, /*@null@*/ double *Us,
            int nV, /*@null@*/ double *Vs, egObject **result)
{
  int      i, stat, ot, mc, len, nmult;
  egObject *obj, *context, *ref;

  *result = NULL;
  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != CURVE) &&
      (object->oclass != SURFACE))   return EGADS_NOTGEOM;
  if  (object->mtype != BSPLINE)     return EGADS_NOTGEOM;
  if  (object->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(object))        return EGADS_CNTXTHRD;

  context = EG_context(object);
  gp_Trsf form = gp_Trsf();

  if (object->oclass == CURVE) {

    egadsCurve *pcurve = (egadsCurve *) object->blind;
    Handle(Geom_Curve)    hCurve = pcurve->handle;
    Handle(Geom_Geometry) nGeom  = hCurve->Transformed(form);
    Handle(Geom_Curve)    nCurve = Handle(Geom_Curve)::DownCast(nGeom);
    Handle(Geom_BSplineCurve)
                        hBSpline = Handle(Geom_BSplineCurve)::DownCast(nCurve);

    if (nU <= 0) {
      try {
        hBSpline->IncreaseDegree(3);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Geometry Creation Error Elev (EG_addKnots)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Warning: Geometry Creation Error Elev (EG_addKnots)!\n");
        return EGADS_GEOMERR;
      }
    } else {
      for (len = i = 1; i < nU; i++)
        if (fabs(Us[i]-Us[i-1]) > KNACC) len++;
      TColStd_Array1OfReal    aKnots(1, len);
      TColStd_Array1OfInteger aMults(1, len);
      aKnots(1) = Us[0];
      aMults(1) = nmult = 1;
      for (len = i = 1; i < nU; i++)
        if (fabs(Us[i]-Us[i-1]) > KNACC) {
          len++;
          aKnots(len) = Us[i];
          aMults(len) = nmult = 1;
        } else {
          nmult++;
          aMults(len) = nmult;
        }
      try {
        hBSpline->IncreaseDegree(3);
        hBSpline->InsertKnots(aKnots, aMults, KNACC, Standard_True);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Geometry Creation Error (EG_addKnots)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Warning: Geometry Creation Error (EG_addKnots)!\n");
        return EGADS_GEOMERR;
      }
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_addKnots)!\n", stat);
      return stat;
    }
    obj->oclass        = CURVE;
    obj->mtype         = BSPLINE;
    egadsCurve *pcurvn = new egadsCurve;
    pcurvn->handle     = hBSpline;
    pcurvn->ref        = NULL;
    pcurvn->topFlg     = 0;
    pcurvn->header     = NULL;
    pcurvn->data       = NULL;
    pcurvn->data_dot   = NULL;
    obj->blind         = pcurvn;
    EG_getGeometry(obj, &ot, &mc, &ref, &pcurvn->header, &pcurvn->data);
    EG_getGeometryLen(obj, &mc, &pcurvn->dataLen);
    pcurvn->trange[0]  = hBSpline->FirstParameter();
    pcurvn->trange[1]  = hBSpline->LastParameter();
    EG_referenceObject(obj, context);
    *result = obj;

  } else {

    egadsSurface *psurface = (egadsSurface *) object->blind;
    Handle(Geom_Surface)  hSurf = psurface->handle;
    Handle(Geom_Geometry) nGeom = hSurf->Transformed(form);
    Handle(Geom_Surface)  nSurf = Handle(Geom_Surface)::DownCast(nGeom);
    Handle(Geom_BSplineSurface)
                       hBSpline = Handle(Geom_BSplineSurface)::DownCast(nSurf);

    if ((nU <= 0) && (nV <= 0)) {
      try {
        hBSpline->IncreaseDegree(3, 3);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Geometry Creation Error Elev (EG_addKnots)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Warning: Geometry Creation Error Elev (EG_addKnots)!\n");
        return EGADS_GEOMERR;
      }
    }

    if ((nU > 0) && (Us != NULL)) {
      for (len = i = 1; i < nU; i++)
        if (fabs(Us[i]-Us[i-1]) > KNACC) len++;
      TColStd_Array1OfReal    aKnots(1, len);
      TColStd_Array1OfInteger aMults(1, len);
      aKnots(1) = Us[0];
      aMults(1) = nmult = 1;
      for (len = i = 1; i < nU; i++)
        if (fabs(Us[i]-Us[i-1]) > KNACC) {
          len++;
          aKnots(len) = Us[i];
          aMults(len) = nmult = 1;
        } else {
          nmult++;
          aMults(len) = nmult;
        }
      try {
        hBSpline->IncreaseDegree(3, 3);
        hBSpline->InsertUKnots(aKnots, aMults, KNACC, Standard_True);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Geometry Creation Error U (EG_addKnots)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Warning: Geometry Creation Error U (EG_addKnots)!\n");
        return EGADS_GEOMERR;
      }
    }

    if ((nV > 0) && (Vs != NULL)) {
      for (len = i = 1; i < nV; i++)
        if (fabs(Vs[i]-Vs[i-1]) > KNACC) len++;
      TColStd_Array1OfReal    aKnots(1, len);
      TColStd_Array1OfInteger aMults(1, len);
      aKnots(1) = Vs[0];
      aMults(1) = nmult = 1;
      for (len = i = 1; i < nV; i++)
        if (fabs(Vs[i]-Vs[i-1]) > KNACC) {
          len++;
          aKnots(len) = Vs[i];
          aMults(len) = nmult = 1;
        } else {
          nmult++;
          aMults(len) = nmult;
        }
      try {
        if ((nU <= 0) || (Us == NULL)) hBSpline->IncreaseDegree(3, 3);
        hBSpline->InsertVKnots(aKnots, aMults, KNACC, Standard_True);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Geometry Creation Error V (EG_addKnots)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Warning: Geometry Creation Error V (EG_addKnots)!\n");
        return EGADS_GEOMERR;
      }
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Surface = %d (EG_addKnots)!\n", stat);
      return stat;
    }
    obj->oclass         = SURFACE;
    obj->mtype          = BSPLINE;
    egadsSurface *psurf = new egadsSurface;
    psurf->handle       = hBSpline;
    psurf->ref          = NULL;
    psurf->topFlg       = 0;
    psurf->header       = NULL;
    psurf->data         = NULL;
    psurf->data_dot     = NULL;
    obj->blind          = psurf;
    EG_getGeometry(obj, &ot, &mc, &ref, &psurf->header, &psurf->data);
    EG_getGeometryLen(obj, &mc, &psurf->dataLen);
    hBSpline->Bounds(psurf->urange[0], psurf->urange[1],
                     psurf->vrange[0], psurf->vrange[1]);
    EG_referenceObject(obj, context);
    *result = obj;
  }

  return EGADS_SUCCESS;
}


int
EG_mapSequen(egObject *src, egObject *dst, egObject **result)
{
  int      i, j, hit, outLevel, stat;
  egObject *context, *obj;

  if  (src == NULL)                return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if ((src->oclass != PCURVE)  && (src->oclass != CURVE) &&
      (src->oclass != SURFACE))    return EGADS_NOTGEOM;
  if  (src->blind == NULL)         return EGADS_NODATA;
  if  (EG_sameThread(src))         return EGADS_CNTXTHRD;
  if  (dst == NULL)                return EGADS_NULLOBJ;
  if  (dst->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if  (dst->oclass != src->oclass) return EGADS_NOTGEOM;
  if  (dst->mtype  != src->mtype)  return EGADS_GEOMERR;
  if  (dst->blind == NULL)         return EGADS_NODATA;
  if  (EG_sameThread(dst))         return EGADS_CNTXTHRD;

  *result  = NULL;
  outLevel = EG_outLevel(dst);
  context  = EG_context(dst);

  if (src->oclass == PCURVE) {

    if (src->mtype != BSPLINE) return EGADS_GEOMERR;
    egadsPCurve *ppcurv            = (egadsPCurve *) src->blind;
    Handle(Geom2d_Curve)    hCurve = ppcurv->handle;
    Handle(Geom2d_BSplineCurve)
                          hBSpline = Handle(Geom2d_BSplineCurve)::DownCast(hCurve);
    egadsPCurve *ppcurvd           = (egadsPCurve *) dst->blind;
    Handle(Geom2d_Curve)    hCurvt = ppcurvd->handle;
    gp_Trsf2d                 form = gp_Trsf2d();
    Handle(Geom2d_Geometry)  nGeom = hCurvt->Transformed(form);
    Handle(Geom2d_Curve)    hCurvd = Handle(Geom2d_Curve)::DownCast(nGeom);
    Handle(Geom2d_BSplineCurve)
                          hBSplind = Handle(Geom2d_BSplineCurve)::DownCast(hCurvd);
    if (hBSpline->NbKnots() != hBSplind->NbKnots()) return EGADS_INDEXERR;
    int len = hBSpline->NbKnots();
    TColStd_Array1OfInteger mults(1, len);
    TColStd_Array1OfInteger multd(1, len);
    hBSpline->Multiplicities(mults);
    hBSplind->Multiplicities(multd);
    for (i = 1; i <= len; i++)
      if (mults(i) != multd(i)) return EGADS_RANGERR;

    TColStd_Array1OfReal knots(1, len);
    TColStd_Array1OfReal knotd(1, len);
    hBSpline->Knots(knots);
    hBSplind->Knots(knotd);
    hit = 0;
    for (i = 2; i < len; i++) {
      double scaledKnot = (knots(i)-knots(1))/(knots(len)-knots(1));
      scaledKnot *= knotd(len)-knotd(1);
      scaledKnot += knotd(1);
      for (j = 2; j < len; j++)
        if (fabs(scaledKnot-knotd(j)) <= 0.5*KNDIFF) break;
      if (j != len) continue;
      hBSplind->InsertKnot(scaledKnot);
      hit++;
      if (outLevel > 1)
        printf("   inserting knot = %lf (%lf)\n", scaledKnot, knots(i));
    }
#ifdef KNOTREMOVE
    if (hit != 0) {
      int nknot = hBSplind->NbKnots();
      for (j = nknot-1; j >= 2; j--) {
        double knot = hBSplind->Knot(j);
        for (i = 2; i < len; i++) {
          double scaledKnot = (knots(i)-knots(1))/(knots(len)-knots(1));
          scaledKnot *= knotd(len)-knotd(1);
          scaledKnot += knotd(1);
          if (fabs(scaledKnot-knot) <= 0.5*KNDIFF) break;
        }
        if (i != len) continue;
        if (outLevel > 1)
          printf("   removing  knot = %lf\n", hBSplind->Knot(j));
        if (!hBSplind->RemoveKnot(j, 0, 1.0))
          if (outLevel > 0)
            printf(" EGADS Warning: Cannot remove knot %lf (EG_mapSequen)!",
                   hBSplind->Knot(j));
      }
    }
#endif
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: PCurve makeObject = %d (EG_mapSequen)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completePCurve(obj, hCurvd);

  } else if (src->oclass == CURVE) {

    if (src->mtype != BSPLINE) return EGADS_GEOMERR;
    egadsCurve *pcurve           = (egadsCurve *) src->blind;
    Handle(Geom_Curve)    hCurve = pcurve->handle;
    Handle(Geom_BSplineCurve)
                        hBSpline = Handle(Geom_BSplineCurve)::DownCast(hCurve);
    egadsCurve *pcurvd           = (egadsCurve *) dst->blind;
    Handle(Geom_Curve)    hCurvt = pcurvd->handle;
    gp_Trsf                 form = gp_Trsf();
    Handle(Geom_Geometry)  nGeom = hCurvt->Transformed(form);
    Handle(Geom_Curve)    hCurvd = Handle(Geom_Curve)::DownCast(nGeom);
    Handle(Geom_BSplineCurve)
                        hBSplind = Handle(Geom_BSplineCurve)::DownCast(hCurvd);
    if (hBSpline->NbKnots() != hBSplind->NbKnots()) return EGADS_INDEXERR;
    int len = hBSpline->NbKnots();
    TColStd_Array1OfInteger mults(1, len);
    TColStd_Array1OfInteger multd(1, len);
    hBSpline->Multiplicities(mults);
    hBSplind->Multiplicities(multd);
    for (i = 1; i <= len; i++)
      if (mults(i) != multd(i)) return EGADS_RANGERR;
/*  GeomAbs_BSplKnotDistribution kDist = hBSplind->KnotDistribution();
    if (kDist == GeomAbs_NonUniform)      printf(" NonUniform!\n");
    if (kDist == GeomAbs_Uniform)         printf(" Uniform!\n");
    if (kDist == GeomAbs_QuasiUniform)    printf(" QuasiUniform!\n");
    if (kDist == GeomAbs_PiecewiseBezier) printf(" PiecewiseBezier!\n");  */

    TColStd_Array1OfReal knots(1, len);
    TColStd_Array1OfReal knotd(1, len);
    hBSpline->Knots(knots);
    hBSplind->Knots(knotd);
    hit = 0;
    for (i = 2; i < len; i++) {
      double scaledKnot = (knots(i)-knots(1))/(knots(len)-knots(1));
      scaledKnot *= knotd(len)-knotd(1);
      scaledKnot += knotd(1);
      for (j = 2; j < len; j++)
        if (fabs(scaledKnot-knotd(j)) <= 0.5*KNDIFF) break;
      if (j != len) continue;
      hBSplind->InsertKnot(scaledKnot);
      hit++;
      if (outLevel > 1)
        printf("   inserting knot = %lf (%lf)\n", scaledKnot, knots(i));
    }
#ifdef KNOTREMOVE
    if (hit != 0) {
      int nknot = hBSplind->NbKnots();
      for (j = nknot-1; j >= 2; j--) {
        double knot = hBSplind->Knot(j);
        for (i = 2; i < len; i++) {
          double scaledKnot = (knots(i)-knots(1))/(knots(len)-knots(1));
          scaledKnot *= knotd(len)-knotd(1);
          scaledKnot += knotd(1);
          if (fabs(scaledKnot-knot) <= 0.5*KNDIFF) break;
        }
        if (i != len) continue;
        if (outLevel > 1)
          printf("   removing  knot = %lf\n", hBSplind->Knot(j));
        if (!hBSplind->RemoveKnot(j, 0, 1.0))
          if (outLevel > 0)
            printf(" EGADS Warning: Cannot remove knot %lf (EG_mapSequen)!",
                   hBSplind->Knot(j));
      }
    }
#endif
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Curve makeObject = %d (EG_mapSequen)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeCurve(obj, hCurvd);

  } else {

    if (src->mtype != BSPLINE) return EGADS_GEOMERR;
    egadsSurface *psurface       = (egadsSurface *) src->blind;
    Handle(Geom_Surface)  hSurf  = psurface->handle;
    Handle(Geom_BSplineSurface)
                        hBSpline = Handle(Geom_BSplineSurface)::DownCast(hSurf);
    egadsSurface *psurfacd       = (egadsSurface *) dst->blind;
    Handle(Geom_Surface)  hSurft = psurfacd->handle;
    gp_Trsf                 form = gp_Trsf();
    Handle(Geom_Geometry)  nGeom = hSurft->Transformed(form);
    Handle(Geom_Surface)  hSurfd = Handle(Geom_Surface)::DownCast(nGeom);
    Handle(Geom_BSplineSurface)
                        hBSplind = Handle(Geom_BSplineSurface)::DownCast(hSurfd);
    if (hBSpline->NbUKnots() != hBSplind->NbUKnots()) return EGADS_INDEXERR;
    if (hBSpline->NbVKnots() != hBSplind->NbVKnots()) return EGADS_INDEXERR;
    int uLen = hBSpline->NbUKnots();
    int vLen = hBSpline->NbVKnots();
    TColStd_Array1OfInteger uMults(1, uLen);
    TColStd_Array1OfInteger uMultd(1, uLen);
    TColStd_Array1OfInteger vMults(1, vLen);
    TColStd_Array1OfInteger vMultd(1, vLen);
    hBSpline->UMultiplicities(uMults);
    hBSplind->UMultiplicities(uMultd);
    for (i = 1; i <= uLen; i++)
      if (uMults(i) != uMultd(i)) return EGADS_RANGERR;
    hBSpline->VMultiplicities(vMults);
    hBSplind->VMultiplicities(vMultd);
    for (i = 1; i <= vLen; i++)
      if (vMults(i) != vMultd(i)) return EGADS_RANGERR;

    TColStd_Array1OfReal uKnots(1, uLen);
    TColStd_Array1OfReal vKnots(1, vLen);
    TColStd_Array1OfReal uKnotd(1, uLen);
    TColStd_Array1OfReal vKnotd(1, vLen);
    hBSpline->UKnots(uKnots);
    hBSpline->VKnots(vKnots);
    hBSplind->UKnots(uKnotd);
    hBSplind->VKnots(vKnotd);

    // u knots
    hit = 0;
    for (i = 2; i < uLen; i++) {
      double scaledKnot = (uKnots(i)-uKnots(1))/(uKnots(uLen)-uKnots(1));
      scaledKnot *= uKnotd(uLen)-uKnotd(1);
      scaledKnot += uKnotd(1);
      for (j = 2; j < uLen; j++)
        if (fabs(scaledKnot-uKnotd(j)) <= 0.5*KNDIFF) break;
      if (j != uLen) continue;
      hBSplind->InsertUKnot(scaledKnot, 1, 0.5*KNDIFF);
      hit++;
      if (outLevel > 1)
        printf("   inserting u knot = %lf (%lf)\n", scaledKnot, uKnots(i));
    }
#ifdef KNOTREMOVE
    if (hit != 0) {
      int nknot = hBSplind->NbUKnots();
      for (j = nknot-1; j >= 2; j--) {
        double knot = hBSplind->UKnot(j);
        for (i = 2; i < uLen; i++) {
          double scaledKnot = (uKnots(i)-uKnots(1))/(uKnots(uLen)-uKnots(1));
          scaledKnot *= uKnotd(uLen)-uKnotd(1);
          scaledKnot += uKnotd(1);
          if (fabs(scaledKnot-knot) <= 0.5*KNDIFF) break;
        }
        if (i != uLen) continue;
        if (outLevel > 1)
          printf("   removing  u knot = %lf\n", hBSplind->UKnot(j));
        if (!hBSplind->RemoveUKnot(j, 0, 1.0))
          if (outLevel > 0)
            printf(" EGADS Warning: Cannot remove U knot %lf (EG_mapSequen)!",
                   hBSplind->UKnot(j));
      }
    }
#endif

    // v knots
    hit = 0;
    for (i = 2; i < vLen; i++) {
      double scaledKnot = (vKnots(i)-vKnots(1))/(vKnots(vLen)-vKnots(1));
      scaledKnot *= vKnotd(vLen)-vKnotd(1);
      scaledKnot += vKnotd(1);
      for (j = 2; j < vLen; j++)
        if (fabs(scaledKnot-vKnotd(j)) <= 0.5*KNDIFF) break;
      if (j != vLen) continue;
      hBSplind->InsertVKnot(scaledKnot, 1, 0.5*KNDIFF);
      hit++;
      if (outLevel > 1)
        printf("   inserting v knot = %lf (%lf)\n", scaledKnot, vKnots(i));
    }
#ifdef KNOTREMOVE
    if (hit != 0) {
      int nknot = hBSplind->NbVKnots();
      for (j = nknot-1; j >= 2; j--) {
        double knot = hBSplind->VKnot(j);
        for (i = 2; i < vLen; i++) {
          double scaledKnot = (vKnots(i)-vKnots(1))/(vKnots(vLen)-vKnots(1));
          scaledKnot *= vKnotd(vLen)-vKnotd(1);
          scaledKnot += vKnotd(1);
          if (fabs(scaledKnot-knot) <= 0.5*KNDIFF) break;
        }
        if (i != vLen) continue;
        if (outLevel > 1)
          printf("   removing  v knot = %lf\n", hBSplind->VKnot(j));
        if (!hBSplind->RemoveVKnot(j, 0, 1.0))
          if (outLevel > 0)
            printf(" EGADS Warning: Cannot remove V knot %lf (EG_mapSequen)!",
                   hBSplind->VKnot(j));
      }
    }
#endif
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface makeObject = %d (EG_mapSequen)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeSurf(obj, hSurfd);

  }

  *result = obj;
  return EGADS_SUCCESS;
}


void
EG_mapTessTs(egTess1D src, egTess1D dst)
{
  double t;

  if (src.npts != dst.npts) {
    printf(" EGADS Warning: Len Mismatch src = %d, dst = %d (EG_mapTessTs)!\n",
           src.npts, dst.npts);
    return;
  }
  egadsEdge *pedgs = (egadsEdge *) src.obj->blind;
  if (pedgs == NULL) {
    printf(" EGADS Warning: NULL src Edge Object (EG_mapTessTs)!\n");
    return;
  }
  egObject  *curvs = pedgs->curve;
  if (curvs == NULL) {
    printf(" EGADS Warning: No curve Object for src Edge (EG_mapTessTs)!\n");
    return;
  }
  egadsCurve *pcurvs = (egadsCurve *) curvs->blind;
  if (pcurvs == NULL) {
    printf(" EGADS Warning: No curve Data for src Edge (EG_mapTessTs)!\n");
    return;
  }
  egadsEdge *pedgd = (egadsEdge *) dst.obj->blind;
  if (pedgd == NULL) {
    printf(" EGADS Warning: NULL dst Edge Object (EG_mapTessTs)!\n");
    return;
  }
  egObject  *curvd = pedgd->curve;
  if (curvd == NULL) {
    printf(" EGADS Warning: No curve Object for dst Edge (EG_mapTessTs)!\n");
    return;
  }
  egadsCurve *pcurvd = (egadsCurve *) curvd->blind;
  if (pcurvd == NULL) {
    printf(" EGADS Warning: No curve Data for dst Edge (EG_mapTessTs)!\n");
    return;
  }

  int n = src.npts;
  GeomAdaptor_Curve ACsrc(pcurvs->handle);
  GeomAdaptor_Curve ACdst(pcurvd->handle);
  double slen = GCPnts_AbscissaPoint::Length(ACsrc, src.t[0], src.t[n-1]);
  double dlen = GCPnts_AbscissaPoint::Length(ACdst, dst.t[0], dst.t[n-1]);

  //
  // have the relative arcLengths in the destination match the source
  for (int i = 1; i < n-1; i++) {
    double srcAlen = GCPnts_AbscissaPoint::Length(ACsrc, src.t[0], src.t[i]);
    double tgtAlen = dlen*srcAlen/slen;
    try {
      GCPnts_AbscissaPoint AP(ACdst, tgtAlen, dst.t[0], dst.t[i]);
      if (!AP.IsDone()) continue;
      t = AP.Parameter();
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Warning: GCPnts_AbscissaPoint (EG_mapTessTs)!\n");
      printf("                %s\n", e.GetMessageString());
      continue;
    }
    catch (...) {
      printf(" EGADS Warning: GCPnts_AbscissaPoint (EG_mapTessTs)!\n");
      continue;
    }
/*  printf("    %d:  %lf %lf\n", i+1, dst.t[i], t);  */
    gp_Pnt P0;
    pcurvd->handle->D0(t, P0);
    dst.t[i]       = t;
    dst.xyz[3*i  ] = P0.X();
    dst.xyz[3*i+1] = P0.Y();
    dst.xyz[3*i+2] = P0.Z();
  }
}


int
EG_relPosTs(egObject *geom, int n, const double *rel, double *ts, double *xyzs)
{
  egadsCurve *pcurve;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != CURVE)  &&
      (geom->oclass != EDGE))      return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;

  if (geom->oclass == CURVE) {

    pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;

  } else {

    if (geom->mtype == DEGENERATE) return EGADS_DEGEN;
    egadsEdge *pedge  = (egadsEdge *) geom->blind;
    if (pedge  == NULL) return EGADS_NULLOBJ;
    egObject  *curvo  = pedge->curve;
    if (curvo  == NULL) return EGADS_NULLOBJ;
    pcurve = (egadsCurve *) curvo->blind;

  }
  if (pcurve == NULL) return EGADS_NULLOBJ;
  GeomAdaptor_Curve AC(pcurve->handle);
  double alen = GCPnts_AbscissaPoint::Length(AC, ts[0], ts[n-1]);
  if (alen == 0.0) {
    printf(" EGADS Error: ArcLength of Segment is Zero (EG_relPosTs)!\n");
    return EGADS_GEOMERR;
  }
  if (rel == NULL) {
    for (int i = 1; i < n-1; i++) {
      GCPnts_AbscissaPoint AP(AC, i*alen/(n-1), ts[0]);
      if (!AP.IsDone()) continue;
      double t    = AP.Parameter();
      gp_Pnt P0;
      pcurve->handle->D0(t, P0);
      ts[i]       = t;
      xyzs[3*i  ] = P0.X();
      xyzs[3*i+1] = P0.Y();
      xyzs[3*i+2] = P0.Z();
    }
  } else {
    for (int i = 1; i < n-1; i++) {
      GCPnts_AbscissaPoint AP(AC, rel[i-1]*alen, ts[0]);
      if (!AP.IsDone()) continue;
      double t    = AP.Parameter();
/*    printf("    %d:  %lf %lf\n", i, ts[i], t);  */
      gp_Pnt P0;
      pcurve->handle->D0(t, P0);
      ts[i]       = t;
      xyzs[3*i  ] = P0.X();
      xyzs[3*i+1] = P0.Y();
      xyzs[3*i+2] = P0.Z();
    }
  }

  return EGADS_SUCCESS;
}
