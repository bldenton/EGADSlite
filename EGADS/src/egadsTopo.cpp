/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Topology Functions
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
#include "egadsClasses.h"

#define OCC_SOLIDS

#define UVTOL    1.e-4

/* structures */

  typedef struct {
    int    nFace;		/* the number of Faces in the Edge */
    int    seq;                 /* the sequence number */
    int    *fIndices;           /* the Face indices (nFace in len) */
    char   *ID;                 /* the ID including seq # */
    double CG[4];               /* Center of Gravity, then Length */
  } edgeID;


  extern "C" int  EG_destroyTopology( egObject *topo );
  extern "C" int  EG_copyObject( const egObject *object, /*@null@*/ void *ptr,
                                       egObject **copy );
  extern "C" int  EG_fullAttrs( const egObject *obj );
  extern "C" int  EG_isSame( const egObject *obj1, const egObject *obj2 );
  extern "C" int  EG_attributeRet( const egObject *obj, const char *name,
                                   int *atype, int *len,
                                   /*@null@*/ const int **ints,
                                   /*@null@*/ const double **reals,
                                   /*@null@*/ const char **str );
  extern "C" int  EG_attributeAdd( egObject *obj, const char *name, int atype,
                                   int len, /*@null@*/ const int    *ints,
                                            /*@null@*/ const double *reals,
                                            /*@null@*/ const char   *str );
  extern "C" int  EG_mapKnots( egObject *src, egObject *dst, egObject **result );

  extern "C" int  EG_getGeometry( const egObject *geom, int *oclass, int *mtype,
                                  egObject **refGeom, /*@null@*/ int **ivec,
                                  /*@null@*/ double **rvec );
  extern "C" int  EG_otherCurve( const egObject *surface, const egObject *curve,
                                 double tol, egObject **newcurve );
  extern "C" int  EG_evaluate( const egObject *geom,
                               /*@null@*/ const double *param, double *result );
  extern "C" int  EG_invEvaluate( const egObject *obj, double *xyz,
                                  double *param, double *results );

  extern "C" int  EG_tolerance( const egObject *topo, double *tol );
  extern "C" int  EG_getTolerance( const egObject *topo, double *tol );
  extern "C" int  EG_setTolerance( const egObject *topo, double  tol );
  extern "C" int  EG_delSmallEdges( const egObject *body,
                                    double tol, egObject **newBody );
  extern "C" int  EG_getTopology( const egObject *topo, egObject **geom,
                                  int *oclass, int *type,
                                  /*@null@*/ double *limits, int *nChildren,
                                  egObject ***children, int **senses );
  extern "C" int  EG_makeTopology( egObject *context, /*@null@*/ egObject *geom,
                                   int oclass, int mtype, /*@null@*/ double *limits,
                                   int nChildren, /*@null@*/ egObject **children,
                                   /*@null@*/ int *senses, egObject **topo );
  extern "C" int  EG_makeLoop( int nedge, egObject **edges,
                               /*@null@*/ egObject *geom, double toler,
                               egObject **result );
  extern "C" int  EG_makeFace( egObject *object, int mtype,
                               /*@null@*/ const double *limits, egObject **face );
  extern "C" int  EG_getPlane( const egObject *object, egObject **plane );
  extern "C" int  EG_getArea( egObject *object, /*@null@*/ const double *limits,
                              double *area );
  extern "C" int  EG_getUVbox( const egObject *face, const egObject *loop,
                               double *box );
  extern "C" int  EG_getUVinfo( egObject *face, egObject *loop, double *box,
                                double *area );
  extern "C" int  EG_getBodyTopos( const egObject *body, /*@null@*/ egObject *src,
                                   int oclass, int *ntopo,
                                   /*@null@*/ egObject ***topos );
  extern "C" int  EG_indexBodyTopo( const egObject *body, const egObject *src );
  extern "C" int  EG_objectBodyTopo( const egObject *body, int oclass, int index,
                                     egObject **obj );
  extern "C" int  EG_sameBodyTopo( const egObject *bod1, const egObject *bod2 );
  extern "C" int  EG_makeSolidBody( egObject *context, int stype,
                                    const double *rvec, egObject **body );
  extern "C" int  EG_makeSolidBody_dot( egObject *body, int stype,
                                        const double *rvec,
                                        const double *rvec_dot );
  extern "C" int  EG_getBoundingBox( const egObject *topo, double *box );
  extern "C" int  EG_getMassProperties( const egObject *topo,
                                        /*@null@*/ double *props );
  extern "C" int  EG_isEquivalent( const egObject *topo1, const egObject *topo2 );
  extern "C" int  EG_isPlanar( const egObject *topo );
  extern "C" int  EG_getEdgeUV( const egObject *face, const egObject *edge,
                                int sense, double t, double *result );
  extern "C" int  EG_getEdgeUVs( const egObject *face, const egObject *edge,
                                 int sense, int nt, const double *ts,
                                 double *uvs );
  extern "C" int  EG_getEdgeUVeval( const egObject *face, const egObject *edge,
                                    int sense, double t, double *result );
  extern "C" int  EG_invEdgeUV( const egObject *face, const egObject *edge,
                                int sense, double *uv, double *t, double *uvs );
  extern "C" int  EG_getPCurve( const egObject *face, const egObject *edge,
                                int sense, int *mtype, int **header,
                                double **data );
  extern "C" int  EG_getBody( const egObject *topo, egObject **result );
  extern "C" int  EG_inTopology( const egObject *topo, const double *xyz );
  extern "C" int  EG_inFace( const egObject *face, const double *uv );
  extern "C" int  EG_inFaceOCC( const egObject *face, double tol,
                                const double *uv );
  extern "C" int  EG_inFaceAlt( const egObject *face, const double *uv );
  extern "C" int  EG_sewFaces( int nobj, const egObject **objs, double toler,
                               int opt, egObject **result );
  extern "C" int  EG_replaceFaces( const egObject *body,  int nobj,
                                         egObject **objs, egObject **result );
  extern "C" int  EG_mapBody( const egObject *sBody, const egObject *dBody,
                              const char *fAttr, egObject **mBody );
  extern "C" int  EG_copyGeometry_dot( const egObject *obj,
                                       /*@null@*/ const double *xform,
                                       /*@null@*/ const double *xform_dot,
                                       egObject *copy );
  extern "C" int  EG_hasGeometry_dot( const egObject *obj );

  extern     int  EG_attriBodyDup( const egObject *src, egObject *dst );
  extern     void EG_getGeometryLen( const egObject *geom, int *ni, int *nr );
  extern     void EG_completePCurve( egObject *g, Handle(Geom2d_Curve) &hCurv );
  extern     void EG_completeCurve(  egObject *g, Handle(Geom_Curve)   &hCurv );
  extern     void EG_completeSurf(   egObject *g, Handle(Geom_Surface) &hSurf );
  extern     int  EG_addStrAttr( egObject *obj, const char *name,
                                 const char *str );
  extern     int  EG_inFaceX( const egObject *face, const double *uva,
                              /*@null@*/ double *pt, /*@null@*/ double *uvx );

  extern "C" int  EG_makeSolidBox( egObject *context, const double *data,
                                   egObject **body );
  extern "C" int  EG_makeSolidBox_dot( const double *data, const double *data_dot,
                                       egObject *body );
  extern "C" int  EG_makeSolidSphere( egObject *context, int stypx,
                                      const double *data, egObject **body );
  extern "C" int  EG_makeSolidSphere_dot( int stypx, const double *data,
                                          const double *data_dot,
                                          egObject *body );
  extern "C" int  EG_makeSolidCone( egObject *context, int stypx,
                                    const double *data, egObject **body );
  extern "C" int  EG_makeSolidCone_dot( int stypx, const double *data,
                                        const double *data_dot, egObject *body );
  extern "C" int  EG_makeSolidCylinder( egObject *context, int stypx,
                                        const double *data, egObject **body );
  extern "C" int  EG_makeSolidCylinder_dot( int stypx, const double *data,
                                            const double *data_dot,
                                            egObject *body );
  extern "C" int  EG_makeSolidTorus( egObject *context, int stypx,
                                     const double *data, egObject **body );
  extern "C" int  EG_makeSolidTorus_dot( int stypx, const double *data,
                                         const double *data_dot,
                                         egObject *body );


static void
EG_cleanMaps(egadsMap *map)
{
  if (map->objs == NULL) return;
  EG_free(map->objs);
  map->objs = NULL;
}


static void
EG_printStatus(BRepCheck_Status cStatus)
{
  if (cStatus == BRepCheck_NoError) return;

  if (cStatus == BRepCheck_InvalidPointOnCurve) {
    printf(" EGADS Fault: Invalid Point On Curve\n");
  } else if (cStatus == BRepCheck_InvalidPointOnCurveOnSurface) {
    printf(" EGADS Fault: Invalid Point On Curve On Surface\n");
  } else if (cStatus == BRepCheck_InvalidPointOnSurface) {
    printf(" EGADS Fault: Invalid Point On Surface\n");
  } else if (cStatus == BRepCheck_No3DCurve) {
    printf(" EGADS Fault: No 3D Curve\n");
  } else if (cStatus == BRepCheck_Multiple3DCurve) {
    printf(" EGADS Fault: Multiple 3D Curves\n");
  } else if (cStatus == BRepCheck_Invalid3DCurve) {
    printf(" EGADS Fault: Invalid 3D Curve\n");
  } else if (cStatus == BRepCheck_NoCurveOnSurface) {
    printf(" EGADS Fault: No Curve On Surface\n");
  } else if (cStatus == BRepCheck_InvalidCurveOnSurface) {
    printf(" EGADS Fault: Invalid Curve On Surface\n");
  } else if (cStatus == BRepCheck_InvalidCurveOnClosedSurface) {
    printf(" EGADS Fault: Invalid Curve On Closed Surface\n");
  } else if (cStatus == BRepCheck_InvalidSameRangeFlag) {
    printf(" EGADS Fault: Invalid SameRange Flag\n");
  } else if (cStatus == BRepCheck_InvalidSameParameterFlag) {
    printf(" EGADS Fault: Invalid Same Parameter Flag\n");
  } else if (cStatus == BRepCheck_InvalidDegeneratedFlag) {
    printf(" EGADS Fault: Invalid Degenerated Flag\n");
  } else if (cStatus == BRepCheck_FreeEdge) {
    printf(" EGADS Fault: Free Edge\n");
  } else if (cStatus == BRepCheck_InvalidMultiConnexity) {
    printf(" EGADS Fault: Invalid Multi Connexity\n");
  } else if (cStatus == BRepCheck_InvalidRange) {
    printf(" EGADS Fault: Invalid Range\n");
  } else if (cStatus == BRepCheck_EmptyWire) {
    printf(" EGADS Fault: Empty Wire\n");
  } else if (cStatus == BRepCheck_RedundantEdge) {
    printf(" EGADS Fault: Redundant Edge\n");
  } else if (cStatus == BRepCheck_SelfIntersectingWire) {
    printf(" EGADS Fault: Self Intersecting Wire\n");
  } else if (cStatus == BRepCheck_NoSurface) {
    printf(" EGADS Fault: No Surface\n");
  } else if (cStatus == BRepCheck_InvalidWire) {
    printf(" EGADS Fault: Invalid Wire\n");
  } else if (cStatus == BRepCheck_RedundantWire) {
    printf(" EGADS Fault: Redundant Wire\n");
  } else if (cStatus == BRepCheck_IntersectingWires) {
    printf(" EGADS Fault: Intersecting Wires\n");
  } else if (cStatus == BRepCheck_InvalidImbricationOfWires) {
    printf(" EGADS Fault: Invalid Imbrication Of Wires\n");
  } else if (cStatus == BRepCheck_EmptyShell) {
    printf(" EGADS Fault: Empty Shell\n");
  } else if (cStatus == BRepCheck_RedundantFace) {
    printf(" EGADS Fault: Redundant Face\n");
  } else if (cStatus == BRepCheck_UnorientableShape) {
    printf(" EGADS Fault: Unorientable Shape\n");
  } else if (cStatus == BRepCheck_NotClosed) {
    printf(" EGADS Fault: Not Closed\n");
  } else if (cStatus == BRepCheck_NotConnected) {
    printf(" EGADS Fault: Not Connected\n");
  } else if (cStatus == BRepCheck_SubshapeNotInShape) {
    printf(" EGADS Fault: Subshape Not In Shape\n");
  } else if (cStatus == BRepCheck_BadOrientation) {
    printf(" EGADS Fault: Bad Orientation\n");
  } else if (cStatus == BRepCheck_BadOrientationOfSubshape) {
    printf(" EGADS Fault: Bad Orientation Of Subshape\n");
  } else if (cStatus == BRepCheck_InvalidPolygonOnTriangulation) {
    printf(" EGADS Fault: Invalid Polygon On Triangulation\n");
  } else if (cStatus == BRepCheck_InvalidToleranceValue) {
    printf(" EGADS Fault: Invalid Tolerance Value\n");
  } else if (cStatus == BRepCheck_CheckFail) {
    printf(" EGADS Fault: Check Fail\n");
  } else {
    printf(" EGADS Unknown Fault = %d\n", cStatus);
  }
}


void
EG_checkStatus(const Handle_BRepCheck_Result tResult)
{

  if (tResult.IsNull()) return;
  const BRepCheck_ListOfStatus& tList = tResult->Status();
  if (tList.Extent() <= 0) return;

  BRepCheck_Status cStatus = tList.First();
  EG_printStatus(cStatus);
}


int
EG_destroyTopology(egObject *topo)
{
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;

  if (topo->blind == NULL) return EGADS_SUCCESS;

  if (topo->oclass == MODEL) {
    egadsModel *mshape = (egadsModel *) topo->blind;
    if (mshape->bodies != NULL) {
      for (int i = 0; i < mshape->nbody; i++)
        EG_dereferenceObject(mshape->bodies[i], topo);
      delete [] mshape->bodies;
    }
    mshape->shape.Nullify();
    delete mshape;

  } else if (topo->oclass == BODY) {

    egadsBody *pbody = (egadsBody *) topo->blind;
    if (pbody != NULL) {
      if (topo->mtype == WIREBODY) {
        int nwire = pbody->loops.map.Extent();
        for (int i = 0; i < nwire; i++)
          EG_dereferenceObject(pbody->loops.objs[i], topo);
      } else if (topo->mtype == FACEBODY) {
        int nface = pbody->faces.map.Extent();
        for (int i = 0; i < nface; i++)
          EG_dereferenceObject(pbody->faces.objs[i], topo);
      } else {
        int nshell = pbody->shells.map.Extent();
        for (int i = 0; i < nshell; i++)
          EG_dereferenceObject(pbody->shells.objs[i], topo);
        if (topo->mtype == SOLIDBODY) delete [] pbody->senses;
      }
      EG_cleanMaps(&pbody->shells);
      EG_cleanMaps(&pbody->faces);
      EG_cleanMaps(&pbody->loops);
      EG_cleanMaps(&pbody->edges);
      EG_cleanMaps(&pbody->nodes);
      delete pbody;
    }

  } else if (topo->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    if (pshell != NULL) {
      if (pshell->topFlg == 0) {
        for (int i = 0; i < pshell->nfaces; i++)
          EG_dereferenceObject(pshell->faces[i], topo);
      } else {
        for (int i = 0; i < pshell->nfaces; i++)
          EG_dereferenceTopObj(pshell->faces[i], topo);
      }
      delete [] pshell->faces;
      delete pshell;
    }

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    if (pface != NULL) {
      if (pface->topFlg == 0) {
        for (int i = 0; i < pface->nloops; i++)
          EG_dereferenceObject(pface->loops[i], topo);
        EG_dereferenceObject(pface->surface, topo);
      } else {
        for (int i = 0; i < pface->nloops; i++)
          EG_dereferenceTopObj(pface->loops[i], topo);
        EG_dereferenceTopObj(pface->surface, topo);
      }
      delete [] pface->senses;
      delete [] pface->loops;
      delete pface;
    }

  } else if (topo->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) topo->blind;
    if (ploop != NULL) {
      if (ploop->topFlg == 0) {
        for (int i = 0; i < ploop->nedges; i++) {
          EG_dereferenceObject(ploop->edges[i], topo);
          if (ploop->surface != NULL)
            EG_dereferenceObject(ploop->edges[i+ploop->nedges], topo);
        }
        if (ploop->surface != NULL)
          EG_dereferenceObject(ploop->surface, topo);
      } else {
        for (int i = 0; i < ploop->nedges; i++) {
          EG_dereferenceTopObj(ploop->edges[i], topo);
          if (ploop->surface != NULL)
            EG_dereferenceTopObj(ploop->edges[i+ploop->nedges], topo);
        }
        if (ploop->surface != NULL)
          EG_dereferenceTopObj(ploop->surface, topo);
      }
      delete [] ploop->senses;
      delete [] ploop->edges;
      delete ploop;
    }

  } else if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    if (pedge != NULL) {
      int degen = 0;

      if ((pedge->curve == NULL) &&
          (topo->mtype  == DEGENERATE)) degen = 1;
      if (pedge->topFlg == 0) {
        if (degen == 0) EG_dereferenceObject(pedge->curve, topo);
        EG_dereferenceObject(pedge->nodes[0], topo);
        EG_dereferenceObject(pedge->nodes[1], topo);
      } else {
        if (degen == 0) EG_dereferenceTopObj(pedge->curve, topo);
        EG_dereferenceTopObj(pedge->nodes[0], topo);
        EG_dereferenceTopObj(pedge->nodes[1], topo);
      }
      delete pedge;
    }

  } else {

    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode != NULL) delete pnode;

  }

  return EGADS_SUCCESS;
}


void
EG_splitPeriodics(egadsBody *body)
{
  int          hit    = 0;
  TopoDS_Shape bshape = body->shape;
  TopTools_IndexedMapOfShape MapE;
  TopExp::MapShapes(bshape, TopAbs_EDGE, MapE);
  for (int i = 1; i <= MapE.Extent(); i++) {
    TopoDS_Shape shape = MapE(i);
    TopoDS_Edge  Edge  = TopoDS::Edge(shape);
    if (Edge.Closed()) hit++;
  }
  if (hit == 0) {
    TopTools_IndexedMapOfShape MapF;
    TopExp::MapShapes(bshape, TopAbs_FACE, MapF);
    for (int i = 1; i <= MapF.Extent(); i++) {
      TopoDS_Shape shape = MapF(i);
      TopoDS_Face  Face  = TopoDS::Face(shape);
      BRepAdaptor_Surface aSurf(Face, Standard_True);
      if (aSurf.IsUClosed()) hit++;
      if (aSurf.IsVClosed()) hit++;
    }
  }
  if (hit == 0) return;

  // use the OpenCASCADE method ->

  TopoDS_Shape solid = bshape;
  Handle(ShapeBuild_ReShape) reShape = new ShapeBuild_ReShape();
  ShapeUpgrade_ShapeDivideClosed aShape(bshape);
  aShape.SetNbSplitPoints(1);
  aShape.SetContext(reShape);
  if (aShape.Perform(Standard_False)) {
    solid = reShape->Apply(bshape);
    if (solid.IsNull()) {
      printf(" EGADS Warning: Can't Split Periodics!\n");
      solid = bshape;
    } else {
      BRepCheck_Analyzer fCheck(solid);
      if (!fCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(solid);
        sfs->Perform();
        TopoDS_Shape fixedSolid = sfs->Shape();
        if (fixedSolid.IsNull()) {
          printf(" EGADS Warning: Periodic Split is Invalid!\n");
          solid = bshape;
        } else {
          BRepCheck_Analyzer sfCheck(fixedSolid);
          if (!sfCheck.IsValid()) {
            printf(" EGADS Warning: Periodic Split is InValid!\n");
            solid = bshape;
          } else {
            solid = fixedSolid;
          }
        }
      }
    }
  }

  body->shape = solid;
}


void
EG_splitMultiplicity(egadsBody *body, int outLevel)
{
  TopoDS_Shape bshape = body->shape;

  int hite = 0;
  TopTools_IndexedMapOfShape MapE;
  TopExp::MapShapes(bshape, TopAbs_EDGE, MapE);
  for (int i = 1; i <= MapE.Extent(); i++) {
    Standard_Real t1, t2;
    TopoDS_Shape shape = MapE(i);
    TopoDS_Edge  Edge  = TopoDS::Edge(shape);
    Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
    Handle(Geom_BSplineCurve) hBSpline =
                                    Handle(Geom_BSplineCurve)::DownCast(hCurve);
    if  (hBSpline.IsNull()) continue;
    if  (hBSpline->Continuity() == GeomAbs_C0)  hite++;
/*  if ((hBSpline->Degree() > 1) &&
        (hBSpline->Continuity() == GeomAbs_C1)) hite++;  */
/*  int c = -1;
    if (hBSpline->Continuity() == GeomAbs_C0) c = 0;
    if (hBSpline->Continuity() == GeomAbs_C1) c = 1;
    if (hBSpline->Continuity() == GeomAbs_C2) c = 2;
    if (hBSpline->Continuity() == GeomAbs_C3) c = 3;
    if (hBSpline->Continuity() == GeomAbs_CN) c = 99;
    printf(" Edge %d: degree = %d, c = %d!\n", i, hBSpline->Degree(), c);  */
  }
  int hitf = 0;
  TopTools_IndexedMapOfShape MapF;
  TopExp::MapShapes(bshape, TopAbs_FACE, MapF);
  for (int i = 1; i <= MapF.Extent(); i++) {
    TopoDS_Shape shape = MapF(i);
    TopoDS_Face  Face  = TopoDS::Face(shape);
    Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
    Handle(Geom_BSplineSurface) hBSpline =
                                Handle(Geom_BSplineSurface)::DownCast(hSurface);
    if   (hBSpline.IsNull()) continue;
    if   (hBSpline->Continuity() == GeomAbs_C0)  hitf++;
/*  if (((hBSpline->UDegree() > 1) || (hBSpline->VDegree() > 1)) &&
         (hBSpline->Continuity() == GeomAbs_C1)) hitf++;  */
  }
  if (hite+hitf == 0) return;

  // use the OpenCASCADE method ->

  TopoDS_Shape shape = bshape;
  ShapeUpgrade_ShapeDivideContinuity aShape(bshape);
/*  aShape.SetBoundaryCriterion(GeomAbs_C2);
    aShape.SetPCurveCriterion(GeomAbs_C2);
    aShape.SetSurfaceCriterion(GeomAbs_C2);  */
  if (aShape.Perform()) {
    if (aShape.Status(ShapeExtend_OK)) {
      printf(" EGADS Warning: No Splitting %d %d (EG_splitMultiplicity)!\n",
             hite, hitf);
      return;
    }
    if (outLevel > 1) {
      if (aShape.Status(ShapeExtend_DONE1))
        printf(" EGADS Info: Some Edges Split!\n");
      if (aShape.Status(ShapeExtend_DONE2))
        printf(" EGADS Info: Surfaces were Split!\n");
      if (aShape.Status(ShapeExtend_DONE3))
        printf(" EGADS Info: Surfaces were modified without splitting!\n");
      if (aShape.Status(ShapeExtend_FAIL1))
        printf(" EGADS Info: Some errors occured splitting Wires!\n");
      if (aShape.Status(ShapeExtend_FAIL2))
        printf(" EGADS Info: Faces could not be split!\n");
    }
    if (!aShape.Status(ShapeExtend_DONE)) {
      printf(" EGADS Warning: Not Done (EG_splitMultiplicity)!\n");
      return;
    }
    shape = aShape.Result();
    if (shape.IsNull()) {
      printf(" EGADS Warning: Can't Do Continuity Split!\n");
      shape = bshape;
    } else {
      BRepCheck_Analyzer fCheck(shape);
      if (!fCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          printf(" EGADS Warning: Continuity Split is Invalid!\n");
          shape = bshape;
        } else {
          BRepCheck_Analyzer sfCheck(fixedShape);
          if (!sfCheck.IsValid()) {
            printf(" EGADS Warning: Continuity Split is InValid!\n");
            shape = bshape;
          } else {
            shape = fixedShape;
          }
        }
      }
    }
  } else {
    printf(" EGADS Warning: Perform is False (EG_splitMultiplicity)!\n");
    return;
  }

  body->shape = shape;
}


static int
EG_bodyRecurseGeom_dot(const egObject *obj, egadsBody *pbody)
{
  int stat;
  int i;

  if (obj->oclass == NODE) {

    egadsNode *pnode = (egadsNode *) obj->blind;
    int index = pbody->nodes.map.FindIndex(pnode->node);
    if (index == 0) return EGADS_TOPOERR;
    stat = EG_copyGeometry_dot(obj, NULL, NULL, pbody->nodes.objs[index-1]);
    if (stat != EGADS_SUCCESS) return stat;

  } else if (obj->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) obj->blind;
    int index = pbody->edges.map.FindIndex(pedge->edge);
    if (index != 0) {
      egadsEdge *bedge = (egadsEdge *) pbody->edges.objs[index-1]->blind;
       if (obj->mtype == DEGENERATE) {
        if (pbody->edges.objs[index-1]->mtype != DEGENERATE)
          return EGADS_TOPOERR;
      } else {
        /* copy curve to eliminate redundant copying... */
        stat = EG_copyGeometry_dot(pedge->curve, NULL, NULL, bedge->curve);
        if (stat != EGADS_SUCCESS) return stat;
      }
      /* copy t-range sensitivity */
      bedge->filled = pedge->filled;
      for (i = 0; i < 2; i++) bedge->trange_dot[i] = pedge->trange_dot[i];
    }
    stat = EG_bodyRecurseGeom_dot(pedge->nodes[0], pbody);
    if (stat != EGADS_SUCCESS) return stat;
    if (obj->mtype == TWONODE) {
      stat = EG_bodyRecurseGeom_dot(pedge->nodes[1], pbody);
      if (stat != EGADS_SUCCESS) return stat;
    }

  } else if (obj->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) obj->blind;
    for (int i = 0; i < ploop->nedges; i++) {
      stat = EG_bodyRecurseGeom_dot(ploop->edges[i], pbody);
      if (stat != EGADS_SUCCESS) return stat;
    }

  } else if (obj->oclass == FACE) {

    egadsFace *pface = (egadsFace *) obj->blind;
    int index = pbody->faces.map.FindIndex(pface->face);
    if (index != 0) {
      /* copy surface to eliminate redundant copying... */
      egadsFace *bface = (egadsFace *) pbody->faces.objs[index-1]->blind;
      stat = EG_copyGeometry_dot(pface->surface, NULL, NULL, bface->surface);
      if (stat != EGADS_SUCCESS) return stat;
    }
    for (int i = 0; i < pface->nloops; i++) {
      stat = EG_bodyRecurseGeom_dot(pface->loops[i], pbody);
      if (stat != EGADS_SUCCESS) return stat;
    }

  } else if (obj->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) obj->blind;
    for (int i = 0; i < pshell->nfaces; i++) {
      stat = EG_bodyRecurseGeom_dot(pshell->faces[i], pbody);
      if (stat != EGADS_SUCCESS) return stat;
    }

  }

  return EGADS_SUCCESS;
}


static int
EG_bodyTravGeom_dot(const egObject *obj, egadsBody *pbody)
{
  int stat;

  if (obj->blind == NULL)     return EGADS_NULLOBJ;
  stat = EG_hasGeometry_dot(obj);
  if (stat == EGADS_NOTFOUND) return EGADS_SUCCESS;
  if (stat != EGADS_SUCCESS)  return stat;

  return EG_bodyRecurseGeom_dot(obj, pbody);
}


void
EG_fillPCurves(TopoDS_Face face, egObject *surfo, egObject *loopo,
                                 egObject *topObj)
{
  int           i = 0;
  egObject      *geom;
  Standard_Real f, l;

  egadsLoop *ploop = (egadsLoop *) loopo->blind;
  if (ploop->surface == NULL) return;
  if (ploop->surface != surfo) {
    printf(" EGADS Internal: Loop/Face mismatch on Surface!\n");
    return;
  }

  TopoDS_Wire wire = ploop->loop;
  BRepTools_WireExplorer ExpWE;
  for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
    if (ploop->edges[ploop->nedges+i] != NULL) {
      printf(" EGADS Internal: PCurve already Filled!\n");
      return;
    }
    TopoDS_Shape shape = ExpWE.Current();
    TopoDS_Edge  edge  = TopoDS::Edge(shape);
    if (EG_makeObject(EG_context(surfo),
                      &ploop->edges[ploop->nedges+i]) == EGADS_SUCCESS) {
      geom = ploop->edges[ploop->nedges+i];
      Handle(Geom2d_Curve) hCurve = BRep_Tool::
                                    CurveOnSurface(edge, face, f, l);
      if (hCurve.IsNull()) {
        Handle(Geom_Curve)   hCurv    = BRep_Tool::Curve(edge, f, l);
        Handle(Geom_Surface) hSurface = BRep_Tool::Surface(face);
        double toler = BRep_Tool::Tolerance(edge);
        try {
          hCurve = GeomProjLib::Curve2d(hCurv, f, l, hSurface, toler);
        }
        catch (const Standard_Failure& e) {
          printf(" EGADS Info: Geometry Creation Error (EG_fillPCurves)!\n");
        }
        catch (...) {
          printf(" EGADS Info: Geometry Creation Error (EG_fillPCurves)!\n");
        }
      }
      geom->topObj = topObj;
      EG_completePCurve(geom,  hCurve);
      EG_referenceObject(geom, loopo);
    }
    i++;
  }

}


int
EG_shellClosure(egadsShell *pshell, int mtype)
{
  int ret, i, *hits;

  TopoDS_Shell Shell = pshell->shell;
  TopTools_IndexedMapOfShape MapE;
  TopExp::MapShapes(Shell, TopAbs_EDGE, MapE);
  if (MapE.Extent() == 0) return CLOSED;

  hits = new int[MapE.Extent()];
  for (i = 0; i < MapE.Extent(); i++) hits[i] = 0;

  TopExp_Explorer ExpW;
  for (ExpW.Init(Shell, TopAbs_EDGE); ExpW.More(); ExpW.Next()) {
    TopoDS_Shape shape = ExpW.Current();
    TopoDS_Edge  edge  = TopoDS::Edge(shape);
    if (BRep_Tool::Degenerated(edge)) continue;
    i = MapE.FindIndex(edge);
    if (i == 0) {
      printf(" EGADS Internal: Edge not found (EG_shellClosure)!\n");
      continue;
    }
    hits[i-1]++;
  }

  ret = CLOSED;
  for (i = 0; i < MapE.Extent(); i++)
    if ((hits[i] != 2) && (hits[i] != 0)) ret = OPEN;
  if ((mtype == DEGENERATE) && (ret == OPEN))
    for (i = 0; i < MapE.Extent(); i++)
      printf(" EGADS Info: Edge %d: hits = %d\n", i+1, hits[i]);

  delete [] hits;
  return ret;
}


void
EG_checkLoops(egObject *loop)
{
  int          j, n, stat, atype, alen;
  double       t, d, x0[2], x1[2], data[6];
  const int    *ints;
  const double *reals;
  const char   *str;

  if (loop == NULL)               return;
  if (loop->magicnumber != MAGIC) return;
  if (loop->blind == NULL)        return;
  if (loop->oclass != LOOP)       return;

  egadsLoop *lloop = (egadsLoop *) loop->blind;
  if (lloop->surface == NULL)     return;
  for (n = 0; n < lloop->nedges; n++) {
    if (lloop->edges[n] == NULL) continue;
    egadsEdge *ledge = (egadsEdge *) lloop->edges[n]->blind;
    egObject  *pcobj = lloop->edges[lloop->nedges+n];
    stat  = EG_attributeRet(pcobj, ".Bad", &atype, &alen,
                            &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) continue;
    EG_evaluate(pcobj, &ledge->trange[0], data);
    d     = sqrt(data[2]*data[2] + data[3]*data[3]);
    x0[0] = x0[1] = 0.0;
    if (d != 0.0) {
      x0[0] = data[2]/d;
      x0[1] = data[3]/d;
    }
    for (j = 1; j < 1000; j++) {
      t = ledge->trange[0]+j*(ledge->trange[1]-ledge->trange[0])/999.;
      EG_evaluate(pcobj, &t, data);
      d     = sqrt(data[2]*data[2] + data[3]*data[3]);
      x1[0] = x1[1] = 0.0;
      if (d != 0.0) {
        x1[0] = data[2]/d;
        x1[1] = data[3]/d;
      }
      if (x0[0]*x1[0] + x0[1]*x1[1] < -0.95) {
        stat = EG_addStrAttr(pcobj, ".Bad", "fold");
        if (stat != EGADS_SUCCESS)
          printf(" EGADS Info: EG_addStrAttr fold= %d\n", stat);
      }
      x0[0] = x1[0];
      x0[1] = x1[1];
    }
  }
}


void
EG_fillTopoObjs(egObject *object, egObject *topObj)
{
  int           outLevel, stat, degen = 0;
  egObject      *context;
  Standard_Real t1, t2;
  TopoDS_Vertex V1, V2;

  outLevel = EG_outLevel(object);
  context  = EG_context(object);

  if (object->oclass == EDGE) {

    egObject *geom   = NULL;
    egObject *pn1    = NULL;
    egObject *pn2    = NULL;
    egadsEdge *pedge = (egadsEdge *) object->blind;
    TopoDS_Edge Edge = pedge->edge;
    if (BRep_Tool::Degenerated(Edge)) {
      degen = 1;
    } else {
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
      stat = EG_makeObject(context, &geom);
      if (stat == EGADS_SUCCESS) {
        geom->topObj = topObj;
        EG_completeCurve(geom, hCurve);
      }
    }

    TopExp::Vertices(Edge, V2, V1, Standard_True);
    EG_makeObject(context, &pn1);
    if (pn1 != NULL) {
      egadsNode *pnode   = new egadsNode;
      gp_Pnt pv          = BRep_Tool::Pnt(V1);
      pnode->node        = V1;
      pnode->xyz[0]      = pv.X();
      pnode->xyz[1]      = pv.Y();
      pnode->xyz[2]      = pv.Z();
      pnode->bbox.filled = 0;
      pnode->filled      = 0;
      pnode->xyz_dot[0]  = pnode->xyz_dot[1] = pnode->xyz_dot[2] = 0.0;
      pn1->oclass        = NODE;
      pn1->blind         = pnode;
      pn1->topObj        = topObj;
      BRepCheck_Analyzer v1Check(V1);
      if(!v1Check.IsValid())
        if (outLevel > 0)
        printf(" EGADS Info: Node1 may be invalid (EG_fillTopoObjs)!\n");
    }
    if (V1.IsSame(V2)) {
      object->mtype = ONENODE;
      pn2           = pn1;
    } else {
      object->mtype = TWONODE;
      EG_makeObject(context, &pn2);
      if (pn2 != NULL) {
        egadsNode *pnode   = new egadsNode;
        gp_Pnt pv          = BRep_Tool::Pnt(V2);
        pnode->node        = V2;
        pnode->xyz[0]      = pv.X();
        pnode->xyz[1]      = pv.Y();
        pnode->xyz[2]      = pv.Z();
        pnode->bbox.filled = 0;
        pnode->filled      = 0;
        pnode->xyz_dot[0]  = pnode->xyz_dot[1] = pnode->xyz_dot[2] = 0.0;
        pn2->oclass        = NODE;
        pn2->blind         = pnode;
        pn2->topObj        = topObj;
        BRepCheck_Analyzer v2Check(V2);
        if(!v2Check.IsValid())
          if (outLevel > 0)
          printf(" EGADS Info: Node2 may be invalid (EG_fillTopoObjs)!\n");
      }
    }
    if (Edge.Orientation() != TopAbs_REVERSED) {
      pedge->nodes[0] = pn2;
      pedge->nodes[1] = pn1;
    } else {
      pedge->nodes[0] = pn1;
      pedge->nodes[1] = pn2;
    }

    pedge->curve       = geom;
    pedge->topFlg      = 0;
    pedge->bbox.filled = 0;
    object->topObj     = topObj;
    if (degen == 1) {
      object->mtype    = DEGENERATE;
    } else {
      EG_referenceObject(geom, object);
    }
    EG_referenceObject(pn1,  object);
    EG_referenceObject(pn2,  object);
    BRepCheck_Analyzer eCheck(Edge);
    if (!eCheck.IsValid())
      if (outLevel > 0)
        printf(" EGADS Info: Edge may be invalid (EG_fillTopoObjs)!\n");

  } else if (object->oclass == LOOP) {

    int      *senses = NULL;
    egObject **edgeo = NULL;
    egadsLoop *ploop = (egadsLoop *) object->blind;
    TopoDS_Wire Wire = ploop->loop;
    int            n = 1;
    int           ne = 0;
    int       closed = 0;
    if (Wire.Closed()) closed = 1;
    // more reliable for checking closure of Wires
    TopExp::Vertices(Wire, V1, V2);
    if (!V1.IsNull() && !V2.IsNull())
      if (V1.IsSame(V2)) {
        closed = 1;
      } else {
        closed = 0;
      }
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) ne++;
    if (ploop->surface != NULL) n = 2;

    if (ne > 0) {
      edgeo  = new egObject*[n*ne];
      senses = new int[ne];
    }
    int k = 0;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) {
      if (edgeo == NULL) continue;
      TopoDS_Shape shapW = ExpWE.Current();
      TopoDS_Edge  Edge  = TopoDS::Edge(shapW);
      edgeo[k]           = NULL;
      senses[k]          = 1;
      if (n == 2) edgeo[k+ne] = NULL;
      if (shapW.Orientation() == TopAbs_REVERSED) senses[k] = -1;
      for (int j = 0; j < k; j++) {
        egadsEdge *pedg = (egadsEdge *) edgeo[j]->blind;
        if (Edge.IsSame(pedg->edge)) {
          edgeo[k] = edgeo[j];
          break;
        }
      }
      if (edgeo[k] == NULL) {
        stat = EG_makeObject(context, &edgeo[k]);
        if (stat != EGADS_SUCCESS) continue;
        edgeo[k]->oclass = EDGE;
        egadsEdge *pedge   = new egadsEdge;
        pedge->edge        = Edge;
        pedge->curve       = NULL;
        pedge->nodes[0]    = NULL;
        pedge->nodes[1]    = NULL;
        pedge->bbox.filled = 0;
        BRep_Tool::Range(Edge, pedge->trange[0], pedge->trange[1]);
        pedge->filled      = 0;
        pedge->trange_dot[0] = pedge->trange_dot[1] = 0;
        edgeo[k]->blind  = pedge;
        EG_fillTopoObjs(edgeo[k], topObj);
      }
      EG_referenceObject(edgeo[k], object);
      k++;
    }

    ploop->nedges      = ne;
    ploop->edges       = edgeo;
    ploop->senses      = senses;
    ploop->topFlg      = 0;
    ploop->bbox.filled = 0;
    object->topObj     = topObj;
    object->mtype      = OPEN;
    if (closed == 1) object->mtype = CLOSED;
    BRepCheck_Analyzer wCheck(Wire);
    if (!wCheck.IsValid())
      if (outLevel > 0)
        printf(" EGADS Info: Loop may be invalid (EG_fillTopoObjs)!\n");

  } else {

    printf(" EGADS Internal: Not Implemented (EG_fillTopoObjs)!\n");

  }

}


static void
EG_MapShapes(int outLevel, TopoDS_Shape& shape, TopAbs_ShapeEnum eshape,
             TopTools_IndexedMapOfShape& mapshape, const char *sshape)
{
  mapshape.Clear();
  TopTools_IndexedMapOfShape map;
  TopExp::MapShapes(shape, eshape, map);
  for (int i = 0; i < map.Extent(); i++) {
    if ((map(i+1).Orientation() == TopAbs_INTERNAL) ||
        (map(i+1).Orientation() == TopAbs_EXTERNAL)) continue;
    mapshape.Add(map(i+1));
  }
  if ((outLevel > 0) &&
      (mapshape.Extent() != map.Extent())) {
    printf(" EGADS Info: Ignoring %d internal/external %ss\n",
           map.Extent()-mapshape.Extent(), sshape);
  }
}


int
EG_traverseBody(egObject *context, int i, egObject *bobj,
                egObject *topObj, egadsBody *body, int *nerr)
{
  int          ii, j, k, outLevel, stat, hit = 0, solid = 0;
  TopoDS_Shape shape;
  egObject     *obj, *geom, **curves_obj = NULL, **surfs_obj = NULL;
  TColStd_IndexedMapOfTransient curves_map, surfs_map;
  NCollection_Vector<Handle(Geom_Curve)>   curves_vec;
  NCollection_Vector<Handle(Geom_Surface)> surfs_vec;

  *nerr    = 0;
  outLevel = EG_outLevel(context);
  if (body->shape.ShapeType() == TopAbs_SOLID) solid = 1;

  TopExp_Explorer Exp;
  EG_MapShapes(outLevel, body->shape, TopAbs_VERTEX, body->nodes.map,  "NODE" );
  EG_MapShapes(outLevel, body->shape, TopAbs_EDGE,   body->edges.map,  "EDGE" );
  EG_MapShapes(outLevel, body->shape, TopAbs_WIRE,   body->loops.map,  "LOOP" );
  EG_MapShapes(outLevel, body->shape, TopAbs_FACE,   body->faces.map,  "FACE" );
  EG_MapShapes(outLevel, body->shape, TopAbs_SHELL,  body->shells.map, "SHELL");
  int nNode  = body->nodes.map.Extent();
  int nEdge  = body->edges.map.Extent();
  int nLoop  = body->loops.map.Extent();
  int nFace  = body->faces.map.Extent();
  int nShell = body->shells.map.Extent();
  bobj->oclass = BODY;
  bobj->mtype  = WIREBODY;
  if (nFace > 0) {
    bobj->mtype = FACEBODY;
    if (nShell > 0) {
      bobj->mtype = SHEETBODY;
      if (solid == 1) bobj->mtype = SOLIDBODY;
    }
  }

  if (outLevel > 1)
    printf(" EGADS Info: Shape %d has %d Nodes, %d Edges, %d Loops, %d Faces and %d Shells\n",
           i+1, nNode, nEdge, nLoop, nFace, nShell);

  // check for sensical entity counts

  if (nNode == 0) {
    hit++;
    if (outLevel > 0)
      printf(" EGADS Warning: Shape %d has zero Nodes!\n", i+1);
  }
  if (nEdge == 0) {
    hit++;
    if (outLevel > 0)
      printf(" EGADS Warning: Shape %d has zero Edges!\n", i+1);
  }
  if (nLoop == 0) {
    if (nFace > 0) {
      hit++;
      if (outLevel > 0)
        printf(" EGADS Warning: Shape %d has zero Loops and %d Faces!\n",
               i+1, nFace);
    }
    if (nShell > 0) {
      hit++;
      if (outLevel > 0)
        printf(" EGADS Warning: Shape %d has zero Loops and %d Shells!\n",
               i+1, nShell);
    }
  }
  if (nFace == 0) {
    if (nShell > 0) {
      hit++;
      if (outLevel > 0)
        printf(" EGADS Warning: Shape %d has zero Faces and %d Shells!\n",
               i+1, nShell);
    }
  }

  if (hit > 0) {
    *nerr = hit;
    return EGADS_TOPOERR;
  }

  // allocate ego storage

  if (nNode > 0) {
    body->nodes.objs = (egObject **) EG_alloc(nNode*sizeof(egObject *));
    if (body->nodes.objs == NULL) return EGADS_MALLOC;
    for (j = 0; j < nNode; j++) {
      stat = EG_makeObject(context, &body->nodes.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }
  if (nEdge > 0) {
    body->edges.objs = (egObject **) EG_alloc(nEdge*sizeof(egObject *));
    if (body->edges.objs == NULL) {
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nEdge; j++) {
      stat = EG_makeObject(context, &body->edges.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }
  if (nLoop > 0) {
    body->loops.objs = (egObject **) EG_alloc(nLoop*sizeof(egObject *));
    if (body->loops.objs == NULL) {
      EG_cleanMaps(&body->edges);
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nLoop; j++) {
      stat = EG_makeObject(context, &body->loops.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->loops);
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }
  if (nFace > 0) {
    body->faces.objs = (egObject **) EG_alloc(nFace*sizeof(egObject *));
    if (body->faces.objs == NULL) {
      EG_cleanMaps(&body->loops);
      EG_cleanMaps(&body->edges);
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nFace; j++) {
      stat = EG_makeObject(context, &body->faces.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->faces);
        EG_cleanMaps(&body->loops);
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }
  if (nShell > 0) {
    body->shells.objs = (egObject **) EG_alloc(nShell*sizeof(egObject *));
    if (body->shells.objs == NULL) {
      EG_cleanMaps(&body->faces);
      EG_cleanMaps(&body->loops);
      EG_cleanMaps(&body->edges);
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nShell; j++) {
      stat = EG_makeObject(context, &body->shells.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->shells);
        EG_cleanMaps(&body->faces);
        EG_cleanMaps(&body->loops);
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }

  if (nEdge > 0) {
    Standard_Real t1, t2;

    /* count all unique curves */
    for (j = 0; j < nEdge; j++) {
      shape            = body->edges.map(j+1);
      TopoDS_Edge Edge = TopoDS::Edge(shape);
      Handle(Geom_Curve) hCurve;
      if (BRep_Tool::Degenerated(Edge)) {
        gp_Ax1 A;
        hCurve = new Geom_Line(A); /* dummy curve for indexing into the map */
      } else {
        hCurve = BRep_Tool::Curve(Edge, t1, t2);
      }
      curves_vec.Append(hCurve);
      curves_map.Add(hCurve);
    }
    int nCurv = curves_map.Extent();

    /* allocate curve egos*/
    curves_obj = (egObject **) EG_alloc(nCurv*sizeof(egObject *));
    if (curves_obj == NULL) {
      EG_cleanMaps(&body->shells);
      EG_cleanMaps(&body->faces);
      EG_cleanMaps(&body->loops);
      EG_cleanMaps(&body->edges);
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }

    /* initialize the curve egos */
    for (j = 0; j < nCurv; j++) {
      stat = EG_makeObject(context, &curves_obj[j]);
      if (stat != EGADS_SUCCESS) {
        EG_free(curves_obj);
        EG_cleanMaps(&body->shells);
        EG_cleanMaps(&body->faces);
        EG_cleanMaps(&body->loops);
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }

  if (nFace > 0) {

    /* generate all surface ego's as they are needed for both loops and faces */
    for (j = 0; j < nFace; j++) {
      shape            = body->faces.map(j+1);
      TopoDS_Face Face = TopoDS::Face(shape);
      Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
      surfs_vec.Append(hSurface);
      surfs_map.Add(hSurface);
    }
    int nSurf = surfs_map.Extent();

    /* allocate surface ego's */
    surfs_obj = (egObject **) EG_alloc(nSurf*sizeof(egObject *));
    if (surfs_obj == NULL) {
      EG_free(curves_obj);
      EG_cleanMaps(&body->shells);
      EG_cleanMaps(&body->faces);
      EG_cleanMaps(&body->loops);
      EG_cleanMaps(&body->edges);
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }

    /* construct the surface ego's */
    for (j = 0; j < nSurf; j++) {
      stat = EG_makeObject(context, &surfs_obj[j]);
      if (stat != EGADS_SUCCESS) {
        EG_free(surfs_obj);
        EG_free(curves_obj);
        EG_cleanMaps(&body->shells);
        EG_cleanMaps(&body->faces);
        EG_cleanMaps(&body->loops);
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
      Handle(Geom_Surface) hSurface = Handle(Geom_Surface)::DownCast(surfs_map(j+1));
      surfs_obj[j]->topObj = topObj;
      EG_completeSurf(surfs_obj[j], hSurface);
    }
  }

  // fill our stuff

  for (j = 0; j < nNode; j++) {
    egadsNode *pnode   = new egadsNode;
    obj                = body->nodes.objs[j];
    shape              = body->nodes.map(j+1);
    TopoDS_Vertex Vert = TopoDS::Vertex(shape);
    gp_Pnt pv          = BRep_Tool::Pnt(Vert);
    pnode->node        = Vert;
    pnode->xyz[0]      = pv.X();
    pnode->xyz[1]      = pv.Y();
    pnode->xyz[2]      = pv.Z();
    pnode->bbox.filled = 0;
    pnode->filled      = 0;
    pnode->xyz_dot[0]  = pnode->xyz_dot[1] = pnode->xyz_dot[2] = 0.0;
    obj->oclass        = NODE;
    obj->blind         = pnode;
    obj->topObj        = topObj;
  }

  for (j = 0; j < nEdge; j++) {
    int           n1, n2, degen = 0;
    TopoDS_Vertex V1, V2;
    Standard_Real t1=0, t2=0;

    egadsEdge *pedge = new egadsEdge;
    obj              = body->edges.objs[j];
    shape            = body->edges.map(j+1);
    geom             = curves_obj[curves_map.FindIndex(curves_vec(j))-1];
    TopoDS_Edge Edge = TopoDS::Edge(shape);
    if (geom->oclass == NIL) {
      geom->topObj = topObj;
      if (BRep_Tool::Degenerated(Edge)) {
        degen        = 1;
        geom->oclass = CURVE;
        geom->mtype  = DEGENERATE;
        geom->blind  = NULL;
      } else {
        BRep_Tool::Range(Edge, t1, t2); /* get parameters */
        Handle(Geom_Curve) hCurve = curves_vec(j);
        EG_completeCurve(geom, hCurve);
      }
    }

    TopExp::Vertices(Edge, V2, V1, Standard_True);
    if (Edge.Orientation() != TopAbs_REVERSED) {
      n1 = body->nodes.map.FindIndex(V2);
      n2 = body->nodes.map.FindIndex(V1);
    } else {
      n1 = body->nodes.map.FindIndex(V1);
      n2 = body->nodes.map.FindIndex(V2);
    }
    if (outLevel > 2)
      printf(" Edge %d:  nodes = %d %d  degen = %d (%lf, %lf)\n",
             j+1, n1, n2, degen, t1, t2);

    if (((n1 == 0) || (n2 == 0)) && (degen == 0)) {
      hit++;
      if (outLevel > 0)
        printf(" EGADS Warning: Node(s) not found for Edge (%d %d)!\n", n1, n2);
      Handle(Geom_Curve) hCurve = curves_vec(j);
      if (n1 == 0) {
        gp_Pnt pnt1;
        hCurve->D0(t1, pnt1);
        if (outLevel > 0)
          printf("                beg: %lf %lf %lf",
                 pnt1.X(), pnt1.Y(), pnt1.Z());
        double dist = 1.e308;
        for (ii = 0; ii < nNode; ii++) {
          TopoDS_Shape  shapv = body->nodes.map(ii+1);
          TopoDS_Vertex vertv = TopoDS::Vertex(shapv);
          gp_Pnt pv           = BRep_Tool::Pnt(vertv);
          double d            = sqrt((pnt1.X()-pv.X())*(pnt1.X()-pv.X()) +
                                     (pnt1.Y()-pv.Y())*(pnt1.Y()-pv.Y()) +
                                     (pnt1.Z()-pv.Z())*(pnt1.Z()-pv.Z()));
          if (d >= dist) continue;
          dist = d;
          n1   = ii+1;
        }
        if (outLevel > 0)
          printf(" vert = %d dist = %le\n", n1, dist);
      }
      if (n2 == 0) {
        gp_Pnt pnt2;
        hCurve->D0(t2, pnt2);
        if (outLevel > 0)
          printf("                end: %lf %lf %lf",
                 pnt2.X(), pnt2.Y(), pnt2.Z());
        double dist = 1.e308;
        for (ii = 0; ii < nNode; ii++) {
          TopoDS_Shape  shapv = body->nodes.map(ii+1);
          TopoDS_Vertex vertv = TopoDS::Vertex(shapv);
          gp_Pnt pv           = BRep_Tool::Pnt(vertv);
          double d            = sqrt((pnt2.X()-pv.X())*(pnt2.X()-pv.X()) +
                                     (pnt2.Y()-pv.Y())*(pnt2.Y()-pv.Y()) +
                                     (pnt2.Z()-pv.Z())*(pnt2.Z()-pv.Z()));
          if (d >= dist) continue;
          dist = d;
          n2   = ii+1;
        }
        if (outLevel > 0)
          printf(" vert = %d dist = %le\n", n2, dist);
      }
    } else if ((n1 == 0) && (n2 == 0) && (degen == 1)) {
      hit++;
      if (outLevel > 0)
        printf(" EGADS Warning: Node not found for Degen Edge (%d)!\n", n1);
      Bnd_Box Box;
      double bbox[6];
      BRepBndLib::Add(Edge, Box);
      Box.Get(bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
      gp_Pnt pnt1(bbox[0], bbox[1], bbox[2]);
      if (outLevel > 0)
        printf("                pnt: %lf %lf %lf",
               pnt1.X(), pnt1.Y(), pnt1.Z());
      double dist = 1.e308;
      for (ii = 0; ii < nNode; ii++) {
        TopoDS_Shape  shapv = body->nodes.map(ii+1);
        TopoDS_Vertex vertv = TopoDS::Vertex(shapv);
        gp_Pnt pv           = BRep_Tool::Pnt(vertv);
        double d            = sqrt((pnt1.X()-pv.X())*(pnt1.X()-pv.X()) +
                                   (pnt1.Y()-pv.Y())*(pnt1.Y()-pv.Y()) +
                                   (pnt1.Z()-pv.Z())*(pnt1.Z()-pv.Z()));
        if (d >= dist) continue;
        dist = d;
        n1   = n2 = ii+1;
      }
      if (outLevel > 0)
        printf(" vert = %d dist = %le\n", n1, dist);
    } else if (((n1 == 0) || (n2 == 0)) && (degen == 1)) {
      if (n1 == 0) n1 = n2;
      if (n2 == 0) n2 = n1;
    }
    egObject *pn1   = body->nodes.objs[n1-1];
    egObject *pn2   = body->nodes.objs[n2-1];

    pedge->edge        = Edge;
    pedge->curve       = geom;
    pedge->nodes[0]    = pn1;
    pedge->nodes[1]    = pn2;
    pedge->topFlg      = 0;
    pedge->bbox.filled = 0;
    BRep_Tool::Range(Edge, pedge->trange[0], pedge->trange[1]);
    pedge->filled      = 0;
    pedge->trange_dot[0] = pedge->trange_dot[1] = 0;
    // special catch for old egads files and the use of zero radius circles
    if ((degen == 0) && (geom->mtype == CIRCLE)) {
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
      GeomAdaptor_Curve AC(hCurve);
      double alen = GCPnts_AbscissaPoint::Length(AC, t1, t2);
      if (alen == 0.0) degen = 1;
    }
    obj->oclass        = EDGE;
    obj->blind         = pedge;
    obj->topObj        = topObj;
    obj->mtype         = TWONODE;
    if (n1 == n2) obj->mtype = ONENODE;
    if (degen == 1) {
      obj->mtype = DEGENERATE;
    } else {
      EG_referenceObject(geom, obj);
    }
    EG_referenceObject(pn1, obj);
    EG_referenceObject(pn2, obj);
  }

  for (j = 0; j < nLoop; j++) {
    int      *senses = NULL, closed = 0, ne = 0;
    egObject **edgeo = NULL;

    egadsLoop *ploop = new egadsLoop;
    obj              = body->loops.objs[j];
    shape            = body->loops.map(j+1);
    obj->oclass      = LOOP;
    if (shape.Closed()) closed = 1;
    TopoDS_Wire Wire = TopoDS::Wire(shape);
    // more reliable for checking closure of Wires
    TopoDS_Vertex V1, V2;
    TopExp::Vertices(Wire, V1, V2);
    if (!V1.IsNull() && !V2.IsNull())
      if (V1.IsSame(V2)) {
        closed = 1;
      } else {
        closed = 0;
      }
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) ne++;
    if (outLevel > 2)
      printf(" Loop %d: # edges = %d, closed = %d\n", j+1, ne, closed);
    if ((outLevel > 0) && (ne == 0))
      printf(" EGADS Info: Loop %d has no Edges!\n", j+1);

    // find the Face
    TopoDS_Face Face;
    int         hit;
    for (hit = k = 0; k < nFace; k++) {
      TopoDS_Shape shapf = body->faces.map(k+1);
      Face = TopoDS::Face(shapf);
      TopExp_Explorer ExpW;
      for (ExpW.Init(shapf, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
        TopoDS_Shape shapw = ExpW.Current();
        TopoDS_Wire  fwire = TopoDS::Wire(shapw);
        if (fwire.IsSame(Wire)) {
          hit++;
          break;
        }
      }
      if (hit != 0) break;
    }
    if ((hit == 0) && (outLevel > 0) && (nFace != 0))
      printf(" EGADS Internal: Loop without a Face!\n");
    if (hit != 0) {
      geom = surfs_obj[surfs_map.FindIndex(surfs_vec(k))-1];
      hit = 2;
      if (geom->mtype == PLANE) hit = 1;
    } else {
      hit = 1;
    }
    if (hit == 1) {
      geom = NULL;
    } else {
      EG_referenceObject(geom, obj);
    }

    if (ne > 0) {
      edgeo  = new egObject*[hit*ne];
      senses = new int[ne];
      for (k = 0; k < hit*ne; k++) edgeo[k] = NULL;
    }
    k = 0;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shapW = ExpWE.Current();
      TopoDS_Edge  Edge  = TopoDS::Edge(shapW);
      int          ed    = body->edges.map.FindIndex(Edge);
      senses[k]          = 1;
      if (shapW.Orientation() == TopAbs_REVERSED) senses[k] = -1;
      if (ed != 0) {
        egObject *eobj = body->edges.objs[ed-1];
        edgeo[k]       = eobj;
        EG_referenceObject(eobj, obj);
        if (outLevel > 2)
          printf("        %d  edge = %d   sense = %d\n", k, ed, senses[k]);
        k++;
      } else {
        printf(" EGADS Warning: Edge not found for Loop!\n");
      }
    }
    ploop->loop        = Wire;
    ploop->surface     = geom;
    ploop->nedges      = k;
    ploop->edges       = edgeo;
    ploop->senses      = senses;
    ploop->topFlg      = 0;
    ploop->bbox.filled = 0;
    obj->blind         = ploop;
    obj->topObj        = topObj;
    obj->mtype         = OPEN;
    if (closed == 1) obj->mtype = CLOSED;
    if (bobj->mtype == WIREBODY) EG_referenceObject(obj, bobj);
    EG_checkLoops(obj);
  }

  int *lcnt = (int *) EG_alloc(nLoop*sizeof(int));
  if (lcnt != NULL) for (j = 0; j < nLoop; j++) lcnt[j] = 0;
  for (j = 0; j < nFace; j++) {
    int      *senses = NULL;
    egObject **loopo = NULL;

    egadsFace *pface = new egadsFace;
    obj              = body->faces.objs[j];
    shape            = body->faces.map(j+1);
    geom             = surfs_obj[surfs_map.FindIndex(surfs_vec(j))-1];
    obj->oclass      = FACE;
    TopoDS_Face Face = TopoDS::Face(shape);
    EG_referenceObject(geom, obj);

    int nl = 0;
    TopExp_Explorer ExpW;
    for (ExpW.Init(shape, TopAbs_WIRE); ExpW.More(); ExpW.Next()) nl++;
    if (outLevel > 2)
      printf(" Face %d: # loops = %d\n", j+1, nl);
    TopoDS_Wire oWire = BRepTools::OuterWire(Face);

    if (nl > 0) {
      loopo  = new egObject*[nl];
      senses = new int[nl];
    }
    k = 0;
    for (ExpW.Init(shape, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Shape shapw = ExpW.Current();
      TopoDS_Wire  Wire  = TopoDS::Wire(shapw);
      loopo[k]           = NULL;
      senses[k]          = -1;
      if (Wire.IsSame(oWire)) senses[k] = 1;
      int lp = body->loops.map.FindIndex(Wire);
      if (lp != 0) {
        if (lcnt != NULL) lcnt[lp-1]++;
        loopo[k] = body->loops.objs[lp-1];
        if (lcnt != NULL) {
          if (lcnt[lp-1] == 1) {
            EG_fillPCurves(Face, geom, loopo[k], topObj);
          } else {
            senses[k] *= 2;
            if (outLevel > 0)
              printf(" EGADS Info: Loop %d found in 2 Faces (sense = %d in Face %d)!\n",
                   lp, senses[k], j+1);
          }
        } else {
          EG_fillPCurves(Face, geom, loopo[k], topObj);
        }
        EG_referenceObject(loopo[k], obj);
      } else {
        printf(" EGADS Warning: Loop not found for Face!\n");
      }
      if (outLevel > 2)
        printf("        %d  loop = %d     outer = %d\n", k, lp, senses[k]);
      k++;
    }
    pface->face        = Face;
    pface->surface     = geom;
    pface->nloops      = nl;
    pface->loops       = loopo;
    pface->senses      = senses;
    pface->topFlg      = 0;
    pface->bbox.filled = 0;
    BRepTools::UVBounds(Face, pface->urange[0], pface->urange[1],
                              pface->vrange[0], pface->vrange[1]);
    obj->blind         = pface;
    obj->topObj        = topObj;
    obj->mtype         = SFORWARD;
    if (Face.Orientation() == TopAbs_REVERSED) obj->mtype = SREVERSE;
    if (bobj->mtype == FACEBODY) EG_referenceObject(obj, bobj);
  }
  if (lcnt != NULL) EG_free(lcnt);

  if (nShell > 0) {
    TopoDS_Shell oShell;
    if (solid == 1) {
      TopoDS_Solid Solid = TopoDS::Solid(body->shape);
      oShell = BRepClass3d::OuterShell(Solid);
      body->senses = new int[nShell];
    }

    for (j = 0; j < nShell; j++) {
      egObject   **faceo = NULL;
      egadsShell *pshell = new egadsShell;
      obj                = body->shells.objs[j];
      shape              = body->shells.map(j+1);
      obj->oclass        = SHELL;
      TopoDS_Shell Shell = TopoDS::Shell(shape);
      if (solid == 1) {
        body->senses[j] = -1;
        if (Shell.IsSame(oShell)) body->senses[j] = 1;
      }

      int nf = 0;
      TopExp_Explorer ExpF;
      for (ExpF.Init(shape, TopAbs_FACE); ExpF.More(); ExpF.Next()) {
        TopoDS_Shape shapf = ExpF.Current();
        if (shapf.Orientation() == TopAbs_INTERNAL ||
            shapf.Orientation() == TopAbs_EXTERNAL) continue;
        nf++;
      }

      if (nf > 0) faceo = new egObject*[nf];

      k = 0;
      for (ExpF.Init(shape, TopAbs_FACE); ExpF.More(); ExpF.Next()) {
        if (faceo == NULL) continue;
        TopoDS_Shape shapf = ExpF.Current();
        if (shapf.Orientation() == TopAbs_INTERNAL ||
            shapf.Orientation() == TopAbs_EXTERNAL) continue;
        TopoDS_Face  Face  = TopoDS::Face(shapf);
        faceo[k]           = NULL;
        int fa = body->faces.map.FindIndex(Face);
        if (fa != 0) {
          faceo[k] = body->faces.objs[fa-1];
          EG_referenceObject(faceo[k], obj);
        } else {
          printf(" EGADS Warning: Face not found for Shell!\n");
        }
        if (outLevel > 2)
          printf(" Shell %d/%d: Face = %d\n", k, j+1, fa);
        k++;
      }
      pshell->shell       = Shell;
      pshell->nfaces      = nf;
      pshell->faces       = faceo;
      pshell->topFlg      = 0;
      pshell->bbox.filled = 0;
      obj->blind          = pshell;
      obj->topObj         = topObj;
      obj->mtype          = EG_shellClosure(pshell, 0);
      if (bobj->mtype >= SHEETBODY) EG_referenceObject(obj, bobj);
    }
  }

  EG_free(curves_obj);
  EG_free(surfs_obj);

  *nerr = hit;
  return EGADS_SUCCESS;
}


int
EG_tolerance(const egObject *topo, double *tol)
{
  int    i, stat;
  double toler;

  *tol = 0.0;
  if  (topo == NULL)                                  return EGADS_NULLOBJ;
  if  (topo->magicnumber != MAGIC)                    return EGADS_NOTOBJ;
  if  (topo->blind == NULL)                           return EGADS_NODATA;
  if ((topo->oclass < NODE) || (topo->oclass > BODY)) return EGADS_NOTTOPO;

  if (topo->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) topo->blind;
    *tol = BRep_Tool::Tolerance(pnode->node);
  } else if (topo->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) topo->blind;
    *tol = BRep_Tool::Tolerance(pedge->edge);
    if (topo->mtype == TWONODE) {
      if (pedge->nodes[0] != NULL) {
        stat = EG_tolerance(pedge->nodes[0], &toler);
        if (stat != EGADS_SUCCESS) return stat;
        if (toler > *tol) *tol = toler;
      }
      if (pedge->nodes[1] != NULL) {
        stat = EG_tolerance(pedge->nodes[1], &toler);
        if (stat != EGADS_SUCCESS) return stat;
        if (toler > *tol) *tol = toler;
      }
    } else {
      if (pedge->nodes[0] != NULL) {
        stat = EG_tolerance(pedge->nodes[0], &toler);
        if (stat != EGADS_SUCCESS) return stat;
        if (toler > *tol) *tol = toler;
      }
    }
  } else if (topo->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) topo->blind;
    for (i = 0; i < ploop->nedges; i++) {
      stat = EG_tolerance(ploop->edges[i], &toler);
      if (stat != EGADS_SUCCESS) return stat;
      if (toler > *tol) *tol = toler;
    }
  } else if (topo->oclass == FACE) {
    egadsFace *pface = (egadsFace *) topo->blind;
    *tol = BRep_Tool::Tolerance(pface->face);
    for (i = 0; i < pface->nloops; i++) {
      stat = EG_tolerance(pface->loops[i], &toler);
      if (stat != EGADS_SUCCESS) return stat;
      if (toler > *tol) *tol = toler;
    }
  } else if (topo->oclass == SHELL) {
    egadsShell *pshell = (egadsShell *) topo->blind;
    for (i = 0; i < pshell->nfaces; i++) {
      stat = EG_tolerance(pshell->faces[i], &toler);
      if (stat != EGADS_SUCCESS) return stat;
      if (toler > *tol) *tol = toler;
    }
  } else {
    egadsBody *pbody = (egadsBody *) topo->blind;
    for (i = 0; i < pbody->shells.map.Extent(); i++) {
      stat = EG_tolerance(pbody->shells.objs[i], &toler);
      if (stat != EGADS_SUCCESS) return stat;
      if (toler > *tol) *tol = toler;
    }
  }

  return EGADS_SUCCESS;
}


int
EG_getTolerance(const egObject *topo, double *tol)
{
  *tol = 0.0;
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->oclass <  NODE)       return EGADS_NOTTOPO;
  if (topo->oclass >= MODEL)      return EGADS_NOTTOPO;

  if (topo->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode != NULL) *tol = BRep_Tool::Tolerance(pnode->node);
  } else if (topo->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) topo->blind;
    if (pedge != NULL) *tol = BRep_Tool::Tolerance(pedge->edge);
  } else if (topo->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) topo->blind;
    if (ploop != NULL)
      for (int i = 0; i < ploop->nedges; i++) {
        egadsEdge *pedge = (egadsEdge *) ploop->edges[i]->blind;
        if (pedge == NULL) continue;
        double toler = BRep_Tool::Tolerance(pedge->edge);
        if (toler > *tol) *tol = toler;
      }
  } else if (topo->oclass == FACE) {
    egadsFace *pface = (egadsFace *) topo->blind;
    if (pface != NULL) *tol = BRep_Tool::Tolerance(pface->face);
  } else if (topo->oclass == SHELL) {
    egadsShell *pshell = (egadsShell *) topo->blind;
    if (pshell != NULL)
      for (int i = 0; i < pshell->nfaces; i++) {
        egadsFace *pface = (egadsFace *) pshell->faces[i]->blind;
        if (pface == NULL) continue;
        double toler = BRep_Tool::Tolerance(pface->face);
        if (toler > *tol) *tol = toler;
      }
  } else {
    egadsBody *pbody = (egadsBody *) topo->blind;
    if (pbody != NULL)
      if (topo->mtype == WIREBODY) {
        int nedge = pbody->edges.map.Extent();
        for (int i = 0; i < nedge; i++) {
          TopoDS_Edge Edge = TopoDS::Edge(pbody->edges.map(i+1));
          double toler     = BRep_Tool::Tolerance(Edge);
          if (toler > *tol) *tol = toler;
        }
      } else {
        int nface = pbody->faces.map.Extent();
        for (int i = 0; i < nface; i++) {
          TopoDS_Face Face = TopoDS::Face(pbody->faces.map(i+1));
          double toler     = BRep_Tool::Tolerance(Face);
          if (toler > *tol) *tol = toler;
        }
      }
  }

  return EGADS_SUCCESS;
}


int
EG_setTolerance(const egObject *topo, double tol)
{
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->oclass <  NODE)       return EGADS_NOTTOPO;
  if (topo->oclass >= MODEL)      return EGADS_NOTTOPO;

  ShapeFix_ShapeTolerance sTol;

  if (topo->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode != NULL)  sTol.SetTolerance(pnode->node, tol);
  } else if (topo->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) topo->blind;
    if (pedge != NULL)  sTol.SetTolerance(pedge->edge, tol);
  } else if (topo->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) topo->blind;
    if (ploop != NULL)  sTol.SetTolerance(ploop->loop, tol, TopAbs_EDGE);
  } else if (topo->oclass == FACE) {
    egadsFace *pface = (egadsFace *) topo->blind;
    if (pface != NULL)  sTol.SetTolerance(pface->face, tol);
  } else if (topo->oclass == SHELL) {
    egadsShell *pshell = (egadsShell *) topo->blind;
    if (pshell != NULL) sTol.SetTolerance(pshell->shell, tol, TopAbs_FACE);
  } else {
    egadsBody *pbody = (egadsBody *) topo->blind;
    if (pbody != NULL)  sTol.SetTolerance(pbody->shape, tol, TopAbs_SHELL);
  }

  return EGADS_SUCCESS;
}


int
EG_delSmallEdges(const egObject *body, double tol, egObject **newBody)
{
  int             stat, cnt, nerr = 0;
  egObject        *obj;
  Standard_Real   t1, t2;
  TopExp_Explorer exp0, exp1;

  *newBody = NULL;
  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->oclass != BODY)       return EGADS_NOTBODY;
  if (body->blind == NULL)        return EGADS_NODATA;
  if (EG_sameThread(body))        return EGADS_CNTXTHRD;

  int outLevel       = EG_outLevel(body);
  egObject *context  = EG_context(body);
  egadsBody *pbodx   = (egadsBody *) body->blind;
  TopoDS_Shape shape = pbodx->shape;

  /* note that Edges are found twice */
  for (exp0.Init(shape, TopAbs_EDGE); exp0.More(); exp0.Next()) {
    TopoDS_Edge edge = TopoDS::Edge(exp0.Current());
    if (BRep_Tool::Degenerated(edge)) continue;
    Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
    GeomAdaptor_Curve AC(hCurve);
    if (GCPnts_AbscissaPoint::Length(AC, t1, t2) < tol) nerr++;
  }
  if (nerr == 0) return EGADS_OUTSIDE;

  Handle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape();

  /* adjust the shape so we can use RemoveSmallEdges */
  for (exp0.Init(shape, TopAbs_EDGE); exp0.More(); exp0.Next()) {
    TopoDS_Edge edge = TopoDS::Edge(exp0.Current());
    if (BRep_Tool::Degenerated(edge)) continue;
    Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
    GeomAdaptor_Curve AC(hCurve);
    if (GCPnts_AbscissaPoint::Length(AC, t1, t2) >= tol) continue;
    TopoDS_Vertex v1, v2;
    if (edge.Orientation() == TopAbs_REVERSED) {
      TopExp::Vertices(edge, v2, v1, Standard_True);
    } else {
      TopExp::Vertices(edge, v1, v2, Standard_True);
    }
    /* remove vertex2 and the Edge */
    for (exp1.Init(shape, TopAbs_EDGE); exp1.More(); exp1.Next()) {
      TopoDS_Edge other = TopoDS::Edge(exp1.Current());
      if (BRep_Tool::Degenerated(other)) continue;
      if (other.IsSame(edge)) continue;
      TopoDS_Vertex o1, o2;
      TopExp::Vertices(other, o1, o2, Standard_True);
      if (o1.IsSame(v2) || o2.IsSame(v2)) {
        TopoDS_Edge newedge;
        ShapeBuild_Edge sbe;
        if (o1.IsSame(v2)) {
          newedge = sbe.CopyReplaceVertices(other, v1, o2);
        } else {
          newedge = sbe.CopyReplaceVertices(other, o1, v1);
        }
        rebuild->Replace(other, newedge);
      }
    }
    rebuild->Remove(v2);
    rebuild->Remove(edge);
  }
  shape = rebuild->Apply(shape);
  cnt   = 0;
  for (exp1.Init(shape, TopAbs_EDGE); exp1.More(); exp1.Next()) {
    TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
    if (BRep_Tool::Degenerated(edge)) continue;
    Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
    GeomAdaptor_Curve AC(hCurve);
    if (GCPnts_AbscissaPoint::Length(AC, t1, t2) >= tol) continue;
    TopoDS_Vertex v1, v2;
    if (edge.Orientation() == TopAbs_REVERSED) {
      TopExp::Vertices(edge, v2, v1, Standard_True);
    } else {
      TopExp::Vertices(edge, v1, v2, Standard_True);
    }
    if (v1.IsSame(v2)) cnt++;
  }
  if (cnt != 0) {
    Handle(ShapeBuild_ReShape) rebuilx = new ShapeBuild_ReShape();
    for (exp1.Init(shape, TopAbs_EDGE); exp1.More(); exp1.Next()) {
      TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
      if (BRep_Tool::Degenerated(edge)) continue;
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
      GeomAdaptor_Curve AC(hCurve);
      if (GCPnts_AbscissaPoint::Length(AC, t1, t2) >= tol) continue;
      TopoDS_Vertex v1, v2;
      if (edge.Orientation() == TopAbs_REVERSED) {
        TopExp::Vertices(edge, v2, v1, Standard_True);
      } else {
        TopExp::Vertices(edge, v1, v2, Standard_True);
      }
      if (v1.IsSame(v2)) rebuilx->Remove(edge);
    }
    shape = rebuilx->Apply(shape);
  }
  shape = ShapeFix::RemoveSmallEdges(shape, tol, rebuild);

/*
  Handle(ShapeFix_FixSmallFace) sffsm = new ShapeFix_FixSmallFace();
  sffsm->Init(shape);
  sffsm->SetPrecision(100.*tol);
  sffsm->Perform();
  shape = sffsm->FixShape();
 */
/*
  Handle(ShapeFix_FixSmallFace) sffsm = new ShapeFix_FixSmallFace();
  sffsm->Init(shape);
  sffsm->SetPrecision(100.*tol);
  shape = sffsm->RemoveSmallFaces();
 */

  /* check and fix if necessary */
  BRepCheck_Analyzer sCheck(shape);
  if (!sCheck.IsValid()) {
    Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(shape);
    sfs->Perform();
    TopoDS_Shape fixShape = sfs->Shape();
    if (fixShape.IsNull()) {
      if (outLevel > 0)
      printf(" EGADS Warning: Shape is invalid (EG_delSmallEdges)!\n");
      return EGADS_CONSTERR;
    }
    shape = fixShape;
  }

  /* make the EGADS Object */
  stat  = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
    printf(" EGADS Error: Cannot make Object (EG_delSmallEdges)!\n");
    return stat;
  }
  obj->oclass        = body->oclass;
  obj->mtype         = body->mtype;
  egadsBody *pbody   = new egadsBody;
  pbody->nodes.objs  = NULL;
  pbody->edges.objs  = NULL;
  pbody->loops.objs  = NULL;
  pbody->faces.objs  = NULL;
  pbody->shells.objs = NULL;
  pbody->senses      = NULL;
  pbody->shape       = shape;
  pbody->bbox.filled = 0;
  pbody->massFill    = 0;
  obj->blind         = pbody;
  stat = EG_traverseBody(context, 0, obj, obj, pbody, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbody;
    return stat;
  }
  EG_attriBodyDup(body, obj);
  EG_referenceObject(obj, context);

  *newBody = obj;
  return EGADS_SUCCESS;
}


int
EG_getTopology(const egObject *topo, egObject **geom, int *oclass,
               int *type, /*@null@*/ double *limits, int *nChildren,
               egObject ***children, int **senses)
{
  *geom      = NULL;
  *oclass    = *type = 0;
  *nChildren = 0;
  *children  = NULL;
  *senses    = NULL;
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->oclass < NODE)        return EGADS_NOTTOPO;
  *oclass = topo->oclass;
  *type   = topo->mtype;

  if (topo->oclass == NODE) {

    egadsNode *pnode = (egadsNode *) topo->blind;
    if ((limits != NULL) && (pnode != NULL)) {
      limits[0] = pnode->xyz[0];
      limits[1] = pnode->xyz[1];
      limits[2] = pnode->xyz[2];
    }

  } else if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    if (pedge != NULL) {
      *geom      = pedge->curve;
      *nChildren = 1;
      if (topo->mtype == TWONODE) *nChildren = 2;
      *children = pedge->nodes;
      if (limits != NULL)
        BRep_Tool::Range(pedge->edge, limits[0], limits[1]);
    }

  } else if (topo->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) topo->blind;
    if (ploop != NULL) {
      *geom      = ploop->surface;
      *nChildren = ploop->nedges;
      *children  = ploop->edges;
      *senses    = ploop->senses;
    }

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    if (pface != NULL) {
      *geom      = pface->surface;
      *nChildren = pface->nloops;
      *children  = pface->loops;
      *senses    = pface->senses;
      if (limits != NULL) {
        Standard_Real umin, umax, vmin, vmax;

        BRepTools::UVBounds(pface->face, umin, umax, vmin, vmax);
        limits[0] = umin;
        limits[1] = umax;
        limits[2] = vmin;
        limits[3] = vmax;
      }
    }

  } else if (topo->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    if (pshell != NULL) {
      *nChildren = pshell->nfaces;
      *children  = pshell->faces;
    }

  } else if (topo->oclass == BODY) {

    egadsBody *pbody = (egadsBody *) topo->blind;
    if (pbody != NULL)
      if (topo->mtype == WIREBODY) {
        *nChildren = pbody->loops.map.Extent();
        *children  = pbody->loops.objs;
      } else if (topo->mtype == FACEBODY) {
        *nChildren = pbody->faces.map.Extent();
        *children  = pbody->faces.objs;
      } else {
        *nChildren = pbody->shells.map.Extent();
        *children  = pbody->shells.objs;
        if (topo->mtype == SOLIDBODY) *senses = pbody->senses;
      }

  } else {

    egadsModel *pmodel = (egadsModel *) topo->blind;
    if (pmodel != NULL) {
      *nChildren = pmodel->nbody;
      *children  = pmodel->bodies;
    }

  }

  return EGADS_SUCCESS;
}


static void
EG_makePCurves(TopoDS_Face& face, egObject *surfo, egObject *loopo,
               Standard_Real prec, int flag)
{
  int      i = 0, j, *periodic=NULL;
  egObject *geom;

  egadsLoop *ploop = (egadsLoop *) loopo->blind;
  if (ploop->surface == NULL) return;
  if (ploop->surface->mtype == PLANE) return;
  if (ploop->surface != surfo) {
    printf(" EGADS Internal: Loop/Face mismatch on Surface (EG_makePCurves)!\n");
    return;
  }

  int outLevel     = EG_outLevel(surfo);
  periodic         = (int*)EG_alloc(ploop->nedges*sizeof(int));
  if (periodic == NULL) {
    printf(" EGADS Error: Cannot malloc (EG_makePCurves)!\n");
  }

  /* look for periodic edges in the loop */
  if (periodic != NULL) {
    for (i = 0; i < ploop->nedges; i++) periodic[i] = 0;
    for (i = 0; i < ploop->nedges; i++)
      for (j = i+1; j < ploop->nedges; j++)
        if (ploop->edges[i] == ploop->edges[j]) {
          periodic[i] =  j;
          periodic[j] = -1;
          break;
        }
  }

  TopoDS_Wire wire = ploop->loop;
  BRep_Builder           Builder;
  BRepTools_WireExplorer ExpWE;
  for (i = 0, ExpWE.Init(wire); ExpWE.More(); ExpWE.Next(), i++) {
    geom                = ploop->edges[ploop->nedges+i];
    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    TopoDS_Shape shape  = ExpWE.Current();
    TopoDS_Edge  edge   = TopoDS::Edge(shape);
    if (ppcurv == NULL) continue;
    Handle(Geom2d_Curve) hCurv2d = ppcurv->handle;

    if ((flag == 0) || (outLevel > 2)) {
      TopoDS_Vertex V1, V2;
      Standard_Real t, t1, t2, delta, mdelta;
      gp_Pnt        pnt, pnte, pv1, pv2;
      gp_Pnt2d      uv;

      egadsSurface *psurf = (egadsSurface *) surfo->blind;
      Handle(Geom_Surface) hSurface = psurf->handle;
      if (edge.Orientation() == TopAbs_REVERSED) {
        TopExp::Vertices(edge, V2, V1, Standard_True);
      } else {
        TopExp::Vertices(edge, V1, V2, Standard_True);
      }
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
      pv1 = BRep_Tool::Pnt(V1);
      pv2 = BRep_Tool::Pnt(V2);
      if (outLevel > 2)
        printf(" PCurve #%d: Limits = %lf %lf    prec = %le\n",
               i, t1, t2, prec);
      hCurv2d->D0(t1, uv);
      hSurface->D0(uv.X(), uv.Y(), pnt);
      mdelta = sqrt((pnt.X()-pv1.X())*(pnt.X()-pv1.X()) +
                    (pnt.Y()-pv1.Y())*(pnt.Y()-pv1.Y()) +
                    (pnt.Z()-pv1.Z())*(pnt.Z()-pv1.Z()));
      if (outLevel > 2)
        printf("            delta for 1st Node     = %le  %lf %lf %lf\n",
               mdelta, pv1.X(), pv1.Y(), pv1.Z());
      delta = 0.0;
      for (int j = 1; j < 36; j++) {
        t = t1 + j*(t2-t1)/36.0;
        hCurv2d->D0(t, uv);
        hSurface->D0(uv.X(), uv.Y(), pnt);
        if (BRep_Tool::Degenerated(edge)) {
          pnte = pv1;
        } else {
          hCurve->D0(t, pnte);
        }
        delta += sqrt((pnt.X()-pnte.X())*(pnt.X()-pnte.X()) +
                      (pnt.Y()-pnte.Y())*(pnt.Y()-pnte.Y()) +
                      (pnt.Z()-pnte.Z())*(pnt.Z()-pnte.Z()));
      }
      delta /= 35;
      if (outLevel > 2)
        printf("            ave delta against Edge = %le\n", delta);
      if (delta > mdelta) mdelta = delta;
      hCurv2d->D0(t2, uv);
      hSurface->D0(uv.X(), uv.Y(), pnt);
      delta = sqrt((pnt.X()-pv2.X())*(pnt.X()-pv2.X()) +
                   (pnt.Y()-pv2.Y())*(pnt.Y()-pv2.Y()) +
                   (pnt.Z()-pv2.Z())*(pnt.Z()-pv2.Z()));
      if (outLevel > 2)
        printf("            delta for 2nd Node     = %le  %lf %lf %lf\n",
               delta, pv2.X(), pv2.Y(), pv2.Z());
      if (delta > mdelta) mdelta = delta;
      if ((flag == 0) && (mdelta*1.001 > prec))
        Builder.SameParameter(edge, Standard_False);
    }

    if (periodic != NULL && periodic[i] < 0) continue; /* periodic edge already processed */
    if (periodic != NULL && periodic[i] > 0) {
      egObject *geomR      = ploop->edges[ploop->nedges+periodic[i]];
      egadsPCurve *ppcurvR = (egadsPCurve *) geomR->blind;
      Handle(Geom2d_Curve) hCurv2dR = ppcurvR->handle;

      if (edge.Orientation() == TopAbs_FORWARD)
        Builder.UpdateEdge(edge, hCurv2d, hCurv2dR, face, prec);
      else
        Builder.UpdateEdge(edge, hCurv2dR, hCurv2d, face, prec);
    } else {
      Builder.UpdateEdge(edge, hCurv2d, face, prec);
    }
  }

  EG_free(periodic);
}


static int
EG_examineFace(TopoDS_Face Face, int nLoop, egObject **loops, int outLevel)
{
  TopTools_IndexedMapOfShape MapL;
  TopExp::MapShapes(Face, TopAbs_WIRE, MapL);
  if (MapL.Extent() != nLoop) {
    if (outLevel > 0)
      printf(" EGADS Info: # Loops input = %d  fixed = %d\n",
             nLoop, MapL.Extent());
    return EGADS_TOPOERR;
  }
  if (outLevel <= 1) return EGADS_SUCCESS;

  int cnt = 0;
  for (int i = 1; i <= nLoop; i++) {
    TopoDS_Shape shape = MapL(i);
    TopoDS_Wire  wire  = TopoDS::Wire(shape);
    int hit = 0;
    for (int j = 0; j < nLoop; j++) {
      egadsLoop *ploop = (egadsLoop *) loops[j]->blind;
      if (wire.IsSame(ploop->loop)) {
        hit = -1;
        if (wire.IsEqual(ploop->loop)) hit = 1;
      }
    }
    if (hit == 0) {
      printf(" EGADS Info: Input Loop %d Not Found\n", i);
    } else if (hit == -1) {
      printf(" EGADS Info: Input Loop %d needs to be Reversed!\n", i);
      cnt++;
    }
  }
  if (cnt == nLoop)
    printf("             Or the Face orientation should be Reversed!\n");

  return EGADS_SUCCESS;
}


int
EG_makeTopology(egObject *context, /*@null@*/ egObject *geom,
                int oclass, int mtypex, /*@null@*/ double *limits,
                int nChildren, /*@null@*/ egObject **children,
                /*@null@*/ int *senses, egObject **topo)
{
  int      i, n, stat, mtype, outLevel, nerr;
  egCntxt  *cntx;
  egObject *obj;

  *topo = NULL;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (EG_sameThread(context))        return EGADS_CNTXTHRD;
  cntx = (egCntxt *) context->blind;
  if (cntx == NULL)                  return EGADS_NODATA;
  outLevel = cntx->outLevel;
  mtype    = mtypex;

  if ((oclass < NODE) || (oclass > MODEL)) {
    if (outLevel > 0)
      printf(" EGADS Error: oclass = %d (EG_makeTopology)!\n", oclass);
    return EGADS_NOTTOPO;
  }

  if (oclass == NODE) {

    if (limits == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Node with no Data (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }

    gp_Pnt pnt(limits[0], limits[1], limits[2]);
    TopoDS_Vertex vert = BRepBuilderAPI_MakeVertex(pnt);
    BRepCheck_Analyzer vCheck(vert);
    if (!vCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Node is invalid (EG_makeTopology)!\n");
      return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Node object (EG_makeTopology)!\n");
      return stat;
    }
    egadsNode *pnode   = new egadsNode;
    pnode->node        = vert;
    pnode->xyz[0]      = limits[0];
    pnode->xyz[1]      = limits[1];
    pnode->xyz[2]      = limits[2];
    pnode->bbox.filled = 0;
    pnode->filled      = 0;
    pnode->xyz_dot[0]  = pnode->xyz_dot[1] = pnode->xyz_dot[2] = 0.0;
    obj->oclass        = NODE;
    obj->blind         = pnode;
    obj->topObj        = context;
    EG_referenceObject(obj, context);

  } else if (oclass == EDGE) {
    TopoDS_Vertex V1, V2;
    Standard_Real P1, P2;

    if (limits == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Limits is NULL (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    if (limits[0] >= limits[1]) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge Tmin (%lf) >= Tmax (%lf) (EG_makeTopology)!\n",
               limits[0], limits[1]);
      return EGADS_RANGERR;
    }
    if ((mtype == TWONODE) && (nChildren != 2)) {
      if (outLevel > 0)
        printf(" EGADS Error: TWONODE Edge with %d Nodes (EG_makeTopology)!\n",
               nChildren);
      return EGADS_TOPOERR;
    }
    if ((mtype == ONENODE) && (nChildren != 1)) {
      if (outLevel > 0)
        printf(" EGADS Error: ONENODE Edge with %d Nodes (EG_makeTopology)!\n",
               nChildren);
      return EGADS_TOPOERR;
    }

    if (mtype == DEGENERATE) {
      if (nChildren != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Degen Edge with %d Nodes (EG_makeTopology)!\n",
                 nChildren);
        return EGADS_TOPOERR;
      }
      if (children[0] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Degen Edge w/ Node NULL (EG_makeTopology)!\n");
        return EGADS_NULLOBJ;
      }
      if (children[0]->oclass != NODE) {
        if (outLevel > 0)
          printf(" EGADS Error: Degen Edge w/ nonNode Child (EG_makeTopology)!\n");
        return EGADS_NOTTOPO;
      }
      if (children[0]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Degen Edge w/ NULL Node Child (EG_makeTopology)!\n");
        return EGADS_NODATA;
      }
      if (EG_context(children[0]) != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Degen Edge w/ BAD Child Cntxt (EG_makeTopology)!\n");
        return EGADS_MIXCNTX;
      }

      // make a degenerate Edge
      egadsNode *pnode = (egadsNode *) children[0]->blind;
      V1               = pnode->node;
      P1               = limits[0];
      P2               = limits[1];

      Handle(Geom_Surface)    hSurf;
      Handle(Geom2d_Curve)    hPcurve;
      gp_Pnt                  pv = BRep_Tool::Pnt(V1);
      // construct temporary BSpline Degenerate Surface
      TColgp_Array2OfPnt      aPoles  (1, 5, 1, 2);
      TColStd_Array2OfReal    aWeights(1, 5, 1, 2);
      TColStd_Array1OfReal    uKnots(1, 3);
      TColStd_Array1OfInteger uMults(1, 3);
      TColStd_Array1OfReal    vKnots(1, 2);
      TColStd_Array1OfInteger vMults(1, 2);
      uKnots(1)     = P1;
      uKnots(2)     = 0.5*(P1 + P2);
      uKnots(3)     = P2;
      uMults(1)     = 3;
      uMults(2)     = 2;
      uMults(3)     = 3;
      vKnots(1)     = 0.0;
      vKnots(2)     = 1.0;
      vMults(1)     = 2;
      vMults(2)     = 2;
      aWeights(1,1) = 1.000000000000;
      aWeights(2,1) = 0.707106781187;
      aWeights(3,1) = 1.000000000000;
      aWeights(4,1) = 0.707106781187;
      aWeights(5,1) = 1.000000000000;
      aWeights(1,2) = 1.000000000000;
      aWeights(2,2) = 0.707106781187;
      aWeights(3,2) = 1.000000000000;
      aWeights(4,2) = 0.707106781187;
      aWeights(5,2) = 1.000000000000;
      aPoles(1,1)   = gp_Pnt(pv.X(), pv.Y(), pv.Z());
      aPoles(2,1)   = aPoles(1,1);
      aPoles(3,1)   = aPoles(1,1);
      aPoles(4,1)   = aPoles(1,1);
      aPoles(5,1)   = aPoles(1,1);
      aPoles(1,2)   = gp_Pnt(pv.X()-0.25, pv.Y()+0.00, pv.Z()-1.0);
      aPoles(2,2)   = gp_Pnt(pv.X()-0.25, pv.Y()+0.25, pv.Z()-1.0);
      aPoles(3,2)   = gp_Pnt(pv.X()+0.00, pv.Y()+0.25, pv.Z()-1.0);
      aPoles(4,2)   = gp_Pnt(pv.X()+0.25, pv.Y()+0.25, pv.Z()-1.0);
      aPoles(5,2)   = gp_Pnt(pv.X()+0.25, pv.Y()+0.00, pv.Z()-1.0);
      try {
        hSurf = new Geom_BSplineSurface(aPoles, aWeights, uKnots, vKnots,
                                        uMults, vMults, 2, 1,
                                        Standard_False, Standard_False);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Surface Creation Error (EG_makeTopology)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Warning: Surface Creation Error (EG_makeTopology)!\n");
        return EGADS_GEOMERR;
      }
      // make the PCurve along the Surface's Degenerate side
      gp_Pnt2d pntl(P1,    0.0);
      gp_Dir2d dirl(P2-P1, 0.0);
      try {
        hPcurve = new Geom2d_Line(pntl, dirl);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: PCurve Creation Error (EG_makeTopology)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Warning: PCurve Creation Error (EG_makeTopology)!\n");
        return EGADS_GEOMERR;
      }
      // Build the Degenerate Edge with our placeholders
      BRepBuilderAPI_MakeEdge MEdge(hPcurve, hSurf, V1, V1, P1, P2);
      if (!MEdge.IsDone()) {
        printf(" EGADS Warning: Edge Creation Error (EG_makeTopology)!\n");
        return EGADS_TOPOERR;
      }
      TopoDS_Edge Edge = MEdge.Edge();
      BRep_Builder Builder;
      Builder.Degenerated(Edge, Standard_True);
      if (!BRep_Tool::Degenerated(Edge))
        printf(" EGADS Info: Degenerate Edge NOT Degenerate!\n");
      BRepCheck_Analyzer eCheck(Edge);
      if (!eCheck.IsValid()) {
        if (outLevel > 0)
          printf(" EGADS Info: Degen Edge is invalid (EG_makeTopology)!\n");
        EG_checkStatus(eCheck.Result(Edge));
        return EGADS_CONSTERR;
      }
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Cannot make Degen Edge obj (EG_makeTopology)!\n");
        return stat;
      }
      egadsEdge *pedge   = new egadsEdge;
      pedge->edge        = Edge;
      pedge->curve       = NULL;
      pedge->nodes[0]    = children[0];
      pedge->nodes[1]    = children[0];
      pedge->topFlg      = 1;
      pedge->bbox.filled = 0;
      BRep_Tool::Range(Edge, pedge->trange[0], pedge->trange[1]);
      pedge->filled      = 0;
      pedge->trange_dot[0] = pedge->trange_dot[1] = 0;
      obj->oclass        = EDGE;
      obj->blind         = pedge;
      obj->topObj        = context;
      obj->mtype         = DEGENERATE;
      EG_referenceTopObj(pedge->nodes[0],  obj);
      EG_referenceTopObj(pedge->nodes[1],  obj);
      EG_referenceObject(obj,          context);

      *topo = obj;
      return EGADS_SUCCESS;
    }

    if (geom == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with NULL Geom (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (geom->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with No Geom (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    if (geom->oclass != CURVE) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge Geom not CURVE (EG_makeTopology)!\n");
      return EGADS_NOTGEOM;
    }
    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with NULL Children (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if ((nChildren != 1) && (nChildren != 2)) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with %d Verts (EG_makeTopology)!\n",
               nChildren);
      return EGADS_TOPOERR;
    }
    if (children[0] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with Vert[0] NULL (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (children[0]->oclass != NODE) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with nonNode Child[0] (EG_makeTopology)!\n");
      return EGADS_NOTTOPO;
    }
    if (children[0]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge w/ NULL Node Child[0] (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    if (EG_context(children[0]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge w/ BAD Child[0] Cntxt (EG_makeTopology)!\n");
      return EGADS_MIXCNTX;
    }
    if (nChildren == 2) {
      if (children[1] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge with Vert[1] NULL (EG_makeTopology)!\n");
        return EGADS_NULLOBJ;
      }
      if (children[1]->oclass != NODE) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge w/ nonNode Child[1] (EG_makeTopology)!\n");
        return EGADS_NOTTOPO;
      }
      if (children[1]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge w/ NULL Node Child[1] (EG_makeTopology)!\n");
        return EGADS_NODATA;
      }
      if (EG_context(children[1]) != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge w/ BAD Child[1] Cntxt (EG_makeTopology)!\n");
        return EGADS_MIXCNTX;
      }
    }
    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    egadsNode  *pnode1 = (egadsNode *)  children[0]->blind;
    egadsNode  *pnode2 = (egadsNode *)  children[0]->blind;
    if (nChildren == 2) pnode2 = (egadsNode *) children[1]->blind;

    P1 = limits[0];
    P2 = limits[1];
    V1 = pnode1->node;
    V2 = pnode2->node;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    gp_Pnt pnt1, pnt2, pv1, pv2;
    hCurve->D0(P1, pnt1);
    hCurve->D0(P2, pnt2);
    pv1 = BRep_Tool::Pnt(V1);
    pv2 = BRep_Tool::Pnt(V2);
    if (outLevel > 2) {
      printf(" P1 = %lf %lf %lf  %lf %lf %lf\n", pnt1.X(), pnt1.Y(), pnt1.Z(),
                pnode1->xyz[0], pnode1->xyz[1], pnode1->xyz[2]);
      printf("      vert = %lf %lf %lf\n", pv1.X(), pv1.Y(), pv1.Z());
      printf(" P2 = %lf %lf %lf  %lf %lf %lf\n", pnt2.X(), pnt2.Y(), pnt2.Z(),
              pnode2->xyz[0], pnode2->xyz[1], pnode2->xyz[2]);
      printf("      vert = %lf %lf %lf\n", pv2.X(), pv2.Y(), pv2.Z());
    }
    double delta1 = sqrt((pnt1.X()-pv1.X())*(pnt1.X()-pv1.X()) +
                         (pnt1.Y()-pv1.Y())*(pnt1.Y()-pv1.Y()) +
                         (pnt1.Z()-pv1.Z())*(pnt1.Z()-pv1.Z()));
    double delta2 = sqrt((pnt2.X()-pv2.X())*(pnt2.X()-pv2.X()) +
                         (pnt2.Y()-pv2.Y())*(pnt2.Y()-pv2.Y()) +
                         (pnt2.Z()-pv2.Z())*(pnt2.Z()-pv2.Z()));

    Standard_Real old  = BRepBuilderAPI::Precision();
    Standard_Real prec = old;
    if (outLevel > 2)
      printf("   Limits = %f %lf, Tol = %le %le   %le\n",
             P1, P2, delta1, delta2, old);
    if (delta1*1.001 > prec) prec = 1.001*delta1;
    if (delta2*1.001 > prec) prec = 1.001*delta2;
    BRepBuilderAPI::Precision(prec);
    BRepBuilderAPI_MakeEdge MEdge;
    MEdge.Init(hCurve, V1, V2, P1, P2);
    BRepBuilderAPI::Precision(old);
    if (!MEdge.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem with the Edge (EG_makeTopology)!\n");
       return EGADS_NODATA;
    }
    TopoDS_Edge Edge = MEdge.Edge();

    if (hCurve->IsPeriodic()) {
      Standard_Real cf = hCurve->FirstParameter();
      Standard_Real cl = hCurve->LastParameter();

      /* undo periodic range adjustments if
         the upper specified parameters is in range of the curve*/
      if ( (P2 >= cf) && (P2 <= cl) ) {
        BRep_Builder Builder;
        Builder.Range(Edge,P1,P2);
      }
    }
    /* check the edge */
    BRepCheck_Analyzer eCheck(Edge);
    if (!eCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Edge is invalid (EG_makeTopology)!\n");
      EG_checkStatus(eCheck.Result(Edge));
      return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Edge object (EG_makeTopology)!\n");
      return stat;
    }
    egadsEdge *pedge   = new egadsEdge;
    pedge->edge        = Edge;
    pedge->curve       = geom;
    pedge->nodes[0]    = children[0];
    pedge->nodes[1]    = children[0];
    pedge->topFlg      = 1;
    pedge->bbox.filled = 0;
    BRep_Tool::Range(Edge, pedge->trange[0], pedge->trange[1]);
    pedge->filled      = 0;
    pedge->trange_dot[0] = pedge->trange_dot[1] = 0;
    obj->oclass        = EDGE;
    obj->blind         = pedge;
    obj->topObj        = context;
    obj->mtype         = ONENODE;
    if (nChildren == 2) {
      obj->mtype      = TWONODE;
      pedge->nodes[1] = children[1];
    }
    EG_referenceTopObj(geom,             obj);
    EG_referenceTopObj(pedge->nodes[0],  obj);
    EG_referenceTopObj(pedge->nodes[1],  obj);
    EG_referenceObject(obj,          context);

  } else if (oclass == LOOP) {

    if ((children == NULL) || (senses == NULL)) {
      if (outLevel > 0)
        printf(" EGADS Error: Loop with NULL Input (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Loop with %d Edges (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    n = 1;
    if (geom != NULL) {
      if (geom->oclass != SURFACE) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop Geom not SURFACE (EG_makeTopology)!\n");
        return EGADS_NOTGEOM;
      }
      if (geom->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop with No Geom Data (EG_makeTopology)!\n");
        return EGADS_NODATA;
      }
      if (geom->mtype == PLANE) {
        if (outLevel > 0)
          printf(" EGADS Info: Loop with Planar Surface (EG_makeTopology)!\n");
      } else {
        n = 2;
      }
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop with Edge[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != EDGE) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop w/ nonEdge Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop w/ NULL Edge Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
      if (EG_context(children[i]) != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop w/ BAD Child[%d] Cntxt (EG_makeTopology)!\n",
                 i);
        return EGADS_MIXCNTX;
      }
      if (n == 1) continue;
      if (children[i+nChildren] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop with PCurve[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i+nChildren]->oclass != PCURVE) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop w/ nonPCurve Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i+nChildren]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop w/ NULL PCurve Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
      if (EG_context(children[i+nChildren]) != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge w/ BAD PCurve[%d] Cntxt (EG_makeTopology)!\n",
                 i);
        return EGADS_MIXCNTX;
      }
    }
    BRepBuilderAPI_MakeWire MW;
    for (i = 0; i < nChildren; i++) {
      egadsEdge *pedge = (egadsEdge *) children[i]->blind;
      TopoDS_Edge edge;
      if ((children[i]->mtype == DEGENERATE) && (geom != NULL)) {
        // always replace degenerate edge with the actual surface/pcurve
        egadsNode    *pnode = (egadsNode *)    pedge->nodes[0]->blind;
        egadsSurface *psurf = (egadsSurface *) geom->blind;
        egadsPCurve  *pcurv = (egadsPCurve *)  children[i+nChildren]->blind;
        if (!BRep_Tool::Degenerated(pedge->edge))
          printf(" EGADS Info: Degenerate Edge %d in Loop NOT Degenerate!\n",
                 i+1);

        // Build the Degenerate Edge with the actual pcurve and surface
        BRepBuilderAPI_MakeEdge MEdge(pcurv->handle, psurf->handle,
                                      pnode->node, pnode->node,
                                      pedge->trange[0], pedge->trange[1]);
        if (!MEdge.IsDone()) {
          printf(" EGADS Error: Degen Edge %d in Loop is invalid (EG_makeTopology)!\n",
                 i+1);
          return EGADS_TOPOERR;
        }
        edge = MEdge.Edge();
        BRep_Builder Builder;
        Builder.Degenerated(edge, Standard_True);

        BRepCheck_Analyzer eCheck(edge);
        if (!eCheck.IsValid()) {
          if (outLevel > 0)
            printf(" EGADS Info: Degen Edge %d in Loop is invalid!\n", i+1);
          EG_checkStatus(eCheck.Result(edge));
          return EGADS_CONSTERR;
        }
        if (!BRep_Tool::Degenerated(edge))
          printf(" EGADS Info: Degenerate Edge %d in Loop NOT Degenerate!\n",
                 i+1);
        pedge->edge = edge;
      } else {
        edge = pedge->edge;
      }

      // may only be required for the first Edge, must be in order!
      if (edge.Orientation() == TopAbs_REVERSED) {
        if (senses[i] ==  1) edge.Orientation(TopAbs_FORWARD);
      } else {
        if (senses[i] == -1) edge.Orientation(TopAbs_REVERSED);
      }
/*
      if (edge.Orientation() == TopAbs_REVERSED) {
        if (senses[i] ==  1) {
          TopoDS_Shape shape = edge.Oriented(TopAbs_FORWARD);
          edge = TopoDS::Edge(shape);
        }
      } else {
        if (senses[i] == -1) {
          TopoDS_Shape shape = edge.Oriented(TopAbs_REVERSED);
          edge = TopoDS::Edge(shape);
        }
      }
*/
      try {
        MW.Add(edge);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Cannot Add Edge %d in Wire (EG_makeTopology)!\n",
               i+1);
        printf("                %s\n", e.GetMessageString());
        return EGADS_TOPOERR;
      }
      catch (...) {
        printf(" EGADS Warning: Cannot Add Edge %d in Wire (EG_makeTopology)!\n",
               i+1);
        return EGADS_TOPOERR;
      }
      if (MW.Error()) {
        if (outLevel > 0)
          printf(" EGADS Error: Problem with Edge %d (EG_makeTopology)!\n",
                 i+1);
        return EGADS_NODATA;
      }
    }
    if (!MW.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem with Loop (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    TopoDS_Wire wire = MW.Wire();

    // validate against the senses
    if (outLevel > 2) {
      i = 0;
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
        TopoDS_Shape shape = ExpWE.Current();
        TopoDS_Edge  Edge  = TopoDS::Edge(shape);
        int sense = 1;
        if (shape.Orientation() == TopAbs_REVERSED) sense = -1;
        egadsEdge *pedge = (egadsEdge *) children[i]->blind;
        if (Edge.IsSame(pedge->edge)) {
          printf("  %d: Edges same senses = %d %d\n",
                 i, senses[i], sense);
        } else {
          printf("  %d: Edges NOT the same senses = %d %d\n",
                 i, senses[i], sense);
        }
        i++;
      }
    }

    BRepCheck_Analyzer wCheck(wire);
    if (!wCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Wire is invalid (EG_makeTopology)!\n");
      EG_checkStatus(wCheck.Result(wire));
      return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Loop object (EG_makeTopology)!\n");
      return stat;
    }
    obj->oclass = LOOP;

    egadsLoop *ploop  = new egadsLoop;
    egObject  **edgeo = new egObject*[n*nChildren];
    int       *esense = new int[nChildren];
    int       closed  = 0;
    if (wire.Closed()) closed = 1;
    // more reliable for checking closure of Wires
    TopoDS_Vertex V1, V2;
    TopExp::Vertices(wire, V1, V2);
    if (!V1.IsNull() && !V2.IsNull())
      if (V1.IsSame(V2)) {
        closed = 1;
      } else {
        closed = 0;
      }
    for (i = 0; i < nChildren; i++) {
      edgeo[i]  = children[i];
      esense[i] = senses[i];
      EG_referenceTopObj(children[i], obj);
      if (n == 1) continue;
      edgeo[i+nChildren] = children[i+nChildren];
      EG_referenceTopObj(children[i+nChildren], obj);
    }
    ploop->loop        = wire;
    ploop->surface     = geom;
    ploop->nedges      = nChildren;
    ploop->edges       = edgeo;
    ploop->senses      = esense;
    ploop->topFlg      = 1;
    ploop->bbox.filled = 0;
    obj->blind         = ploop;
    obj->topObj        = context;
    obj->mtype         = OPEN;
    if (closed == 1) obj->mtype = CLOSED;
    EG_referenceObject(obj, context);
    if (geom != NULL)
      if (geom->mtype == PLANE) {
        ploop->surface = NULL;
      } else {
        EG_referenceTopObj(geom, obj);
      }
    if (mtype == CLOSED)
      if ((outLevel > 0) && (closed == 0))
        printf(" EGADS Info: Wire is Open (EG_makeTopology)!\n");
    if (mtype == OPEN)
      if ((outLevel > 0) && (closed == 1))
        printf(" EGADS Info: Wire is Closed (EG_makeTopology)!\n");

  } else if (oclass == FACE) {

    int ignorePCurves = 0;
    if ((mtype == PCURVE) || (mtype == -PCURVE)) {
      if (outLevel > 1)
        printf(" EGADS Info: ignore PCurves when making Face!\n");
      ignorePCurves = 1;
      mtype /= PCURVE;
    }
    if ((mtype != SFORWARD) && (mtype != SREVERSE)) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with MType = %d (EG_makeTopology)!\n",
               mtype);
      return EGADS_RANGERR;
    }
    if (geom == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with NULL Geom (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (geom->oclass != SURFACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Geom not SURFACE (EG_makeTopology)!\n");
      return EGADS_NOTGEOM;
    }
    if (geom->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with No Geom (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    egadsSurface *psurf = (egadsSurface *) geom->blind;
    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with NULL Children (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (senses == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with NULL Senses (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with %d Loops (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Face with Loop[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != LOOP) {
        if (outLevel > 0)
          printf(" EGADS Error: Face w/ nonLoop Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->mtype != CLOSED) {
        if (outLevel > 0)
          printf(" EGADS Error: Face with OPEN Loop[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Face w/ NULL Loop Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
      egadsLoop *ploop = (egadsLoop *) children[i]->blind;
      if ((ploop->surface != geom) && (geom->mtype != PLANE)) {
        if (outLevel > 0)
          printf(" EGADS Error: Face/Loop[%d] Geom Mismatch (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTGEOM;
      }
      if (EG_context(children[i]) != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Face w/ BAD Child[%d] Cntxt (EG_makeTopology)!\n",
                 i);
        return EGADS_MIXCNTX;
      }
    }

    TopoDS_Face face;
    BRepBuilderAPI_MakeFace MFace;
    Standard_Real old = BRepBuilderAPI::Precision();
    Standard_Real prc = old;
    int nTrys = 5;              // # of tol attempts before setting SameParam = F
    for (int itry = 0; itry < nTrys; itry++) {
      if (itry != 0) face.Nullify();
      if (itry == nTrys-1) {    // last attempt -- set tol back
        prc = old;
        BRepBuilderAPI::Precision(old);
      }
//    this can change the underlying surface!
//             a copy prevents this when the surface is shared between Faces
/*    gp_Trsf               form  = gp_Trsf();
      Handle(Geom_Surface)  hSurf = psurf->handle;
      Handle(Geom_Geometry) nGeom = hSurf->Transformed(form);
      Handle(Geom_Surface)  nSurf = Handle(Geom_Surface)::DownCast(nGeom);
 */
      Handle(Geom_Surface)  nSurf = psurf->handle;
      MFace.Init(nSurf, Standard_False, prc);
      for (i = 0; i < nChildren; i++) {
        egadsLoop *ploop = (egadsLoop *) children[i]->blind;
        TopoDS_Wire wire = ploop->loop;
        if (mtype == SREVERSE) wire.Reverse();
        MFace.Add(wire);
        if (MFace.Error()) {
          if (outLevel > 0)
            printf(" EGADS Error: Problem with Loop %d (EG_makeTopology)!\n",
                   i+1);
          BRepBuilderAPI::Precision(old);
          return EGADS_NODATA;
        }
      }
      if (MFace.IsDone()) {
        face = MFace.Face();
        if (mtype == SREVERSE) {
          face.Orientation(TopAbs_REVERSED);
        } else {
          face.Orientation(TopAbs_FORWARD);
        }
        if (ignorePCurves == 0)
          for (i = 0; i < nChildren; i++)
            EG_makePCurves(face, geom, children[i], prc, nTrys-itry-1);
        BRepLib::SameParameter(face);
        BRepCheck_Analyzer oCheck(face);
        if (oCheck.IsValid()) break;
      }
      prc *= 10.0;
      BRepBuilderAPI::Precision(prc);
      if (outLevel > 1)
        printf(" EGADS Info: Adjusting Precision for Face - itry = %d  prec = %lf\n",
               itry, prc);
    }
    BRepBuilderAPI::Precision(old);
    if (!MFace.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem with the Face (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    for (i = 0; i < nChildren; i++) {
      egadsLoop *ploop = (egadsLoop *) children[i]->blind;
      BRepCheck_Wire wireCheck(ploop->loop);
      try {
        wireCheck.InContext(face);
        BRepCheck_ListOfStatus stati;
        stati = wireCheck.Status();
        if (stati.Extent() > 0) {
          BRepCheck_Status cStatus = stati.First();
          if (outLevel > 0)
            EG_printStatus(cStatus);
          if (cStatus != BRepCheck_NoError) return EGADS_GEOMERR;
        }
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Cannot Check Wire Intersections (EG_makeTopology)!\n");
        printf("                %s\n", e.GetMessageString());
/*      return EGADS_TOPOERR;  */
      }
      catch (...) {
        printf(" EGADS Warning: Cannot Check Wire Intersections (EG_makeTopology)!\n");
/*      return EGADS_TOPOERR;  */
      }
    }
    BRepCheck_Analyzer fCheck(face);
    if (!fCheck.IsValid()) {
      // try to fix the fault
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(face);
      sfs->Perform();
      TopoDS_Shape fixedFace = sfs->Shape();
      if (fixedFace.IsNull()) {
        if (outLevel > 0) {
          printf(" EGADS Info: Invalid Face w/ NULL Fix (EG_makeTopology)!\n");
          EG_checkStatus(fCheck.Result(face));
        }
        return  EGADS_CONSTERR;
      }
      BRepCheck_Analyzer fxCheck(fixedFace);
      if (!fxCheck.IsValid()) {
        if (outLevel > 0) {
          printf(" EGADS Info: Face is invalid (EG_makeTopology)!\n");
          EG_checkStatus(fxCheck.Result(fixedFace));
        }
        return  EGADS_CONSTERR;
      }
      face = TopoDS::Face(fixedFace);
      if (outLevel > 0)
        printf(" EGADS Internal: Face has been fixed (EG_makeTopology)!\n");
      stat = EG_examineFace(face, nChildren, children, outLevel);
      if (stat != EGADS_SUCCESS) return stat;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Face object (EG_makeTopology)!\n");
      return stat;
    }
    obj->oclass = FACE;

    egadsFace *pface = new egadsFace;
    egObject **loopo = new egObject*[nChildren];
    int      *lsense = new int[nChildren];
    for (i = 0; i < nChildren; i++) {
      loopo[i]  = children[i];
      lsense[i] = senses[i];
      EG_referenceTopObj(children[i], obj);
    }
    pface->face        = face;
    pface->surface     = geom;
    pface->nloops      = nChildren;
    pface->loops       = loopo;
    pface->senses      = lsense;
    pface->topFlg      = 1;
    pface->bbox.filled = 0;
    BRepTools::UVBounds(face, pface->urange[0], pface->urange[1],
                              pface->vrange[0], pface->vrange[1]);
    obj->blind         = pface;
    obj->topObj        = context;
    obj->mtype         = mtype;
    EG_referenceTopObj(geom, obj);
    EG_referenceObject(obj,  context);

  } else if (oclass == SHELL) {

    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Shell with NULL Input (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Shell with %d Faces (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Shell with Face[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != FACE) {
        if (outLevel > 0)
          printf(" EGADS Error: Shell w/ nonFace Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Shell w/ NULL Face Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
      if (EG_context(children[i]) != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Shell w/ BAD Child[%d] Cntxt (EG_makeTopology)!\n",
                 i);
        return EGADS_MIXCNTX;
      }
    }
    BRep_Builder builder3D;
    TopoDS_Shell shell;
    builder3D.MakeShell(shell);
    for (i = 0; i < nChildren; i++) {
      egadsFace *pface = (egadsFace *) children[i]->blind;
      builder3D.Add(shell, pface->face);
    }
    BRepLib::SameParameter(shell);
    BRepCheck_Analyzer shCheck(shell);
    if (!shCheck.IsValid()) {
      // try to fix the fault
      Handle_ShapeFix_Shape sss = new ShapeFix_Shape(shell);
      // Do not attempt to orient faces as this might break up the shell
      if (mtypex == OPEN) sss->FixShellTool()->FixOrientationMode() = 0;
      sss->Perform();
      TopoDS_Shape fixedShell = sss->Shape();
      if (fixedShell.IsNull()) {
        if (outLevel > 0) {
          printf(" EGADS Info: Invalid Shell w/ NULL Fix (EG_makeTopology)!\n");
          EG_checkStatus(shCheck.Result(shell));
        }
        if (mtype != DEGENERATE) return EGADS_CONSTERR;
      } else {
        BRepCheck_Analyzer shfxCheck(fixedShell);
        if (!shfxCheck.IsValid() || (fixedShell.ShapeType() != TopAbs_SHELL)) {
          if (outLevel > 0) {
            printf(" EGADS Info: Shell is invalid (EG_makeTopology)!\n");
            EG_checkStatus(shfxCheck.Result(fixedShell));
          }
          if (mtype != DEGENERATE) return EGADS_CONSTERR;
        } else {
          shell = TopoDS::Shell(fixedShell);
        }
      }
    }
    if (mtypex == CLOSED) {
      ShapeFix_Shell fixer;
      if (fixer.FixFaceOrientation(shell) == Standard_True) {
        if (fixer.NbShells() == 1) {
          TopoDS_Shell shellfix = fixer.Shell();

          /* FixFaceOrientation might change the face indexing
           * Sort the faces so they are consistent with the input indexing */
          TopoDS_Shell shellsort;
          builder3D.MakeShell(shellsort);

          int nface = 0;
          TopExp_Explorer Exp;
          for (Exp.Init(shell,TopAbs_FACE); Exp.More(); Exp.Next()) {
            TopExp_Explorer ExpFix;
            for (ExpFix.Init(shellfix,TopAbs_FACE); ExpFix.More(); ExpFix.Next()) {
              if (Exp.Current().IsSame(ExpFix.Current())) {
                builder3D.Add(shellsort, ExpFix.Current());
                nface++;
                break;
              }
            }
          }
          if (nface != nChildren) {
            printf(" EGADS Warning: Unable to sort faces! (EG_makeTopology)!\n");
            shell = shellfix;
          } else {
            /* use the sorted shell */
            shell = shellsort;
          }
        }
      }
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Shell object (EG_makeTopology)!\n");
      return stat;
    }
   obj->oclass = SHELL;

    egadsShell *pshell = new egadsShell;
    egObject   **faceo = new egObject*[nChildren];
    for (i = 0; i < nChildren; i++) {
      faceo[i] = children[i];
      EG_referenceTopObj(children[i], obj);
    }
    pshell->shell       = shell;
    pshell->nfaces      = nChildren;
    pshell->faces       = faceo;
    pshell->topFlg      = 1;
    pshell->bbox.filled = 0;
    obj->blind          = pshell;
    obj->topObj         = context;
    obj->mtype          = EG_shellClosure(pshell, mtype);
    EG_referenceObject(obj, context);
    if (mtype == CLOSED)
      if ((outLevel > 0) && (obj->mtype == OPEN))
        printf(" EGADS Info: Shell is Open (EG_makeTopology)!\n");
    if (mtype == OPEN)
      if ((outLevel > 0) && (obj->mtype == CLOSED))
        printf(" EGADS Info: Shell is Closed (EG_makeTopology)!\n");
    if (mtype == DEGENERATE)
      if (obj->mtype == OPEN) {
        printf(" EGADS Info: Shell is Open (EG_makeTopology)!\n");
      } else {
        printf(" EGADS Info: Shell is Closed (EG_makeTopology)!\n");
      }

  } else if (oclass == BODY) {

    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Body with NULL Children (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Body with %d Children (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    if ((mtype < WIREBODY) || (mtype > SOLIDBODY)) {
      if (outLevel > 0)
        printf(" EGADS Error: Body with mtype = %d (EG_makeTopology)!\n",
               mtype);
      return EGADS_RANGERR;
    }
    int                    cclass = SHELL;
    if (mtype == FACEBODY) cclass = FACE;
    if (mtype == WIREBODY) cclass = LOOP;
    if ((mtype != SOLIDBODY) && (mtype != SHEETBODY) && (nChildren != 1)) {
      if (outLevel > 0)
        printf(" EGADS Error: non Solid/Sheet w/ %d children (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Body with child[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != cclass) {
        if (outLevel > 0)
          printf(" EGADS Error: Body w/ invalid Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Body with NULL Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
      if ((children[i]->mtype != CLOSED) && (mtype == SOLIDBODY)) {
        if (outLevel > 0)
          printf(" EGADS Error: Solid w/ nonClosed Shell[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_RANGERR;
      }
      if (EG_context(children[i]) != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Body w/ BAD Child[%d] Cntxt (EG_makeTopology)!\n",
                 i);
        return EGADS_MIXCNTX;
      }
    }
    int          fixedSld = 0;
    TopoDS_Shape shape;
    if (mtype == WIREBODY) {
      egadsLoop *ploop = (egadsLoop *) children[0]->blind;
      shape = ploop->loop;
    } else if (mtype == FACEBODY) {
      egadsFace *pface = (egadsFace *) children[0]->blind;
      shape = pface->face;
    } else if (mtype == SHEETBODY) {
      if (nChildren == 1) {
        egadsShell *pshell = (egadsShell *) children[0]->blind;
        shape = pshell->shell;
      } else {
        BRep_Builder    builder3D;
        TopoDS_Compound compound;
        builder3D.MakeCompound(compound);
        for (i = 0; i < nChildren; i++) {
          egadsShell *pshell = (egadsShell *) children[i]->blind;
          builder3D.Add(compound, pshell->shell);
        }
        shape = compound;
      }
    } else {
      BRep_Builder builder3D;
      TopoDS_Solid solid;
      builder3D.MakeSolid(solid);
      for (i = 0; i < nChildren; i++) {
        egadsShell *pshell = (egadsShell *) children[i]->blind;
        builder3D.Add(solid, pshell->shell);
      }
      try {
        BRepLib::OrientClosedSolid(solid);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Cannot Orient Solid (EG_makeTopology)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_TOPOERR;
      }
      catch (...) {
        printf(" EGADS Warning: Cannot Orient Solid (EG_makeTopology)!\n");
        return EGADS_TOPOERR;
      }
      BRepCheck_Analyzer sCheck(solid);
      if (!sCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(solid);
        sfs->Perform();
        TopoDS_Shape fixedSolid = sfs->Shape();
        if (fixedSolid.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Solid is invalid (EG_makeTopology)!\n");
          solid.Nullify();
          return EGADS_CONSTERR;
        } else {
          if (outLevel > 0)
            printf(" EGADS Warning: Fixing Solid (EG_makeTopology)!\n");
          solid = TopoDS::Solid(fixedSolid);
          fixedSld = 1;
        }
      }
      shape = solid;
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Body object (EG_makeTopology)!\n");
      return stat;
    }
    obj->oclass        = oclass;
    obj->mtype         = mtype;
    egadsBody *pbody   = new egadsBody;
    pbody->nodes.objs  = NULL;
    pbody->edges.objs  = NULL;
    pbody->loops.objs  = NULL;
    pbody->faces.objs  = NULL;
    pbody->shells.objs = NULL;
    pbody->senses      = NULL;
    pbody->shape       = shape;
    pbody->bbox.filled = 0;
    pbody->massFill    = 0;
    obj->blind         = pbody;
    stat = EG_traverseBody(context, 0, obj, obj, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      EG_deleteObject(obj);
      return stat;
    }
    for (i = 0; i < nChildren; i++) {
      EG_attriBodyDup(children[i], obj);
      stat = EG_bodyTravGeom_dot(children[i], pbody);
      if (stat != EGADS_SUCCESS) {
        EG_deleteObject(obj);
        return stat;
      }
    }
    /* patch up Face Attributes for fixed solid */
    if (fixedSld == 1) {
      for (i = 0; i < nChildren; i++) {
        egadsShell *pshell = (egadsShell *) children[i]->blind;
        for (int j = 0; j < pshell->nfaces; j++) {
          egadsFace *pface = (egadsFace *) pshell->faces[j]->blind;
          int index = pbody->faces.map.FindIndex(pface->face);
          if (index != 0) continue;
          for (int k = 0; k < pbody->faces.map.Extent(); k++) {
            stat = EG_isSame(pshell->faces[j], pbody->faces.objs[k]);
            if (stat != EGADS_SUCCESS) continue;
            EG_attributeDup(pshell->faces[j], pbody->faces.objs[k]);
            break;
          }
        }
      }
    }
    EG_referenceObject(obj, context);

  } else {

    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Model with NULL Children (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Model with %d Bodies (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Model with Body[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != BODY) {
        if (outLevel > 0)
          printf(" EGADS Error: Model w/ nonBody Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->topObj != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Model w/ body[%d] reference (EG_makeTopology)!\n",
                 i);
        return EGADS_REFERCE;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Model w/ NULL Body Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
      if (EG_context(children[i]) != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Model w/ BAD Child[%d] Cntxt (EG_makeTopology)!\n",
                 i);
        return EGADS_MIXCNTX;
      }
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Model object (EG_makeTopology)!\n");
      return stat;
    }
    egadsModel *pmodel  = new egadsModel;
    pmodel->bodies      = new egObject*[nChildren];
    pmodel->nbody       = nChildren;
    pmodel->bbox.filled = 0;
    obj->oclass         = MODEL;
    obj->blind          = pmodel;
    TopoDS_Compound compound;
    BRep_Builder    builder3D;
    builder3D.MakeCompound(compound);
    for (i = 0; i < nChildren; i++) {
      pmodel->bodies[i] = children[i];
      egadsBody *pbody  = (egadsBody *) children[i]->blind;
      builder3D.Add(compound, pbody->shape);
      EG_referenceObject(children[i], obj);
      EG_removeCntxtRef(children[i]);
    }
    pmodel->shape = compound;
    EG_referenceObject(obj, context);

  }

  *topo = obj;
  return EGADS_SUCCESS;
}


int
EG_makeLoop(int nedge, egObject **edges, /*@null@*/ egObject *surf,
            double toler, egObject **result)
{
  int      i, first, stat, outLevel, oclass, mtype, mtypc, close, nn, nc, nnew;
  int      hit, *lsenses, *senses;
  double   tol, nodetol, tlim[2], xyz0[3], xyz1[3], xyz2[3], xyz3[3];
  egObject *context, *ref, *curv, **children, **nodes, *fnode, *geom;
  egObject *nnode = NULL, **newedges, *nds[2];

  *result = NULL;
  geom    = surf;
  if (nedge < 1)                          return EGADS_EMPTY;
  if (edges == NULL)                      return EGADS_NULLOBJ;
  for (first = 0; first < nedge; first++)
    if (edges[first] != NULL) break;
  if (first == nedge)                     return EGADS_NULLOBJ;
  if (edges[first]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (surf != NULL) {
    if (surf->oclass != SURFACE)          return EGADS_NOTGEOM;
    if (surf->mtype  == PLANE) geom = NULL;
  }
  if (EG_sameThread(edges[first]))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(edges[first]);
  context  = EG_context(edges[first]);

  tol = nodetol = 0.0;
  for (i = 0; i < nedge; i++) {
    if (edges[i] == NULL) continue;
    if (edges[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not an EGO (EG_makeLoop)!\n", i+1);
      return EGADS_NOTOBJ;
    }
    if (edges[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d has no data (EG_makeLoop)!\n", i+1);
      return EGADS_NODATA;
    }
    if (EG_context(edges[i]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d Context Mismatch (EG_makeLoop)!\n",
               i+1);
      return EGADS_MIXCNTX;
    }
    if (edges[i]->oclass != EDGE) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is Not an EDGE (EG_makeLoop)!\n", i+1);
      return EGADS_NOTTOPO;
    }
    if (edges[i]->mtype == DEGENERATE) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is DEGENERATE (EG_makeLoop)!\n", i+1);
      return EGADS_DEGEN;
    }
    egadsEdge *pedge = (egadsEdge *) edges[i]->blind;
    double toltmp = BRep_Tool::Tolerance(pedge->edge);
    if (toltmp > tol) tol = toltmp;
    BRep_Tool::Range(pedge->edge, tlim[0], tlim[1]);
    egObject  *curvo  = pedge->curve;
    if (curvo  == NULL) continue;
    egadsCurve *pcurve = (egadsCurve *) curvo->blind;
    if (pcurve == NULL) continue;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    GeomAdaptor_Curve AC(hCurve);
    toltmp = GCPnts_AbscissaPoint::Length(AC, tlim[0], tlim[1]);
    if (nodetol == 0.0) {
      nodetol = toltmp;
    } else {
      if (nodetol > toltmp) nodetol = toltmp;
    }
  }
  nodetol *= 0.25;
  if (nodetol > tol) nodetol = tol;
  if (toler   > tol) tol     = toler;
  lsenses = (int *) EG_alloc(nedge*sizeof(int));
  if (lsenses == NULL) return EGADS_MALLOC;
  newedges = (egObject **) EG_alloc(2*nedge*sizeof(egObject *));
  if (newedges == NULL) {
    EG_free(lsenses);
    return EGADS_MALLOC;
  }
  if (outLevel > 1)
    printf(" EG_makeLoop: Nedge = %d  Tolerance = %le (%le)  nodetol = %le\n",
           nedge, tol, toler, nodetol);

  /* set the first edge */
  xyz0[0] = xyz0[1] = xyz0[2] = 0.0;
  stat = EG_getTopology(edges[first], &ref, &oclass, &mtype, tlim, &nn, &nodes,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    EG_free(newedges);
    EG_free(lsenses);
    return stat;
  }
  stat = EG_getTopology(nodes[0], &ref, &oclass, &mtype, xyz0, &nc, &children,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    EG_free(newedges);
    EG_free(lsenses);
    return stat;
  }
  xyz1[0] = xyz0[0];
  xyz1[1] = xyz0[1];
  xyz1[2] = xyz0[2];
  fnode   = nodes[0];
  if (nn != 1) {
    nnode = nodes[1];
    stat  = EG_getTopology(nodes[1], &ref, &oclass, &mtype, xyz1, &nc, &children,
                           &senses);
    if (stat != EGADS_SUCCESS) {
      EG_free(newedges);
      EG_free(lsenses);
      return stat;
    }
  }
  newedges[0]  = edges[first];
  lsenses[0]   = 1;
  edges[first] = NULL;
  nnew         = 1;
  if (nn == 1) {
    if (geom != NULL) {
      stat = EG_otherCurve(geom, ref, tol, &newedges[1]);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
    }
    /* make the loop */
    stat = EG_makeTopology(context, geom, LOOP, CLOSED, NULL, nnew, newedges,
                           lsenses, result);
    EG_free(newedges);
    EG_free(lsenses);
    if (stat == EGADS_SUCCESS)
      for (i = 0; i < nedge; i++)
        if (edges[i] != NULL) stat++;
    return stat;
  }
  close = OPEN;

  /* serach for next edges (forward) */
  do {
    for (hit = i = 0; i < nedge; i++) {
      if (edges[i] == NULL) continue;
      stat = EG_getTopology(edges[i], &curv, &oclass, &mtypc, tlim, &nn, &nodes,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      if (nn == 1) continue;
      stat = EG_getTopology(nodes[0], &ref, &oclass, &mtype, xyz2, &nc,
                            &children, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      stat = EG_getTopology(nodes[1], &ref, &oclass, &mtype, xyz3, &nc,
                            &children, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      if (sqrt((xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0]) +
               (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1]) +
               (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2])) <= nodetol) {
        lsenses[nnew] =  1;
        xyz1[0] = xyz3[0];
        xyz1[1] = xyz3[1];
        xyz1[2] = xyz3[2];
        nds[0]  = nnode;
        nds[1]  = nnode = nodes[1];
      } else if (sqrt((xyz1[0]-xyz3[0])*(xyz1[0]-xyz3[0]) +
                      (xyz1[1]-xyz3[1])*(xyz1[1]-xyz3[1]) +
                      (xyz1[2]-xyz3[2])*(xyz1[2]-xyz3[2])) <= nodetol) {
        lsenses[nnew] = -1;
        xyz1[0] = xyz2[0];
        xyz1[1] = xyz2[1];
        xyz1[2] = xyz2[2];
        nds[1]  = nnode;
        nds[0]  = nnode = nodes[0];
      } else {
        continue;
      }
      hit++;
      /* are we closed -- xyz1 == xyz0? */
      if (sqrt((xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0]) +
               (xyz1[1]-xyz0[1])*(xyz1[1]-xyz0[1]) +
               (xyz1[2]-xyz0[2])*(xyz1[2]-xyz0[2])) <= nodetol) {
        close  = CLOSED;
        if (lsenses[nnew] == 1) {
          nds[1] = fnode;
        } else {
          nds[0] = fnode;
        }
      }
      if ((nodes[0] == nds[0]) && (nodes[1] == nds[1])) {
        newedges[nnew] = edges[i];
      } else {
        if (outLevel > 1)
          printf(" EG_makeLoop: New Edge for %d\n", i);
        stat = EG_makeTopology(context, curv, EDGE, mtypc, tlim, 2, nds,
                               senses, &newedges[nnew]);
        if (stat != EGADS_SUCCESS) {
          EG_free(newedges);
          EG_free(lsenses);
          return stat;
        }
      }
      edges[i] = NULL;
      nnew++;
    }
  } while ((close == OPEN) && (hit != 0));

  if (close == CLOSED) {
    if (geom != NULL) {
      for (i = 0; i < nnew; i++) {
        stat = EG_getTopology(newedges[i], &curv, &oclass, &mtypc, tlim, &nn,
                              &nodes, &senses);
        if (stat != EGADS_SUCCESS) {
          EG_free(newedges);
          EG_free(lsenses);
          return stat;
        }
        stat = EG_otherCurve(geom, curv, tol, &newedges[nnew+i]);
        if (stat != EGADS_SUCCESS) {
          EG_free(newedges);
          EG_free(lsenses);
          return stat;
        }
      }
    }
    /* make the loop */
    stat = EG_makeTopology(context, geom, LOOP, CLOSED, NULL, nnew, newedges,
                           lsenses, result);
    EG_free(newedges);
    EG_free(lsenses);
    if (stat == EGADS_SUCCESS)
      for (i = 0; i < nedge; i++)
        if (edges[i] != NULL) stat++;
    return stat;
  }

  /* start from the first and look back -- will be open */
  for (i = 0; i < nnew/2; i++) {
    nn           =  nnew-i-1;
    ref          =  newedges[i];
    newedges[i]  =  newedges[nn];
    newedges[nn] =  ref;
    nc           =  lsenses[i];
    lsenses[i]   = -lsenses[nn];
    lsenses[nn]  = -nc;
  }
  xyz1[0] = xyz0[0];
  xyz1[1] = xyz0[1];
  xyz1[2] = xyz0[2];
  nnode   = fnode;

  do {
    for (hit = i = 0; i < nedge; i++) {
      if (edges[i] == NULL) continue;
      stat = EG_getTopology(edges[i], &curv, &oclass, &mtypc, tlim, &nn, &nodes,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      if (nn == 1) continue;
      stat = EG_getTopology(nodes[0], &ref, &oclass, &mtype, xyz2, &nc,
                            &children, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      stat = EG_getTopology(nodes[1], &ref, &oclass, &mtype, xyz3, &nc,
                            &children, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      if (sqrt((xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0]) +
               (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1]) +
               (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2])) <= nodetol) {
        lsenses[nnew] =  1;
        xyz1[0] = xyz3[0];
        xyz1[1] = xyz3[1];
        xyz1[2] = xyz3[2];
        nds[0]  = nnode;
        nds[1]  = nnode = nodes[1];
      } else if (sqrt((xyz1[0]-xyz3[0])*(xyz1[0]-xyz3[0]) +
                      (xyz1[1]-xyz3[1])*(xyz1[1]-xyz3[1]) +
                      (xyz1[2]-xyz3[2])*(xyz1[2]-xyz3[2])) <= nodetol) {
        lsenses[nnew] = -1;
        xyz1[0] = xyz2[0];
        xyz1[1] = xyz2[1];
        xyz1[2] = xyz2[2];
        nds[1]  = nnode;
        nds[0]  = nnode = nodes[0];
      } else {
        continue;
      }
      hit++;
      if ((nodes[0] == nds[0]) && (nodes[1] == nds[1])) {
        newedges[nnew] = edges[i];
      } else {
        if (outLevel > 1)
          printf(" EG_makeLoop: New Edge for %d\n", i);
        stat = EG_makeTopology(context, curv, EDGE, mtypc, tlim, 2, nds,
                               senses, &newedges[nnew]);
        if (stat != EGADS_SUCCESS) {
          EG_free(newedges);
          EG_free(lsenses);
          return stat;
        }
      }
      edges[i] = NULL;
      nnew++;
    }
  } while (hit != 0);

  if (geom != NULL) {
    for (i = 0; i < nnew; i++) {
      stat = EG_getTopology(newedges[i], &curv, &oclass, &mtypc, tlim, &nn,
                            &nodes, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      stat = EG_otherCurve(geom, curv, tol, &newedges[nnew+i]);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
    }
  }
  /* make the loop */
  stat = EG_makeTopology(context, geom, LOOP, OPEN, NULL, nnew, newedges,
                         lsenses, result);
  EG_free(newedges);
  EG_free(lsenses);
  if (stat == EGADS_SUCCESS)
    for (i = 0; i < nedge; i++)
      if (edges[i] != NULL) stat++;
  return stat;
}


int
EG_getPlane(const egObject *object, egObject **plane)
{
  int      stat, oclass, mtype;
  double   tol;
  egObject *ref;

  *plane = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass != LOOP)       return EGADS_GEOMERR;
  if (object->blind == NULL)        return EGADS_NODATA;
  egadsLoop *ploop  = (egadsLoop *) object->blind;
  if (ploop->surface != NULL)       return EGADS_GEOMERR;
  if (EG_sameThread(object))        return EGADS_CNTXTHRD;
  int outLevel = EG_outLevel(object);
  egObject *context = EG_context(object);

  // try to fit a plane
  stat = EG_tolerance(object, &tol);
  if (stat != EGADS_SUCCESS) return stat;
  Handle(Geom_Surface) hSurface;
  int nTrys = 4;              // # of tol attempts
  for (int itry = 0; itry < nTrys; itry++) {
    BRepBuilderAPI_FindPlane FPlane(ploop->loop, tol);
    if (FPlane.Found()) {
      hSurface = FPlane.Plane();
      break;
    }
    tol *= 10.0;
    if (outLevel > 1)
      printf(" EGADS Info: Adjusting Prec for getPlane - itry = %d  prec = %le\n",
             itry, tol);
  }
  if (hSurface.IsNull()) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Planar Surface (EG_getPlane)!\n");
    return EGADS_GEOMERR;
  }

  egObject *obj;
  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: make Surface = %d (EG_getPlane)!\n", stat);
    return stat;
  }
  obj->oclass         = SURFACE;
  obj->mtype          = PLANE;
  egadsSurface *psurf = new egadsSurface;
  psurf->handle       = hSurface;
  psurf->ref          = NULL;
  psurf->topFlg       = 1;
  psurf->header       = NULL;
  psurf->data         = NULL;
  psurf->data_dot     = NULL;
  obj->blind          = psurf;
  EG_getGeometry(obj, &oclass, &mtype, &ref, &psurf->header, &psurf->data);
  EG_getGeometryLen(obj, &mtype, &psurf->dataLen);
  hSurface->Bounds(psurf->urange[0], psurf->urange[1],
                   psurf->vrange[0], psurf->vrange[1]);
  EG_referenceObject(obj, context);

  *plane = obj;
  return EGADS_SUCCESS;
}


int
EG_getArea(egObject *object, /*@null@*/ const double *limits, double *area)
{
  double      sense = 1.0;
  TopoDS_Face Face;

  *area = 0.0;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != SURFACE) && (object->oclass != LOOP) &&
      (object->oclass != FACE))     return EGADS_GEOMERR;
  if (object->blind == NULL)        return EGADS_NODATA;
  int outLevel = EG_outLevel(object);

  if (object->oclass == FACE) {

    egadsFace *pface = (egadsFace *) object->blind;
    Face = pface->face;

  } else if (object->oclass == SURFACE) {

    if (limits == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface with NULL Limits (EG_getArea)!\n");
      return EGADS_NODATA;
    }
    egadsSurface *psurf = (egadsSurface *) object->blind;
    Standard_Real tol   = BRepLib::Precision();
    BRepLib_MakeFace MFace(psurf->handle, limits[0], limits[1],
                                          limits[2], limits[3], tol);
    Face = MFace.Face();

  } else {

    egadsLoop *ploop  = (egadsLoop *) object->blind;
    if (ploop->surface == NULL) {

      // try to fit a plane
      Standard_Real tol = Precision::Confusion();
      Handle(Geom_Surface) hSurface;
      int nTrys = 4;              // # of tol attempts
      for (int itry = 0; itry < nTrys; itry++) {
        BRepBuilderAPI_FindPlane FPlane(ploop->loop, tol);
        if (FPlane.Found()) {
          hSurface = FPlane.Plane();
          break;
        }
        tol *= 10.0;
        if (outLevel > 1)
          printf(" EGADS Info: Adjusting Prec for makeFace - itry = %d  prec = %le\n",
                 itry, tol);
      }
      if (hSurface.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Error: Cannot make Planar Surface (EG_getArea)!\n");
        return EGADS_GEOMERR;
      }

      try {
        BRepLib_MakeFace MFace(hSurface, ploop->loop);
        Face = MFace.Face();
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Warning: Cannot makeFace (EG_getArea)!\n");
        printf("                %s\n", e.GetMessageString());
        return EGADS_CONSTERR;
      }
      catch (...) {
        printf(" EGADS Warning: Cannot makeFace (EG_getArea)!\n");
        return EGADS_CONSTERR;
      }

      // did making the Face flip the Loop?
      TopExp_Explorer ExpW;
      for (ExpW.Init(Face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
        TopoDS_Shape shapw = ExpW.Current();
        TopoDS_Wire  wire  = TopoDS::Wire(shapw);
        if (!wire.IsSame(ploop->loop)) continue;
        if (!wire.IsEqual(ploop->loop)) sense = -1.0;
      }

    } else {

      // make a standard Face
      egObject *geom = ploop->surface;
      if (geom->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop had NULL Ref Surface (EG_getArea)!\n");
        return EGADS_NOTGEOM;
      }
      egadsSurface *psurf = (egadsSurface *) geom->blind;
      BRepBuilderAPI_MakeFace MFace;
      Standard_Real old = BRepBuilderAPI::Precision();
      Standard_Real prc = old;
      int nTrys = 5;              // # of tol attempts before setting SameParam = F
      for (int itry = 0; itry < nTrys; itry++) {
        if (itry != 0) Face.Nullify();
        if (itry == nTrys-1) {    // last attempt -- set tol back
          prc = old;
          BRepBuilderAPI::Precision(old);
        }
        MFace.Init(psurf->handle, Standard_False, prc);
        MFace.Add(ploop->loop);
        if (MFace.Error()) {
          if (outLevel > 0)
            printf(" EGADS Error: Problem with Loop (EG_getArea)!\n");
          return EGADS_NODATA;
        }
        if (MFace.IsDone()) {
          Face = MFace.Face();
          EG_makePCurves(Face, ploop->surface, object, prc, nTrys-itry-1);
          BRepLib::SameParameter(Face);
          BRepCheck_Analyzer oCheck(Face);
          if (oCheck.IsValid()) break;
        }
        prc *= 10.0;
        BRepBuilderAPI::Precision(prc);
        if (outLevel > 1)
          printf(" EGADS Info: Adjusting Precision for Face - itry = %d  prec = %lf\n",
                 itry, prc);
      }
      BRepBuilderAPI::Precision(old);
      if (!MFace.IsDone()) {
        if (outLevel > 0)
          printf(" EGADS Error: Problem making the Face (EG_getArea)!\n");
        return EGADS_NODATA;
      }

    }

  }

  BRepGProp    BProps;
  GProp_GProps SProps;
  BProps.SurfaceProperties(Face, SProps);
  *area = sense*SProps.Mass();

  return EGADS_SUCCESS;
}


int
EG_getUVbox(const egObject *face, const egObject *loop, double *box)
{
  int i, outLevel;

  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_TOPOERR;
  if (face->blind == NULL)        return EGADS_NODATA;
  if (loop == NULL)               return EGADS_NULLOBJ;
  if (loop->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (loop->oclass != LOOP)       return EGADS_TOPOERR;
  if (loop->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(face);

  egadsFace *pface = (egadsFace *) face->blind;
  for (i = 0; i < pface->nloops; i++)
    if (loop == pface->loops[i]) break;
  if (i == pface->nloops) {
    if (outLevel > 0)
      printf(" EGADS Error: Loop not in Face (EG_getUVbox)!\n");
    return EGADS_TOPOERR;
  }
  egadsLoop *ploop = (egadsLoop *) loop->blind;
  BRepTools::UVBounds(pface->face, ploop->loop, box[0], box[1], box[2], box[3]);

  return EGADS_SUCCESS;
}


int
EG_getUVinfo(egObject *surface, egObject *loop, double *box, double *area)
{
  if (surface == NULL)               return EGADS_NULLOBJ;
  if (surface->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (surface->oclass != SURFACE)    return EGADS_TOPOERR;
  if (surface->blind == NULL)        return EGADS_NODATA;
  if (loop == NULL)                  return EGADS_NULLOBJ;
  if (loop->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (loop->oclass != LOOP)          return EGADS_TOPOERR;
  if (loop->blind == NULL)           return EGADS_NODATA;

  *area               = 0.0;
  box[0] = box[1]     = box[2] = box[3] = 0.0;
  int outLevel        = EG_outLevel(surface);
  egadsSurface *psurf = (egadsSurface *) surface->blind;
  egadsLoop *ploop    = (egadsLoop *) loop->blind;
  TopoDS_Wire wire    = ploop->loop;

  BRepBuilderAPI_MakeFace MFace;
  Standard_Real old = BRepBuilderAPI::Precision();
  Standard_Real prc = old;
  int nTrys = 5;
  for (int itry = 0; itry < nTrys; itry++) {
    Handle(Geom_Surface) nSurf = psurf->handle;
    MFace.Init(nSurf, Standard_False, prc);
    MFace.Add(wire);
    if (MFace.Error()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem with Loop (EG_getUVinfo)!\n");
      BRepBuilderAPI::Precision(old);
      return EGADS_NODATA;
    }
    if (MFace.IsDone()) break;

    prc *= 10.0;
    BRepBuilderAPI::Precision(prc);
    if (outLevel > 1)
      printf(" EGADS Info: Adjusting Precision for Face - itry = %d  prec = %lf\n",
             itry, prc);
  }
  BRepBuilderAPI::Precision(old);
  if (!MFace.IsDone()) {
    if (outLevel > 0)
      printf(" EGADS Error: Problem with the Face (EG_getUVinfo)!\n");
    return EGADS_NODATA;
  }
  TopoDS_Face Face = MFace.Face();
  EG_makePCurves(Face, surface, loop, prc, 0);
  BRepLib::SameParameter(Face);

  BRepCheck_Analyzer fCheck(Face);
  if (!fCheck.IsValid()) {
    Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(Face);
    sfs->Perform();
    TopoDS_Shape fixedFace = sfs->Shape();
    if (fixedFace.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Info: Fixing Face Failed (EG_getUVinfo)!\n");
      return EGADS_CONSTERR;
    } else {
      BRepCheck_Analyzer sfCheck(fixedFace);
      if (!sfCheck.IsValid()) {
        if (outLevel > 0) {
          printf(" EGADS Info:  Face is invalid (EG_getUVinfo)!\n");
          EG_checkStatus(sfCheck.Result(fixedFace));
        }
        return EGADS_CONSTERR;
      } else {
        Face = TopoDS::Face(fixedFace);
      }
    }
  }

  BRepTools::UVBounds(Face, wire, box[0], box[1], box[2], box[3]);
  Handle(ShapeExtend_WireData) seWire = new ShapeExtend_WireData(ploop->loop);
  *area = ShapeAnalysis::TotCross2D(seWire, Face);

  return EGADS_SUCCESS;
}


int
EG_makeFace(egObject *object, int mtype,
            /*@null@*/ const double *limits, egObject **face)
{
  int      stat, outLevel, *mark, nl = 1;
  egObject *obj, *context, *loop = NULL, *geom = NULL;

  *face = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != SURFACE) && (object->oclass != LOOP) &&
      (object->oclass != FACE))     return EGADS_GEOMERR;
  if (object->blind == NULL)        return EGADS_NODATA;
  if (EG_sameThread(object))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);

  if (object->oclass != FACE)
    if ((mtype != SFORWARD) && (mtype != SREVERSE)) {
      if (outLevel > 0)
        printf(" EGADS Error: Mtype = %d (EG_makeFace)!\n", mtype);
      return EGADS_TOPOERR;
    }
  if (object->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) object->blind;
    if (ploop->surface != NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Loop had Ref Surface (EG_makeFace)!\n");
      return EGADS_NOTGEOM;
    }
  } else {
    if (limits == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface/Face with NULL Limits (EG_makeFace)!\n");
      return EGADS_NODATA;
    }
  }

  TopoDS_Face Face;
  if (object->oclass == SURFACE) {

    egadsSurface *psurf = (egadsSurface *) object->blind;
    Handle(Geom_Surface) hSurf = psurf->handle;
    Standard_Real tol = BRepLib::Precision();
    BRepLib_MakeFace MFace(hSurf, limits[0], limits[1],
                                  limits[2], limits[3], tol);
    Face = MFace.Face();
    if (mtype == SREVERSE) {
      Face.Orientation(TopAbs_REVERSED);
    } else {
      Face.Orientation(TopAbs_FORWARD);
    }
    BRepLib::SameParameter(Face);
    BRepCheck_Analyzer fCheck(Face);
    if (!fCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Face may be invalid (EG_makeFace)!\n");
      EG_checkStatus(fCheck.Result(Face));
      return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Face object (EG_makeFace)!\n");
      return stat;
    }
    obj->oclass = FACE;
    geom        = object;
    TopExp_Explorer ExpW;
    for (ExpW.Init(Face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Shape shapw = ExpW.Current();
      TopoDS_Wire  Wire  = TopoDS::Wire(shapw);
      stat = EG_makeObject(context, &loop);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Cannot make Surface object (EG_makeFace)!\n");
        obj->oclass = NIL;
        EG_deleteObject(obj);
        return stat;
      }
      loop->oclass       = LOOP;
      egadsLoop *ploop   = new egadsLoop;
      loop->blind        = ploop;
      ploop->loop        = Wire;
      ploop->nedges      = 0;
      ploop->edges       = NULL;
      ploop->senses      = NULL;
      ploop->topFlg      = 0;
      ploop->surface     = NULL;
      ploop->bbox.filled = 0;
      if (object->mtype != PLANE) {
        ploop->surface = geom;
        EG_referenceObject(geom, loop);
      }
      EG_fillTopoObjs(loop, obj);
      EG_fillPCurves(Face, geom, loop, obj);
      EG_checkLoops(loop);
      break;
    }

  } else if (object->oclass == FACE) {

    egadsFace *pfasrc = (egadsFace *) object->blind;
    geom = pfasrc->surface;
    if (geom->mtype != PLANE) {
      if (outLevel > 0)
        printf(" EGADS Error: Face is NOT Planar (EG_makeFace)!\n");
      return EGADS_GEOMERR;
    }
    if (pfasrc->nloops != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: Face has multiple Loops (EG_makeFace)!\n");
      return EGADS_GEOMERR;
    }
    BRepOffsetAPI_MakeOffset mOffset(pfasrc->face);
    mOffset.Perform(-limits[0]);
    mOffset.Build();
    TopoDS_Shape shapf = mOffset.Shape();

    int nWire = 0;
    TopExp_Explorer Exp;
    if (shapf.ShapeType() == TopAbs_WIRE) {
      TopoDS_Wire wire = TopoDS::Wire(shapf);
      BRepCheck_Analyzer sCheck(wire);
      if (!sCheck.IsValid()) {
        if (outLevel > 0)
          printf(" EGADS Info: Offset Wire is NOT valid (EG_makeFace)!\n");
        return EGADS_CONSTERR;
      }
      nWire++;
//    printf(" Is a Wire!\n");
    } else if (shapf.ShapeType() == TopAbs_COMPOUND) {
      for (Exp.Init(shapf, TopAbs_WIRE); Exp.More(); Exp.Next()) {
        TopoDS_Wire wire = TopoDS::Wire(Exp.Current());
        BRepCheck_Analyzer wCheck(wire);
        if (!wCheck.IsValid()) {
          if (outLevel > 0)
            printf(" EGADS Info: Offset Wire is NOT valid (EG_makeFace)!\n");
          return EGADS_CONSTERR;
        }
        nWire++;
      }
//    printf(" Is a Compound -- with %d Wires!\n", nWire);
    } else {
      printf(" EGADS Info: What is this shape (EG_makeFace)?!\n");
    }
    if (nWire == 0) {                             // not really an error
      *face = object;                             //  just no hole!
      return EGADS_OUTSIDE;
    }
    TopoDS_Wire *wires = new TopoDS_Wire[nWire];
    if (shapf.ShapeType() == TopAbs_WIRE) {
      wires[0] = TopoDS::Wire(shapf);
    } else {
      int i = 0;
      for (Exp.Init(shapf, TopAbs_WIRE); Exp.More(); Exp.Next()) {
        wires[i] = TopoDS::Wire(Exp.Current());
        i++;
      }
    }

    // fillet the new Wires (holes)?
    if (limits[1] > 0.0) {
      BRepFilletAPI_MakeFillet2d f2d;
      egadsSurface *psurf = (egadsSurface *) geom->blind;
      for (int i = 0; i < nWire; i++) {
        int OK = 1;
        BRepLib_MakeFace MFace(psurf->handle, wires[i]);
        TopoDS_Face Facex = MFace.Face();

        TopTools_IndexedMapOfShape vmap;
        TopExp::MapShapes(wires[i], TopAbs_VERTEX, vmap);
        if (vmap.Extent() == 0) continue;
        mark = (int *) EG_alloc(vmap.Extent()*sizeof(int));
        if (mark == NULL) continue;

        f2d.Init(Facex);
        for (int j = 0; j < vmap.Extent(); j++) {
          TopoDS_Shape  shapv = vmap(j+1);
          TopoDS_Vertex vert  = TopoDS::Vertex(shapv);
          f2d.AddFillet(vert, limits[1]);
          mark[j] = 0;
          ChFi2d_ConstructionError err = f2d.Status();
/*
          printf(" Vertex %d of %d ->", j+1, vmap.Extent());
          if (err == ChFi2d_NotPlanar)            printf(" status = NotPlanar!\n");
          if (err == ChFi2d_NoFace)               printf(" status = NoFace!\n");
          if (err == ChFi2d_InitialisationError)  printf(" status = InitErr!\n");
          if (err == ChFi2d_ParametersError)      printf(" status = ParamErr!\n");
          if (err == ChFi2d_Ready)                printf(" status = Ready!\n");
          if (err == ChFi2d_IsDone)               printf(" status = OK!\n");
          if (err == ChFi2d_ComputationError)     printf(" status = ComputErr!\n");
          if (err == ChFi2d_ConnexionError)       printf(" status = Connection!\n");
          if (err == ChFi2d_TangencyError)        printf(" status = TangenErr!\n");
          if (err == ChFi2d_FirstEdgeDegenerated) printf(" status = FirstDegen!\n");
          if (err == ChFi2d_LastEdgeDegenerated)  printf(" status = LastDegen!\n");
          if (err == ChFi2d_BothEdgesDegenerated) printf(" status = BothDegen!\n");
          if (err == ChFi2d_NotAuthorized)        printf(" status = NotAuthor!\n");
*/
          if (err == ChFi2d_IsDone) mark[j] = 1;
        }

        f2d.Init(Facex);
        int cnt = 0;
        for (int j = 0; j < vmap.Extent(); j++) {
          if (mark[j] == 0) continue;
          TopoDS_Shape  shapv = vmap(j+1);
          TopoDS_Vertex vert  = TopoDS::Vertex(shapv);
          f2d.AddFillet(vert, limits[1]);
          cnt++;
        }
        EG_free(mark);
        if (cnt == 0) continue;
        try {
          TopoDS_Shape shape = f2d.Shape();
//        if (shape.ShapeType() == TopAbs_FACE) printf(" filletted shape: Face!\n");
          Face = TopoDS::Face(shape);
        }
        catch (const Standard_Failure& e) {
          printf(" EGADS Error: Cannot get Filletted Face (EG_makeFace)!\n");
          printf("                %s\n", e.GetMessageString());
/*        delete [] wires;
          return EGADS_CONSTERR;  */
          OK = 0;
        }
        catch (...) {
          printf(" EGADS Error: Cannot get Filletted Face (EG_makeFace)!\n");
/*        delete [] wires;
          return EGADS_CONSTERR;  */
          OK = 0;
        }

        // update the wires with the filletted version
        if (OK == 1) {
          TopTools_IndexedMapOfShape wmap;
          TopExp::MapShapes(Face, TopAbs_WIRE, wmap);
          if (wmap.Extent() != 1)
            if (outLevel > 0)
              printf(" EGADS Warning: Fillet %d #Loops = %d (EG_makeFace)!\n",
                     i+1, wmap.Extent());
          wires[i] = TopoDS::Wire(wmap(1));
        }
      }
    }

    // make the new Face
    BRepBuilderAPI_MakeFace mFace(pfasrc->face);
    for (int i = 0; i < nWire; i++) {
      if (object->mtype == SFORWARD) wires[i].Reverse();   // reverse holes
      mFace.Add(wires[i]);
    }
    if (!mFace.IsDone()) {
      delete [] wires;
      printf(" EGADS Warning: BRepBuilderAPI_MakeFace Error (EG_makeFace)!\n");
      return EGADS_GEOMERR;
    }
    Face = mFace.Face();

    BRepCheck_Analyzer fCheck(Face);
    if (!fCheck.IsValid()) {
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(Face);
      sfs->Perform();
      TopoDS_Shape fixedFace = sfs->Shape();
      if (fixedFace.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Info: Fixing Face Failed (EG_makeFace)!\n");
        delete [] wires;
        return EGADS_CONSTERR;
      } else {
        BRepCheck_Analyzer sfCheck(fixedFace);
        if (!sfCheck.IsValid()) {
          if (outLevel > 0)
            printf(" EGADS Info: Face may be invalid (EG_makeFace)!\n");
          delete [] wires;
          return EGADS_CONSTERR;
        } else {
          Face = TopoDS::Face(fixedFace);
        }
      }
    }

    // make the EGADS objects
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      delete [] wires;
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Face object (EG_makeFace)!\n");
      return stat;
    }
    obj->oclass      = FACE;
    egObject **loopo = new egObject*[nWire+1];
    int *senses      = new int[nWire+1];
    loopo[0]         = pfasrc->loops[0];
    senses[0]        = 1;
    for (int i = 0; i < nWire; i++) {
      stat = EG_makeObject(context, &loop);
      if (stat != EGADS_SUCCESS) {
        delete [] loopo;
        delete [] senses;
        delete [] wires;
        if (outLevel > 0)
          printf(" EGADS Error: Cannot make Loop object (EG_makeFace)!\n");
        obj->oclass = NIL;
        EG_deleteObject(obj);
        return stat;
      }
      loop->oclass       = LOOP;
      egadsLoop *ploop   = new egadsLoop;
      loop->blind        = ploop;
      ploop->loop        = wires[i];
      ploop->nedges      = 0;
      ploop->edges       = NULL;
      ploop->senses      = NULL;
      ploop->topFlg      = 0;
      ploop->surface     = NULL;
      ploop->bbox.filled = 0;
      EG_fillTopoObjs(loop, obj);
      EG_fillPCurves(Face, geom, loop, obj);
      EG_checkLoops(loop);
      loopo[i+1]       = loop;
      senses[i+1]      = -1;
    }
    delete [] wires;

    egadsFace *pface   = new egadsFace;
    pface->face        = Face;
    pface->surface     = geom;
    pface->nloops      = nWire+1;
    pface->loops       = loopo;
    pface->senses      = senses;
    pface->topFlg      = 0;
    pface->bbox.filled = 0;
    BRepTools::UVBounds(Face, pface->urange[0], pface->urange[1],
                              pface->vrange[0], pface->vrange[1]);
    obj->blind       = pface;
    obj->mtype       = object->mtype;

    EG_referenceObject(geom, obj);
    for (int i = 0; i <= nWire; i++) EG_referenceObject(loopo[i], obj);
    EG_referenceObject(obj,  context);
    EG_attributeDup(object,  obj);

    *face = obj;
    return EGADS_SUCCESS;

  } else {

    egadsLoop *ploop  = (egadsLoop *) object->blind;
    Standard_Real tol = Precision::Confusion();
    Handle(Geom_Surface) hSurface;
    int nTrys = 4;              // # of tol attempts
    for (int itry = 0; itry < nTrys; itry++) {
      BRepBuilderAPI_FindPlane FPlane(ploop->loop, tol);
      if (FPlane.Found()) {
        hSurface = FPlane.Plane();
        break;
      }
      tol *= 10.0;
      if (outLevel > 1)
        printf(" EGADS Info: Adjusting Prec for makeFace - itry = %d  prec = %le\n",
               itry, tol);
    }
    if (hSurface.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Planar Surface (EG_makeFace)!\n");
      return EGADS_GEOMERR;
    }

    try {
      BRepLib_MakeFace MFace(hSurface, ploop->loop);
      Face = MFace.Face();
      if (mtype == SREVERSE) {
        Face.Orientation(TopAbs_REVERSED);
      } else {
        Face.Orientation(TopAbs_FORWARD);
      }
      BRepLib::SameParameter(Face);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Warning: Cannot makeFace (EG_makeFace)!\n");
      printf("                %s\n", e.GetMessageString());
      return EGADS_CONSTERR;
    }
    catch (...) {
      printf(" EGADS Warning: Cannot makeFace (EG_makeFace)!\n");
      return EGADS_CONSTERR;
    }

/*  BRepCheck_Analyzer fCheck(Face);
    if (!fCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Face may be invalid (EG_makeFace)!\n");
      return EGADS_CONSTERR;
    }  */
    BRepCheck_Wire wireCheck(ploop->loop);
    wireCheck.InContext(Face);
    BRepCheck_ListOfStatus stati;
    stati = wireCheck.Status();
    if (stati.Extent() > 0) {
      BRepCheck_Status cStatus = stati.First();
      if (outLevel > 0)
        EG_printStatus(cStatus);
      if (cStatus != BRepCheck_NoError) return EGADS_GEOMERR;
    }

    BRepCheck_Analyzer fCheck(Face);
    if (!fCheck.IsValid()) {
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(Face);
      sfs->Perform();
      TopoDS_Shape fixedFace = sfs->Shape();
      if (fixedFace.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Info: Fixing Face Failed (EG_makeFace)!\n");
        return EGADS_CONSTERR;
      } else {
        BRepCheck_Analyzer sfCheck(fixedFace);
        if (!sfCheck.IsValid()) {
          if (outLevel > 0)
            printf(" EGADS Info: Face may be invalid (EG_makeFace)!\n");
          return EGADS_CONSTERR;
        } else {
          Face = TopoDS::Face(fixedFace);
          if (outLevel > 0)
            printf(" EGADS Internal: Face has been fixed (EG_makeFace)!\n");
        }
      }
    }
    stat = EG_examineFace(Face, 1, &object, outLevel);
    if (stat != EGADS_SUCCESS) return stat;

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Face object (EG_makeFace)!\n");
      return stat;
    }
    obj->oclass = FACE;

    loop = object;
    stat = EG_makeObject(context, &geom);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Surface object (EG_makeFace)!\n");
      obj->oclass = NIL;
      EG_deleteObject(obj);
      return stat;
    }
    geom->topObj = obj;
    EG_completeSurf(geom, hSurface);

  }

  egadsFace *pface   = new egadsFace;
  egObject **loopo   = new egObject*[nl];
  int *senses        = new int[nl];
  loopo[0]           = loop;
  senses[0]          = 1;
  pface->face        = Face;
  pface->surface     = geom;
  pface->nloops      = nl;
  pface->loops       = loopo;
  pface->senses      = senses;
  pface->topFlg      = 0;
  pface->bbox.filled = 0;
  BRepTools::UVBounds(Face, pface->urange[0], pface->urange[1],
                            pface->vrange[0], pface->vrange[1]);
  obj->blind         = pface;
  obj->mtype         = mtype;

  EG_referenceObject(geom, obj);
  EG_referenceObject(loop, obj);
  EG_referenceObject(obj,  context);

  *face = obj;
  return EGADS_SUCCESS;
}


int
EG_getBodyTopos(const egObject *body, /*@null@*/ egObject *src,
                int oclass, int *ntopo, /*@null@*/ egObject ***topos)
{
  int      outLevel, i, n, index;
  egadsMap *map;
  egObject **objs;

  *ntopo = 0;
  if (topos != NULL) *topos = NULL;
  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->oclass != BODY)       return EGADS_NOTBODY;
  if (body->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(body);

  if ((oclass < NODE) || (oclass > SHELL)) {
    if (outLevel > 0)
      printf(" EGADS Error: oclass = %d (EG_getBodyTopos)!\n",
             oclass);
    return EGADS_NOTTOPO;
  }

  egadsBody *pbody = (egadsBody *) body->blind;
  if (oclass == NODE) {
    map = &pbody->nodes;
  } else if (oclass == EDGE) {
    map = &pbody->edges;
  } else if (oclass == LOOP) {
    map = &pbody->loops;
  } else if (oclass == FACE) {
    map = &pbody->faces;
  } else {
    map = &pbody->shells;
  }

  if (src == NULL) {

    n = map->map.Extent();
    if (n == 0) return EGADS_SUCCESS;
    if (topos == NULL) {
      *ntopo = n;
      return EGADS_SUCCESS;
    }
    objs = (egObject **) EG_alloc(n*sizeof(egObject *));
    if (objs == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Malloc oclass = %d, n = %d (EG_getBodyTopos)!\n",
               oclass, n);
      return EGADS_MALLOC;
    }
    for (i = 0; i < n; i++) objs[i] = map->objs[i];

  } else {

    if (src->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: src not an EGO (EG_getBodyTopos)!\n");
      return EGADS_NOTOBJ;
    }
    if ((src->oclass < NODE) || (src->oclass > SHELL)) {
      if (outLevel > 0)
        printf(" EGADS Error: src not a Topo (EG_getBodyTopos)!\n");
      return EGADS_NOTTOPO;
    }
    if (src->oclass == oclass) {
      if (outLevel > 0)
        printf(" EGADS Error: src Topo is oclass (EG_getBodyTopos)!\n");
      return EGADS_TOPOERR;
    }
    if (EG_context(body) != EG_context(src)) {
      if (outLevel > 0)
        printf(" EGADS Error: Context mismatch (EG_getBodyTopos)!\n");
      return EGADS_MIXCNTX;
    }
    if (src->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL src pointer (EG_getBodyTopos)!\n");
      return EGADS_NODATA;
    }

    TopoDS_Shape     shape;
    TopAbs_ShapeEnum senum;
    if (src->oclass == NODE) {
      egadsNode *pnode = (egadsNode *) src->blind;
      shape = pnode->node;
      senum = TopAbs_VERTEX;
    } else if (src->oclass == EDGE) {
      egadsEdge *pedge = (egadsEdge *) src->blind;
      shape = pedge->edge;
      senum = TopAbs_EDGE;
    } else if (src->oclass == LOOP) {
      egadsLoop *ploop = (egadsLoop *) src->blind;
      shape = ploop->loop;
      senum = TopAbs_WIRE;
    } else if (src->oclass == FACE) {
      egadsFace *pface = (egadsFace *) src->blind;
      shape = pface->face;
      senum = TopAbs_FACE;
    } else {
      egadsShell *pshell = (egadsShell *) src->blind;
      shape = pshell->shell;
      senum = TopAbs_SHELL;
    }

    if (src->oclass > oclass) {

      // look down the tree (get sub-shapes)
      if (oclass == NODE) {
        senum = TopAbs_VERTEX;
      } else if (oclass == EDGE) {
        senum = TopAbs_EDGE;
      } else if (oclass == LOOP) {
        senum = TopAbs_WIRE;
      } else if (oclass == FACE) {
        senum = TopAbs_FACE;
      } else {
        senum = TopAbs_SHELL;
      }
      TopTools_IndexedMapOfShape smap;
      TopExp::MapShapes(shape, senum, smap);
      n = smap.Extent();
      if (n == 0) return EGADS_SUCCESS;
      if (topos == NULL) {
        *ntopo = n;
        return EGADS_SUCCESS;
      }
      objs = (egObject **) EG_alloc(n*sizeof(egObject *));
      if (objs == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc oclass = %d, n = %d (EG_getBodyTopos)!\n",
                 oclass, n);
        return EGADS_MALLOC;
      }
      for (i = 0; i < n; i++) {
        objs[i] = NULL;
        TopoDS_Shape shapo = smap(i+1);
        index = map->map.FindIndex(shapo);
        if (index == 0) {
          if (outLevel > 0)
            printf(" EGADS Warning: %d/%d NotFound oclass = %d (EG_getBodyTopos)!\n",
                   i+1, n, oclass);
        } else {
          objs[i] = map->objs[index-1];
        }
      }

    } else {

      // look up (get super-shapes)
      for (n = i = 0; i < map->map.Extent(); i++) {
        TopoDS_Shape shapo = map->map(i+1);
        TopTools_IndexedMapOfShape smap;
        TopExp::MapShapes(shapo, senum, smap);
        if (smap.FindIndex(shape) != 0) n++;
      }
      if (n == 0) return EGADS_SUCCESS;
      if (topos == NULL) {
        *ntopo = n;
        return EGADS_SUCCESS;
      }
      objs = (egObject **) EG_alloc(n*sizeof(egObject *));
      if (objs == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc oclass = %d, N = %d (EG_getBodyTopos)!\n",
                 oclass, n);
        return EGADS_MALLOC;
      }
      for (n = i = 0; i < map->map.Extent(); i++) {
        TopoDS_Shape shapo = map->map(i+1);
        TopTools_IndexedMapOfShape smap;
        TopExp::MapShapes(shapo, senum, smap);
        if (smap.FindIndex(shape) != 0) {
          objs[n] = map->objs[i];
          n++;
        }
      }
    }
  }

  *ntopo = n;
  *topos = objs;

  return EGADS_SUCCESS;
}


int
EG_indexBodyTopo(const egObject *body, const egObject *src)
{
  int outLevel, index;

  if (src  == NULL)               return EGADS_NULLOBJ;
  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->oclass != BODY)       return EGADS_NOTBODY;
  if (body->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(body);

  if (src->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: src not an EGO (EG_indexBodyTopo)!\n");
    return EGADS_NOTOBJ;
  }
  if ((src->oclass < NODE) || (src->oclass > SHELL)) {
    if (outLevel > 0)
      printf(" EGADS Error: src not a Topo (EG_indexBodyTopo)!\n");
    return EGADS_NOTTOPO;
  }
  if (EG_context(body) != EG_context(src)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_indexBodyTopo)!\n");
    return EGADS_MIXCNTX;
  }
  if (src->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL src pointer (EG_indexBodyTopo)!\n");
    return EGADS_NODATA;
  }

  egadsBody *pbody = (egadsBody *) body->blind;
  if (src->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) src->blind;
    index = pbody->nodes.map.FindIndex(pnode->node);
  } else if (src->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) src->blind;
    index = pbody->edges.map.FindIndex(pedge->edge);
  } else if (src->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) src->blind;
    index = pbody->loops.map.FindIndex(ploop->loop);
  } else if (src->oclass == FACE) {
    egadsFace *pface = (egadsFace *) src->blind;
    index = pbody->faces.map.FindIndex(pface->face);
  } else {
    egadsShell *pshell = (egadsShell *) src->blind;
    index = pbody->shells.map.FindIndex(pshell->shell);
  }

  if (index == 0) index = EGADS_NOTFOUND;
  return index;
}


int
EG_objectBodyTopo(const egObject *body, int oclass, int index, egObject **obj)
{
  egadsMap *map;

  if  (body == NULL)                       return EGADS_NULLOBJ;
  if  (body->magicnumber != MAGIC)         return EGADS_NOTOBJ;
  if  (body->oclass != BODY)               return EGADS_NOTBODY;
  if  (body->blind == NULL)                return EGADS_NODATA;
  if ((oclass < NODE) || (oclass > SHELL)) return EGADS_NOTTOPO;
  if  (index <= 0)                         return EGADS_INDEXERR;

  egadsBody *pbody = (egadsBody *) body->blind;
  if (oclass == NODE) {
    map = &pbody->nodes;
  } else if (oclass == EDGE) {
    map = &pbody->edges;
  } else if (oclass == LOOP) {
    map = &pbody->loops;
  } else if (oclass == FACE) {
    map = &pbody->faces;
  } else {
    map = &pbody->shells;
  }
  if (index > map->map.Extent()) return EGADS_INDEXERR;

  *obj = map->objs[index-1];
  return EGADS_SUCCESS;
}


int
EG_sameBodyTopo(const egObject *bod1, const egObject *bod2)
{
  int      i, j, outLevel, ind1, ind2, nc, err = 0;
  egObject *obj1, *obj2;

  if (bod1 == NULL)               return EGADS_NULLOBJ;
  if (bod1->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (bod1->oclass != BODY)       return EGADS_NOTBODY;
  if (bod1->blind == NULL)        return EGADS_NODATA;
  if (bod2 == NULL)               return EGADS_NULLOBJ;
  if (bod2->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (bod2->oclass != BODY)       return EGADS_NOTBODY;
  if (bod2->blind == NULL)        return EGADS_NODATA;
  if (bod1->mtype != bod2->mtype) return EGADS_TOPOERR;
  egadsBody *pbod1 = (egadsBody *) bod1->blind;
  egadsBody *pbod2 = (egadsBody *) bod2->blind;
  outLevel = EG_outLevel(bod1);

  /* sizes must match */
  if ((bod1->mtype == SHEETBODY) || (bod1->mtype == SOLIDBODY))
    if (pbod1->shells.map.Extent() != pbod2->shells.map.Extent()) {
      if (outLevel > 1)
        printf(" EGADS Warning: #Shells = %d %d (EG_sameBodyTopo)!\n",
               pbod1->shells.map.Extent(), pbod2->shells.map.Extent());
      return EGADS_TOPOCNT;
    }

  if (bod1->mtype != WIREBODY)
    if (pbod1->faces.map.Extent() != pbod2->faces.map.Extent()) {
      if (outLevel > 1)
        printf(" EGADS Warning: #Faces = %d %d (EG_sameBodyTopo)!\n",
               pbod1->faces.map.Extent(), pbod2->faces.map.Extent());
      return EGADS_TOPOCNT;
    }

  if (pbod1->loops.map.Extent() != pbod2->loops.map.Extent()) {
    if (outLevel > 1)
      printf(" EGADS Warning: #Loops = %d %d (EG_sameBodyTopo)!\n",
             pbod1->loops.map.Extent(), pbod2->loops.map.Extent());
    return EGADS_TOPOCNT;
  }
  if (pbod1->edges.map.Extent() != pbod2->edges.map.Extent()) {
    if (outLevel > 1)
      printf(" EGADS Warning: #Edges = %d %d (EG_sameBodyTopo)!\n",
             pbod1->edges.map.Extent(), pbod2->edges.map.Extent());
    return EGADS_TOPOCNT;
  }
  if (pbod1->nodes.map.Extent() != pbod2->nodes.map.Extent()) {
    if (outLevel > 1)
      printf(" EGADS Warning: #Nodes = %d %d (EG_sameBodyTopo)!\n",
             pbod1->nodes.map.Extent(), pbod2->nodes.map.Extent());
    return EGADS_TOPOCNT;
  }

  /* look at children */

  for (i = 0; i < pbod1->edges.map.Extent(); i++) {
    obj1 = pbod1->edges.objs[i];
    obj2 = pbod2->edges.objs[i];
    if (obj1->mtype != obj2->mtype) {
      if (outLevel > 1)
        printf(" EGADS Warning: Edge %d - Types = %d %d (EG_sameBodyTopo)!\n",
               i+1, obj1->mtype, obj2->mtype);
      err++;
      continue;
    }
    egadsEdge *pedg1 = (egadsEdge *) obj1->blind;
    egadsEdge *pedg2 = (egadsEdge *) obj2->blind;
                                nc = 1;
    if (obj1->mtype == TWONODE) nc = 2;
    for (j = 0; j < nc; j++) {
      egObject  *nod1  = pedg1->nodes[j];
      egObject  *nod2  = pedg2->nodes[j];
      egadsNode *pnod1 = (egadsNode *) nod1->blind;
      egadsNode *pnod2 = (egadsNode *) nod2->blind;
      ind1 = pbod1->nodes.map.FindIndex(pnod1->node);
      ind2 = pbod2->nodes.map.FindIndex(pnod2->node);
      if (ind1 != ind2) {
        if (outLevel > 1)
          printf(" EGADS Warning: Edge %d - nIndices = %d %d (EG_sameBodyTopo)!\n",
                 i+1, ind1, ind2);
        err++;
      }
    }
  }
  if (err != 0) return EGADS_TOPOERR;

  for (i = 0; i < pbod1->loops.map.Extent(); i++) {
    obj1 = pbod1->loops.objs[i];
    obj2 = pbod2->loops.objs[i];
    if (obj1->mtype != obj2->mtype) {
      if (outLevel > 1)
        printf(" EGADS Warning: Loop %d - Types = %d %d (EG_sameBodyTopo)!\n",
               i+1, obj1->mtype, obj2->mtype);
      return EGADS_TOPOERR;
    }
    egadsLoop *ploo1 = (egadsLoop *) obj1->blind;
    egadsLoop *ploo2 = (egadsLoop *) obj2->blind;
    nc = ploo1->nedges;
    if (nc != ploo2->nedges) {
      if (outLevel > 1)
        printf(" EGADS Warning: Loop %d - nEdges = %d %d (EG_sameBodyTopo)!\n",
               i+1, nc, ploo2->nedges);
      return EGADS_TOPOERR;
    }
    for (j = 0; j < nc; j++) {
      egObject  *edg1  = ploo1->edges[j];
      egObject  *edg2  = ploo2->edges[j];
      egadsEdge *pedg1 = (egadsEdge *) edg1->blind;
      egadsEdge *pedg2 = (egadsEdge *) edg2->blind;
      ind1 = pbod1->edges.map.FindIndex(pedg1->edge);
      ind2 = pbod2->edges.map.FindIndex(pedg2->edge);
      if (ind1 != ind2) {
        if (outLevel > 1)
          printf(" EGADS Warning: Loop %d - eIndices = %d %d (EG_sameBodyTopo)!\n",
                 i+1, ind1, ind2);
        err++;
      }
    }
  }
  if (err != 0) return EGADS_TOPOERR;
  if (bod1->mtype == WIREBODY) return EGADS_SUCCESS;

  for (i = 0; i < pbod1->faces.map.Extent(); i++) {
    obj1 = pbod1->faces.objs[i];
    obj2 = pbod2->faces.objs[i];
    if (obj1->mtype != obj2->mtype) {
      if (outLevel > 1)
        printf(" EGADS Warning: Face %d - Types = %d %d (EG_sameBodyTopo)!\n",
               i+1, obj1->mtype, obj2->mtype);
      return EGADS_TOPOERR;
    }
    egadsFace *pfac1 = (egadsFace *) obj1->blind;
    egadsFace *pfac2 = (egadsFace *) obj2->blind;
    nc = pfac1->nloops;
    if (nc != pfac2->nloops) {
      if (outLevel > 1)
        printf(" EGADS Warning: Face %d - nLoops = %d %d (EG_sameBodyTopo)!\n",
               i+1, nc, pfac2->nloops);
      return EGADS_TOPOERR;
    }
    for (j = 0; j < nc; j++) {
      egObject  *loo1  = pfac1->loops[j];
      egObject  *loo2  = pfac2->loops[j];
      egadsLoop *ploo1 = (egadsLoop *) loo1->blind;
      egadsLoop *ploo2 = (egadsLoop *) loo2->blind;
      ind1 = pbod1->loops.map.FindIndex(ploo1->loop);
      ind2 = pbod2->loops.map.FindIndex(ploo2->loop);
      if (ind1 != ind2) {
        if (outLevel > 1)
          printf(" EGADS Warning: Face %d - lIndices = %d %d (EG_sameBodyTopo)!\n",
                 i+1, ind1, ind2);
        err++;
      }
    }
  }
  if (err != 0) return EGADS_TOPOERR;
  if (bod1->mtype == FACEBODY) return EGADS_SUCCESS;

  for (i = 0; i < pbod1->shells.map.Extent(); i++) {
    obj1 = pbod1->shells.objs[i];
    obj2 = pbod2->shells.objs[i];
    if (obj1->mtype != obj2->mtype) {
      if (outLevel > 1)
        printf(" EGADS Warning: Shell %d - Types = %d %d (EG_sameBodyTopo)!\n",
               i+1, obj1->mtype, obj2->mtype);
      return EGADS_TOPOERR;
    }
    egadsShell *pshl1 = (egadsShell *) obj1->blind;
    egadsShell *pshl2 = (egadsShell *) obj2->blind;
    nc = pshl1->nfaces;
    if (nc != pshl2->nfaces) {
      if (outLevel > 1)
        printf(" EGADS Warning: Shell %d - nFaces = %d %d (EG_sameBodyTopo)!\n",
               i+1, nc, pshl2->nfaces);
      return EGADS_TOPOERR;
    }
    for (j = 0; j < nc; j++) {
      egObject  *fac1  = pshl1->faces[j];
      egObject  *fac2  = pshl2->faces[j];
      egadsFace *pfac1 = (egadsFace *) fac1->blind;
      egadsFace *pfac2 = (egadsFace *) fac2->blind;
      ind1 = pbod1->faces.map.FindIndex(pfac1->face);
      ind2 = pbod2->faces.map.FindIndex(pfac2->face);
      if (ind1 != ind2) {
        if (outLevel > 1)
          printf(" EGADS Warning: Shell %d - fIndices = %d %d (EG_sameBodyTopo)!\n",
                 i+1, ind1, ind2);
        err++;
      }
    }
  }
  if (err != 0) return EGADS_TOPOERR;

  return EGADS_SUCCESS;
}


int
EG_makeSolidBody(egObject *context, int stypx, const double *data,
                 egObject **body)
{
  int           outLevel, stype, stat, nerr;
  egObject      *obj;
  TopoDS_Shape  solid;
#ifdef OCC_SOLIDS
  Standard_Real height;
#endif

  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (EG_sameThread(context))        return EGADS_CNTXTHRD;
  if (data == NULL)                  return EGADS_NODATA;
  outLevel = EG_outLevel(context);
  stype    = abs(stypx);

  if ((stype < BOX) || (stype > TORUS)) {
    if (outLevel > 0)
      printf(" EGADS Error: stype = %d (EG_makeSolidBody)!\n", stype);
    return EGADS_RANGERR;
  }

  switch (stype) {

  /* box = 1 */
  case BOX:
#ifdef OCC_SOLIDS
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeBox(gp_Pnt(data[0], data[1], data[2]),
                                data[3], data[4], data[5]);
    break;
#else
    return EG_makeSolidBox(context, data, body);
#endif

  /* sphere = 2 */
  case SPHERE:
#ifdef OCC_SOLIDS
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeSphere(gp_Pnt(data[0], data[1], data[2]),
                                   data[3]);
    break;
#else
    return EG_makeSolidSphere(context, stypx, data, body);
#endif

  /* cone = 3 */
  case CONE:
#ifdef OCC_SOLIDS
    height = sqrt( (data[3]-data[0])*(data[3]-data[0]) +
                   (data[4]-data[1])*(data[4]-data[1]) +
                   (data[5]-data[2])*(data[5]-data[2]) );
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(data[0], data[1], data[2]),
                                 gp_Dir(data[3]-data[0], data[4]-data[1],
                                        data[5]-data[2])),
                                 0.0, data[6], height);
    break;
#else
    return EG_makeSolidCone(context, stypx, data, body);
#endif

  /* cylinder = 4 */
  case CYLINDER:
#ifdef OCC_SOLIDS
    height = sqrt( (data[3]-data[0])*(data[3]-data[0]) +
                   (data[4]-data[1])*(data[4]-data[1]) +
                   (data[5]-data[2])*(data[5]-data[2]) );
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(data[0], data[1], data[2]),
                                     gp_Dir(data[3]-data[0], data[4]-data[1],
                                            data[5]-data[2])),
                                     data[6], height);
    break;
#else
    return EG_makeSolidCylinder(context, stypx, data, body);
#endif

  /* torus = 5 */
  case TORUS:
#ifdef OCC_SOLIDS
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(data[0], data[1], data[2]),
                                  gp_Dir(data[3], data[4], data[5])),
                                  data[6], data[7]);
#else
    return EG_makeSolidTorus(context, stypx, data, body);
#endif
    break;
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_makeSolidBody)!\n");
    return stat;
  }
  egadsBody *pbody   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = SOLIDBODY;
  pbody->nodes.objs  = NULL;
  pbody->edges.objs  = NULL;
  pbody->loops.objs  = NULL;
  pbody->faces.objs  = NULL;
  pbody->shells.objs = NULL;
  pbody->senses      = NULL;
  pbody->shape       = solid;
  pbody->bbox.filled = 0;
  pbody->massFill    = 0;
  obj->blind         = pbody;
  if (stypx > BOX) EG_splitPeriodics(pbody);
  stat = EG_traverseBody(context, 0, obj, obj, pbody, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbody;
    return stat;
  }

  EG_referenceObject(obj, context);
  *body = obj;
  return EGADS_SUCCESS;
}


int
EG_makeSolidBody_dot(egObject *body, int stypx,
                     const double *data, const double *data_dot)
{
  int           outLevel, stype;

  if (body == NULL)                  return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (body->oclass != BODY)          return EGADS_NOTBODY;
  if (body->mtype  != SOLIDBODY)     return EGADS_NOTBODY;
  if (data == NULL)                  return EGADS_NODATA;
  if (data_dot == NULL)              return EGADS_NODATA;
  if (EG_sameThread(body))           return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(body);
  stype    = abs(stypx);

  if ((stype < BOX) || (stype > TORUS)) {
    if (outLevel > 0)
      printf(" EGADS Error: stype = %d (EG_makeSolidBody_dot)!\n", stype);
    return EGADS_RANGERR;
  }

  switch (stype) {

  /* box = 1 */
  case BOX:
    return EG_makeSolidBox_dot(data, data_dot, body);

  /* sphere = 2 */
  case SPHERE:
    return EG_makeSolidSphere_dot(stypx, data, data_dot, body);

  /* cone = 3 */
  case CONE:
    return EG_makeSolidCone_dot(stypx, data, data_dot, body);

  /* cylinder = 4 */
  case CYLINDER:
    return EG_makeSolidCylinder_dot(stypx, data, data_dot, body);

  /* torus = 5 */
  case TORUS:
    return EG_makeSolidTorus_dot(stypx, data, data_dot, body);
  }

  if (outLevel > 0)
    printf(" EGADS Error: stype = %d (EG_makeSolidBody_dot)!\n", stype);
  return EGADS_RANGERR;
}


int
EG_getBoundingBox(const egObject *topo, double *bbox)
{
  int          i, n;
  egObject     *obj;
  Bnd_Box      Box;
  egadsBox     *ebox;
  TopoDS_Shape shape;

  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->oclass < NODE)        return EGADS_NOTTOPO;
  if (topo->oclass > MODEL)       return EGADS_NOTTOPO;
  if (topo->blind == NULL)        return EGADS_NODATA;

  /* are we cached? */
  if (topo->oclass == NODE) {

    egadsNode *pnode = (egadsNode *) topo->blind;
    ebox  = &pnode->bbox;
    shape = pnode->node;

  } else if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    ebox  = &pedge->bbox;
    shape = pedge->edge;

  } else if (topo->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) topo->blind;
    ebox  = &ploop->bbox;
    shape = ploop->loop;

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    ebox  = &pface->bbox;
    shape = pface->face;

  } else if (topo->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    ebox  = &pshell->bbox;
    shape = pshell->shell;

  } else if (topo->oclass == BODY) {

    egadsBody *pbody = (egadsBody *) topo->blind;
    ebox  = &pbody->bbox;
    shape = pbody->shape;

  } else {

    egadsModel *pmodel = (egadsModel *) topo->blind;
    ebox  = &pmodel->bbox;
    shape = pmodel->shape;

  }
  if (ebox->filled != 0) {
    bbox[0] = ebox->box[0];
    bbox[1] = ebox->box[1];
    bbox[2] = ebox->box[2];
    bbox[3] = ebox->box[3];
    bbox[4] = ebox->box[4];
    bbox[5] = ebox->box[5];
    return EGADS_SUCCESS;
  }

  /* do we have a single degenerate Edge? */
  TopTools_IndexedMapOfShape nmap;
  TopExp::MapShapes(shape, TopAbs_VERTEX, nmap);
  if (nmap.Extent() == 1) {
    TopTools_IndexedMapOfShape emap;
    TopExp::MapShapes(shape, TopAbs_EDGE, emap);
    if (emap.Extent() == 1) {
      TopoDS_Edge edge = TopoDS::Edge(emap(1));
      if (BRep_Tool::Degenerated(edge)) {
        TopoDS_Vertex vert = TopoDS::Vertex(nmap(1));
        gp_Pnt pv = BRep_Tool::Pnt(vert);
        /* fill the cache */
        ebox->filled = 1;
        ebox->box[0] = pv.X();
        ebox->box[1] = pv.Y();
        ebox->box[2] = pv.Z();
        ebox->box[3] = pv.X();
        ebox->box[4] = pv.Y();
        ebox->box[5] = pv.Z();
        /* fill the return value */
        bbox[0] = ebox->box[0];
        bbox[1] = ebox->box[1];
        bbox[2] = ebox->box[2];
        bbox[3] = ebox->box[3];
        bbox[4] = ebox->box[4];
        bbox[5] = ebox->box[5];
        return EGADS_SUCCESS;
      }
    }
  }

  /* no -- lets compute the bounding box */
  try {

    if (topo->oclass == NODE) {

      egadsNode *pnode = (egadsNode *) topo->blind;
#if CASVER >= 730
      BRepBndLib::AddOptimal(pnode->node, Box, Standard_False);
#else
      BRepBndLib::Add(pnode->node, Box);
#endif

    } else if (topo->oclass == EDGE) {

      egadsEdge *pedge = (egadsEdge *) topo->blind;
#if CASVER >= 730
      BRepBndLib::AddOptimal(pedge->edge, Box, Standard_False);
#else
      BRepBndLib::Add(pedge->edge, Box);
#endif

    } else if (topo->oclass == LOOP) {

      egadsLoop *ploop = (egadsLoop *) topo->blind;
#if CASVER >= 730
      BRepBndLib::AddOptimal(ploop->loop, Box, Standard_False);
#else
      BRepBndLib::Add(ploop->loop, Box);
#endif

    } else if (topo->oclass == FACE) {

      egadsFace *pface = (egadsFace *) topo->blind;
#if CASVER >= 730
      BRepBndLib::AddOptimal(pface->face, Box, Standard_False);
#else
      BRepBndLib::AddClose(pface->face, Box);
#endif

    } else if (topo->oclass == SHELL) {

      egadsShell *pshell = (egadsShell *) topo->blind;
#if CASVER >= 730
      BRepBndLib::AddOptimal(pshell->shell, Box, Standard_False);
#else
      BRepBndLib::Add(pshell->shell, Box);
#endif

    } else if (topo->oclass == BODY) {

      egadsBody *pbody = (egadsBody *) topo->blind;
#if CASVER >= 730
      BRepBndLib::AddOptimal(pbody->shape, Box, Standard_False);
#else
      BRepBndLib::Add(pbody->shape, Box);
#endif

    } else {

      egadsModel *pmodel = (egadsModel *) topo->blind;
      n = pmodel->nbody;
      for (i = 0; i < n; i++) {
        obj = pmodel->bodies[i];
        if (obj == NULL) continue;
        egadsBody *pbody = (egadsBody *) obj->blind;
        if (pbody == NULL) continue;
 #if CASVER >= 730
        BRepBndLib::AddOptimal(pbody->shape, Box, Standard_False);
#else
        BRepBndLib::Add(pbody->shape, Box);
#endif
      }

    }

    Box.Get(bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Warning: BoundingBox failure (EG_getBoundingBox)!\n");
    printf("                %s\n", e.GetMessageString());
    return EGADS_TOPOERR;
  }
  catch (...) {
    printf(" EGADS Warning: BoundingBox failure (EG_getBoundingBox)!\n");
    return EGADS_TOPOERR;
  }

  /* save it away */
  ebox->filled = 1;
  ebox->box[0] = bbox[0];
  ebox->box[1] = bbox[1];
  ebox->box[2] = bbox[2];
  ebox->box[3] = bbox[3];
  ebox->box[4] = bbox[4];
  ebox->box[5] = bbox[5];

  return EGADS_SUCCESS;
}


int
EG_getMassProperties(const egObject *topo, /*@null@*/ double *data)
{
  int           i;
  gp_Pnt        CofG, pv;
  gp_Mat        Inert;
  BRepGProp     BProps;
  GProp_GProps  SProps, VProps;
  TopoDS_Shape  shape;
  TopoDS_Vertex vert;

  if  (topo == NULL)               return EGADS_NULLOBJ;
  if  (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((topo->oclass < EDGE) ||
      (topo->oclass > BODY))       return EGADS_NOTTOPO;
  if  (topo->blind == NULL)        return EGADS_NODATA;
  if ((topo->oclass == BODY) && (data == NULL)) {
    egadsBody *pbody = (egadsBody *) topo->blind;
    BRepTools::Clean(pbody->shape);
    return EGADS_SUCCESS;
  }
  if  (data == NULL)               return EGADS_NONAME;
  for (i = 0; i < 14; i++) data[i] = 0.0;

  /* are we degenerate (a single Node)? */
  if (topo->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) topo->blind;
    shape = pedge->edge;
  } else if (topo->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) topo->blind;
    shape = ploop->loop;
  } else if (topo->oclass == FACE) {
    egadsFace *pface = (egadsFace *) topo->blind;
    shape = pface->face;
  } else if (topo->oclass == SHELL) {
    egadsShell *pshell = (egadsShell *) topo->blind;
    shape = pshell->shell;
  } else {
    egadsBody *pbody = (egadsBody *) topo->blind;
    shape = pbody->shape;
    if (pbody->massFill != 0) {
      for (i = 0; i < 14; i++) data[i] = pbody->massProp[i];
      return EGADS_SUCCESS;
    }
  }
  TopTools_IndexedMapOfShape nmap;
  TopExp::MapShapes(shape, TopAbs_VERTEX, nmap);
  if (nmap.Extent() == 1) {
    TopTools_IndexedMapOfShape emap;
    TopExp::MapShapes(shape, TopAbs_EDGE, emap);
    if (emap.Extent() == 1) {
      TopoDS_Edge edge = TopoDS::Edge(emap(1));
      if (BRep_Tool::Degenerated(edge)) {
        vert    = TopoDS::Vertex(nmap(1));
        pv      = BRep_Tool::Pnt(vert);
        data[2] = pv.X();
        data[3] = pv.Y();
        data[4] = pv.Z();
        return EGADS_SUCCESS;
      }
    }
  }

  /* use the appropriate dimensional methods */
  if ((topo->oclass == EDGE) || (topo->oclass == LOOP)) {

    BProps.LinearProperties(shape, SProps);
    CofG  = SProps.CentreOfMass();
    Inert = SProps.MatrixOfInertia();

  } else if ((topo->oclass == FACE) || (topo->oclass == SHELL)) {

    BProps.SurfaceProperties(shape, SProps);
    CofG  = SProps.CentreOfMass();
    Inert = SProps.MatrixOfInertia();

  } else {

    if (topo->mtype == SOLIDBODY) {
      BProps.SurfaceProperties(shape, SProps);
      BProps.VolumeProperties( shape, VProps);
      CofG    = VProps.CentreOfMass();
      Inert   = VProps.MatrixOfInertia();
      data[0] = VProps.Mass();
    } else if (topo->mtype == WIREBODY) {
      BProps.LinearProperties(shape, SProps);
      CofG    = SProps.CentreOfMass();
      Inert   = SProps.MatrixOfInertia();
    } else {
      BProps.SurfaceProperties(shape, SProps);
      CofG    = SProps.CentreOfMass();
      Inert   = SProps.MatrixOfInertia();
    }

  }

  data[ 1] = SProps.Mass();
  data[ 2] = CofG.X();
  data[ 3] = CofG.Y();
  data[ 4] = CofG.Z();
  data[ 5] = Inert.Value(1,1);
  data[ 6] = Inert.Value(1,2);
  data[ 7] = Inert.Value(1,3);
  data[ 8] = Inert.Value(2,1);
  data[ 9] = Inert.Value(2,2);
  data[10] = Inert.Value(2,3);
  data[11] = Inert.Value(3,1);
  data[12] = Inert.Value(3,2);
  data[13] = Inert.Value(3,3);
  
  if (topo->oclass == BODY) {
    egadsBody *pbody = (egadsBody *) topo->blind;
    for (i = 0; i < 14; i++) pbody->massProp[i] = data[i];
    pbody->massFill = 1;
  }

  return EGADS_SUCCESS;
}


int
EG_isEquivalent(const egObject *topo1, const egObject *topo2)
{
  int          i, stat;
  TopoDS_Shape shape1, shape2;

  if (topo1 == topo2)                 return EGADS_SUCCESS;
  if (topo1 == NULL)                  return EGADS_NULLOBJ;
  if (topo1->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (topo1->oclass < NODE)           return EGADS_NOTTOPO;
  if (topo1->blind == NULL)           return EGADS_NODATA;
  if (topo2 == NULL)                  return EGADS_NULLOBJ;
  if (topo2->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (topo2->oclass != topo1->oclass) return EGADS_NOTTOPO;
  if (topo2->blind == NULL)           return EGADS_NODATA;

  if (topo1->oclass == NODE) {

    egadsNode *pnode1 = (egadsNode *) topo1->blind;
    egadsNode *pnode2 = (egadsNode *) topo2->blind;
    shape1            = pnode1->node;
    shape2            = pnode2->node;

  } else if (topo1->oclass == EDGE) {

    egadsEdge *pedge1 = (egadsEdge *) topo1->blind;
    egadsEdge *pedge2 = (egadsEdge *) topo2->blind;
    shape1            = pedge1->edge;
    shape2            = pedge2->edge;

  } else if (topo1->oclass == LOOP) {

    egadsLoop *ploop1 = (egadsLoop *) topo1->blind;
    egadsLoop *ploop2 = (egadsLoop *) topo2->blind;
    shape1            = ploop1->loop;
    shape2            = ploop2->loop;

  } else if (topo1->oclass == FACE) {

    egadsFace *pface1 = (egadsFace *) topo1->blind;
    egadsFace *pface2 = (egadsFace *) topo2->blind;
    shape1            = pface1->face;
    shape2            = pface2->face;

  } else if (topo1->oclass == SHELL) {

    egadsShell *pshell1 = (egadsShell *) topo1->blind;
    egadsShell *pshell2 = (egadsShell *) topo2->blind;
    shape1              = pshell1->shell;
    shape2              = pshell2->shell;

  } else if (topo1->oclass == BODY) {

    egadsBody *pbody1 = (egadsBody *) topo1->blind;
    egadsBody *pbody2 = (egadsBody *) topo2->blind;
    shape1            = pbody1->shape;
    shape2            = pbody2->shape;

  } else {

    egadsModel *pmodel1 = (egadsModel *) topo1->blind;
    egadsModel *pmodel2 = (egadsModel *) topo2->blind;
    shape1              = pmodel1->shape;
    shape2              = pmodel2->shape;
  }
  if (shape1.IsSame(shape2)) return EGADS_SUCCESS;

  /* try to match via geometry */

  if (topo1->oclass == NODE) {

    return EG_isSame(topo1, topo2);

  } else if (topo1->oclass == EDGE) {

    stat = EG_isSame(topo1, topo2);
    if (stat != EGADS_SUCCESS) return stat;
    egadsEdge *pedge1 = (egadsEdge *) topo1->blind;
    egadsEdge *pedge2 = (egadsEdge *) topo2->blind;
    stat = EG_isSame(pedge1->nodes[0], pedge2->nodes[0]);
    if (stat != EGADS_SUCCESS) return stat;
    return EG_isSame(pedge1->nodes[1], pedge2->nodes[1]);

  } else if (topo1->oclass == LOOP) {

    egadsLoop *ploop1 = (egadsLoop *) topo1->blind;
    egadsLoop *ploop2 = (egadsLoop *) topo2->blind;

    if (ploop1->nedges != ploop2->nedges) return EGADS_OUTSIDE;
    for (i = 0; i < ploop1->nedges; i++) {
      if (ploop1->senses[i] != ploop2->senses[i]) return EGADS_OUTSIDE;
      stat = EG_isEquivalent(ploop1->edges[i], ploop2->edges[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
    return EGADS_SUCCESS;

  } else if (topo1->oclass == FACE) {

    egadsFace *pface1 = (egadsFace *) topo1->blind;
    egadsFace *pface2 = (egadsFace *) topo2->blind;

    stat = EG_isSame(topo1, topo2);
    if (stat != EGADS_SUCCESS) return stat;
    if (pface1->nloops != pface2->nloops) return EGADS_OUTSIDE;
    for (i = 0; i < pface1->nloops; i++) {
      stat = EG_isEquivalent(pface1->loops[i], pface2->loops[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
    return EGADS_SUCCESS;

  } else if (topo1->oclass == SHELL) {

    egadsShell *pshell1 = (egadsShell *) topo1->blind;
    egadsShell *pshell2 = (egadsShell *) topo2->blind;

    if (pshell1->nfaces != pshell2->nfaces) return EGADS_OUTSIDE;
    for (i = 0; i < pshell1->nfaces; i++) {
      stat = EG_isEquivalent(pshell1->faces[i], pshell2->faces[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
    return EGADS_SUCCESS;

  } else if (topo1->oclass == BODY) {

    egadsBody *pbody1 = (egadsBody *) topo1->blind;
    egadsBody *pbody2 = (egadsBody *) topo2->blind;

    if (pbody1->shells.map.Extent() != pbody2->shells.map.Extent())
      return EGADS_OUTSIDE;
    if (pbody1->faces.map.Extent()  != pbody2->faces.map.Extent())
      return EGADS_OUTSIDE;
    if (pbody1->loops.map.Extent()  != pbody2->loops.map.Extent())
      return EGADS_OUTSIDE;
    if (pbody1->edges.map.Extent()  != pbody2->edges.map.Extent())
      return EGADS_OUTSIDE;
    if (pbody1->nodes.map.Extent()  != pbody2->nodes.map.Extent())
      return EGADS_OUTSIDE;

    if (pbody1->shells.map.Extent() != 0) {
      for (i = 0; i < pbody1->shells.map.Extent(); i++) {
        stat = EG_isEquivalent(pbody1->shells.objs[i], pbody2->shells.objs[i]);
        if (stat != EGADS_SUCCESS) return stat;
      }
    } else if (pbody1->faces.map.Extent() != 0) {
      for (i = 0; i < pbody1->faces.map.Extent(); i++) {
        stat = EG_isEquivalent(pbody1->faces.objs[i],  pbody2->faces.objs[i]);
        if (stat != EGADS_SUCCESS) return stat;
      }
    } else {
      for (i = 0; i < pbody1->edges.map.Extent(); i++) {
        stat = EG_isEquivalent(pbody1->edges.objs[i],  pbody2->edges.objs[i]);
        if (stat != EGADS_SUCCESS) return stat;
      }
    }
    return EGADS_SUCCESS;

  }

  return EGADS_OUTSIDE;
}


int
EG_isPlanar(const egObject *topo)
{
  TopoDS_Shape shape;

  if (topo == NULL)                  return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (topo->oclass <= NODE)          return EGADS_NOTTOPO;
  if (topo->blind == NULL)           return EGADS_NODATA;

  if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    shape            = pedge->edge;

  } else if (topo->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) topo->blind;
    shape            = ploop->loop;

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    shape            = pface->face;

  } else if (topo->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    shape              = pshell->shell;

  } else if (topo->oclass == BODY) {

    egadsBody *pbody = (egadsBody *) topo->blind;
    shape            = pbody->shape;

  } else {

    egadsModel *pmodel = (egadsModel *) topo->blind;
    shape              = pmodel->shape;

  }

  // note: single linear edges return not planar!
  BRepBuilderAPI_FindPlane planar(shape);
  if (planar.Found()) return EGADS_SUCCESS;
  return EGADS_OUTSIDE;
}


int
EG_getEdgeUV(const egObject *face, const egObject *topo, int sense, double t,
             double *uv)
{
  int             outLevel, found, stat;
  double          result[6];
  egObject        *pcurv = NULL;
  gp_Pnt2d        P2d;
  TopoDS_Edge     edge;
  TopExp_Explorer ExpW;

  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(face);

  if (topo == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Reference (EG_getEdgeUV)!\n");
    return EGADS_NULLOBJ;
  }
  if (topo->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: topo not an EGO (EG_getEdgeUV)!\n");
    return EGADS_NOTOBJ;
  }
  if ((topo->oclass != EDGE) && (topo->oclass != NODE)) {
    if (outLevel > 0)
      printf(" EGADS Error: Not an Edge (EG_getEdgeUV)!\n");
    return EGADS_NOTTOPO;
  }
  if ((sense < -2) || (sense > 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: Sense = %d (EG_getEdgeUV)!\n", sense);
    return EGADS_RANGERR;
  }
  if (EG_context(face) != EG_context(topo)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_getEdgeUV)!\n");
    return EGADS_MIXCNTX;
  }
  egadsFace *pface = (egadsFace *) face->blind;
  if (pface == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face pointer (EG_getEdgeUV)!\n");
    return EGADS_NODATA;
  }

  // undocumented -- topo is a node!
  if (topo->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Node pointer (EG_getEdgeUV)!\n");
      return EGADS_NODATA;
    }
    TopoDS_Vertex vert = TopoDS::Vertex(pnode->node);

    double tt = 0.0;
    found     = 0;
    for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Wire wire = TopoDS::Wire(ExpW.Current());
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
        edge = TopoDS::Edge(ExpWE.Current());
        if (BRep_Tool::Degenerated(edge)) continue;
        TopoDS_Vertex V1, V2;
        if (edge.Orientation() == TopAbs_REVERSED) {
          TopExp::Vertices(edge, V2, V1, Standard_True);
        } else {
          TopExp::Vertices(edge, V1, V2, Standard_True);
        }
        if ((vert.IsSame(V1)) || (vert.IsSame(V2))) {
          double t1, t2;
          Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
          if (vert.IsSame(V1)) {
            tt = t1;
          } else {
            tt = t2;
          }
          found = 1;
          break;
        }
      }
      if (found != 0) break;
    }
    if (found == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Node not in Face (EG_getEdgeUV)!\n");
      return EGADS_NOTFOUND;
    }

    BRepAdaptor_Curve2d Curve2d(edge, pface->face);
    Curve2d.D0(tt, P2d);
    uv[0] = P2d.X();
    uv[1] = P2d.Y();

    return EGADS_SUCCESS;
  }

  // topo is an edge
  egadsEdge *pedge = (egadsEdge *) topo->blind;
  if (pedge == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Edge pointer (EG_getEdgeUV)!\n");
    return EGADS_NODATA;
  }
  edge  = pedge->edge;
  found = 0;
  for (int iloop = 0; iloop < pface->nloops; iloop++) {
    if (pface->loops[iloop] == NULL) continue;
    egadsLoop *ploop = (egadsLoop *) pface->loops[iloop]->blind;
    if (ploop == NULL) continue;
    for (int iedge = 0; iedge < ploop->nedges; iedge++)
      if (ploop->edges[iedge] == topo) {
        if (ploop->surface != NULL) pcurv = ploop->edges[iedge+ploop->nedges];
        found++;
      }
  }
  if (found == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge/Sense not in Face (EG_getEdgeUV)!\n");
    return EGADS_NOTFOUND;
  }

  // check if we can use sense of zero
  if (sense == 0)
    if (found > 1) {
      if (outLevel > 1)
        printf(" EGADS Warning: Edge in Face twice & sense=0 (EG_getEdgeUV)!\n");
      return EGADS_TOPOERR;
    }

  // we are only in the Face once
  if (found == 1)
    if (pcurv != NULL) {
      stat = EG_evaluate(pcurv, &t, result);
      if (stat == EGADS_SUCCESS) {
        uv[0] = result[0];
        uv[1] = result[1];
      }
      return stat;
    } else {
      BRepAdaptor_Curve2d Curve2d(edge, pface->face);
      Curve2d.D0(t, P2d);
      uv[0] = P2d.X();
      uv[1] = P2d.Y();
      return EGADS_SUCCESS;
    }

  // the Edge in the Face more than once -- find it!
  found = 0;
  for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
    TopoDS_Shape shapw = ExpW.Current();
    TopoDS_Wire  wire  = TopoDS::Wire(shapw);
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shape = ExpWE.Current();
      TopoDS_Edge  wedge = TopoDS::Edge(shape);
      if (!wedge.IsSame(edge)) continue;
      if (sense < 0) {
        if (shape.Orientation() == TopAbs_REVERSED) found++;
      } else {
        if (shape.Orientation() != TopAbs_REVERSED) found++;
      }
      if (found == 0) continue;
      edge = wedge;
      break;
    }
    if (found != 0) break;
  }

  if (!BRep_Tool::SameRange(edge)) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge & PCurve not SameRange (EG_getEdgeUV)!\n");
    return EGADS_GEOMERR;
  }

  // get and evaluate the pcurve
  BRepAdaptor_Curve2d Curve2d(edge, pface->face);
  Curve2d.D0(t, P2d);
  uv[0] = P2d.X();
  uv[1] = P2d.Y();

  return EGADS_SUCCESS;
}


int
EG_getEdgeUVs(const egObject *face, const egObject *edge, int sense, int nts,
              const double *t, double *uv)
{
  int             i, outLevel, found;
  double          result[6];
  egObject        *pcurv = NULL;
  gp_Pnt2d        P2d;
  TopoDS_Edge     edgo;
  TopExp_Explorer ExpW;

  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(face);

  if (edge == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Reference (EG_getEdgeUVs)!\n");
    return EGADS_NULLOBJ;
  }
  if (edge->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: topo not an EGO (EG_getEdgeUVs)!\n");
    return EGADS_NOTOBJ;
  }
  if (edge->oclass != EDGE) {
    if (outLevel > 0)
      printf(" EGADS Error: Not an Edge (EG_getEdgeUVs)!\n");
    return EGADS_NOTTOPO;
  }
  if ((sense < -2) || (sense > 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: Sense = %d (EG_getEdgeUVs)!\n", sense);
    return EGADS_RANGERR;
  }
  if (EG_context(face) != EG_context(edge)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_getEdgeUVs)!\n");
    return EGADS_MIXCNTX;
  }
  egadsFace *pface = (egadsFace *) face->blind;
  if (pface == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face pointer (EG_getEdgeUVs)!\n");
    return EGADS_NODATA;
  }

  egadsEdge *pedge = (egadsEdge *) edge->blind;
  if (pedge == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Edge pointer (EG_getEdgeUVs)!\n");
    return EGADS_NODATA;
  }
  edgo  = pedge->edge;
  found = 0;
  for (int iloop = 0; iloop < pface->nloops; iloop++) {
    if (pface->loops[iloop] == NULL) continue;
    egadsLoop *ploop = (egadsLoop *) pface->loops[iloop]->blind;
    if (ploop == NULL) continue;
    for (int iedge = 0; iedge < ploop->nedges; iedge++)
      if (ploop->edges[iedge] == edge) {
        if (ploop->surface != NULL) pcurv = ploop->edges[iedge+ploop->nedges];
        found++;
      }
  }
  if (found == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge/Sense not in Face (EG_getEdgeUVs)!\n");
    return EGADS_NOTFOUND;
  }

  // check if we can use sense of zero
  if (sense == 0)
    if (found > 1) {
      if (outLevel > 1)
        printf(" EGADS Warning: Edge in Face twice & sense=0 (EG_getEdgeUVs)!\n");
      return EGADS_TOPOERR;
    }

  // we are only in the Face once
  if (found == 1) {
    if (pcurv != NULL) {
      for (i = 0; i < nts; i++) {
        result[0] = result[1] = 0.0;
        EG_evaluate(pcurv, &t[i], result);
        uv[2*i  ] = result[0];
        uv[2*i+1] = result[1];
      }
    } else {
      BRepAdaptor_Curve2d Curve2d(edgo, pface->face);
      for (i = 0; i < nts; i++) {
        Curve2d.D0(t[i], P2d);
        uv[2*i  ] = P2d.X();
        uv[2*i+1] = P2d.Y();
      }
    }
    return EGADS_SUCCESS;
  }

  // the Edge in the Face more than once -- find it!
  found = 0;
  for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
    TopoDS_Shape shapw = ExpW.Current();
    TopoDS_Wire  wire  = TopoDS::Wire(shapw);
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shape = ExpWE.Current();
      TopoDS_Edge  wedge = TopoDS::Edge(shape);
      if (!wedge.IsSame(edgo)) continue;
      if (sense < 0) {
        if (shape.Orientation() == TopAbs_REVERSED) found++;
      } else {
        if (shape.Orientation() != TopAbs_REVERSED) found++;
      }
      if (found == 0) continue;
      edgo = wedge;
      break;
    }
    if (found != 0) break;
  }

  if (!BRep_Tool::SameRange(edgo)) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge & PCurve not SameRange (EG_getEdgeUVs)!\n");
    return EGADS_GEOMERR;
  }

  // get and evaluate the pcurve
  BRepAdaptor_Curve2d Curve2d(edgo, pface->face);
  for (i = 0; i < nts; i++) {
    Curve2d.D0(t[i], P2d);
    uv[2*i  ] = P2d.X();
    uv[2*i+1] = P2d.Y();
  }

  return EGADS_SUCCESS;
}


int
EG_getEdgeUVeval(const egObject *face, const egObject *topo, int sense,
                 double t, double *result)
{
  int             outLevel, found;
  egObject        *pcurv = NULL;
  gp_Pnt2d        P2d;
  gp_Vec2d        V12d, V22d;
  TopoDS_Edge     edge;
  TopExp_Explorer ExpW;

  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(face);

  if (topo == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Reference (EG_getEdgeUVeval)!\n");
    return EGADS_NULLOBJ;
  }
  if (topo->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: topo not an EGO (EG_getEdgeUVeval)!\n");
    return EGADS_NOTOBJ;
  }
  if ((topo->oclass != EDGE) && (topo->oclass != NODE)) {
    if (outLevel > 0)
      printf(" EGADS Error: Not an Edge (EG_getEdgeUVeval)!\n");
    return EGADS_NOTTOPO;
  }
  if ((sense < -2) || (sense > 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: Sense = %d (EG_getEdgeUVeval)!\n", sense);
    return EGADS_RANGERR;
  }
  if (EG_context(face) != EG_context(topo)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_getEdgeUVeval)!\n");
    return EGADS_MIXCNTX;
  }
  egadsFace *pface = (egadsFace *) face->blind;
  if (pface == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face pointer (EG_getEdgeUVeval)!\n");
    return EGADS_NODATA;
  }

  // undocumented -- topo is a node!
  if (topo->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Node pointer (EG_getEdgeUVeval)!\n");
      return EGADS_NODATA;
    }
    TopoDS_Vertex vert = TopoDS::Vertex(pnode->node);

    double tt = 0.0;
    found     = 0;
    for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Wire wire = TopoDS::Wire(ExpW.Current());
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
        edge = TopoDS::Edge(ExpWE.Current());
        if (BRep_Tool::Degenerated(edge)) continue;
        TopoDS_Vertex V1, V2;
        if (edge.Orientation() == TopAbs_REVERSED) {
          TopExp::Vertices(edge, V2, V1, Standard_True);
        } else {
          TopExp::Vertices(edge, V1, V2, Standard_True);
        }
        if ((vert.IsSame(V1)) || (vert.IsSame(V2))) {
          double t1, t2;
          Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
          if (vert.IsSame(V1)) {
            tt = t1;
          } else {
            tt = t2;
          }
          found = 1;
          break;
        }
      }
      if (found != 0) break;
    }
    if (found == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Node not in Face (EG_getEdgeUVeval)!\n");
      return EGADS_NOTFOUND;
    }

    BRepAdaptor_Curve2d Curve2d(edge, pface->face);
    Curve2d.D2(tt, P2d, V12d, V22d);
    result[0] = P2d.X();
    result[1] = P2d.Y();
    result[2] = V12d.X();
    result[3] = V12d.Y();
    result[4] = V22d.X();
    result[5] = V22d.Y();

    return EGADS_SUCCESS;
  }

  // topo is an edge
  egadsEdge *pedge = (egadsEdge *) topo->blind;
  if (pedge == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Edge pointer (EG_getEdgeUVeval)!\n");
    return EGADS_NODATA;
  }
  edge  = pedge->edge;
  found = 0;
  for (int iloop = 0; iloop < pface->nloops; iloop++) {
    if (pface->loops[iloop] == NULL) continue;
    egadsLoop *ploop = (egadsLoop *) pface->loops[iloop]->blind;
    if (ploop == NULL) continue;
    for (int iedge = 0; iedge < ploop->nedges; iedge++)
      if (ploop->edges[iedge] == topo) {
        if (ploop->surface != NULL) pcurv = ploop->edges[iedge+ploop->nedges];
        found++;
      }
  }
  if (found == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge/Sense not in Face (EG_getEdgeUVeval)!\n");
    return EGADS_NOTFOUND;
  }

  // we are only in the Face once
  if (found == 1)
    if (pcurv != NULL) {
      return  EG_evaluate(pcurv, &t, result);
    } else {
      BRepAdaptor_Curve2d Curve2d(edge, pface->face);
      Curve2d.D2(t, P2d, V12d, V22d);
      result[0] = P2d.X();
      result[1] = P2d.Y();
      result[2] = V12d.X();
      result[3] = V12d.Y();
      result[4] = V22d.X();
      result[5] = V22d.Y();
      return EGADS_SUCCESS;
    }

  // the Edge in the Face more than once -- find it!
  found = 0;
  for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
    TopoDS_Shape shapw = ExpW.Current();
    TopoDS_Wire  wire  = TopoDS::Wire(shapw);
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shape = ExpWE.Current();
      TopoDS_Edge  wedge = TopoDS::Edge(shape);
      if (!wedge.IsSame(edge)) continue;
      if (sense < 0) {
        if (shape.Orientation() == TopAbs_REVERSED) found++;
      } else {
        if (shape.Orientation() != TopAbs_REVERSED) found++;
      }
      if (found == 0) continue;
      edge = wedge;
      break;
    }
    if (found != 0) break;
  }

  if (!BRep_Tool::SameRange(edge)) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge & PCurve not SameRange (EG_getEdgeUVeval)!\n");
    return EGADS_GEOMERR;
  }

  // get and evaluate the pcurve
  BRepAdaptor_Curve2d Curve2d(edge, pface->face);
  Curve2d.D2(t, P2d, V12d, V22d);
  result[0] = P2d.X();
  result[1] = P2d.Y();
  result[2] = V12d.X();
  result[3] = V12d.Y();
  result[4] = V22d.X();
  result[5] = V22d.Y();

  return EGADS_SUCCESS;
}


int
EG_invEdgeUV(const egObject *face, const egObject *topo, int sense, double *uv,
             double *t, double *uvs)
{
  int             i, outLevel, found;
  double          a, b, tx, tmin, tmax, pw[2];
  gp_Pnt2d        pnt;
  gp_Vec2d        t1, t2;
  TopoDS_Edge     edge;
  TopExp_Explorer ExpW;
  static double   ratios[5] = {0.02, 0.25, 0.5, 0.75, 0.98};

  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(face);

  if (topo == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Reference (EG_invEdgeUV)!\n");
    return EGADS_NULLOBJ;
  }
  if (topo->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: topo not an EGO (EG_invEdgeUV)!\n");
    return EGADS_NOTOBJ;
  }
  if (topo->oclass != EDGE) {
    if (outLevel > 0)
      printf(" EGADS Error: Not an Edge (EG_invEdgeUV)!\n");
    return EGADS_NOTTOPO;
  }
  if ((sense < -2) || (sense > 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: Sense = %d (EG_invEdgeUV)!\n", sense);
    return EGADS_RANGERR;
  }
  if (EG_context(face) != EG_context(topo)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_invEdgeUV)!\n");
    return EGADS_MIXCNTX;
  }
  egadsFace *pface = (egadsFace *) face->blind;
  if (pface == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face pointer (EG_invEdgeUV)!\n");
    return EGADS_NODATA;
  }

  egadsEdge *pedge = (egadsEdge *) topo->blind;
  if (pedge == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Edge pointer (EG_invEdgeUV)!\n");
    return EGADS_NODATA;
  }
  edge  = pedge->edge;
  found = 0;
  for (int iloop = 0; iloop < pface->nloops; iloop++) {
    if (pface->loops[iloop] == NULL) continue;
    egadsLoop *ploop = (egadsLoop *) pface->loops[iloop]->blind;
    if (ploop == NULL) continue;
    for (int iedge = 0; iedge < ploop->nedges; iedge++)
      if (ploop->edges[iedge] == topo) found++;
  }
  if (found == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge/Sense not in Face (EG_invEdgeUV)!\n");
    return EGADS_NOTFOUND;
  }

  // is edge in the face more than once?
  if (found > 1) {
    found = 0;
    for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Shape shapw = ExpW.Current();
      TopoDS_Wire  wire  = TopoDS::Wire(shapw);
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
        TopoDS_Shape shape = ExpWE.Current();
        TopoDS_Edge  wedge = TopoDS::Edge(shape);
        if (!wedge.IsSame(edge)) continue;
        if (sense < 0) {
          if (shape.Orientation() == TopAbs_REVERSED) found++;
        } else {
          if (shape.Orientation() != TopAbs_REVERSED) found++;
        }
        if (found == 0) continue;
        edge = wedge;
        break;
      }
      if (found != 0) break;
    }
  }

  BRep_Tool::Range(edge, tmin, tmax);
  BRepAdaptor_Curve2d Curve2d(edge, pface->face);
  b = 0.0;
  for (i = 0; i < 5; i++) {
    tx = (1.0-ratios[i])*tmin + ratios[i]*tmax;
    Curve2d.D0(tx, pnt);
    a = (pnt.X()-uv[0])*(pnt.X()-uv[0]) +
        (pnt.Y()-uv[1])*(pnt.Y()-uv[1]);
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

  for (i = 0; i < 20; i++) {
    if ((*t < tmin) || (*t > tmax)) break;
    Curve2d.D2(*t, pnt, t1, t2);
    pw[0] = pnt.X() - uv[0];
    pw[1] = pnt.Y() - uv[1];
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

  Curve2d.D0(*t, pnt);
  uvs[0] = pnt.X();
  uvs[1] = pnt.Y();

  return EGADS_SUCCESS;
}


int
EG_getPCurve(const egObject *face, const egObject *topo, int sense, int *mtype,
             int **ivec, double **rvec)
{
  int             i, j, len, outLevel, found, *ints = NULL;
  double          *data = NULL;
  TopoDS_Edge     edge;
  TopExp_Explorer ExpW;

  *mtype = 0;
  *ivec  = NULL;
  *rvec  = NULL;
  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(face);

  if (topo == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Reference (EG_getPCurve)!\n");
    return EGADS_NULLOBJ;
  }
  if (topo->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: topo not an EGO (EG_getPCurve)!\n");
    return EGADS_NOTOBJ;
  }
  if (topo->oclass != EDGE) {
    if (outLevel > 0)
      printf(" EGADS Error: Not an Edge (EG_getPCurve)!\n");
    return EGADS_NOTTOPO;
  }
  if ((sense < -2) || (sense > 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: Sense = %d (EG_getPCurve)!\n", sense);
    return EGADS_RANGERR;
  }
  if (EG_context(face) != EG_context(topo)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_getPCurve)!\n");
    return EGADS_MIXCNTX;
  }
  egadsFace *pface = (egadsFace *) face->blind;
  if (pface == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face pointer (EG_getPCurve)!\n");
    return EGADS_NODATA;
  }

  egadsEdge *pedge = (egadsEdge *) topo->blind;
  if (pedge == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Edge pointer (EG_getPCurve)!\n");
    return EGADS_NODATA;
  }
  edge  = pedge->edge;
  found = 0;
  for (int iloop = 0; iloop < pface->nloops; iloop++) {
    if (pface->loops[iloop] == NULL) continue;
    egadsLoop *ploop = (egadsLoop *) pface->loops[iloop]->blind;
    if (ploop == NULL) continue;
    for (int iedge = 0; iedge < ploop->nedges; iedge++)
      if (ploop->edges[iedge] == topo) found++;
  }
  if (found == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge/Sense not in Face (EG_getPCurve)!\n");
    return EGADS_NOTFOUND;
  }

  // is edge in the face more than once?
  if (found > 1) {
    found = 0;
    for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Shape shapw = ExpW.Current();
      TopoDS_Wire  wire  = TopoDS::Wire(shapw);
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
        TopoDS_Shape shape = ExpWE.Current();
        TopoDS_Edge  wedge = TopoDS::Edge(shape);
        if (!wedge.IsSame(edge)) continue;
        if (sense < 0) {
          if (shape.Orientation() == TopAbs_REVERSED) found++;
        } else {
          if (shape.Orientation() != TopAbs_REVERSED) found++;
        }
        if (found == 0) continue;
        edge = wedge;
        break;
      }
      if (found != 0) break;
    }
  }

  BRepAdaptor_Curve2d Curve2d(edge, pface->face);
  GeomAbs_CurveType ctype = Curve2d.GetType();
  if (ctype == GeomAbs_OtherCurve) {
    if (outLevel > 0)
      printf(" EGADS Error: PCurve is type OtherCurve (EG_getPCurve)!\n");
    return EGADS_GEOMERR;
  }
  Handle(Geom2d_Curve) hCurve = Curve2d.Curve();

  // return the data associated with the PCurve

  switch (ctype) {

    case GeomAbs_Line:
      data = (double *) EG_alloc(4*sizeof(double));
      if (data == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc on PLine (EG_getPCurve)!\n");
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
        *mtype  = LINE;
      }
      break;

    case GeomAbs_Circle:
      data = (double *) EG_alloc(7*sizeof(double));
      if (data == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc on PCircle (EG_getPCurve)!\n");
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
        *mtype  = CIRCLE;
      }
      break;

    case GeomAbs_Ellipse:
      data = (double *) EG_alloc(8*sizeof(double));
      if (data == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc on PEllipse (EG_getPCurve)!\n");
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
        *mtype  = ELLIPSE;
      }
      break;

    case GeomAbs_Parabola:
      data = (double *) EG_alloc(7*sizeof(double));
      if (data == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc on PParabola (EG_getPCurve)!\n");
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
        *mtype  = PARABOLA;
      }
      break;

    case GeomAbs_Hyperbola:
      data = (double *) EG_alloc(8*sizeof(double));
      if (data == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc on PHyperbola (EG_getPCurve)!\n");
        return EGADS_MALLOC;
      } else {
        Handle(Geom2d_Hyperbola)
          hHypr = Handle(Geom2d_Hyperbola)::DownCast(hCurve);
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
        *mtype  = HYPERBOLA;
      }
      break;

    case GeomAbs_BezierCurve:
      ints = (int *) EG_alloc(3*sizeof(int));
      if (ints == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc on PBezier Curve (EG_getPCurve)!\n");
        return EGADS_MALLOC;
      } else {
        Handle(Geom2d_BezierCurve)
          hBezier = Handle(Geom2d_BezierCurve)::DownCast(hCurve);
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
            printf(" EGADS Error: Malloc on PBezier Data (EG_getPCurve)!\n");
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
        *ivec  = ints;
        *rvec  = data;
        *mtype = BEZIER;
      }
      break;

    case GeomAbs_BSplineCurve:
      ints = (int *) EG_alloc(4*sizeof(int));
      if (ints == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc on PBSpline Curve (EG_getPCurve)!\n");
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
            printf(" EGADS Error: Malloc PBSpline Data (EG_getPCurve)!\n");
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
          for (i = 1; i <= ints[2]; i++,len++) data[len] = hBSpline->Weight(i);
        *ivec  = ints;
        *rvec  = data;
        *mtype = BSPLINE;
      }
      break;

    default:
      if (outLevel > 0)
        printf(" EGADS Error: PCurve is unknown type (EG_getPCurve)!\n");
      return EGADS_GEOMERR;
  }

  return EGADS_SUCCESS;
}


int
EG_getBody(const egObject *obj, egObject **body)
{
  int i, index;

  *body = NULL;
  if (obj == NULL)                  return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (obj->blind == NULL)           return EGADS_NODATA;
  if ((obj->oclass < NODE) ||
      (obj->oclass > SHELL))        return EGADS_NOTTOPO;
  egObject *topObj = obj->topObj;
  if (topObj == NULL)               return EGADS_NULLOBJ;
  if (topObj->magicnumber != MAGIC) return EGADS_NOTOBJ;

  if (topObj->oclass == BODY) {
    *body = topObj;
  } else if (topObj->oclass == MODEL) {
    egadsModel *pmodel = (egadsModel *) topObj->blind;
    if (pmodel != NULL) {
      for (i = 0; i < pmodel->nbody; i++) {
        egObject *bod = pmodel->bodies[i];
        if (bod == NULL) continue;
        egadsBody *pbody = (egadsBody *) bod->blind;
        if (obj->oclass == NODE) {
          egadsNode *pnode = (egadsNode *) obj->blind;
          index = pbody->nodes.map.FindIndex(pnode->node);
        } else if (obj->oclass == EDGE) {
          egadsEdge *pedge = (egadsEdge *) obj->blind;
          index = pbody->edges.map.FindIndex(pedge->edge);
        } else if (obj->oclass == LOOP) {
          egadsLoop *ploop = (egadsLoop *) obj->blind;
          index = pbody->loops.map.FindIndex(ploop->loop);
        } else if (obj->oclass == FACE) {
          egadsFace *pface = (egadsFace *) obj->blind;
          index = pbody->faces.map.FindIndex(pface->face);
        } else {
          egadsShell *pshell = (egadsShell *) obj->blind;
          index = pbody->shells.map.FindIndex(pshell->shell);
        }
        if (index != 0) {
          *body = bod;
          break;
        }
      }
    }
  }

  return EGADS_SUCCESS;
}


int
EG_inTopology(const egObject *topo, const double *xyz)
{
  int           stat, outLevel;
  Standard_Real tol, t, u, v, range[2];

  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(topo);

  gp_Pnt pnt(xyz[0], xyz[1], xyz[2]);
  stat = EG_tolerance(topo, &tol);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_tolerance = %d (EG_inTopology)!\n", stat);
    return stat;
  }

  if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    egObject  *curv  = pedge->curve;
    if (curv == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No curve Object for Edge (EG_inTopology)!\n");
      return EGADS_NULLOBJ;
    }
    egadsCurve *pcurve = (egadsCurve *) curv->blind;
    if (pcurve == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No curve Data for Edge (EG_inTopology)!\n");
      return EGADS_NODATA;
    }
    Handle(Geom_Curve) hCurve = pcurve->handle;
    GeomAPI_ProjectPointOnCurve projPnt(pnt, hCurve);
    if (projPnt.NbPoints() == 0) {
      if (outLevel > 0)
        printf(" EGADS Warning: No projection on Curve (EG_inTopology)!\n");
      return EGADS_NOTFOUND;
    }
    if (projPnt.LowerDistance() > tol) return EGADS_OUTSIDE;
    t = projPnt.LowerDistanceParameter();
    BRep_Tool::Range(pedge->edge, range[0], range[1]);
    if ((t < range[0]) || (t > range[1])) return EGADS_OUTSIDE;
    return EGADS_SUCCESS;

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    egObject  *surf  = pface->surface;
    if (surf == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No Surf Object for Face (EG_inTopology)!\n");
      return EGADS_NULLOBJ;
    }
    egadsSurface *psurf = (egadsSurface *) surf->blind;
    if (psurf == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No Surf Data for Face (EG_inTopology)!\n");
      return EGADS_NODATA;
    }
    GeomAPI_ProjectPointOnSurf projPnt(pnt, psurf->handle);
    if (!projPnt.IsDone()) {
      printf(" EGADS Warning: GeomAPI_ProjectPointOnSurf (EG_inTopology)!\n");
      return EGADS_GEOMERR;
    }
    if (projPnt.LowerDistance() > tol) return EGADS_OUTSIDE;
    projPnt.LowerDistanceParameters(u, v);
    gp_Pnt2d pnt2d(u, v);
    BRepClass_FaceClassifier aClassifier(pface->face, pnt2d, tol);
    if (aClassifier.State() == TopAbs_OUT) return EGADS_OUTSIDE;
    return EGADS_SUCCESS;

  } else if ((topo->oclass == SHELL) && (topo->mtype == CLOSED)) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    TopoDS_Solid solid;
    BRep_Builder builder3D;
    builder3D.MakeSolid(solid);
    builder3D.Add(solid, pshell->shell);
    try {
      BRepLib::OrientClosedSolid(solid);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Warning: Cannot Orient Solid (EG_inTopology)!\n");
      printf("                %s\n", e.GetMessageString());
      return EGADS_TOPOERR;
    }
    catch (...) {
      printf(" EGADS Warning: Cannot Orient Solid (EG_inTopology)!\n");
      return EGADS_TOPOERR;
    }
    BRepClass3d_SolidClassifier sClassifier(solid, pnt, tol);
    if (sClassifier.State() == TopAbs_OUT) return EGADS_OUTSIDE;
    return EGADS_SUCCESS;

  } else if ((topo->oclass == BODY) && (topo->mtype == SOLIDBODY)) {

    egadsBody   *pbody = (egadsBody *) topo->blind;
    TopoDS_Solid solid = TopoDS::Solid(pbody->shape);
    BRepClass3d_SolidClassifier sClassifier(solid, pnt, tol);
    if (sClassifier.State() == TopAbs_OUT) return EGADS_OUTSIDE;
    return EGADS_SUCCESS;

  }

  return EGADS_NOTTOPO;
}


int
EG_inFace(const egObject *face, const double *uv)
{
  return EG_inFaceX(face, uv, NULL, NULL);
}


int
EG_inFaceOCC(const egObject *face, double tol, const double *uv)
{
  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;

  egadsFace  *pface = (egadsFace *) face->blind;
  gp_Pnt2d pnt2d(uv[0], uv[1]);
  BRepClass_FaceClassifier aClassifier(pface->face, pnt2d, tol);
  if (aClassifier.State() == TopAbs_OUT) return EGADS_OUTSIDE;

  return EGADS_SUCCESS;
}


int
EG_inFaceAlt(const egObject *face, const double *uv)
{
  int           i, mtype, sense, sen, outLevel;
  double        tol, uvtol, dist, a, b, d, dd, tx, t, ts;
  double        tmin, tmax, umin, umax, vmin, vmax, pw[2], uvs[2];
  gp_Pnt2d      pnt;
  gp_Vec2d      t1, t2;
  gp_Pnt        P0, E0;
  static double ratios[5] = {0.02, 0.25, 0.5, 0.75, 0.98};

  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(face);

  egadsFace *pface = (egadsFace *) face->blind;
  BRepTools::UVBounds(pface->face, umin, umax, vmin, vmax);
  uvtol = umax - umin;
  if (vmax - vmin > uvtol) uvtol = vmax - vmin;
  uvtol *= 1.e-10;
  if (uv[0]+uvtol < umin) return EGADS_OUTSIDE;
  if (uv[0]-uvtol > umax) return EGADS_OUTSIDE;
  if (uv[1]+uvtol < vmin) return EGADS_OUTSIDE;
  if (uv[1]-uvtol > vmax) return EGADS_OUTSIDE;

  egObject *surf = pface->surface;
  if (surf == NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: No Surf Object for Face (EG_inFaceAlt)!\n");
    return EGADS_NULLOBJ;
  }
  egadsSurface *psurf = (egadsSurface *) surf->blind;
  if (psurf == NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: No Surf Data for Face (EG_inFaceAlt)!\n");
    return EGADS_NODATA;
  }
  Handle(Geom_Surface) hSurface = psurf->handle;
  hSurface->D0(uv[0], uv[1], P0);
  mtype = SFORWARD;
  if (pface->face.Orientation() == TopAbs_REVERSED) mtype = SREVERSE;

  dist  = 1.e300;
  ts    = 0.0;
  sense = 1;
  TopoDS_Edge edge;
  TopExp_Explorer ExpW;
  for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
    TopoDS_Shape shapw = ExpW.Current();
    TopoDS_Wire  wire  = TopoDS::Wire(shapw);
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shape = ExpWE.Current();
      TopoDS_Edge  wedge = TopoDS::Edge(shape);
//    if (BRep_Tool::Degenerated(wedge)) continue;
      sen = SFORWARD;
      if (shape.Orientation() == TopAbs_REVERSED) sen = SREVERSE;
      BRep_Tool::Range(wedge, tmin, tmax);
      BRepAdaptor_Curve2d Curve2d(wedge, pface->face);
      b = 0.0;
      t = tmin;
      for (i = 0; i < 5; i++) {
        tx = (1.0-ratios[i])*tmin + ratios[i]*tmax;
        Curve2d.D0(tx, pnt);
        a = (pnt.X()-uv[0])*(pnt.X()-uv[0]) +
            (pnt.Y()-uv[1])*(pnt.Y()-uv[1]);
        if (i == 0) {
          t = tx;
          b = a;
        } else {
          if (a < b) {
            t = tx;
            b = a;
          }
        }
      }

      for (i = 0; i < 20; i++) {
        if ((t < tmin) || (t > tmax)) break;
        Curve2d.D2(t, pnt, t1, t2);
        pw[0] = pnt.X() - uv[0];
        pw[1] = pnt.Y() - uv[1];
        b     = -( pw[0]*t1.X() +  pw[1]*t1.Y());
        a     =  (t1.X()*t1.X() + t1.Y()*t1.Y()) +
                 ( pw[0]*t2.X() +  pw[1]*t2.Y());
        if (a == 0.0) break;
        b /= a;
        if (fabs(b) < 1.e-10*(tmax-tmin)) break;
        t += b;
      }
      if (t < tmin) t = tmin;
      if (t > tmax) t = tmax;

      Curve2d.D1(t, pnt, t1);
      d = (pnt.X()-uv[0])*(pnt.X()-uv[0]) + (pnt.Y()-uv[1])*(pnt.Y()-uv[1]);
      if (d < uvtol*uvtol) return EGADS_SUCCESS;
      if (d >= dist) continue;
      if (!BRep_Tool::Degenerated(wedge)) {
        tol = BRep_Tool::Tolerance(wedge);
        Handle(Geom_Curve) hCurve = BRep_Tool::Curve(wedge, tmin, tmax);
        hCurve->D0(t, E0);
        dist = sqrt((P0.X()-E0.X())*(P0.X()-E0.X()) +
                    (P0.Y()-E0.Y())*(P0.Y()-E0.Y()) +
                    (P0.Z()-E0.Z())*(P0.Z()-E0.Z()));
        if (dist <= tol) return EGADS_SUCCESS;
      }

//    if ((t == tmin) || (t == tmax)) {
        uvs[0] = uv[0] - pnt.X();
        uvs[1] = uv[1] - pnt.Y();
        dd = sqrt(uvs[0]*uvs[0] + uvs[1]*uvs[1]);
        if (dd != 0.0) {
          uvs[0] /= dd;
          uvs[1] /= dd;
        }
        pw[0] = t1.X();
        pw[1] = t1.Y();
        dd = sqrt(pw[0]*pw[0] + pw[1]*pw[1]);
        if (dd != 0.0) {
          pw[0] /= dd;
          pw[1] /= dd;
        }
        dd = uvs[0]*pw[1] - uvs[1]*pw[0];
        if (dd == 0.0) d = 1.1e300;
//    }
      if (d < dist) {
        dist  = d;
        edge  = wedge;
        ts    = t;
        sense = sen;
      }
    }
  }
  if (edge.IsNull()) {
    if (outLevel > 0)
      printf(" No Edge Found in EG_inFaceAlt!\n");
    return EGADS_NOTFOUND;
  }

  BRepAdaptor_Curve2d curve2d(edge, pface->face);
  curve2d.D1(ts, pnt, t1);
  uvs[0] = uv[0] - pnt.X();
  uvs[1] = uv[1] - pnt.Y();
  d = sqrt(uvs[0]*uvs[0] + uvs[1]*uvs[1]);
  if (d != 0.0) {
    uvs[0] /= d;
    uvs[1] /= d;
  }
  pw[0] = t1.X();
  pw[1] = t1.Y();
  d = sqrt(pw[0]*pw[0] + pw[1]*pw[1]);
  if (d != 0.0) {
    pw[0] /= d;
    pw[1] /= d;
  }
  d = uvs[0]*pw[1] - uvs[1]*pw[0];

  if (d*sense*mtype > 0.0) return EGADS_OUTSIDE;
  return EGADS_SUCCESS;
}


static int
EG_getRecurFrag(egObject *object, int oclass, int *ntopo, egObject ***topos)
{
  int      stat, i, oc, mt, nobj, *sens;
  double   limits[4];
  egObject *geom, **objs, **tmp;
  
  /* save away our objects */
  if (object->oclass == oclass) {
    tmp = *topos;
    for (i = 0; i < *ntopo; i++)
      if (object == tmp[i]) return EGADS_SUCCESS;

    if (*ntopo == 0) {
      tmp = (egObject **) EG_alloc(sizeof(egObject *));
    } else {
      tmp = (egObject **) EG_reall(*topos, (*ntopo+1)*sizeof(egObject *));
    }
    if (tmp == NULL) return EGADS_MALLOC;
    tmp[*ntopo] = object;
    *ntopo     += 1;
    *topos      = tmp;
    return EGADS_SUCCESS;
  }
  
  /* examine the children */
  stat = EG_getTopology(object, &geom, &oc, &mt, limits, &nobj, &objs, &sens);
  if (stat != EGADS_SUCCESS) return stat;

  for (i = 0; i < nobj; i++) {
    stat = EG_getRecurFrag(objs[i], oclass, ntopo, topos);
    if (stat != EGADS_SUCCESS) return stat;
  }
  
  return EGADS_SUCCESS;
}


static int
EG_getTopoFrags(const egObject *object, int oclass,
                int *ntopo, egObject ***topos)
{
  *ntopo = 0;
  *topos = NULL;
  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass < LOOP) ||
      (object->oclass > SHELL))      return EGADS_NOTTOPO;
  if  (object->blind == NULL)        return EGADS_NODATA;
  if ((oclass != NODE) &&
      (oclass != EDGE))              return EGADS_INDEXERR;
  
  return EG_getRecurFrag((egObject *) object, oclass, ntopo, topos);
}


int
EG_sewFaces(int nobj, const egObject **objs, double toler, int opt,
            egObject **result)
{
  int      i, j, k, n, no, outLevel, stat, nerr, fullAttr, *amap;
  double   tol, toltmp, tmin, tmax, trange[2], xyz[3];
  egObject *context, *omodel, *geom, **obs;

  *result = NULL;
  if (nobj <= 1)                     return EGADS_EMPTY;
  if (objs == NULL)                  return EGADS_NULLOBJ;
  if (objs[0] == NULL)               return EGADS_NULLOBJ;
  if (objs[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(objs[0]))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(objs[0]);
  context  = EG_context(objs[0]);
  fullAttr = EG_fullAttrs(objs[0]);

  tol = 0.0;
  for (i = 0; i < nobj; i++) {
    if (objs[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Object %d (EG_sewFaces)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (objs[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not an EGO (EG_sewFaces)!\n", i+1);
      return EGADS_NOTOBJ;
    }
    if (objs[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d has no data (EG_sewFaces)!\n", i+1);
      return EGADS_NODATA;
    }
    if (EG_context(objs[i]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d Context Mismatch (EG_sewFaces)!\n",
               i+1);
      return EGADS_MIXCNTX;
    }
    if (objs[i]->oclass == BODY) {
      if (objs[i]->mtype == WIREBODY) {
        if (outLevel > 0)
          printf(" EGADS Error: Object %d is a WireBody (EG_sewFaces)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
      egadsBody *pbody = (egadsBody *) objs[i]->blind;
      int nface = pbody->faces.map.Extent();
      for (j = 0; j < nface; j++) {
        stat = EG_tolerance(pbody->faces.objs[j], &toltmp);
        if (stat != EGADS_SUCCESS) return stat;
        if (toltmp > tol) tol = toltmp;
      }
    } else if (objs[i]->oclass == SHELL) {
      egadsShell *pshell = (egadsShell *) objs[i]->blind;
      for (j = 0; j < pshell->nfaces; j++) {
        stat = EG_tolerance(pshell->faces[j], &toltmp);
        if (stat != EGADS_SUCCESS) return stat;
        if (toltmp > tol) tol = toltmp;
      }
    } else if (objs[i]->oclass == FACE) {
      stat = EG_tolerance(objs[i], &toltmp);
      if (stat != EGADS_SUCCESS) return stat;
      if (toltmp > tol) tol = toltmp;
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d Does Not have Faces (EG_sewFaces)!\n",
               i+1);
      return EGADS_NOTTOPO;
    }
  }
  if (toler > tol) tol = toler;

  Standard_Boolean flag = Standard_False;
  if (opt == 1)    flag = Standard_True;
  BRepBuilderAPI_Sewing sew(tol, Standard_True, Standard_True, Standard_True,
                            flag);
  for (i = 0; i < nobj; i++) {
    TopoDS_Shape shape;
    if (objs[i]->oclass == BODY) {
      egadsBody *pbody = (egadsBody *) objs[i]->blind;
      shape = pbody->shape;
    } else if (objs[i]->oclass == SHELL) {
      egadsShell *pshell = (egadsShell *) objs[i]->blind;
      shape = pshell->shell;
    } else {
      egadsFace *pface = (egadsFace *) objs[i]->blind;
      shape = pface->face;
    }
    sew.Add(shape);
  }
  sew.Perform();
  TopoDS_Shape sewShape = sew.SewedShape();
  BRepCheck_Analyzer fCheck(sewShape);
  if (!fCheck.IsValid()) {
    Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(sewShape);
#if CASVER >= 700
    if (opt == 1) sfs->FixShellTool()->SetNonManifoldFlag(Standard_True);
#endif
    sfs->Perform();
    TopoDS_Shape fixedShape = sfs->Shape();
    if (fixedShape.IsNull()) {
      printf(" EGADS Warning: fixed Sew is Invalid!\n");
    } else {
      BRepCheck_Analyzer sfCheck(fixedShape);
      if (!sfCheck.IsValid()) {
        printf(" EGADS Warning: Fixed Sew is InValid!\n");
      } else {
        sewShape = fixedShape;
      }
    }
  }

  // map the faces & check orientation for non-manifold results
  TopTools_IndexedMapOfShape smap;
  TopExp::MapShapes(sewShape, TopAbs_FACE, smap);
  amap = NULL;
  n    = smap.Extent();
  if (n > 0) {
    amap = (int *) EG_alloc(3*n*sizeof(int));
    if (amap == NULL) {
      printf(" EGADS Warning: Allocation for %d fMap (EG_sewFaces)!\n", n);
    } else {
      int *flip = &amap[2*n];
      for (j = 0; j < 3*n; j++) amap[j] = 0;
      for (i = 0; i < nobj; i++) {
        if (objs[i]->oclass == BODY) {
          egadsBody *pbodf = (egadsBody *) objs[i]->blind;
          int nface = pbodf->faces.map.Extent();
          for (j = 0; j < nface; j++) {
            TopoDS_Face faci =
                     TopoDS::Face(sew.ModifiedSubShape(pbodf->faces.map(j+1)));
            k = smap.FindIndex(faci) - 1;
            if (k < 0) continue;
            amap[2*k]   = i+1;
            amap[2*k+1] = j+1;
            if (opt == 1) {
/*            printf(" OR %d = %d   %d %d\n", k+1, objs[i]->mtype,
                     faci.Orientation(), smap(k+1).Orientation());  */
              if (faci.Orientation() != smap(k+1).Orientation()) flip[k] = 1;
            }
          }
        } else if (objs[i]->oclass == SHELL) {
          egadsShell *pshelf = (egadsShell *) objs[i]->blind;
          for (j = 0; j < pshelf->nfaces; j++) {
            egadsFace *pfacf = (egadsFace *) pshelf->faces[j]->blind;
            if (pfacf == NULL) continue;
            TopoDS_Face faci = TopoDS::Face(sew.ModifiedSubShape(pfacf->face));
            k = smap.FindIndex(faci) - 1;
            if (k < 0) continue;
            amap[2*k]   = i+1;
            amap[2*k+1] = j+1;
            if (opt == 1) {
/*            printf(" OR %d = %d   %d %d\n", k+1, objs[i]->mtype,
                     faci.Orientation(), smap(k+1).Orientation());  */
              if (faci.Orientation() != smap(k+1).Orientation()) flip[k] = 1;
            }
          }
        } else {
          egadsFace *pface = (egadsFace *) objs[i]->blind;
          TopoDS_Face faci = TopoDS::Face(sew.Modified(pface->face));
          k = smap.FindIndex(faci) - 1;
          if (k < 0) continue;
          amap[2*k]   = i+1;
          amap[2*k+1] = 1;
          if (opt == 1) {
/*          printf(" OR %d = %d   %d %d\n", k+1, objs[i]->mtype,
                   faci.Orientation(), smap(k+1).Orientation());  */
            if (faci.Orientation() != smap(k+1).Orientation()) flip[k] = 1;
          }
        }
      }

      // reorient flipped faces (non-manifold)
      if (opt == 1) {
        int cnt = 0;
        for (k = 0; k < n; k++) if (flip[k] == 1) cnt++;
        if (cnt != 0) {
          BRepTools_ReShape reshape;
          for (k = 0; k < n; k++)
            if (flip[k] == 1)
              reshape.Replace(smap(k+1), smap(k+1).Reversed());
          TopoDS_Shape newShape = reshape.Apply(sewShape, TopAbs_FACE);
          BRepCheck_Analyzer fCheck(newShape);
          if (!fCheck.IsValid()) {
            Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
            sfs->Perform();
            TopoDS_Shape fixedShape = sfs->Shape();
            if (fixedShape.IsNull()) {
              printf(" EGADS Warning: reshaped Sew is Invalid!\n");
            } else {
              BRepCheck_Analyzer sfCheck(fixedShape);
              if (!sfCheck.IsValid()) {
                printf(" EGADS Warning: Reshaped Sew is InValid!\n");
              } else {
                newShape = fixedShape;
              }
            }
          }
          sewShape = newShape;
        }
      }

    }
  }

  // check for promoting sheets to solids
  TopoDS_Compound newShapes;
  BRep_Builder    builder3D;
  builder3D.MakeCompound(newShapes);

  TopExp_Explorer Exp;
  for (Exp.Init(sewShape, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) {
    TopoDS_Shape shape = Exp.Current();
    TopoDS_Shell shell = TopoDS::Shell(shape);
    TopoDS_Solid solid;
    builder3D.MakeSolid(solid);
    builder3D.Add(solid, shell);
    BRepCheck_Analyzer sCheck(solid);
    if (!sCheck.IsValid()) {
      solid.Nullify();
    } else {
      try {
        BRepGProp    BProps;
        GProp_GProps VProps;

/*      BRepLib::OrientClosedSolid(solid);  */
        BProps.VolumeProperties(solid, VProps);
        if (VProps.Mass() < 0.0) solid.Reverse();
      }
      catch (const Standard_Failure& e) {
        solid.Nullify();
      }
      catch (...) {
        solid.Nullify();
      }
    }
    if (!solid.IsNull()) {
      builder3D.Add(newShapes, solid);
    } else {
      builder3D.Add(newShapes, shape);
    }
  }
  for (Exp.Init(sewShape, TopAbs_FACE, TopAbs_SHELL);
       Exp.More(); Exp.Next())
    builder3D.Add(newShapes, Exp.Current());
  for (Exp.Init(sewShape, TopAbs_SOLID); Exp.More(); Exp.Next())
    builder3D.Add(newShapes, Exp.Current());

  sewShape = newShapes;

  // count our bodies
  int nFace  = 0;
  int nSheet = 0;
  int nSolid = 0;
  for (Exp.Init(sewShape, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) nFace++;
  for (Exp.Init(sewShape, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) nSheet++;
  for (Exp.Init(sewShape, TopAbs_SOLID); Exp.More(); Exp.Next()) nSolid++;

  if (outLevel > 1)
    printf(" EGADS Info: Sewn Result has %d Solids, %d Sheets and %d Faces\n",
           nSolid, nSheet, nFace);

  int nBody = nFace+nSheet+nSolid;
  if (nBody == 0) {
    if (amap != NULL) EG_free(amap);
    sewShape.Nullify();
    if (outLevel > 0)
      printf(" EGADS Warning: Nothing found in Result (EG_sewFaces)!\n");
    return EGADS_NODATA;
  }

  egadsModel *mshape  = new egadsModel;
  mshape->shape       = sewShape;
  mshape->nbody       = nBody;
  mshape->bodies      = new egObject*[nBody];
  mshape->bbox.filled = 0;
  for (i = 0; i < nBody; i++) {
    stat = EG_makeObject(context, &mshape->bodies[i]);
    if (stat != EGADS_SUCCESS) {
      for (int j = 0; j < i; j++) {
        egObject  *obj   = mshape->bodies[j];
        egadsBody *pbody = (egadsBody *) obj->blind;
        delete pbody;
        EG_deleteObject(mshape->bodies[j]);
      }
      delete [] mshape->bodies;
      delete mshape;
      if (amap != NULL) EG_free(amap);
      return stat;
    }
    egObject  *pobj    = mshape->bodies[i];
    egadsBody *pbody   = new egadsBody;
    pbody->nodes.objs  = NULL;
    pbody->edges.objs  = NULL;
    pbody->loops.objs  = NULL;
    pbody->faces.objs  = NULL;
    pbody->shells.objs = NULL;
    pbody->senses      = NULL;
    pbody->bbox.filled = 0;
    pbody->massFill    = 0;
    pobj->blind        = pbody;
  }

  i = 0;
  for (Exp.Init(mshape->shape, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
  }
  for (Exp.Init(mshape->shape, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
  }
  for (Exp.Init(mshape->shape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
/*
    TopoDS_Shape shape = Exp.Current();
    TopoDS_Solid solid = TopoDS::Solid(shape);
    BRepCheck_Analyzer sCheck(solid);
    if (!sCheck.IsValid()) {
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(solid);
      sfs->Perform();
      TopoDS_Shape fixedSolid = sfs->Shape();
      if (fixedSolid.IsNull()) {
        printf(" EGADS Warning: Cannot Fix Solid (EG_sewFaces)!\n");
      } else {
        if (outLevel > 0)
          printf(" EGADS Warning: Fixing Solid (EG_sewFaces)!\n");
        solid = TopoDS::Solid(fixedSolid);
      }
    }
    pbody->shape = solid;
 */
  }

  stat = EG_makeObject(context, &omodel);
  if (stat != EGADS_SUCCESS) {
    sewShape.Nullify();
    for (i = 0; i < nBody; i++) {
      egObject  *obj   = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) obj->blind;
      delete pbody;
      EG_deleteObject(mshape->bodies[i]);
    }
    delete [] mshape->bodies;
    delete mshape;
    if (amap != NULL) EG_free(amap);
    return stat;
  }
  omodel->oclass = MODEL;
  omodel->blind  = mshape;
  EG_referenceObject(omodel, context);

  for (i = 0; i < nBody; i++) {
    egObject  *pobj  = mshape->bodies[i];
    egadsBody *pbody = (egadsBody *) pobj->blind;
    pobj->topObj     = omodel;
    stat = EG_traverseBody(context, i, pobj, omodel, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      mshape->nbody = i;
      EG_destroyTopology(omodel);
      if (amap != NULL) EG_free(amap);
      return stat;
    }
  }

  // set the attributes
  if (amap != NULL) {
    TopTools_IndexedMapOfShape mmap;
    TopExp::MapShapes(sewShape, TopAbs_FACE, mmap);
    for (i = 0; i < nBody; i++) {
      egObject  *pobj  = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) pobj->blind;
      for (j = 0; j < pbody->faces.map.Extent(); j++) {
        TopoDS_Face facb = TopoDS::Face(pbody->faces.map(j+1));
        int index = mmap.FindIndex(facb) - 1;
        if (index < 0) continue;
//      printf(" body %d: face %d -> %d %d\n", i+1, j+1, index, amap[2*index]);
        k     = amap[2*index+1] - 1;
        index = amap[2*index  ] - 1;
        if (index < 0) continue;
        if (objs[index]->oclass == BODY) {
          egadsBody *pbodf = (egadsBody *) objs[index]->blind;
          if ((k < 0) || (k >= pbodf->faces.map.Extent())) continue;
          EG_attributeDup(pbodf->faces.objs[k], pbody->faces.objs[j]);
        } else if (objs[index]->oclass == SHELL) {
          egadsShell *pshelf = (egadsShell *) objs[index]->blind;
          if ((k < 0) || (k >= pshelf->nfaces)) continue;
          EG_attributeDup(pshelf->faces[k], pbody->faces.objs[j]);
        } else {
          EG_attributeDup(objs[index], pbody->faces.objs[j]);
        }
      }
    }
    EG_free(amap);
  }
  
  /* deal with extended attribution */
  if (fullAttr == 1) {

    for (i = 0; i < nobj; i++) {
      if (objs[i]->oclass == BODY) {
        stat = EG_getBodyTopos(objs[i], NULL, EDGE, &no, &obs);
      } else {
        stat = EG_getTopoFrags(objs[i], EDGE, &no, &obs);
      }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 1)
          printf(" EGADS Warning: Getting Edge objects returns = %d\n", stat);
        continue;
      }
      for (int ib = 0; ib < nBody; ib++) {
        egObject  *pobj  = mshape->bodies[ib];
        egadsBody *pbody = (egadsBody *) pobj->blind;
        for (int k = 0; k < no; k++) {
          if (obs[k]->mtype == DEGENERATE) continue;
          egadsEdge *pedge = (egadsEdge *) obs[k]->blind;
          trange[0] = pedge->trange[0];
          trange[1] = pedge->trange[1];
          geom      = pedge->curve;
          for (j = 0; j < pbody->edges.map.Extent(); j++) {
            if (pbody->edges.objs[j]->mtype == DEGENERATE) continue;
            if (EG_isSame(obs[k], pbody->edges.objs[j]) == EGADS_SUCCESS) {
              egadsEdge *pedgd = (egadsEdge *) pbody->edges.objs[j]->blind;
              egadsNode *pnod0 = (egadsNode *) pedgd->nodes[0]->blind;
              egadsNode *pnod1 = (egadsNode *) pedgd->nodes[1]->blind;
              stat = EG_invEvaluate(geom, pnod0->xyz, &tmin, xyz);
              if (stat != EGADS_SUCCESS) continue;
              stat = EG_invEvaluate(geom, pnod1->xyz, &tmax, xyz);
              if (stat != EGADS_SUCCESS) continue;
              if ((tmin+UVTOL < trange[0]) || (tmax-UVTOL > trange[1])) continue;
              EG_attributeDup(obs[k], pbody->edges.objs[j]);
            }
          }
        }
      }
      EG_free(obs);
    }
    
    for (i = 0; i < nobj; i++) {
      if (objs[i]->oclass == BODY) {
        stat = EG_getBodyTopos(objs[i], NULL, NODE, &no, &obs);
      } else {
        stat = EG_getTopoFrags(objs[i], NODE, &no, &obs);
      }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 1)
          printf(" EGADS Warning: Getting Node objects returns = %d\n", stat);
        continue;
      }
      for (int ib = 0; ib < nBody; ib++) {
        egObject  *pobj  = mshape->bodies[ib];
        egadsBody *pbody = (egadsBody *) pobj->blind;
        for (int k = 0; k < no; k++)
          for (j = 0; j < pbody->nodes.map.Extent(); j++)
            if (EG_isSame(obs[k], pbody->nodes.objs[j]) == EGADS_SUCCESS)
              EG_attributeDup(obs[k], pbody->nodes.objs[j]);
      }
      EG_free(obs);
    }
    
  }

  *result = omodel;
  return EGADS_SUCCESS;
}


int
EG_replaceFaces(const egObject *body, int nobj, egObject **objs,
                      egObject **result)
{
  int      i, j, k, outLevel, stat, mtype, nerr, fullAttr, no;
  double   tmin, tmax, trange[2], xyz[3];
  egObject *context, *obj, *src, *geom, **obs;

  *result = NULL;
  if  (body == NULL)               return EGADS_NULLOBJ;
  if  (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (body->oclass != BODY)       return EGADS_NOTBODY;
  if  (body->blind == NULL)        return EGADS_NODATA;
  if ((body->mtype != SHEETBODY) &&
      (body->mtype != SOLIDBODY))  return EGADS_TOPOERR;
  if  (nobj < 1)                   return EGADS_EMPTY;
  if  (objs == NULL)               return EGADS_NULLOBJ;
  if  (EG_sameThread(body))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(body);
  context  = EG_context(body);
  fullAttr = EG_fullAttrs(body);

  // check the input objects
  egadsBody *pbody = (egadsBody *) body->blind;
  if (body->mtype == SOLIDBODY)
    if (pbody->shells.map.Extent() != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: SolidBody with > 1 Shell (EG_replaceFaces)!\n");
      return EGADS_TOPOERR;
    }
  for (i = 0; i < nobj; i++) {
    if (objs[2*i  ] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Src Object %d (EG_replaceFaces)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (objs[2*i  ]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Src Object %d is not an EGO (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTOBJ;
    }
    if (objs[2*i  ]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Src Object %d has no data (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NODATA;
    }
    if (objs[2*i  ]->oclass != FACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Src Object %d is Not a Face (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTTOPO;
    }
    egadsFace *pface = (egadsFace *) objs[2*i  ]->blind;
    if (pbody->faces.map.FindIndex(pface->face) == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is NOT in Body (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTBODY;
    }

    if (objs[2*i+1] == NULL) continue;
    if (objs[2*i+1]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Repl Object %d is not an EGO (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTOBJ;
    }
    if (objs[2*i+1]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Repl Object %d has no data (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NODATA;
    }
    if (EG_context(objs[2*i+1]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Repl Obj %d Context Mismatch (EG_replaceFaces)!\n",
               i+1);
      return EGADS_MIXCNTX;
    }
    if (objs[2*i+1]->oclass != FACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Repl Object %d is Not a Face (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTTOPO;
    }
  }

  // make the changes
  BRepTools_ReShape reshape;
  for (i = 0; i < nobj; i++) {
    egadsFace *pface = (egadsFace *) objs[2*i  ]->blind;
    if (objs[2*i+1] == NULL) {
      reshape.Remove(pface->face);
    } else {
      egadsFace *pfacn = (egadsFace *) objs[2*i+1]->blind;
      reshape.Replace(pface->face, pfacn->face);
    }
  }

  // get mappings for attributes
  int nNull = 0;
  int *amap = new int[pbody->faces.map.Extent()];
  for (i = 1; i <= pbody->faces.map.Extent(); i++) {
    TopoDS_Face face = TopoDS::Face(pbody->faces.map(i));
    TopoDS_Face facn = TopoDS::Face(reshape.Value(face));
    if (facn.IsNull()) {
      amap[i-1] = 0;
      nNull++;
    } else if (facn.IsSame(face)) {
      amap[i-1] = i;
    } else {
      for (j = 0; j < nobj; j++) {
        if (objs[2*j+1] == NULL) continue;
        egadsFace *pfacn = (egadsFace *) objs[2*j+1]->blind;
        if (facn.IsSame(pfacn->face)) {
          amap[i-1] = -(j+1);
          break;
        }
      }
      if (j == nobj) {
        if (outLevel > 0)
          printf(" EGADS Error: New Face not Found (EG_replaceFaces)!\n");
        delete [] amap;
        return EGADS_NOTFOUND;
      }
    }
  }

  // apply the changes
  TopoDS_Shape newShape;
  if (body->mtype == SOLIDBODY) {
    egadsShell *pshell = (egadsShell *) pbody->shells.objs[0]->blind;
    newShape = reshape.Apply(pshell->shell, TopAbs_FACE);
  } else {
    newShape = reshape.Apply(pbody->shape, TopAbs_FACE);
  }

  // check our result -- sheet body
  if (newShape.ShapeType() != TopAbs_SHELL) {
    if (outLevel > 0)
      printf(" EGADS Error: Result Not a Sheet (EG_replaceFaces)!\n");
    delete [] amap;
    return EGADS_CONSTERR;
  }
  mtype = SHEETBODY;

  // promote to a solid?
  if ((body->mtype == SOLIDBODY) && (nNull == 0)) {
    TopoDS_Shell shell = TopoDS::Shell(newShape);
    BRep_Builder builder3D;
    TopoDS_Solid solid;
    builder3D.MakeSolid(solid);
    builder3D.Add(solid, shell);
    try {
      BRepLib::OrientClosedSolid(solid);
    }
    catch (const Standard_Failure& e) {
/*    printf(" EGADS Warning: Cannot Orient Solid (EG_replaceFaces)!\n");
      printf("                %s\n", e.GetMessageString()); */
      solid.Nullify();
    }
    catch (...) {
//    printf(" EGADS Warning: Cannot Orient Solid (EG_replaceFaces)!\n");
      solid.Nullify();
    }
    if (!solid.IsNull()) {
      BRepCheck_Analyzer sCheck(solid);
      if (!sCheck.IsValid()) {
      printf(" EGADS Warning: Invalid Solid (EG_replaceFaces)!\n");
        solid.Nullify();
      }
    }
    if (!solid.IsNull()) {
      newShape = solid;
      mtype    = SOLIDBODY;
    }
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_replaceFaces)!\n");
    delete [] amap;
    return stat;
  }
  egadsBody *pbodn   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = mtype;
  pbodn->nodes.objs  = NULL;
  pbodn->edges.objs  = NULL;
  pbodn->loops.objs  = NULL;
  pbodn->faces.objs  = NULL;
  pbodn->shells.objs = NULL;
  pbodn->senses      = NULL;
  pbodn->shape       = newShape;
  pbodn->bbox.filled = 0;
  pbodn->massFill    = 0;
  obj->blind         = pbodn;
  stat = EG_traverseBody(context, 0, obj, obj, pbodn, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete [] amap;
    delete pbodn;
    return stat;
  }

  // assign attributes
  if (fullAttr == 1) EG_attriBodyDup(body, obj);
  if (pbody->faces.map.Extent() != pbodn->faces.map.Extent()+nNull) {
    if (outLevel > 0)
      printf(" EGADS Warning: Attribute length Mismatch (EG_replaceFaces)!\n");
  } else {
    for (j = i = 0; i < pbody->faces.map.Extent(); i++) {
      if (amap[i] == 0) continue;
      if (amap[i] >  0) {
        src   = pbody->faces.objs[amap[i]-1];
      } else {
        int k = -amap[i]-1;
        src   = objs[2*k+1];
      }
      if (fullAttr == 1) EG_attributeDel(pbodn->faces.objs[j], NULL);
      EG_attributeDup(src, pbodn->faces.objs[j]);
      j++;
    }
  }
  delete [] amap;
  
  /* deal with extended attribution */
  if (fullAttr == 1) {

    for (i = 0; i < nobj; i++) {
      if (objs[2*i+1] == NULL) continue;
      stat = EG_getTopoFrags(objs[2*i+1], EDGE, &no, &obs);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 1)
          printf(" EGADS Warning: Getting Edge objects returns = %d\n", stat);
        continue;
      }
      for (k = 0; k < no; k++) {
        if (obs[k]->mtype == DEGENERATE) continue;
        egadsEdge *pedge = (egadsEdge *) obs[k]->blind;
        trange[0] = pedge->trange[0];
        trange[1] = pedge->trange[1];
        geom      = pedge->curve;
        for (j = 0; j < pbodn->edges.map.Extent(); j++) {
          if (pbodn->edges.objs[j]->mtype == DEGENERATE) continue;
          if (EG_isSame(obs[k], pbodn->edges.objs[j]) == EGADS_SUCCESS) {
            egadsEdge *pedgd = (egadsEdge *) pbodn->edges.objs[j]->blind;
            egadsNode *pnod0 = (egadsNode *) pedgd->nodes[0]->blind;
            egadsNode *pnod1 = (egadsNode *) pedgd->nodes[1]->blind;
            stat = EG_invEvaluate(geom, pnod0->xyz, &tmin, xyz);
            if (stat != EGADS_SUCCESS) continue;
            stat = EG_invEvaluate(geom, pnod1->xyz, &tmax, xyz);
            if (stat != EGADS_SUCCESS) continue;
            if ((tmin+UVTOL < trange[0]) || (tmax-UVTOL > trange[1])) continue;
            EG_attributeDup(obs[k], pbodn->edges.objs[j]);
          }
        }
      }
      EG_free(obs);
    }
    
    for (i = 0; i < nobj; i++) {
      if (objs[2*i+1] == NULL) continue;
      stat = EG_getTopoFrags(objs[2*i+1], NODE, &no, &obs);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 1)
          printf(" EGADS Warning: Getting Node objects returns = %d\n", stat);
        continue;
      }
      for (k = 0; k < no; k++)
        for (j = 0; j < pbodn->nodes.map.Extent(); j++)
          if (EG_isSame(obs[k], pbodn->nodes.objs[j]) == EGADS_SUCCESS)
            EG_attributeDup(obs[k], pbodn->nodes.objs[j]);
      EG_free(obs);
    }
    
  }

  EG_referenceObject(obj, context);
  *result = obj;

  return EGADS_SUCCESS;
}


static int
EG_edgeCmp(double *cg0, double *cg1)
{
  if (fabs(cg0[0]-cg1[0]) > 1.e-6) {
    if (cg0[0] < cg1[0]) return 1;
    return 0;
  }
  if (fabs(cg0[1]-cg1[1]) > 1.e-6) {
    if (cg0[1] < cg1[1]) return 1;
    return 0;
  }
  if (fabs(cg0[2]-cg1[2]) > 1.e-6) {
    if (cg0[2] < cg1[2]) return 1;
    return 0;
  }

  if (cg0[3] < cg1[3]) return 1;
  return 0;
}


static int
EG_getEdgeIDs(egadsBody *pbody, const char *fAttr, edgeID **IDs)
{
  int          i, j, k, m, n, len, index, cnt, hit, stat, aType, aLen;
  char         **facAttr, line[1025];
  const int    *ints;
  const char   *str;
  const double *reals;
  edgeID       *edgeIDs;
  gp_Pnt       CofG;
  BRepGProp    BProps;
  GProp_GProps SProps;

  *IDs    = NULL;
  len     = pbody->edges.map.Extent();
  edgeIDs = (edgeID *) EG_alloc(len*sizeof(edgeID));
  if (edgeIDs == NULL) return EGADS_MALLOC;
  for (i = 0; i < len; i++) {
    edgeIDs[i].nFace    = 0;
    edgeIDs[i].seq      = 0;
    edgeIDs[i].fIndices = NULL;
    edgeIDs[i].ID       = NULL;
  }

  // find # of Faces per Edge
  for (i = 1; i <= pbody->faces.map.Extent(); i++) {
    TopTools_IndexedMapOfShape smap;
    TopExp::MapShapes(pbody->faces.map(i), TopAbs_EDGE, smap);
    n = smap.Extent();
    for (j = 1; j <= n; j++) {
      index = pbody->edges.map.FindIndex(smap(j));
      if (index == 0) {
        printf(" EGADS Internal: Face %d -- cannot find edge %d!\n", i, j);
        continue;
      }
      edgeIDs[index-1].nFace++;
    }
  }

  // allocate our space per Edge
  for (i = 0; i < len; i++) {
    if (edgeIDs[i].nFace == 0) continue;
    edgeIDs[i].fIndices = (int *) EG_alloc((edgeIDs[i].nFace+1)*sizeof(int));
    if (edgeIDs[i].fIndices == NULL) {
      for (j = 0; j < i; j++)
        if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
      EG_free(edgeIDs);
      return EGADS_MALLOC;
    }
    for (j = 0; j < edgeIDs[i].nFace+1; j++) edgeIDs[i].fIndices[j] = 0;
    edgeIDs[i].nFace = 0;
  }

  // store Face indices
  for (i = 1; i <= pbody->faces.map.Extent(); i++) {
    TopTools_IndexedMapOfShape smap;
    TopExp::MapShapes(pbody->faces.map(i), TopAbs_EDGE, smap);
    n = smap.Extent();
    for (j = 1; j <= n; j++) {
      index = pbody->edges.map.FindIndex(smap(j));
      edgeIDs[index-1].fIndices[edgeIDs[index-1].nFace] = i;
      edgeIDs[index-1].nFace++;
    }
  }

  // look for duplicates
  for (i = 0; i < len-1; i++) {
    cnt = 0;
    if (edgeIDs[i].fIndices[edgeIDs[i].nFace] != 0) continue;
    for (j = i+1; j < len; j++) {
      if (edgeIDs[i].nFace != edgeIDs[j].nFace) continue;
      for (k = 0; k < edgeIDs[i].nFace; k++)
        if (edgeIDs[i].fIndices[k] != edgeIDs[j].fIndices[k]) break;
      if (k == edgeIDs[i].nFace) {
        if (edgeIDs[i].fIndices[k] == 0) {
          BProps.LinearProperties(pbody->edges.map(i+1), SProps);
          CofG = SProps.CentreOfMass();
          edgeIDs[i].seq         = 1;
          edgeIDs[i].fIndices[k] = -i-1;
          edgeIDs[i].CG[0]       = CofG.X();
          edgeIDs[i].CG[1]       = CofG.Y();
          edgeIDs[i].CG[2]       = CofG.Z();
          edgeIDs[i].CG[3]       = SProps.Mass();
/*        printf("  %d: %lf %lf %lf %lf\n", i+1, CofG.X(), CofG.Y(), CofG.Z(),
                 SProps.Mass());  */
        }
        BProps.LinearProperties(pbody->edges.map(j+1), SProps);
        CofG = SProps.CentreOfMass();
        edgeIDs[j].seq         = cnt+2;
        edgeIDs[j].fIndices[k] = -i-1;
        edgeIDs[j].CG[0]       = CofG.X();
        edgeIDs[j].CG[1]       = CofG.Y();
        edgeIDs[j].CG[2]       = CofG.Z();
        edgeIDs[j].CG[3]       = SProps.Mass();
/*      printf("  %d: %lf %lf %lf %lf\n", j+1, CofG.X(), CofG.Y(), CofG.Z(),
               SProps.Mass());  */
        cnt++;
      }
    }
    if (cnt != 0) {
      do {
        hit = 0;
        for (m = 1; m <= cnt; m++) {
          for (j = 0; j < len; j++)
            if ((edgeIDs[j].fIndices[edgeIDs[j].nFace] == -i-1) &&
                (edgeIDs[j].seq == m)) break;
          if (j == len) {
            printf("  Not found ERROR  %d!\n", m);
            continue;
          }
          for (k = 0; k < len; k++)
            if ((edgeIDs[k].fIndices[edgeIDs[k].nFace] == -i-1) &&
                (edgeIDs[k].seq == m+1)) break;
          if (k == len) {
            printf("  Not found ERROR  %d!\n", m+1);
            continue;
          }
          if (EG_edgeCmp(edgeIDs[j].CG, edgeIDs[k].CG) == 1) continue;
          edgeIDs[k].seq = m;
          edgeIDs[j].seq = m+1;
          hit++;
        }
      } while (hit != 0);
    }
  }

  // convert FaceIDs to strings

  facAttr = (char **) EG_alloc(pbody->faces.map.Extent()*sizeof(char *));
  if (facAttr == NULL) {
    for (j = 0; j < len; j++)
      if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
    EG_free(edgeIDs);
    return EGADS_MALLOC;
  }
  for (i = 0; i < pbody->faces.map.Extent(); i++) facAttr[i] = NULL;
  for (i = 0; i < pbody->faces.map.Extent(); i++) {
    stat = EG_attributeRet(pbody->faces.objs[i], fAttr, &aType, &aLen,
                           &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Face Attr = %d for Face %d (EG_mapBody)!\n",
             stat, i+1);
      for (j = 0; j < i; j++) EG_free(facAttr[i]);
      EG_free(facAttr);
      for (j = 0; j < len; j++)
        if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
      EG_free(edgeIDs);
      return stat;
    }
    if (aType == ATTRSTRING) {
      facAttr[i] = EG_strdup(str);
    } else if (aType == ATTRREAL) {
      sprintf(line, "%20.13le", reals[0]);
      for (j = 1; j < aLen; j++) {
        hit = strlen(line);
        if (hit+21 > 1024) {
          printf(" EGADS Warning: FaceID %d truncated to %d ints (EG_mapBody)!\n",
                 i+1, j);
          break;
        }
        sprintf(&line[hit], ":%20.13le", reals[j]);
      }
      facAttr[i] = EG_strdup(line);
    } else if (aType == ATTRINT) {
      sprintf(line, "%d", ints[0]);
      for (j = 1; j < aLen; j++) {
        hit = strlen(line);
        if (hit+10 > 1024) {
          printf(" EGADS Warning: FaceID %d truncated to %d ints (EG_mapBody)!\n",
                 i+1, j);
          break;
        }
        sprintf(&line[hit], ":%d", ints[j]);
      }
      facAttr[i] = EG_strdup(line);
    }
  }

  // do we have all of the Face strings?
  for (i = 0; i < pbody->faces.map.Extent(); i++)
    if (facAttr[i] == NULL) {
      printf(" EGADS Error: FaceID %d converted to NULL (EG_mapBody)!\n",
             i+1);
      for (i = 0; i < pbody->faces.map.Extent(); i++) EG_free(facAttr[i]);
      EG_free(facAttr);
      for (j = 0; j < len; j++)
        if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
      EG_free(edgeIDs);
      return EGADS_MALLOC;
    }

  // make the Edge ID strings from the Faces

  for (i = 0; i < len; i++) {
    // order the Faces based on strcmps
    do {
      hit = 0;
      for (m = 0; m < edgeIDs[i].nFace-1; m++) {
        if (strcmp(facAttr[edgeIDs[i].fIndices[m  ]-1],
                   facAttr[edgeIDs[i].fIndices[m+1]-1]) <= 0) continue;
        k                        = edgeIDs[i].fIndices[m  ];
        edgeIDs[i].fIndices[m  ] = edgeIDs[i].fIndices[m+1];
        edgeIDs[i].fIndices[m+1] = k;
        hit++;
      }
    } while (hit != 0);
    for (cnt = j = 0; j < edgeIDs[i].nFace; j++)
      cnt += strlen(facAttr[edgeIDs[i].fIndices[j]-1]) + 1;
    edgeIDs[i].ID = (char *) EG_alloc((cnt+11)*sizeof(char));
    if (edgeIDs[i].ID == NULL) continue;
    // concatinate the ordered Face strings
    for (cnt = j = 0; j < edgeIDs[i].nFace; j++) {
      sprintf(&edgeIDs[i].ID[cnt], "%s|", facAttr[edgeIDs[i].fIndices[j]-1]);
      cnt += strlen(facAttr[edgeIDs[i].fIndices[j]-1]) + 1;
    }
    // append the sequence number
    cnt--;
    sprintf(&edgeIDs[i].ID[cnt], "\\%d", edgeIDs[i].seq);
  }
  for (i = 0; i < pbody->faces.map.Extent(); i++) EG_free(facAttr[i]);
  EG_free(facAttr);

  // do we have them all?
  for (i = 0; i < pbody->edges.map.Extent(); i++) {
/*  printf("  Edge %d: %d Faces ->", i+1, edgeIDs[i].nFace);
    for (j = 0; j < edgeIDs[i].nFace; j++)
      printf("  %d", edgeIDs[i].fIndices[j]);
    printf(" : %d\n", edgeIDs[i].seq);  */
    if (edgeIDs[i].ID == NULL) {
      printf(" EGADS Error: edgeID %d is NULL (EG_mapBody)!\n", i+1);
      for (j = 0; j < len; j++) {
        if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
        if (edgeIDs[j].ID       != NULL) EG_free(edgeIDs[j].ID);
      }
      EG_free(edgeIDs);
      return EGADS_MALLOC;
    }
    if (edgeIDs[i].fIndices != NULL) EG_free(edgeIDs[i].fIndices);
/*  printf(" Edge %d: %s\n", i+1, edgeIDs[i].ID);  */
  }

  *IDs = edgeIDs;
  return EGADS_SUCCESS;
}


static int
EG_sameAttrs(int aTypes, int lens, const int *ints, const double *reals,
             const char *strs, int aTyped, int lend, const int *intd,
             const double *reald, const char *strd)
{
  int i;

  if (aTypes != aTyped) return EGADS_OUTSIDE;

  if (aTypes == ATTRINT) {
    if (lens != lend) return EGADS_OUTSIDE;
    for (i = 0; i < lens; i++)
      if (ints[i] != intd[i]) return EGADS_OUTSIDE;
  } else if ((aTypes == ATTRREAL) || (aTypes == ATTRCSYS)) {
    if (lens != lend) return EGADS_OUTSIDE;
    for (i = 0; i < lens; i++)
      if (reals[i] != reald[i]) return EGADS_OUTSIDE;
  } else if (aTypes == ATTRSTRING) {
    if (strcmp(strs,strd) != 0) return EGADS_OUTSIDE;
  }

  return EGADS_SUCCESS;
}


static int
EG_bodyMapping(const egObject  *sBody, const egObject *dBody, const char *fAttr,
               int **pnMap, int **peMap, int **pfMap)
{
  int          i, j, in0, in1, jn0, jn1, outLevel, stat;
  int          aTypes, aTyped, lens, lend;
  int          *nMap = NULL, *eMap = NULL, *fMap = NULL;
  edgeID       *edgeIDs, *edgeIDd;
  const int    *ints,  *intd;
  const char   *strs,  *strd;
  const double *reals, *reald;

  egadsBody *pbods = (egadsBody *) sBody->blind;
  egadsBody *pbodd = (egadsBody *) dBody->blind;
  outLevel = EG_outLevel(sBody);

  if (fAttr == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face Attribute (EG_mapBody)!\n");
    return EGADS_NONAME;
  }

  // Faces first
  fMap = (int *) EG_alloc(pbods->faces.map.Extent()*sizeof(int));
  if (fMap == NULL) return EGADS_MALLOC;
  for (i = 0; i < pbods->faces.map.Extent(); i++) fMap[i] = 0;
  for (i = 0; i < pbods->faces.map.Extent(); i++) {
    stat = EG_attributeRet(pbods->faces.objs[i], fAttr, &aTypes, &lens,
                           &ints, &reals, &strs);
    if (stat == EGADS_NOTFOUND) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Attr Not in Face %d (EG_mapBody)!\n",
               i+1);
      EG_free(fMap);
      return stat;
    }
    for (j = 0; j < pbodd->faces.map.Extent(); j++) {
      stat = EG_attributeRet(pbodd->faces.objs[j], fAttr, &aTyped, &lend,
                             &intd, &reald, &strd);
      if (stat != EGADS_SUCCESS) continue;
      if (EG_sameAttrs(aTypes, lens, ints, reals, strs,
                       aTyped, lend, intd, reald, strd) == EGADS_SUCCESS) {
        if (fMap[i] == 0) {
          fMap[i] = j+1;
        } else {
          if (outLevel > 0)
            printf(" EGADS Error: Face Attr Multip in %d (EG_mapBody)!\n",
                   i+1);
          EG_free(fMap);
          return EGADS_ATTRERR;
        }
      }
    }
  }
/*
  for (i = 0; i < pbods->faces.map.Extent(); i++)
    printf("  fMap %d: %d\n", i+1, fMap[i]);
 */
  for (i = 0; i < pbods->faces.map.Extent(); i++)
    if (fMap[i] == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: No Mapping for src Face %d (EG_mapBody)!\n",
               i+1);
      EG_free(fMap);
      return EGADS_ATTRERR;
    }

  // now Edges & Nodes
  eMap = (int *) EG_alloc(pbods->edges.map.Extent()*sizeof(int));
  if (eMap == NULL) {
    EG_free(fMap);
    return EGADS_MALLOC;
  }
  for (i = 0; i < pbods->edges.map.Extent(); i++) eMap[i] = 0;
  nMap = (int *) EG_alloc(pbods->nodes.map.Extent()*sizeof(int));
  if (nMap == NULL) {
    EG_free(eMap);
    EG_free(fMap);
    return EGADS_MALLOC;
  }
  for (i = 0; i < pbods->nodes.map.Extent(); i++) nMap[i] = 0;

  // make the Edge IDs
  stat = EG_getEdgeIDs(pbods, fAttr, &edgeIDs);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getEdgeIDs src ret = %d (EG_mapBody)!\n", stat);
    EG_free(nMap);
    EG_free(eMap);
    EG_free(fMap);
    return stat;
  }
  stat = EG_getEdgeIDs(pbodd, fAttr, &edgeIDd);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getEdgeIDs dst ret = %d (EG_mapBody)!\n", stat);
    for (i = 0; i < pbods->edges.map.Extent(); i++)
      if (edgeIDs[i].ID != NULL) EG_free(edgeIDs[i].ID);
    EG_free(edgeIDs);
    EG_free(nMap);
    EG_free(eMap);
    EG_free(fMap);
    return stat;
  }

  for (i = 0; i < pbods->edges.map.Extent(); i++) {
    egadsEdge     *pedge = (egadsEdge *) pbods->edges.objs[i]->blind;
    egadsNode     *pnode = (egadsNode *) pedge->nodes[0]->blind;
    TopoDS_Vertex vert   = TopoDS::Vertex(pnode->node);
    in0 = in1            = pbods->nodes.map.FindIndex(vert);
    if (in0 <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot find N0 in Edge %d (EG_mapBody)!\n",
               i+1);
      for (j = 0; j < pbods->edges.map.Extent(); j++)
        if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
      EG_free(edgeIDs);
      for (j = 0; j < pbodd->edges.map.Extent(); j++)
        if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
      EG_free(edgeIDd);
      EG_free(nMap);
      EG_free(eMap);
      EG_free(fMap);
      return EGADS_TOPOERR;
    }
    if (pbods->edges.objs[i]->mtype == TWONODE) {
      pnode = (egadsNode *) pedge->nodes[1]->blind;
      vert  = TopoDS::Vertex(pnode->node);
      in1   = pbods->nodes.map.FindIndex(vert);
      if (in1 <= 0) {
        if (outLevel > 0)
          printf(" EGADS Error: Cannot find N1 in Edge %d (EG_mapBody)!\n",
                 i+1);
        for (j = 0; j < pbods->edges.map.Extent(); j++)
          if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
        EG_free(edgeIDs);
        for (j = 0; j < pbodd->edges.map.Extent(); j++)
          if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
        EG_free(edgeIDd);
        EG_free(nMap);
        EG_free(eMap);
        EG_free(fMap);
        return EGADS_TOPOERR;
      }
    }
    for (j = 0; j < pbodd->edges.map.Extent(); j++) {
      if (strcmp(edgeIDs[i].ID, edgeIDd[j].ID) == 0) {
        if (eMap[i] == 0) {
          if (pbods->edges.objs[i]->mtype != pbodd->edges.objs[j]->mtype) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge Type Mismatch %d %d (EG_mapBody)!\n",
                     i+1, j+1);
            for (j = 0; j < pbods->edges.map.Extent(); j++)
              if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
            EG_free(edgeIDs);
            for (j = 0; j < pbodd->edges.map.Extent(); j++)
              if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
            EG_free(edgeIDd);
            EG_free(nMap);
            EG_free(eMap);
            EG_free(fMap);
            return EGADS_TOPOERR;
          }
          eMap[i]   = j+1;
          pedge     = (egadsEdge *) pbodd->edges.objs[j]->blind;
          pnode     = (egadsNode *) pedge->nodes[0]->blind;
          vert      = TopoDS::Vertex(pnode->node);
          jn0 = jn1 = pbodd->nodes.map.FindIndex(vert);
          if (jn0 <= 0) {
            if (outLevel > 0)
              printf(" EGADS Error: Cannot find N0 dstE %d (EG_mapBody)!\n",
                     j+1);
            for (j = 0; j < pbods->edges.map.Extent(); j++)
              if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
            EG_free(edgeIDs);
            for (j = 0; j < pbodd->edges.map.Extent(); j++)
              if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
            EG_free(edgeIDd);
            EG_free(nMap);
            EG_free(eMap);
            EG_free(fMap);
            return EGADS_TOPOERR;
          }
          if (pbodd->edges.objs[j]->mtype == TWONODE) {
            pnode = (egadsNode *) pedge->nodes[1]->blind;
            vert  = TopoDS::Vertex(pnode->node);
            jn1   = pbodd->nodes.map.FindIndex(vert);
            if (jn1 <= 0) {
              if (outLevel > 0)
                printf(" EGADS Error: Cannot find N1 dstE %d (EG_mapBody)!\n",
                       j+1);
              for (j = 0; j < pbods->edges.map.Extent(); j++)
                if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
              EG_free(edgeIDs);
              for (j = 0; j < pbodd->edges.map.Extent(); j++)
                if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
              EG_free(edgeIDd);
              EG_free(nMap);
              EG_free(eMap);
              EG_free(fMap);
              return EGADS_TOPOERR;
            }
          }
          if (nMap[in0-1] == 0) {
            nMap[in0-1] = jn0;
          } else if (nMap[in0-1] != jn0) {
            if (outLevel > 0)
              printf(" EGADS Error: Node attr Multip %d: %d %d (EG_mapBody)!\n",
                     in0, nMap[in0-1], jn0);
            for (j = 0; j < pbods->edges.map.Extent(); j++)
              if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
            EG_free(edgeIDs);
            for (j = 0; j < pbodd->edges.map.Extent(); j++)
              if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
            EG_free(edgeIDd);
            EG_free(nMap);
            EG_free(eMap);
            EG_free(fMap);
            return EGADS_ATTRERR;
          }
          if (nMap[in1-1] == 0) {
            nMap[in1-1] = jn1;
          } else if (nMap[in1-1] != jn1) {
            if (outLevel > 0)
              printf(" EGADS Error: Node Attr Multip %d: %d %d (EG_mapBody)!\n",
                     in1, nMap[in1-1], jn1);
            for (j = 0; j < pbods->edges.map.Extent(); j++)
              if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
            EG_free(edgeIDs);
            for (j = 0; j < pbodd->edges.map.Extent(); j++)
              if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
            EG_free(edgeIDd);
            EG_free(nMap);
            EG_free(eMap);
            EG_free(fMap);
            return EGADS_ATTRERR;
          }
        } else {
          if (outLevel > 0)
            printf(" EGADS Error: Edge Attr Multip %d: %d %d (EG_mapBody)!\n",
                   i+1, eMap[i], j+1);
          for (j = 0; j < pbods->edges.map.Extent(); j++)
            if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
          EG_free(edgeIDs);
          for (j = 0; j < pbodd->edges.map.Extent(); j++)
            if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
          EG_free(edgeIDd);
          EG_free(nMap);
          EG_free(eMap);
          EG_free(fMap);
          return EGADS_ATTRERR;
        }
      }
    }
  }
/*
  for (i = 0; i < pbods->edges.map.Extent(); i++) {
    j = eMap[i] - 1;
    printf("  eMap %d: %d  %d %d\n", i+1, eMap[i],
    pbods->edges.objs[i]->mtype, pbodd->edges.objs[j]->mtype);
  }
 */
  for (i = 0; i < pbods->edges.map.Extent(); i++)
    if (eMap[i] == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: No Mapping for src Edge %d (EG_mapBody)!\n",
               i+1);
      for (j = 0; j < pbods->edges.map.Extent(); j++)
        if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
      EG_free(edgeIDs);
      for (j = 0; j < pbodd->edges.map.Extent(); j++)
        if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
      EG_free(edgeIDd);
      EG_free(nMap);
      EG_free(eMap);
      EG_free(fMap);
      return EGADS_ATTRERR;
    }
/*
  for (i = 0; i < pbods->nodes.map.Extent(); i++)
    printf("  nMap %d: %d\n", i+1, nMap[i]);
 */
  for (i = 0; i < pbods->nodes.map.Extent(); i++)
    if (nMap[i] == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: No Mapping for src Node %d (EG_mapBody)!\n",
               i+1);
      for (j = 0; j < pbods->edges.map.Extent(); j++)
        if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
      EG_free(edgeIDs);
      for (j = 0; j < pbodd->edges.map.Extent(); j++)
        if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
      EG_free(edgeIDd);
      EG_free(nMap);
      EG_free(eMap);
      EG_free(fMap);
      return EGADS_ATTRERR;
    }

  for (j = 0; j < pbods->edges.map.Extent(); j++)
    if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
  EG_free(edgeIDs);
  for (j = 0; j < pbodd->edges.map.Extent(); j++)
    if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
  EG_free(edgeIDd);

  *pnMap = nMap;
  *peMap = eMap;
  *pfMap = fMap;
  return EGADS_SUCCESS;
}


int
EG_mapBody(const egObject  *sBody, const egObject *dBody, const char *fAttr,
                 egObject **mBody)
{
  int outLevel, stat, *nMap = NULL, *eMap = NULL, *fMap = NULL, exacTopo = 0;

  *mBody = NULL;
  if  (sBody == NULL)                return EGADS_NULLOBJ;
  if  (sBody->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if  (sBody->oclass != BODY)        return EGADS_NOTBODY;
  if  (sBody->blind == NULL)         return EGADS_NODATA;
  if ((sBody->mtype != FACEBODY)  &&
      (sBody->mtype != SHEETBODY) &&
      (sBody->mtype != SOLIDBODY))   return EGADS_TOPOERR;
  if  (EG_sameThread(sBody))         return EGADS_CNTXTHRD;
  if  (dBody == NULL)                return EGADS_NULLOBJ;
  if  (dBody->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if  (dBody->oclass != BODY)        return EGADS_NOTBODY;
  if  (dBody->blind == NULL)         return EGADS_NODATA;
  if  (dBody->mtype != sBody->mtype) return EGADS_TOPOERR;
  if  (EG_sameThread(dBody))         return EGADS_CNTXTHRD;
  egadsBody *pbodd = (egadsBody *) dBody->blind;
  outLevel = EG_outLevel(sBody);

  stat = EG_sameBodyTopo(sBody, dBody);
  if (stat == EGADS_TOPOCNT) return stat;
  if (stat == EGADS_SUCCESS) exacTopo = 1;
  if (sBody->mtype == FACEBODY) {
    if (exacTopo != 1) return EGADS_TOPOERR;
    /* fix spline surface? */
    goto fixBSplSur;
  }

  // differing topology -- make mappings
  if (exacTopo == 0) {
    stat = EG_bodyMapping(sBody, dBody, fAttr, &nMap, &eMap, &fMap);
    if (stat != EGADS_SUCCESS) return stat;
  }

  // look at BSpline surfaces for potential Face updates

fixBSplSur:
  // can we just use the destination Body?
  if ((*mBody == NULL) && (exacTopo == 1)) return EGADS_SUCCESS;

  // make a new body if not already done
  if (*mBody == NULL) {
    stat = EG_copyObject(dBody, NULL, mBody);
    if (stat != EGADS_SUCCESS) {
      if (nMap != NULL) EG_free(nMap);
      if (eMap != NULL) EG_free(eMap);
      if (fMap != NULL) EG_free(fMap);
      return stat;
    }
  }

  // put the mapping attributes on the body
  if (nMap != NULL) {
    stat = EG_attributeAdd(*mBody, ".nMap", ATTRINT, pbodd->nodes.map.Extent(),
                           nMap, NULL, NULL);
    EG_free(nMap);
    if (stat != EGADS_SUCCESS) {
      EG_deleteObject(*mBody);
      *mBody = NULL;
      if (outLevel > 0)
        printf(" EGADS Warning: Attribute Node Mapping is %d (EG_mapBody)!\n",
               stat);
      EG_free(eMap);
      EG_free(fMap);
      return stat;
    }
  }
  if (eMap != NULL) {
    stat = EG_attributeAdd(*mBody, ".eMap", ATTRINT, pbodd->edges.map.Extent(),
                           eMap, NULL, NULL);
    EG_free(eMap);
    if (stat != EGADS_SUCCESS) {
      EG_deleteObject(*mBody);
      *mBody = NULL;
      if (outLevel > 0)
        printf(" EGADS Warning: Attribute Edge Mapping is %d (EG_mapBody)!\n",
               stat);
      EG_free(fMap);
      return stat;
    }
  }
  if (fMap != NULL) {
    stat = EG_attributeAdd(*mBody, ".fMap", ATTRINT, pbodd->faces.map.Extent(),
                           fMap, NULL, NULL);
    EG_free(fMap);
    if (stat != EGADS_SUCCESS) {
      EG_deleteObject(*mBody);
      *mBody = NULL;
      if (outLevel > 0)
        printf(" EGADS Warning: Attribute Face Mapping is %d (EG_mapBody)!\n",
               stat);
      return stat;
    }
  }

  return EGADS_SUCCESS;
}
