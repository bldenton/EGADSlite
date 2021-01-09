/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             High-Level Functions
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
#include "egadsStack.h"

#define OCC_EXTRUDE
#define OCC_ROTATE

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)          (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


  typedef struct {
    int sense;                  /* sense use for the loop construction */
    int index;                  /* index in loop */
    int lIndex;                 /* loop index */
  } loopInfo;


  extern "C" int  EG_destroyTopology( egObject *topo );
  extern "C" int  EG_fullAttrs( const egObject *obj );
  extern "C" int  EG_sewFaces( int nobj, const egObject **objs, double toler,
                               int flag, egObject **result );
  extern "C" int  EG_getBodyTopos( const egObject *body, /*@null@*/ egObject *src,
                                   int oclass, int *ntopo, egObject ***topos );
  extern "C" int  EG_attributeNum( const egObject *obj, int *num );
  extern "C" int  EG_attributeDup( const egObject *src, ego dst );
  extern "C" int  EG_attributeAdd( egObject *obj, const char *name, int type,
                                   int len, /*@null@*/ const int    *ints,
                                            /*@null@*/ const double *reals,
                                            /*@null@*/ const char   *str );
  extern "C" int  EG_attributeDel( egObject *obj, /*@null@*/ const char *name );
  extern "C" int  EG_attributeRet( const egObject *obj, const char *name,
                                   int *atype, int *len,
                                   /*@null@*/ const int    **ints,
                                   /*@null@*/ const double **reals,
                                   /*@null@*/ const char   **str );
  extern "C" int  EG_isSame( const egObject *geom1, const egObject *geom2 );
  extern "C" int  EG_isEquivalent( const egObject *top1, const egObject *top2 );
  extern "C" int  EG_getRange( const egObject *geom, double *range, int *per );
  extern "C" int  EG_evaluate( const egObject *geom,
                               /*@null@*/ const double *param, double *results );
             int  EG_evaluate( const egObject *geom,
                               /*@null@*/ const SurrealS<1> *param,
                               SurrealS<1> *result );
  extern "C" int  EG_tolerance( const egObject *topo, double *tol );
  extern "C" int  EG_getTopology( const egObject *topo, egObject **geom,
                                  int *oclass, int *type,
                                  /*@null@*/ double *limits,
                                  int *nChildren, egObject ***children,
                                  int **sense );
  extern "C" int  EG_copyObject( const egObject *object, /*@null@*/ void *oform,
                                 egObject **copy );
  extern "C" int  EG_splitBody( const egObject *body, int nedge,
                                const egObject **facEdg, egObject **result );

  extern "C" int  EG_generalBoolean( egObject *src, egObject *tool, int oper,
                                     double tol, egObject **model );
  extern "C" int  EG_fuseSheets( const egObject *src, const egObject *tool,
                                 egObject **sheet );
  extern "C" int  EG_solidBoolean( const egObject *src, const egObject *tool,
                                   int oper, egObject **model );
  extern "C" int  EG_intersection( const egObject *src, const egObject *tool,
                                   int *nedge, /*@null@*/ egObject ***facEdg,
                                   egObject **model );
  extern "C" int  EG_imprintBody( const egObject *src, int nedge,
                                  const egObject **facEdg, egObject **result );
  extern "C" int  EG_filletBody( const egObject *src, int nedge,
                                 const egObject **edges, double radius,
                                       egObject **result, int **faceMap );
  extern "C" int  EG_chamferBody( const egObject *src, int nedge,
                                  const egObject **edges, const egObject **faces,
                                  double dis1, double dis2, egObject **result,
                                  int **faceMap );
  extern "C" int  EG_offsetEdge( const egObject *src,  const egObject *edge,
                                 const egObject *face, double offset,
                                 egObject **result );
  extern "C" int  EG_hollowBody( const egObject *src, int nface,
                                 const egObject **faces, double offset, int join,
                                       egObject **result, int **faceMap );
  extern "C" int  EG_extrude( const egObject *src, double dist,
                              const double *dir, egObject **result );
  extern "C" int  EG_extrude_dot( egObject *body, const egObject *src,
                                  double dist, double dist_dot,
                                  const double *dir, const double *dir_dot );

  extern "C" int  EG_rotate( const egObject *src, double angle,
                             const double *axis, egObject **result );
  extern "C" int  EG_sweep( const egObject *src, const egObject *spine, int mode,
                                  egObject **result );
  extern "C" int  EG_loft( int nsec, const egObject **secs, int opt,
                                           egObject **result );

  extern     int  EG_traverseBody( egObject *context, int i, egObject *bobj,
                                   egObject *topObj, egadsBody *body,
                                   int *nerr );
  extern     int  EG_attriBodyDup( const egObject *src, egObject *dst );
  extern     void EG_completePCurve( egObject *g, Handle(Geom2d_Curve) &hCurv );
  extern     void EG_completeSurf(   egObject *g, Handle(Geom_Surface) &hSurf );

  extern "C" int  EG_makeGeometry( egObject *context, int oclass, int mtype,
                                   /*@null@*/ egObject *refGeo, const int *ivec,
                                   const double *rvec, egObject **geom );
  extern "C" int  EG_makeTopology( egObject *context, /*@null@*/ egObject *geom,
                                   int oclass, int mtype, /*@null@*/ double *limits,
                                   int nChildren, /*@null@*/ egObject **children,
                                   /*@null@*/ int *senses, egObject **topo );

  extern "C" int  EG_makeTransform( egObject *context, const double *xform,
                                    egObject **oform );
  extern "C" int  EG_getGeometry( const egObject *geom, int *oclass, int *type,
                                  egObject **rGeom, /*@null@*/ int **ivec,
                                  /*@null@*/ double **rvec );
  extern     int  EG_getGeometry( const ego geom, int *oclass, int *mtype,
                                  ego *refGeom, int **ivec, SurrealS<1> **rvec );

  extern "C" int  EG_hasGeometry_dot( const egObject *geom );
             int  EG_setGeometry_dot( egObject *geom, int oclass, int mtype,
                                      /*@null@*/ const int *ivec,
                                      SurrealS<1> *data_dot );
             int  EG_copyGeometry_dot( const egObject *obj,
                                       /*@null@*/ const SurrealS<1> *xform,
                                       ego copy );


static void
EG_attrFixup(int outLevel, egObject *body, const egObject *src)
{
  int      i, j, na, stat, nface, nfsrc, iper, ipers;
  double   uvbox[4], uvbsrc[4], fuzz=1.0e-4;
  egObject **faces, **fsrcs;

  stat = EG_getBodyTopos(body, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Internal EG_attrFixup: EG_getBodyTopos body = %d!\n", stat);
    return;
  }
  stat = EG_getBodyTopos(src, NULL, FACE, &nfsrc, &fsrcs);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Internal EG_attrFixup: EG_getBodyTopos src  = %d!\n", stat);
    EG_free(faces);
    return;
  }

  for (j = 0; j < nface; j++) {
    stat = EG_attributeNum(faces[j], &na);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Internal EG_attrFixup: EG_attributeNum %d = %d!\n",
             j+1, stat);
      continue;
    }
    if (na != 0) continue;

    for (i = 0; i < nfsrc; i++) {
      /* compare geometry */
      stat = EG_isSame(faces[j], fsrcs[i]);
      if (stat != EGADS_SUCCESS) continue;
      /* compare uvbox */
      stat = EG_getRange(fsrcs[i], uvbsrc, &ipers);
      if (stat != EGADS_SUCCESS) continue;
      stat = EG_getRange(faces[j], uvbox,  &iper);
      if (stat != EGADS_SUCCESS) continue;
      if (uvbox[0] <  uvbsrc[0]-fuzz) continue;
      if (uvbox[1] >  uvbsrc[1]+fuzz) continue;
      if (uvbox[2] <  uvbsrc[2]-fuzz) continue;
      if (uvbox[3] >  uvbsrc[3]+fuzz) continue;
      EG_attributeNum(fsrcs[i], &na);
      if (na == 0) break;
      if (outLevel > 1)
        printf(" EG_attrFixup -- Face %d: Fixup with %d's atrributes!\n",
               j+1, i+1);
      stat = EG_attributeDup(fsrcs[i], faces[j]);
      if (stat != EGADS_SUCCESS)
        printf(" EGADS Internal EG_attrFixup: EG_attributeDup = %d!\n", stat);
      break;
    }
  }

  EG_free(fsrcs);
  EG_free(faces);
}


#if CASVER >= 700
static int
EG_matchGeneral(BRepAlgoAPI_BuilderAlgo& BSO, egObject *src,
                egObject *tool, TopoDS_Shape result,
                egObject ***emapping, egObject ***fmapping)
{
  int         i, j, index, nbody, fullAttrs;
  egObject    **emap = NULL, **fmap = NULL, **bodies;
  TopoDS_Edge edge, genedge;
  TopoDS_Face face, genface;

  *emapping = NULL;
  *fmapping = NULL;
  fullAttrs = EG_fullAttrs(src);

  TopTools_IndexedMapOfShape rmape, rmap;
  TopExp::MapShapes(result, TopAbs_EDGE, rmape);
  if (fullAttrs != 0) {
    if (rmape.Extent() != 0) {
      emap = (egObject **) EG_alloc(rmape.Extent()*sizeof(egObject *));
      if (emap == NULL) return EGADS_MALLOC;
    }
  }
  TopExp::MapShapes(result, TopAbs_FACE, rmap);
  if (rmap.Extent() != 0) {
    fmap = (egObject **) EG_alloc(rmap.Extent()*sizeof(egObject *));
    if (fmap == NULL) {
      if (emap != NULL) EG_free(emap);
      return EGADS_MALLOC;
    }
  }
  if ((emap == NULL) && (fmap == NULL)) return EGADS_CONSTERR;

  if (emap != NULL)
    for (i = 0; i < rmape.Extent(); i++) emap[i] = NULL;
  if (fmap != NULL)
    for (i = 0; i < rmap.Extent();  i++) fmap[i] = NULL;

  if (src->oclass == MODEL) {
    egadsModel *pmdl = (egadsModel *) src->blind;
    nbody            = pmdl->nbody;
    bodies           = pmdl->bodies;
  } else {
    nbody            = 1;
    bodies           = &src;
  }

  for (i = 0; i < nbody; i++) {
    if (bodies[i]        == NULL) return EGADS_NULLOBJ;
    if (bodies[i]->blind == NULL) return EGADS_NODATA;
    egadsBody *pbody = (egadsBody *) bodies[i]->blind;
    if (emap != NULL)
      for (j = 1; j <= pbody->edges.map.Extent(); j++) {
        edge = TopoDS::Edge(pbody->edges.map(j));
        if (BSO.IsDeleted(edge)) continue;
        const TopTools_ListOfShape& listEdges = BSO.Modified(edge);
        if (listEdges.Extent() > 0) {
          /* modified Edges */
          TopTools_ListIteratorOfListOfShape it(listEdges);
          for (; it.More(); it.Next()) {
            genedge = TopoDS::Edge(it.Value());
            index   = rmape.FindIndex(genedge);
            if (index > 0) emap[index-1] = pbody->edges.objs[j-1];
          }
        }
      }
    if (fmap != NULL)
      for (j = 1; j <= pbody->faces.map.Extent(); j++) {
        face = TopoDS::Face(pbody->faces.map(j));
        if (BSO.IsDeleted(face)) continue;
        const TopTools_ListOfShape& listFaces = BSO.Modified(face);
        if (listFaces.Extent() > 0) {
          /* modified Faces */
          TopTools_ListIteratorOfListOfShape it(listFaces);
          for (; it.More(); it.Next()) {
            genface = TopoDS::Face(it.Value());
            index   = rmap.FindIndex(genface);
            if (index > 0) fmap[index-1] = pbody->faces.objs[j-1];
          }
        }
      }
    }

  if (tool->oclass == MODEL) {
    egadsModel *pmdl = (egadsModel *) tool->blind;
    nbody            = pmdl->nbody;
    bodies           = pmdl->bodies;
  } else {
    nbody  = 1;
    bodies = &tool;
  }

  for (i = 0; i < nbody; i++) {
    if (bodies[i]        == NULL) return EGADS_NULLOBJ;
    if (bodies[i]->blind == NULL) return EGADS_NODATA;
    egadsBody *pbody = (egadsBody *) bodies[i]->blind;
    if (emap != NULL)
      for (j = 1; j <= pbody->edges.map.Extent(); j++) {
        edge = TopoDS::Edge(pbody->edges.map(j));
        if (BSO.IsDeleted(edge)) continue;
        const TopTools_ListOfShape& listEdges = BSO.Modified(edge);
        if (listEdges.Extent() > 0) {
          /* modified Edges */
          TopTools_ListIteratorOfListOfShape it(listEdges);
          for (; it.More(); it.Next()) {
            genedge = TopoDS::Edge(it.Value());
            index   = rmape.FindIndex(genedge);
            if (index > 0) emap[index-1] = pbody->edges.objs[j-1];
          }
        }
      }
    if (fmap != NULL)
      for (j = 1; j <= pbody->faces.map.Extent(); j++) {
        face = TopoDS::Face(pbody->faces.map(j));
        if (BSO.IsDeleted(face)) continue;
        const TopTools_ListOfShape& listFaces = BSO.Modified(face);
        if (listFaces.Extent() > 0) {
          /* modified Faces */
          TopTools_ListIteratorOfListOfShape it(listFaces);
          for (; it.More(); it.Next()) {
            genface = TopoDS::Face(it.Value());
            index   = rmap.FindIndex(genface);
            if (index > 0) fmap[index-1] = pbody->faces.objs[j-1];
          }
        }
      }
  }

  *emapping = emap;
  *fmapping = fmap;
  return EGADS_SUCCESS;
}
#endif


static void
EG_matchMdlFace(BRepAlgoAPI_BooleanOperation& BSO, TopoDS_Shape src,
                int iface, TopoDS_Shape tool, TopoDS_Shape result,
                int fullAttrs, int **emapping, int **fmapping)
{
  int                        i, j, nf, ne, *emap = NULL, *map = NULL;
  TopoDS_Edge                edge, genedge;
  TopoDS_Face                face, genface;
  TopTools_IndexedMapOfShape rmap, smap, tmap, rmape, smape, tmape;

  *emapping = NULL;
  *fmapping = NULL;
  if (fullAttrs != 0) {
    TopExp::MapShapes(result, TopAbs_EDGE, rmape);
    TopExp::MapShapes(src,    TopAbs_EDGE, smape);
    ne = rmape.Extent();
  } else {
    ne = 0;
  }
  TopExp::MapShapes(result, TopAbs_FACE, rmap);
  TopExp::MapShapes(src,    TopAbs_FACE, smap);
  nf = rmap.Extent();
  if ((nf == 0) && (ne == 0)) return;

  if (ne != 0) {
    emap = (int *) EG_alloc(ne*sizeof(int));
    if (emap == NULL) return;
    for (i = 0; i < ne; i++) emap[i] = 0;
    *emapping = emap;
  }

  if (nf != 0) {
    map = (int *) EG_alloc(nf*sizeof(int));
    if (map == NULL) return;
    for (i = 0; i < nf; i++) map[i] = 0;
    *fmapping = map;
  }

  if (ne != 0)
    for (i = 1; i <= smape.Extent(); i++) {
      edge = TopoDS::Edge(smape(i));
      if (BSO.IsDeleted(edge)) continue;
      const TopTools_ListOfShape& listEdges = BSO.Modified(edge);
      if (listEdges.Extent() > 0) {
        /* modified Edges */
        TopTools_ListIteratorOfListOfShape it(listEdges);
        for (; it.More(); it.Next()) {
          genedge = TopoDS::Edge(it.Value());
          j = rmape.FindIndex(genedge);
          if (j > 0) emap[j-1] = i;
        }
      } else {
        j = rmape.FindIndex(edge);
        if (j > 0) emap[j-1] = i;
      }
    }

  if (nf != 0)
    for (i = 1; i <= smap.Extent(); i++) {
      face = TopoDS::Face(smap(i));
      if (BSO.IsDeleted(face)) continue;
      const TopTools_ListOfShape& listFaces = BSO.Modified(face);
      if (listFaces.Extent() > 0) {
        /* modified Faces */
        TopTools_ListIteratorOfListOfShape it(listFaces);
        for (; it.More(); it.Next()) {
          genface = TopoDS::Face(it.Value());
          j = rmap.FindIndex(genface);
          if (j > 0) map[j-1] = i;
        }
      } else {
        j = rmap.FindIndex(face);
        if (j > 0) map[j-1] = i;
      }
    }

  if (iface == 0) {

    if (ne != 0) {
      TopExp::MapShapes(tool, TopAbs_EDGE, tmape);
      for (i = 1; i <= tmape.Extent(); i++) {
        edge = TopoDS::Edge(tmape(i));
        if (BSO.IsDeleted(edge)) continue;
        const TopTools_ListOfShape& listEdges = BSO.Modified(edge);
        if (listEdges.Extent() > 0) {
          /* modified Edges */
          TopTools_ListIteratorOfListOfShape it(listEdges);
          for (; it.More(); it.Next()) {
            genedge = TopoDS::Edge(it.Value());
            j = rmape.FindIndex(genedge);
            if (j > 0) emap[j-1] = -i;
          }
        }
      }
    }

    if (nf != 0) {
      TopExp::MapShapes(tool, TopAbs_FACE, tmap);
      for (i = 1; i <= tmap.Extent(); i++) {
        face = TopoDS::Face(tmap(i));
        if (BSO.IsDeleted(face)) continue;
        const TopTools_ListOfShape& listFaces = BSO.Modified(face);
        if (listFaces.Extent() > 0) {
          /* modified Faces */
          TopTools_ListIteratorOfListOfShape it(listFaces);
          for (; it.More(); it.Next()) {
            genface = TopoDS::Face(it.Value());
            j = rmap.FindIndex(genface);
            if (j > 0) map[j-1] = -i;
          }
        }
      }
    }

  } else {

    /* NOTE: No persistence of Edges! */
    if (nf != 0) {
      face = TopoDS::Face(tool);
      if (!BSO.IsDeleted(face)) {
        const TopTools_ListOfShape& listFaces = BSO.Modified(face);
        if (listFaces.Extent() > 0) {
          /* modified Faces */
          TopTools_ListIteratorOfListOfShape it(listFaces);
          for (; it.More(); it.Next()) {
            genface = TopoDS::Face(it.Value());
            j = rmap.FindIndex(genface);
            if (j > 0) map[j-1] = -1;
          }
        }
      }
    }

  }

}


static void
EG_matchFaces(BRepAlgoAPI_BooleanOperation& BSO, const egObject *src,
              const egObject *tool, TopAbs_ShapeEnum type, TopoDS_Shape& result,
              int ***emapping, int ***fmapping)
{
  int             i, j, k, ns, nface, fullAttrs, **map, **emap = NULL;
  egadsBody       *pbods, *pbodt;
  const egObject  *oface = NULL;
  TopoDS_Shape    solid;
  TopoDS_Edge     edge, genedge;
  TopoDS_Face     face, genface;
  TopExp_Explorer Exp;

  *emapping = NULL;
  *fmapping = NULL;
  fullAttrs = EG_fullAttrs(src);

  pbods     = (egadsBody *) src->blind;
  if (((tool->oclass == FACE) ||
      ((tool->oclass == BODY) && (tool->mtype == FACEBODY)))) {
    if (tool->oclass == FACE) {
      pbodt = NULL;
      oface = tool;
    } else {
      pbodt = (egadsBody *) tool->blind;
      oface = pbodt->faces.objs[0];
    }
  } else {
    pbodt = (egadsBody *) tool->blind;
  }

  ns = 0;
  for (Exp.Init(result, type); Exp.More(); Exp.Next()) ns++;
  if (ns == 0) return;
  map = (int **) EG_alloc(ns*sizeof(int *));
  if (map == NULL) return;
  for (i = 0; i < ns; i++) map[i] = NULL;
  if (fullAttrs != 0) {
    emap = (int **) EG_alloc(ns*sizeof(int *));
    if (emap == NULL) {
      EG_free(map);
      return;
    }
    for (i = 0; i < ns; i++) emap[i] = NULL;
  }

  k = 0;
  for (Exp.Init(result, type); Exp.More(); Exp.Next()) {
    solid = Exp.Current();
    TopTools_IndexedMapOfShape MapE, MapF;
    TopExp::MapShapes(solid, TopAbs_EDGE, MapE);
    if (fullAttrs != 0) {
      int nedge = MapE.Extent();
      if (nedge > 0) emap[k] = (int *) EG_alloc(nedge*sizeof(int));
      if (emap[k] != NULL)
        for (j = 0; j < nedge; j++) emap[k][j] = 0;
    }
    TopExp::MapShapes(solid, TopAbs_FACE, MapF);
    nface = MapF.Extent();
    if (nface > 0) map[k] = (int *) EG_alloc(nface*sizeof(int));
    if (map[k] != NULL)
      for (j = 0; j < nface; j++) map[k][j] = 0;
    k++;
  }

  *emapping = emap;
  *fmapping =  map;

  /* look at source shape */

  if (fullAttrs != 0) {
    for (i = 1; i <= pbods->edges.map.Extent(); i++) {
      edge = TopoDS::Edge(pbods->edges.map(i));
      if (BSO.IsDeleted(edge)) continue;
      const TopTools_ListOfShape& listEdges = BSO.Modified(edge);
      if (listEdges.Extent() > 0) {
        /* modified Edges */
        TopTools_ListIteratorOfListOfShape it(listEdges);
        for (; it.More(); it.Next()) {
          genedge = TopoDS::Edge(it.Value());
          k       = 0;
          for (Exp.Init(result, type); Exp.More(); Exp.Next()) {
            solid = Exp.Current();
            TopTools_IndexedMapOfShape MapE;
            TopExp::MapShapes(solid, TopAbs_EDGE, MapE);
            if (emap[k] != NULL) {
              j = MapE.FindIndex(genedge);
              if (j > 0) emap[k][j-1] = i;
            }
            k++;
          }
        }
      }
    }
  }

  for (i = 1; i <= pbods->faces.map.Extent(); i++) {
    face = TopoDS::Face(pbods->faces.map(i));
    if (BSO.IsDeleted(face)) continue;
    const TopTools_ListOfShape& listFaces = BSO.Modified(face);
    if (listFaces.Extent() > 0) {
      /* modified Faces */
      TopTools_ListIteratorOfListOfShape it(listFaces);
      for (; it.More(); it.Next()) {
        genface = TopoDS::Face(it.Value());
        k       = 0;
        for (Exp.Init(result, type); Exp.More(); Exp.Next()) {
          solid = Exp.Current();
          TopTools_IndexedMapOfShape MapF;
          TopExp::MapShapes(solid, TopAbs_FACE, MapF);
          if (map[k] != NULL) {
            j = MapF.FindIndex(genface);
            if (j > 0) map[k][j-1] = i;
          }
          k++;
        }
      }
    }
  }

  /* look at tool shape */

  if (fullAttrs != 0)
    if (tool->oclass != FACE) {
      /* NOTE: No Edge persistence for Face but OK for all others */
      for (i = 1; i <= pbodt->edges.map.Extent(); i++) {
        edge = TopoDS::Edge(pbodt->edges.map(i));
        if (BSO.IsDeleted(edge)) continue;
        const TopTools_ListOfShape& listEdges = BSO.Modified(edge);
        if (listEdges.Extent() > 0) {
          /* modified Edges */
          TopTools_ListIteratorOfListOfShape it(listEdges);
          for (; it.More(); it.Next()) {
            genedge = TopoDS::Edge(it.Value());
            k       = 0;
            for (Exp.Init(result, type); Exp.More(); Exp.Next()) {
              solid = Exp.Current();
              TopTools_IndexedMapOfShape MapE;
              TopExp::MapShapes(solid, TopAbs_EDGE, MapE);
              if (emap[k] != NULL) {
                j = MapE.FindIndex(genedge);
                if (j > 0) emap[k][j-1] = -i;
              }
              k++;
            }
          }
        }
      }
    }

  if (oface == NULL) {

    for (i = 1; i <= pbodt->faces.map.Extent(); i++) {
      face = TopoDS::Face(pbodt->faces.map(i));
      if (BSO.IsDeleted(face)) continue;
      const TopTools_ListOfShape& listFaces = BSO.Modified(face);
      if (listFaces.Extent() > 0) {
        /* modified Faces */
        TopTools_ListIteratorOfListOfShape it(listFaces);
        for (; it.More(); it.Next()) {
          genface = TopoDS::Face(it.Value());
          k       = 0;
          for (Exp.Init(result, type); Exp.More(); Exp.Next()) {
            solid = Exp.Current();
            TopTools_IndexedMapOfShape MapF;
            TopExp::MapShapes(solid, TopAbs_FACE, MapF);
            if (map[k] != NULL) {
              j = MapF.FindIndex(genface);
              if (j > 0) map[k][j-1] = -i;
            }
            k++;
          }
        }
      }
    }

  } else {

    egadsFace *pface = (egadsFace *) oface->blind;
    face = pface->face;
    k    = 0;
    for (Exp.Init(result, type); Exp.More(); Exp.Next()) {
      solid = Exp.Current();
      TopTools_IndexedMapOfShape MapF;
      TopExp::MapShapes(solid, TopAbs_FACE, MapF);
      if (map[k] != NULL) {
        j = MapF.FindIndex(genface);
        if (j > 0) map[k][j-1] = -1;
      }
      k++;
    }
    if (!BSO.IsDeleted(face)) {
      const TopTools_ListOfShape& listFaces = BSO.Modified(face);
      if (listFaces.Extent() > 0) {
        /* modified Faces */
        TopTools_ListIteratorOfListOfShape it(listFaces);
        for (; it.More(); it.Next()) {
          genface = TopoDS::Face(it.Value());
          k       = 0;
          for (Exp.Init(result, type); Exp.More(); Exp.Next()) {
            solid = Exp.Current();
            TopTools_IndexedMapOfShape MapF;
            TopExp::MapShapes(solid, TopAbs_FACE, MapF);
            if (map[k] != NULL) {
              j = MapF.FindIndex(genface);
              if (j > 0) map[k][j-1] = -1;
            }
            k++;
          }
        }
      }
    }
  }

}


int
EG_generalBoolean(egObject *src, egObject *tool, int oper, double tol,
                  egObject **model)
{
#if CASVER < 700
  printf(" EGADS Error: Not implemented for this Rev of OCC (EG_generalBoolean)!\n");
  return EGADS_NOTFOUND;
#else
  int          i, j, index, nerr, outLevel, stat;
  double       toler = Precision::Confusion();
  egObject     *context, *omodel, **emap = NULL, **fmap = NULL;
  TopoDS_Shape result;

  *model = NULL;
  if  (src == NULL)               return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (EG_sameThread(src))        return EGADS_CNTXTHRD;
  if ((src->oclass != MODEL) && (src->oclass != BODY))
                                  return EGADS_NOTBODY;
  if  (src->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);
  if (tol > 0.0) toler = tol;
  if ((oper != SUBTRACTION) && (oper != INTERSECTION) &&
      (oper != FUSION)      && (oper != SPLITTER)) {
    if (outLevel > 0)
      printf(" EGADS Error: BAD Operator = %d (EG_generalBoolean)!\n",
             oper);
    return EGADS_RANGERR;
  }
  if (tool == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Tool (EG_generalBoolean)!\n");
    return EGADS_NULLOBJ;
  }
  if (tool->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Tool is not an EGO (EG_generalBoolean)!\n");
    return EGADS_NOTOBJ;
  }
  if (tool->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Tool has no data (EG_generalBoolean)!\n");
    return EGADS_NODATA;
  }
  if (EG_sameThread(tool)) {
    if (outLevel > 0)
      printf(" EGADS Error: Bad Thread on Tool (EG_generalBoolean)!\n");
    return EGADS_CNTXTHRD;
  }
  if (EG_context(tool) != context) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_generalBoolean)!\n");
    return EGADS_MIXCNTX;
  }
  if ((tool->oclass != MODEL) && (tool->oclass != BODY)) {
    if (outLevel > 0)
      printf(" EGADS Error: Tool is not a BODY or MODEL (EG_generalBoolean)!\n");
    return EGADS_NOTBODY;
  }

  TopTools_ListOfShape sList, tList;
  if (src->oclass == MODEL) {
    egadsModel *pmdl = (egadsModel *) src->blind;
    sList.Append(TopoDS::Compound(pmdl->shape));
  } else {
    egadsBody *pbody = (egadsBody *) src->blind;
    if (src->mtype == SOLIDBODY) {
      sList.Append(TopoDS::Solid(pbody->shape));
    } else {
      sList.Append(pbody->shape);
    }
  }
  if (tool->oclass == MODEL) {
    egadsModel *pmdl = (egadsModel *) tool->blind;
    tList.Append(TopoDS::Compound(pmdl->shape));
  } else {
    egadsBody *pbody = (egadsBody *) tool->blind;
    if (tool->mtype == SOLIDBODY) {
      tList.Append(TopoDS::Solid(pbody->shape));
    } else {
      tList.Append(pbody->shape);
    }
  }

  if (oper == INTERSECTION) {

    try {
      BRepAlgoAPI_Common BSO;
      BSO.SetArguments(sList);
      BSO.SetTools(tList);
//    BSO.SetGlue(BOPAlgo_GlueFull);
      BSO.SetFuzzyValue(toler);
      BSO.SetNonDestructive(Standard_True);
      BSO.SetUseOBB(Standard_True);
      BSO.Build();
#if !(defined(__APPLE__) && !defined(__clang__))
      if (BSO.HasErrors()) {
        BSO.DumpErrors(std::cout);
        return EGADS_GEOMERR;
      }
#endif
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Intersection (EG_generalBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      stat   = EG_matchGeneral(BSO, src, tool, result, &emap, &fmap);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Error: SBO Intersection Exception (EG_generalBoolean)!\n");
      printf("              %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Error: SBO Intersection Exception (EG_generalBoolean)!\n");
      return EGADS_GEOMERR;
    }

  } else if (oper == SUBTRACTION) {

    try {
      BRepAlgoAPI_Cut BSO;
      BSO.SetArguments(sList);
      BSO.SetTools(tList);
//    BSO.SetGlue(BOPAlgo_GlueFull);
      BSO.SetFuzzyValue(toler);
      BSO.SetNonDestructive(Standard_True);
      BSO.SetUseOBB(Standard_True);
      BSO.Build();
#if !(defined(__APPLE__) && !defined(__clang__))
      if (BSO.HasErrors()) {
        BSO.DumpErrors(std::cout);
        return EGADS_GEOMERR;
      }
#endif
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Subtraction (EG_generalBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      stat   = EG_matchGeneral(BSO, src, tool, result, &emap, &fmap);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Error: SBO Subtraction Exception (EG_generalBoolean)!\n");
      printf("              %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Error: SBO Subtraction Exception (EG_generalBoolean)!\n");
      return EGADS_GEOMERR;
    }

  } else if (oper == FUSION) {

    try {
      BRepAlgoAPI_Fuse BSO;
      BSO.SetArguments(sList);
      BSO.SetTools(tList);
//    BSO.SetGlue(BOPAlgo_GlueFull);
      BSO.SetFuzzyValue(toler);
      BSO.SetNonDestructive(Standard_True);
      BSO.SetUseOBB(Standard_True);
      BSO.Build();
#if !(defined(__APPLE__) && !defined(__clang__))
      if (BSO.HasErrors()) {
        BSO.DumpErrors(std::cout);
        return EGADS_GEOMERR;
      }
#endif
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Fusion (EG_generalBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      stat   = EG_matchGeneral(BSO, src, tool, result, &emap, &fmap);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Error: SBO Fusion Exception (EG_generalBoolean)!\n");
      printf("              %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Error: SBO Fusion Exception (EG_generalBoolean)!\n");
      return EGADS_GEOMERR;
    }

  } else {

    try {
      BRepAlgoAPI_Splitter BSO;
      BSO.SetArguments(sList);
      BSO.SetTools(tList);
//    BSO.SetGlue(BOPAlgo_GlueFull);
      BSO.SetFuzzyValue(toler);
      BSO.SetNonDestructive(Standard_True);
      BSO.SetUseOBB(Standard_True);
      BSO.Build();
#if !(defined(__APPLE__) && !defined(__clang__))
      if (BSO.HasErrors()) {
        BSO.DumpErrors(std::cout);
        return EGADS_GEOMERR;
      }
#endif
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Split (EG_generalBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      stat   = EG_matchGeneral(BSO, src, tool, result, &emap, &fmap);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Error: SBO Split Exception (EG_generalBoolean)!\n");
      printf("              %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Error: SBO Split Exception (EG_generalBoolean)!\n");
      return EGADS_GEOMERR;
    }

  }
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Attribute Mapping failed = %d (EG_generalBoolean)!\n",
             stat);
    return EGADS_CONSTERR;
  }

  // parse the result
  int nWire  = 0;
  int nFace  = 0;
  int nSheet = 0;
  int nSolid = 0;
  TopExp_Explorer Exp;
  for (Exp.Init(result, TopAbs_WIRE,  TopAbs_FACE);
       Exp.More(); Exp.Next()) nWire++;
  for (Exp.Init(result, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) nFace++;
  for (Exp.Init(result, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) nSheet++;
  for (Exp.Init(result, TopAbs_SOLID); Exp.More(); Exp.Next()) nSolid++;
  if (outLevel > 1)
    printf(" Info: result has %d Solids, %d Sheets, %d Faces and %d Wires\n",
           nSolid, nSheet, nFace, nWire);

  int nBody = nWire+nFace+nSheet+nSolid;
  if (nBody == 0) {
    result.Nullify();
    if (outLevel > 0)
      printf(" EGADS Warning: Nothing found in result (EG_generalBoolean)!\n");
    if (emap != NULL) EG_free(emap);
    if (fmap != NULL) EG_free(fmap);
    return EGADS_NODATA;
  }

  egadsModel *mshape  = new egadsModel;
  mshape->shape       = result;
  mshape->nbody       = nBody;
  mshape->bbox.filled = 0;
  mshape->bodies      = new egObject*[nBody];
  for (i = 0; i < nBody; i++) {
    stat = EG_makeObject(context, &mshape->bodies[i]);
    if (stat != EGADS_SUCCESS) {
      for (j = 0; j < i; j++) {
        egObject  *obj   = mshape->bodies[j];
        egadsBody *pbody = (egadsBody *) obj->blind;
        delete pbody;
        EG_deleteObject(mshape->bodies[j]);
      }
      delete [] mshape->bodies;
      delete mshape;
      if (emap != NULL) EG_free(emap);
      if (fmap != NULL) EG_free(fmap);
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
  for (Exp.Init(mshape->shape, TopAbs_WIRE,  TopAbs_FACE);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
  }
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
  }
  for (i = 0; i < nBody; i++) {
    egObject  *pobj  = mshape->bodies[i];
    egadsBody *pbody = (egadsBody *) pobj->blind;
    BRepCheck_Analyzer sCheck(pbody->shape);
    if (!sCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Error: Result %d/%d is inValid (EG_solidBoolean)!\n",
               i+1, nBody);
      for (j = 0; j < nBody; j++) {
        egObject  *obj   = mshape->bodies[j];
        egadsBody *pbody = (egadsBody *) obj->blind;
        delete pbody;
        EG_deleteObject(mshape->bodies[j]);
      }
      delete [] mshape->bodies;
      delete mshape;
      if (emap != NULL) EG_free(emap);
      if (fmap != NULL) EG_free(fmap);
      return EGADS_GEOMERR;
    }
  }

  stat = EG_makeObject(context, &omodel);
  if (stat != EGADS_SUCCESS) {
    result.Nullify();
    for (i = 0; i < nBody; i++) {
      egObject  *obj   = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) obj->blind;
      delete pbody;
      EG_deleteObject(mshape->bodies[i]);
    }
    delete [] mshape->bodies;
    delete mshape;
    if (emap != NULL) EG_free(emap);
    if (fmap != NULL) EG_free(fmap);
    return stat;
  }
  omodel->oclass = MODEL;
  omodel->blind  = mshape;
  TopTools_IndexedMapOfShape rmape, rmap;
  TopExp::MapShapes(result, TopAbs_EDGE, rmape);
  TopExp::MapShapes(result, TopAbs_FACE, rmap);

  // set up the bodies
  for (i = 0; i < nBody; i++) {
    egObject  *pobj  = mshape->bodies[i];
    egadsBody *pbody = (egadsBody *) pobj->blind;
    pobj->topObj     = omodel;
    stat = EG_traverseBody(context, i, pobj, omodel, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      mshape->nbody = i;
      EG_destroyTopology(omodel);
      if (emap != NULL) EG_free(emap);
      if (fmap != NULL) EG_free(fmap);
      return stat;
    }
    // copy the attributes
    EG_attriBodyDup(src,  pobj);
    EG_attriBodyDup(tool, pobj);
    if (emap != NULL)
      for (j = 0; j < pbody->edges.map.Extent(); j++) {
        TopoDS_Edge dsedge = TopoDS::Edge(pbody->edges.map(j+1));
        index = rmape.FindIndex(dsedge);
        if (index == 0) continue;
        if (emap[index-1] == NULL) continue;
        EG_attributeDup(emap[index-1], pbody->edges.objs[j]);
      }
    if (fmap != NULL)
      for (j = 0; j < pbody->faces.map.Extent(); j++) {
        TopoDS_Face dsface = TopoDS::Face(pbody->faces.map(j+1));
        index = rmap.FindIndex(dsface);
        if (index == 0) continue;
        if (fmap[index-1] == NULL) continue;
        EG_attributeDup(fmap[index-1], pbody->faces.objs[j]);
      }
  }
  if (emap != NULL) EG_free(emap);
  if (fmap != NULL) EG_free(fmap);
  EG_referenceObject(omodel, context);

  *model = omodel;
  return EGADS_SUCCESS;
#endif
}


static int
EG_modelBoolean(const egObject *src, const egObject *tool, int oper,
                      egObject **model)
{
  int             i, j, k, stat, outLevel, iface, index, nerr, fullAttrs;
  int             *emap = NULL, *fmap = NULL;
  egObject        *context, *omodel;
  TopoDS_Shape    result;
  const egObject  *face = NULL;
  TopExp_Explorer Exp;

  if (src->blind == NULL) return EGADS_NODATA;
  outLevel  = EG_outLevel(src);
  context   = EG_context(src);
  fullAttrs = EG_fullAttrs(src);

  if ((oper != INTERSECTION) && (oper != FUSION) && (oper != SUBTRACTION)) {
    if (outLevel > 0)
      printf(" EGADS Error: BAD Operator = %d (EG_solidBoolean)!\n",
             oper);
    return EGADS_RANGERR;
  }
  if (tool == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Tool (EG_solidBoolean)!\n");
    return EGADS_NULLOBJ;
  }
  if (tool->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Tool is not an EGO (EG_solidBoolean)!\n");
    return EGADS_NOTOBJ;
  }
  if (tool->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Tool has no data (EG_solidBoolean)!\n");
    return EGADS_NODATA;
  }
  if (EG_sameThread(tool)) {
    if (outLevel > 0)
      printf(" EGADS Error: Bad Thread on Tool (EG_solidBoolean)!\n");
    return EGADS_CNTXTHRD;
  }
  if (EG_context(tool) != context) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_solidBoolean)!\n");
    return EGADS_MIXCNTX;
  }
  egadsModel      *pmdl = (egadsModel *) src->blind;
  TopoDS_Compound ssrc  = TopoDS::Compound(pmdl->shape);

  if (oper == FUSION) {

    if  ((tool->oclass != FACE) &&
        ((tool->oclass != BODY) || (tool->mtype != FACEBODY))) {
      printf(" EGADS Error: Face Tool is wrong type (EG_solidBoolean)!\n");
      return EGADS_NOTTOPO;
    }
    if (tool->oclass == FACE) {
      face = tool;
    } else {
      egadsBody *pbodf = (egadsBody *) tool->blind;
      face = pbodf->faces.objs[0];
    }
    if (face == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Face Tool (EG_solidBoolean)!\n");
      return EGADS_NULLOBJ;
    }
    if (face->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Tool is not an EGO (EG_solidBoolean)!\n");
      return EGADS_NOTOBJ;
    }
    if (face->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Tool has no data (EG_solidBoolean)!\n");
      return EGADS_NODATA;
    }
    egadsFace    *pface = (egadsFace *) face->blind;
    TopoDS_Shape stool  = pface->face;
    try {
      BRepAlgoAPI_Fuse BSO(ssrc, stool);
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Fusion (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      stat   = 0;
      TopExp_Explorer Expe;
      for (Expe.Init(result, TopAbs_EDGE); Expe.More(); Expe.Next()) {
        TopoDS_Vertex V1, V2;
        TopoDS_Edge   Edge = TopoDS::Edge(Expe.Current());
        if (BRep_Tool::Degenerated(Edge)) continue;
        TopExp::Vertices(Edge, V2, V1, Standard_True);
        if (V2.IsNull() && V1.IsNull()) stat++;
      }
      if (stat != 0) {
        // extend the tool Face and try again
        Handle(Geom_Surface) hSurf = BRep_Tool::Surface(pface->face);
        BRepLib_MakeFace MFace(hSurf, Standard_True);
        TopoDS_Face eFace = MFace.Face();
        // get the intersection Edge(s)
        BRepAlgoAPI_Section Sec(ssrc, eFace, Standard_False);
        Sec.ComputePCurveOn1(Standard_True);
        Sec.Approximation(Standard_True);
        Sec.Build();
        if (!Sec.IsDone()) {
          if (outLevel > 0)
            printf(" EGADS Error: Can't Section (EG_solidBoolean)!\n");
          return EGADS_GEOMERR;
        }
        TopoDS_Shape scribe = Sec.Shape();
        // scribe the complete intersection
        BRepFeat_SplitShape Split(ssrc);
        for (Expe.Init(scribe, TopAbs_EDGE); Expe.More(); Expe.Next()) {
          TopoDS_Edge Edge = TopoDS::Edge(Expe.Current());
          TopoDS_Face Face;
          if (!Sec.HasAncestorFaceOn1(Edge, Face)) continue;
          Split.Add(Edge, Face);
        }
        Split.Build();
        if (!Split.IsDone()) {
          if (outLevel > 0)
            printf(" EGADS Error: Can't Split (EG_solidBoolean)!\n");
          return EGADS_GEOMERR;
        }
        TopoDS_Shape newShape = Split.Shape();

        // map the Edges & Faces for future union attribution
        TopTools_IndexedMapOfShape splmape, smape, splmap, smap;
        TopExp::MapShapes(newShape, TopAbs_EDGE, splmape);
        TopExp::MapShapes(ssrc,     TopAbs_EDGE, smape);
        TopExp::MapShapes(newShape, TopAbs_FACE, splmap);
        TopExp::MapShapes(ssrc,     TopAbs_FACE, smap);
        int *spltabe = NULL, *spltab = NULL;

        if (fullAttrs != 0) {
          if (splmape.Extent() != 0) {
            spltabe = (int *) EG_alloc(splmape.Extent()*sizeof(int));
            if (spltabe == NULL) {
              if (outLevel > 0)
                printf(" EGADS Error: Malloc Edge SpliTable (EG_solidBoolean)!\n");
              return EGADS_MALLOC;
            }
            for (i = 0; i < splmape.Extent(); i++) spltabe[i] = 0;
            for (i = 1; i <= smape.Extent(); i++) {
              TopoDS_Edge dsedge = TopoDS::Edge(smape(i));
              const TopTools_ListOfShape& listEdges = Split.Modified(dsedge);
              if (listEdges.Extent() > 0) {
                /* modified Edges */
                TopTools_ListIteratorOfListOfShape it(listEdges);
                for (; it.More(); it.Next()) {
                  TopoDS_Edge genedge = TopoDS::Edge(it.Value());
                  index = splmape.FindIndex(genedge);
                  if (index <= 0) continue;
                  spltabe[index-1] = i;
                }
              } else {
                index = splmape.FindIndex(dsedge);
                if (index <= 0) continue;
                spltabe[index-1] = i;
              }
            }
          }
        }

        if (splmap.Extent() != 0) {
          spltab = (int *) EG_alloc(splmap.Extent()*sizeof(int));
          if (spltab == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc of Split Table (EG_solidBoolean)!\n");
            if (spltabe != NULL) EG_free(spltabe);
            return EGADS_MALLOC;
          }
          for (i = 0; i < splmap.Extent(); i++) spltab[i] = 0;
          for (i = 1; i <= smap.Extent(); i++) {
            TopoDS_Face dsface = TopoDS::Face(smap(i));
            const TopTools_ListOfShape& listFaces = Split.Modified(dsface);
            if (listFaces.Extent() > 0) {
              /* modified Faces */
              TopTools_ListIteratorOfListOfShape it(listFaces);
              for (; it.More(); it.Next()) {
                TopoDS_Face genface = TopoDS::Face(it.Value());
                index = splmap.FindIndex(genface);
                if (index <= 0) continue;
                spltab[index-1] = i;
              }
            } else {
              index = splmap.FindIndex(dsface);
              if (index <= 0) continue;
              spltab[index-1] = i;
            }
          }
          for (i = 0; i < splmap.Extent(); i++)
            if (spltab[i] == 0) {
              printf(" EGADS Error: Mapping Failed (EG_solidBoolean)!\n");
              if (spltabe != NULL) EG_free(spltabe);
              EG_free(spltab);
              return EGADS_GEOMERR;
            }
        }

        // redo the union
        BRepAlgoAPI_Fuse BSO(newShape, stool);
        if (!BSO.IsDone()) {
          printf(" EGADS Error: Can't do SBO Fusion (EG_solidBoolean)!\n");
          if (spltabe != NULL) EG_free(spltabe);
          if (spltab  != NULL) EG_free(spltab);
          return EGADS_GEOMERR;
        }
        result = BSO.Shape();
        iface  = 1;
        EG_matchMdlFace(BSO, newShape, iface, stool, result, fullAttrs,
                        &emap, &fmap);
#if CASVER >= 710
        TopTools_IndexedMapOfShape rmap, rmape;
        TopExp::MapShapes(result, TopAbs_EDGE, rmape);
        TopExp::MapShapes(result, TopAbs_FACE, rmap);
        if (fmap != NULL) {
          for (i = 0; i < rmap.Extent(); i++)
            if (fmap[i] == 0) {
              TopoDS_Face rface = TopoDS::Face(rmap(i+1));
              if (rface.IsSame(stool)) fmap[i] = -1;
            }
        }
        BRep_Builder aBB;
        TopoDS_Iterator aIt;
#if CASVER >= 730
        TopTools_ListOfShape aLCB;
        TopTools_ListIteratorOfListOfShape aItLCB;
#else
        BOPCol_ListOfShape aLCB;
        BOPCol_ListIteratorOfListOfShape aItLCB;
#endif

        // find all connected Faces and put them in Shells
        BOPTools_AlgoTools::MakeConnexityBlocks(result, TopAbs_EDGE, TopAbs_FACE,
                                                aLCB);
        TopoDS_Shape aRC;
        BOPTools_AlgoTools::MakeContainer(TopAbs_COMPOUND, aRC);
        int nbdy = 0;
        for (aItLCB.Initialize(aLCB); aItLCB.More(); aItLCB.Next()) {
          TopoDS_Shape aRCB;
          BOPTools_AlgoTools::MakeContainer(TopAbs_SHELL, aRCB);
          const TopoDS_Shape& aCB = aItLCB.Value();
          for (aIt.Initialize(aCB); aIt.More(); aIt.Next())
            aBB.Add(aRCB, aIt.Value());
          nbdy++;
          aBB.Add(aRC, aRCB);
        }
        result = aRC;

        // remap the Edges
        if (emap != NULL) {
          TopTools_IndexedMapOfShape smape;
          TopExp::MapShapes(result, TopAbs_EDGE, smape);
          if (rmape.Extent() != smape.Extent()) {
            printf(" EGADS Warning: Fuse/Shell # Edge Mismatch %d %d!\n",
                   rmape.Extent(), smape.Extent());
          } else {
            int *tmap = (int *) EG_alloc(rmape.Extent()*sizeof(int));
            if (tmap == NULL) {
              printf(" EGADS Warning: Fuse/Shell Edge mapping allocation %d!\n",
                     rmape.Extent());
            } else {
              for (j = 0; j < rmape.Extent(); j++) tmap[j] = 0;
              for (i = 0; i < rmape.Extent(); i++) {
                TopoDS_Edge redge = TopoDS::Edge(rmape(i+1));
                j = smape.FindIndex(redge) - 1;
//              printf(" %d/%d: Edge maps to %d\n", i+1, rmape.Extent(), j+1);
                if (j >= 0) tmap[j] = emap[i];
              }
              EG_free(emap);
              emap = tmap;
            }
          }
        }

        // remap the Faces
        if (fmap != NULL) {
          TopTools_IndexedMapOfShape smap;
          TopExp::MapShapes(result, TopAbs_FACE, smap);
          if (rmap.Extent() != smap.Extent()) {
            printf(" EGADS Warning: Fuse/Shell # Face Mismatch %d %d!\n",
                   rmap.Extent(), smap.Extent());
          } else {
            int *tmap = (int *) EG_alloc(rmap.Extent()*sizeof(int));
            if (tmap == NULL) {
              printf(" EGADS Warning: Fuse/Shell mapping allocation %d!\n",
                     rmap.Extent());
            } else {
              for (j = 0; j < rmap.Extent(); j++) tmap[j] = 0;
              for (i = 0; i < rmap.Extent(); i++) {
                TopoDS_Face rface = TopoDS::Face(rmap(i+1));
                j = smap.FindIndex(rface) - 1;
//              printf(" %d/%d:  maps to %d\n", i+1, rmap.Extent(), j+1);
                if (j >= 0) tmap[j] = fmap[i];
              }
              EG_free(fmap);
              fmap = tmap;
            }
          }
        }
#endif
        // patch up the Edge map
        if (spltabe != NULL) {
          TopTools_IndexedMapOfShape rmape;
          TopExp::MapShapes(result, TopAbs_EDGE, rmape);
          for (i = 0; i < rmape.Extent(); i++)
            if (emap[i] > 0)
              if (spltabe[emap[i]-1] == 0) {
                if (outLevel > 1)
                  printf(" EGADS Warning: Mapping maybe lost for Edge %d!\n",
                         i+1);
                emap[i] = 0;
              } else {
                emap[i] = spltabe[emap[i]-1];
                continue;
              }
          EG_free(spltabe);
        }
        // patch up the Face map
        if (spltab != NULL) {
          TopTools_IndexedMapOfShape rmap;
          TopExp::MapShapes(result, TopAbs_FACE, rmap);
          for (i = 0; i < rmap.Extent(); i++)
            if (fmap[i] > 0) {
              fmap[i] = spltab[fmap[i]-1];
              continue;
            }
          EG_free(spltab);
        }
      } else {
        iface = 1;
        EG_matchMdlFace(BSO, ssrc, iface, stool, result, fullAttrs,
                        &emap, &fmap);
#if CASVER >= 710
        TopTools_IndexedMapOfShape rmap, rmape;
        TopExp::MapShapes(result, TopAbs_EDGE, rmape);
        TopExp::MapShapes(result, TopAbs_FACE, rmap);
        if (fmap != NULL) {
          for (i = 0; i < rmap.Extent(); i++)
            if (fmap[i] == 0) {
              TopoDS_Face rface = TopoDS::Face(rmap(i+1));
              if (rface.IsSame(stool)) fmap[i] = -1;
            }
        }
        BRep_Builder aBB;
        TopoDS_Iterator aIt;
#if CASVER >= 730
        TopTools_ListOfShape aLCB;
        TopTools_ListIteratorOfListOfShape aItLCB;
#else
        BOPCol_ListOfShape aLCB;
        BOPCol_ListIteratorOfListOfShape aItLCB;
#endif

        // find all connected Faces and put them in Shells
        BOPTools_AlgoTools::MakeConnexityBlocks(result, TopAbs_EDGE, TopAbs_FACE,
                                                aLCB);
        TopoDS_Shape aRC;
        BOPTools_AlgoTools::MakeContainer(TopAbs_COMPOUND, aRC);
        int nbdy = 0;
        for (aItLCB.Initialize(aLCB); aItLCB.More(); aItLCB.Next()) {
          TopoDS_Shape aRCB;
          BOPTools_AlgoTools::MakeContainer(TopAbs_SHELL, aRCB);
          const TopoDS_Shape& aCB = aItLCB.Value();
          for (aIt.Initialize(aCB); aIt.More(); aIt.Next())
            aBB.Add(aRCB, aIt.Value());
          nbdy++;
          aBB.Add(aRC, aRCB);
        }
        result = aRC;

        // remap the Edges
        if (emap != NULL) {
          TopTools_IndexedMapOfShape smape;
          TopExp::MapShapes(result, TopAbs_EDGE, smape);
          if (rmape.Extent() != smape.Extent()) {
            printf(" EGADS Warning: Fuse/Shell # Edge Mismatch %d %d!\n",
                   rmape.Extent(), smape.Extent());
          } else {
            int *tmap = (int *) EG_alloc(rmape.Extent()*sizeof(int));
            if (tmap == NULL) {
              printf(" EGADS Warning: Fuse/Shell Edge mapping allocation %d!\n",
                     rmape.Extent());
            } else {
              for (j = 0; j < rmape.Extent(); j++) tmap[j] = 0;
              for (i = 0; i < rmape.Extent(); i++) {
                TopoDS_Edge redge = TopoDS::Edge(rmape(i+1));
                j = smape.FindIndex(redge) - 1;
//              printf(" %d/%d: Edge maps to %d\n", i+1, rmape.Extent(), j+1);
                if (j >= 0) tmap[j] = emap[i];
              }
              EG_free(emap);
              emap = tmap;
            }
          }
        }

        // remap the Faces
        if (fmap != NULL) {
          TopTools_IndexedMapOfShape smap;
          TopExp::MapShapes(result, TopAbs_FACE, smap);
          if (rmap.Extent() != smap.Extent()) {
            printf(" EGADS Warning: Fuse/Shell # Face Mismatch %d %d!\n",
                   rmap.Extent(), smap.Extent());
          } else {
            int *tmap = (int *) EG_alloc(rmap.Extent()*sizeof(int));
            if (tmap == NULL) {
              printf(" EGADS Warning: Fuse/Shell mapping allocation %d!\n",
                     rmap.Extent());
            } else {
              for (j = 0; j < rmap.Extent(); j++) tmap[j] = 0;
              for (i = 0; i < rmap.Extent(); i++) {
                TopoDS_Face rface = TopoDS::Face(rmap(i+1));
                j = smap.FindIndex(rface) - 1;
//              printf(" %d/%d:  maps to %d\n", i+1, rmap.Extent(), j+1);
                if (j >= 0) tmap[j] = fmap[i];
              }
              EG_free(fmap);
              fmap = tmap;
            }
          }
        }
#endif
      }
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Error: SBO Fusion Exception (EG_solidBoolean)!\n");
      printf("              %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Error: SBO Fusion Exception (EG_solidBoolean)!\n");
      return EGADS_GEOMERR;
    }

  } else if (oper == SUBTRACTION) {

    if (tool->oclass != BODY) {
      if (outLevel > 0)
        printf(" EGADS Error: Tool is not a Body (EG_solidBoolean)!\n");
      return EGADS_NOTBODY;
    }
    if (tool->mtype != SOLIDBODY) {
      if (outLevel > 0)
        printf(" EGADS Error: Tool is not a Solid Body (EG_solidBoolean)!\n");
      return EGADS_NOTTOPO;
    }
    egadsBody    *pbods = (egadsBody *) tool->blind;
    TopoDS_Solid stool  = TopoDS::Solid(pbods->shape);
    try {
      BRepAlgoAPI_Cut BSO(ssrc, stool);
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Subtraction (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      iface  = 0;
      EG_matchMdlFace(BSO, ssrc, iface, stool, result, fullAttrs, &emap, &fmap);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Error: SBO Subtraction Exception (EG_solidBoolean)!\n");
      printf("              %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Error: SBO Subtraction Exception (EG_solidBoolean)!\n");
      return EGADS_GEOMERR;
    }

  } else {

    if (tool->oclass != BODY) {
      if (outLevel > 0)
        printf(" EGADS Error: Tool is not a Body (EG_solidBoolean)!\n");
      return EGADS_NOTBODY;
    }
    if (tool->mtype != SOLIDBODY) {
      if (outLevel > 0)
        printf(" EGADS Error: Tool is not a Solid Body (EG_solidBoolean)!\n");
      return EGADS_NOTTOPO;
    }
    egadsBody    *pbods = (egadsBody *) tool->blind;
    TopoDS_Solid stool  = TopoDS::Solid(pbods->shape);
    try {
      BRepAlgoAPI_Common BSO(ssrc, stool);
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Intersection (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      iface  = 0;
      EG_matchMdlFace(BSO, ssrc, iface, stool, result, fullAttrs, &emap, &fmap);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Error: SBO Inters Exception (EG_solidBoolean)!\n");
      printf("              %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Error: SBO Inters Exception (EG_solidBoolean)!\n");
      return EGADS_GEOMERR;
    }

  }

  int nWire  = 0;
  int nFace  = 0;
  int nSheet = 0;
  int nSolid = 0;
  for (Exp.Init(result, TopAbs_WIRE,  TopAbs_FACE);
       Exp.More(); Exp.Next()) nWire++;
  for (Exp.Init(result, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) nFace++;
  for (Exp.Init(result, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) nSheet++;
  for (Exp.Init(result, TopAbs_SOLID); Exp.More(); Exp.Next()) nSolid++;
  if (outLevel > 1)
    printf(" Info: result has %d Solids, %d Sheets, %d Faces and %d Wires\n",
           nSolid, nSheet, nFace, nWire);

  int nBody = nWire+nFace+nSheet+nSolid;
  if (nBody == 0) {
    result.Nullify();
    if (outLevel > 0)
      printf(" EGADS Warning: Nothing found in result (EG_solidBoolean)!\n");
    return EGADS_NODATA;
  }

  egadsModel *mshape  = new egadsModel;
  mshape->shape       = result;
  mshape->nbody       = nBody;
  mshape->bbox.filled = 0;
  mshape->bodies      = new egObject*[nBody];
  for (i = 0; i < nBody; i++) {
    stat = EG_makeObject(context, &mshape->bodies[i]);
    if (stat != EGADS_SUCCESS) {
      for (j = 0; j < i; j++) {
        egObject  *obj   = mshape->bodies[j];
        egadsBody *pbody = (egadsBody *) obj->blind;
        delete pbody;
        EG_deleteObject(mshape->bodies[j]);
      }
      delete [] mshape->bodies;
      delete mshape;
      if (emap != NULL) EG_free(emap);
      if (fmap != NULL) EG_free(fmap);
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
  for (Exp.Init(mshape->shape, TopAbs_WIRE,  TopAbs_FACE);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
  }
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
  }
  for (i = 0; i < nBody; i++) {
    egObject  *pobj  = mshape->bodies[i];
    egadsBody *pbody = (egadsBody *) pobj->blind;
    BRepCheck_Analyzer sCheck(pbody->shape);
    if (!sCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Error: Result %d/%d is inValid (EG_solidBoolean)!\n",
               i+1, nBody);
      for (j = 0; j < nBody; j++) {
        egObject  *obj   = mshape->bodies[j];
        egadsBody *pbody = (egadsBody *) obj->blind;
        delete pbody;
        EG_deleteObject(mshape->bodies[j]);
      }
      delete [] mshape->bodies;
      delete mshape;
      if (emap != NULL) EG_free(emap);
      if (fmap != NULL) EG_free(fmap);
      return EGADS_GEOMERR;
    }
  }

  stat = EG_makeObject(context, &omodel);
  if (stat != EGADS_SUCCESS) {
    result.Nullify();
    for (i = 0; i < nBody; i++) {
      egObject  *obj   = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) obj->blind;
      delete pbody;
      EG_deleteObject(mshape->bodies[i]);
    }
    delete [] mshape->bodies;
    delete mshape;
    if (emap != NULL) EG_free(emap);
    if (fmap != NULL) EG_free(fmap);
    return stat;
  }
  omodel->oclass = MODEL;
  omodel->blind  = mshape;
  TopTools_IndexedMapOfShape smap, rmap, smape, rmape;
  TopExp::MapShapes(ssrc,   TopAbs_EDGE, smape);
  TopExp::MapShapes(result, TopAbs_EDGE, rmape);
  TopExp::MapShapes(ssrc,   TopAbs_FACE, smap);
  TopExp::MapShapes(result, TopAbs_FACE, rmap);

  for (i = 0; i < nBody; i++) {
    egObject  *pobj  = mshape->bodies[i];
    egadsBody *pbody = (egadsBody *) pobj->blind;
    pobj->topObj     = omodel;
    stat = EG_traverseBody(context, i, pobj, omodel, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      mshape->nbody = i;
      EG_destroyTopology(omodel);
      if (emap != NULL) EG_free(emap);
      if (fmap != NULL) EG_free(fmap);
      return stat;
    }
    EG_attriBodyDup(src, pobj);
    if (iface == 0) EG_attriBodyDup(tool, pobj);
    if (emap != NULL) {
      // fill in the attributes from cut Edges
      for (j = 0; j < pbody->edges.map.Extent(); j++) {
        TopoDS_Edge dsedge = TopoDS::Edge(pbody->edges.map(j+1));
        index = rmape.FindIndex(dsedge);
        if (index == 0) continue;
        index = emap[index-1];
        if (index == 0) continue;
        if (outLevel > 2)
          printf(" %d:  Edge mapping[%d] = %d\n", i, j, index);
        if (index > 0) {
          for (k = 0; k < pmdl->nbody; k++) {
            egObject  *bsrc  = pmdl->bodies[k];
            egadsBody *pbods = (egadsBody *) bsrc->blind;
            int ind = pbods->edges.map.FindIndex(smape(index));
            if (ind == 0) continue;
            EG_attributeDup(pbods->edges.objs[ind-1], pbody->edges.objs[j]);
            break;
          }
        } else {
          if (iface == 0) {
            /* only from Body tools */
            egadsBody *pbodt = (egadsBody *) tool->blind;
            EG_attributeDup(pbodt->edges.objs[-index-1], pbody->edges.objs[j]);
          }
        }
      }
    }
    if (fmap != NULL) {
      // fill in the attributes from cut Faces
      for (j = 0; j < pbody->faces.map.Extent(); j++) {
        TopoDS_Face dsface = TopoDS::Face(pbody->faces.map(j+1));
        index = rmap.FindIndex(dsface);
        if (index == 0) continue;
        index = fmap[index-1];
        if (index == 0) continue;
        if (outLevel > 2)
          printf(" %d:  Face mapping[%d] = %d\n", i, j, index);
        if (index > 0) {
          for (k = 0; k < pmdl->nbody; k++) {
            egObject  *bsrc  = pmdl->bodies[k];
            egadsBody *pbods = (egadsBody *) bsrc->blind;
            int ind = pbods->faces.map.FindIndex(smap(index));
            if (ind == 0) continue;
            EG_attributeDup(pbods->faces.objs[ind-1], pbody->faces.objs[j]);
            break;
          }
        } else {
          if (iface == 0) {
            egadsBody *pbodt = (egadsBody *) tool->blind;
            EG_attributeDup(pbodt->faces.objs[-index-1], pbody->faces.objs[j]);
          } else {
            EG_attributeDup(face, pbody->faces.objs[j]);
          }
        }
      }
    }
  }
  if (emap != NULL) EG_free(emap);
  if (fmap != NULL) EG_free(fmap);
  EG_referenceObject(omodel, context);

  *model = omodel;
  return EGADS_SUCCESS;
}


static int
EG_unionMatch(const egObject *src, const egObject *tool, egObject **model)
{
  int    i, j, k, nsf, ntf, ne, *fmap, *emap;
  double sBox[6], tBox[6];

  int outLevel     = EG_outLevel(src);
  egadsBody *psrc  = (egadsBody *) src->blind;
  egadsBody *ptool = (egadsBody *) tool->blind;
  nsf  = psrc->faces.map.Extent();
  ntf  = ptool->faces.map.Extent();
  k    = ne = psrc->edges.map.Extent();
  if (ntf > ne) k = ntf;
  fmap = (int *) EG_alloc(nsf*sizeof(int));
  emap = (int *) EG_alloc(  k*sizeof(int));
  if ((fmap == NULL) || (emap == NULL)) {
    if (emap != NULL) EG_free(emap);
    if (fmap != NULL) EG_free(fmap);
    return EGADS_MALLOC;
  }

  // find matching Edges by BBox, centroid & arclength
  for (i = 0; i < ne; i++) {
    emap[i] = 0;
    TopoDS_Edge sedge = TopoDS::Edge(psrc->edges.map(i+1));
    if (BRep_Tool::Degenerated(sedge)) {
      emap[i] = -1;
      continue;
    }
    double stol = BRep_Tool::Tolerance(sedge);
    Bnd_Box sbox;
    BRepBndLib::Add(sedge, sbox);
    sbox.Get(sBox[0], sBox[1], sBox[2], sBox[3], sBox[4], sBox[5]);
    BRepGProp    BProps;
    GProp_GProps sProps;
    BProps.LinearProperties(sedge, sProps);
    gp_Pnt sCG = sProps.CentreOfMass();
    for (j = 0; j < ptool->edges.map.Extent(); j++) {
      TopoDS_Edge tedge = TopoDS::Edge(ptool->edges.map(j+1));
      if (BRep_Tool::Degenerated(tedge)) continue;
      double ttol = BRep_Tool::Tolerance(tedge);
      if (stol > ttol) ttol = stol;
      Bnd_Box tbox;
      BRepBndLib::Add(tedge, tbox);
      tbox.Get(tBox[0], tBox[1], tBox[2], tBox[3], tBox[4], tBox[5]);
      if ((fabs(sBox[0]-tBox[0]) <= ttol) && (fabs(sBox[1]-tBox[1]) <= ttol) &&
          (fabs(sBox[2]-tBox[2]) <= ttol) && (fabs(sBox[3]-tBox[3]) <= ttol) &&
          (fabs(sBox[4]-tBox[4]) <= ttol) && (fabs(sBox[5]-tBox[5]) <= ttol)) {
        GProp_GProps tProps;
        BProps.LinearProperties(tedge, tProps);
        gp_Pnt tCG = tProps.CentreOfMass();
        if ((fabs(sCG.X() - tCG.X()) <= ttol) &&
            (fabs(sCG.Y() - tCG.Y()) <= ttol) &&
            (fabs(sCG.Z() - tCG.Z()) <= ttol)) {
          double toler = tProps.Mass();
          if (toler < 2.0) toler = 2.0;
          if (fabs(sProps.Mass()-tProps.Mass()) <= toler*ttol) {
            if (outLevel > 2)
              printf(" EGADS Info: Edges pass = %d/%d  %lf/%lf\n",
                     i, j, sProps.Mass(), tProps.Mass());
            emap[i] = j+1;
            break;
          }
        }
      }
    }
  }

  // find matching Faces by Edges, BBox, centroid & surface area
  for (i = 0; i < nsf; i++) {
    fmap[i] = 0;
    TopoDS_Face sface = TopoDS::Face(psrc->faces.map(i+1));
    // look at the Edges
    TopTools_IndexedMapOfShape smap;
    TopExp::MapShapes(sface, TopAbs_EDGE, smap);
    for (k = j = 1; j <= smap.Extent(); j++) {
      int index = psrc->edges.map.FindIndex(smap(j));
      if (index == 0) {
        k = 0;
        break;
      }
      if (emap[index-1] == 0) {
        k = 0;
        break;
      }
    }
    if (k == 0) continue;
    // ckeck BBox
    double stol = BRep_Tool::Tolerance(sface);
    Bnd_Box sbox;
    BRepBndLib::Add(sface, sbox);
    sbox.Get(sBox[0], sBox[1], sBox[2], sBox[3], sBox[4], sBox[5]);
    BRepGProp    BProps;
    GProp_GProps sProps;
    BProps.SurfaceProperties(sface, sProps);
    gp_Pnt sCG = sProps.CentreOfMass();
    for (j = 0; j < ptool->faces.map.Extent(); j++) {
      TopoDS_Face tface = TopoDS::Face(ptool->faces.map(j+1));
      double ttol = BRep_Tool::Tolerance(tface);
      if (stol > ttol) ttol = stol;
      Bnd_Box tbox;
      BRepBndLib::Add(tface, tbox);
      tbox.Get(tBox[0], tBox[1], tBox[2], tBox[3], tBox[4], tBox[5]);
      if ((fabs(sBox[0]-tBox[0]) <= ttol) && (fabs(sBox[1]-tBox[1]) <= ttol) &&
          (fabs(sBox[2]-tBox[2]) <= ttol) && (fabs(sBox[3]-tBox[3]) <= ttol) &&
          (fabs(sBox[4]-tBox[4]) <= ttol) && (fabs(sBox[5]-tBox[5]) <= ttol)) {
        GProp_GProps tProps;
        BProps.SurfaceProperties(tface, tProps);
        gp_Pnt tCG = tProps.CentreOfMass();
        // center of mass
        if ((fabs(sCG.X() - tCG.X()) <= ttol) &&
            (fabs(sCG.Y() - tCG.Y()) <= ttol) &&
            (fabs(sCG.Z() - tCG.Z()) <= ttol)) {
          // now surface area
          double toler = tProps.Mass()*ttol;
          if (toler < 4.0) toler = 4.0;
          if (fabs(sProps.Mass()-tProps.Mass()) <= toler*ttol) {
            if (outLevel > 2)
              printf(" EGADS Info: Faces pass = %d/%d  %lf/%lf\n",
                     i, j, sProps.Mass(), tProps.Mass());
            fmap[i] = j+1;
            break;
          }
        }
      }
    }
  }

  // any matched Faces?
  for (i = 0; i < ntf; i++) emap[i] = 0;
  for (j = i = 0; i < nsf; i++) {
    if (fmap[i] == 0) continue;
    emap[fmap[i]-1] = i+1;
    j++;
  }
  if (j == 0) {
    EG_free(emap);
    EG_free(fmap);
    return EGADS_EMPTY;
  }

  // take the faces that don't match up and sew them together
  k = nsf + ntf - 2*j;
  const egObject **faces = (const egObject **) EG_alloc(k*sizeof(egObject *));
  if (faces == NULL) {
    EG_free(emap);
    EG_free(fmap);
    return EGADS_MALLOC;
  }
  for (k = i = 0; i < nsf; i++) {
    if (fmap[i] != 0) continue;
    faces[k] = psrc->faces.objs[i];
    k++;
  }
  for (i = 0; i < ntf; i++) {
    if (emap[i] != 0) continue;
    faces[k] = ptool->faces.objs[i];
    k++;
  }
  EG_free(emap);
  EG_free(fmap);

  int stat = EG_sewFaces(k, faces, 0.0, 0, model);
  if (outLevel > 1)
    printf(" EGADS Info: EG_sewFaces = %d (EG_unionMatch)!\n", stat);
  EG_free(faces);

  return stat;
}


int
EG_fuseSheets(const egObject *srcx, const egObject *toolx, egObject **sheet)
{
  int            j, index, stat, outLevel, nShell, nerr, oclass, mtype, eei = 0;
  int            *senses, **emap = NULL, **fmap = NULL;
  egObject       *context, *obj, *newSrc, *model, *newModel, **bodies;
  TopoDS_Shape   result;
  const egObject *src, *tool;
/*
  double         tol, toler;
*/

  *sheet = NULL;
  src    = srcx;
  tool   = toolx;
  if  (src == NULL)                 return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if  (src->oclass != BODY)         return EGADS_NOTBODY;
  if ((src->mtype != SHEETBODY) &&
      (src->mtype != FACEBODY))     return EGADS_NOTTOPO;
  if  (src->blind == NULL)          return EGADS_NODATA;
  if  (EG_sameThread(src))          return EGADS_CNTXTHRD;
  if  (tool == NULL)                return EGADS_NULLOBJ;
  if  (tool->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if  (tool->oclass != BODY)        return EGADS_NOTBODY;
  if ((tool->mtype != SHEETBODY) &&
      (tool->mtype != FACEBODY))    return EGADS_NOTTOPO;
  if  (tool->blind == NULL)         return EGADS_NODATA;
  if  (EG_sameThread(tool))         return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);
  if  (context != EG_context(tool)) return EGADS_MIXCNTX;

  egadsBody *pbods = (egadsBody *) src->blind;
  egadsBody *pbodt = (egadsBody *) tool->blind;
  if (pbods->faces.map.Extent() == 1) {
    /* swap source and tool so tool will be the single Face */
    src   = toolx;
    tool  = srcx;
    pbods = (egadsBody *) src->blind;
    pbodt = (egadsBody *) tool->blind;
  }
  int nft = pbodt->faces.map.Extent();
/*
  stat = EG_tolerance(src, &toler);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_tolerance src = %d (EG_fuseSheets)!\n", stat);
    return stat;
  }
  stat = EG_tolerance(tool, &tol);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_tolerance tool = %d (EG_fuseSheets)!\n", stat);
    return stat;
  }
  if (tol > toler) toler = tol;
*/

  try {

#if CASVER >= 710
    TopTools_ListOfShape sList, tList;
    sList.Append(pbods->shape);
    tList.Append(pbodt->shape);

    BRepAlgoAPI_Fuse BSO;
    BSO.SetArguments(sList);
    BSO.SetTools(tList);
//  BSO.SetFuzzyValue(toler);
    BSO.SetNonDestructive(Standard_True);
    BSO.SetUseOBB(Standard_True);
    BSO.Build();
    TopoDS_Shape resultx = BSO.Shape();

    BRep_Builder    aBB;
    TopoDS_Iterator aIt;
#if CASVER >= 730
    TopTools_ListOfShape aLCB;
    TopTools_ListIteratorOfListOfShape aItLCB;
#else
    BOPCol_ListOfShape aLCB;
    BOPCol_ListIteratorOfListOfShape aItLCB;
#endif

    // find all connected Faces and put them in Shells
    BOPTools_AlgoTools::MakeConnexityBlocks(resultx, TopAbs_EDGE, TopAbs_FACE,
                                            aLCB);
    nShell = 0;
    for (aItLCB.Initialize(aLCB); aItLCB.More(); aItLCB.Next()) {
      TopoDS_Shape aRC;
      BOPTools_AlgoTools::MakeContainer(TopAbs_SHELL, aRC);
      const TopoDS_Shape& aCB = aItLCB.Value();
      for (aIt.Initialize(aCB); aIt.More(); aIt.Next())
        aBB.Add(aRC, aIt.Value());
      result = aRC;
      nShell++;
    }
#else
    /* perform the fuse operation */
    BRepAlgoAPI_Fuse BSO(pbods->shape, pbodt->shape);
    if (!BSO.IsDone()) {
      printf(" EGADS Error: Can't do SBO Fusion (EG_fuseSheets)!\n");
      return EGADS_GEOMERR;
    }
    TopoDS_Shape resultx = BSO.Shape();

    nShell = 0;
    TopExp_Explorer Exp;
    for (Exp.Init(resultx, TopAbs_SHELL); Exp.More(); Exp.Next()) {
      result = Exp.Current();
      nShell++;
    }
#endif
    if (nShell != 1) {
      if (outLevel > 0)
        printf(" EGADS Warning: Fuse has %d Fragments (EG_fuseSheets)!\n",
               nShell);
      return EGADS_GEOMERR;
    }

    /* did we fail but with a single Face as the tool?
       if so try special scribe option found in modelBoolean */
    if (nft == 1) {
      TopTools_IndexedMapOfShape emap;
      TopExp::MapShapes(result, TopAbs_EDGE, emap);
      for (eei = j = 0; j < emap.Extent(); j++)
        if ((emap(j+1).Orientation() == TopAbs_INTERNAL) ||
            (emap(j+1).Orientation() == TopAbs_EXTERNAL)) eei++;

      if (eei != 0) {
        stat = EG_copyObject(src, NULL, &newSrc);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Cannot copyObject %d (EG_fuseSheets)!\n",
                   stat);
          return stat;
        }
        stat = EG_makeTopology(context, NULL, MODEL, 0, NULL, 1, &newSrc, NULL,
                               &model);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Cannot make Model %d (EG_fuseSheets)!\n",
                   stat);
          return stat;
        }
        if (tool->mtype == SHEETBODY) {
          stat = EG_modelBoolean(model, pbodt->faces.objs[0], FUSION, &newModel);
        } else {
          stat = EG_modelBoolean(model, tool, FUSION, &newModel);
        }
        EG_deleteObject(model);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_modelBoolean = %d (EG_fuseSheets)!\n",
                   stat);
          return stat;
        }
        stat = EG_getTopology(newModel, &obj, &oclass, &mtype, NULL,
                              &j, &bodies, &senses);
        if (j != 1) {
          if (outLevel > 0)
            printf(" EGADS Warning: Fuse has %d Bodies (EG_fuseSheets)!\n", j);
          EG_deleteObject(newModel);
          return EGADS_GEOMERR;
        }
        if (bodies[0]->mtype != SHEETBODY) {
          if (outLevel > 0)
            printf(" EGADS Warning: Result is not a Sheet (EG_fuseSheets)!\n");
          EG_deleteObject(newModel);
          return EGADS_GEOMERR;
        }
        stat = EG_copyObject(bodies[0], NULL, &obj);
        EG_deleteObject(newModel);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Warning: Copy result = %d (EG_fuseSheets)!\n", stat);
          return EGADS_GEOMERR;
        }

        *sheet = obj;
        return EGADS_SUCCESS;
      }
    }

    BRepCheck_Analyzer rCheck(result);
    if (!rCheck.IsValid()) {
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(result);
#if CASVER >= 700
      sfs->FixShellTool()->SetNonManifoldFlag(Standard_True);
#endif
      sfs->Perform();
      TopoDS_Shape fixedShape = sfs->Shape();
      if (fixedShape.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Error: Resultant Body is invalid (EG_fuseSheets)!\n");
        return EGADS_GEOMERR;
      } else {
        BRepCheck_Analyzer sfCheck(fixedShape);
        if (!sfCheck.IsValid()) {
          if (outLevel > 0)
            printf(" EGADS Error: Fixed Body is invalid (EG_fuseSheets)!\n");
          return EGADS_GEOMERR;
        } else {
          result = fixedShape;
        }
      }
    }
    EG_matchFaces(BSO, src, tool, TopAbs_SHELL, result, &emap, &fmap);
  }
  catch (...) {
    printf(" EGADS Error: Construction Exception (EG_fuseSheets)!\n");
    return EGADS_CONSTERR;
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_fuseSheets)!\n");
    return stat;
  }
  egadsBody *pbodr   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = SHEETBODY;
  pbodr->nodes.objs  = NULL;
  pbodr->edges.objs  = NULL;
  pbodr->loops.objs  = NULL;
  pbodr->faces.objs  = NULL;
  pbodr->shells.objs = NULL;
  pbodr->senses      = NULL;
  pbodr->shape       = result;
  pbodr->bbox.filled = 0;
  pbodr->massFill    = 0;
  obj->blind         = pbodr;
  stat = EG_traverseBody(context, 0, obj, obj, pbodr, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbodr;
    return stat;
  }
  EG_referenceObject(obj, context);

  EG_attriBodyDup(src,  obj);
  EG_attriBodyDup(tool, obj);
  if (emap != NULL) {
    // fill in the attributes from cut Edges
    for (j = 0; j < pbodr->edges.map.Extent(); j++) {
      index = emap[0][j];
      if (index == 0) continue;
      if (outLevel > 2)
        printf(" %d:  Edge mapping[%d] = %d\n", 0, j, index);
      if (index > 0) {
        EG_attributeDup(pbods->edges.objs[ index-1], pbodr->edges.objs[j]);
      } else {
        EG_attributeDup(pbodt->edges.objs[-index-1], pbodr->edges.objs[j]);
      }
    }
    EG_free(emap[0]);
    EG_free(emap);
  }
  if (fmap != NULL) {
    // fill in the attributes from cut Faces
    for (j = 0; j < pbodr->faces.map.Extent(); j++) {
      index = fmap[0][j];
      if (index == 0) continue;
      if (outLevel > 2)
        printf(" %d:  Face mapping[%d] = %d\n", 0, j, index);
      if (index > 0) {
        EG_attributeDup(pbods->faces.objs[ index-1], pbodr->faces.objs[j]);
      } else {
        EG_attributeDup(pbodt->faces.objs[-index-1], pbodr->faces.objs[j]);
      }
    }
    EG_free(fmap[0]);
    EG_free(fmap);
  }

  *sheet = obj;
  return EGADS_SUCCESS;
}


int
EG_solidBoolean(const egObject *src, const egObject *tool, int oper,
                      egObject **model)
{
  int             i, j, outLevel, index, stat, hit, nerr, rev = 0;
  int             **emap = NULL, **fmap = NULL;
  egObject        *context, *omodel;
  const egObject  *face  = NULL;
  egadsBody       *pbodt = NULL;
  TopoDS_Shape    stool, result;
  TopExp_Explorer Exp;

  *model = NULL;
  if (src == NULL)               return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(src))        return EGADS_CNTXTHRD;
  if (src->oclass == MODEL)      return EG_modelBoolean(src, tool, oper, model);
  if (src->oclass != BODY)       return EGADS_NOTBODY;
  if (src->mtype != SOLIDBODY)   return EGADS_NOTTOPO;
  if (src->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);

  if ((oper != SUBTRACTION) && (oper != INTERSECTION) &&
      (oper != FUSION)) {
    if (outLevel > 0)
      printf(" EGADS Error: BAD Operator = %d (EG_solidBoolean)!\n",
             oper);
    return EGADS_RANGERR;
  }
  if (tool == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Tool (EG_solidBoolean)!\n");
    return EGADS_NULLOBJ;
  }
  if (tool->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Tool is not an EGO (EG_solidBoolean)!\n");
    return EGADS_NOTOBJ;
  }
  if (tool->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Tool has no data (EG_solidBoolean)!\n");
    return EGADS_NODATA;
  }
  if (EG_sameThread(tool)) {
    if (outLevel > 0)
      printf(" EGADS Error: Bad Thread on Tool (EG_solidBoolean)!\n");
    return EGADS_CNTXTHRD;
  }
  if (EG_context(tool) != context) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_solidBoolean)!\n");
    return EGADS_MIXCNTX;
  }
  if ((oper == SUBTRACTION) && ((tool->oclass == FACE) ||
      ((tool->oclass == BODY) && (tool->mtype == FACEBODY)))) {
    if (tool->oclass == FACE) {
      face = tool;
    } else {
      egadsBody *pbodf = (egadsBody *) tool->blind;
      face = pbodf->faces.objs[0];
    }
    if (face == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Face Tool (EG_solidBoolean)!\n");
      return EGADS_NULLOBJ;
    }
    if (face->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Tool is not an EGO (EG_solidBoolean)!\n");
      return EGADS_NOTOBJ;
    }
    if (face->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Tool has no data (EG_solidBoolean)!\n");
      return EGADS_NODATA;
    }
  } else {
    if (tool->oclass != BODY) {
      if (outLevel > 0)
        printf(" EGADS Error: Tool is not a Body (EG_solidBoolean)!\n");
      return EGADS_NOTBODY;
    }
  }

  int checkAgain    = 1;
  egadsBody *pbods  = (egadsBody *) src->blind;
  TopoDS_Solid ssrc = TopoDS::Solid(pbods->shape);
  if (face == NULL) {
    pbodt = (egadsBody *) tool->blind;
    stool = pbodt->shape;
  } else {
    egadsFace *pface = (egadsFace *) face->blind;
    stool = pface->face;
  }

  if (oper == INTERSECTION) {

    try {
      BRepAlgoAPI_Common BSO(ssrc, stool);
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Intersection (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      checkAgain = 0;
      for (Exp.Init(result, TopAbs_SOLID); Exp.More(); Exp.Next()) {
        TopoDS_Solid solid = TopoDS::Solid(Exp.Current());
        BRepCheck_Analyzer sCheck(solid);
        if (!sCheck.IsValid()) {
          checkAgain = 1;
          if (outLevel > 1)
            printf(" EGADS Info: Inters Failed -- try reverse (EG_solidBoolean)!\n");
          // intersection is associtive -- see if it works the other way
          rev = 1;
          BRepAlgoAPI_Common BSO(stool, ssrc);
          if (!BSO.IsDone()) {
            printf(" EGADS Error: Can't do SBO Inters (EG_solidBoolean)!\n");
            return EGADS_GEOMERR;
          }
          result = BSO.Shape();
          EG_matchFaces(BSO, src, tool, TopAbs_SOLID, result, &emap, &fmap);
          break;
        }
      }
      if (rev == 0) EG_matchFaces(BSO, src, tool, TopAbs_SOLID, result,
                                  &emap, &fmap);
    }
    catch (...) {
      if (rev == 1) {
        printf(" EGADS Error: SBO Inters Exception (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
      rev = -1;
    }

    if (rev == -1) {
      rev = 1;
      if (outLevel > 1)
        printf(" EGADS Info: Inters Aborted -- try reverse (EG_solidBoolean)!\n");
      try {
        BRepAlgoAPI_Common BSO(stool, ssrc);
        if (!BSO.IsDone()) {
          printf(" EGADS Error: Can't do SBO Inters (EG_solidBoolean)!\n");
          return EGADS_GEOMERR;
        }
        result = BSO.Shape();
        EG_matchFaces(BSO, src, tool, TopAbs_SOLID, result, &emap, &fmap);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Error: SBO Inters Exception (EG_solidBoolean)!\n");
        printf("              %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Error: SBO Inters Exception (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
    }

  } else if (oper == SUBTRACTION) {

    try {
      BRepAlgoAPI_Cut BSO(ssrc, stool);
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Subtraction (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      EG_matchFaces(BSO, src, tool, TopAbs_SOLID, result, &emap, &fmap);
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Error: SBO Subraction Exception (EG_solidBoolean)!\n");
      printf("              %s\n", e.GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...) {
      printf(" EGADS Error: SBO Subraction Exception (EG_solidBoolean)!\n");
      return EGADS_GEOMERR;
    }

  } else {

    try {
      BRepAlgoAPI_Fuse BSO(ssrc, stool);
      if (!BSO.IsDone()) {
        printf(" EGADS Error: Can't do SBO Fusion (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
      result = BSO.Shape();
      i = 0;
      for (Exp.Init(result, TopAbs_SOLID); Exp.More(); Exp.Next()) i++;
      if (i == 2) {
        // try explicit matching
        stat = EG_unionMatch(src, tool, model);
        if (stat == EGADS_SUCCESS) return stat;
      }
      if ((i != 1) && (outLevel > 0))
        printf(" EGADS Warning: Union has %d Bodies (EG_solidBoolean)!\n", i);
      checkAgain = 0;
      for (Exp.Init(result, TopAbs_SOLID); Exp.More(); Exp.Next()) {
        TopoDS_Solid solid = TopoDS::Solid(Exp.Current());
        BRepCheck_Analyzer sCheck(solid);
        if (!sCheck.IsValid()) {
          checkAgain = 1;
          if (outLevel > 1)
            printf(" EGADS Info: Union Failed -- try reverse!\n");
          // union is associtive -- see if it works the other way
          rev = 1;
          BRepAlgoAPI_Fuse BSO(stool, ssrc);
          if (!BSO.IsDone()) {
            printf(" EGADS Error: Can't do SBO Union (EG_solidBoolean)!\n");
            return EGADS_GEOMERR;
          }
          result = BSO.Shape();
          EG_matchFaces(BSO, src, tool, TopAbs_SOLID, result, &emap, &fmap);
          break;
        }
      }
      if (rev == 0) EG_matchFaces(BSO, src, tool, TopAbs_SOLID, result,
                                  &emap, &fmap);
    }
    catch (...) {
      if (rev == 1) {
        printf(" EGADS Error: SBO Fusion Exception (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
      rev = -1;
    }

    if (rev == -1) {
      rev = 1;
      if (outLevel > 1)
        printf(" EGADS Info: Union Aborted -- try reverse!\n");
      try {
        BRepAlgoAPI_Fuse BSO(stool, ssrc);
        if (!BSO.IsDone()) {
          printf(" EGADS Error: Can't do SBO Fusion (EG_solidBoolean)!\n");
          return EGADS_GEOMERR;
        }
        result = BSO.Shape();
        EG_matchFaces(BSO, src, tool, TopAbs_SOLID, result, &emap, &fmap);
      }
      catch (const Standard_Failure& e) {
        printf(" EGADS Error: SBO Fusion Exception (EG_solidBoolean)!\n");
        printf("              %s\n", e.GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...) {
        printf(" EGADS Error: SBO Fusion Exception (EG_solidBoolean)!\n");
        return EGADS_GEOMERR;
      }
    }

  }

  int nBody = 0;
  for (Exp.Init(result, TopAbs_SOLID); Exp.More(); Exp.Next()) nBody++;
  if ((oper == FUSION) && (nBody == 2)) {
    stat = EG_unionMatch(src, tool, model);
    if (emap != NULL) {
      for (j = 0; j < nBody; j++) EG_free(emap[j]);
      EG_free(emap);
    }
    if (fmap != NULL) {
      for (j = 0; j < nBody; j++) EG_free(fmap[j]);
      EG_free(fmap);
    }
    return stat;
  }
  if (outLevel > 1)
    printf("   Boolean Solid Oper result has #%d solids!\n", nBody);
  if (nBody == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL SBO Result (EG_solidBoolean)!\n");
    return EGADS_NOTFOUND;
  }

  // vaild solids?
  if (checkAgain == 1) {
    i = 0;
    for (Exp.Init(result, TopAbs_SOLID); Exp.More(); Exp.Next()) {
      TopoDS_Solid solid = TopoDS::Solid(Exp.Current());
      BRepCheck_Analyzer sCheck(solid);
      if (!sCheck.IsValid()) {
        if (outLevel > 0)
          printf(" EGADS Warning: Solid %d/%d is invalid (EG_solidBoolean)!\n",
                 i+1, nBody);
        if (emap != NULL) {
          for (j = 0; j < nBody; j++) EG_free(emap[j]);
          EG_free(emap);
        }
        if (fmap != NULL) {
          for (j = 0; j < nBody; j++) EG_free(fmap[j]);
          EG_free(fmap);
        }
        return EGADS_CONSTERR;
      }
      i++;
    }
  }

  egadsModel *mshape  = new egadsModel;
  mshape->shape       = result;
  mshape->nbody       = nBody;
  mshape->bbox.filled = 0;
  mshape->bodies      = new egObject*[nBody];
  for (i = 0; i < nBody; i++) {
    stat = EG_makeObject(context, &mshape->bodies[i]);
    if (stat != EGADS_SUCCESS) {
      for (j = 0; j < i; j++) {
        egObject  *obj   = mshape->bodies[j];
        egadsBody *pbody = (egadsBody *) obj->blind;
        delete pbody;
        EG_deleteObject(mshape->bodies[j]);
      }
      delete [] mshape->bodies;
      delete mshape;
      if (emap != NULL) {
        for (j = 0; j < nBody; j++) EG_free(emap[j]);
        EG_free(emap);
      }
      if (fmap != NULL) {
        for (j = 0; j < nBody; j++) EG_free(fmap[j]);
        EG_free(fmap);
      }
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
  for (Exp.Init(mshape->shape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
  }

  stat = EG_makeObject(context, &omodel);
  if (stat != EGADS_SUCCESS) {
    result.Nullify();
    for (i = 0; i < nBody; i++) {
      egObject  *obj   = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) obj->blind;
      delete pbody;
      EG_deleteObject(mshape->bodies[i]);
    }
    delete [] mshape->bodies;
    delete mshape;
    if (emap != NULL) {
      for (j = 0; j < nBody; j++) EG_free(emap[j]);
      EG_free(emap);
    }
    if (fmap != NULL) {
      for (j = 0; j < nBody; j++) EG_free(fmap[j]);
      EG_free(fmap);
    }
    return stat;
  }
  omodel->oclass = MODEL;
  omodel->blind  = mshape;
  EG_referenceObject(omodel, context);

  for (hit = i = 0; i < nBody; i++) {
    egObject  *pobj  = mshape->bodies[i];
    egadsBody *pbody = (egadsBody *) pobj->blind;
    pobj->topObj     = omodel;
    stat = EG_traverseBody(context, i, pobj, omodel, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      EG_deleteObject(omodel);
      if (emap != NULL) {
        for (j = 0; j < nBody; j++) EG_free(emap[j]);
        EG_free(emap);
      }
      if (fmap != NULL) {
        for (j = 0; j < nBody; j++) EG_free(fmap[j]);
        EG_free(fmap);
      }
      return stat;
    }
    hit += nerr;

    EG_attriBodyDup(src,  pobj);
    if (face == NULL) EG_attriBodyDup(tool, pobj);
    if (emap != NULL) {
      // fill in the attributes from cut Edges
      for (j = 0; j < pbody->edges.map.Extent(); j++) {
        index = emap[i][j];
        if (index == 0) continue;
        if (outLevel > 2)
          printf(" %d:  Edge mapping[%d] = %d\n", i, j, index);
        if (index > 0) {
          EG_attributeDup(pbods->edges.objs[index-1], pbody->edges.objs[j]);
        } else {
          if (face == NULL)
            EG_attributeDup(pbodt->edges.objs[-index-1], pbody->edges.objs[j]);
        }
      }
      EG_free(emap[i]);
    }
    if (fmap != NULL) {
      // fill in the attributes from cut Faces
      for (j = 0; j < pbody->faces.map.Extent(); j++) {
        index = fmap[i][j];
        if (index == 0) continue;
        if (outLevel > 2)
          printf(" %d:  Face mapping[%d] = %d\n", i, j, index);
        if (index > 0) {
          EG_attributeDup(pbods->faces.objs[index-1], pbody->faces.objs[j]);
        } else {
          if (face == NULL) {
            EG_attributeDup(pbodt->faces.objs[-index-1], pbody->faces.objs[j]);
          } else {
            EG_attributeDup(face, pbody->faces.objs[j]);
          }
        }
      }
      EG_free(fmap[i]);
    }
  }
  if (emap != NULL) EG_free(emap);
  if (fmap != NULL) EG_free(fmap);

  if (hit != 0) {
    EG_deleteObject(omodel);
    return EGADS_TOPOERR;
  }

  *model = omodel;
  return EGADS_SUCCESS;
}


int
EG_intersection(const egObject *src, const egObject *tool, int *nEdge,
                /*@null@*/ egObject ***facEdg, egObject **model)
{
  int             i, j, n, stat, outLevel, nloop, sense, index, found, nerr;
  int             plane = 1;
  egObject        *context = NULL, *omodel = NULL, *geom = NULL, **list = NULL;
  const egObject  *face  = NULL;
  egadsFace       *pface = NULL;
  TopoDS_Shape    s1, result;
  TopoDS_Vertex   V1, V2, Vs, lV1, lV2;
  TopExp_Explorer Exp;

  *nEdge = 0;
  *model = NULL;
  if (facEdg != NULL) *facEdg = NULL;
  if  (src == NULL)               return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (src->oclass != BODY)       return EGADS_NOTBODY;
  if ((src->mtype != SOLIDBODY) && (src->mtype != SHEETBODY) &&
      (src->mtype != FACEBODY))   return EGADS_NOTTOPO;
  if  (src->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(src))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);

  if (tool == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Tool (EG_intersection)!\n");
    return EGADS_NULLOBJ;
  }
  if (tool->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Tool is not an EGO (EG_intersection)!\n");
    return EGADS_NOTOBJ;
  }
  if (tool->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Tool has no data (EG_intersection)!\n");
    return EGADS_NODATA;
  }
  if (EG_sameThread(tool)) {
    if (outLevel > 0)
      printf(" EGADS Error: Bad Thread on Tool (EG_intersection)!\n");
    return EGADS_CNTXTHRD;
  }
  if (EG_context(tool) != context) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_intersection)!\n");
    return EGADS_MIXCNTX;
  }

  if (tool->oclass == BODY) {
    egadsBody *pbodf = (egadsBody *) tool->blind;
    if (tool->mtype == FACEBODY) {
      face = pbodf->faces.objs[0];
    } else if (tool->mtype == SHEETBODY) {
      if (pbodf->faces.map.Extent() == 1) face = pbodf->faces.objs[0];
    } else if (tool->mtype == WIREBODY) {
      if (outLevel > 0)
        printf(" EGADS Error: Tool is a Wire Body (EG_intersection)!\n");
      return EGADS_NOTTOPO;
    }
  } else {
    if (tool->oclass != FACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Tool is not a Face or Body (EG_intersection)!\n");
      return EGADS_NOTBODY;
    }
    face = tool;
  }
  if (face != NULL) {
    /* a single Face in tool */
    pface = (egadsFace *) face->blind;
    geom  = pface->surface;
    if (geom->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Tool Surface is NULL (EG_intersection)!\n");
      return EGADS_NOTGEOM;
    }
    if (geom->mtype != PLANE) plane = 0;
    s1 = pface->face;
  } else {
    /* multiple Faces in tool */
    egadsBody *pbodf = (egadsBody *) tool->blind;
    s1 = pbodf->shape;
  }

  egadsBody *pbody = (egadsBody *) src->blind;
  BRepAlgoAPI_Section Sec(s1, pbody->shape, Standard_False);
  Sec.ComputePCurveOn1(Standard_True);
  Sec.ComputePCurveOn2(Standard_True);
  Sec.Approximation(Standard_True);
  Sec.Build();
  if (!Sec.IsDone()) {
    if (outLevel > 0)
      printf(" EGADS Error: Can't Section (EG_intersection)!\n");
    return EGADS_GEOMERR;
  }
  result = Sec.Shape();

  TopTools_IndexedMapOfShape MapE;
  TopExp::MapShapes(result, TopAbs_EDGE, MapE);
  int nedge = MapE.Extent();
  if (nedge == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: No Intersection (EG_intersection)!\n");
    return EGADS_CONSTERR;
  }

  // find the Loops
  loopInfo *info = new loopInfo[nedge];
  for (i = 0; i <  nedge; i++) info[i].lIndex = info[i].sense = 0;
  for (nloop = i = 1; i <= nedge; i++) {
    if (info[i-1].lIndex != 0) continue;
    TopoDS_Shape shape = MapE(i);
    TopoDS_Edge  Edge  = TopoDS::Edge(shape);
    TopExp::Vertices(Edge, V2, V1, Standard_True);
    sense = -1;
    if (Edge.Orientation() != TopAbs_REVERSED) {
      sense = 1;
      Vs    = V2;
      V2    = V1;
      V1    = Vs;
    }
    Vs    = V1;
    index = 0;
    info[i-1].lIndex = nloop;
    info[i-1].index  = index;
    info[i-1].sense  = sense;
    while (!Vs.IsSame(V2)) {
      for (j = 1; j <= nedge; j++) {
        if (info[j-1].lIndex != 0) continue;
        TopoDS_Edge lEdge = TopoDS::Edge(MapE(j));
        TopExp::Vertices(lEdge, lV1, lV2, Standard_True);
        if (V2.IsSame(lV1)) {
          index++;
          sense = 1;
          if (Edge.Orientation() == TopAbs_REVERSED) sense = -1;
          info[j-1].lIndex = nloop;
          info[j-1].index  = index;
          info[j-1].sense  = sense;
          V2 = lV2;
          break;
        } else if (V2.IsSame(lV2)) {
          index++;
          sense = -1;
          if (Edge.Orientation() == TopAbs_REVERSED) sense = 1;
          info[j-1].lIndex = nloop;
          info[j-1].index  = index;
          info[j-1].sense  = sense;
          V2 = lV1;
          break;
        }
      }
      if (j > nedge) {
        /* we are open -- check the other direction */
        TopExp::Vertices(Edge, V2, V1, Standard_True);
        if (Edge.Orientation() != TopAbs_FORWARD) {
          Vs = V2;
          V2 = V1;
          V1 = Vs;
        }
        j = 1;
        while (j <= nedge) {
          for (j = 1; j <= nedge; j++) {
            if (info[j-1].lIndex != 0) continue;
            TopoDS_Edge lEdge = TopoDS::Edge(MapE(j));
            TopExp::Vertices(lEdge, lV1, lV2, Standard_True);
            if (V2.IsSame(lV1)) {
              index++;
              sense = 1;
              if (Edge.Orientation() == TopAbs_FORWARD) sense = -1;
              info[j-1].lIndex = nloop;
              info[j-1].index  = index;
              info[j-1].sense  = sense;
              V2 = lV2;
              break;
            } else if (V2.IsSame(lV2)) {
              index++;
              sense = -1;
              if (Edge.Orientation() == TopAbs_FORWARD) sense = 1;
              info[j-1].lIndex = nloop;
              info[j-1].index  = index;
              info[j-1].sense  = sense;
              V2 = lV1;
              break;
            }
          }
        }
        break;
      }
    }
    nloop++;
  }
  nloop--;

  // create the EGADS objects for the WireBodies
  egObject **wireo = new egObject*[nloop];
  for (i = 0; i < nloop; i++) {
    stat = EG_makeObject(context, &wireo[i]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Loop object (EG_intersection)!\n");
      for (j = 0; j < i; j++) EG_deleteObject(wireo[j]);
      delete [] info;
      return stat;
    }
  }

  // make the OCC Wires and then the WireBodies
  TopoDS_Compound compound;
  BRep_Builder    builder3D;
  builder3D.MakeCompound(compound);
  for (i = 0; i < nloop; i++) {
    BRepBuilderAPI_MakeWire MW;
    index = 0;
    do {
      for (found = j = 0; j < nedge; j++)
        if (info[j].index == index)
          if (info[j].lIndex == i+1) {
            found = j+1;
            break;
          }
      if (found == 0) continue;
      TopoDS_Shape shape = MapE(found);
      TopoDS_Edge  Edge  = TopoDS::Edge(shape);
      if (Edge.Orientation() == TopAbs_REVERSED) {
        if (info[found-1].sense ==  1) Edge.Orientation(TopAbs_FORWARD);
      } else {
        if (info[found-1].sense == -1) Edge.Orientation(TopAbs_REVERSED);
      }
      MW.Add(Edge);
      if (MW.Error()) {
        if (outLevel > 0)
          printf(" EGADS Error: Problem with Edge %d (EG_intersection)!\n",
                 found);
        for (j = 0; j < nloop; j++) EG_deleteObject(wireo[j]);
        delete [] info;
        return EGADS_NODATA;
      }
      index++;
    } while (found != 0);
    if (!MW.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem with Loop (EG_intersection)!\n");
      for (j = 0; j < nloop; j++) EG_deleteObject(wireo[j]);
      delete [] info;
      return EGADS_NODATA;
    }
    TopoDS_Wire Wire = MW.Wire();
    builder3D.Add(compound, Wire);
    if (outLevel > 1)
      printf(" Wire %d made with %d edges!\n", i+1, index);

    egadsBody *pbodw   = new egadsBody;
    wireo[i]->oclass   = BODY;
    wireo[i]->mtype    = WIREBODY;
    pbodw->nodes.objs  = NULL;
    pbodw->edges.objs  = NULL;
    pbodw->loops.objs  = NULL;
    pbodw->faces.objs  = NULL;
    pbodw->shells.objs = NULL;
    pbodw->senses      = NULL;
    pbodw->shape       = Wire;
    pbodw->bbox.filled = 0;
    pbodw->massFill    = 0;
    wireo[i]->blind    = pbodw;
    stat = EG_traverseBody(context, i, wireo[i], wireo[i], pbodw, &nerr);
    if (stat != EGADS_SUCCESS) {
      delete pbodw;
      for (j = 0; j < nloop; j++) EG_deleteObject(wireo[j]);
      delete [] info;
      return stat;
    }
  }
  delete [] info;

  // fix up the WireBodies for PCurves on tool single Face

  if (plane == 0) {
    egadsSurface *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurface = psurf->handle;
    for (i = 0; i < nloop; i++) {
      egObject  *bobj  = wireo[i];
      egadsBody *pbodw = (egadsBody *) bobj->blind;
      egObject  *lobj  = pbodw->loops.objs[0];
      egadsLoop *ploop = (egadsLoop *) lobj->blind;
      egObject  *sobj  = NULL;
      egObject **edgeo = new egObject*[2*ploop->nedges];
      for (j = 0; j < ploop->nedges; j++) {
        edgeo[j] = ploop->edges[j];
        edgeo[j+ploop->nedges] = NULL;
      }
      delete [] ploop->edges;
      ploop->edges = edgeo;
      stat = EG_makeObject(context, &sobj);
      if (stat != EGADS_SUCCESS) continue;
      sobj->topObj = bobj;
      EG_completeSurf(sobj, hSurface);
      EG_referenceObject(sobj, lobj);
      ploop->surface = sobj;
      for (j = 0; j < ploop->nedges; j++) {
        egObject  *eobj  = ploop->edges[j];
        egadsEdge *pedge = (egadsEdge *) eobj->blind;
        TopoDS_Edge Edge = pedge->edge;
        Standard_Real f, l;
        Handle(Geom2d_Curve) hPCurv = BRep_Tool::CurveOnSurface(Edge, pface->face,
                                                                f, l);
        if (hPCurv.IsNull()) {
          Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, f, l);
          double toler = BRep_Tool::Tolerance(Edge);
          try {
            hPCurv = GeomProjLib::Curve2d(hCurve, f, l, hSurface, toler);
          }
          catch (const Standard_Failure& e) {
            continue;
          }
          catch (...) {
            continue;
          }
        }
        stat = EG_makeObject(context, &edgeo[j+ploop->nedges]);
        if (stat != EGADS_SUCCESS) continue;
        edgeo[j+ploop->nedges]->topObj = bobj;
        EG_completePCurve(edgeo[j+ploop->nedges],  hPCurv);
        EG_referenceObject(edgeo[j+ploop->nedges], lobj);
      }
    }
  }

  // Attach the tool Face Attributes to the new Edges

  for (i = 0; i < nloop; i++) {
    egObject  *bobj  = wireo[i];
    egadsBody *pbodw = (egadsBody *) bobj->blind;
    for (j = 0; j < pbodw->edges.map.Extent(); j++) {
      if (tool->oclass == BODY) {
        egadsBody *pbodf = (egadsBody *) tool->blind;
        TopoDS_Edge Edge = TopoDS::Edge(pbodw->edges.map(j+1));
        TopoDS_Face Face;
        if (Sec.HasAncestorFaceOn1(Edge, Face)) {
          index = pbodf->faces.map.FindIndex(Face);
          if (index <= 0) continue;
          EG_attributeDup(pbodf->faces.objs[index-1], pbodw->edges.objs[j]);
        }
      } else {
        EG_attributeDup(tool, pbodw->edges.objs[j]);
      }
    }
  }

  // fill in the Face/Edge pairs (if requested)
  if (facEdg != NULL) {
    *facEdg = list = (egObject **) EG_alloc(2*nedge*sizeof(egObject *));
    n = 0;
    if (list != NULL)
      for (i = 0; i < nloop; i++) {
        egObject  *bobj  = wireo[i];
        egadsBody *pbodw = (egadsBody *) bobj->blind;
        for (j = 0; j < pbodw->edges.map.Extent(); j++) {
          TopoDS_Edge Edge = TopoDS::Edge(pbodw->edges.map(j+1));
          TopoDS_Face Face;
          if (Sec.HasAncestorFaceOn2(Edge, Face)) {
            index = pbody->faces.map.FindIndex(Face);
            if (index <= 0) continue;
            list[n  ] = pbody->faces.objs[index-1];
            list[n+1] = pbodw->edges.objs[j];
            n += 2;
          }
        }
      }
    *nEdge = n/2;
  }

  // make the EGADS model

  egadsModel *mshape  = new egadsModel;
  mshape->shape       = compound;
  mshape->nbody       = nloop;
  mshape->bodies      = wireo;
  mshape->bbox.filled = 0;
  stat = EG_makeObject(context, &omodel);
  if (stat != EGADS_SUCCESS) {
    compound.Nullify();
    for (i = 0; i < nloop; i++) {
      egObject  *obj   = mshape->bodies[i];
      egadsBody *pbodw = (egadsBody *) obj->blind;
      delete pbodw;
      EG_deleteObject(mshape->bodies[i]);
    }
    delete [] mshape->bodies;
    delete mshape;
    return stat;
  }
  omodel->oclass = MODEL;
  omodel->blind  = mshape;
  for (i = 0; i < nloop; i++) {
    EG_referenceObject(mshape->bodies[i], omodel);
    EG_removeCntxtRef(mshape->bodies[i]);
  }
  EG_referenceObject(omodel, context);

  *model = omodel;
  return EGADS_SUCCESS;
}


int
EG_imprintBody(const egObject *src, int nedge, const egObject **facEdg,
                     egObject **result)
{
  int      i, outLevel, stat, nerr;
  egObject *context, *obj;

  *result = NULL;
  if  (src == NULL)               return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (src->oclass != BODY)       return EGADS_NOTBODY;
  if ((src->mtype != SOLIDBODY) && (src->mtype != SHEETBODY) &&
      (src->mtype != FACEBODY))   return EGADS_NOTTOPO;
  if  (src->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(src))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);

  if (nedge <= 0) {
    if (outLevel > 0)
      printf(" EGADS Warning: No Edges (EG_imprintBody)!\n");
    return EG_copyObject(src, NULL, result);
  }
  if (facEdg == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL facEdg Pointer (EG_imprintBody)!\n");
    return EGADS_NULLOBJ;
  }
  for (i = 0; i < 2*nedge; i++) {
    if (facEdg[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Object %d (EG_imprintBody)!\n", i/2+1);
      return EGADS_NULLOBJ;
    }
    if (facEdg[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not an EGO (EG_imprintBody)!\n",
               i/2+1);
      return EGADS_NOTOBJ;
    }
    if (facEdg[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d has no data (EG_imprintBody)!\n",
               i/2+1);
      return EGADS_NODATA;
    }
    if (EG_context(facEdg[i]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is wrong Context (EG_imprintBody)!\n",
               i/2+1);
      return EGADS_MIXCNTX;
    }
  }
  egadsBody *pbody = (egadsBody *) src->blind;
  for (i = 0; i < nedge; i++) {
    if (facEdg[2*i  ]->oclass != FACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not FACE (EG_imprintBody)!\n", i);
      return EGADS_NOTTOPO;
    }
    egadsFace *pface = (egadsFace *) facEdg[2*i]->blind;
    if (pbody->faces.map.FindIndex(pface->face) <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d is not in Body (EG_imprintBody)!\n", i);
      return EGADS_NOTBODY;
    }
    if ((facEdg[2*i+1]->oclass != EDGE) && (facEdg[2*i+1]->oclass != LOOP)) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not EDGE/LOOP (EG_imprintBody)!\n",
               i);
      return EGADS_NOTTOPO;
    }
    if (facEdg[2*i+1]->oclass == EDGE) {
      egadsEdge *pedge = (egadsEdge *) facEdg[2*i+1]->blind;
      if (pbody->edges.map.FindIndex(pedge->edge) > 0) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d is in Body (EG_imprintBody)!\n", i);
        return EGADS_NOTBODY;
      }
    }
  }

  /* try our own code */
  stat = EG_splitBody(src, nedge, facEdg, result);
  if ((stat == EGADS_SUCCESS) && (*result != NULL)) return EGADS_SUCCESS;

  /* use OpenCASCADE */
  if (outLevel > 0)
    printf(" EGADS Info: splitBody = %d\n", stat);
  if (result != NULL) {
    EG_deleteObject(*result);
    *result = NULL;
  }

  TopoDS_Shape newShape;
  BRepFeat_SplitShape Split(pbody->shape);
  try {
    for (i = 0; i < nedge; i++) {
      egadsFace *pface = (egadsFace *) facEdg[2*i  ]->blind;
      if (facEdg[2*i+1]->oclass == EDGE) {
        egadsEdge *pedge = (egadsEdge *) facEdg[2*i+1]->blind;
        Split.Add(pedge->edge, pface->face);
      } else {
        egadsLoop *ploop = (egadsLoop *) facEdg[2*i+1]->blind;
        Split.Add(ploop->loop, pface->face);
      }
    }
    Split.Build();
    if (!Split.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Can't Split (EG_imprintBody)!\n");
      return EGADS_GEOMERR;
    }
    newShape = Split.Shape();
    BRepCheck_Analyzer sCheck(newShape);
    if (!sCheck.IsValid()) {
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
      sfs->Perform();
      TopoDS_Shape fixedShape = sfs->Shape();
      if (fixedShape.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Error: Resultant Body is invalid (EG_imprintBody)!\n");
        return EGADS_GEOMERR;
      } else {
        BRepCheck_Analyzer sfCheck(fixedShape);
        if (!sfCheck.IsValid()) {
          printf(" EGADS Error: Fixed Body is invalid (EG_imprintBody)!\n");
          return EGADS_GEOMERR;
        } else {
          newShape = fixedShape;
        }
      }
    }
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Warning: Split Construction Error (EG_imprintBody)!\n");
    printf("                %s\n", e.GetMessageString());
    return EGADS_CONSTERR;
  }
  catch (...) {
    printf(" EGADS Warning: Split Construction Error (EG_imprintBody)!\n");
    return EGADS_CONSTERR;
  }
  TopTools_IndexedMapOfShape mapE;
  TopExp::MapShapes(newShape, TopAbs_EDGE, mapE);
  if (mapE.Extent() <= pbody->edges.map.Extent()) {
    if (outLevel > 0)
      printf(" EGADS Error: # new Edges = %d, # on src = %d (EG_imprintBody)!\n",
             mapE.Extent(), pbody->edges.map.Extent());
    return EGADS_CONSTERR;
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_imprintBody)!\n");
    return stat;
  }
  egadsBody *pbods   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = src->mtype;
  pbods->nodes.objs  = NULL;
  pbods->edges.objs  = NULL;
  pbods->loops.objs  = NULL;
  pbods->faces.objs  = NULL;
  pbods->shells.objs = NULL;
  pbods->senses      = NULL;
  pbods->shape       = newShape;
  pbods->bbox.filled = 0;
  pbods->massFill    = 0;
  obj->blind         = pbods;
  if (src->mtype == FACEBODY) {
    TopTools_IndexedMapOfShape mapF;
    TopExp::MapShapes(newShape, TopAbs_FACE, mapF);
    if (mapF.Extent() > 1) obj->mtype = SHEETBODY;
  }
  stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbods;
    return stat;
  }

  // map the Attributes
  EG_attriBodyDup(src, obj);
  for (i = 0; i < nedge; i++) {
    egadsFace *pface = (egadsFace *) facEdg[2*i]->blind;
    const TopTools_ListOfShape& listFaces = Split.Modified(pface->face);
    if (listFaces.Extent() > 0) {
      /* modified Faces */
      TopTools_ListIteratorOfListOfShape it(listFaces);
      for (; it.More(); it.Next()) {
        TopoDS_Face genface = TopoDS::Face(it.Value());
        int index = pbods->faces.map.FindIndex(genface);
        if (index <= 0) continue;
        EG_attributeDup(facEdg[2*i], pbods->faces.objs[index-1]);
      }
    }
  }
  // mapping does not always work here -- so
  EG_attrFixup(outLevel, obj, src);

  EG_referenceObject(obj, context);
  *result = obj;
  return EGADS_SUCCESS;
}


int
EG_filletBody(const egObject *src, int nedge, const egObject **edges,
              double radius, egObject **result, int **faceMap)
{
  int      i, k, outLevel, stat, nerr, fullAttrs, *fmap;
  egObject *context, *obj, *face, *edge, *node;

  *result = NULL;
  if (faceMap != NULL) *faceMap = NULL;
  if  (src == NULL)               return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (src->oclass != BODY)       return EGADS_NOTBODY;
  if ((src->mtype != SOLIDBODY) &&
      (src->mtype != SHEETBODY))  return EGADS_NOTTOPO;
  if  (src->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(src))        return EGADS_CNTXTHRD;
  outLevel  = EG_outLevel(src);
  context   = EG_context(src);
  fullAttrs = EG_fullAttrs(src);

  if (nedge <= 0) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edges (EG_filletBody)!\n");
    return EGADS_NODATA;
  }
  if (edges == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Edge Pointer (EG_filletBody)!\n");
    return EGADS_NULLOBJ;
  }
  egadsBody *pbody = (egadsBody *) src->blind;
  for (k = i = 0; i < nedge; i++) {
    if (edges[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Edge Object %d (EG_filletBody)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (edges[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge Object %d is not an EGO (EG_filletBody)!\n",
               i+1);
      return EGADS_NOTOBJ;
    }
    if (edges[i]->oclass != EDGE) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not EDGE (EG_filletBody)!\n", i+1);
      return EGADS_NOTTOPO;
    }
    if (edges[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge Object %d has no data (EG_filletBody)!\n",
               i+1);
      return EGADS_NODATA;
    }
    egadsEdge *pedge = (egadsEdge *) edges[i]->blind;
    if (pbody->edges.map.FindIndex(pedge->edge) == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge %d is NOT in Body (EG_filletBody)!\n",
               i+1);
      return EGADS_NOTBODY;
    }
    if (edges[i]->mtype != DEGENERATE) k++;
  }
  if (k == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: No nonDegenerate Edges (EG_filletBody)!\n");
    return EGADS_NODATA;
  }
  BRepCheck_Analyzer check(pbody->shape);
  if (!check.IsValid()) {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid Input Body (EG_filletBody)!\n");
    return EGADS_TOPOERR;
  }

  // fillet the body
  BRepFilletAPI_MakeFillet fillet(pbody->shape);
  for (i = 0; i < nedge; i++) {
    if (edges[i]->mtype == DEGENERATE) continue;
    egadsEdge *pedge = (egadsEdge *) edges[i]->blind;
    fillet.Add(radius, pedge->edge);
  }

  try {
    fillet.Build();
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Error: Fillet Exception (EG_filletBody)!\n");
    printf("              %s\n", e.GetMessageString());
    return EGADS_GEOMERR;
  }
  catch (...) {
    printf(" EGADS Error: Fillet Exception (EG_filletBody)!\n");
    return EGADS_GEOMERR;
  }

  if (!fillet.IsDone()) {
    if (outLevel > 0)
      printf(" EGADS Error: Can't Fillet (EG_filletBody)!\n");
    return EGADS_GEOMERR;
  }
  TopoDS_Shape newShape = fillet.Shape();
  BRepCheck_Analyzer fCheck(newShape);
  if (!fCheck.IsValid()) {
    Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
    sfs->Perform();
    TopoDS_Shape fixedShape = sfs->Shape();
    if (fixedShape.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Error: Filleted Body is invalid (EG_filletBody)!\n");
      return EGADS_GEOMERR;
    } else {
      BRepCheck_Analyzer sfCheck(fixedShape);
      if (!sfCheck.IsValid()) {
        printf(" EGADS Error: Fixed Body is invalid (EG_filletBody)!\n");
        return EGADS_GEOMERR;
      } else {
        newShape = fixedShape;
      }
    }
  }

  // make sure we have the correct result!
  if (newShape.ShapeType() == TopAbs_COMPOUND) {
    int nshell = 0, nsolid = 0;
    TopExp_Explorer Exp;
    for (Exp.Init(newShape, TopAbs_SHELL, TopAbs_SOLID);
         Exp.More(); Exp.Next()) nshell++;
    for (Exp.Init(newShape, TopAbs_SOLID); Exp.More(); Exp.Next()) nsolid++;
    if (nshell+nsolid != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: Number of Results = %d (EG_filletBody)!\n",
               nshell+nsolid);
      return EGADS_CONSTERR;
    }
    if (nshell == 1) {
      Exp.Init(newShape, TopAbs_SHELL, TopAbs_SOLID);
      newShape = Exp.Current();
      }
    if (nsolid == 1) {
      Exp.Init(newShape, TopAbs_SOLID);
      newShape = Exp.Current();
    }
  }
  if ((newShape.ShapeType() != TopAbs_SOLID) &&
      (newShape.ShapeType() != TopAbs_SHELL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Incorrect Result (EG_filletBody)!\n");
    return EGADS_CONSTERR;
  }
  if ((src->mtype == SOLIDBODY) &&
      (newShape.ShapeType() != TopAbs_SOLID)) {
    if (outLevel > 0)
      printf(" EGADS Error: Result Not a Solid (EG_filletBody)!\n");
    return EGADS_CONSTERR;
  }
  if ((src->mtype == SHEETBODY) &&
      (newShape.ShapeType() != TopAbs_SHELL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Result Not a Sheet (EG_filletBody)!\n");
    return EGADS_CONSTERR;
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_filletBody)!\n");
    return stat;
  }
  egadsBody *pbods   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = src->mtype;
  pbods->nodes.objs  = NULL;
  pbods->edges.objs  = NULL;
  pbods->loops.objs  = NULL;
  pbods->faces.objs  = NULL;
  pbods->shells.objs = NULL;
  pbods->senses      = NULL;
  pbods->shape       = newShape;
  pbods->bbox.filled = 0;
  pbods->massFill    = 0;
  obj->blind         = pbods;
  stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbods;
    return stat;
  }

  // map the Attributes
  EG_attriBodyDup(src, obj);
  if (fullAttrs != 0) {
    for (i = 0; i < pbody->edges.map.Extent(); i++) {
      edge = pbody->edges.objs[i];
      egadsEdge *pedge = (egadsEdge *) edge->blind;
      const TopTools_ListOfShape& listEdges = fillet.Modified(pedge->edge);
      if (listEdges.Extent() > 0) {
        /* modified Edges */
        TopTools_ListIteratorOfListOfShape it(listEdges);
        for (; it.More(); it.Next()) {
          TopoDS_Edge genedge = TopoDS::Edge(it.Value());
          int index = pbods->edges.map.FindIndex(genedge);
          if (index > 0) EG_attributeDup(edge, pbods->edges.objs[index-1]);
        }
      }
    }
  }
  
  for (i = 0; i < pbody->faces.map.Extent(); i++) {
    face = pbody->faces.objs[i];
    egadsFace *pface = (egadsFace *) face->blind;
    const TopTools_ListOfShape& listFaces = fillet.Modified(pface->face);
    if (listFaces.Extent() > 0) {
      /* modified Faces */
      TopTools_ListIteratorOfListOfShape it(listFaces);
      for (; it.More(); it.Next()) {
        TopoDS_Face genface = TopoDS::Face(it.Value());
        int index = pbods->faces.map.FindIndex(genface);
        if (index > 0) EG_attributeDup(face, pbods->faces.objs[index-1]);
      }
    }
  }

  // Face map info: source, modified, generated
  if (faceMap != NULL) {
    fmap = (int *) EG_alloc(2*pbods->faces.map.Extent()*sizeof(int));
    if (fmap != NULL) {
      for (i = 0; i < 2*pbods->faces.map.Extent(); i++) fmap[i] = 0;
      for (i = 0; i <   pbody->faces.map.Extent(); i++) {
        face = pbody->faces.objs[i];
        egadsFace *pface = (egadsFace *) face->blind;
        int index = pbods->faces.map.FindIndex(pface->face);
        if (index > 0) {
          fmap[2*index-2] = FACEDUP;
          fmap[2*index-1] = i+1;
        }
      }
      for (i = 0; i < pbody->faces.map.Extent(); i++) {
        face = pbody->faces.objs[i];
        egadsFace *pface = (egadsFace *) face->blind;
        const TopTools_ListOfShape& listFaces = fillet.Modified(pface->face);
        if (listFaces.Extent() > 0) {
          /* modified faces */
          TopTools_ListIteratorOfListOfShape it(listFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) {
              fmap[2*index-2] = FACECUT;
              fmap[2*index-1] = i+1;
            }
          }
        }
      }
      for (i = 0; i < pbody->edges.map.Extent(); i++) {
        edge = pbody->edges.objs[i];
        egadsEdge *pedge = (egadsEdge *) edge->blind;
        const TopTools_ListOfShape& genFaces = fillet.Generated(pedge->edge);
        if (genFaces.Extent() > 0) {
          /* generated faces from Edges */
          TopTools_ListIteratorOfListOfShape it(genFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) {
              fmap[2*index-2] = EDGEOFF;
              fmap[2*index-1] = i+1;
            }
          }
        }
      }
      for (i = 0; i < pbody->nodes.map.Extent(); i++) {
        node = pbody->nodes.objs[i];
        egadsNode *pnode = (egadsNode *) node->blind;
        const TopTools_ListOfShape& genFaces = fillet.Generated(pnode->node);
        if (genFaces.Extent() > 0) {
          /* generated faces from Nodes */
          TopTools_ListIteratorOfListOfShape it(genFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) {
              fmap[2*index-2] = NODEOFF;
              fmap[2*index-1] = i+1;
            }
          }
        }
      }
      *faceMap = fmap;
    }
  }

  EG_referenceObject(obj, context);
  *result = obj;
  return EGADS_SUCCESS;
}


int
EG_chamferBody(const egObject *src, int nedge, const egObject **edges,
               const egObject **faces, double dis1, double dis2,
               egObject **result, int **faceMap)
{
  int      i, k, outLevel, stat, nerr, fullAttrs, *fmap;
  egObject *context, *obj, *face, *edge, *node;

  *result = NULL;
  if (faceMap != NULL) *faceMap = NULL;
  if  (src == NULL)               return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (src->oclass != BODY)       return EGADS_NOTBODY;
  if ((src->mtype != SOLIDBODY) &&
      (src->mtype != SHEETBODY))  return EGADS_NOTTOPO;
  if  (src->blind == NULL)        return EGADS_NODATA;
  if  (EG_sameThread(src))        return EGADS_CNTXTHRD;
  outLevel  = EG_outLevel(src);
  context   = EG_context(src);
  fullAttrs = EG_fullAttrs(src);

  if (nedge <= 0) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edges (EG_chamferBody)!\n");
    return EGADS_NODATA;
  }
  if (edges == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Edge Pointer (EG_chamferBody)!\n");
    return EGADS_NULLOBJ;
  }
  if (faces == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face Pointer (EG_chamferBody)!\n");
    return EGADS_NULLOBJ;
  }
  egadsBody *pbody = (egadsBody *) src->blind;
  for (k = i = 0; i < nedge; i++) {

    if (edges[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Edge Object %d (EG_chamferBody)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (edges[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge Object %d is not an EGO (EG_chamferBody)!\n",
               i+1);
      return EGADS_NOTOBJ;
    }
    if (edges[i]->oclass != EDGE) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not EDGE (EG_chamferBody)!\n", i+1);
      return EGADS_NOTTOPO;
    }
    if (edges[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge Object %d has no data (EG_chamferBody)!\n",
               i+1);
      return EGADS_NODATA;
    }
    egadsEdge *pedge = (egadsEdge *) edges[i]->blind;
    if (pbody->edges.map.FindIndex(pedge->edge) == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge %d is NOT in Body (EG_chamferBody)!\n",
               i+1);
      return EGADS_NOTBODY;
    }
    if (edges[i]->mtype != DEGENERATE) k++;

    if (faces[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Face Object %d (EG_chamferBody)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (faces[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Object %d is not an EGO (EG_chamferBody)!\n",
               i+1);
      return EGADS_NOTOBJ;
    }
    if (faces[i]->oclass != FACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not FACE (EG_chamferBody)!\n", i+1);
      return EGADS_NOTTOPO;
    }
    if (faces[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Object %d has no data (EG_chamferBody)!\n",
               i+1);
      return EGADS_NODATA;
    }
    egadsFace *pface = (egadsFace *) faces[i]->blind;
    if (pbody->faces.map.FindIndex(pface->face) == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d is NOT in Body (EG_chamferBody)!\n",
               i+1);
      return EGADS_NOTBODY;
    }
  }
  if (k == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: No nonDegenerate Edges (EG_chamferBody)!\n");
    return EGADS_NODATA;
  }
  BRepCheck_Analyzer check(pbody->shape);
  if (!check.IsValid()) {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid Input Body (EG_chamferBody)!\n");
    return EGADS_TOPOERR;
  }

  // chamfer the body
  BRepFilletAPI_MakeChamfer chamfer(pbody->shape);
  for (i = 0; i < nedge; i++) {
    if (edges[i]->mtype == DEGENERATE) continue;
    egadsEdge *pedge = (egadsEdge *) edges[i]->blind;
    egadsFace *pface = (egadsFace *) faces[i]->blind;
    chamfer.Add(dis1, dis2, pedge->edge, pface->face);
  }

  try {
    chamfer.Build();
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Error: Chamfer Exception (EG_chamferBody)!\n");
    printf("              %s\n", e.GetMessageString());
    return EGADS_GEOMERR;
  }
  catch (...) {
    printf(" EGADS Error: Chamfer Exception (EG_chamferBody)!\n");
    return EGADS_GEOMERR;
  }

  if (!chamfer.IsDone()) {
    if (outLevel > 0)
      printf(" EGADS Error: Can't Chamfer (EG_chamferBody)!\n");
    return EGADS_GEOMERR;
  }
  TopoDS_Shape newShape = chamfer.Shape();
  BRepCheck_Analyzer fCheck(newShape);
  if (!fCheck.IsValid()) {
    Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
    sfs->Perform();
    TopoDS_Shape fixedShape = sfs->Shape();
    if (fixedShape.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Error: Chamfered Body is invalid (EG_chamferBody)!\n");
      return EGADS_GEOMERR;
    } else {
      BRepCheck_Analyzer sfCheck(fixedShape);
      if (!sfCheck.IsValid()) {
        printf(" EGADS Error: Fixed Body is invalid (EG_chamferBody)!\n");
        return EGADS_GEOMERR;
      } else {
        newShape = fixedShape;
      }
    }
  }

  // make sure we have the correct result!
  if (newShape.ShapeType() == TopAbs_COMPOUND) {
    int nshell = 0, nsolid = 0;
    TopExp_Explorer Exp;
    for (Exp.Init(newShape, TopAbs_SHELL, TopAbs_SOLID);
         Exp.More(); Exp.Next()) nshell++;
    for (Exp.Init(newShape, TopAbs_SOLID); Exp.More(); Exp.Next()) nsolid++;
    if (nshell+nsolid != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: Number of Results = %d (EG_chamferBody)!\n",
               nshell+nsolid);
      return EGADS_CONSTERR;
    }
    if (nshell == 1) {
      Exp.Init(newShape, TopAbs_SHELL, TopAbs_SOLID);
      newShape = Exp.Current();
      }
    if (nsolid == 1) {
      Exp.Init(newShape, TopAbs_SOLID);
      newShape = Exp.Current();
    }
  }
  if ((newShape.ShapeType() != TopAbs_SOLID) &&
      (newShape.ShapeType() != TopAbs_SHELL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Incorrect Result (EG_chamferBody)!\n");
    return EGADS_CONSTERR;
  }
  if ((src->mtype == SOLIDBODY) &&
      (newShape.ShapeType() != TopAbs_SOLID)) {
    if (outLevel > 0)
      printf(" EGADS Error: Result Not a Solid (EG_chamferBody)!\n");
    return EGADS_CONSTERR;
  }
  if ((src->mtype == SHEETBODY) &&
      (newShape.ShapeType() != TopAbs_SHELL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Result Not a Sheet (EG_chamferBody)!\n");
    return EGADS_CONSTERR;
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_chamferBody)!\n");
    return stat;
  }
  egadsBody *pbods   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = src->mtype;
  pbods->nodes.objs  = NULL;
  pbods->edges.objs  = NULL;
  pbods->loops.objs  = NULL;
  pbods->faces.objs  = NULL;
  pbods->shells.objs = NULL;
  pbods->senses      = NULL;
  pbods->shape       = newShape;
  pbods->bbox.filled = 0;
  pbods->massFill    = 0;
  obj->blind         = pbods;
  stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbods;
    return stat;
  }

  // map the Attributes
  EG_attriBodyDup(src, obj);
  if (fullAttrs != 0) {
    for (i = 0; i < pbody->edges.map.Extent(); i++) {
      edge = pbody->edges.objs[i];
      egadsEdge *pedge = (egadsEdge *) edge->blind;
      const TopTools_ListOfShape& listEdges = chamfer.Modified(pedge->edge);
      if (listEdges.Extent() > 0) {
        /* modified Edges */
        TopTools_ListIteratorOfListOfShape it(listEdges);
        for (; it.More(); it.Next()) {
          TopoDS_Edge genedge = TopoDS::Edge(it.Value());
          int index = pbods->edges.map.FindIndex(genedge);
          if (index > 0) EG_attributeDup(edge, pbods->edges.objs[index-1]);
        }
      }
    }
  }

  for (i = 0; i < pbody->faces.map.Extent(); i++) {
    face = pbody->faces.objs[i];
    egadsFace *pface = (egadsFace *) face->blind;
    const TopTools_ListOfShape& listFaces = chamfer.Modified(pface->face);
    if (listFaces.Extent() > 0) {
      /* modified Faces */
      TopTools_ListIteratorOfListOfShape it(listFaces);
      for (; it.More(); it.Next()) {
        TopoDS_Face genface = TopoDS::Face(it.Value());
        int index = pbods->faces.map.FindIndex(genface);
        if (index > 0) EG_attributeDup(face, pbods->faces.objs[index-1]);
      }
    }
  }

  // Face map info: source, modified, generated
  if (faceMap != NULL) {
    fmap = (int *) EG_alloc(2*pbods->faces.map.Extent()*sizeof(int));
    if (fmap != NULL) {
      for (i = 0; i < 2*pbods->faces.map.Extent(); i++) fmap[i] = 0;
      for (i = 0; i <   pbody->faces.map.Extent(); i++) {
        face = pbody->faces.objs[i];
        egadsFace *pface = (egadsFace *) face->blind;
        int index = pbods->faces.map.FindIndex(pface->face);
        if (index > 0) {
          fmap[2*index-2] = FACEDUP;
          fmap[2*index-1] = i+1;
        }
      }
      for (i = 0; i < pbody->faces.map.Extent(); i++) {
        face = pbody->faces.objs[i];
        egadsFace *pface = (egadsFace *) face->blind;
        const TopTools_ListOfShape& listFaces = chamfer.Modified(pface->face);
        if (listFaces.Extent() > 0) {
          /* modified Faces */
          TopTools_ListIteratorOfListOfShape it(listFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) {
              fmap[2*index-2] = FACECUT;
              fmap[2*index-1] = i+1;
            }
          }
        }
      }
      for (i = 0; i < pbody->edges.map.Extent(); i++) {
        edge = pbody->edges.objs[i];
        egadsEdge *pedge = (egadsEdge *) edge->blind;
        const TopTools_ListOfShape& genFaces = chamfer.Generated(pedge->edge);
        if (genFaces.Extent() > 0) {
          /* generated faces from Edges */
          TopTools_ListIteratorOfListOfShape it(genFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) {
              fmap[2*index-2] = EDGEOFF;
              fmap[2*index-1] = i+1;
            }
          }
        }
      }
      for (i = 0; i < pbody->nodes.map.Extent(); i++) {
        node = pbody->nodes.objs[i];
        egadsNode *pnode = (egadsNode *) node->blind;
        const TopTools_ListOfShape& genFaces = chamfer.Generated(pnode->node);
        if (genFaces.Extent() > 0) {
          /* generated faces from Nodes */
          TopTools_ListIteratorOfListOfShape it(genFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) {
              fmap[2*index-2] = NODEOFF;
              fmap[2*index-1] = i+1;
            }
          }
        }
      }
      *faceMap = fmap;
    }
  }

  EG_referenceObject(obj, context);
  *result = obj;
  return EGADS_SUCCESS;
}


#ifdef OCC_EXTRUDE
int
EG_extrude(const egObject *src, double dist, const double *dir,
           egObject **result)
{
  int          outLevel, stat, mtype, nerr;
  double       d, vec[3];
  egObject     *context, *obj, *edge;
  TopoDS_Shape shape, newShape;

  *result = NULL;
  if (src == NULL)               return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (src->blind == NULL)        return EGADS_NODATA;
  if (EG_sameThread(src))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);

  mtype = SOLIDBODY;
  if  (src->oclass == BODY) {
    if ((src->mtype == WIREBODY) || (src->mtype == FACEBODY)) {
      egadsBody *pbody = (egadsBody *) src->blind;
      shape = pbody->shape;
      if (src->mtype == WIREBODY) mtype = SHEETBODY;
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Body src must be Wire or Face (EG_extrude)!\n");
      return EGADS_NOTTOPO;
    }
  } else if (src->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) src->blind;
    shape = ploop->loop;
    mtype = SHEETBODY;
  } else if (src->oclass == FACE) {
    egadsFace *pface = (egadsFace *) src->blind;
    shape = pface->face;
  } else {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid src type (EG_extrude)!\n");
    return EGADS_NOTTOPO;
  }

  d = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  if (d == 0.0) {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid Direction (EG_extrude)!\n");
    return EGADS_GEOMERR;
  }
  vec[0] = dist*dir[0]/d;
  vec[1] = dist*dir[1]/d;
  vec[2] = dist*dir[2]/d;
  BRepPrimAPI_MakePrism prism(shape, gp_Vec(vec[0],vec[1],vec[2]));
  try {
    newShape = prism.Shape();
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Error: MakePrism Exception (EG_extrude)!\n");
    printf("              %s\n", e.GetMessageString());
    return EGADS_GEOMERR;
  }
  catch (...) {
    printf(" EGADS Error: MakePrism Exception (EG_extrude)!\n");
    return EGADS_GEOMERR;
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_extrude)!\n");
    return stat;
  }
  egadsBody *pbods   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = mtype;
  pbods->nodes.objs  = NULL;
  pbods->edges.objs  = NULL;
  pbods->loops.objs  = NULL;
  pbods->faces.objs  = NULL;
  pbods->shells.objs = NULL;
  pbods->senses      = NULL;
  pbods->shape       = newShape;
  pbods->bbox.filled = 0;
  pbods->massFill    = 0;
  obj->blind         = pbods;
  stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
  if (stat != EGADS_SUCCESS) {
    EG_deleteObject(obj);
    return stat;
  }

  /* copy attributes from the generating Edges */
  if  (src->oclass == BODY) {
    egadsBody *pbody = (egadsBody *) src->blind;
    for (int i = 0; i < pbody->edges.map.Extent(); i++) {
      edge = pbody->edges.objs[i];
      egadsEdge *pedge = (egadsEdge *) edge->blind;
      const TopTools_ListOfShape& genFaces = prism.Generated(pedge->edge);
      if (genFaces.Extent() > 0) {
        /* generated Faces from Edges */
        TopTools_ListIteratorOfListOfShape it(genFaces);
        for (; it.More(); it.Next()) {
          TopoDS_Shape qFace = it.Value();
          if (qFace.IsNull()) continue;
          if (qFace.ShapeType() != TopAbs_FACE) continue;
          TopoDS_Face genface = TopoDS::Face(qFace);
          int index = pbods->faces.map.FindIndex(genface);
          if (index > 0) EG_attributeDup(edge, pbods->faces.objs[index-1]);
        }
      }
    }
  } else if (src->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) src->blind;
    for (int i = 0; i < ploop->nedges; i++) {
      edge = ploop->edges[i];
      egadsEdge *pedge = (egadsEdge *) edge->blind;
      const TopTools_ListOfShape& genFaces = prism.Generated(pedge->edge);
      if (genFaces.Extent() > 0) {
        /* generated Faces from Edges */
        TopTools_ListIteratorOfListOfShape it(genFaces);
        for (; it.More(); it.Next()) {
          TopoDS_Shape qFace = it.Value();
          if (qFace.IsNull()) continue;
          if (qFace.ShapeType() != TopAbs_FACE) continue;
          TopoDS_Face genface = TopoDS::Face(qFace);
          int index = pbods->faces.map.FindIndex(genface);
          if (index > 0) EG_attributeDup(edge, pbods->faces.objs[index-1]);
        }
      }
    }
  } else {
    egadsFace *pface = (egadsFace *) src->blind;
    for (int j = 0; j < pface->nloops; j++) {
      egadsLoop *ploop = (egadsLoop *) pface->loops[j]->blind;
      for (int i = 0; i < ploop->nedges; i++) {
        edge = ploop->edges[i];
        egadsEdge *pedge = (egadsEdge *) edge->blind;
        const TopTools_ListOfShape& genFaces = prism.Generated(pedge->edge);
        if (genFaces.Extent() > 0) {
          /* generated Faces from Edges */
          TopTools_ListIteratorOfListOfShape it(genFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.IsNull()) continue;
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) EG_attributeDup(edge, pbods->faces.objs[index-1]);
          }
        }
      }
    }
  }

  EG_referenceObject(obj, context);
  *result = obj;
  return EGADS_SUCCESS;
}


int
EG_extrude_dot(egObject *body, const egObject *src,
               double dist, double dist_dot,
               const double *dir, const double *dir_dot)
{
  printf(" EGADS Internal Error: Sensitivities not available for OCC extrude! (EG_extrude_dot)\n");
  return EGADS_GEOMERR;
}

#else

// no-name namespace makes functions private to this file
namespace {

double value(double val)
{
  return val;
}


template <int N>
double value(SurrealS<N> val)
{
  return val.value();
}


template<class T>
int
EG_extrude_plane(egObject *edge, const T* dir, int *mtype, T *data)
{
  int stat = EGADS_SUCCESS;
  T mZ;

  /* The plane is reversed to be consistent with OCC */

  *mtype = -PLANE;

  if (edge->oclass != EDGE) return EGADS_TOPOERR;

  egadsEdge *pedge = (egadsEdge*)edge->blind;

  egObject *ecurve = pedge->curve;

  if (ecurve->oclass != CURVE) return EGADS_GEOMERR;

  egadsCurve *pcurve = (egadsCurve *) ecurve->blind;
  Handle(Geom_Curve) hCurve = pcurve->handle;

  /*
   * from OCC GeomAdaptor_SurfaceOfLinearExtrusion::Plane
   */

  gp_Pnt P;
  T D1u[3], newZ[3];
  double UFirst = pedge->trange[0];
  double ULast  = pedge->trange[1];
  if (Precision::IsNegativeInfinite(UFirst) &&
      Precision::IsPositiveInfinite(ULast)) {
    UFirst = -100.;
    ULast  = 100.;
  }
  else if (Precision::IsNegativeInfinite(UFirst)) {
    UFirst = ULast - 200.;
  }
  else if (Precision::IsPositiveInfinite(ULast)) {
    ULast = UFirst + 200.;
  }
  double deltau = (ULast-UFirst)/20.;
  for (int i = 1; i <= 21; i++) {
    T prm = UFirst + (i-1)*deltau;
    stat = EG_evaluate(edge, &prm, data);
    if (stat != EGADS_SUCCESS) return stat;

    D1u[0] = data[3];
    D1u[1] = data[4];
    D1u[2] = data[5];
    CROSS(newZ, D1u, dir);
    mZ = sqrt(DOT(newZ, newZ));
    if (mZ > 1.e-12) break;
  }
//data[0]           /* center from evaluate */
//data[1]
//data[2]
  data[3] = D1u[0]; /* x-axis */
  data[4] = D1u[1];
  data[5] = D1u[2];
  data[6] = -dir[0]; /* y-axis */
  data[7] = -dir[1];
  data[8] = -dir[2];

  gp_Dir Dir(value(dir[0]), value(dir[1]), value(dir[2]));
  gp_Ax3 Ax3(P,gp_Dir(value(newZ[0]), value(newZ[1]), value(newZ[2])),
               gp_Dir(value( D1u[0]), value( D1u[2]), value( D1u[2])));
  if (Dir.Dot(Ax3.YDirection()) < 0.0){
    data[6] = dir[0]; /* y-axis */
    data[7] = dir[1];
    data[8] = dir[2];
    *mtype = -(*mtype);
  }

  return stat;
}


template<class T>
int
EG_extrude_surf(egObject *edge, int *mtype, const T* dir, T *data)
{
  T vec[3];

  if (edge->oclass != EDGE) return EGADS_TOPOERR;

  egObject *ecurve = ((egadsEdge*)edge->blind)->curve;

  if (ecurve->oclass != CURVE) return EGADS_GEOMERR;

  egadsCurve *pcurve = (egadsCurve *) ecurve->blind;
  Handle(Geom_Curve) hCurve = pcurve->handle;

  /*
   * from OCC BRepSweep_Translation::MakeEmptyFace
   *  "extruded surfaces are inverted correspondingly to the topology, so reverse."
   *
   * However, reversing the direction makes the topology unnecessarily convoluted
   *
   */
  vec[0] = dir[0];
  vec[1] = dir[1];
  vec[2] = dir[2];

  gp_Dir Dir(value(vec[0]), value(vec[1]), value(vec[2]));

  if (ecurve->mtype == LINE) {
    Handle(Geom_Line) hLine = Handle(Geom_Line)::DownCast(hCurve);

    gp_Dir D = hLine->Lin().Direction();
    if (!Dir.IsParallel( D, Precision::Angular()))
      return EG_extrude_plane(edge, vec, mtype, data);
  }

  /* create the Extruded surface */
  *mtype  = EXTRUSION;
  data[0] = vec[0]; /* direction */
  data[1] = vec[1];
  data[2] = vec[2];

  return EGADS_SUCCESS;
}


template<class T>
int
EG_extrude_line(egObject* node, const T* dir, T *data)
{
  int status = EGADS_SUCCESS;
  T xyz[3];

  status = EG_evaluate(node, NULL, xyz);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* create the Line (point and direction) */
  data[0] = xyz[0];
  data[1] = xyz[1];
  data[2] = xyz[2];
  data[3] = dir[0];
  data[4] = dir[1];
  data[5] = dir[2];

cleanup:
  return status;
}


int
EG_extrude_loop(objStack *stack, const egObject *loop0, const egObject *loop1,
                const double dist, const double* dir, egObject **faces)
{
  int      stat = EGADS_SUCCESS;
  int      i, mtype, lsens[1] = {SFORWARD}, *ssens;
  int      esens[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE};
  double   data[18], tdata[2] = {0, dist};
  egObject **surfs=NULL, *line, **edges=NULL, *context, *enodes[2];
  egObject *eedges[8], *eloop, *eface;

  if (loop0->oclass != LOOP) return EGADS_TOPOERR;
  if (loop1->oclass != LOOP) return EGADS_TOPOERR;

  context  = EG_context(loop0);

  egadsLoop *ploop0 = (egadsLoop *) loop0->blind;
  egadsLoop *ploop1 = (egadsLoop *) loop1->blind;

  edges = (egObject **) EG_alloc((ploop0->nedges+1)*sizeof(egObject *));
  surfs = (egObject **) EG_alloc( ploop0->nedges   *sizeof(egObject *));
  ssens = (int *)       EG_alloc( ploop0->nedges   *sizeof(int));

  /* construct new geometry and new edges */
  for (i = 0; i < ploop0->nedges; i++) {
    egadsEdge *pedge0 = (egadsEdge *) ploop0->edges[i]->blind;
    egadsEdge *pedge1 = (egadsEdge *) ploop1->edges[i]->blind;

    /* make the surface */
    stat = EG_extrude_surf(ploop0->edges[i], &mtype, dir, data);
    if (stat != EGADS_SUCCESS) goto cleanup;

    ssens[i] = mtype/abs(mtype);
    mtype = abs(mtype);

    if (mtype == EXTRUSION)
      stat = EG_makeGeometry(context, SURFACE, mtype, pedge0->curve, NULL,
                             data, &surfs[i]);
    else
      stat = EG_makeGeometry(context, SURFACE, mtype, NULL, NULL,
                             data, &surfs[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, surfs[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (ploop0->senses[i] == SFORWARD) {
      enodes[0] = pedge0->nodes[0];
      enodes[1] = pedge1->nodes[0];
    } else {
      enodes[0] = pedge0->nodes[1];
      enodes[1] = pedge1->nodes[1];
    }

    /* make a line */
    stat = EG_extrude_line(enodes[0], dir, data);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                           data, &line);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, line);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, line, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, edges[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  if (loop0->mtype == OPEN) {
    /* create the last line/edge for open loop */
    egadsEdge *pedge0 = (egadsEdge *) ploop0->edges[i-1]->blind;
    egadsEdge *pedge1 = (egadsEdge *) ploop1->edges[i-1]->blind;

    if (ploop0->senses[i-1] == SFORWARD) {
      enodes[0] = pedge0->nodes[1];
      enodes[1] = pedge1->nodes[1];
    } else {
      enodes[0] = pedge0->nodes[0];
      enodes[1] = pedge1->nodes[0];
    }

    stat = EG_extrude_line(enodes[0], dir, data);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                           data, &line);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, line);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, line, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, edges[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
  } else {
    /* repeat edge for closed loop */
    edges[i] = edges[0];
  }

  /* construct loops/faces */
  for (i = 0; i < ploop0->nedges; i++) {

    egadsEdge *pedge0 = (egadsEdge *) ploop0->edges[i]->blind;

    eedges[0] = ploop0->senses[i] == SFORWARD ? edges[i] : edges[i+1];
    eedges[1] = ploop0->edges[i];
    eedges[2] = ploop0->senses[i] == SFORWARD ? edges[i+1] : edges[i];
    eedges[3] = ploop1->edges[i];

    if (surfs[i]->mtype == PLANE) {

      /* make the loop */
      stat = EG_makeTopology(context, NULL, LOOP, CLOSED,
                             NULL, 4, eedges, esens, &eloop);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(stack, eloop);
      if (stat != EGADS_SUCCESS) goto cleanup;

    } else {
      tdata[0] = pedge0->trange[0];
      tdata[1] = pedge0->trange[1];

      /* create P-curves */
      data[0] = tdata[0]; data[1] = 0.;    /* u == UMIN */
      data[2] = 0.;       data[3] = 1.;
      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                             &eedges[4+0]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(stack, eedges[4+0]);
      if (stat != EGADS_SUCCESS) goto cleanup;

      data[0] = 0.;       data[1] = 0.;    /* v == VMIN */
      data[2] = 1.;       data[3] = 0.;
      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                             &eedges[4+1]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(stack, eedges[4+1]);
       if (stat != EGADS_SUCCESS) goto cleanup;

      data[0] = tdata[1]; data[1] = 0.;    /* u == UMAX */
      data[2] = 0.;       data[3] = 1.;
      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                             &eedges[4+2]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(stack, eedges[4+2]);
       if (stat != EGADS_SUCCESS) goto cleanup;

      data[0] = 0.;       data[1] = dist;  /* v == VMAX */
      data[2] = 1.;       data[3] = 0.;
      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                             &eedges[4+3]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(stack, eedges[4+3]);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* make the loop */
      stat = EG_makeTopology(context, surfs[i], LOOP, CLOSED,
                             NULL, 4, eedges, esens, &eloop);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(stack, eloop);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    /* make the face */
    stat = EG_makeTopology(context, surfs[i], FACE, ssens[i],
                           NULL, 1, &eloop, lsens, &eface);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, eface);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* copy generating edge attributes to the face */
    EG_attributeDup(ploop0->edges[i], eface);

    /* set the output */
    faces[i] = eface;
  }

cleanup:
  EG_free(edges);
  EG_free(surfs);
  EG_free(ssens);

  return stat;
}


int
EG_extrude_loop_dot(const egObject *loop0, const SurrealS<1>& dist,
                    const SurrealS<1> *dir, egObject **efaces)
{
  int         stat = EGADS_SUCCESS;
  int         i, j, k, mtype;
  int         iedge[4];
  SurrealS<1> data[18], tdata[2] = {0, dist};
  egObject    *eline, **eedges=NULL, *enodes[2], *esurf, *eloop;

  if (loop0->oclass != LOOP) return EGADS_TOPOERR;

  SurrealS<1> mat[12] = {1.00, 0.00, 0.00, dist*dir[0],
                         0.00, 1.00, 0.00, dist*dir[1],
                         0.00, 0.00, 1.00, dist*dir[2]};

  egadsLoop *ploop0 = (egadsLoop*) loop0->blind;

  /* construct new geometry and new edges */
  for (i = 0; i < ploop0->nedges; i++) {
    egadsEdge *pedge0 = (egadsEdge *) ploop0->edges[i]->blind;

    /* get the surface and loop for the face */
    egadsFace *pface  = (egadsFace *) efaces[i]->blind;
    if (pface->nloops != 1) {
      stat = EGADS_TOPOERR;
      goto cleanup;
    }
    eloop = pface->loops[0];
    esurf = pface->surface;

    /* get the edges from the loop */
    egadsLoop *ploop = (egadsLoop *) eloop->blind;
    if (ploop->nedges != 4) {
      stat = EGADS_TOPOERR;
      goto cleanup;
    }
    eedges = ploop->edges;

    /* get the surface data */
    stat = EG_extrude_surf(ploop0->edges[i], &mtype, dir, data);
    if (stat != EGADS_SUCCESS) goto cleanup;
    mtype = abs(mtype);

    /* set the sensitivity for the surface */
    stat = EG_setGeometry_dot(esurf, SURFACE, mtype, NULL, data);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (abs(mtype) == EXTRUSION) {
      /* copy the curve velocity */
      egadsSurface *psurf = (egadsSurface *) esurf->blind;
      stat = EG_copyGeometry_dot(pedge0->curve, NULL, psurf->ref);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    if (ploop0->senses[i] == SFORWARD) {
      enodes[0] = pedge0->nodes[0];
      enodes[1] = pedge0->nodes[1];
    } else {
      enodes[0] = pedge0->nodes[1];
      enodes[1] = pedge0->nodes[0];
    }

    /* reconstruct the original loop ordering in EG_extrude_loop */
    iedge[1] = -1;
    for (j = 0; j < 4; j++) {
      if (pedge0->edge.IsSame(((egadsEdge*)eedges[j]->blind)->edge)) {
        iedge[1] = j;
        break;
      }
    }
    if (iedge[1] == -1) {
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    {
      egadsNode *pnode0 = (egadsNode *) enodes[0]->blind;
      iedge[0] = -1;
      for (j = 0; j < 4; j++) {
        if (j == iedge[1]) continue;
        egadsEdge *pedge = (egadsEdge *) eedges[j]->blind;
        if (pnode0->node.IsSame(((egadsNode *) pedge->nodes[0]->blind)->node)) {
          iedge[0] = j;
          break;
        }
      }
      if (iedge[0] == -1) {
        stat = EGADS_TOPOERR;
        goto cleanup;
      }
    }

    {
      egadsNode *pnode0 = (egadsNode *) enodes[1]->blind;
      iedge[2] = -1;
      for (j = 0; j < 4; j++) {
        if (j == iedge[0]) continue;
        if (j == iedge[1]) continue;
        egadsEdge *pedge = (egadsEdge *) eedges[j]->blind;
        if (pnode0->node.IsSame(((egadsNode *) pedge->nodes[0]->blind)->node)) {
          iedge[2] = j;
          break;
        }
      }
      if (iedge[2] == -1) {
        stat = EGADS_TOPOERR;
        goto cleanup;
      }
    }

    {
      iedge[3] = -1;
      egadsEdge *pedgeLine = (egadsEdge *) eedges[iedge[0]]->blind;
      egadsNode *pnode     = (egadsNode *) pedgeLine->nodes[1]->blind;
      k = ploop0->senses[i] == SFORWARD ? 0 : 1;
      for (j = 0; j < 4; j++) {
        if (j == iedge[0]) continue;
        if (j == iedge[1]) continue;
        if (j == iedge[2]) continue;
        egadsEdge *pedge = (egadsEdge *) eedges[j]->blind;
        if (pnode->node.IsSame(((egadsNode *) pedge->nodes[k]->blind)->node)) {
          iedge[3] = j;
          break;
        }
      }
      if (iedge[3] == -1) {
        stat = EGADS_TOPOERR;
        goto cleanup;
      }
    }

    /* get the line data */
    stat = EG_extrude_line(enodes[0], dir, data);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* get the line from the edge */
    egadsEdge *pedge = (egadsEdge *) eedges[iedge[0]]->blind;
    eline = pedge->curve;
    if (eline->mtype != LINE) {
      stat = EGADS_GEOMERR;
      goto cleanup;
    }

    /* set the sensitivity for the line, nodes and t-range */
    stat = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_copyGeometry_dot(enodes[0], NULL, pedge->nodes[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_copyGeometry_dot(enodes[0], mat, pedge->nodes[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_setGeometry_dot(eedges[iedge[0]], EDGE, TWONODE, NULL, tdata);
    if (stat != EGADS_SUCCESS) goto cleanup;


    /* copy sensitivities from the src curve (not nodes) */
    stat = EG_copyGeometry_dot(ploop0->edges[i], NULL, eedges[iedge[1]]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_copyGeometry_dot(ploop0->edges[i],  mat, eedges[iedge[3]]);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  if (loop0->mtype == OPEN) {
    /* set the last line/edge for open loop */

    /* get the line data */
    stat = EG_extrude_line(enodes[1], dir, data);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* get the line from the edge */
    egadsEdge *pedge = (egadsEdge *) eedges[iedge[2]]->blind;
    eline = pedge->curve;
    if (eline->mtype != LINE) {
      stat = EGADS_GEOMERR;
      goto cleanup;
    }

    /* set the sensitivity for the line and nodes */
    stat = EG_setGeometry_dot(eline, CURVE, LINE, NULL, data);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_copyGeometry_dot(enodes[1], NULL, pedge->nodes[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_copyGeometry_dot(enodes[1], mat, pedge->nodes[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_setGeometry_dot(eedges[iedge[2]], EDGE, TWONODE, NULL, tdata);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

cleanup:
  return stat;
}
} // no-name namespace


int
EG_extrude(const egObject *src, double dist, const double *dir,
           egObject **result)
{
  int            outLevel, stat, i, mtype, nfaces, nloops, offset;
  double         d, vec[3];
  egObject       *context, *ref, **loops1=NULL, **faces, *efaces[2];
  egObject       *shell, *body, *exform, *ecopy;
  const egObject **loops0=NULL;
  objStack       stack;

  nfaces = 2;

  *result = NULL;
  if (src == NULL)               return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (src->blind == NULL)        return EGADS_NODATA;
  if (EG_sameThread(src))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);

  if (dist == 0.0) {
    if (outLevel > 0)
      printf(" EGADS Error: Zero Distance (EG_extrude)!\n");
    return EGADS_GEOMERR;
  }

  d = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  if (d == 0.0) {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid Direction (EG_extrude)!\n");
    return EGADS_GEOMERR;
  }
  vec[0] = dir[0]/d;
  vec[1] = dir[1]/d;
  vec[2] = dir[2]/d;

  if (dist < 0) {
    dist = -dist;
    vec[0] = -vec[0];
    vec[1] = -vec[1];
    vec[2] = -vec[2];
  }

  double mat[12] = {1.00, 0.00, 0.00, dist*vec[0],
                    0.00, 1.00, 0.00, dist*vec[1],
                    0.00, 0.00, 1.00, dist*vec[2]};

  /* create stack for gracefully cleaning up objects */
  stat = EG_stackInit(&stack);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_makeTransform(context, mat, &exform);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, exform);
  if (stat != EGADS_SUCCESS) goto cleanup;

  mtype = SOLIDBODY;
  if  (src->oclass == BODY) {
    if ((src->mtype == WIREBODY) || (src->mtype == FACEBODY)) {
      egadsBody *pbody = (egadsBody *) src->blind;
      nloops = pbody->loops.map.Extent();
      loops0 = const_cast<const egObject**>(pbody->loops.objs);

      if (src->mtype == WIREBODY) {
        mtype  = SHEETBODY;
        nfaces = 0;
        stat   = EG_copyObject(loops0[0], exform, &ecopy);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat   = EG_stackPush(&stack, ecopy);
        if (stat != EGADS_SUCCESS) goto cleanup;

        loops1 = &ecopy;
      } else {
        efaces[0] = pbody->faces.objs[0];
        stat = EG_copyObject(efaces[0], exform, &efaces[1]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, efaces[1]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        egadsFace *pface = (egadsFace *) efaces[1]->blind;
        loops1 = pface->loops;
      }
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Body src must be Wire or Face (EG_extrude)!\n");
      stat = EGADS_NOTTOPO;
      goto cleanup;
    }
  } else if (src->oclass == LOOP) {
    nloops = 1;
    loops0 = &src;
    mtype  = SHEETBODY;
    nfaces = 0;

    stat = EG_copyObject(src, exform, &ecopy);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, ecopy);
    if (stat != EGADS_SUCCESS) goto cleanup;

    loops1 = &ecopy;
  } else if (src->oclass == FACE) {
    egadsFace *pface = (egadsFace *) src->blind;
    nloops = pface->nloops;
    loops0 = const_cast<const egObject**>(pface->loops);

    efaces[0] = const_cast<egObject*>(src);
    stat = EG_copyObject(efaces[0], exform, &efaces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, efaces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    pface  = (egadsFace *) efaces[1]->blind;
    loops1 = pface->loops;
  } else {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid src type (EG_extrude)!\n");
    stat = EGADS_NOTTOPO;
    goto cleanup;
  }


  for (i = 0; i < nloops; i++) {
    egadsLoop *ploop = (egadsLoop *) loops0[i]->blind;
    nfaces += ploop->nedges;
  }
  faces = (egObject **) EG_alloc(nfaces*sizeof(egObject *));

  offset = 0;
  for (i = 0; i < nloops; i++) {
    stat = EG_extrude_loop(&stack, loops0[i], loops1[i],
                           dist, vec, faces + offset);
    if (stat != EGADS_SUCCESS) goto cleanup;

    offset += ((egadsLoop *) loops0[i]->blind)->nedges;
  }

  if (mtype == SOLIDBODY) {
    faces[offset+0] = efaces[0];
    faces[offset+1] = efaces[1];
  }

  /* put all the faces together */
  stat = EG_makeTopology(context, NULL, SHELL, mtype == SOLIDBODY ? CLOSED : OPEN,
                         NULL, nfaces, faces, NULL, &shell);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeTopo Shell = %d (EG_extrude)!\n", stat);
    goto cleanup;
  }
  stat = EG_stackPush(&stack, shell);
  if (stat != EGADS_SUCCESS) goto cleanup;

  body = NULL;
  stat = EG_makeTopology(context, NULL, BODY, mtype, NULL, 1,
                         &shell, NULL, &body);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeTopo Body = %d (EG_extrude)!\n", stat);
    if (body != NULL) EG_deleteObject(body);
    goto cleanup;
  }

  stat = EGADS_SUCCESS;
  *result = body;

cleanup:
  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    int i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (EG_extrude)!\n", i);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);
  EG_free(faces);

  return stat;
}


int
EG_extrude_dot(egObject *body, const egObject *src,
               double dist, double dist_dot,
               const double *dir, const double *dir_dot)
{
  int            outLevel, stat, i, mtype, nfaces, nloops, offset;
  SurrealS<1>    d, vecS[3], dirS[3], distS;
  egObject       **faces;
  const egObject **loops0=NULL, *srcFace;
  egadsBody      *pbody;

  nfaces = 2;

  if (src == NULL)                  return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (src->blind == NULL)           return EGADS_NODATA;
  if (EG_sameThread(src))           return EGADS_CNTXTHRD;

  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->blind == NULL)        return EGADS_NODATA;
  if (body->oclass != BODY)       return EGADS_NOTBODY;
  if (EG_sameThread(body))        return EGADS_CNTXTHRD;

  outLevel = EG_outLevel(src);

  if (EG_hasGeometry_dot(src) != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: src without data_dot (EG_extrude_dot)!\n");
    return EGADS_NODATA;
  }

  if (dist == 0.0) {
    if (outLevel > 0)
      printf(" EGADS Error: Zero Distance (EG_extrude_dot)!\n");
    return EGADS_GEOMERR;
  }

  distS.value() = dist;
  distS.deriv() = dist_dot;

  for (i = 0; i < 3; i++) {
    dirS[i].value() = dir[i];
    dirS[i].deriv() = dir_dot[i];
  }

  d = sqrt(dirS[0]*dirS[0] + dirS[1]*dirS[1] + dirS[2]*dirS[2]);
  if (d == 0.0) {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid Direction (EG_extrude_dot)!\n");
    return EGADS_GEOMERR;
  }
  vecS[0] = dirS[0]/d;
  vecS[1] = dirS[1]/d;
  vecS[2] = dirS[2]/d;

  if (distS < 0) {
    distS = -distS;
    vecS[0] = -vecS[0];
    vecS[1] = -vecS[1];
    vecS[2] = -vecS[2];
  }

  SurrealS<1> mat[12] = {1.00, 0.00, 0.00, distS*vecS[0],
                         0.00, 1.00, 0.00, distS*vecS[1],
                         0.00, 0.00, 1.00, distS*vecS[2]};


  mtype = SOLIDBODY;
  if  (src->oclass == BODY) {
    if ((src->mtype == WIREBODY) || (src->mtype == FACEBODY)) {
      egadsBody *pbody = (egadsBody *) src->blind;
      nloops = pbody->loops.map.Extent();
      loops0 = const_cast<const egObject**>(pbody->loops.objs);

      if (src->mtype == WIREBODY) {
        mtype = SHEETBODY;
        nfaces = 0;
      }
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: src must be Wire or Face (EG_extrude_dot)!\n");
      stat = EGADS_NOTTOPO;
      goto cleanup;
    }
  } else if (src->oclass == LOOP) {
    nloops = 1;
    loops0 = &src;
    mtype = SHEETBODY;
    nfaces = 0;
  } else if (src->oclass == FACE) {
    egadsFace *pface = (egadsFace *) src->blind;
    nloops = pface->nloops;
    loops0 = const_cast<const egObject**>(pface->loops);
  } else {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid src type (EG_extrude_dot)!\n");
    stat = EGADS_NOTTOPO;
    goto cleanup;
  }

  if (body->mtype != mtype) {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid result type (EG_extrude_dot)!\n");
    stat = EGADS_NOTTOPO;
    goto cleanup;
  }

  for (i = 0; i < nloops; i++) {
    egadsLoop *ploop = (egadsLoop *) loops0[i]->blind;
    nfaces += ploop->nedges;
  }

  /* get the Faces from the Body */
  pbody = (egadsBody *) body->blind;
  if (pbody->faces.map.Extent() != nfaces) {
    if (outLevel > 0)
      printf(" EGADS Error: Body src is not consistent with result body (EG_extrude_dot)!\n");
    stat = EGADS_TOPOERR;
    goto cleanup;
  }
  faces = pbody->faces.objs;

  offset = 0;
  for (i = 0; i < nloops; i++) {
    stat = EG_extrude_loop_dot(loops0[i], distS, vecS, faces + offset);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Body src is not consistent with result body (EG_extrude_dot)!\n");
      goto cleanup;
    }

    offset += ((egadsLoop *) loops0[i]->blind)->nedges;
  }

  if (mtype == SOLIDBODY) {

    if  (src->oclass == BODY) {
      egadsBody *pbody = (egadsBody *) src->blind;
      srcFace = pbody->faces.objs[0];
    } else {
      srcFace = src;
    }

    egadsFace *psrcFace = (egadsFace *) srcFace->blind;

    stat = EG_copyGeometry_dot(psrcFace->surface, NULL,
                               ((egadsFace*)faces[offset+0]->blind)->surface);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_copyGeometry_dot(psrcFace->surface,  mat,
                               ((egadsFace*)faces[offset+1]->blind)->surface);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  stat = EGADS_SUCCESS;

cleanup:
  return stat;
}
#endif


#ifdef OCC_ROTATE
int
EG_rotate(const egObject *src, double angle, const double *axis,
          egObject **result)
{
  int          outLevel, stat, mtype, nerr;
  double       ang;
  egObject     *context, *obj, *edge;
  TopoDS_Shape shape, newShape;

  *result = NULL;
  if (src == NULL)               return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (src->blind == NULL)        return EGADS_NODATA;
  if (EG_sameThread(src))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);

  mtype = SOLIDBODY;
  if  (src->oclass == BODY) {
    if ((src->mtype == WIREBODY) || (src->mtype == FACEBODY)) {
      egadsBody *pbody = (egadsBody *) src->blind;
      shape = pbody->shape;
      if (src->mtype == WIREBODY) mtype = SHEETBODY;
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Body src must be Wire or Face (EG_rotate)!\n");
      return EGADS_NOTTOPO;
    }
  } else if (src->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) src->blind;
    shape = ploop->loop;
    mtype = SHEETBODY;
  } else if (src->oclass == FACE) {
    egadsFace *pface = (egadsFace *) src->blind;
    shape = pface->face;
  } else {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid src type (EG_rotate)!\n");
    return EGADS_NOTTOPO;
  }

  gp_Pnt pnt(axis[0], axis[1], axis[2]);
  gp_Dir dir(axis[3], axis[4], axis[5]);
  gp_Ax1 axi(pnt, dir);
  ang = 360.0;
  if ((angle > -360.0) && (angle < 360.0)) ang = angle;
  BRepPrimAPI_MakeRevol revol(shape, axi, ang*PI/180.0);
  try {
    newShape = revol.Shape();
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Error: MakeRevol Exception (EG_rotate)!\n");
    printf("              %s\n", e.GetMessageString());
    return EGADS_GEOMERR;
  }
  catch (...) {
    printf(" EGADS Error: MakeRevol Exception (EG_rotate)!\n");
    return EGADS_GEOMERR;
  }

  if (mtype == SOLIDBODY) {
    GProp_GProps VProps;
    BRepGProp    BProps;

    BProps.VolumeProperties(newShape, VProps);
    if (VProps.Mass() < 0.0) {
      if (outLevel > 0)
        printf(" EGADS Info: Volume = %lf reversing (EG_rotate)!\n",
               VProps.Mass());
      newShape.Reverse();
    }
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_rotate)!\n");
    return stat;
  }
  egadsBody *pbods   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = mtype;
  pbods->nodes.objs  = NULL;
  pbods->edges.objs  = NULL;
  pbods->loops.objs  = NULL;
  pbods->faces.objs  = NULL;
  pbods->shells.objs = NULL;
  pbods->senses      = NULL;
  pbods->shape       = newShape;
  pbods->bbox.filled = 0;
  pbods->massFill    = 0;
  obj->blind         = pbods;
  stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbods;
    return stat;
  }

  /* copy attributes from the generating Edges */
  if (src->oclass == BODY) {
    egadsBody *pbody = (egadsBody *) src->blind;
    for (int i = 0; i < pbody->edges.map.Extent(); i++) {
      edge = pbody->edges.objs[i];
      egadsEdge *pedge = (egadsEdge *) edge->blind;
      const TopTools_ListOfShape& genFaces = revol.Generated(pedge->edge);
      if (genFaces.Extent() > 0) {
        /* generated Faces from Edges */
        TopTools_ListIteratorOfListOfShape it(genFaces);
        for (; it.More(); it.Next()) {
          TopoDS_Shape qFace = it.Value();
          if (qFace.ShapeType() != TopAbs_FACE) continue;
          TopoDS_Face genface = TopoDS::Face(qFace);
          int index = pbods->faces.map.FindIndex(genface);
          if (index > 0) EG_attributeDup(edge, pbods->faces.objs[index-1]);
        }
      }
    }
  } else if (src->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) src->blind;
    for (int i = 0; i < ploop->nedges; i++) {
      edge = ploop->edges[i];
      egadsEdge *pedge = (egadsEdge *) edge->blind;
      const TopTools_ListOfShape& genFaces = revol.Generated(pedge->edge);
      if (genFaces.Extent() > 0) {
        /* generated Faces from Edges */
        TopTools_ListIteratorOfListOfShape it(genFaces);
        for (; it.More(); it.Next()) {
          TopoDS_Shape qFace = it.Value();
          if (qFace.ShapeType() != TopAbs_FACE) continue;
          TopoDS_Face genface = TopoDS::Face(qFace);
          int index = pbods->faces.map.FindIndex(genface);
          if (index > 0) EG_attributeDup(edge, pbods->faces.objs[index-1]);
        }
      }
    }
  } else {
    egadsFace *pface = (egadsFace *) src->blind;
    for (int j = 0; j < pface->nloops; j++) {
      egadsLoop *ploop = (egadsLoop *) pface->loops[j]->blind;
      for (int i = 0; i < ploop->nedges; i++) {
        edge = ploop->edges[i];
        egadsEdge *pedge = (egadsEdge *) edge->blind;
        const TopTools_ListOfShape& genFaces = revol.Generated(pedge->edge);
        if (genFaces.Extent() > 0) {
          /* generated Faces from Edges */
          TopTools_ListIteratorOfListOfShape it(genFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) EG_attributeDup(edge, pbods->faces.objs[index-1]);
          }
        }
      }
    }
  }

  EG_referenceObject(obj, context);
  *result = obj;
  return EGADS_SUCCESS;
}
#else

extern "C" int  EG_saveModel( const egObject *model, const char *name );

// no-name namespace makes functions private to this file
namespace {

template<class T>
void
EG_rotationMatrix(const T& angle, const T *axis, T *mat)
{
  /* unit axis */
  T ux = axis[3];
  T uy = axis[4];
  T uz = axis[5];

  T mag = sqrt(ux*ux + uy*uy + uz*uz);

  ux /= mag;
  uy /= mag;
  uz /= mag;

  T cosa = cos(angle);
  T sina = sin(angle);

  /* rotation matrix about axis (not including the point */
  T R[3][3];

  R[0][0] = cosa + ux*ux*(1 - cosa);
  R[0][1] = ux*uy*(1 - cosa) - uz*sina;
  R[0][2] = ux*uz*(1 - cosa) + uy*sina;

  R[1][0] = uy*ux*(1 - cosa) + uz*sina;
  R[1][1] = cosa + uy*uy*(1 - cosa);
  R[1][2] = uy*uz*(1 - cosa) - ux*sina;

  R[2][0] = uz*ux*(1 - cosa) - uy*sina;
  R[2][1] = uz*uy*(1 - cosa) + ux*sina;
  R[2][2] = cosa + uz*uz*(1 - cosa);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      mat[4*i + j] = R[i][j];

  /* set the translation to include the point */
  mat[4*0 + 3] = axis[0];
  mat[4*1 + 3] = axis[1];
  mat[4*2 + 3] = axis[2];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      mat[4*i + 3] -= R[i][j]*axis[j];
}

#if 0
template<class T>
void
EG_Ax2(T *dirz, T *dirx, T *diry)
{
  T A, Aabs, B, Babs, C, Cabs;

  /* taken from gp_Ax2 constructor
   *
   * This logic introduces discontinuities for finite differencing though...
   */
  A = dirz[0];
  B = dirz[1];
  C = dirz[2];
  Aabs = A;
  if (Aabs < 0) Aabs = - Aabs;
  Babs = B;
  if (Babs < 0) Babs = - Babs;
  Cabs = C;
  if (Cabs < 0) Cabs = - Cabs;

  if      ( Babs <= Aabs && Babs <= Cabs) {
    if (Aabs > Cabs) { dirx[0] = -C; dirx[1] = 0.; dirx[2] =  A; }
    else             { dirx[0] =  C; dirx[1] = 0.; dirx[2] = -A; }
  }
  else if ( Aabs <= Babs && Aabs <= Cabs) {
    if (Babs > Cabs) { dirx[0] = 0.; dirx[1] = -C; dirx[2] =  B; }
    else             { dirx[0] = 0.; dirx[1] =  C; dirx[2] = -B; }
  }
  else {
    if (Aabs > Babs) { dirx[0] = -B; dirx[1] =  A; dirx[2] = 0.; }
    else             { dirx[0] =  B; dirx[1] = -A; dirx[2] = 0.; }
  }

  CROSS(diry, dirz, dirx);
}

template<class T>
int
EG_rotate_surf(egObject *edge, int *mtype, const T *axis, T *data)
{
  if (edge->oclass != EDGE) return EGADS_TOPOERR;

  egObject *ecurve = ((egadsEdge*)edge->blind)->curve;

  if (ecurve->oclass != CURVE) return EGADS_GEOMERR;

  egadsCurve *pcurve = (egadsCurve *) ecurve->blind;
  Handle(Geom_Curve) hCurve = pcurve->handle;

  gp_Dir Dir(value(vec[0]), value(vec[1]), value(vec[2]));

  if (ecurve->mtype == LINE) {
    Handle(Geom_Line) hLine = Handle(Geom_Line)::DownCast(hCurve);

    gp_Dir D = hLine->Lin().Direction();
    if (!Dir.IsParallel( D, Precision::Angular()))
      return EG_extrude_plane(edge, vec, mtype, data);
  }

    /* create the Revolved surface */
  *mtype = REVOLUTION;
  data[0] = axis[0]; /* center */
  data[1] = axis[1];
  data[2] = axis[2];
  data[3] = axis[3]; /* direction */
  data[4] = axis[4];
  data[5] = axis[5];

  return EGADS_SUCCESS;
}
#endif

template<class T>
int
EG_rotate_circle(egObject* node, const T *axis, T *data)
{
  int status = EGADS_SUCCESS;
  T xyz[3], r, dirx[3], diry[3], vec[3], cent[3], dot;

  status = EG_evaluate(node, NULL, xyz);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* find the center on the axis perpendicular to the node */
  vec[0] = xyz[0] - axis[0];
  vec[1] = xyz[1] - axis[1];
  vec[2] = xyz[2] - axis[2];

  dot = DOT(vec, (axis+3));

  cent[0] = axis[0] + dot*axis[3];
  cent[1] = axis[1] + dot*axis[4];
  cent[2] = axis[2] + dot*axis[5];

  /* compute the x-axis and radius */
  dirx[0] = xyz[0] - cent[0];
  dirx[1] = xyz[1] - cent[1];
  dirx[2] = xyz[2] - cent[2];

  CROSS(diry, (axis+3), dirx);

  r = sqrt(DOT(dirx, dirx));

  /* the Circle data */
  data[0] = cent[0]; /* center */
  data[1] = cent[1];
  data[2] = cent[2];
  data[3] = dirx[0]; /* x-axis */
  data[4] = dirx[1];
  data[5] = dirx[2];
  data[6] = diry[0]; /* y-axis */
  data[7] = diry[1];
  data[8] = diry[2];
  data[9] = r;

cleanup:
  return status;
}


int
EG_rotate_loop(objStack *stack, const egObject *loop0, const egObject *loop1,
               const double angle, const double *axis, egObject **faces)
{
  int      stat = EGADS_SUCCESS;
  int      i, lsens[1] = {SFORWARD}, *ssens;
  int      esens[4] = {SREVERSE, SFORWARD, SFORWARD, SREVERSE};
  double   data[18], tdata[2];
  egObject **surfs=NULL, *circle, **edges=NULL, *context, *enodes[2];
  egObject *eedges[8], *eloop, *eface, *trimmed;

  if (angle > 0) {
    tdata[0] = 0;
    tdata[1] = angle;
  } else {
    tdata[0] = angle;
    tdata[1] = 0;
  }

  if (loop0->oclass != LOOP) return EGADS_TOPOERR;
  if (loop1->oclass != LOOP) return EGADS_TOPOERR;

  context  = EG_context(loop0);

  egadsLoop *ploop0 = (egadsLoop *) loop0->blind;
  egadsLoop *ploop1 = (egadsLoop *) loop1->blind;

  edges = (egObject **) EG_alloc((ploop0->nedges+1)*sizeof(egObject *));
  surfs = (egObject **) EG_alloc( ploop0->nedges   *sizeof(egObject *));
  ssens = (int *)       EG_alloc( ploop0->nedges   *sizeof(int));

  /* construct new geometry and new edges */
  for (i = 0; i < ploop0->nedges; i++) {
    egadsEdge *pedge0 = (egadsEdge *) ploop0->edges[i]->blind;
    egadsEdge *pedge1 = (egadsEdge *) ploop1->edges[i]->blind;

    ssens[i] = SFORWARD;

    /* trim the curve to be consistent with OCC */
    stat = EG_makeGeometry(context, CURVE, TRIMMED, pedge0->curve, NULL,
                           pedge0->trange, &trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the surface */
    stat = EG_makeGeometry(context, SURFACE, REVOLUTION, trimmed, NULL,
                           axis, &surfs[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, surfs[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (ploop0->senses[i] == SFORWARD) {
      enodes[0] = pedge0->nodes[0];
      enodes[1] = pedge1->nodes[0];
    } else {
      enodes[0] = pedge0->nodes[1];
      enodes[1] = pedge1->nodes[1];
    }

    /* degeneracy check from BRepSweep_Rotation::MakeEmptyDirectingEdge */
    gp_Pnt P = BRep_Tool::Pnt( ((egadsNode *)enodes[0]->blind)->node );
    gp_Vec V(gp_Dir(axis[3], axis[4], axis[5]));
    gp_Pnt O(axis[0], axis[1], axis[2]);
    O.Translate(V.Dot(gp_Vec(O,P)) * V);
    if (O.IsEqual(P,Precision::Confusion())) {
      /* make a degenerated edge */
      stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                             tdata, 1, enodes, NULL, &edges[i]);
    } else {
      /* make a circle */
      stat = EG_rotate_circle(enodes[0], axis, data);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                             data, &circle);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(stack, circle);
      if (stat != EGADS_SUCCESS) goto cleanup;

#if 0
      {
        egObject *loop, *body, *model, *edge;
        tdata[1] = 2*PI;
        stat = EG_makeTopology(context, circle, EDGE, ONENODE,
                               tdata, 1, enodes, NULL, &edge);
        tdata[1] = angle;

        stat = EG_makeTopology(context, NULL, LOOP, CLOSED,
                               NULL, 1, &edge, esens, &loop);

        stat = EG_makeTopology(context, NULL, BODY, WIREBODY,
                               NULL, 1, &loop, NULL, &body);

        stat = EG_makeTopology(context, NULL, MODEL, 0,
                               NULL, 1, &body, NULL, &model);

        char filename[42];
        sprintf(filename, "circle_%d.egads", i);
        EG_saveModel(model, filename);
      }
#endif

      /* make the edge on the circle */
      if (enodes[0] == enodes[1]) {
        stat = EG_makeTopology(context, circle, EDGE, ONENODE,
                               tdata, 1, enodes, NULL, &edges[i]);
      } else {
        stat = EG_makeTopology(context, circle, EDGE, TWONODE,
                               tdata, 2, enodes, NULL, &edges[i]);
      }
    }
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, edges[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  if (loop0->mtype == OPEN) {
    /* create the last line/edge for open loop */
    egadsEdge *pedge0 = (egadsEdge *) ploop0->edges[i-1]->blind;
    egadsEdge *pedge1 = (egadsEdge *) ploop1->edges[i-1]->blind;

    if (ploop0->senses[i-1] == SFORWARD) {
      enodes[0] = pedge0->nodes[1];
      enodes[1] = pedge1->nodes[1];
    } else {
      enodes[0] = pedge0->nodes[0];
      enodes[1] = pedge1->nodes[0];
    }

    /* degeneracy check from BRepSweep_Rotation::MakeEmptyDirectingEdge */
    gp_Pnt P = BRep_Tool::Pnt( ((egadsNode *)enodes[0]->blind)->node );
    gp_Vec V(gp_Dir(axis[3], axis[4], axis[5]));
    gp_Pnt O(axis[0], axis[1], axis[2]);
    O.Translate(V.Dot(gp_Vec(O,P)) * V);
    if (O.IsEqual(P,Precision::Confusion())) {
      /* make a degenerated edge */
      stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                             tdata, 1, enodes, NULL, &edges[i]);
    } else {
      /* make a circle */
      stat = EG_rotate_circle(enodes[0], axis, data);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                             data, &circle);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(stack, circle);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* make the edge on the circle */
      if (enodes[0] == enodes[1]) {
        stat = EG_makeTopology(context, circle, EDGE, ONENODE,
                               tdata, 1, enodes, NULL, &edges[i]);
      } else {
        stat = EG_makeTopology(context, circle, EDGE, TWONODE,
                               tdata, 2, enodes, NULL, &edges[i]);
      }
    }
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, edges[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
  } else {
    /* repeat edge for closed loop */
    edges[i] = edges[0];
  }

  /* construct loops/faces */
  for (i = 0; i < ploop0->nedges; i++) {

    egadsEdge *pedge0 = (egadsEdge *) ploop0->edges[i]->blind;

    eedges[0] = ploop0->edges[i];
    eedges[1] = ploop0->senses[i] == SFORWARD ? edges[i] : edges[i+1];
    eedges[2] = ploop1->edges[i];
    eedges[3] = ploop0->senses[i] == SFORWARD ? edges[i+1] : edges[i];

    /* check if the central axis must be removed */
    if ((eedges[1]->mtype == DEGENERATE) &&
        (eedges[3]->mtype == DEGENERATE)) {
      egadsEdge *pedge = (egadsEdge *) eedges[0]->blind;

      if (pedge->curve->mtype == LINE) {
        faces[i] = NULL;
        continue;
      }
    }

#if 0
      {
        egObject *loop, *body, *model;

        stat = EG_makeTopology(context, NULL, LOOP, OPEN,
                               NULL, 4, eedges, esens, &loop);

        stat = EG_makeTopology(context, NULL, BODY, WIREBODY,
                               NULL, 1, &loop, NULL, &body);

        stat = EG_makeTopology(context, NULL, MODEL, 0,
                               NULL, 1, &body, NULL, &model);

        char filename[42];
        sprintf(filename, "loop_%d.egads", i);
        EG_saveModel(model, filename);
      }
#endif
    /* create P-curves */
    data[0] = tdata[0]; data[1] = 0.;    /* u == UMIN */
    data[2] = 0.;       data[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                           &eedges[4+0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, eedges[4+0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    data[0] = 0.;       data[1] = pedge0->trange[0];  /* v == VMIN */
    data[2] = 1.;       data[3] = 0.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                           &eedges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, eedges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    data[0] = tdata[1]; data[1] = 0.;    /* u == UMAX */
    data[2] = 0.;       data[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                           &eedges[4+2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, eedges[4+2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    data[0] = 0.;       data[1] = pedge0->trange[1];  /* v == VMAX */
    data[2] = 1.;       data[3] = 0.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                           &eedges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, eedges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the loop */
    stat = EG_makeTopology(context, surfs[i], LOOP, CLOSED,
                           NULL, 4, eedges, esens, &eloop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, eloop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the face */
    stat = EG_makeTopology(context, surfs[i], FACE, ssens[i],
                           NULL, 1, &eloop, lsens, &eface);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, eface);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* copy generating edge attributes to the face */
    EG_attributeDup(ploop0->edges[i], eface);

    /* set the output */
    faces[i] = eface;

#if 0
      {
        egObject *body, *model;

        stat = EG_makeTopology(context, NULL, BODY, FACEBODY,
                               NULL, 1, &eface, NULL, &body);

        stat = EG_makeTopology(context, NULL, MODEL, 0,
                               NULL, 1, &body, NULL, &model);

        char filename[42];
        sprintf(filename, "face_%d.egads", i);
        EG_saveModel(model, filename);
      }
#endif
  }

cleanup:
  EG_free(edges);
  EG_free(surfs);
  EG_free(ssens);

  return stat;
}


int
EG_rotate_copy(objStack *stack, egObject *src, egObject *exform,
               egObject ***loops0_out, egObject ***loops1_out, egObject **efaces)
{
  int      stat, i, j, nloops, found=0;
  egObject *context, *copy, *surf1;
  egObject **loops0=NULL, **loops1=NULL, **loops1_edges=NULL;

  context  = EG_context(src);

  stat = EG_copyObject(src, exform, &copy);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(stack, copy);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (src->oclass == FACE) {
    egadsFace *pface0 = (egadsFace *) src->blind;
    egadsFace *pface1 = (egadsFace *) copy->blind;

    loops1 = (egObject **) EG_alloc(pface1->nloops*sizeof(egObject *));

    for (i = 0; i < pface0->nloops; i++) {

      egadsLoop *ploop0 = (egadsLoop *) pface0->loops[i]->blind;
      egadsLoop *ploop1 = (egadsLoop *) pface1->loops[i]->blind;

      loops1_edges = (egObject **) EG_alloc(ploop1->nedges*sizeof(egObject *));

      /* find equivalent edges */
      for (j = 0; j < ploop0->nedges; j++) {

        loops1_edges[j] = ploop1->edges[j];

        if (EG_isEquivalent(ploop0->edges[j], loops1_edges[j]) == EGADS_SUCCESS) {
          loops1_edges[j] = ploop0->edges[j];
          found = 1;
        }
      }

      /* make the loop */
      surf1 = pface1->surface->mtype == PLANE ? NULL : pface1->surface;
      stat = EG_makeTopology(context, surf1, LOOP, CLOSED,
                             NULL, ploop1->nedges, loops1_edges, ploop1->senses, &loops1[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(stack, loops1[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;

      EG_free(loops1_edges); loops1_edges = NULL;
    }

    /* make the face */
    stat = EG_makeTopology(context, pface1->surface, FACE, copy->mtype,
                           NULL, pface1->nloops, loops1, pface1->senses, &copy);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(stack, copy);
    if (stat != EGADS_SUCCESS) goto cleanup;


    EG_free(loops0); loops0 = NULL;
    EG_free(loops1); loops1 = NULL;

    *loops0_out = pface0->loops;
    *loops1_out = pface1->loops;
    efaces[0] = src;
    efaces[1] = copy;
  }

cleanup:
  EG_free(loops0);
  EG_free(loops1);
  EG_free(loops1_edges);

  return stat;
}

} // no-name namespace


int
EG_rotate(const egObject *src, double angle, const double *axis,
          egObject **result)
{
  int      outLevel, stat, i, mtype, nfaces, nloops, offset;
  double   mat[12];
  egObject *context, *ref, **faces, *efaces[2];
  egObject *shell, *body, *exform, *ecopy;
  egObject **loops0=NULL, **loops1=NULL;
  objStack stack;

  nfaces = 2;

  *result = NULL;
  if (src == NULL)               return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (src->blind == NULL)        return EGADS_NODATA;
  if (EG_sameThread(src))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);

  angle *= PI/180.0;
  if ((angle < -2*PI) && (angle > 2*PI)) angle = 2*PI;
  if ((fabs(angle-2*PI) < Precision::Angular()) ||
      (fabs(angle+2*PI) < Precision::Angular())) angle = 2*PI;

  if (fabs(angle) < Precision::Angular()) {
    if (outLevel > 0)
      printf(" EGADS Error: Angle (%le) below tolerance (%le) (EG_rotate)!\n",
             angle, Precision::Angular());
    return EGADS_GEOMERR;
  }

  /* create stack for gracefully cleaning up objects */
  stat = EG_stackInit(&stack);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* get the rotation matrix */
  EG_rotationMatrix(angle, axis, mat);

  if ((fabs(angle-2*PI) < Precision::Angular()) ||
      (fabs(angle+2*PI) < Precision::Angular())) {
    exform = NULL;
  } else {
    stat = EG_makeTransform(context, mat, &exform);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, exform);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  mtype = SOLIDBODY;
  if  (src->oclass == BODY) {
    if ((src->mtype == WIREBODY) || (src->mtype == FACEBODY)) {
      egadsBody *pbody = (egadsBody *) src->blind;
      nloops = pbody->loops.map.Extent();
      loops0 = pbody->loops.objs;

      if ((fabs(angle-2*PI) < Precision::Angular()) ||
          (fabs(angle+2*PI) < Precision::Angular())) {
        nfaces = 0;
        loops1 = loops0;
        if ((src->mtype == WIREBODY) && (loops0[0]->mtype == OPEN))
          mtype  = SHEETBODY;
      } else {
        if (src->mtype == WIREBODY) {
          mtype  = SHEETBODY;
          nfaces = 0;
          stat   = EG_copyObject(loops0[0], exform, &ecopy);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat   = EG_stackPush(&stack, ecopy);
          if (stat != EGADS_SUCCESS) goto cleanup;

          loops1 = &ecopy;
        } else {
          efaces[0] = pbody->faces.objs[0];
          stat = EG_copyObject(efaces[0], exform, &efaces[1]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, efaces[1]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          egadsFace *pface = (egadsFace *) efaces[1]->blind;
          loops1 = pface->loops;
        }
      }
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Body src must be Wire or Face (EG_rotate)!\n");
      return EGADS_NOTTOPO;
    }
  } else if (src->oclass == LOOP) {
    nloops = 1;
    loops0 = const_cast<egObject**>(&src);
    mtype  = SHEETBODY;
    nfaces = 0;

    if ((fabs(angle-2*PI) < Precision::Angular()) ||
        (fabs(angle+2*PI) < Precision::Angular())) {
      nfaces = 0;
      loops1 = loops0;
      if (loops0[0]->mtype == OPEN)
        mtype  = SHEETBODY;
    } else {
      stat = EG_copyObject(src, exform, &ecopy);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, ecopy);
      if (stat != EGADS_SUCCESS) goto cleanup;

      loops1 = &ecopy;
    }
  } else if (src->oclass == FACE) {

    egadsFace *pface = (egadsFace *) src->blind;
    nloops = pface->nloops;
    loops0 = pface->loops;

    if ((fabs(angle-2*PI) < Precision::Angular()) ||
        (fabs(angle+2*PI) < Precision::Angular())) {
      nfaces = 0;
      loops1 = loops0;
    } else {
//      efaces[0] = const_cast<egObject*>(src);
//      stat = EG_copyObject(efaces[0], exform, &efaces[1]);
//      if (stat != EGADS_SUCCESS) goto cleanup;
//      stat = EG_stackPush(&stack, efaces[1]);
//      if (stat != EGADS_SUCCESS) goto cleanup;
//
//      pface  = (egadsFace *) efaces[1]->blind;
//      loops1 = pface->loops;

      stat = EG_rotate_copy(&stack, const_cast<egObject*>(src), exform, &loops0,
                            &loops1, efaces);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }
  } else {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid src type (EG_rotate)!\n");
    return EGADS_NOTTOPO;
  }


  for (i = 0; i < nloops; i++) {
    egadsLoop *ploop = (egadsLoop *) loops0[i]->blind;
    nfaces += ploop->nedges;
  }
  faces = (egObject **) EG_alloc(nfaces*sizeof(egObject *));

  offset = 0;
  for (i = 0; i < nloops; i++) {
    stat = EG_rotate_loop(&stack, loops0[i], loops1[i],
                           angle, axis, faces + offset);
    if (stat != EGADS_SUCCESS) goto cleanup;

    offset += ((egadsLoop *) loops0[i]->blind)->nedges;
  }

  if (!((fabs(angle-2*PI) < Precision::Angular()) ||
        (fabs(angle+2*PI) < Precision::Angular()))) {
    if (mtype == SOLIDBODY) {
      faces[offset+0] = efaces[0];
      faces[offset+1] = efaces[1];
    }
  }

  /* remove possible NULL faces */
  offset = nfaces;
  nfaces = 0;
  for (i = 0; i < offset; i++) {
    if (faces[i] == NULL) {
      mtype = SOLIDBODY;
      continue;
    }
    faces[nfaces++] = faces[i];
  }

#if 1
  for (i = 0; i < nfaces; i++) {
    egObject *body, *model;

    stat = EG_makeTopology(context, NULL, BODY, FACEBODY,
                           NULL, 1, &faces[i], NULL, &body);

    stat = EG_makeTopology(context, NULL, MODEL, 0,
                           NULL, 1, &body, NULL, &model);

    char filename[42];
    sprintf(filename, "face_%d.egads", i);
    EG_saveModel(model, filename);
  }
#endif

  /* put all the faces together */
  stat = EG_makeTopology(context, NULL, SHELL, mtype == SOLIDBODY ? CLOSED : OPEN,
                         NULL, nfaces, faces, NULL, &shell);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeTopo Shell = %d (EG_rotate)!\n", stat);
    goto cleanup;
  }
  stat = EG_stackPush(&stack, shell);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (shell->mtype == CLOSED) mtype = SOLIDBODY;

  body = NULL;
  stat = EG_makeTopology(context, NULL, BODY, mtype, NULL, 1,
                         &shell, NULL, &body);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeTopo Body = %d (EG_rotate)!\n", stat);
    if (body != NULL) EG_deleteObject(body);
    goto cleanup;
  }


  stat = EGADS_SUCCESS;
  *result = body;

cleanup:
  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    int i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (EG_extrude)!\n", i);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);
  EG_free(faces);

  return stat;
}
#endif


int
EG_offsetEdge(const egObject *src, const egObject *edge, const egObject *face,
              double offset, egObject **result)
{
  int          stat, outLevel, nerr;
  egObject     *context, *obj;
  TopoDS_Shape newFace, chamferFace, newEdge, newShape;

  *result = NULL;
  if (src == NULL)               return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(src))        return EGADS_CNTXTHRD;
  if (src->oclass != BODY)       return EGADS_NOTBODY;
  if (src->mtype != SOLIDBODY)   return EGADS_NOTTOPO;
  if (src->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);

  egadsBody *pbody = (egadsBody *) src->blind;
  BRepCheck_Analyzer check(pbody->shape);
  if (!check.IsValid()) {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid Input Body (EG_offsetEdge)!\n");
    return EGADS_TOPOERR;
  }

  if (face->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Face Object is not an EGO (EG_offsetEdge)!\n");
    return EGADS_NOTOBJ;
  }
  if (face->oclass != FACE) {
    if (outLevel > 0)
      printf(" EGADS Error: Object is not FACE (EG_offsetEdge)!\n");
    return EGADS_NOTTOPO;
  }
  if (face->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face Object has no data (EG_offsetEdge)!\n");
    return EGADS_NODATA;
  }
  egadsFace *pface = (egadsFace *) face->blind;
  if (pbody->faces.map.FindIndex(pface->face) == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Face is NOT in Body (EG_offsetEdge)!\n");
    return EGADS_NOTBODY;
  }

  if (edge->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge Object is not an EGO (EG_offsetEdge)!\n");
    return EGADS_NOTOBJ;
  }
  if (edge->oclass != EDGE) {
    if (outLevel > 0)
      printf(" EGADS Error: Object is not EDGE (EG_offsetEdge)!\n");
    return EGADS_NOTTOPO;
  }
  if (edge->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge Object has no data (EG_offsetEdge)!\n");
    return EGADS_NODATA;
  }
  if (edge->mtype == DEGENERATE) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge Object is Degenerate (EG_offsetEdge)!\n");
    return EGADS_DEGEN;
  }
  egadsEdge *pedge = (egadsEdge *) edge->blind;
  if (pbody->edges.map.FindIndex(pedge->edge) == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge is NOT in Body (EG_offsetEdge)!\n");
    return EGADS_NOTBODY;
  }

  // chamfer the body
  BRepFilletAPI_MakeChamfer chamfer(pbody->shape);
  chamfer.Add(offset, offset, pedge->edge, pface->face);

  try {
    chamfer.Build();
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Error: Chamfer Exception (EG_offsetEdge)!\n");
    printf("              %s\n", e.GetMessageString());
    return EGADS_GEOMERR;
  }
  catch (...) {
    printf(" EGADS Error: Chamfer Exception (EG_offsetEdge)!\n");
    return EGADS_GEOMERR;
  }

  if (!chamfer.IsDone()) {
    if (outLevel > 0)
      printf(" EGADS Error: Can't Chamfer (EG_offsetEdge)!\n");
    return EGADS_GEOMERR;
  }
  if (chamfer.NbContours() != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Number of Countour != 1 (EG_offsetEdge)!\n");
    return EGADS_CONSTERR;
  }

  newShape = chamfer.Shape();
  BRepCheck_Analyzer fCheck(newShape);
  if (!fCheck.IsValid()) {
    Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
    sfs->Perform();
    TopoDS_Shape fixedShape = sfs->Shape();
    if (fixedShape.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Error: Chamfered Body is invalid (EG_offsetEdge)!\n");
      return EGADS_GEOMERR;
    } else {
      BRepCheck_Analyzer sfCheck(fixedShape);
      if (!sfCheck.IsValid()) {
        if (outLevel > 0)
          printf(" EGADS Error: Fixed Body is invalid (EG_offsetEdge)!\n");
        return EGADS_GEOMERR;
      } else {
        newShape = fixedShape;
      }
    }
  }

  TopTools_IndexedMapOfShape emap, fmap, efmap;
  TopExp::MapShapes(newShape,    TopAbs_EDGE,  emap);
  TopExp::MapShapes(newShape,    TopAbs_FACE,  fmap);
  TopExp::MapShapes(pface->face, TopAbs_EDGE, efmap);

  const TopTools_ListOfShape& listFaces = chamfer.Modified(pface->face);
  if (listFaces.Extent() > 0) {
    /* modified Faces */
    TopTools_ListIteratorOfListOfShape it(listFaces);
    for (; it.More(); it.Next()) {
      TopoDS_Face modface = TopoDS::Face(it.Value());
      int index = fmap.FindIndex(modface);
      if (outLevel > 1)
        printf(" Old Face       = %d   New Face = %d [%d]\n",
               pbody->faces.map.FindIndex(pface->face), index, fmap.Extent());
      if (index == 0) continue;
      newFace = modface;
    }
  }

  for (int i = 1; i < efmap.Extent(); i++) {
    const TopTools_ListOfShape& listShapes = chamfer.Generated(efmap(i));
    if (listFaces.Extent() > 0) {
      TopTools_ListIteratorOfListOfShape it(listShapes);
      for (; it.More(); it.Next()) {
        TopoDS_Shape qFace = it.Value();
        if (qFace.ShapeType() != TopAbs_FACE) continue;
        if (outLevel > 1)
          printf(" Chamfered Face = %d\n", fmap.FindIndex(qFace));
        if (fmap.FindIndex(qFace) == 0) continue;
        chamferFace = qFace;
      }
    }
  }
  if (newFace.IsNull() || chamferFace.IsNull()) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot find Faces (EG_offsetEdge)!\n");
    return EGADS_CONSTERR;
  }

  // get the new Edge
  TopTools_IndexedMapOfShape cmap, nmap;
  TopExp::MapShapes(chamferFace, TopAbs_EDGE, cmap);
  TopExp::MapShapes(    newFace, TopAbs_EDGE, nmap);
  for (int i = 1; i <= nmap.Extent(); i++) {
    int index = cmap.FindIndex(nmap(i));
    if (index == 0) continue;
    if (outLevel > 1)
      printf(" Shared Edge    = %d\n", emap.FindIndex(nmap(i)));
    newEdge = nmap(i);
  }
  if (newEdge.IsNull()) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot find New Edge (EG_offsetEdge)!\n");
    return EGADS_CONSTERR;
  }

  BRepFeat_SplitShape Split(pbody->shape);
  try {
    Split.Add(TopoDS::Edge(newEdge), pface->face);
    Split.Build();
    if (!Split.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Can't Split (EG_offsetEdge)!\n");
      return EGADS_GEOMERR;
    }
    newShape = Split.Shape();
    BRepCheck_Analyzer sCheck(newShape);
    if (!sCheck.IsValid()) {
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
      sfs->Perform();
      TopoDS_Shape fixedShape = sfs->Shape();
      if (fixedShape.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Error: Resultant Body is invalid (EG_offsetEdge)!\n");
        return EGADS_GEOMERR;
      } else {
        BRepCheck_Analyzer sfCheck(fixedShape);
        if (!sfCheck.IsValid()) {
          printf(" EGADS Error: Fixed Body is invalid (EG_offsetEdge)!\n");
          return EGADS_GEOMERR;
        } else {
          newShape = fixedShape;
        }
      }
    }
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Warning: Split Construction Error (EG_offsetEdge)!\n");
    printf("                %s\n", e.GetMessageString());
    return EGADS_CONSTERR;
  }
  catch (...) {
    printf(" EGADS Warning: Split Construction Error (EG_offsetEdge)!\n");
    return EGADS_CONSTERR;
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_offsetEdge)!\n");
    return stat;
  }
  egadsBody *pbods   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = src->mtype;
  pbods->nodes.objs  = NULL;
  pbods->edges.objs  = NULL;
  pbods->loops.objs  = NULL;
  pbods->faces.objs  = NULL;
  pbods->shells.objs = NULL;
  pbods->senses      = NULL;
  pbods->shape       = newShape;
  pbods->bbox.filled = 0;
  pbods->massFill    = 0;
  obj->blind         = pbods;
  stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbods;
    return stat;
  }

  // map the Attributes
  EG_attriBodyDup(src, obj);
  const TopTools_ListOfShape& splitFaces = Split.Modified(pface->face);
  if (splitFaces.Extent() > 0) {
    /* modified Faces */
    TopTools_ListIteratorOfListOfShape it(splitFaces);
    for (; it.More(); it.Next()) {
      TopoDS_Face modface = TopoDS::Face(it.Value());
      int index = pbods->faces.map.FindIndex(modface);
      if (index > 0) EG_attributeDup(face, pbods->faces.objs[index-1]);
    }
  }
  // mapping does not always work here -- so
  EG_attrFixup(outLevel, obj, src);

  EG_referenceObject(obj, context);
  *result = obj;
  return EGADS_SUCCESS;
}


static int
EG_hollowFace(const egObject *object, int nedge, const egObject **edges,
              double offset, int joined, egObject **result, int **faceMap)
{
  int            i, j, k, m, stat, outLevel, iper, nfaces, nedx;
  double         d, range[4], uv[2], data[18], dir[3], *x0, *x1;
  egObject       *solid, **edgex, **faces, *hollow;
  const int      *ints;
  const double   *reals;
  const char     *str;
  const egObject **xfaces;

  *result  = NULL;
  outLevel = EG_outLevel(object);

  /* find direction to extrude */
  stat = EG_getRange(object, range, &iper);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getRange = %d (EG_hollowFace)!\n", stat);
    return stat;
  }
  uv[0] = 0.5*(range[0] + range[1]);
  uv[1] = 0.5*(range[2] + range[3]);
  stat  = EG_evaluate(object, uv, data);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_evaluate = %d (EG_hollowFace)!\n", stat);
    return stat;
  }
  x0 = &data[3];
  x1 = &data[6];
  CROSS(dir, x0, x1);
  d  = -sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2])*object->mtype;
  if (d == 0.0) {
    if (outLevel > 0)
      printf(" EGADS Error: Degenerate direction (EG_hollowFace)!\n");
    return EGADS_DEGEN;
  }
  dir[0] /= d;
  dir[1] /= d;
  dir[2] /= d;
  d       = 1.5*fabs(offset);

  /* extrude the Face to make a Solid */
  stat    = EG_extrude(object, d, dir, &solid);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_extrude = %d (EG_hollowFace)!\n", stat);
    return stat;
  }

  /* get the new Faces and make our new list */
  stat = EG_getBodyTopos(solid, NULL, FACE, &nfaces, &faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getBodyTopos S = %d (EG_hollowFace)!\n", stat);
    return stat;
  }

  /* get the Edges (now Faces in the Solid) */
  xfaces = (const egObject **) EG_alloc((nedge+1)*sizeof(egObject *));
  if (xfaces == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc xface (EG_hollowFace)!\n");
    EG_free(faces);
    return EGADS_MALLOC;
  }

  for (m = i = 0; i < nfaces; i++) {
    stat = EG_getBodyTopos(solid, faces[i], EDGE, &nedx, &edgex);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getBodyTopos E = %d (EG_hollowFace)!\n", stat);
      EG_free(xfaces);
      EG_free(faces);
      return stat;
    }
    for (j = 0; j < nedx; j++)
      for (k = 0; k < nedge; k++)
        if (EG_isEquivalent(edges[k], edgex[j]) == EGADS_SUCCESS) {
          if (outLevel > 1)
            printf(" Face %d/%d: obj %d isEquivalent Edge\n", i+1, nfaces, k+1);
          xfaces[m] = faces[i];
          m++;
          goto eloop;
        }
eloop:
    EG_free(edgex);
  }
  k = nfaces - 2;
  if (m == 0) {
    xfaces[m] = faces[k];
    m++;
  }

  /* mark our Face */
  stat = EG_attributeAdd(faces[k], ".hollow", ATTRSTRING, 7, NULL, NULL,
                         "target");
  EG_free(faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_attributeAdd = %d (EG_hollowFace)!\n", stat);
    EG_free(xfaces);
    return stat;
  }
  if (m != nedge+1) {
    if (outLevel > 0)
      printf(" EGADS Error: nFaces = %d (%d) (EG_hollowFace)!\n", m, nedge+1);
    EG_free(xfaces);
    return EGADS_INDEXERR;
  }

  /* do the solid-based hollow */
  stat = EG_hollowBody(solid, nedge+1, xfaces, offset, joined, &hollow, NULL);
  EG_deleteObject(solid);
  EG_free(xfaces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_hollowBody = %d (EG_hollowFace)!\n", stat);
    return stat;
  }

  /* pull Face from the new hollowed body */
  stat = EG_getBodyTopos(hollow, NULL, FACE, &nfaces, &faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getBodyTopos H = %d (EG_hollowFace)!\n", stat);
    EG_deleteObject(hollow);
    return stat;
  }
  for (i = 0; i < nfaces; i++) {
    stat = EG_attributeRet(faces[i], ".hollow", &j, &m, &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) continue;
    if (outLevel > 1)
      printf(" EG_hollowFace: Attribute Found in Face %d of %d!\n", i, nfaces);
    if (*result != NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Multiple Resultant Faces (EG_hollowFace)!\n");
      EG_free(faces);
      EG_deleteObject(*result);
      EG_deleteObject(hollow);
      *result = NULL;
      return EGADS_CONSTERR;
    }
    stat = EG_copyObject(faces[i], NULL, result);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_copyObject = %d (EG_hollowFace)!\n", stat);
    }
  }
  EG_free(faces);
  EG_deleteObject(hollow);
  if (*result == NULL) return EGADS_NOTFOUND;

  /* remove the attribute used to find the Face */
  EG_attributeDel(*result, ".hollow");

  return EGADS_SUCCESS;
}


int
EG_hollowBody(const egObject *src, int nface, const egObject **faces,
              double offset, int joined, egObject **result, int **faceMap)
{
  int      i, outLevel, stat, nerr, fullAttrs, *fmap;
  double   tol = Precision::Confusion();
  egObject *context, *obj, *face, *edge, *node;

  *result = NULL;
  if (faceMap != NULL) *faceMap = NULL;
  if  (src == NULL)                return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if  (EG_sameThread(src))         return EGADS_CNTXTHRD;
  if  (src->oclass == FACE)
    return EG_hollowFace(src, nface, faces, offset, joined, result, faceMap);
  if  (src->oclass != BODY)        return EGADS_NOTBODY;
  if ((src->mtype != SOLIDBODY) && (src->mtype != SHEETBODY) &&
      (src->mtype != FACEBODY))    return EGADS_NOTTOPO;
  if  (src->blind == NULL)         return EGADS_NODATA;
  outLevel  = EG_outLevel(src);
  context   = EG_context(src);
  fullAttrs = EG_fullAttrs(src);

  if (src->mtype != SOLIDBODY) {
    if (nface != 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Removing Faces on a nonSolid (EG_hollowBody)!\n");
      return EGADS_NODATA;
    }
  }
  if (nface < 0) {
    if (outLevel > 0)
      printf(" EGADS Error: No Faces (EG_hollowBody)!\n");
    return EGADS_NODATA;
  }
  if ((faces == NULL) && (nface > 0)) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face Pointer (EG_hollowBody)!\n");
    return EGADS_NULLOBJ;
  }
  TopTools_ListOfShape aList;
  egadsBody *pbody = (egadsBody *) src->blind;
  for (i = 0; i < nface; i++) {
    if (faces[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Face Object %d (EG_hollowBody)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (faces[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Object %d is not an EGO (EG_hollowBody)!\n",
               i+1);
      return EGADS_NOTOBJ;
    }
    if (faces[i]->oclass != FACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not FACE (EG_hollowBody)!\n", i+1);
      return EGADS_NOTTOPO;
    }
    if (faces[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Object %d has no data (EG_hollowBody)!\n",
               i+1);
      return EGADS_NODATA;
    }
    egadsFace *pface = (egadsFace *) faces[i]->blind;
    if (pbody->faces.map.FindIndex(pface->face) == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d is NOT in Body (EG_hollowBody)!\n",
               i+1);
      return EGADS_NOTBODY;
    }
    aList.Append(pface->face);
    if (tol < BRep_Tool::Tolerance(pface->face))
      tol = BRep_Tool::Tolerance(pface->face);
  }
  BRepCheck_Analyzer check(pbody->shape);
  if (!check.IsValid()) {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid Input Body (EG_hollowBody)!\n");
    return EGADS_TOPOERR;
  }

  GeomAbs_JoinType join = GeomAbs_Arc;
  if (joined == 1) join = GeomAbs_Intersection;

  if ((nface == 0) && (src->mtype != SHEETBODY)) {
    // offset the body
    TopoDS_Shape newShape;
    try {
#if CASVER >= 720
      BRepOffsetAPI_MakeOffsetShape offshape;
      offshape.PerformByJoin(pbody->shape, offset, tol,
                             BRepOffset_Skin, Standard_False,
                             Standard_False, join);
#else
      BRepOffsetAPI_MakeOffsetShape offshape(pbody->shape, offset, tol,
                                             BRepOffset_Skin, Standard_False,
                                             Standard_False, join);
#endif
      newShape = offshape.Shape();
      BRepCheck_Analyzer fCheck(newShape);
      if (!fCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Error: Offset Body is invalid (EG_hollowBody)!\n");
          return EGADS_CONSTERR;
        } else {
          BRepCheck_Analyzer sfCheck(fixedShape);
          if (!sfCheck.IsValid()) {
            printf(" EGADS Error: Offset Fixed Body is invalid (EG_hollowBody)!\n");
            return EGADS_CONSTERR;
          } else {
            newShape = fixedShape;
          }
        }
      }

      if (src->mtype == FACEBODY) {
        if (newShape.ShapeType() == TopAbs_SHELL) {
          int nface = 0;
          TopExp_Explorer Exp;
          for (Exp.Init(newShape, TopAbs_FACE); Exp.More(); Exp.Next()) nface++;
          if (nface == 1) {
            Exp.Init(newShape, TopAbs_FACE);
            newShape = Exp.Current();
          }
        }
        if (newShape.ShapeType() != TopAbs_FACE) {
          if (outLevel > 0)
            printf(" EGADS Error: Offset Result Not a Face (EG_hollowBody)!\n");
          return EGADS_CONSTERR;
        }
      } else {
        int nsolid = 0;
        if (newShape.ShapeType() == TopAbs_COMPOUND) {
          TopExp_Explorer Exp;
          for (Exp.Init(newShape, TopAbs_SOLID); Exp.More(); Exp.Next()) nsolid++;
          if (nsolid == 1) {
            Exp.Init(newShape, TopAbs_SOLID);
            newShape = Exp.Current();
          }
        }
        if (newShape.ShapeType() != TopAbs_SOLID) {
          if (outLevel > 0) {
            printf(" EGADS Error: Offset Result Not a Solid (EG_hollowBody)!\n");
            if (newShape.ShapeType() == TopAbs_COMPOUND) {
              printf("              Type = COMPOUND; nSolid = %d\n", nsolid);
              if (nsolid == 0) {
                TopExp_Explorer Exp;
                for (Exp.Init(newShape, TopAbs_SHELL); Exp.More(); Exp.Next())
                  nsolid++;
                printf("                               nShell = %d\n", nsolid);
              }
              if (nsolid == 0) {
                TopExp_Explorer Exp;
                for (Exp.Init(newShape, TopAbs_FACE); Exp.More(); Exp.Next())
                  nsolid++;
                printf("                               nFace  = %d\n", nsolid);
              }
            }
            if (newShape.ShapeType() == TopAbs_COMPSOLID)
              printf("              Type = COMPSOLID!\n");
            if (newShape.ShapeType() == TopAbs_SHELL)
              printf("              Type = SHELL!\n");
            if (newShape.ShapeType() == TopAbs_FACE)
              printf("              Type = FACE!\n");
            if (newShape.ShapeType() == TopAbs_WIRE)
              printf("              Type = WIRE!\n");
            if (newShape.ShapeType() == TopAbs_EDGE)
              printf("              Type = EDGE!\n");
            if (newShape.ShapeType() == TopAbs_VERTEX)
              printf("              Type = VERTEX!\n");
            if (newShape.ShapeType() == TopAbs_SHAPE)
              printf("              Type = SHAPE!\n");
          }
          return EGADS_CONSTERR;
        }
      }

      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Cannot Make Body object (EG_hollowBody)!\n");
        return stat;
      }
      egadsBody *pbods   = new egadsBody;
      obj->oclass        = BODY;
      obj->mtype         = src->mtype;
      pbods->nodes.objs  = NULL;
      pbods->edges.objs  = NULL;
      pbods->loops.objs  = NULL;
      pbods->faces.objs  = NULL;
      pbods->shells.objs = NULL;
      pbods->senses      = NULL;
      pbods->shape       = newShape;
      pbods->bbox.filled = 0;
      pbods->massFill    = 0;
      obj->blind         = pbods;
      stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
      if (stat != EGADS_SUCCESS) {
        delete pbods;
        return stat;
      }

      // Face map info: source, modified, generated
      if (faceMap != NULL) {
        fmap = (int *) EG_alloc(2*pbods->faces.map.Extent()*sizeof(int));
        if (fmap != NULL) {
          for (i = 0; i < 2*pbods->faces.map.Extent(); i++) fmap[i] = 0;
          for (i = 0; i <   pbody->faces.map.Extent(); i++) {
            face = pbody->faces.objs[i];
            egadsFace *pface = (egadsFace *) face->blind;
            int index = pbods->faces.map.FindIndex(pface->face);
            if (index > 0) {
              fmap[2*index-2] = FACEDUP;
              fmap[2*index-1] = i+1;
            }
          }
          for (i = 0; i < pbody->faces.map.Extent(); i++) {
            face = pbody->faces.objs[i];
            egadsFace *pface = (egadsFace *) face->blind;
            const TopTools_ListOfShape& listFaces = offshape.Modified(pface->face);
            if (listFaces.Extent() > 0) {
              /* modified Faces */
              TopTools_ListIteratorOfListOfShape it(listFaces);
              for (; it.More(); it.Next()) {
                TopoDS_Shape qFace = it.Value();
                if (qFace.ShapeType() != TopAbs_FACE) continue;
                TopoDS_Face genface = TopoDS::Face(qFace);
                int index = pbods->faces.map.FindIndex(genface);
                if (index > 0) {
                  fmap[2*index-2] = FACECUT;
                  fmap[2*index-1] = i+1;
                }
              }
            }
          }
          for (i = 0; i < pbody->faces.map.Extent(); i++) {
            face = pbody->faces.objs[i];
            egadsFace *pface = (egadsFace *) face->blind;
            const TopTools_ListOfShape& listFaces = offshape.Generated(pface->face);
            if (listFaces.Extent() > 0) {
              /* generated Faces from Faces */
              TopTools_ListIteratorOfListOfShape it(listFaces);
              for (; it.More(); it.Next()) {
                TopoDS_Shape qFace = it.Value();
                if (qFace.ShapeType() != TopAbs_FACE) continue;
                TopoDS_Face genface = TopoDS::Face(qFace);
                int index = pbods->faces.map.FindIndex(genface);
                if (index > 0) {
                  fmap[2*index-2] = FACEOFF;
                  fmap[2*index-1] = i+1;
                }
              }
            }
          }
          for (i = 0; i < pbody->edges.map.Extent(); i++) {
            edge = pbody->edges.objs[i];
            egadsEdge *pedge = (egadsEdge *) edge->blind;
            const TopTools_ListOfShape& genFaces = offshape.Generated(pedge->edge);
            if (genFaces.Extent() > 0) {
              /* generated Faces from Edges */
              TopTools_ListIteratorOfListOfShape it(genFaces);
              for (; it.More(); it.Next()) {
                TopoDS_Shape qFace = it.Value();
                if (qFace.ShapeType() != TopAbs_FACE) continue;
                TopoDS_Face genface = TopoDS::Face(qFace);
                int index = pbods->faces.map.FindIndex(genface);
                if (index > 0) {
                  fmap[2*index-2] = EDGEOFF;
                  fmap[2*index-1] = i+1;
                }
              }
            }
          }
          for (i = 0; i < pbody->nodes.map.Extent(); i++) {
            node = pbody->nodes.objs[i];
            egadsNode *pnode = (egadsNode *) node->blind;
            const TopTools_ListOfShape& genFaces = offshape.Generated(pnode->node);
            if (genFaces.Extent() > 0) {
              /* generated Faces from Nodes */
              TopTools_ListIteratorOfListOfShape it(genFaces);
              for (; it.More(); it.Next()) {
                TopoDS_Shape qFace = it.Value();
                if (qFace.ShapeType() != TopAbs_FACE) continue;
                TopoDS_Face genface = TopoDS::Face(qFace);
                int index = pbods->faces.map.FindIndex(genface);
                if (index > 0) {
                  fmap[2*index-2] = NODEOFF;
                  fmap[2*index-1] = i+1;
                }
              }
            }
          }

          /* check the fmap for all non-zero*/
          int count = 0;
          for (i = 0; i < 2*pbods->faces.map.Extent(); i++) {
            if (fmap[i] == 0) count++;
          }
          if (count != 0) {
            printf(" EGADS Error: Invalid faceMap (EG_hollowBody)!\n");
            EG_deleteObject(obj);
            EG_free(fmap);
            return EGADS_CONSTERR;
          }

          *faceMap = fmap;
        }
      }
    }
    catch (const Standard_Failure& e) {
      printf(" EGADS Error: MakeOffsetShape Exception (EG_hollowBody)!\n");
      printf("              %s\n", e.GetMessageString());
      return EGADS_CONSTERR;
    }
    catch (...) {
      printf(" EGADS Error: MakeOffsetShape Exception (EG_hollowBody)!\n");
      return EGADS_CONSTERR;
    }

    EG_referenceObject(obj, context);
    *result = obj;
    return EGADS_SUCCESS;
  }

  // hollow the solid body / "solidize" the sheet body
  try {
#if CASVER >= 720
    BRepOffsetAPI_MakeThickSolid hollow;
    if (src->mtype == SHEETBODY) {
      hollow.MakeThickSolidBySimple(pbody->shape, -offset);
    } else {
      hollow.MakeThickSolidByJoin(pbody->shape, aList, -offset, tol,
                                  BRepOffset_Skin, Standard_False,
                                  Standard_False, join);
    }
#else
    if (src->mtype == SHEETBODY) return EGADS_NOTTOPO;
    BRepOffsetAPI_MakeThickSolid hollow(pbody->shape, aList, -offset, tol,
                                        BRepOffset_Skin, Standard_False,
                                        Standard_False, join);
#endif
    hollow.Build();
    TopoDS_Shape newShape = hollow.Shape();
    if (src->mtype == SHEETBODY) {
      BRepGProp    BProps;
      GProp_GProps VProps;
      BProps.VolumeProperties(newShape, VProps);
      if (VProps.Mass() < 0.0) newShape.Reverse();
    }
    BRepCheck_Analyzer fCheck(newShape);
    if (!fCheck.IsValid()) {
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
      sfs->Perform();
      TopoDS_Shape fixedShape = sfs->Shape();
      if (fixedShape.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Error: Hollowed Body is invalid (EG_hollowBody)!\n");
        return EGADS_CONSTERR;
      } else {
        BRepCheck_Analyzer sfCheck(fixedShape);
        if (!sfCheck.IsValid()) {
          printf(" EGADS Error: Fixed Body is invalid (EG_hollowBody)!\n");
          return EGADS_CONSTERR;
        } else {
          newShape = fixedShape;
        }
      }
    }

    // make sure we have the correct result!
    if (newShape.ShapeType() == TopAbs_COMPOUND) {
      int nsolid = 0;
      TopExp_Explorer Exp;
      for (Exp.Init(newShape, TopAbs_SOLID); Exp.More(); Exp.Next()) nsolid++;
      if (nsolid == 1) {
        Exp.Init(newShape, TopAbs_SOLID);
        newShape = Exp.Current();
      }
    }
    if (newShape.ShapeType() != TopAbs_SOLID) {
      if (outLevel > 0)
        printf(" EGADS Error: Result Not a Solid (EG_hollowBody)!\n");
      return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Body object (EG_hollowBody)!\n");
      return stat;
    }
    egadsBody *pbods   = new egadsBody;
    obj->oclass        = BODY;
    obj->mtype         = SOLIDBODY;
    pbods->nodes.objs  = NULL;
    pbods->edges.objs  = NULL;
    pbods->loops.objs  = NULL;
    pbods->faces.objs  = NULL;
    pbods->shells.objs = NULL;
    pbods->senses      = NULL;
    pbods->shape       = newShape;
    pbods->bbox.filled = 0;
    pbods->massFill    = 0;
    obj->blind         = pbods;
    stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
    if (stat != EGADS_SUCCESS) {
      delete pbods;
      return stat;
    }

    // map the Attributes
    EG_attriBodyDup(src, obj);
    if (fullAttrs != 0) {
      for (i = 0; i < pbody->edges.map.Extent(); i++) {
        ego edge = pbody->edges.objs[i];
        egadsEdge *pedge = (egadsEdge *) edge->blind;
        const TopTools_ListOfShape& listEdges = hollow.Modified(pedge->edge);
        if (listEdges.Extent() > 0) {
          /* modified Edges */
          TopTools_ListIteratorOfListOfShape it(listEdges);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qEdge = it.Value();
            if (qEdge.ShapeType() != TopAbs_EDGE) continue;
            TopoDS_Edge genedge = TopoDS::Edge(qEdge);
            int index = pbods->edges.map.FindIndex(genedge);
            if (index > 0) EG_attributeDup(edge, pbods->edges.objs[index-1]);
          }
        }
      }
    }

    for (i = 0; i < pbody->faces.map.Extent(); i++) {
      ego face = pbody->faces.objs[i];
      egadsFace *pface = (egadsFace *) face->blind;
      const TopTools_ListOfShape& listFaces = hollow.Modified(pface->face);
      if (listFaces.Extent() > 0) {
        /* modified Faces */
        TopTools_ListIteratorOfListOfShape it(listFaces);
        for (; it.More(); it.Next()) {
          TopoDS_Shape qFace = it.Value();
          if (qFace.ShapeType() != TopAbs_FACE) continue;
          TopoDS_Face genface = TopoDS::Face(qFace);
          int index = pbods->faces.map.FindIndex(genface);
          if (index > 0) EG_attributeDup(face, pbods->faces.objs[index-1]);
        }
      }
    }

    // Face map info: source, modified, generated
    if (faceMap != NULL) {
      fmap = (int *) EG_alloc(2*pbods->faces.map.Extent()*sizeof(int));
      if (fmap != NULL) {
        for (i = 0; i < 2*pbods->faces.map.Extent(); i++) fmap[i] = 0;
        for (i = 0; i <   pbody->faces.map.Extent(); i++) {
          face = pbody->faces.objs[i];
          egadsFace *pface = (egadsFace *) face->blind;
          int index = pbods->faces.map.FindIndex(pface->face);
          if (index > 0) {
            fmap[2*index-2] = FACEDUP;
            fmap[2*index-1] = i+1;
          }
        }
        for (i = 0; i < pbody->faces.map.Extent(); i++) {
          face = pbody->faces.objs[i];
          egadsFace *pface = (egadsFace *) face->blind;
          const TopTools_ListOfShape& listFaces = hollow.Modified(pface->face);
          if (listFaces.Extent() > 0) {
            /* modified Faces */
            TopTools_ListIteratorOfListOfShape it(listFaces);
            for (; it.More(); it.Next()) {
              TopoDS_Shape qFace = it.Value();
              if (qFace.ShapeType() != TopAbs_FACE) continue;
              TopoDS_Face genface = TopoDS::Face(qFace);
              int index = pbods->faces.map.FindIndex(genface);
              if (index > 0) {
                fmap[2*index-2] = FACECUT;
                fmap[2*index-1] = i+1;
              }
            }
          }
        }
        for (i = 0; i < pbody->faces.map.Extent(); i++) {
          face = pbody->faces.objs[i];
          egadsFace *pface = (egadsFace *) face->blind;
          const TopTools_ListOfShape& listFaces = hollow.Generated(pface->face);
          if (listFaces.Extent() > 0) {
            /* generated Faces from Faces */
            TopTools_ListIteratorOfListOfShape it(listFaces);
            for (; it.More(); it.Next()) {
              TopoDS_Shape qFace = it.Value();
              if (qFace.ShapeType() != TopAbs_FACE) continue;
              TopoDS_Face genface = TopoDS::Face(qFace);
              int index = pbods->faces.map.FindIndex(genface);
              if (index > 0) {
                fmap[2*index-2] = FACEOFF;
                fmap[2*index-1] = i+1;
              }
            }
          }
        }
        for (i = 0; i < pbody->edges.map.Extent(); i++) {
          edge = pbody->edges.objs[i];
          egadsEdge *pedge = (egadsEdge *) edge->blind;
          const TopTools_ListOfShape& genFaces = hollow.Generated(pedge->edge);
          if (genFaces.Extent() > 0) {
            /* generated Faces from Edges */
            TopTools_ListIteratorOfListOfShape it(genFaces);
            for (; it.More(); it.Next()) {
              TopoDS_Shape qFace = it.Value();
              if (qFace.ShapeType() != TopAbs_FACE) continue;
              TopoDS_Face genface = TopoDS::Face(qFace);
              int index = pbods->faces.map.FindIndex(genface);
              if (index > 0) {
                fmap[2*index-2] = EDGEOFF;
                fmap[2*index-1] = i+1;
              }
            }
          }
        }
        for (i = 0; i < pbody->nodes.map.Extent(); i++) {
          node = pbody->nodes.objs[i];
          egadsNode *pnode = (egadsNode *) node->blind;
          const TopTools_ListOfShape& genFaces = hollow.Generated(pnode->node);
          if (genFaces.Extent() > 0) {
            /* generated Faces from Nodes */
            TopTools_ListIteratorOfListOfShape it(genFaces);
            for (; it.More(); it.Next()) {
              TopoDS_Shape qFace = it.Value();
              if (qFace.ShapeType() != TopAbs_FACE) continue;
              TopoDS_Face genface = TopoDS::Face(qFace);
              int index = pbods->faces.map.FindIndex(genface);
              if (index > 0) {
                fmap[2*index-2] = NODEOFF;
                fmap[2*index-1] = i+1;
              }
            }
          }
        }

        /* check the fmap for all non-zero*/
        int count = 0;
        for (i = 0; i < 2*pbods->faces.map.Extent(); i++) {
          if (fmap[i] == 0) count++;
        }
        if (count != 0) {
          printf(" EGADS Error: Invalid faceMap (EG_hollowBody)!\n");
          EG_deleteObject(obj);
          EG_free(fmap);
          return EGADS_CONSTERR;
        }

        *faceMap = fmap;
      }
    }

  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Error: MakeThickSolid Exception (EG_hollowBody)!\n");
    printf("              %s\n", e.GetMessageString());
    return EGADS_CONSTERR;
  }
  catch (...) {
    printf(" EGADS Error: MakeThickSolid Exception (EG_hollowBody)!\n");
    return EGADS_CONSTERR;
  }

  EG_referenceObject(obj, context);
  *result = obj;
  return EGADS_SUCCESS;
}


int
EG_sweep(const egObject *src, const egObject *spine, int mode,
               egObject **result)
{
  int          outLevel, stat, mtype, nerr, nseg = 1;
  egObject     *context, *obj, *edge;
  TopoDS_Shape shape, newShape;
  TopoDS_Wire  wire;

  *result = NULL;
  if (src == NULL)               return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (src->blind == NULL)        return EGADS_NODATA;
  if (EG_sameThread(src))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(src);
  context  = EG_context(src);

  if (spine == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL edge Reference (EG_sweep)!\n");
    return EGADS_NULLOBJ;
  }
  if (spine->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: edge not an EGO (EG_sweep)!\n");
    return EGADS_NOTOBJ;
  }
  if ((spine->oclass != EDGE) && (spine->oclass != LOOP) &&
      (spine->oclass != BODY)) {
    if (outLevel > 0)
      printf(" EGADS Error: Not an Edge/LOOP/WIREBODY (EG_sweep)!\n");
    return EGADS_NOTTOPO;
  }
  if (context != EG_context(spine)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_sweep)!\n");
    return EGADS_MIXCNTX;
  }
  if ((mode < GeomFill_IsCorrectedFrenet) ||
      (mode > GeomFill_IsDiscreteTrihedron)) {
    if (outLevel > 0)
      printf(" EGADS Error: Mode = %d [%d - %d] (EG_sweep)!\n", mode,
             GeomFill_IsCorrectedFrenet, GeomFill_IsDiscreteTrihedron);
    return EGADS_INDEXERR;
  }

  // get spine shape
  if (spine->oclass == EDGE) {
    BRepBuilderAPI_MakeWire MW;
    egadsEdge *pedge = (egadsEdge *) spine->blind;
    TopoDS_Edge edge = pedge->edge;
    edge.Orientation(TopAbs_FORWARD);
    MW.Add(edge);
    if (MW.Error()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem adding Edge (EG_sweep)!\n");
      return EGADS_NODATA;
    }
    if (!MW.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem with Loop (EG_sweep)!\n");
      return EGADS_NODATA;
    }
    wire = MW.Wire();
  } else if (spine->oclass == LOOP) {
    egadsLoop *ploos = (egadsLoop *) spine->blind;
    wire = ploos->loop;
    nseg = ploos->nedges;
  } else {
    if (spine->mtype != WIREBODY) {
      if (outLevel > 0)
        printf(" EGADS Error: Body not a WIREBODY (EG_sweep)!\n");
      return EGADS_NOTTOPO;
    }
    egadsBody *pbods = (egadsBody *) spine->blind;
    if (pbods->loops.map.Extent() != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: WIREBODY with %d Loops (EG_sweep)!\n",
               pbods->loops.map.Extent());
      return EGADS_NOTTOPO;
    }
    egadsLoop *ploos = (egadsLoop *) pbods->loops.objs[0]->blind;
    wire = ploos->loop;
    nseg = ploos->nedges;
  }

  // get sweeped shape
  mtype = SOLIDBODY;
  if  (src->oclass == BODY) {
    if ((src->mtype == WIREBODY) || (src->mtype == FACEBODY)) {
      egadsBody *pbody = (egadsBody *) src->blind;
      shape = pbody->shape;
      if (src->mtype == WIREBODY) mtype = SHEETBODY;
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Body src must be Wire or Face (EG_sweep)!\n");
      return EGADS_NOTTOPO;
    }
  } else if (src->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) src->blind;
    shape = ploop->loop;
    mtype = SHEETBODY;
  } else if (src->oclass == FACE) {
    egadsFace *pface = (egadsFace *) src->blind;
    shape = pface->face;
  } else {
    if (outLevel > 0)
      printf(" EGADS Error: Invalid src type (EG_sweep)!\n");
    return EGADS_NOTTOPO;
  }

  // make the sweeped shape
  try {
    GeomFill_Trihedron aMode = (GeomFill_Trihedron) mode;
    Standard_Boolean   fac1  = Standard_False;
    if (nseg != 1)     fac1  = Standard_True;
    BRepOffsetAPI_MakePipe pipe(wire, shape, aMode, fac1);
    newShape = pipe.Shape();
    if (mtype == SOLIDBODY) {
      if (newShape.ShapeType() != TopAbs_SOLID) {
        if (outLevel > 0)
          printf(" EGADS Error: Sweep Result Not a Solid (EG_sweep)!\n");
        return EGADS_CONSTERR;
      }
    } else {
      if (newShape.ShapeType() != TopAbs_SHELL) {
        if (outLevel > 0)
          printf(" EGADS Error: Sweep Result Not a Shell (EG_sweep)!\n");
        return EGADS_CONSTERR;
      }
    }

    // are we OK?
    BRepCheck_Analyzer check(newShape);
    if (!check.IsValid()) {
      // try to fix the fault
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
      sfs->Perform();
      TopoDS_Shape fixedShape = sfs->Shape();
      if (fixedShape.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Info: Invalid Shape w/ NULL Fix (EG_sweep)!\n");
        return EGADS_CONSTERR;
      }
      BRepCheck_Analyzer fxCheck(fixedShape);
      if (!fxCheck.IsValid()) {
        if (outLevel > 0)
          printf(" EGADS Info: Result is invalid (EG_sweep)!\n");
        return EGADS_CONSTERR;
      }
      newShape = fixedShape;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Body object (EG_sweep)!\n");
      return stat;
    }
    egadsBody *pbods   = new egadsBody;
    obj->oclass        = BODY;
    obj->mtype         = mtype;
    pbods->nodes.objs  = NULL;
    pbods->edges.objs  = NULL;
    pbods->loops.objs  = NULL;
    pbods->faces.objs  = NULL;
    pbods->shells.objs = NULL;
    pbods->senses      = NULL;
    pbods->shape       = newShape;
    pbods->bbox.filled = 0;
    pbods->massFill    = 0;
    obj->blind         = pbods;
    stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
    if (stat != EGADS_SUCCESS) {
      delete pbods;
      return stat;
    }

    // Note: end-caps not handled

    /* copy attributes from the generating Edges */
    if (src->oclass == BODY) {
      egadsBody *pbody = (egadsBody *) src->blind;
      for (int i = 0; i < pbody->edges.map.Extent(); i++) {
        edge = pbody->edges.objs[i];
        egadsEdge *pedge = (egadsEdge *) edge->blind;
#if CASVER >= 700 // Not sure this is the right version...
        const TopTools_ListOfShape& genFaces = pipe.Generated(pedge->edge);
#else
        const TopTools_ListOfShape& genFaces = static_cast<BRepPrimAPI_MakeSweep&>(pipe).Generated(pedge->edge);
#endif
        if (genFaces.Extent() > 0) {
          /* generated Faces from Edges */
          TopTools_ListIteratorOfListOfShape it(genFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) EG_attributeDup(edge, pbods->faces.objs[index-1]);
          }
        }
      }
    } else if (src->oclass == LOOP) {
      egadsLoop *ploop = (egadsLoop *) src->blind;
      for (int i = 0; i < ploop->nedges; i++) {
        edge = ploop->edges[i];
        egadsEdge *pedge = (egadsEdge *) edge->blind;
#if CASVER >= 700 // Not sure this is the right version...
        const TopTools_ListOfShape& genFaces = pipe.Generated(pedge->edge);
#else
        const TopTools_ListOfShape& genFaces = static_cast<BRepPrimAPI_MakeSweep&>(pipe).Generated(pedge->edge);
#endif
        if (genFaces.Extent() > 0) {
          /* generated Faces from Edges */
          TopTools_ListIteratorOfListOfShape it(genFaces);
          for (; it.More(); it.Next()) {
            TopoDS_Shape qFace = it.Value();
            if (qFace.ShapeType() != TopAbs_FACE) continue;
            TopoDS_Face genface = TopoDS::Face(qFace);
            int index = pbods->faces.map.FindIndex(genface);
            if (index > 0) EG_attributeDup(edge, pbods->faces.objs[index-1]);
          }
        }
      }
    } else {
      egadsFace *pface = (egadsFace *) src->blind;
      for (int j = 0; j < pface->nloops; j++) {
        egadsLoop *ploop = (egadsLoop *) pface->loops[j]->blind;
        for (int i = 0; i < ploop->nedges; i++) {
          edge = ploop->edges[i];
          egadsEdge *pedge = (egadsEdge *) edge->blind;
#if CASVER >= 700 // Not sure this is the right version...
          const TopTools_ListOfShape& genFaces = pipe.Generated(pedge->edge);
#else
          const TopTools_ListOfShape& genFaces = static_cast<BRepPrimAPI_MakeSweep&>(pipe).Generated(pedge->edge);
#endif
          if (genFaces.Extent() > 0) {
            /* generated Faces from Edges */
            TopTools_ListIteratorOfListOfShape it(genFaces);
            for (; it.More(); it.Next()) {
              TopoDS_Shape qFace = it.Value();
              if (qFace.ShapeType() != TopAbs_FACE) continue;
              TopoDS_Face genface = TopoDS::Face(qFace);
              int index = pbods->faces.map.FindIndex(genface);
              if (index > 0) EG_attributeDup(edge, pbods->faces.objs[index-1]);
            }
          }
        }
      }
    }

    EG_referenceObject(obj, context);
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Error: BRepOffsetAPI_MakePipe Exception (EG_sweep)!\n");
    printf("              %s\n", e.GetMessageString());
    return EGADS_GEOMERR;
  }
  catch (...) {
    printf(" EGADS Error: BRepOffsetAPI_MakePipe Exception (EG_sweep)!\n");
    return EGADS_GEOMERR;
  }

  *result = obj;
  return EGADS_SUCCESS;
}


int
EG_loft(int nsec, const egObject **secs, int opt, egObject **result)
{
  int      i, outLevel, stat, nerr;
  egObject *context, *obj;

  *result = NULL;
  if (nsec <= 1)                     return EGADS_EMPTY;
  if (secs == NULL)                  return EGADS_NULLOBJ;
  if (secs[0] == NULL)               return EGADS_NULLOBJ;
  if (secs[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(secs[0]))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(secs[0]);
  context  = EG_context(secs[0]);

  for (i = 0; i < nsec; i++) {
    if (secs[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Section Object %d (EG_loft)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (secs[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is not an EGO (EG_loft)!\n", i+1);
      return EGADS_NOTOBJ;
    }
    if (secs[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d has no data (EG_loft)!\n", i+1);
      return EGADS_NODATA;
    }
    if (EG_context(secs[i]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d Context Mismatch (EG_loft)!\n", i+1);
      return EGADS_MIXCNTX;
    }
    if (secs[i]->oclass == NODE) {
      if ((i != 0) && (i != nsec-1)) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Node and not Bound (EG_loft)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
    } else if (secs[i]->oclass == BODY) {
      if (secs[i]->mtype != WIREBODY) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Not a WireBody (EG_loft)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
    } else if (secs[i]->oclass != LOOP) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is Not a Loop (EG_loft)!\n", i+1);
      return EGADS_NOTTOPO;
    }
  }

  Standard_Boolean isSolid = ((opt&1) == 1);
  Standard_Boolean isRuled = ((opt&2) == 2);
  TopoDS_Shape     newShape;
  try {
    BRepOffsetAPI_ThruSections loft(isSolid, isRuled);
    for (i = 0; i < nsec; i++)
      if (secs[i]->oclass == NODE) {
        egadsNode *pnode = (egadsNode *) secs[i]->blind;
        loft.AddVertex(pnode->node);
      } else if (secs[i]->oclass == BODY) {
        egadsBody *pbody = (egadsBody *) secs[i]->blind;
        TopoDS_Shape bshape = pbody->shape;
        TopoDS_Wire wire = TopoDS::Wire(bshape);
        loft.AddWire(wire);
      } else {
        egadsLoop *ploop = (egadsLoop *) secs[i]->blind;
        loft.AddWire(ploop->loop);
      }
    loft.Build();
    if (!loft.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Can't Loft (EG_loft)!\n");
      return EGADS_GEOMERR;
    }
    newShape = loft.Shape();
  }
  catch (const Standard_Failure& e) {
    printf(" EGADS Error: ThruSections Exception (EG_loft)!\n");
    printf("              %s\n", e.GetMessageString());
    return EGADS_GEOMERR;
  }
  catch (...) {
    printf(" EGADS Error: ThruSections Exception (EG_loft)!\n");
    return EGADS_GEOMERR;
  }

  // are we OK?
  BRepCheck_Analyzer check(newShape);
  if (!check.IsValid()) {
    // try to fix the fault
    Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(newShape);
    sfs->Perform();
    TopoDS_Shape fixedShape = sfs->Shape();
    if (fixedShape.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Info: Invalid Shape w/ NULL Fix (EG_loft)!\n");
      return EGADS_CONSTERR;
    }
    BRepCheck_Analyzer fxCheck(fixedShape);
    if (!fxCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Result is invalid (EG_loft)!\n");
      return EGADS_CONSTERR;
    }
    newShape = fixedShape;
  }
  if (newShape.ShapeType() == TopAbs_COMPOUND) {
    if (outLevel > 0)
      printf(" EGADS Error: Result is Compound (EG_loft)!\n");
    return EGADS_CONSTERR;
  }
  if (isSolid) {
    if (newShape.ShapeType() != TopAbs_SOLID) {
      if (outLevel > 0)
        printf(" EGADS Error: Result Not a Solid (EG_loft)!\n");
      return EGADS_CONSTERR;
    }
  } else {
    if (newShape.ShapeType() != TopAbs_SHELL) {
      if (outLevel > 0)
        printf(" EGADS Error: Result Not a Shell (EG_loft)!\n");
      return EGADS_CONSTERR;
    }
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_loft)!\n");
    return stat;
  }
  egadsBody *pbods   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = SHEETBODY;
  pbods->nodes.objs  = NULL;
  pbods->edges.objs  = NULL;
  pbods->loops.objs  = NULL;
  pbods->faces.objs  = NULL;
  pbods->shells.objs = NULL;
  pbods->senses      = NULL;
  pbods->shape       = newShape;
  pbods->bbox.filled = 0;
  pbods->massFill    = 0;
  obj->blind         = pbods;
  if (isSolid) obj->mtype = SOLIDBODY;
  stat = EG_traverseBody(context, 0, obj, obj, pbods, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbods;
    return stat;
  }

  // do we want to do anything with attributes?

  EG_referenceObject(obj, context);
  *result = obj;
  return EGADS_SUCCESS;
}
