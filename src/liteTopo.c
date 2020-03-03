/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Lite Topology Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "liteClasses.h"


#define PARAMACC         1.0e-4         /* parameter accuracy */
#define CROSS(a,b,c)       a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                           a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                           a[2] = (b[0]*c[1]) - (b[1]*c[0])


  extern int EG_evaluate( const egObject *geom, const double *param,
                          double *result );
  extern int EG_invEvaluate( const egObject *geom, const double *xyz,
                             double *param, double *result );
  extern int EG_invEvaLimits( const egObject *geom, /*@null@*/ const double *lim,
                             const double *xyz, double *param, double *result );
  extern int EG_getRange( const egObject *geom, double *range, int *periodic );
  extern int EG_inFaceX( const egObject *face, const double *uva,
                         /*@null@*/ double *pt, /*@null@*/ double *uvx );


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
    liteNode *pnode;
    pnode = (liteNode *) topo->blind;
    if ((limits != NULL) && (pnode != NULL)) {
      limits[0] = pnode->xyz[0];
      limits[1] = pnode->xyz[1];
      limits[2] = pnode->xyz[2];
    }
  } else if (topo->oclass == EDGE) {
    liteEdge *pedge;
    pedge = (liteEdge *) topo->blind;
    if (pedge != NULL) {
      *geom      = pedge->curve;
      *nChildren = 1;
      if (topo->mtype == TWONODE) *nChildren = 2;
      *children = pedge->nodes;
      if (limits != NULL) {
        limits[0] = pedge->trange[0];
        limits[1] = pedge->trange[1];
      }
    }
  } else if (topo->oclass == LOOP) {
    liteLoop *ploop;
    ploop = (liteLoop *) topo->blind;
    if (ploop != NULL) {
      *geom      = ploop->surface;
      *nChildren = ploop->nedges;
      *children  = ploop->edges;
      *senses    = ploop->senses;
    }
  } else if (topo->oclass == FACE) {
    liteFace *pface;
    pface = (liteFace *) topo->blind;
    if (pface != NULL) {
      *geom       = pface->surface;
      *nChildren  = pface->nloops;
      *children   = pface->loops;
      *senses     = pface->senses;
      if (limits != NULL) {
        limits[0] = pface->urange[0];
        limits[1] = pface->urange[1];
        limits[2] = pface->vrange[0];
        limits[3] = pface->vrange[1];
      }
    }
  } else if (topo->oclass == SHELL) {
    liteShell *pshell;
    pshell = (liteShell *) topo->blind;
    if (pshell != NULL) {
      *nChildren = pshell->nfaces;
      *children  = pshell->faces;
    }
  } else if (topo->oclass == BODY) {
    liteBody *pbody;
    pbody = (liteBody *) topo->blind;
    if (pbody != NULL)
      if (topo->mtype == WIREBODY) {
        *nChildren = pbody->loops.nobjs;
        *children  = pbody->loops.objs;
      } else if (topo->mtype == FACEBODY) {
        *nChildren = pbody->faces.nobjs;
        *children  = pbody->faces.objs;
      } else {
        *nChildren = pbody->shells.nobjs;
        *children  = pbody->shells.objs;
        if (topo->mtype == SOLIDBODY) *senses = pbody->senses;
      }
  } else {
    liteModel *pmodel;
    pmodel = (liteModel *) topo->blind;
    if (pmodel != NULL) {
      *nChildren = pmodel->nbody;
      *children  = pmodel->bodies;
    }
  }
  
  return EGADS_SUCCESS;
}


static int
EG_contained(egObject *obj, egObject *src)
{
  int i, stat;
  
  if (src->oclass == EDGE) {
    liteEdge *pedge;
    pedge = (liteEdge *) src->blind;
    for (i = 0; i < 2; i++)
      if (pedge->nodes[i] == obj) return EGADS_SUCCESS;
  } else if (src->oclass == LOOP) {
    liteLoop *ploop;
    ploop = (liteLoop *) src->blind;
    if (obj->oclass == EDGE) {
      for (i = 0; i < ploop->nedges; i++)
        if (ploop->edges[i] == obj) return EGADS_SUCCESS;
    } else {
      for (i = 0; i < ploop->nedges; i++) {
        stat = EG_contained(obj, ploop->edges[i]);
        if (stat == EGADS_SUCCESS) return stat;
      }
    }
  } else if (src->oclass == FACE) {
    liteFace *pface;
    pface = (liteFace *) src->blind;
    if (obj->oclass == LOOP) {
      for (i = 0; i < pface->nloops; i++)
        if (pface->loops[i] == obj) return EGADS_SUCCESS;
    } else {
      for (i = 0; i < pface->nloops; i++) {
        stat = EG_contained(obj, pface->loops[i]);
        if (stat == EGADS_SUCCESS) return stat;
      }
    }
  } else if (src->oclass == SHELL) {
    liteShell *pshell;
    pshell = (liteShell *) src->blind;
    if (obj->oclass == FACE) {
      for (i = 0; i < pshell->nfaces; i++)
        if (pshell->faces[i] == obj) return EGADS_SUCCESS;
    } else {
      for (i = 0; i < pshell->nfaces; i++) {
        stat = EG_contained(obj, pshell->faces[i]);
        if (stat == EGADS_SUCCESS) return stat;
      }
    }
  }
  
  return EGADS_OUTSIDE;
}


int
EG_getTolerance(const egObject *topo, double *toler)
{
  if  (topo == NULL)               return EGADS_NULLOBJ;
  if  (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (topo->blind == NULL)        return EGADS_NODATA;
  if ((topo->oclass != NODE) && (topo->oclass != EDGE) &&
      (topo->oclass != FACE))      return EGADS_NOTTOPO;
  
  if (topo->oclass == NODE) {
    liteNode *pnode;
    pnode  = (liteNode *) topo->blind;
    *toler = pnode->tol;
  } else if (topo->oclass == EDGE) {
    liteEdge *pedge;
    pedge = (liteEdge *) topo->blind;
    *toler = pedge->tol;
  } else {
    liteFace *pface;
    pface = (liteFace *) topo->blind;
    *toler = pface->tol;
  }
  
  return EGADS_SUCCESS;
}


int
EG_tolerance(const egObject *topo, double *toler)
{
  int    i, stat;
  double tol;

  if  (topo == NULL)                                  return EGADS_NULLOBJ;
  if  (topo->magicnumber != MAGIC)                    return EGADS_NOTOBJ;
  if  (topo->blind == NULL)                           return EGADS_NODATA;
  if ((topo->oclass < NODE) || (topo->oclass > BODY)) return EGADS_NOTTOPO;
  
  *toler = 0.0;
  if (topo->oclass == NODE) {
    liteNode *pnode;
    pnode  = (liteNode *) topo->blind;
    *toler = pnode->tol;
  } else if (topo->oclass == EDGE) {
    liteEdge *pedge;
    pedge = (liteEdge *) topo->blind;
    *toler = pedge->tol;
    stat   = EG_tolerance(pedge->nodes[0], &tol);
    if (stat != EGADS_SUCCESS) return stat;
    if (tol > *toler) *toler = tol;
    stat   = EG_tolerance(pedge->nodes[1], &tol);
    if (stat != EGADS_SUCCESS) return stat;
    if (tol > *toler) *toler = tol;
  } else if (topo->oclass == LOOP) {
    liteLoop *ploop;
    ploop = (liteLoop *) topo->blind;
    for (i = 0; i < ploop->nedges; i++) {
      stat = EG_tolerance(ploop->edges[i], &tol);
      if (stat != EGADS_SUCCESS) return stat;
      if (tol > *toler) *toler = tol;
    }
  } else if (topo->oclass == FACE) {
    liteFace *pface;
    pface = (liteFace *) topo->blind;
    *toler = pface->tol;
    for (i = 0; i < pface->nloops; i++) {
      stat = EG_tolerance(pface->loops[i], &tol);
      if (stat != EGADS_SUCCESS) return stat;
      if (tol > *toler) *toler = tol;
    }
  } else if (topo->oclass == SHELL) {
    liteShell *pshell;
    pshell = (liteShell *) topo->blind;
    for (i = 0; i < pshell->nfaces; i++) {
      stat = EG_tolerance(pshell->faces[i], &tol);
      if (stat != EGADS_SUCCESS) return stat;
      if (tol > *toler) *toler = tol;
    }
  } else {
    liteBody *pbody;
    pbody = (liteBody *) topo->blind;
    for (i = 0; i < pbody->shells.nobjs; i++) {
      stat = EG_tolerance(pbody->shells.objs[i], &tol);
      if (stat != EGADS_SUCCESS) return stat;
      if (tol > *toler) *toler = tol;
    }
  }

  return EGADS_SUCCESS;
}


int
EG_getBodyTopos(const egObject *body, /*@null@*/ egObject *src,
                int oclass, int *ntopo, /*@null@*/ egObject ***topos)
{
  int      i, n;
  egObject **objs;
  liteBody *pbody;
  liteMap  map;

  *ntopo = 0;
  if (topos != NULL) *topos = NULL;
  if  (body == NULL)                       return EGADS_NULLOBJ;
  if  (body->magicnumber != MAGIC)         return EGADS_NOTOBJ;
  if  (body->oclass != BODY)               return EGADS_NOTBODY;
  if  (body->blind == NULL)                return EGADS_NODATA;
  if ((oclass < NODE) || (oclass > SHELL)) return EGADS_NOTTOPO;
  
  pbody = (liteBody *) body->blind;
  if (oclass == NODE) {
    map = pbody->nodes;
  } else if (oclass == EDGE) {
    map = pbody->edges;
  } else if (oclass == LOOP) {
    map = pbody->loops;
  } else if (oclass == FACE) {
    map = pbody->faces;
  } else {
    map = pbody->shells;
  }
  
  /* all body objects */
  if (src == NULL) {
    if (map.nobjs == 0) return EGADS_SUCCESS;
    if (topos == NULL) {
      *ntopo = map.nobjs;
      return EGADS_SUCCESS;
    }
    objs = (egObject **) EG_alloc(map.nobjs*sizeof(egObject *));
    if (objs == NULL) return EGADS_MALLOC;
    for (i = 0; i < map.nobjs; i++) objs[i] = map.objs[i];
    *ntopo = map.nobjs;
    *topos = objs;
    return EGADS_SUCCESS;
  }

  /* selective body objects */
  if  (src->magicnumber != MAGIC)                    return EGADS_NOTOBJ;
  if ((src->oclass < NODE) || (src->oclass > SHELL)) return EGADS_NOTTOPO;
  if  (src->oclass == oclass)                        return EGADS_TOPOERR;
  if  (src->blind == NULL)                           return EGADS_NODATA;
  
  /* look down the tree */
  if (src->oclass > oclass) {
    for (n = i = 0; i < map.nobjs; i++)
      if (EG_contained(map.objs[i], src) == EGADS_SUCCESS) n++;
    if (n == 0) return EGADS_SUCCESS;
    if (topos == NULL) {
      *ntopo = n;
      return EGADS_SUCCESS;
    }
    objs = (egObject **) EG_alloc(n*sizeof(egObject *));
    if (objs == NULL) return EGADS_MALLOC;
    for (n = i = 0; i < map.nobjs; i++)
      if (EG_contained(map.objs[i], src) == EGADS_SUCCESS) {
        objs[n] = map.objs[i];
        n++;
      }
    *ntopo = n;
    *topos = objs;
    return EGADS_SUCCESS;
  }
  
  /* look up the tree */
  for (n = i = 0; i < map.nobjs; i++)
    if (EG_contained(src, map.objs[i]) == EGADS_SUCCESS) n++;
  if (n == 0) return EGADS_SUCCESS;
  if (topos == NULL) {
    *ntopo = n;
    return EGADS_SUCCESS;
  }
  objs = (egObject **) EG_alloc(n*sizeof(egObject *));
  if (objs == NULL) return EGADS_MALLOC;
  for (n = i = 0; i < map.nobjs; i++)
    if (EG_contained(src, map.objs[i]) == EGADS_SUCCESS) {
      objs[n] = map.objs[i];
      n++;
    }
  *ntopo = n;
  *topos = objs;
  return EGADS_SUCCESS;
}


int
EG_indexBodyTopo(const egObject *body, const egObject *src)
{
  int      i;
  liteBody *pbody;
  
  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->oclass != BODY)       return EGADS_NOTBODY;
  if (body->blind == NULL)        return EGADS_NODATA;
  if (src->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if ((src->oclass < NODE) || (src->oclass > SHELL))
                                  return EGADS_NOTTOPO;
  if (src->blind == NULL)         return EGADS_NODATA;

  pbody = (liteBody *) body->blind;
  if (src->oclass == NODE) {
    for (i = 0; i < pbody->nodes.nobjs; i++)
      if (src == pbody->nodes.objs[i]) return i+1;
  } else if (src->oclass == EDGE) {
    for (i = 0; i < pbody->edges.nobjs; i++)
      if (src == pbody->edges.objs[i]) return i+1;
  } else if (src->oclass == LOOP) {
    for (i = 0; i < pbody->loops.nobjs; i++)
      if (src == pbody->loops.objs[i]) return i+1;
  } else if (src->oclass == FACE) {
    for (i = 0; i < pbody->faces.nobjs; i++)
      if (src == pbody->faces.objs[i]) return i+1;
  } else {
    for (i = 0; i < pbody->shells.nobjs; i++)
      if (src == pbody->shells.objs[i]) return i+1;
  }

  return EGADS_NOTFOUND;
}


int
EG_objectBodyTopo(const egObject *body, int oclass, int index, egObject **obj)
{
  liteBody *pbody;
  liteMap  map;
  
  if  (body == NULL)                       return EGADS_NULLOBJ;
  if  (body->magicnumber != MAGIC)         return EGADS_NOTOBJ;
  if  (body->oclass != BODY)               return EGADS_NOTBODY;
  if  (body->blind == NULL)                return EGADS_NODATA;
  if ((oclass < NODE) || (oclass > SHELL)) return EGADS_NOTTOPO;
  if  (index <= 0)                         return EGADS_INDEXERR;
  
  pbody = (liteBody *) body->blind;
  if (oclass == NODE) {
    map = pbody->nodes;
  } else if (oclass == EDGE) {
    map = pbody->edges;
  } else if (oclass == LOOP) {
    map = pbody->loops;
  } else if (oclass == FACE) {
    map = pbody->faces;
  } else {
    map = pbody->shells;
  }
  if (index > map.nobjs) return EGADS_INDEXERR;
  
  *obj = map.objs[index-1];
  return EGADS_SUCCESS;
}


int
EG_getBoundingBox(const egObject *topo, double *bbox)
{
  int i;
  
  bbox[0] = bbox[1] = bbox[2] = bbox[3] = bbox[4] = bbox[5] = 0.0;
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->oclass < NODE)        return EGADS_NOTTOPO;
  if (topo->blind == NULL)        return EGADS_NODATA;
  
  if (topo->oclass == NODE) {
    liteNode *pnode;
    pnode = (liteNode *) topo->blind;
    if (pnode != NULL) {
      bbox[0] = bbox[3] = pnode->xyz[0];
      bbox[1] = bbox[4] = pnode->xyz[1];
      bbox[2] = bbox[5] = pnode->xyz[2];
    }
  } else if (topo->oclass == EDGE) {
    liteEdge *pedge;
    pedge = (liteEdge *) topo->blind;
    if (pedge != NULL)
      for (i = 0; i < 6; i++) bbox[i] = pedge->bbox[i];
  } else if (topo->oclass == LOOP) {
    liteLoop *ploop;
    ploop = (liteLoop *) topo->blind;
    if (ploop != NULL)
      for (i = 0; i < 6; i++) bbox[i] = ploop->bbox[i];
  } else if (topo->oclass == FACE) {
    liteFace *pface;
    pface = (liteFace *) topo->blind;
    if (pface != NULL)
      for (i = 0; i < 6; i++) bbox[i] = pface->bbox[i];
  } else if (topo->oclass == SHELL) {
    liteShell *pshell;
    pshell = (liteShell *) topo->blind;
    if (pshell != NULL)
      for (i = 0; i < 6; i++) bbox[i] = pshell->bbox[i];
  } else if (topo->oclass == BODY) {
    liteBody *pbody;
    pbody = (liteBody *) topo->blind;
    if (pbody != NULL)
      for (i = 0; i < 6; i++) bbox[i] = pbody->bbox[i];
  } else {
    liteModel *pmodel;
    pmodel = (liteModel *) topo->blind;
    if (pmodel != NULL)
      for (i = 0; i < 6; i++) bbox[i] = pmodel->bbox[i];
  }
  
  return EGADS_SUCCESS;
}


int
EG_getEdgeUV(const egObject *face, const egObject *edge, int sense, double t,
             double *result)
{
  int      i, j, stat, found;
  double   data[9], xyz[3];
  egObject *surface, *loop, *pcurve;
  liteLoop *ploop;
  liteFace *pface;
  
  result[0] = result[1] = 0.0;
  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  if (edge == NULL)               return EGADS_NULLOBJ;
  if (edge->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (edge->oclass != EDGE)       return EGADS_NOTTOPO;
  if (edge->blind == NULL)        return EGADS_NODATA;
  if ((sense != -1) && (sense != 0) &&
      (sense !=  1))              return EGADS_RANGERR;
  pface   = (liteFace *) face->blind;
  surface = pface->surface;
  if (surface == NULL)            return EGADS_NULLOBJ;

  /* special code for Planar Faces -- no PCurves */
  if (surface->mtype == PLANE) {
    stat = EG_evaluate(edge, &t, data);
    if (stat != EGADS_SUCCESS) return stat;
    return EG_invEvaluate(pface->surface, data, result, xyz);
  }
  
  /* check if we can use sense of zero */
  if (sense == 0) {
    found = 0;
    for (i = 0; i < pface->nloops; i++) {
      if (pface->loops[i] == NULL) continue;
      ploop = (liteLoop *) pface->loops[i]->blind;
      if (ploop == NULL) continue;
      for (j = 0; j < ploop->nedges; j++)
        if (ploop->edges[j] == edge) found++;
    }
    if (found > 1) {
      printf(" EGADS Warning: Edge in Face twice & sense=0 (EG_getEdgeUV)!\n");
      return EGADS_TOPOERR;
    }
  }
    
  /* find Edge/Sense pair in Face */
  for (i = 0; i < pface->nloops; i++) {
    loop  = pface->loops[i];
    if (loop == NULL) continue;
    ploop = (liteLoop *) loop->blind;
    if (ploop == NULL) continue;
    for (j = 0; j < ploop->nedges; j++)
      if (ploop->edges[j] == edge) {
        if (sense != 0)
          if (sense != ploop->senses[j]) continue;
        pcurve = ploop->edges[j+ploop->nedges];
        stat   = EG_evaluate(pcurve, &t, data);
        if (stat == EGADS_SUCCESS) {
          result[0] = data[0];
          result[1] = data[1];
        }
        return stat;
      }
  }
  
  return EGADS_NOTFOUND;
}


int
EG_getEdgeUVs(const egObject *face, const egObject *edge, int sense, int nt,
              const double *t, double *uvs)
{
  int      i, j, k, stat;
  double   data[9], xyz[3];
  egObject *surface, *loop, *pcurve;
  liteLoop *ploop;
  liteFace *pface;
  
  for (k = 0; k < 2*nt; k++) uvs[k] = 0.0;
  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  if (edge == NULL)               return EGADS_NULLOBJ;
  if (edge->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (edge->oclass != EDGE)       return EGADS_NOTTOPO;
  if (edge->blind == NULL)        return EGADS_NODATA;
  if ((sense != -1) && (sense != 0) &&
      (sense !=  1))              return EGADS_RANGERR;
  pface   = (liteFace *) face->blind;
  surface = pface->surface;
  if (surface == NULL)            return EGADS_NULLOBJ;
  
  /* special code for Planar Faces -- no PCurves */
  if (surface->mtype == PLANE) {
    for (k = 0; k < nt; k++) {
      stat = EG_evaluate(edge, &t[k], data);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_invEvaluate(pface->surface, data, &uvs[2*k], xyz);
      if (stat != EGADS_SUCCESS) return stat;
    }
    return EGADS_SUCCESS;
  }
  
  /* find Edge/Sense pair in Face */
  for (i = 0; i < pface->nloops; i++) {
    loop  = pface->loops[i];
    if (loop == NULL) continue;
    ploop = (liteLoop *) loop->blind;
    if (ploop == NULL) continue;
    for (j = 0; j < ploop->nedges; j++)
      if (ploop->edges[j] == edge) {
        if (sense != 0)
          if (sense != ploop->senses[j]) continue;
        pcurve = ploop->edges[j+ploop->nedges];
        for (k = 0; k < nt; k++) {
          stat = EG_evaluate(pcurve, &t[k], data);
          if (stat != EGADS_SUCCESS) return stat;
          uvs[2*k  ] = data[0];
          uvs[2*k+1] = data[1];
        }
        return EGADS_SUCCESS;
      }
  }
  
  return EGADS_NOTFOUND;
}


int
EG_getEdgeUVeval(const egObject *face, const egObject *edge, int sense,
                 double t, double *result)
{
  int      i, j, stat;
  double   data[9], xyz[3], eval[18];
  egObject *surface, *loop, *pcurve;
  liteLoop *ploop;
  liteFace *pface;
  
  result[0] = result[1] = result[2] = result[3] = result[4] = result[5] = 0.0;
  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  if (edge == NULL)               return EGADS_NULLOBJ;
  if (edge->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (edge->oclass != EDGE)       return EGADS_NOTTOPO;
  if (edge->blind == NULL)        return EGADS_NODATA;
  if (edge->mtype == DEGENERATE)  return EGADS_DEGEN;
  if ((sense != -1) && (sense != 0) &&
      (sense !=  1))              return EGADS_RANGERR;
  pface   = (liteFace *) face->blind;
  surface = pface->surface;
  if (surface == NULL)            return EGADS_NULLOBJ;
  
  /* special code for Planar Faces -- no PCurves */
  if (surface->mtype == PLANE) {
    stat = EG_evaluate(edge, &t, data);
    if (stat != EGADS_SUCCESS) return stat;
    /* this may not be correct! */
    stat = EG_invEvaluate(pface->surface, data, result, xyz);
    if (stat != EGADS_SUCCESS) return stat;
    stat = EG_evaluate(face, result, eval);
    if (stat != EGADS_SUCCESS) return stat;
    result[2] = (data[3]*eval[3] + data[4]*eval[4] + data[5]*eval[5]) /
                (eval[3]*eval[3] + eval[4]*eval[4] + eval[5]*eval[5]);
    result[3] = (data[3]*eval[6] + data[4]*eval[7] + data[5]*eval[8]) /
                (eval[6]*eval[6] + eval[7]*eval[7] + eval[8]*eval[8]);
  }
  
  /* find Edge/Sense pair in Face */
  for (i = 0; i < pface->nloops; i++) {
    loop  = pface->loops[i];
    if (loop == NULL) continue;
    ploop = (liteLoop *) loop->blind;
    if (ploop == NULL) continue;
    for (j = 0; j < ploop->nedges; j++)
      if (ploop->edges[j] == edge) {
        if (sense != 0)
          if (sense != ploop->senses[j]) continue;
        pcurve = ploop->edges[j+ploop->nedges];
        return EG_evaluate(pcurve, &t, result);
      }
  }
  
  return EGADS_NOTFOUND;
}


int
EG_getBody(const egObject *obj, egObject **body)
{
  int       i;
  liteModel *pmodel;
  egObject  *topObj, *bod;
  
  *body = NULL;
  if (obj == NULL)                  return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (obj->blind == NULL)           return EGADS_NODATA;
  if ((obj->oclass < NODE) ||
      (obj->oclass > SHELL))        return EGADS_NOTTOPO;
  topObj = obj->topObj;
  if (topObj == NULL)               return EGADS_NULLOBJ;
  if (topObj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  
  if (topObj->oclass == BODY) {
    *body = topObj;
  } else if (topObj->oclass == MODEL) {
    pmodel = (liteModel *) topObj->blind;
    if (pmodel != NULL) {
      for (i = 0; i < pmodel->nbody; i++) {
        bod = pmodel->bodies[i];
        if (bod == NULL) continue;
        if (EG_indexBodyTopo(bod, obj) <= 0) continue;
        *body = bod;
        break;
      }
    }
  }
  
  return EGADS_SUCCESS;
}


int
EG_inFace(const egObject *face, const double *uv)
{
  return EG_inFaceX(face, uv, NULL, NULL);
}


int
EG_inTopology(const egObject *topo, const double *xyz)
{
  int       i, j, stat;
  double    d, dist, param[2], uv[2], coord[3], dir[3], norm[3], data[18];
  liteEdge  *pedge;
  liteFace  *pface;
  liteShell *pshell;
  liteBody  *pbody;
  egObject  *face;
  
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->blind == NULL)        return EGADS_NODATA;
  
  if (topo->oclass == EDGE) {
    pedge = (liteEdge *) topo->blind;
    stat  = EG_invEvaluate(pedge->curve, xyz, param, coord);
    if (stat != EGADS_SUCCESS) return stat;
    dist  = sqrt((xyz[0]-coord[0])*(xyz[0]-coord[0]) +
                 (xyz[1]-coord[1])*(xyz[1]-coord[1]) +
                 (xyz[2]-coord[2])*(xyz[2]-coord[2]));
    if (dist > pedge->tol) return EGADS_OUTSIDE;
    if ((param[0] < pedge->trange[0]) || (param[0] > pedge->trange[1]))
      return EGADS_OUTSIDE;
  } else if (topo->oclass == FACE) {
    pface = (liteFace *) topo->blind;
    stat  = EG_invEvaluate(pface->surface, xyz, param, coord);
    if (stat != EGADS_SUCCESS) return stat;
    dist  = sqrt((xyz[0]-coord[0])*(xyz[0]-coord[0]) +
                 (xyz[1]-coord[1])*(xyz[1]-coord[1]) +
                 (xyz[2]-coord[2])*(xyz[2]-coord[2]));
    if (dist > pface->tol) return EGADS_OUTSIDE;
    return EG_inFace(topo, param);
  } else if ((topo->oclass == SHELL) && (topo->mtype == CLOSED)) {
    dist   = 1.e308;
    face   = NULL;
    pshell = (liteShell *) topo->blind;
    for (i = 0; i < pshell->nfaces; i++) {
      pface = (liteFace *) pshell->faces[i]->blind;
      if (pface == NULL) continue;
      stat  = EG_invEvaluate(pshell->faces[i], xyz, param, coord);
      if (stat != EGADS_SUCCESS) return stat;
      d = sqrt((xyz[0]-coord[0])*(xyz[0]-coord[0]) +
               (xyz[1]-coord[1])*(xyz[1]-coord[1]) +
               (xyz[2]-coord[2])*(xyz[2]-coord[2]));
      if (d < dist) {
        if (d <= pface->tol) return EGADS_SUCCESS;
        dist  = d;
        face  = pshell->faces[i];
        uv[0] = param[0];
        uv[1] = param[1];
      }
    }
    if (face == NULL) return EGADS_OUTSIDE;
    stat = EG_evaluate(face, uv, data);
    if (stat != EGADS_SUCCESS) return stat;
    coord[0] = data[3];
    coord[1] = data[4];
    coord[2] = data[5];
    dir[0]   = data[6];
    dir[1]   = data[7];
    dir[2]   = data[8];
    CROSS(norm, coord, dir);
    if (face->mtype == SREVERSE) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
      norm[2] = -norm[2];
    }
    d        = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    if (d == 0.0) return EGADS_DEGEN;
    norm[0] /= d;
    norm[1] /= d;
    norm[2] /= d;
    dir[0]   = xyz[0] - data[0];
    dir[1]   = xyz[1] - data[1];
    dir[2]   = xyz[2] - data[2];
    d        = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    dir[0]  /= d;
    dir[1]  /= d;
    dir[2]  /= d;
    if (dir[0]*norm[0]+dir[1]*norm[1]+dir[2]*norm[2] > 0.0) return EGADS_OUTSIDE;
    return EGADS_SUCCESS;
    
  } else if ((topo->oclass == BODY) && (topo->mtype == SOLIDBODY)) {
    pbody = (liteBody *) topo->blind;
    dist  = 1.e308;
    face  = NULL;
    for (j = 0; j < pbody->shells.nobjs; j++) {
      pshell = (liteShell *) pbody->shells.objs[j]->blind;
      if (pshell == NULL) continue;
      for (i = 0; i < pshell->nfaces; i++) {
        pface = (liteFace *) pshell->faces[i]->blind;
        if (pface == NULL) continue;
        stat  = EG_invEvaluate(pshell->faces[i], xyz, param, coord);
        if (stat == EGADS_DEGEN)   continue;
        if (stat != EGADS_SUCCESS) return stat;
        d = sqrt((xyz[0]-coord[0])*(xyz[0]-coord[0]) +
                 (xyz[1]-coord[1])*(xyz[1]-coord[1]) +
                 (xyz[2]-coord[2])*(xyz[2]-coord[2]));
        if (d < dist) {
          if (d <= pface->tol) return EGADS_SUCCESS;
          dist  = d;
          face  = pshell->faces[i];
          uv[0] = param[0];
          uv[1] = param[1];
        }
      }
    }
    if (face == NULL) return EGADS_OUTSIDE;
    stat = EG_evaluate(face, uv, data);
    if (stat != EGADS_SUCCESS) return stat;
    coord[0] = data[3];
    coord[1] = data[4];
    coord[2] = data[5];
    dir[0]   = data[6];
    dir[1]   = data[7];
    dir[2]   = data[8];
    CROSS(norm, coord, dir);
    if (face->mtype == SREVERSE) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
      norm[2] = -norm[2];
    }
    d        = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    if (d == 0.0) return EGADS_DEGEN;
    norm[0] /= d;
    norm[1] /= d;
    norm[2] /= d;
    
    dir[0]   = xyz[0] - data[0];
    dir[1]   = xyz[1] - data[1];
    dir[2]   = xyz[2] - data[2];
    d        = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    dir[0]  /= d;
    dir[1]  /= d;
    dir[2]  /= d;
    if (dir[0]*norm[0]+dir[1]*norm[1]+dir[2]*norm[2] > 0.0) return EGADS_OUTSIDE;
    return EGADS_SUCCESS;
  }
  
  return EGADS_NOTTOPO;
}
