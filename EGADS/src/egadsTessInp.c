/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Tessellation Input Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __CUDACC__
#include "liteString.h"
#else
#include <string.h>
#endif

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsTris.h"
#include "emp.h"

#include "regQuads.h"


  typedef struct {
    void     *mutex;              /* the mutex or NULL for single thread */
    long     master;              /* master thread ID */
    int      end;                 /* end of loop */
    int      index;               /* current loop index */
    egTessel *ntess;              /* tessellation structure */
  /*@dependent@*/
    bodyQuad *bodydata;           /* the quad storage */
  } EMPquad;


  typedef struct {
    int vert2;                    /* the second vert */
    int nvert;                    /* the new vert */
    int next;                     /* the next in the list for the first vert */
  } midside;


#define CRXSS(a,b,c)      c[0] = ((a)[1]*(b)[2]) - ((a)[2]*(b)[1]);\
                          c[1] = ((a)[2]*(b)[0]) - ((a)[0]*(b)[2]);\
                          c[2] = ((a)[0]*(b)[1]) - ((a)[1]*(b)[0])

#define FACE_NORMAL_FOLD_CHECK
#define REGULAR
#define NOTFILLED       -1
#define EPS10           1.e-10

#ifdef DEBUG
#define REPORT
#endif

#ifdef REPORT
#include <time.h>
#endif

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif
#ifdef __PROTO_H_AND_D__
#undef __PROTO_H_AND_D__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#define __PROTO_H_AND_D__   extern "C" __host__ __device__
#else
#define __HOST_AND_DEVICE__
#define __PROTO_H_AND_D__ extern
#endif

__PROTO_H_AND_D__ void EG_cleanupTess( egTessel *btess );
__PROTO_H_AND_D__ void EG_cleanupTessMaps( egTessel *btess );
__PROTO_H_AND_D__ void EG_makeConnect( int k1, int k2, int *tri, int *kedge,
                                       int *ntable, connect *etable, int face );
__PROTO_H_AND_D__ int  EG_fillArea( int ncontours, const int *cntr,
                                    const double *vertices, int *triangles,
                                    int *n_fig8, int pass, fillArea *fa );

__PROTO_H_AND_D__ int  EG_getTopology( const egObject *topo, egObject **geom,
                                       int *oclas, int *type,
                                       /*@null@*/ double *limits, int *nChild,
                                       egObject ***children, int **senses );
__PROTO_H_AND_D__ int  EG_getBodyTopos( const egObject *body,
                                        /*@null@*/ egObject *src, int oclass,
                                        int *nto, /*@null@*/ egObject ***topo );
__PROTO_H_AND_D__ int  EG_indexBodyTopo( const egObject *body,
                                         const egObject *src );
__PROTO_H_AND_D__ int  EG_evaluate( const egObject *geom,
                                    /*@null@*/ const double *prm, double *dat );
__PROTO_H_AND_D__ int  EG_getRange( const egObject *geom, double *range,
                                    int *pflag );
__PROTO_H_AND_D__ int  EG_invEvaluate( const egObject *geom, double *xyz,
                                       double *param, double *result );
__PROTO_H_AND_D__ int  EG_invEvaluateGuess( const egObject *geom, double *xyz,
                                            double *param, double *result );
__PROTO_H_AND_D__ int  EG_getEdgeUV( const egObject *face, const egObject *edge,
                                     int sense, double t, double *result );
__PROTO_H_AND_D__ int  EG_getEdgeUVeval( const egObject *face,
                                         const egObject *topo, int sense,
                                         double t, double *result);
__PROTO_H_AND_D__ int  EG_getTessEdge( const egObject *tess, int indx, int *len,
                                       const double **xyz, const double **t );
__PROTO_H_AND_D__ int  EG_getTessFace( const egObject *tess, int indx, int *len,
                                       const double **xyz, const double **uv,
                                       const int **ptype, const int **pindex,
                                       int *ntri, const int **tris,
                                                  const int **tric );
__PROTO_H_AND_D__ int  EG_attributeRet( const egObject *obj, const char *name,
                                        int *type, int *len,
                                        /*@null@*/ const int **ints,
                                        /*@null@*/ const double **reals,
                                        /*@null@*/ const char   **str );
#ifndef LITE
           extern int  EG_attributeDel( egObject *obj,
                                        /*@null@*/ const char *name );
           extern int  EG_attributeAdd( egObject *obj, const char *name,
                                        int type, int len,
                                        /*@null@*/ const int    *ints,
                                        /*@null@*/ const double *reals,
                                        /*@null@*/ const char   *str );
#endif


__HOST_AND_DEVICE__ int
EG_openTessBody(egObject *tess)
{
  egTessel *btess;
  egObject *obj;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (EG_sameThread(tess))          return EGADS_CNTXTHRD;
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_openTessBody)!\n");
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    printf(" EGADS Error: NULL Source Object (EG_openTessBody)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    printf(" EGADS Error: Source Not an Object (EG_openTessBody)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    printf(" EGADS Error: Source Not Body (EG_openTessBody)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->done != 1) return EGADS_TESSTATE;

  /* set open state and clean up any local/global mappings */
  btess->done = 0;
  EG_cleanupTessMaps(btess);

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_initTessBody(egObject *object, egObject **tess)
{
  int      i, j, k, n, stat, outLevel, nedge, nloop, nface, oclass, mtype;
  int      nnode, lor, ndum, *senses, *lsense, *finds;
  double   limits[4];
  egTessel *btess;
  egObject *ttess, *context, *geom, **faces, **loops, **edges, **nodes, **dum;

  *tess = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass != BODY)       return EGADS_NOTBODY;
  if (EG_sameThread(object))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);
  if (context == NULL)              return EGADS_NULLOBJ;

  stat = EG_getBodyTopos(object, NULL, EDGE, &nedge, &edges);
  if (stat != EGADS_SUCCESS) return stat;
  stat = EG_getBodyTopos(object, NULL, FACE, &nface, &faces);
  if (stat  != EGADS_SUCCESS) return stat;

  btess = (egTessel *) EG_alloc(sizeof(egTessel));
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Blind Malloc (EG_initTessBody)!\n");
    EG_free(faces);
    EG_free(edges);
    return EGADS_MALLOC;
  }
  btess->src       = object;
  btess->xyzs      = NULL;
  btess->tess1d    = NULL;
  btess->tess2d    = NULL;
  btess->globals   = NULL;
  btess->nGlobal   = 0;
  btess->nEdge     = nedge;
  btess->nFace     = nface;
  btess->nu        = 0;
  btess->nv        = 0;
  btess->done      = 0;
  btess->params[0] = 0.0;
  btess->params[1] = 0.0;
  btess->params[2] = 0.0;
  for (i = 0; i < MTESSPARAM; i++) btess->tparam[i] = 0.0;

  btess->tess1d = (egTess1D *) EG_alloc(nedge*sizeof(egTess1D));
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Alloc %d Edges (EG_initTessBody)!\n", nedge);
    EG_free(faces);
    EG_free(edges);
    EG_free(btess);
    return EGADS_MALLOC;
  }
  for (j = 0; j < nedge; j++) {
    btess->tess1d[j].obj            = edges[j];
    btess->tess1d[j].faces[0].index = 0;
    btess->tess1d[j].faces[0].nface = 0;
    btess->tess1d[j].faces[0].faces = NULL;
    btess->tess1d[j].faces[0].tric  = NULL;
    btess->tess1d[j].faces[1].index = 0;
    btess->tess1d[j].faces[1].nface = 0;
    btess->tess1d[j].faces[1].faces = NULL;
    btess->tess1d[j].faces[1].tric  = NULL;
    btess->tess1d[j].nodes[0]       = 0;
    btess->tess1d[j].nodes[1]       = 0;
    btess->tess1d[j].xyz            = NULL;
    btess->tess1d[j].t              = NULL;
    btess->tess1d[j].global         = NULL;
    btess->tess1d[j].npts           = 0;
  }
  EG_free(edges);
  for (j = 0; j < nedge; j++) {
    stat = EG_getTopology(btess->tess1d[j].obj, &geom, &oclass, &mtype, limits,
                          &nnode, &nodes, &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTopology = %d (EG_initTessBody)!\n", stat);
      EG_free(faces);
      EG_cleanupTess(btess);
      EG_free(btess);
      return EGADS_MALLOC;
    }
    stat = EG_indexBodyTopo(object, nodes[0]);
    if (stat <= EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_indexBodyTopo0 = %d (EG_initTessBody)!\n",
               stat);
      EG_free(faces);
      EG_cleanupTess(btess);
      EG_free(btess);
      return EGADS_MALLOC;
    }
    btess->tess1d[j].nodes[0] = btess->tess1d[j].nodes[1] =  stat;
    if (mtype == DEGENERATE)    btess->tess1d[j].nodes[1] = -stat;
    if (nnode == 1) continue;
    stat = EG_indexBodyTopo(object, nodes[1]);
    if (stat < EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_indexBodyTopo1 = %d (EG_initTessBody)!\n",
               stat);
      EG_free(faces);
      EG_cleanupTess(btess);
      EG_free(btess);
      return EGADS_MALLOC;
    }
    btess->tess1d[j].nodes[1] = stat;
  }

  if (nface != 0) {
    /* get the Edge Face indices */
    for (i = 0; i < nface; i++) {
      stat = EG_getTopology(faces[i], &geom, &oclass, &mtype, limits,
                            &nloop, &loops, &lsense);
      if (stat != EGADS_SUCCESS) continue;
      for (j = 0; j < nloop; j++) {
        lor = 1;
        if ((lsense[j] == 2) || (lsense[j] == -2)) lor = -1;
        stat = EG_getTopology(loops[j], &geom, &oclass, &mtype, limits,
                              &ndum, &dum, &senses);
        if (stat != EGADS_SUCCESS) continue;
        for (k = 0; k < ndum; k++) {
          n = EG_indexBodyTopo(object, dum[k]);
          if (n <= EGADS_SUCCESS) continue;
          if (senses[k]*lor < 0) {
            if (btess->tess1d[n-1].faces[0].nface != 0) {
              if (btess->tess1d[n-1].faces[0].nface == 1) {
                btess->tess1d[n-1].faces[0].faces = (int *) EG_alloc(2*sizeof(int));
                if (btess->tess1d[n-1].faces[0].faces == NULL) {
                  if (outLevel > 0)
                    printf(" EGADS Error: Alloc (-) Edge %d (EG_initTessBody)!\n",
                           n);
                  EG_free(faces);
                  EG_cleanupTess(btess);
                  EG_free(btess);
                  return EGADS_MALLOC;
                }
                btess->tess1d[n-1].faces[0].faces[0] = btess->tess1d[n-1].faces[0].index;
                btess->tess1d[n-1].faces[0].faces[1] = i+1;
              } else {
                finds = (int *) EG_reall( btess->tess1d[n-1].faces[0].faces,
                                         (btess->tess1d[n-1].faces[0].nface+1)*sizeof(int));
                if (finds == NULL) {
                  if (outLevel > 0)
                    printf(" EGADS Error: ReAlloc (-) Edge %d (EG_initTessBody)!\n",
                           n);
                  EG_free(faces);
                  EG_cleanupTess(btess);
                  EG_free(btess);
                  return EGADS_MALLOC;
                }
                finds[btess->tess1d[n-1].faces[0].nface] = i+1;
                btess->tess1d[n-1].faces[0].faces = finds;
              }
            }
            btess->tess1d[n-1].faces[0].index = i+1;
            btess->tess1d[n-1].faces[0].nface++;
          } else {
            if (btess->tess1d[n-1].faces[1].nface != 0) {
              if (btess->tess1d[n-1].faces[1].nface == 1) {
                btess->tess1d[n-1].faces[1].faces = (int *) EG_alloc(2*sizeof(int));
                if (btess->tess1d[n-1].faces[1].faces == NULL) {
                  if (outLevel > 0)
                    printf(" EGADS Error: Alloc (+) Edge %d (EG_initTessBody)!\n",
                           n);
                  EG_free(faces);
                  EG_cleanupTess(btess);
                  EG_free(btess);
                  return EGADS_MALLOC;
                }
                btess->tess1d[n-1].faces[1].faces[0] = btess->tess1d[n-1].faces[1].index;
                btess->tess1d[n-1].faces[1].faces[1] = i+1;
              } else {
                finds = (int *) EG_reall( btess->tess1d[n-1].faces[1].faces,
                                         (btess->tess1d[n-1].faces[1].nface+1)*sizeof(int));
                if (finds == NULL) {
                  if (outLevel > 0)
                    printf(" EGADS Error: ReAlloc (+) Edge %d (EG_initTessBody)!\n",
                           n);
                  EG_free(faces);
                  EG_cleanupTess(btess);
                  EG_free(btess);
                  return EGADS_MALLOC;
                }
                finds[btess->tess1d[n-1].faces[1].nface] = i+1;
                btess->tess1d[n-1].faces[1].faces = finds;
              }
            }
            btess->tess1d[n-1].faces[1].index = i+1;
            btess->tess1d[n-1].faces[1].nface++;
          }
        }
      }
    }
    EG_free(faces);
    
    btess->tess2d = (egTess2D *) EG_alloc(2*nface*sizeof(egTess2D));
    if (btess->tess2d == NULL) {
      printf(" EGADS Error: Alloc %d Faces (EG_initTessBody)!\n", nface);
      EG_cleanupTess(btess);
      EG_free(btess);
      return EGADS_MALLOC;
    }
    for (j = 0; j < 2*nface; j++) {
      btess->tess2d[j].mKnots = NULL;
      btess->tess2d[j].xyz    = NULL;
      btess->tess2d[j].uv     = NULL;
      btess->tess2d[j].global = NULL;
      btess->tess2d[j].ptype  = NULL;
      btess->tess2d[j].pindex = NULL;
      btess->tess2d[j].bary   = NULL;
      btess->tess2d[j].frame  = NULL;
      btess->tess2d[j].frlps  = NULL;
      btess->tess2d[j].tris   = NULL;
      btess->tess2d[j].tric   = NULL;
      btess->tess2d[j].patch  = NULL;
      btess->tess2d[j].npts   = 0;
      btess->tess2d[j].nframe = 0;
      btess->tess2d[j].nfrlps = 0;
      btess->tess2d[j].ntris  = 0;
      btess->tess2d[j].npatch = 0;
      btess->tess2d[j].tfi    = 0;
    }
  }

  stat = EG_makeObject(context, &ttess);
  if (stat != EGADS_SUCCESS) {
    EG_cleanupTess(btess);
    EG_free(btess);
    return stat;
  }
  ttess->oclass = TESSELLATION;
  ttess->blind  = btess;
  EG_referenceObject(ttess,  context);
  EG_referenceTopObj(object, ttess);
  *tess = ttess;

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_computeTessMap(egTessel *btess, int outLevel)
{
  int i, j, k, n, npts, pt, pi, nNode, *inode;

  if (btess->nGlobal !=    0) return EGADS_EXISTS;
  if (btess->globals != NULL) return EGADS_EXISTS;

  /* special case of a degenerate WIREBODY -- NodeBody */
  if ((btess->nEdge == 1) && (btess->tess1d[0].obj->mtype == DEGENERATE) &&
      (btess->tess1d[0].nodes[0] == 1)) {
    btess->xyzs = (double *) EG_alloc(3*sizeof(double));
    if (btess->xyzs == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation of 1 Nodes (EG_computeTessMap)!\n");
      return EGADS_MALLOC;
    }
    btess->xyzs[0] = btess->tess1d[0].xyz[0];
    btess->xyzs[1] = btess->tess1d[0].xyz[1];
    btess->xyzs[2] = btess->tess1d[0].xyz[2];

    btess->tess1d[0].global = (int *) EG_alloc(2*sizeof(int));
    if (btess->tess1d[0].global == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation of 1 Global (EG_computeTessMap)!\n");
      return EGADS_MALLOC;
    }
    btess->tess1d[0].global[0] = 1;
    btess->tess1d[0].global[1] = 1;

    btess->globals = (int *) EG_alloc(2*sizeof(int));
    if (btess->globals == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation of 1 Nodes (EG_computeTessMap)!\n");
      return EGADS_MALLOC;
    }
    btess->globals[0] = 0;
    btess->globals[1] = 1;
    btess->nGlobal    = 1;
    return EGADS_SUCCESS;
  }

  /* get Node and Edge sizes */
  for (nNode = npts = i = 0; i < btess->nEdge; i++) {
    if (btess->tess1d[i].obj        ==       NULL) continue;
    if (btess->tess1d[i].obj->mtype == DEGENERATE) continue;
    if (btess->tess1d[i].nodes[1]   <           0) continue;
    npts += btess->tess1d[i].npts;
    if (nNode < btess->tess1d[i].nodes[0]) nNode = btess->tess1d[i].nodes[0];
    if (nNode < btess->tess1d[i].nodes[1]) nNode = btess->tess1d[i].nodes[1];
  }
  inode = (int *) EG_alloc(nNode*sizeof(int));
  if (inode == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation of %d Nodes (EG_computeTessMap)!\n",
             nNode);
    return EGADS_MALLOC;
  }
  btess->xyzs = (double *) EG_alloc(3*nNode*sizeof(double));
  if (btess->xyzs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation Of %d Nodes (EG_computeTessMap)!\n",
             nNode);
    EG_free(inode);
    return EGADS_MALLOC;
  }

  for (i = 0; i < nNode; i++) {
    inode[i]           = 0;
    btess->xyzs[3*i  ] = 0.0;
    btess->xyzs[3*i+1] = 0.0;
    btess->xyzs[3*i+2] = 0.0;
  }
  for (i = 0; i < btess->nEdge; i++) {
    if (btess->tess1d[i].obj        ==       NULL) continue;
    if (btess->tess1d[i].obj->mtype == DEGENERATE) continue;
    if (btess->tess1d[i].nodes[1]   <           0) continue;
    j = btess->tess1d[i].nodes[0]-1;
    btess->xyzs[3*j  ] = btess->tess1d[i].xyz[0];
    btess->xyzs[3*j+1] = btess->tess1d[i].xyz[1];
    btess->xyzs[3*j+2] = btess->tess1d[i].xyz[2];
    inode[j]++;
    j = btess->tess1d[i].nodes[1]-1;
    k = btess->tess1d[i].npts;
    btess->xyzs[3*j  ] = btess->tess1d[i].xyz[3*k-3];
    btess->xyzs[3*j+1] = btess->tess1d[i].xyz[3*k-2];
    btess->xyzs[3*j+2] = btess->tess1d[i].xyz[3*k-1];
    inode[j]++;
    btess->tess1d[i].global = (int *) EG_alloc(k*sizeof(int));
    if (btess->tess1d[i].global == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: %d Allocation of %d Global (EG_computeTessMap)!\n",
               i+1, k);
      EG_cleanupTessMaps(btess);
      EG_free(inode);
      return EGADS_MALLOC;
    }
    for (j = 0; j < k; j++) btess->tess1d[i].global[j] = 0;
  }

  if (btess->nFace == 0) {

    /* deal with wirebodies */

    for (i = 0; i < nNode; i++)
      if (inode[i] != 0) npts -= inode[i]-1;

    btess->globals = (int *) EG_alloc(2*npts*sizeof(int));
    if (btess->globals == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation of %d Globals (EG_computeTessMap)!\n",
               npts);
      EG_cleanupTessMaps(btess);
      EG_free(inode);
      return EGADS_MALLOC;
    }

    for (k = i = 0; i < btess->nEdge; i++) {
      if (btess->tess1d[i].obj        ==       NULL) continue;
      if (btess->tess1d[i].obj->mtype == DEGENERATE) continue;
      if (btess->tess1d[i].nodes[1]   <           0) continue;
      n = btess->tess1d[i].nodes[0] - 1;
      if (inode[n] > 0) {
        btess->globals[2*k  ]      = 0;
        btess->globals[2*k+1]      = n+1;
        k++;
        btess->tess1d[i].global[0] = k;
        inode[n] = -k;
      } else {
        btess->tess1d[i].global[0] = -inode[n];
      }

      for (j = 1; j < btess->tess1d[i].npts-1; j++) {
        btess->globals[2*k  ]      = j+1;
        btess->globals[2*k+1]      = i+1;
        k++;
        btess->tess1d[i].global[j] = k;
      }

      n = btess->tess1d[i].nodes[1] - 1;
      j = btess->tess1d[i].npts     - 1;
      if (inode[n] > 0) {
        btess->globals[2*k  ]      = 0;
        btess->globals[2*k+1]      = n+1;
        k++;
        btess->tess1d[i].global[j] = k;
        inode[n] = -k;
      } else {
        btess->tess1d[i].global[j] = -inode[n];
      }

    }

    EG_free(inode);
    btess->nGlobal = npts;
    return EGADS_SUCCESS;
  }

  for (k = i = 0; i < btess->nFace; i++) {
    n = btess->tess2d[i].npts;
    if (n == 0) continue;
    btess->tess2d[i].global = (int *) EG_alloc(n*sizeof(int));
    if (btess->tess2d[i].global == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: %d Allocation of %d Globals (EG_computeTessMap)!\n",
               i+1, n);
      EG_cleanupTessMaps(btess);
      EG_free(inode);
      return EGADS_MALLOC;
    }
    for (j = 0; j < n; j++) btess->tess2d[i].global[j] = 0;
    k += n;
  }
  btess->globals = (int *) EG_alloc(2*k*sizeof(int));
  if (btess->globals == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation of %d Globals (EG_computeTessMap)!\n",
             k);
    EG_cleanupTessMaps(btess);
    EG_free(inode);
    return EGADS_MALLOC;
  }

  for (i = 0; i < nNode; i++) inode[i] = 0;
  for (k = i = 0; i < btess->nFace; i++) {
    n = btess->tess2d[i].npts;
    for (j = 0; j < n; j++) {
      pt = btess->tess2d[i].ptype[j];
      pi = btess->tess2d[i].pindex[j];
      if (pt == 0) {
        if (inode[pi-1] == 0) {
          btess->globals[2*k  ]            = pt;
          btess->globals[2*k+1]            = pi;
          k++;
          btess->tess2d[i].global[j]       = k;
          inode[pi-1]                      = k;
        } else {
          btess->tess2d[i].global[j]       = inode[pi-1];
        }
      } else if (pt > 0) {
        if (btess->tess1d[pi-1].global[pt-1] == 0) {
          btess->globals[2*k  ]            = pt;
          btess->globals[2*k+1]            = pi;
          k++;
          btess->tess2d[i].global[j]       = k;
          btess->tess1d[pi-1].global[pt-1] = k;
        } else {
          btess->tess2d[i].global[j]       = btess->tess1d[pi-1].global[pt-1];
        }
      } else {
        btess->globals[2*k  ]              = -j-1;
        btess->globals[2*k+1]              =  i+1;
        k++;
        btess->tess2d[i].global[j]         =  k;
      }
    }
  }

  /* patch up beginning and end of Edges */
  for (i = 0; i < btess->nEdge; i++) {
    if (btess->tess1d[i].obj        ==       NULL) continue;
    if (btess->tess1d[i].obj->mtype == DEGENERATE) continue;
    if (btess->tess1d[i].nodes[1]   <           0) continue;
    n = btess->tess1d[i].nodes[0] - 1;
    if (inode[n] == 0) {
      btess->globals[2*k  ]      = 0;
      btess->globals[2*k+1]      = n+1;
      k++;
      inode[n]                   = k;
    } else {
      btess->tess1d[i].global[0] = inode[n];
    }
    n = btess->tess1d[i].nodes[1] - 1;
    j = btess->tess1d[i].npts     - 1;
    if (inode[n] == 0) {
      btess->globals[2*k  ]      = 0;
      btess->globals[2*k+1]      = n+1;
      k++;
      inode[n]                   = k;
    } else {
      btess->tess1d[i].global[j] = inode[n];
    }
  }

  EG_free(inode);
  btess->nGlobal = k;

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_statusTessBody(egObject *tess, egObject **body, int *state, int *npts)
{
  int          i, j, k, stat, outLevel, atype, alen;
  egTessel     *btess;
  egObject     *obj;
  const int    *ints;
  const double *reals;
  const char   *str;

  *body  = NULL;
  *state = *npts = 0;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_statusTessBody)!\n");
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_statusTessBody)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_statusTessBody)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_statusTessBody)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_statusTessBody)!\n");
    return EGADS_NODATA;
  }
  if ((btess->tess2d == NULL) && (btess->nFace != 0)) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_statusTessBody)!\n");
    return EGADS_NODATA;
  }
  *body  = obj;
  *state = btess->done;
  *npts  = btess->nGlobal;

  if (btess->done == 0) {

    /* are we done? */
    for (i = 0; i < btess->nEdge; i++) {
      if (btess->tess1d[i].nodes[0] == -btess->tess1d[i].nodes[1]) continue;
      if (btess->tess1d[i].npts == 0) return EGADS_OUTSIDE;
    }
    for (i = 0; i < btess->nFace; i++)
      if (btess->tess2d[i].npts == 0) return EGADS_OUTSIDE;

    stat = EG_attributeRet(tess, ".mixed", &atype, &alen, &ints, &reals, &str);
    if (stat == EGADS_SUCCESS) {
      if ((alen != btess->nFace) || (atype != ATTRINT)) {
#ifndef LITE
        stat = EG_attributeDel(tess, ".mixed");
        if ((stat != EGADS_SUCCESS) && (outLevel > 0))
          printf(" EGADS Error: Deleting Attribute %d (EG_statusTessBody)!\n",
                 stat);
#endif
      } else {
        for (j = i = 0; i < btess->nFace; i++) {
          if (2*ints[i] == btess->tess2d[i].ntris) j++;
          if (2*ints[i] >  btess->tess2d[i].ntris) {
            printf(" EGADS Error: %d 2*nQuads (%d) > nTris (%d) (EG_statusTessBody)!\n",
                   i+1, ints[i], btess->tess2d[i].ntris);
            j = -1;
            break;
          }
          /* check pairs */
          for (k = btess->tess2d[i].ntris-2*ints[i]; k < btess->tess2d[i].ntris;
               k += 2)
            if ((btess->tess2d[i].tris[3*k  ] != btess->tess2d[i].tris[3*k+3]) ||
                (btess->tess2d[i].tris[3*k+2] != btess->tess2d[i].tris[3*k+4])) {
              printf(" EGADS Error: %d %d Bad Pair %d %d %d - %d %d %d!\n", i+1,
                     k+1, btess->tess2d[i].tris[3*k], btess->tess2d[i].tris[3*k+1],
                     btess->tess2d[i].tris[3*k+2], btess->tess2d[i].tris[3*k+3],
                     btess->tess2d[i].tris[3*k+4], btess->tess2d[i].tris[3*k+5]);
              j = -1;
              break;
            }
          if (j == -1) break;
        }
        if (j == -1) {
#ifndef LITE
          stat = EG_attributeDel(tess, ".mixed");
          if ((stat != EGADS_SUCCESS) && (outLevel > 0))
            printf(" EGADS Error: Deleting Attribute %d (EG_statusTessBody)!\n",
                   stat);
#endif
        } else if (j == btess->nFace) {
#ifndef LITE
          stat = EG_attributeAdd(tess, ".tessType", ATTRSTRING, 4, NULL, NULL,
                                 "Quad");
          if (stat != EGADS_SUCCESS)
            if (outLevel > 0)
              printf(" EGADS Warning: EG_attributeAdd = %d (EG_statusTessBody)!\n",
                     stat);
#endif
        } else {
#ifndef LITE
          stat = EG_attributeAdd(tess, ".tessType", ATTRSTRING, 5, NULL, NULL,
                                 "Mixed");
          if (stat != EGADS_SUCCESS)
            if (outLevel > 0)
              printf(" EGADS Warning: EG_attributeAdd = %d (EG_statusTessBody)!\n",
                     stat);
#endif
        }
      }
    }

    *state = btess->done = 1;
  }

  if (btess->globals != NULL) return EGADS_SUCCESS;

  /* compute the mappings and return the number of global vertices */
  stat = EG_computeTessMap(btess, outLevel);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_computeTessMap = %d (EG_statusTessBody)!\n",
              stat);
    return stat;
  }
  *npts = btess->nGlobal;

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_setTessEdge(const egObject *tess, int index, int len, const double *xyz,
               const double *t)
{
  int      i, j, k, n, stat, outLevel, oclass, mtype, nnode, *senses;
  double   xyz0[3], xyz1[3], trange[2], *xyzs, *ts;
  egTessel *btess;
  egObject *obj, **nodes, *geom, **objs;

  if  (tess == NULL)                 return EGADS_NULLOBJ;
  if  (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if  (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if  (EG_sameThread(tess))          return EGADS_CNTXTHRD;
  if  (len <= 1)                     return EGADS_NODATA;
  if ((xyz == NULL) || (t == NULL))  return EGADS_NODATA;
  outLevel = EG_outLevel(tess);

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_setTessEdge)!\n");
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_setTessEdge)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_setTessEdge)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_setTessEdge)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_setTessEdge)!\n");
    return EGADS_NODATA;
  }
  if (btess->done == 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Complete Tessellation (EG_setTessEdge)!\n");
    return EGADS_EXISTS;
  }
  if ((index < 1) || (index > btess->nEdge)) {
    if (outLevel > 0)
      printf(" EGADS Error: Index = %d [1-%d] (EG_setTessEdge)!\n",
             index, btess->nEdge);
    return EGADS_INDEXERR;
  }

  for (i = 1; i < len; i++) {
    if (t[i] > t[i-1]) continue;
    printf(" EGADS Error: ts not in order %d %lf  %d %lf (EG_setTessEdge)!\n",
           i-1, t[i-1], i, t[i]);
    return EGADS_RANGERR;
  }

  /* are any of our Faces already set? */
  stat = EG_getBodyTopos(obj, btess->tess1d[index-1].obj, FACE, &n, &objs);
  if (stat  != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge %d - EG_getBodyTopos = %d (EG_setTessEdge)!\n",
             index, stat);
    return stat;
  }
  if ((n != 0) && (objs != NULL))
    for (i = 0; i < n; i++) {
      j = EG_indexBodyTopo(obj, objs[i]);
      if (j <= EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d - indexTopoBody %d = %d (EG_setTessEdge)!\n",
                 index, i+1, j);
        EG_free(objs);
        return j;
      }
      if (btess->tess2d[j-1].npts != 0) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d - Face %d set (EG_setTessEdge)!\n",
                 index, j);
        EG_free(objs);
        return EGADS_EXISTS;
      }
    }
  if (objs != NULL) EG_free(objs);

  /* get the bounding information */
  stat = EG_getTopology(btess->tess1d[index-1].obj, &geom, &oclass, &mtype,
                        trange, &nnode, &nodes, &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getTopology = %d for Edge %d (EG_setTessEdge)!\n",
             stat, index);
    return stat;
  }

  stat = EG_getTopology(nodes[0], &geom, &oclass, &mtype, xyz0, &i, &objs,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getTopology = %d - Edge %d/0 (EG_setTessEdge)!\n",
             stat, index);
    return stat;
  }
  j = nnode - 1;
  stat = EG_getTopology(nodes[j], &geom, &oclass, &mtype, xyz1, &i, &objs,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getTopology = %d - Edge %d/1 (EG_setTessEdge)!\n",
             stat, index);
    return stat;
  }

  /* allocate the data */
  xyzs = (double *) EG_alloc(3*len*sizeof(double));
  if (xyzs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocating %d Coordinates (EG_setTessEdge)!\n",
             len);
    return EGADS_MALLOC;
  }
  ts = (double *) EG_alloc(len*sizeof(double));
  if (ts == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocating %d Parameters (EG_setTessEdge)!\n",
             len);
    EG_free(xyzs);
    return EGADS_MALLOC;
  }
#ifdef SETNODE
  if ((xyz[0] != xyz0[0]) || (xyz[1] != xyz0[1]) || (xyz[2] != xyz0[2])) {
    printf(" EGADS Warning: Edge %d- Node %d  XYZ misMatch (EG_setTessEdge)!\n",
           index, EG_indexBodyTopo(obj, nodes[0]));
  }
  ts[0]   = trange[0];
  xyzs[0] = xyz0[0];
  xyzs[1] = xyz0[1];
  xyzs[2] = xyz0[2];
  for (i = 1; i < len-1; i++) {
    ts[i]       = t[i];
    xyzs[3*i  ] = xyz[3*i  ];
    xyzs[3*i+1] = xyz[3*i+1];
    xyzs[3*i+2] = xyz[3*i+2];
  }
  if ((xyz[3*len-3] != xyz1[0]) || (xyz[3*len-2] != xyz1[1]) ||
      (xyz[3*len-1] != xyz1[2])) {
    printf(" EGADS Warning: Edge %d+ Node %d  XYZ misMatch (EG_setTessEdge)!\n",
           index, EG_indexBodyTopo(obj, nodes[j]));
  }
  ts[len-1]     = trange[1];
  xyzs[3*len-3] = xyz1[0];
  xyzs[3*len-2] = xyz1[1];
  xyzs[3*len-1] = xyz1[2];
#else
  for (i = 0; i < len; i++) {
    ts[i]       = t[i];
    xyzs[3*i  ] = xyz[3*i  ];
    xyzs[3*i+1] = xyz[3*i+1];
    xyzs[3*i+2] = xyz[3*i+2];
  }
#endif

  /* set the data */
  if (btess->tess1d[index-1].xyz != NULL) EG_free(btess->tess1d[index-1].xyz);
  if (btess->tess1d[index-1].t   != NULL) EG_free(btess->tess1d[index-1].t);

  btess->tess1d[index-1].npts = len;
  btess->tess1d[index-1].xyz  = xyzs;
  btess->tess1d[index-1].t    = ts;

  if (n > 0) {
    if (btess->tess1d[index-1].faces[0].tric == NULL)
      btess->tess1d[index-1].faces[0].tric = (int *)
                                             EG_alloc((n*(len-1))*sizeof(int));
    if (btess->tess1d[index-1].faces[0].tric == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: Alloc %d Tric- Edge %d (EG_setTessEdge)!\n",
               len, index);
    } else {
      for (i = 0; i < len-1; i++)
        for (k = 0; k < n; k++)
          btess->tess1d[index-1].faces[0].tric[i*n+k] = 0;
    }

    if (btess->tess1d[index-1].faces[1].tric == NULL)
      btess->tess1d[index-1].faces[1].tric = (int *)
                                             EG_alloc((n*(len-1))*sizeof(int));
    if (btess->tess1d[index-1].faces[1].tric == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: Alloc %d Tric+ Edge %d (EG_setTessEdge)!\n",
               len, index);
    } else {
      for (i = 0; i < len-1; i++)
        for (k = 0; k < n; k++)
          btess->tess1d[index-1].faces[1].tric[i*n+k] = 0;
    }
  }

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
findPoint(double *range, double *uvp, int ptype, int pindex,
          int len, int *table, const double *uv, double *uvx)
{
  int    i, j;
  double du, dv;

  for (i = 0; i < len; i++)
    if ((ptype == table[2*i  ]) && (pindex == table[2*i+1])) {
      du = fabs(uv[2*i  ] - uvp[0])/(range[1]-range[0]);
      dv = fabs(uv[2*i+1] - uvp[1])/(range[3]-range[2]);
      if ((du < 0.25) && (dv < 0.25)) {
        uvx[0] = uv[2*i  ];
        uvx[1] = uv[2*i+1];
        return i;
      }
    }

  /* Degenerate Node? -- update and let pass */
  if (ptype == 0) {
    for (j = i = 0; i < len; i++)
      if ((ptype == table[2*i  ]) && (pindex == table[2*i+1])) j++;
    if (j == 1)
      for (i = 0; i < len; i++)
        if ((ptype == table[2*i  ]) && (pindex == table[2*i+1])) {
/*        printf(" EGADS Info: Differing UV @ Node %d  %lf %lf -- %lf %lf!\n",
                 table[2*i+1], uvp[0], uvp[1], uv[2*i  ], uv[2*i+1]);  */
          du     = fabs(uv[2*i  ] - uvp[0])/(range[1]-range[0]);
          dv     = fabs(uv[2*i+1] - uvp[1])/(range[3]-range[2]);
          uvx[0] = uv[2*i  ];
          uvx[1] = uv[2*i+1];
          if (du >= 0.25) uvx[0] = uvp[0];
          if (dv >= 0.25) uvx[1] = uvp[1];
          return i;
        }
  }

  return EGADS_NOTFOUND;
}


__HOST_AND_DEVICE__ static int
makeNeighbors(int f, int nverts, int ntri, int *tris, int *tric,
              int nseg, triSeg *segs)
{
  int     *ntab, nside, j;
  connect *etab;

  ntab = (int *) EG_alloc(nverts*sizeof(int));
  if (ntab == NULL) {
    printf(" EGADS Error: Vert Table Malloc (EG_setTessFace)!\n");
    return EGADS_MALLOC;
  }
  etab = (connect *) EG_alloc(ntri*3*sizeof(connect));
  if (etab == NULL) {
    printf(" EGADS Error: Edge Table Malloc (EG_setTessFace)!\n");
    EG_free(ntab);
    return EGADS_MALLOC;
  }

  nside = -1;
  for (j = 0; j < nverts; j++) ntab[j] = NOTFILLED;
  for (j = 0; j < ntri;  j++) {
    EG_makeConnect(tris[3*j+1], tris[3*j+2], &tric[3*j  ], &nside,ntab,etab, f);
    EG_makeConnect(tris[3*j  ], tris[3*j+2], &tric[3*j+1], &nside,ntab,etab, f);
    EG_makeConnect(tris[3*j  ], tris[3*j+1], &tric[3*j+2], &nside,ntab,etab, f);
  }

  for (j = 0; j < nseg; j++)
    EG_makeConnect(segs[j].indices[0], segs[j].indices[1], &segs[j].neighbor,
                   &nside, ntab, etab, f);

  /* report any unconnected triangle sides */
  for (j = 0; j <= nside; j++) {
    if (etab[j].tri == NULL) continue;
    printf(" EGADS Info: Face %d, Unconnected Side %d %d = %d\n",
           f, etab[j].node1+1, etab[j].node2+1, *etab[j].tri);
    *etab[j].tri = 0;
  }

  EG_free(etab);
  EG_free(ntab);
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_setTessFace(const egObject *tess, int index, int len, const double *xyz,
               const double *uv, int ntri, const int *tris)
{
  int      i, j, k, m, n, hit, iedge, outLevel, stat, nedge, *table, *map;
  int      oclass, mtype, nloop, np, sen, ori, lor, pt, pi, *senses, *lsenses;
  int      ntot, st, nseg, ntrix, nf8, nd, mm, mp;
  int      *frlps, *frame, *ptype, *pindex, *trix, *tric, *sns;
  double   smallu, smallv;
  double   range[4], trange[2], uvm[2], uvp[2], uvx[2], *uvs, *xyzs, *intEdg;
  triSeg   *segs;
  fillArea fast;
  egTessel *btess;
  egObject *obj, *geom, *face, **faces, **loops, **edges, **nds;
  static int    sides[3][2] = {{1,2}, {2,0}, {0,1}};
  static double scl[3][2]   = {{1.0, 1.0},  {10.0, 1.0},  {0.1, 10.0}};

  if  (tess == NULL)                 return EGADS_NULLOBJ;
  if  (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if  (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if  (EG_sameThread(tess))          return EGADS_CNTXTHRD;
  if ((len <= 1) || (ntri < 1))      return EGADS_NODATA;
  if ((xyz == NULL) || (uv == NULL)) return EGADS_NODATA;
  if  (tris == NULL)                 return EGADS_NODATA;
  outLevel = EG_outLevel(tess);

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_setTessFace)!\n");
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_setTessFace)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_setTessFace)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_setTessFace)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess2d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_setTessFace)!\n");
    return EGADS_NODATA;
  }
  if (btess->done == 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Complete Tessellation (EG_setTessFace)!\n");
    return EGADS_EXISTS;
  }
  if ((index < 1) || (index > btess->nFace)) {
    if (outLevel > 0)
      printf(" EGADS Error: Index = %d [1-%d] (EG_setTessFace)!\n",
             index, btess->nFace);
    return EGADS_INDEXERR;
  }

  /* check triangle indices */
  for (i = 0; i < ntri; i++) {
    for (j = 0; j < 3; j++) {
      if ((tris[3*i+j] < 1) || (tris[3*i+j] > len)) {
        printf(" EGADS Error: Face %d - tris %d/%d = %d [1-%d] (EG_setTessFace)!\n",
               index, i+1, j, tris[3*i+j], len);
        return EGADS_INDEXERR;
      }
    }
    if ((tris[3*i  ] == tris[3*i+1]) || (tris[3*i  ] == tris[3*i+2]) ||
        (tris[3*i+1] == tris[3*i+2])) {
      printf(" EGADS Error: Face %d - tris %d is degenerate = %d %d %d (EG_setTessFace)!\n",
             index, i+1, tris[3*i  ], tris[3*i+1], tris[3*i+2]);
      return EGADS_INDEXERR;
    }
  }

  /* get our Face object */
  stat = EG_getBodyTopos(obj, NULL, FACE, &i, &faces);
  if (stat  != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - EG_getBodyTopos Faces = %d (EG_setTessFace)!\n",
             index, stat);
    return stat;
  }
  face = faces[index-1];
  EG_free(faces);

  /* make sure we have all of the edge tessellations */
  stat = EG_getBodyTopos(obj, face, EDGE, &nedge, &edges);
  if (stat  != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - EG_getBodyTopos Edges = %d (EG_setTessFace)!\n",
             index, stat);
    return stat;
  }
  for (i = 0; i < nedge; i++) {
    iedge = EG_indexBodyTopo(obj, edges[i]);
    if (iedge <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d - EG_indexTopoBody Edge %d = %d (EG_setTessFace)!\n",
               index, i+1, iedge);
      return iedge;
    }
    if (edges[i]->mtype == DEGENERATE) continue;
    if (btess->tess1d[iedge-1].nodes[0] == -btess->tess1d[iedge-1].nodes[1])
      continue;
    if (btess->tess1d[iedge-1].npts > 0) continue;
    printf(" EGADS Error: Face %d - No Tessellation for Edge %d (EG_setTessFace)!\n",
           index, iedge);
    EG_free(edges);
    return EGADS_NOTFOUND;
  }

  /* setup ptype/pindex */
  table = (int *) EG_alloc(3*len*sizeof(int));
  if (table == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d Vert Table (EG_setTessFace)!\n",
             index, len);
    EG_free(edges);
    return EGADS_MALLOC;
  }
  map = &table[2*len];
  for (i = 0; i < len; i++)
    table[2*i  ] = table[2*i+1] = map[i] = -1;
  for (k = 0; k < len; k++) {
    for (hit = i = 0; i < nedge; i++) {
      iedge = EG_indexBodyTopo(obj, edges[i]);
      if (edges[i]->mtype == DEGENERATE) continue;
      if (btess->tess1d[iedge-1].nodes[0] == -btess->tess1d[iedge-1].nodes[1])
        continue;
      for (j = 0; j < btess->tess1d[iedge-1].npts; j++) {
        if (xyz[3*k  ] != btess->tess1d[iedge-1].xyz[3*j  ]) continue;
        if (xyz[3*k+1] != btess->tess1d[iedge-1].xyz[3*j+1]) continue;
        if (xyz[3*k+2] != btess->tess1d[iedge-1].xyz[3*j+2]) continue;
        table[2*k  ] = j+1;
        table[2*k+1] = iedge;
        if (j == 0) {
          table[2*k  ] = 0;
          table[2*k+1] = btess->tess1d[iedge-1].nodes[0];
        }
        if (j == btess->tess1d[iedge-1].npts-1) {
          table[2*k  ] = 0;
          table[2*k+1] = btess->tess1d[iedge-1].nodes[1];
        }
        hit++;
        break;
      }
      if (hit != 0) break;
    }
  }
  EG_free(edges);

  /* reorder based on loops */
  stat = EG_getTopology(face, &geom, &oclass, &ori, range, &nloop, &loops,
                        &lsenses);
  if (stat  != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - EG_getTopology = %d (EG_setTessFace)!\n",
             index, stat);
    EG_free(table);
    return stat;
  }
  smallu = 0.00005*(range[1] - range[0]);
  smallv = 0.00005*(range[3] - range[2]);

  /* get total number of points in all of the loops */
  for (ntot = i = 0; i < nloop; i++) {
    stat = EG_getTopology(loops[i], &geom, &oclass, &mtype, NULL, &nedge,
                          &edges, &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d - EG_getTopology Loop %d = %d (EG_setTessFace)!\n",
               index, i+1, stat);
      EG_free(table);
      return stat;
    }
    for (j = 0; j < nedge; j++) {
      iedge = EG_indexBodyTopo(obj, edges[j]);
      if (edges[j]->mtype == DEGENERATE) continue;
      if (btess->tess1d[iedge-1].nodes[0] == -btess->tess1d[iedge-1].nodes[1])
        continue;
      ntot += btess->tess1d[iedge-1].npts-1;
    }
  }
  ntrix = ntot-2 + 2*(nloop-1);
  segs  = (triSeg *) EG_alloc(ntot*sizeof(triSeg));
  if (segs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d Segs (EG_setTessFace)!\n",
             index, ntot);
    EG_free(table);
    return EGADS_MALLOC;
  }
  uvs = (double *) EG_alloc((2*ntot+2)*sizeof(double));
  if (uvs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d uvs (EG_setTessFace)!\n",
             index, ntot);
    EG_free(segs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  uvs[0] = uvs[1] = 0.0;

  /* find the Edge vertices */
  frlps = (int *) EG_alloc(nloop*sizeof(int));
  if (frlps == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d Loops (EG_setTessFace)!\n",
             index, nloop);
    EG_free(uvs);
    EG_free(segs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  for (np = i = 0; i < nloop; i++) {
    st   = np;
    stat = EG_getTopology(loops[i], &geom, &oclass, &mtype, NULL, &nedge,
                          &edges, &senses);
    if (stat != EGADS_SUCCESS) continue;
    lor = 1;
    if ((lsenses[i] == 2) || (lsenses[i] == -2)) lor = -1;
    n = 0;
    if (ori*lor == SREVERSE) n = nedge-1;
    for (j = 0; j < nedge; j++, n += ori*lor) {
      iedge = EG_indexBodyTopo(obj, edges[n]);
      if (edges[n]->mtype == DEGENERATE) continue;
      if (btess->tess1d[iedge-1].nodes[0] == -btess->tess1d[iedge-1].nodes[1])
        continue;
      stat = EG_getTopology(edges[n], &geom, &oclass, &mtype, trange, &nd,
                            &nds, &sns);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EG_getTopology = %d  for Face = %d, Edge = %d\n",
               stat, index, iedge);
        EG_free(frlps);
        EG_free(uvs);
        EG_free(segs);
        EG_free(table);
        return stat;
      }
      sen = senses[n]*ori*lor;

      /* internal Edge? */
      intEdg = NULL;
      for (m = 0; m < nedge; m++) {
        if (m == n) continue;
        if (iedge == EG_indexBodyTopo(obj, edges[m])) {
          uvm[0] = uvm[1] = -1.0;
          uvp[0] = uvp[1] =  1.0;
          EG_getEdgeUV(face, edges[n], -1, 0.5*(trange[0]+trange[1]), uvm);
          EG_getEdgeUV(face, edges[n],  1, 0.5*(trange[0]+trange[1]), uvp);
          if ((uvm[0] == uvp[0]) && (uvm[1] == uvp[1]) && (intEdg == NULL)) {
            printf(" EGADS Info: ");
            printf("Face #%d -> Edge #%d (%d) Internally in Loop %d %d, sen = %d!\n",
                   index, iedge, nedge, n+1, m+1, sen);
            intEdg = (double *)
                     EG_alloc(4*btess->tess1d[iedge-1].npts*sizeof(double));
            if (intEdg == NULL) {
              printf(" EGADS Internal: Cannot Allocate %d intEgdes!\n",
                     btess->tess1d[iedge-1].npts);
              continue;
            }
            for (m = 0; m < btess->tess1d[iedge-1].npts; m++) {
              stat = EG_getEdgeUV(face, edges[n], senses[n]*lor,
                                  btess->tess1d[iedge-1].t[m], &intEdg[4*m]);
              if (stat != EGADS_SUCCESS) {
                printf(" EGADS Error: getEdgeUV! = %d  for Face %d, Edge = %d\n",
                       stat, index, iedge);
                EG_free(intEdg);
                EG_free(frlps);
                EG_free(uvs);
                EG_free(segs);
                EG_free(table);
                return stat;
              }
            }
            for (m = 0; m < btess->tess1d[iedge-1].npts; m++) {
              mm = m - 1;
              mp = m + 1;
              if (mm <  0) mm = 0;
              if (mp >= btess->tess1d[iedge-1].npts)
                mp = btess->tess1d[iedge-1].npts - 1;
              uvm[0] = intEdg[4*mp  ] - intEdg[4*mm  ];
              uvm[1] = intEdg[4*mp+1] - intEdg[4*mm+1];
              uvp[0] = atan2(-uvm[1], uvm[0]);
              intEdg[4*m+2] = sen*smallu*sin(uvp[0]);
              intEdg[4*m+3] = sen*smallv*cos(uvp[0]);
            }
          }
        }
      }

      if (sen == 1) {
        for (m = 0; m < btess->tess1d[iedge-1].npts-1; m++, np++) {
          pt = m+1;
          pi = iedge;
          if (m == 0) {
            pt = 0;
            pi = btess->tess1d[iedge-1].nodes[0];
          }
          uvp[0] = uvp[1] = 0.0;
          stat   = EG_getEdgeUV(face, edges[n], senses[n],
                                btess->tess1d[iedge-1].t[m], uvp);
          if (stat != EGADS_SUCCESS)
            printf(" EGADS Internal: Face %d - EdgeUV+ = %d (EG_setTessFace)!\n",
                   index, stat);
          stat = findPoint(range, uvp, pt, pi, len, table, uv, uvx);
          if (stat < EGADS_SUCCESS) {
            printf(" EGADS Error: Face %d - FindPt+ %d/%d = %d (EG_setTessFace)!\n",
                   index, pt, pi, stat);
            if (intEdg != NULL) EG_free(intEdg);
            EG_free(frlps);
            EG_free(uvs);
            EG_free(segs);
            EG_free(table);
            return stat;
          }
          map[stat]           =  np;
          if (intEdg == NULL) {
            uvs[2*np+2]       =  uvx[0];
            uvs[2*np+3]       =  uvx[1];
          } else {
            uvs[2*np+2]       =  intEdg[4*m  ] + intEdg[4*m+2];
            uvs[2*np+3]       =  intEdg[4*m+1] + intEdg[4*m+3];
          }
          segs[np].indices[0] =  np+1;
          segs[np].indices[1] =  np+2;
          segs[np].neighbor   = -iedge;
          segs[np].edge       =  senses[n]*lor*iedge;
          segs[np].index      =  m+1;
        }
      } else {
        for (m = btess->tess1d[iedge-1].npts-1; m > 0; m--, np++) {
          pt = m+1;
          pi = iedge;
          if (m == btess->tess1d[iedge-1].npts-1) {
            pt = 0;
            pi = btess->tess1d[iedge-1].nodes[1];
          }
          uvm[0] = uvm[1] = 0.0;
          stat   = EG_getEdgeUV(face, edges[n], senses[n],
                                btess->tess1d[iedge-1].t[m], uvm);
          if (stat != EGADS_SUCCESS)
            printf(" EGADS Internal: Face %d - EdgeUV- = %d (EG_setTessFace)!\n",
                   index, stat);
          stat = findPoint(range, uvm, pt, pi, len, table, uv, uvx);
          if (stat < EGADS_SUCCESS) {
            printf(" EGADS Error: Face %d - FindPt- %d/%d = %d (EG_setTessFace)!\n",
                   index, pt, pi, stat);
            if (intEdg != NULL) EG_free(intEdg);
            EG_free(frlps);
            EG_free(uvs);
            EG_free(segs);
            EG_free(table);
            return stat;
          }
          map[stat]           =  np;
          if (intEdg == NULL) {
            uvs[2*np+2]       =  uvx[0];
            uvs[2*np+3]       =  uvx[1];
          } else {
            uvs[2*np+2]       =  intEdg[4*m  ] + intEdg[4*m+2];
            uvs[2*np+3]       =  intEdg[4*m+1] + intEdg[4*m+3];
          }
          segs[np].indices[0] =  np+1;
          segs[np].indices[1] =  np+2;
          segs[np].neighbor   = -iedge;
          segs[np].edge       =  senses[n]*lor*iedge;
          segs[np].index      =  m+1;
        }
      }
      if (intEdg != NULL) EG_free(intEdg);
    }
    if (np > 0) segs[np-1].indices[1] = st+1;
    frlps[i] = np - st;
  }
  nseg = np;

  /* fill up 2D tess structure */

  for (i = 0; i < len; i++)
    if (map[i] == -1) {
      map[i] = np;
      np++;
    }

  /* get the frame */
  fast.pts   = NULL;
  fast.segs  = NULL;
  fast.front = NULL;
  frame = (int *) EG_alloc(3*ntrix*sizeof(int));
  if (frame == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d frame (EG_setTessFace)!\n",
             index, ntrix);
    EG_free(frlps);
    EG_free(uvs);
    EG_free(segs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  n = EG_fillArea(nloop, frlps, uvs, frame, &nf8, 0, &fast);
  /* adjust for figure 8 configurations */
  if (nf8 != 0) {
    printf(" EGADS Warning: Face %d -> Found %d figure 8's!\n", index, nf8);
    for (i = 0; i < nf8; i++) if (n+2*i == ntrix) ntrix = n;
  }
  if (n != ntrix) {
    range[0] = range[2] = uvs[2];
    range[1] = range[3] = uvs[3];
    for (i = 2; i <= ntot; i++) {
      if (uvs[2*i  ] < range[0]) range[0] = uvs[2*i  ];
      if (uvs[2*i+1] < range[1]) range[1] = uvs[2*i+1];
      if (uvs[2*i  ] > range[2]) range[2] = uvs[2*i  ];
      if (uvs[2*i+1] > range[3]) range[3] = uvs[2*i+1];
    }
    for (i = 1; i <= ntot; i++) {
      uvs[2*i  ] = (uvs[2*i  ]-range[0])/(range[2]-range[0]);
      uvs[2*i+1] = (uvs[2*i+1]-range[1])/(range[3]-range[1]);
    }
    for (j = 0; j < 3; j++) {
      for (i = 1; i <= ntot; i++) {
        uvs[2*i  ] *= scl[j][0];
        uvs[2*i+1] *= scl[j][1];
      }
      n = EG_fillArea(nloop, frlps, uvs, frame, &nf8, 1, &fast);
      printf(" EGADS Internal: Face %d -> Renormalizing %d, ntris = %d (%d)!\n",
             index, j, ntrix, n);
      if (n == ntrix) break;
    }
  }
  if (fast.segs  != NULL) EG_free(fast.segs);
  if (fast.pts   != NULL) EG_free(fast.pts);
  if (fast.front != NULL) EG_free(fast.front);
  EG_free(uvs);
  if (n != ntrix) {
    printf(" EGADS Error: Face %d - Can't Triangulate Frame (EG_setTessFace)!\n",
           index);
    EG_free(frame);
    EG_free(frlps);
    EG_free(segs);
    EG_free(table);
    return EGADS_DEGEN;
  }

  /* set the triangle data */
  trix = (int *) EG_alloc(3*ntri*sizeof(int));
  tric = (int *) EG_alloc(3*ntri*sizeof(int));
  if ((trix == NULL) || (tric == NULL)) {
    if (trix != NULL) EG_free(trix);
    if (tric != NULL) EG_free(tric);
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d tris (EG_setTessFace)!\n",
             index, ntri);
    EG_free(frame);
    EG_free(frlps);
    EG_free(segs);
    EG_free(table);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ntri; i++) {
    trix[3*i  ] = map[tris[3*i  ]-1] + 1;
    trix[3*i+1] = map[tris[3*i+1]-1] + 1;
    trix[3*i+2] = map[tris[3*i+2]-1] + 1;
    tric[3*i  ] = i+1;
    tric[3*i+1] = i+1;
    tric[3*i+2] = i+1;
  }
  stat = makeNeighbors(index, len, ntri, trix, tric, nseg, segs);
  EG_free(segs);
  if (stat != EGADS_SUCCESS) {
    EG_free(trix);
    EG_free(tric);
    EG_free(frame);
    EG_free(frlps);
    EG_free(segs);
    EG_free(table);
    return stat;
  }

  /* set the reordered vertices */
  ptype  = (int *)    EG_alloc(  len*sizeof(int));
  pindex = (int *)    EG_alloc(  len*sizeof(int));
  uvs    = (double *) EG_alloc(2*len*sizeof(double));
  xyzs   = (double *) EG_alloc(3*len*sizeof(double));
  if ((ptype == NULL) || (pindex == NULL) || (xyzs == NULL) || (uvs == NULL)) {
    if (ptype  != NULL) EG_free(ptype);
    if (pindex != NULL) EG_free(pindex);
    if (uvs    != NULL) EG_free(uvs);
    if (xyzs   != NULL) EG_free(xyzs);
    if (outLevel > 0)
      printf(" EGADS Error: Face %d - Allocating %d verts (EG_setTessFace)!\n",
             index, len);
    EG_free(trix);
    EG_free(tric);
    EG_free(frame);
    EG_free(frlps);
    EG_free(table);
    return EGADS_MALLOC;
  }
  for (j = 0; j < len; j++) {
    i           = map[j];
    ptype[i]    = table[2*j  ];
    pindex[i]   = table[2*j+1];
    uvs[2*i  ]  = uv[2*j  ];
    uvs[2*i+1]  = uv[2*j+1];
    xyzs[3*i  ] = xyz[3*j  ];
    xyzs[3*i+1] = xyz[3*j+1];
    xyzs[3*i+2] = xyz[3*j+2];
  }
  EG_free(table);
  
  /* are we OK with the Frame? */
  for (i = 0; i < ntri; i++)
    for (j = 0; j < 3; j++)
      if (tric[3*i+j] == 0) {
        printf(" Face %d: tri = %d, side = %d -- No Neighbor!\n",
               index, i+1, j);
      } else if (tric[3*i+j] < 0) {
        mm = trix[3*i+sides[j][0]]-1;
        mp = trix[3*i+sides[j][1]]-1;
        if (ptype[mm] < 0)
          printf(" Face %d: Edge = %d  verts = %d %d -> Not in Frame-!\n",
                 index, -tric[3*i+j], mm+1, mp+1);
        if (ptype[mp] < 0)
          printf(" Face %d: Edge = %d  verts = %d %d -> Not in Frame+!\n",
                 index, -tric[3*i+j], mm+1, mp+1);
      }

  /* update the Face pointers */
  if (btess->tess2d[index-1].xyz    != NULL)
    EG_free(btess->tess2d[index-1].xyz);
  if (btess->tess2d[index-1].uv     != NULL)
    EG_free(btess->tess2d[index-1].uv);
  if (btess->tess2d[index-1].ptype  != NULL)
    EG_free(btess->tess2d[index-1].ptype);
  if (btess->tess2d[index-1].pindex != NULL)
    EG_free(btess->tess2d[index-1].pindex);
  if (btess->tess2d[index-1].bary   != NULL)
    EG_free(btess->tess2d[index-1].bary);
  if (btess->tess2d[index-1].frame  != NULL)
    EG_free(btess->tess2d[index-1].frame);
  if (btess->tess2d[index-1].frlps  != NULL)
    EG_free(btess->tess2d[index-1].frlps);
  if (btess->tess2d[index-1].tris   != NULL)
    EG_free(btess->tess2d[index-1].tris);
  if (btess->tess2d[index-1].tric   != NULL)
    EG_free(btess->tess2d[index-1].tric);
  btess->tess2d[index-1].npts   = len;
  btess->tess2d[index-1].xyz    = xyzs;
  btess->tess2d[index-1].uv     = uvs;
  btess->tess2d[index-1].ptype  = ptype;
  btess->tess2d[index-1].pindex = pindex;
  btess->tess2d[index-1].ntris  = ntri;
  btess->tess2d[index-1].tris   = trix;
  btess->tess2d[index-1].tric   = tric;
  btess->tess2d[index-1].bary   = NULL;
  btess->tess2d[index-1].nframe = ntrix;
  btess->tess2d[index-1].frame  = frame;
  btess->tess2d[index-1].frlps  = frlps;
  btess->tess2d[index-1].nfrlps = nloop;

  for (i = 1; i < nloop; i++)
    btess->tess2d[index-1].frlps[i] += btess->tess2d[index-1].frlps[i-1];

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_localToGlobal(const egObject *tess, int index, int local, int *global)
{
  int      stat;
  egTessel *btess;

  if  (tess == NULL)                          return EGADS_NULLOBJ;
  if  (tess->magicnumber != MAGIC)            return EGADS_NOTOBJ;
  if  (tess->oclass != TESSELLATION)          return EGADS_NOTTESS;
  if  (index == 0)                            return EGADS_INDEXERR;
  if  (local <  1)                            return EGADS_RANGERR;
  btess = (egTessel *) tess->blind;
  if (btess == NULL)                          return EGADS_NOTFOUND;
  if (btess->done == 0)                       return EGADS_TESSTATE;
  if ((index < 0) && (-index > btess->nEdge)) return EGADS_INDEXERR;
  if ((index > 0) && ( index > btess->nFace)) return EGADS_INDEXERR;

  if (btess->globals == NULL) {
    stat = EG_computeTessMap(btess, EG_outLevel(tess));
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_computeTessMap = %d (EG_localToGlobal)!\n", stat);
      return stat;
    }
  }

  if (index < 0) {
    if (btess->tess1d[-index-1].global == NULL) return EGADS_DEGEN;
    if (local > btess->tess1d[-index-1].npts) return EGADS_RANGERR;
    *global = btess->tess1d[-index-1].global[local-1];
  } else {
    if (btess->tess2d[ index-1].global == NULL) return EGADS_DEGEN;
    if (local > btess->tess2d[ index-1].npts) return EGADS_RANGERR;
    *global = btess->tess2d[ index-1].global[local-1];
  }

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_getGlobal(const egObject *tess, int global, int *ptype, int *pindex,
             /*@null@*/ double *xyz)
{
  int      i, j, stat;
  egTessel *btess;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  btess = (egTessel *) tess->blind;
  if (btess == NULL)                return EGADS_NOTFOUND;
  if (btess->done == 0)             return EGADS_TESSTATE;

  if (btess->globals == NULL) {
    stat = EG_computeTessMap(btess, EG_outLevel(tess));
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_computeTessMap = %d (EG_getGlobal)!\n", stat);
      return stat;
    }
  }
  if ((global < 1) || (global > btess->nGlobal)) return EGADS_INDEXERR;

  *ptype  = i = btess->globals[2*global-2];
  *pindex = j = btess->globals[2*global-1];
  if (xyz == NULL) return EGADS_SUCCESS;

  if (i == 0) {
    xyz[0] = btess->xyzs[3*j-3];
    xyz[1] = btess->xyzs[3*j-2];
    xyz[2] = btess->xyzs[3*j-1];
  } else if (i > 0) {
    xyz[0] = btess->tess1d[j-1].xyz[3*i-3];
    xyz[1] = btess->tess1d[j-1].xyz[3*i-2];
    xyz[2] = btess->tess1d[j-1].xyz[3*i-1];
  } else {
    xyz[0] = btess->tess2d[j-1].xyz[-3*i-3];
    xyz[1] = btess->tess2d[j-1].xyz[-3*i-2];
    xyz[2] = btess->tess2d[j-1].xyz[-3*i-1];
  }

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_findMidSide(int i1, int i2, int *table, midside *mid)
{
  int index, last;

  if (table[i1] == NOTFILLED) return 0;

  index = last = table[i1];
  while (index != NOTFILLED) {
    if (mid[index].vert2 == i2) return mid[index].nvert;
    last  = index;
    index = mid[index].next;
  }

  return -last-1;
}


__HOST_AND_DEVICE__ void
EG_getInterior(const ego face, double *xyz, double *uv)
{
  int    count, stat, diverge = 0;
  double a00, a10, a11, b0, b1, det, dist, ldist, frac = 1.0;
  double dx[3], du[2], range[4], result[18];

  stat = EG_getRange(face, range, &count);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Warning: EG_getRange = %d (EG_getInterior)!\n", stat);
    return;
  }
  if ((uv[0] < range[0]) || (uv[0] > range[1]) || (uv[1] < range[2]) ||
      (uv[1] > range[3])) {
    printf(" EGADS Warning: Out of Range (EG_getInterior)!\n");
    return;
  }

  /* newton iteration */
  ldist = 0.0;
  for (count = 0; count < 15; count++) {
    stat = EG_evaluate(face, uv, result);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Warning: EG_evaluate = %d (EG_getInterior)!\n", stat);
      return;
    }
    dx[0] = result[0] - xyz[0];
    dx[1] = result[1] - xyz[1];
    dx[2] = result[2] - xyz[2];
    dist  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    if (dist < EPS10) return;

    b0  =    -dx[0]*result[ 3] -     dx[1]*result[ 4] -     dx[2]*result[ 5];
    b1  =    -dx[0]*result[ 6] -     dx[1]*result[ 7] -     dx[2]*result[ 8];
    a00 = result[3]*result[ 3] + result[4]*result[ 4] + result[5]*result[ 5] +
              dx[0]*result[ 9] +     dx[1]*result[10] +     dx[2]*result[11];
    a10 = result[3]*result[ 6] + result[4]*result[ 7] + result[5]*result[ 8] +
              dx[0]*result[12] +     dx[1]*result[13] +     dx[2]*result[14];
    a11 = result[6]*result[ 6] + result[7]*result[ 7] + result[8]*result[ 8] +
              dx[0]*result[15] +     dx[1]*result[16] +     dx[2]*result[17];

    det = a00*a11 - a10*a10;
    if (det == 0.0) {
      printf(" EGADS Warning: det = zero (EG_getInterior)!\n");
      return;
    }
    det   = 1.0/det;
    du[0] = det*(b0*a11 - b1*a10);
    du[1] = det*(b1*a00 - b0*a10);
    if ((fabs(du[0]) < EPS10) && (fabs(du[1]) < EPS10)) return;
    if (count != 0)
      if (ldist < dist) {
        printf(" EG_getInterior -- diverge %d  %le %le!\n",
               count+1, du[0], du[1]);
/*      return;  */
        diverge++;
      }
    ldist  = dist;

    /* partial update? */
    if (uv[0]+frac*du[0] < range[0]) frac = (uv[0]-range[0])/fabs(du[0]);
    if (uv[0]+frac*du[0] > range[1]) frac = (range[1]-uv[0])/fabs(du[0]);
    if (uv[1]+frac*du[1] < range[2]) frac = (uv[1]-range[2])/fabs(du[1]);
    if (uv[1]+frac*du[1] > range[3]) frac = (range[3]-uv[1])/fabs(du[1]);

    uv[0] += frac*du[0];
    uv[1] += frac*du[1];
    if (frac != 1.0) {
      if (diverge != 0) printf(" EG_getInterior: boundary exit %d\n", diverge);
      return;
    }
  }

  printf(" EG_getInterior: not converged %d!\n", diverge);
}


__HOST_AND_DEVICE__ int
EG_baryInsert(ego face, double w1, double w2, double w3, const double *uvA,
              const double *uvB, const double *uvC, double *uvOUT)
{
  int    i, it, nT = 100, nS = 20, stat, fold = 0, periodic;
  double pA[18], pB[18], pC[18], OA[3], OB[3], OC[3], pIT[18], J[16], JI[16];
  double crA[3]  = {0.,0.,0.}, crB[3] = {0.,0.,0.}, crC[3] = {0.,0.,0.};
  double duA[3]  = {0.,0.,0.}, duB[3] = {0.,0.,0.}, duC[3] = {0.,0.,0.};
  double dvA[3]  = {0.,0.,0.}, dvB[3] = {0.,0.,0.}, dvC[3] = {0.,0.,0.};
  double crIT[3] = {0.,0.,0.}, AC[3], CB[3], BA[3], range[4], uvIT[4], uv[4];
  double grad[3], delta[4] = {0.,0.,0.,0.}, L[4] = {0.,0.,0.,0.};
  double a, b, c, s, da = 0, db = 0.0, dc = 0.0, det, ka, kc, delnrm = 1.0;
  double rsdnrm = 0.0, x3 = 0.0, tol = EPS10;
#ifdef FACE_NORMAL_FOLD_CHECK
  int    mtype;
#endif
#ifdef DEBUG
  double d, e1 = 0.0, e2 = 0.0;
#endif
#ifdef DEBUGG
  double x0 = 0.0, x1 = 0.0;
#endif

  /* Initial guess uv = w1*uvA + w2*uvB + w3*uvC */
  uvIT[0] = w1 * uvA[0] + w2 * uvB[0] + w3 * uvC[0];
  uvIT[1] = w1 * uvA[1] + w2 * uvB[1] + w3 * uvC[1];
  uvIT[2] = uvIT[3] = 1.0;
  ka      = (w1 * w1) / (w2 * w2);
  kc      = (w3 * w3) / (w2 * w2);

  stat    = EG_evaluate(face, uvA , pA);
  stat   += EG_evaluate(face, uvB , pB);
  stat   += EG_evaluate(face, uvC , pC);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_baryInsert: EG_evaluate triangle vertices stat %d !!\n", stat);
    return stat;
  }

#ifdef DEBUG
  EG_evaluate(face, uvIT , pIT);
  OA[0] = pA[0] - pIT[0];
  OA[1] = pA[1] - pIT[1];
  OA[2] = pA[2] - pIT[2];

  OB[0] = pB[0] - pIT[0];
  OB[1] = pB[1] - pIT[1];
  OB[2] = pB[2] - pIT[2];

  OC[0] = pC[0] - pIT[0];
  OC[1] = pC[1] - pIT[1];
  OC[2] = pC[2] - pIT[2];
  /* get areas */
  CRXSS(OB, OC, crA);
  a = sqrt(DOT(crA, crA));
  CRXSS(OC, OA, crB);
  b = sqrt(DOT(crB, crB));
  CRXSS(OA, OB, crC);
  c = sqrt(DOT(crC, crC));
  d = a + b + c;
  printf(" Initial areas TOT %lf  a %lf  b  %lf  c %lf\n", d, a, b, c);
  printf(" ratios %lf %lf %lf (expected %lf %lf %lf)\n", a/d, b/d, c/d,
         w1, w2, w3);
#endif

  CB[0] = pB[0] - pC[0]; CB[1] = pB[1] - pC[1]; CB[2] = pB[2] - pC[2];
  AC[0] = pC[0] - pA[0]; AC[1] = pC[1] - pA[1]; AC[2] = pC[2] - pA[2];
  BA[0] = pA[0] - pB[0]; BA[1] = pA[1] - pB[1]; BA[2] = pA[2] - pB[2];

#ifdef FACE_NORMAL_FOLD_CHECK
  /* get mtype=SFORWARD or mtype=SREVERSE for the face to get topology normal */
  mtype = face->mtype;
#else
  /* use the original triangle normal as the reference */
  CRXSS(AC, BA, crIT);
#endif

  for (it = 0; it < nT; it++) {

    /* perform a line search such that the dot of all triangles is positive
     * also check that the residual decreases from the previous iteration */
    s = 1.0;
    for (i = 0; i < nS; i++) {
      uv[0] = uvIT[0] - s*delta[0];
      uv[1] = uvIT[1] - s*delta[1];
      uv[2] = uvIT[2] - s*delta[2];
      uv[3] = uvIT[3] - s*delta[3];

      stat = EG_evaluate(face, uv, pIT);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_baryInsert: EG_evaluate = %d\n", stat);
        s /= 2.0;
        continue;
      }

      /* Get vectors */
      OA[0] = pA[0] - pIT[0]; OA[1] = pA[1] - pIT[1]; OA[2] = pA[2] - pIT[2];
      OB[0] = pB[0] - pIT[0]; OB[1] = pB[1] - pIT[1]; OB[2] = pB[2] - pIT[2];
      OC[0] = pC[0] - pIT[0]; OC[1] = pC[1] - pIT[1]; OC[2] = pC[2] - pIT[2];
      /* get areas */
      CRXSS(OB, OC, crA);
      CRXSS(OC, OA, crB);
      CRXSS(OA, OB, crC);

#ifdef FACE_NORMAL_FOLD_CHECK
      /* get the normal vector at the proposed point */
      CRXSS(pIT+3, pIT+6, crIT);
      crIT[0] *= mtype;
      crIT[1] *= mtype;
      crIT[2] *= mtype;
#endif

      /* if any of triangles are folded, shorten the line search
       * but only if a non-folded solution was found */
      if ((it > 0) &&
          ((DOT(crA,crIT) < 0.0) || (DOT(crB,crIT) < 0.0) || (DOT(crC,crIT) < 0.0))) {
        s /= 2.0;
        continue;
      } else if (it == 0) {
        /* check if the initial guess is folded */
        if ((DOT(crA,crIT) < 0.0) || (DOT(crB,crIT) < 0.0) || (DOT(crC,crIT) < 0.0)) {
          if (fold == 1) {
            /* already tried better guess with inverse evaluate... */
            printf(" EG_baryInsert: Initial Fold!\n");
#ifdef DEBUGM
            printf("VARIABLES=X, Y, Z, U, V\n");
            printf("ZONE N=4, E=3, ET=TRIANGLE F=FEPOINT\n");
            printf("%f, %f, %f, %f, %f\n", pA[0], pA[1], pA[2], uvA[0], uvA[1]);
            printf("%f, %f, %f, %f, %f\n", pB[0], pB[1], pB[2], uvB[0], uvB[1]);
            printf("%f, %f, %f, %f, %f\n", pC[0], pC[1], pC[2], uvC[0], uvC[1]);
            printf("%f, %f, %f, %f, %f\n", pIT[0], pIT[1], pIT[2], uvIT[0], uvIT[1]);
            printf("1 2 4\n");
            printf("2 3 4\n");
            printf("3 1 4\n");
            printf("ZONE N=3, E=1, ET=TRIANGLE F=FEPOINT\n");
            printf("%f, %f, %f, %f, %f\n", pA[0], pA[1], pA[2], uvA[0], uvA[1]);
            printf("%f, %f, %f, %f, %f\n", pB[0], pB[1], pB[2], uvB[0], uvB[1]);
            printf("%f, %f, %f, %f, %f\n", pC[0], pC[1], pC[2], uvC[0], uvC[1]);
            printf("1 2 3\n");
#endif
            return EGADS_TESSTATE;
          }

          /* use inverse evaluate to get a better initial guess */
          pIT[0] = w1 * pA[0] + w2 * pB[0] + w3 * pC[0];
          pIT[1] = w1 * pA[1] + w2 * pB[1] + w3 * pC[1];
          pIT[2] = w1 * pA[2] + w2 * pB[2] + w3 * pC[2];
          stat = EG_invEvaluateGuess(face, pIT, uvIT, J);

          /* check for success and in range */
          EG_getRange(face, range, &periodic);
          if ((stat != EGADS_SUCCESS) ||
              (uvIT[0] < range[0]) || (uvIT[0] > range[1]) ||
              (uvIT[1] < range[2]) || (uvIT[1] > range[3])) {

            if (stat != EGADS_SUCCESS) {
              printf(" EG_baryInsert: EG_invEvaluateGuess = %d\n", stat);
            } else {
              printf(" EG_baryInsert: EG_invEvaluateGuess out of range!\n");
            }

            /* one last hope for a better initial guess */
            stat = EG_invEvaluate(face, pIT, uvIT, J);
            if (stat != EGADS_SUCCESS) {
              printf(" EG_baryInsert: EG_invEvaluate = %d\n", stat);
              return EGADS_TESSTATE;
            }
          }
          fold = 1;
          continue; /* try again */
        }
      }

      da = 1.0 + uv[2];
      db = 1.0 - ka * uv[2] - kc * uv[3];
      dc = 1.0 + uv[3];
      a  = DOT(crA, crA);
      b  = DOT(crB, crB);
      c  = DOT(crC, crC);
      /* Evaluate fn */
      /* a = w1 (a + b + c)
       * b = w2 (a + b + c)
       * c = w3 (a + b + c)
       * a + b + c = b / w2
       * a = w1/w2 * b; c = w3 / w2 *  b */
      grad[0] = pIT[3]; grad[1] = pIT[4]; grad[2] = pIT[5];
      CRXSS(grad, CB, duA);
      CRXSS(grad, AC, duB);
      CRXSS(grad, BA, duC);
      grad[0] = pIT[6]; grad[1] = pIT[7]; grad[2] = pIT[8];
      CRXSS(grad, CB, dvA);
      CRXSS(grad, AC, dvB);
      CRXSS(grad, BA, dvC);
      L[0] = DOT(duA, crA) * da + DOT(duB, crB) * db +
             DOT(duC, crC) * dc;
      L[1] = DOT(dvA, crA) * da + DOT(dvB, crB) * db +
             DOT(dvC, crC) * dc;
      L[2] = 0.5 * (a - ka * b);
      L[3] = 0.5 * (c - kc * b);

      /* check for convergence on the residual */
      rsdnrm  = sqrt(DOT4(L, L));

      if (it == 0 || rsdnrm < x3) {
        x3 = rsdnrm; /* save off the last resdual and exit */
        break;
      }  else {
        s /= 2.0; /* try again if the residual grew */
        continue;
      }
    }
    if (i == nS) printf(" EG_baryInsert: LineSearch Failed!\n");

#ifdef DEBUGG
    if      (it == 0) x0 = rsdnrm;
    else if (it == 1) x1 = rsdnrm;
    else {
      e1 = fabs(x1 / x0);
      e2 = fabs(rsdnrm / x1);
      x0 = x1;
      x1 = rsdnrm;
    }
    
    printf(" Xn = (%lf %lf %lf %lf )\n", uvIT[0], uvIT[1], uvIT[2], uvIT[3]);
    printf(" L = %lf %lf %lf %lf DELTA %lf %lf %lf %lf SIZE  %1.8le < %1.8le\n",
           L[0], L[1], L[2], L[3], delta[0], delta[1], delta[2], delta[3],
           rsdnrm, tol);
#endif

    /* update the solution */
    uvIT[0] = uv[0];
    uvIT[1] = uv[1];
    uvIT[2] = uv[2];
    uvIT[3] = uv[3];

    if (i      == nS) break;  /* line search failed */
    if (rsdnrm < tol) break;  /* converged! */

    J[ 2] = DOT(duA, crA) - DOT(duB, crB) * ka;
    J[ 3] = DOT(duC, crC) - DOT(duB, crB) * kc;
    J[ 6] = DOT(dvA, crA) - DOT(dvB, crB) * ka;
    J[ 7] = DOT(dvC, crC) - DOT(dvB, crB) * kc;
    J[ 8] = J[ 2];
    J[12] = J[ 3];
    J[ 9] = J[ 6];
    J[13] = J[ 7];
    J[10] = J[11] = 0.0;
    J[14] = J[15] = 0.0;
    J[ 0] = DOT(duA, duA) * da + DOT(duB, duB) * db +
            DOT(duC, duC) * dc;
    J[ 1] = DOT(duA, dvA) * da + DOT(duB, dvB) * db +
            DOT(duC, dvC) * dc;
    J[ 5] = DOT(dvA, dvA) * da + DOT(dvB, dvB) * db +
            DOT(dvC, dvC) * dc;
    grad[0] = pIT[9]; grad[1] = pIT[10]; grad[2] = pIT[11];
    CRXSS(grad, CB, duA);
    CRXSS(grad, AC, duB);
    CRXSS(grad, BA, duC);
    J[0]   += DOT(duA, crA) * da + DOT(duB, crB) * db +
              DOT(duC, crC) * dc;
    grad[0] = pIT[12]; grad[1] = pIT[13]; grad[2] = pIT[14];
    CRXSS(grad, CB, duA);
    CRXSS(grad, AC, duB);
    CRXSS(grad, BA, duC);
    J[1]   += DOT(duA, crA) * da + DOT(duB, crB) * db +
              DOT(duC, crC) * dc;
    grad[0] = pIT[15]; grad[1] = pIT[16]; grad[2] = pIT[17];
    CRXSS(grad, CB, duA);
    CRXSS(grad, AC, duB);
    CRXSS(grad, BA, duC);
    J[5] += DOT(duA, crA) * da + DOT(duB, crB) * db +
            DOT(duC, crC) * dc;
    J[4]  = J[1];
#ifdef DEBUGG
    printf("\n JACOBIAN MATRIX 4 x 4 \n");
    printf(" %lf %lf %lf %lf\n", J[ 0], J[ 1], J[ 2], J[ 3]);
    printf(" %lf %lf %lf %lf\n", J[ 4], J[ 5], J[ 6], J[ 7]);
    printf(" %lf %lf %lf %lf\n", J[ 8], J[ 9], J[10], J[11]);
    printf(" %lf %lf %lf %lf\n", J[12], J[13], J[14], J[15]);
    printf(" ----------------------\n");
#endif
    /* Directly invert 4x4 matrix */
    JI[0 ] = J[6] * J[11] * J[13] - J[ 7] * J[10] * J[13] + J[ 7] * J[ 9] * J[14] -
             J[5] * J[11] * J[14] - J[ 6] * J[ 9] * J[15] + J[ 5] * J[10] * J[15];
    JI[1 ] = J[3] * J[10] * J[13] - J[ 2] * J[11] * J[13] - J[ 3] * J[ 9] * J[14] +
             J[1] * J[11] * J[14] + J[ 2] * J[ 9] * J[15] - J[ 1] * J[10] * J[15];
    JI[2 ] = J[2] * J[ 7] * J[13] - J[ 3] * J[ 6] * J[13] + J[ 3] * J[ 5] * J[14] -
             J[1] * J[ 7] * J[14] - J[ 2] * J[ 5] * J[15] + J[ 1] * J[ 6] * J[15];
    JI[3 ] = J[3] * J[ 6] * J[ 9] - J[ 2] * J[ 7] * J[ 9] - J[ 3] * J[ 5] * J[10] +
             J[1] * J[ 7] * J[10] + J[ 2] * J[ 5] * J[11] - J[ 1] * J[ 6] * J[11];
    JI[4 ] = J[7] * J[10] * J[12] - J[ 6] * J[11] * J[12] - J[ 7] * J[ 8] * J[14] +
             J[4] * J[11] * J[14] + J[ 6] * J[ 8] * J[15] - J[ 4] * J[10] * J[15];
    JI[5 ] = J[2] * J[11] * J[12] - J[ 3] * J[10] * J[12] + J[ 3] * J[ 8] * J[14] -
             J[0] * J[11] * J[14] - J[ 2] * J[ 8] * J[15] + J[ 0] * J[10] * J[15];
    JI[6 ] = J[3] * J[ 6] * J[12] - J[ 2] * J[ 7] * J[12] - J[ 3] * J[4] * J[14] +
             J[0] * J[ 7] * J[14] + J[ 2] * J[ 4] * J[15] - J[ 0] * J[6] * J[15];
    JI[7 ] = J[2] * J[ 7] * J[ 8] - J[ 3] * J[ 6] * J[ 8] + J[ 3] * J[4] * J[10] -
             J[0] * J[ 7] * J[10] - J[ 2] * J[ 4] * J[11] + J[ 0] * J[6] * J[11];
    JI[8 ] = J[5] * J[11] * J[12] - J[ 7] * J[ 9] * J[12] + J[ 7] * J[8] * J[13] -
             J[4] * J[11] * J[13] - J[ 5] * J[ 8] * J[15] + J[ 4] * J[9] * J[15];
    JI[9 ] = J[3] * J[ 9] * J[12] - J[ 1] * J[11] * J[12] - J[ 3] * J[8] * J[13] +
             J[0] * J[11] * J[13] + J[ 1] * J[ 8] * J[15] - J[ 0] * J[9] * J[15];
    JI[10] = J[1] * J[ 7] * J[12] - J[ 3] * J[ 5] * J[12] + J[ 3] * J[4] * J[13] -
             J[0] * J[ 7] * J[13] - J[ 1] * J[ 4] * J[15] + J[ 0] * J[5] * J[15];
    JI[11] = J[3] * J[ 5] * J[ 8] - J[ 1] * J[ 7] * J[ 8] - J[ 3] * J[4] * J[ 9] +
             J[0] * J[ 7] * J[ 9] + J[ 1] * J[ 4] * J[11] - J[ 0] * J[5] * J[11];
    JI[12] = J[6] * J[ 9] * J[12] - J[ 5] * J[10] * J[12] - J[ 6] * J[8] * J[13] +
             J[4] * J[10] * J[13] + J[ 5] * J[ 8] * J[14] - J[ 4] * J[9] * J[14];
    JI[13] = J[1] * J[10] * J[12] - J[ 2] * J[ 9] * J[12] + J[ 2] * J[8] * J[13] -
             J[0] * J[10] * J[13] - J[ 1] * J[ 8] * J[14] + J[ 0] * J[9] * J[14];
    JI[14] = J[2] * J[ 5] * J[12] - J[ 1] * J[ 6] * J[12] - J[ 2] * J[4] * J[13] +
             J[0] * J[ 6] * J[13] + J[ 1] * J[ 4] * J[14] - J[ 0] * J[5] * J[14];
    JI[15] = J[1] * J[ 6] * J[ 8] - J[ 2] * J[ 5] * J[ 8] + J[ 2] * J[4] * J[ 9] -
             J[0] * J[ 6] * J[ 9] - J[ 1] * J[ 4] * J[10] + J[ 0] * J[5] * J[10];
    det    = J[3] * J[ 6] * J[ 9] * J[12] - J[ 2] * J[ 7] * J[ 9] * J[12] -
             J[3] * J[ 5] * J[10] * J[12] + J[ 1] * J[ 7] * J[10] * J[12] +
             J[2] * J[ 5] * J[11] * J[12] - J[ 1] * J[ 6] * J[11] * J[12] -
             J[3] * J[ 6] * J[ 8] * J[13] + J[ 2] * J[ 7] * J[ 8] * J[13] +
             J[3] * J[ 4] * J[10] * J[13] - J[ 0] * J[ 7] * J[10] * J[13] -
             J[2] * J[ 4] * J[11] * J[13] + J[ 0] * J[ 6] * J[11] * J[13] +
             J[3] * J[ 5] * J[ 8] * J[14] - J[ 1] * J[ 7] * J[ 8] * J[14] -
             J[3] * J[ 4] * J[ 9] * J[14] + J[ 0] * J[ 7] * J[ 9] * J[14] +
             J[1] * J[ 4] * J[11] * J[14] - J[ 0] * J[ 5] * J[11] * J[14] -
             J[2] * J[ 5] * J[ 8] * J[15] + J[ 1] * J[ 6] * J[ 8] * J[15] +
             J[2] * J[ 4] * J[ 9] * J[15] - J[ 0] * J[ 6] * J[ 9] * J[15] -
             J[1] * J[ 4] * J[10] * J[15] + J[ 0] * J[ 5] * J[10] * J[15];
    if (det == 0.0) break;
    det    = 1.0 / det;
#ifdef DEBUGG
    printf("\n INVERSE JACOBIAN MATRIX 4 x 4 DET %lf ORI DET %1.12le\n",
           det, 1.0/det);
    printf(" %lf %lf %lf %lf\n",JI[0]*det, JI[1]*det, JI[2]*det, JI[3]*det);
    printf(" %lf %lf %lf %lf\n",JI[4]*det, JI[5]*det, JI[6]*det, JI[7]*det);
    printf(" %lf %lf %lf %lf\n",JI[8]*det, JI[9]*det, JI[10]*det, JI[11]*det);
    printf(" %lf %lf %lf %lf\n",JI[12]*det, JI[13]*det, JI[14]*det, JI[15]*det);
    printf(" ----------------------\n");
#endif
    delta[0]  = det * (JI[0 ] * L[0] + JI[1 ] * L[1] + JI[2 ] * L[2] + JI[ 3] * L[3]);
    delta[1]  = det * (JI[4 ] * L[0] + JI[5 ] * L[1] + JI[6 ] * L[2] + JI[ 7] * L[3]);
    delta[2]  = det * (JI[8 ] * L[0] + JI[9 ] * L[1] + JI[10] * L[2] + JI[11] * L[3]);
    delta[3]  = det * (JI[12] * L[0] + JI[13] * L[1] + JI[14] * L[2] + JI[15] * L[3]);

    /* check for convergence on the parameter update */
    delnrm  = sqrt(DOT4(delta, delta));
    if (delnrm < tol) break;  /* converged! */
  }
  if (rsdnrm >= tol && delnrm >= tol) {
    printf(" EG_baryInsert: not converged -- residual %1.2le delta %1.2le (%1.2le)\n", rsdnrm, delnrm, tol);
  }

#ifdef DEBUG
  printf(" \n\n Found point %lf %lf %lf %lf %lf\n",
         pIT[0], pIT[1], pIT[2], uvIT[0], uvIT[1]);
  if (e1 > 0.0 && e2 > 0.0)
    printf(" Newton iterations %d Rate %lf\n", it, log(e2)/log(e1));
  OA[0] = pA[0] - pIT[0];
  OA[1] = pA[1] - pIT[1];
  OA[2] = pA[2] - pIT[2];

  OB[0] = pB[0] - pIT[0];
  OB[1] = pB[1] - pIT[1];
  OB[2] = pB[2] - pIT[2];

  OC[0] = pC[0] - pIT[0];
  OC[1] = pC[1] - pIT[1];
  OC[2] = pC[2] - pIT[2];
  /* get areas */
  CRXSS(OB, OC, crA);
  a = sqrt(DOT(crA, crA));
  CRXSS(OC, OA, crB);
  b = sqrt(DOT(crB, crB));
  CRXSS(OA, OB, crC);
  c = sqrt(DOT(crC, crC));
  d = a + b + c;
  printf(" Final areas TOT %lf  a %lf  b  %lf  c %lf\n", d, a, b, c);
  printf(" ratios %lf %lf %lf (expected %lf %lf %lf)\n",
         a/d, b/d, c/d, w1, w2, w3);
#endif

  uvOUT[0] = uvIT[0];
  uvOUT[1] = uvIT[1];

  if ((DOT(crA,crIT) < 0.0) || (DOT(crB,crIT) < 0.0) || (DOT(crC,crIT) < 0.0)) {
#ifdef DEBUGM
    printf("VARIABLES=X, Y, Z, U, V\n");
    printf("ZONE N=4, E=3, ET=TRIANGLE F=FEPOINT\n");
    printf("%f, %f, %f, %f, %f\n", pA[0], pA[1], pA[2], uvA[0], uvA[1]);
    printf("%f, %f, %f, %f, %f\n", pB[0], pB[1], pB[2], uvB[0], uvB[1]);
    printf("%f, %f, %f, %f, %f\n", pC[0], pC[1], pC[2], uvC[0], uvC[1]);
    printf("%f, %f, %f, %f, %f\n", pIT[0], pIT[1], pIT[2], uvIT[0], uvIT[1]);
    printf("1 2 4\n");
    printf("2 3 4\n");
    printf("3 1 4\n");
#endif
    return EGADS_TESSTATE;
  }

  return EGADS_SUCCESS;
}


/* P such that splits uv0 and uv2 in fact1, 1-fact1 and
               splits uv1 and uv3 in fact2, 1-fact2 and
               with minimum distance
*/
__HOST_AND_DEVICE__ void
EG_minArc4(const ego face, double fact1, double fact2, const double *uv0,
           const double *uv1, const double *uv2, const double *uv3,
           double *uvOUT)
{
  int    i, it, stat, nT = 100;
  double x2, l[4], p0[18], p2[18], p1[18], p3[18], pIT[18], J[16], JI[16], L[4];
  double Pu0[3], Pu1[3], Pu2[3], Pu3[3], Pv0[3], Pv1[3], Pv2[3], Pv3[3];
  double delta[4], det, k02, k13, d0, d1, d2, d3, aux0[3], aux1[3], aux2[3];
  double aux3[3], r0[3], r1[3], r2[3], r3[3], uvIT[4], range[4], pOUT[18];
 #ifdef DEBUG
  double e1 = 0.0, e2 = 0.0, lt02, lt13, x1, x0;
#endif

  stat = EG_getRange(face, range, &it);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_minArc4: EG_getRange = %d!\n", stat);
    return;
  }

  /* Initial guess uv = 0.5 (uv0 + uv1) */
  uvIT[0]  = 0.25*(uv0[0] + uv1[0] + uv2[0] + uv3[0]);
  uvIT[1]  = 0.25*(uv0[1] + uv1[1] + uv2[1] + uv3[1]);
  uvOUT[0] = uvIT[0];
  uvOUT[1] = uvIT[1];
  uvIT[2]  = uvIT[3] = 0.0;
  i        = EG_evaluate(face, uv0, p0);
  i       += EG_evaluate(face, uv1, p1);
  i       += EG_evaluate(face, uv2, p2);
  i       += EG_evaluate(face, uv3, p3);
  i       += EG_evaluate(face, uvIT, pIT);
#ifdef DEBUG
  printf("%lf %lf %lf 1\n",  p0[0],  p0[1],  p0[2]);
  printf("%lf %lf %lf 2\n",  p1[0],  p1[1],  p1[2]);
  printf("%lf %lf %lf 3\n",  p2[0],  p2[1],  p2[2]);
  printf("%lf %lf %lf 4\n",  p3[0],  p3[1],  p3[2]);
  printf("%lf %lf %lf 4\n", pIT[0], pIT[1], pIT[2]);
#endif
  if (i != EGADS_SUCCESS) {
    printf(" EG_minArc4: EG_evaluate %d!\n", i);
    return;
  }
  
#ifdef DEBUG
  r0[0] = pIT[0] - p0[0]; r0[1] = pIT[1] - p0[1]; r0[2] = pIT[2] - p0[2];
  r1[0] = pIT[0] - p1[0]; r1[1] = pIT[1] - p1[1]; r1[2] = pIT[2] - p1[2];
  r2[0] = pIT[0] - p2[0]; r2[1] = pIT[1] - p2[1]; r2[2] = pIT[2] - p2[2];
  r3[0] = pIT[0] - p3[0]; r3[1] = pIT[1] - p3[1]; r3[2] = pIT[2] - p3[2];
  l[0]  = DOT(r0, r0);    l[1]  = DOT(r1, r1);
  l[2]  = DOT(r2, r2);    l[3]  = DOT(r3, r3);
  l[0]  = sqrt(l[0]);
  l[1]  = sqrt(l[1]);
  l[2]  = sqrt(l[2]);
  l[3]  = sqrt(l[3]);
  lt02  = l[0] + l[2];
  printf(" Initial Guess As distances sqrt(l0) = %lf sqrt (l1) %lf\n",
         l[0], l[2]);
  printf(" ratios %lf %lf (EXPECTED %lf %lf)\n",
         l[0] / lt02, l[2]/lt02, fact1, (1.0 - fact1));
  printf(" As distances sqrt(l0) = %lf sqrt (l1) %lf\n",
           l[1], l[3]);
  printf(" ratios %lf %lf (EXPECTED %lf %lf)\n",
         l[1] / lt02, l[3]/lt02, fact2, (1.0 - fact2));
  printf(" --------------------------------------------------- \n");
#endif
  
  k02 = (fact1*fact1) / ((1.0 - fact1)*(1.0 - fact1)) ;
  k13 = (fact2*fact2) / ((1.0 - fact2)*(1.0 - fact2));
  for (it   = 0; it < nT; it++) {
      r0[0] = pIT[0] - p0[0]; r0[1] = pIT[1] - p0[1]; r0[2] = pIT[2] - p0[2];
      r1[0] = pIT[0] - p1[0]; r1[1] = pIT[1] - p1[1]; r1[2] = pIT[2] - p1[2];

      r2[0] = pIT[0] - p2[0]; r2[1] = pIT[1] - p2[1]; r2[2] = pIT[2] - p2[2];
      r3[0] = pIT[0] - p3[0]; r3[1] = pIT[1] - p3[1]; r3[2] = pIT[2] - p3[2];

      l [0] = DOT(r0, r0); l [1] = DOT(r1, r1);
      l [2] = DOT(r2, r2); l [3] = DOT(r3, r3);
#ifdef DEBUGG
        printf("%lf %lf %lf 1\n",  p0[0],  p0[1],  p0[2]);
        printf("%lf %lf %lf 2\n",  p1[0],  p1[1],  p1[2]);
        printf("%lf %lf %lf 3\n",  p2[0],  p2[1],  p2[2]);
        printf("%lf %lf %lf 4\n",  p3[0],  p3[1],  p3[2]);
        printf("%lf %lf %lf 4\n", pIT[0], pIT[1], pIT[2]);
#endif
      d0 = 1.0 + uvIT[2];
      d1 = 1.0 + uvIT[3];
      d2 = 1.0 - k02*uvIT[2];
      d3 = 1.0 - k13*uvIT[3];
      Pu0[0] = pIT[3]; Pu0[1] = pIT[4]; Pu0[2] = pIT[5];
      Pu1[0] = pIT[3]; Pu1[1] = pIT[4]; Pu1[2] = pIT[5];
      Pu2[0] = pIT[3]; Pu2[1] = pIT[4]; Pu2[2] = pIT[5];
      Pu3[0] = pIT[3]; Pu3[1] = pIT[4]; Pu3[2] = pIT[5];

      Pv0[0] = pIT[6]; Pv0[1] = pIT[7]; Pv0[2] = pIT[8];
      Pv1[0] = pIT[6]; Pv1[1] = pIT[7]; Pv1[2] = pIT[8];
      Pv2[0] = pIT[6]; Pv2[1] = pIT[7]; Pv2[2] = pIT[8];
      Pv3[0] = pIT[6]; Pv3[1] = pIT[7]; Pv3[2] = pIT[8];

      L[0]   = d0*DOT(Pu0, r0) + d1*DOT(Pu1, r1) +
               d2*DOT(Pu2, r2) + d3*DOT(Pu3, r3);
      L[1]   = d0*DOT(Pv0, r0) + d1*DOT(Pv1, r1) +
               d2*DOT(Pv2, r2) + d3*DOT(Pv3, r3);
      L[2]   = 0.5*(l[0] - k02*l[2]);
      L[3]   = 0.5*(l[1] - k13*l[3]);

      aux0[0] = pIT[9]; aux0[1] = pIT[10]; aux0[2] = pIT[11];
      aux1[0] = pIT[9]; aux1[1] = pIT[10]; aux1[2] = pIT[11];
      aux2[0] = pIT[9]; aux2[1] = pIT[10]; aux2[2] = pIT[11];
      aux3[0] = pIT[9]; aux3[1] = pIT[10]; aux3[2] = pIT[11];
      J[ 0]   = d0*(DOT(Pu0, Pu0) + DOT(aux0, r0)) +
                d1*(DOT(Pu1, Pu1) + DOT(aux1, r1)) +
                d2*(DOT(Pu2, Pu2) + DOT(aux2, r2)) +
                d3*(DOT(Pu3, Pu3) + DOT(aux3, r3)) ;
      aux0[0] = pIT[12]; aux0[1] = pIT[13]; aux0[2] = pIT[14];
      aux1[0] = pIT[12]; aux1[1] = pIT[13]; aux1[2] = pIT[14];
      aux2[0] = pIT[12]; aux2[1] = pIT[13]; aux2[2] = pIT[14];
      aux3[0] = pIT[12]; aux3[1] = pIT[13]; aux3[2] = pIT[14];
      J[ 1]   = d0*(DOT(Pu0, Pv0) + DOT(aux0, r0)) +
                d1*(DOT(Pu1, Pv1) + DOT(aux1, r1)) +
                d2*(DOT(Pu2, Pv2) + DOT(aux2, r2)) +
                d3*(DOT(Pu3, Pv3) + DOT(aux3, r3)) ;
      aux0[0] = pIT[15]; aux0[1] = pIT[16]; aux0[2] = pIT[17];
      aux1[0] = pIT[15]; aux1[1] = pIT[16]; aux1[2] = pIT[17];
      aux2[0] = pIT[15]; aux2[1] = pIT[16]; aux2[2] = pIT[17];
      aux3[0] = pIT[15]; aux3[1] = pIT[16]; aux3[2] = pIT[17];
      J[ 5]   = d0*(DOT(Pv0, Pv0) + DOT(aux0, r0)) +
                d1*(DOT(Pv1, Pv1) + DOT(aux1, r1)) +
                d2*(DOT(Pv2, Pv2) + DOT(aux2, r2)) +
                d3*(DOT(Pv3, Pv3) + DOT(aux3, r3)) ;
      J [ 2]  = DOT(Pu0, r0) - k02*DOT(Pu2, r2);
      J [ 3]  = DOT(Pu1, r1) - k13*DOT(Pu3, r3);
      J [ 6]  = DOT(Pv0, r0) - k02*DOT(Pv2, r2);
      J [ 7]  = DOT(Pv1, r1) - k13*DOT(Pv3, r3);
      J [ 4]  = J [ 1];
      J [ 8]  = J [ 2];
      J [ 9]  = J [ 6];
      J [12]  = J [ 3];
      J [13]  = J [ 7];
      J [10]  = J [11] = J [14] = J [15] = 0.0;
      JI[ 0] = J[6]*J[11]*J[13] - J[ 7]*J[10]*J[13] + J[ 7]*J[ 9]*J[14] -
               J[5]*J[11]*J[14] - J[ 6]*J[ 9]*J[15] + J[ 5]*J[10]*J[15];
      JI[ 1] = J[3]*J[10]*J[13] - J[ 2]*J[11]*J[13] - J[ 3]*J[ 9]*J[14] +
               J[1]*J[11]*J[14] + J[ 2]*J[ 9]*J[15] - J[ 1]*J[10]*J[15];
      JI[ 2] = J[2]*J[ 7]*J[13] - J[ 3]*J[ 6]*J[13] + J[ 3]*J[ 5]*J[14] -
               J[1]*J[ 7]*J[14] - J[ 2]*J[ 5]*J[15] + J[ 1]*J[ 6]*J[15];
      JI[ 3] = J[3]*J[ 6]*J[ 9] - J[ 2]*J[ 7]*J[ 9] - J[ 3]*J[ 5]*J[10] +
               J[1]*J[ 7]*J[10] + J[ 2]*J[ 5]*J[11] - J[ 1]*J[ 6]*J[11];
      JI[ 4] = J[7]*J[10]*J[12] - J[ 6]*J[11]*J[12] - J[ 7]*J[ 8]*J[14] +
               J[4]*J[11]*J[14] + J[ 6]*J[ 8]*J[15] - J[ 4]*J[10]*J[15];
      JI[ 5] = J[2]*J[11]*J[12] - J[ 3]*J[10]*J[12] + J[ 3]*J[ 8]*J[14] -
               J[0]*J[11]*J[14] - J[ 2]*J[ 8]*J[15] + J[ 0]*J[10]*J[15];
      JI[ 6] = J[3]*J[ 6]*J[12] - J[ 2]*J[ 7]*J[12] - J[ 3]*J[ 4]*J[14] +
               J[0]*J[ 7]*J[14] + J[ 2]*J[ 4]*J[15] - J[ 0]*J[ 6]*J[15];
      JI[ 7] = J[2]*J[ 7]*J[ 8] - J[ 3]*J[ 6]*J[ 8] + J[ 3]*J[ 4]*J[10] -
               J[0]*J[ 7]*J[10] - J[ 2]*J[ 4]*J[11] + J[ 0]*J[ 6]*J[11];
      JI[ 8] = J[5]*J[11]*J[12] - J[ 7]*J[ 9]*J[12] + J[ 7]*J[ 8]*J[13] -
               J[4]*J[11]*J[13] - J[ 5]*J[ 8]*J[15] + J[ 4]*J[ 9]*J[15];
      JI[ 9] = J[3]*J[ 9]*J[12] - J[ 1]*J[11]*J[12] - J[ 3]*J[ 8]*J[13] +
               J[0]*J[11]*J[13] + J[ 1]*J[ 8]*J[15] - J[ 0]*J[ 9]*J[15];
      JI[10] = J[1]*J[ 7]*J[12] - J[ 3]*J[ 5]*J[12] + J[ 3]*J[ 4]*J[13] -
               J[0]*J[ 7]*J[13] - J[ 1]*J[ 4]*J[15] + J[ 0]*J[ 5]*J[15];
      JI[11] = J[3]*J[ 5]*J[ 8] - J[ 1]*J[ 7]*J[ 8] - J[ 3]*J[ 4]*J[ 9] +
               J[0]*J[ 7]*J[ 9] + J[ 1]*J[ 4]*J[11] - J[ 0]*J[ 5]*J[11];
      JI[12] = J[6]*J[ 9]*J[12] - J[ 5]*J[10]*J[12] - J[ 6]*J[ 8]*J[13] +
               J[4]*J[10]*J[13] + J[ 5]*J[ 8]*J[14] - J[ 4]*J[ 9]*J[14];
      JI[13] = J[1]*J[10]*J[12] - J[ 2]*J[ 9]*J[12] + J[ 2]*J[ 8]*J[13] -
               J[0]*J[10]*J[13] - J[ 1]*J[ 8]*J[14] + J[ 0]*J[ 9]*J[14];
      JI[14] = J[2]*J[ 5]*J[12] - J[ 1]*J[ 6]*J[12] - J[ 2]*J[ 4]*J[13] +
               J[0]*J[ 6]*J[13] + J[ 1]*J[ 4]*J[14] - J[ 0]*J[ 5]*J[14];
      JI[15] = J[1]*J[ 6]*J[ 8] - J[ 2]*J[ 5]*J[ 8] + J[ 2]*J[ 4]*J[ 9] -
               J[0]*J[ 6]*J[ 9] - J[ 1]*J[ 4]*J[10] + J[ 0]*J[ 5]*J[10];
      det    = J[3]*J[ 6]*J[ 9]*J[12] - J[ 2]*J[ 7]*J[ 9]*J[12] -
               J[3]*J[ 5]*J[10]*J[12] + J[ 1]*J[ 7]*J[10]*J[12] +
               J[2]*J[ 5]*J[11]*J[12] - J[ 1]*J[ 6]*J[11]*J[12] -
               J[3]*J[ 6]*J[ 8]*J[13] + J[ 2]*J[ 7]*J[ 8]*J[13] +
               J[3]*J[ 4]*J[10]*J[13] - J[ 0]*J[ 7]*J[10]*J[13] -
               J[2]*J[ 4]*J[11]*J[13] + J[ 0]*J[ 6]*J[11]*J[13] +
               J[3]*J[ 5]*J[ 8]*J[14] - J[ 1]*J[ 7]*J[ 8]*J[14] -
               J[3]*J[ 4]*J[ 9]*J[14] + J[ 0]*J[ 7]*J[ 9]*J[14] +
               J[1]*J[ 4]*J[11]*J[14] - J[ 0]*J[ 5]*J[11]*J[14] -
               J[2]*J[ 5]*J[ 8]*J[15] + J[ 1]*J[ 6]*J[ 8]*J[15] +
               J[2]*J[ 4]*J[ 9]*J[15] - J[ 0]*J[ 6]*J[ 9]*J[15] -
               J[1]*J[ 4]*J[10]*J[15] + J[ 0]*J[ 5]*J[10]*J[15];
     det     = 1.0 / det;
#ifdef DEBUGG
    printf(" Jacobian   \n");
    printf("   %lf %lf %lf %lf \n", J[ 0], J[ 1], J[ 2], J[ 3]);
    printf("   %lf %lf %lf %lf \n", J[ 4], J[ 5], J[ 6], J[ 7]);
    printf("   %lf %lf %lf %lf \n", J[ 8], J[ 9], J[10], J[11]);
    printf("   %lf %lf %lf %lf \n", J[12], J[13], J[14], J[15]);
    printf("\n INVERSE JACOBIAN MATRIX ORI DET %1.12e\n", det);
    printf(" %lf %lf %lf %lf\n", JI[ 0]*det, JI[ 1]*det, JI[ 2]*det, JI[ 3]*det);
    printf(" %lf %lf %lf %lf\n", JI[ 4]*det, JI[ 5]*det, JI[ 6]*det, JI[ 7]*det);
    printf(" %lf %lf %lf %lf\n", JI[ 8]*det, JI[ 9]*det, JI[10]*det, JI[11]*det);
    printf(" %lf %lf %lf %lf\n", JI[12]*det, JI[13]*det, JI[14]*det, JI[15]*det);
    printf(" ----------------------\n");
    printf(" L %lf %lf %lf %lf \n", L[0], L[1], L[2], L[3]);
#endif
    delta[0] = det*(JI[0 ]*L[0] + JI[1 ]*L[1] + JI[2 ]*L[2] + JI[ 3]*L[3]);
    delta[1] = det*(JI[4 ]*L[0] + JI[5 ]*L[1] + JI[6 ]*L[2] + JI[ 7]*L[3]);
    delta[2] = det*(JI[8 ]*L[0] + JI[9 ]*L[1] + JI[10]*L[2] + JI[11]*L[3]);
    delta[3] = det*(JI[12]*L[0] + JI[13]*L[1] + JI[14]*L[2] + JI[15]*L[3]);
    uvIT[0] -= delta[0];
    uvIT[1] -= delta[1];
    uvIT[2] -= delta[2];
    uvIT[3] -= delta[3];
    if (uvIT[0] < range[0] || uvIT[0] > range[1] ||
        uvIT[1] < range[2] || uvIT[1] > range[3]) {
      pIT[0]  = 0.25*(p0[0] + p1[0] + p2[0] + p3[0]);
      pIT[1]  = 0.25*(p0[1] + p1[1] + p2[1] + p3[1]);
      pIT[2]  = 0.25*(p0[2] + p1[2] + p2[2] + p3[2]);
      uvIT[0] = uvOUT[0]; uvIT[1] = uvOUT[1];
      stat    = EG_invEvaluateGuess(face, pIT, uvIT, pOUT);
      if (stat != EGADS_SUCCESS ||
          uvIT[0] < range[0] || uvIT[0] > range[1] ||
          uvIT[1] < range[2] || uvIT[1] > range[3]) {
        stat      = EG_invEvaluate(face, pIT, uvIT, pOUT);
        if (stat != EGADS_SUCCESS ||
            uvIT[0] < range[0] || uvIT[0] > range[1] ||
            uvIT[1] < range[2] || uvIT[1] > range[3]) return;
        uvOUT[0] = uvIT[0];
        uvOUT[1] = uvIT[1];
        /* we are hopelessly out of range... */
        printf(" EG_minArc4: Out of UVbox!\n");
        return;
      }
    }
    x2 = sqrt(DOT4(delta, delta));
#ifdef DEBUG
    if      (it == 0) x0 = x2;
    else if (it == 1) x1 = x2;
    else {
      e1 = fabs(x1 / x0);
      e2 = fabs(x2 / x1);
      x0 = x1;
      x1 = x2;
    }
#endif
    i = EG_evaluate(face, uvIT, pIT);
#ifdef DEBUGG
    printf(" NEW POINT %lf %lf %lf\n", pIT[0], pIT[1], pIT[2]);
    printf(" du  %lf %lf %lf  dv %lf %lf %lf\n",
           pIT[3], pIT[4], pIT[5], pIT[6], pIT[7], pIT[8]);
    printf(" duu %lf %lf %lf  duv %lf %lf %lf  dvv %lf %lf %lf\n",
           pIT[ 9], pIT[10], pIT[11], pIT[12], pIT[13], pIT[14],
           pIT[15], pIT[16], pIT[17]);
#endif
    if (i != EGADS_SUCCESS) printf(" EG_minArc4: EG_evaluate = %d!\n", i);
    if (i != EGADS_SUCCESS || x2 < EPS10) break;
  }
  if (it == nT)
    printf(" EG_minArc4: Not Converged!\n");

#ifdef DEBUG
  printf("\n\n --------------- REPORT MIN ARC 4 --------------------------- \n");
  if (i != EGADS_SUCCESS) printf(" EG_minArc4: EG_evaluate %d !!\n", i);
  printf(" IT %d DELTA SIZE %1.2e < %1.2e n", it, x2, EPS10);
  if (e1 > EPS10 && e2 > EPS10)
      printf(" CONVERGENCE RATE %lf \n", log(e2) / log(e1));
  r0[0] = pIT[0] - p0[0]; r0[1] = pIT[1] - p0[1]; r0[2] = pIT[2] - p0[2];
  r1[0] = pIT[0] - p1[0]; r1[1] = pIT[1] - p1[1]; r1[2] = pIT[2] - p1[2];
  r2[0] = pIT[0] - p2[0]; r2[1] = pIT[1] - p2[1]; r2[2] = pIT[2] - p2[2];
  r3[0] = pIT[0] - p3[0]; r3[1] = pIT[1] - p3[1]; r3[2] = pIT[2] - p3[2];
  l [0] = DOT(r0, r0); l[1] = DOT(r1, r1);
  l [2] = DOT(r2, r2); l[3] = DOT(r3, r3);
  l [0] = sqrt(l[0]); l[1]  = sqrt(l[1]);
  l [2] = sqrt(l[2]); l[3]  = sqrt(l[3]);
   lt02 = l[0] + l[2];
   lt13 = l[1] + l[3];
  printf(" FINAL l02 = %lf l0 = %lf l2 = %lf ratios %lf | %lf (EXP %lf  %lf)\n",
         lt02, l[0], l[2], l[0] / (lt02*lt02), l[2]/(lt02*lt02),
         fact1*fact1, (1.0 - fact1)*(1.0 - fact1));
  printf(" FINAL l13 = %lf l1 = %lf l3 = %lf ratios %lf | %lf (EXP %lf  %lf)\n",
         lt13, l[1], l[3], l[1] / (lt13*lt13), l[3]/(lt13*lt13),
         fact2*fact2, (1.0 - fact2)*(1.0 - fact2) );
  if (fabs(l[0] / lt02 - fact1       )  > 1.0e-08 ||
      fabs(l[2] / lt02 - (1.0 - fact1)) > 1.0e-08 ||
      fabs(l[1] / lt13 - fact2       )  > 1.0e-08 ||
      fabs(l[3] / lt13 - (1.0 - fact2)) > 1.0e-08 ) {
    printf(" EG_minArc4: Did NOT work!!!!!\n");
    return;
  }
#endif

 uvOUT[0] = uvIT[0];
 uvOUT[1] = uvIT[1];
 return;
}


__HOST_AND_DEVICE__ void
EG_getEdgepoint(const ego edge, double w, double tm, double tp, double *tOUT)
{
  int    nT = 50, stat, i, it, nS = 20;
  double pIT[18], pm[18], pp[18], dt = 0.0, vt0[3], vt1[3];
  double r0[3] = {0.,0.,0.}, r1[3] = {0.,0.,0.};
  double k, xyz[3], res = 1.0, t, tIT, l0, l1, f, ff;
  double s, nEPS = 1.e-10;
#ifdef DEBUG
  double lt, x0 = 0.0, x1 = 0.0, x2 = 0.0, e1 = -1.0, e2 = -1.0;
#endif
  
  k     = (w * w) / ((1.0 - w) * (1.0 - w));
  tIT   = tp * w + tm * (1.0 - w);
  *tOUT = tIT;
  stat  = EG_evaluate(edge, &tm,  pm);
  stat += EG_evaluate(edge, &tp,  pp);
  stat += EG_evaluate(edge, &tIT, pIT);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getEdgePoint before Newton: EG_evaluate error sum = %d\n ",
           stat);
    return;
  }
  
#ifdef DEBUG
  r0[0] = pIT[0] - pm[0]; r1[0] = pIT[0] - pp[0];
  r0[1] = pIT[1] - pm[1]; r1[1] = pIT[1] - pp[1];
  r0[2] = pIT[2] - pm[2]; r1[2] = pIT[2] - pp[2];
  l0    = sqrt(DOT(r0, r0));
  l1    = sqrt(DOT(r1, r1));
  lt    = l0 + l1;
  printf("\n Initial l0 = %lf l1 = %lf split edge in fracs %lf %lf (expected %lf %lf)\n",
         l0, l1, l0 / lt, l1 / lt, w, 1.0 - w);
#endif
  
  for (it = 0 ; it < nT; it++) {
    /* perform a line search such that residual decreases from the previous iteration */
    s = 1.0; f = 0.0;
    for (i = 0; i < nS; i++) {
      t    = tIT + s * dt;
      stat = EG_evaluate(edge, &t, pIT);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_getEdgepoint: EG_evaluate = %d\n", stat);
        s /= 2.0;
        continue;
      }
      if (t < tm || t > tp ) {
        if (it > 0) {
          s /= 2.0;
          continue;
        }
        else if (it == 0) {
#ifdef DEBUG
          printf(" t %lf OUT OF RANGE: [%lf %lf]\n", t, tm, tp);
#endif
          if (i == 0) {
            pIT[0]    = w * pp[0] + (1.0 - w) * pm[0];
            pIT[1]    = w * pp[1] + (1.0 - w) * pm[1];
            pIT[2]    = w * pp[2] + (1.0 - w) * pm[2];
            stat      = EG_invEvaluateGuess(edge, pIT, &tIT, xyz);
            if (stat != EGADS_SUCCESS ||
                tIT < tm || tIT > tp ) {
              printf(" EG_getEdgepoint: EG_invEvaluateGuess %d out-of-range!\n",
                     stat);
              /* we are hopelessly out of range... */
              return;
            }
            continue; /* try again */
          } else {
            /* already tried better guess with inverse evaluate... */
            printf(" EG_getEdgepoint: Initial out-of-range!\n");
            return;
          }
        }
      }
      // distance vector
      r0 [0] = pIT[0] - pm[0]; r1[0]  = pIT[0] - pp[0];
      r0 [1] = pIT[1] - pm[1]; r1[1]  = pIT[1] - pp[1];
      r0 [2] = pIT[2] - pm[2]; r1[2]  = pIT[2] - pp[2];
      l0     = DOT(r0, r0);
      l1     = DOT(r1, r1);
      f      = 0.5 * (l0 - k * l1);
      if (it == 0 || fabs(f) < res) {
        res = fabs(f); /* save off the last residual and exit */
        break;
      }  else {
        s /= 2.0; /* try again if the residual grew */
        continue;
      }
    }
    tIT = t;
    if (i == nS) printf(" EG_getEdgepoint: LineSearch Failed!\n");
    if (i == nS || res < nEPS) break;
    // l0  = DOT(p - pm, p - pm) l1 = DOT(p - pp, p - pp)
    vt0[0] = pIT[3]; vt1[0] = pIT[3];
    vt0[1] = pIT[4]; vt1[1] = pIT[4];
    vt0[2] = pIT[5]; vt1[2] = pIT[5];
    ff     = DOT(vt0, r0) - k * DOT(vt1, r1);
    // dl0 = DOT(dt, p - pm) dl1 = DOT(dt - pp, p - pp)
    dt     = - (f / ff);
#ifdef DEBUG
    x2     = fabs(dt);
    if      (it == 0) x0 = x2;
    else if (it == 1) x1 = x2;
    else {
      e1 = fabs(x1 / x0);
      e2 = fabs(x2 / x1);
      x0 = x1;
      x1 = x2;
    }
#endif
    if (fabs(dt) < nEPS) break;
  }
  
#ifdef DEBUG
  if (e1 >= 0.0 && e2 >= 0.0)
    printf("Newton iterations %d convergence %lf (e1 %lf e2 %lf)\n",
           it, log(e2) / log(e1), e1, e2);
  else printf("Newton iterations %d convergence NA\n", it);
  r0[0] = pIT[0] - pm[0]; r0[1] = pIT[1] - pm[1]; r0[2] = pIT[2] - pm[2];
  r1[0] = pIT[0] - pp[0]; r1[1] = pIT[1] - pp[1]; r1[2] = pIT[2] - pp[2];
  l0    = sqrt(DOT(r0, r0));
  l1    = sqrt(DOT(r1, r1));
  lt    = l0 + l1;
  printf(" FINAL l0 = %lf l1 = %lf split edge in fracs %lf %lf (expected %lf %lf)\n",
         l0, l1, l0 / lt, l1 / lt, w, 1.0 - w);
#endif
  
  if (res >= nEPS && fabs(dt) >= nEPS)
    printf(" EG_getEdgepoint: not converged -- residual %1.2le delta %1.2le (%1.2le)!\n",
           res, fabs(dt), nEPS);
  
  /* found a solution out of range -- report it (for now) */
  if (tIT < tm || tIT > tp )
    printf(" EG_getEdgepoint: solution out-of-range!\n");
  
  *tOUT = tIT;
}


/* Find the point P = S(u,v) midpoint of uv0 and uv1 that minimizes the distance
   between min 1/2 ( l_0 ^ 2 + l_1 ^ 2) + lambda * (l_1 - l_0) = L(u, v, lambda)
   grad(L) = (L_1, L_2, L_3) = (0, 0, 0) -->
      solution using Newton x_n+1 = x_n + delta_n (3x3 system)

   A new point * is inserted along the edge formed be point M-P.

   If L or R are present, dot product between the triangle normals

     L-M-x
     M-R-x
     R-P-x
     P-L-x

      P
    / | \
   /  |  \
  L   x   R
   \  |  /
    \ | /
      M

   is guaranteed to be positive.

*/
__HOST_AND_DEVICE__ void
EG_getSidepoint(const ego face, double fact, const double *uvm,
                const double *uvp, /*@null@*/ const double *uvl,
                /*@null@*/ const double *uvr, double *uvOUT)
{
  int    stat, i, it, nT = 50, nS = 20, fold = 0, foldL, foldR;
  double b, s, dlu0 = 0.0, dlu1 = 0.0, dlv0 = 0.0, dlv1 = 0.0, ddl0, ddl1;
  double detJ, rsdnrm = 1.0, x3 = 0.0, delnrm = 1.0, ctt;
  double pM[18], pP[18], pL[18], pR[18], pIT[18], J[3][3], ATJ[3][3];
  double range[4], r0[3] = {0.,0.,0.}, r1[3] = {0.,0.,0.}, l[2] = {0.,0.};
  double delta[3] = {0.,0.,0.}, L[3] = {0.,0.,0.}, uv[3], uvIT[3], xyz[3];
  double MR[3], RP[3], PL[3], LM[3], RO[3], MO[3], LO[3], PO[3];
  double crMR[3], crRP[3], crPL[3], crLM[3], crITL[3], crITR[3];
#ifdef FACE_NORMAL_FOLD_CHECK
  int    mtype;
#endif
#ifdef DEBUG
  double e1 = 0.0, e2 = 0.0, x0, x1;
#endif

  stat = EG_getRange(face, range, &it);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getSidepoint: EG_egetRange = %d!\n ", stat);
    return;
  }
  /* Initial guess uv  */
  uvIT[0]  = fact*uvp[0]  + (1.0-fact)*uvm[0];
  uvIT[1]  = fact*uvp[1]  + (1.0-fact)*uvm[1];
  uvIT[2]  = 0.0;
  uvOUT[0] = uvIT[0];
  uvOUT[1] = uvIT[1];
  stat     = EG_evaluate(face, uvm, pM);
  stat    += EG_evaluate(face, uvp, pP);
  if (uvl != NULL) stat += EG_evaluate(face, uvl, pL);
  if (uvr != NULL) stat += EG_evaluate(face, uvr, pR);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getSidepoint: EG_evaluate = %d!\n ", stat);
    return;
  }
#ifdef DEBUG
  EG_evaluate(face, uvIT, pIT);
  r0[0] = pIT[0] - pM[0];
  r0[1] = pIT[1] - pM[1];
  r0[2] = pIT[2] - pM[2];

  r1[0] = pIT[0] - pP[0];
  r1[1] = pIT[1] - pP[1];
  r1[2] = pIT[2] - pP[2];

  l[0]  = DOT(r0, r0);
  l[1]  = DOT(r1, r1);
  b     = sqrt(l[0]) + sqrt(l[1]);
  printf(" \n\n --------------------------------------------------- \n");
  printf(" INITIAL GUESS: TOTAL ARC %lf l0 = ||pIT - pM||^2 = %lf l1 = ||pIT - pP||^2 = %lf\n",
         b, l[0], l[1]);
  printf(" ratios %lf %lf (EXPECTED %lf %lf)\n",
         l[0]/(b*b), l[1]/(b*b), fact*fact, (1.0-fact)*(1.0-fact));
  l[0]  = sqrt(l[0]);
  l[1]  = sqrt(l[1]);
  printf(" As distances sqrt(l0) = ||pIT - pM|| = %lf sqrt (l1) = ||pIT - pP|| = %lf\n",
         l[0], l[1]);
  printf(" ratios %lf %lf (EXPECTED %lf %lf)\n", l[0]/b, l[1]/b, fact, (1.0-fact));
  printf(" --------------------------------------------------- \n");
#endif

  /* segments connected to L */
  if (uvl != NULL) {
    PL[0] = pL[0] - pP[0]; PL[1] = pL[1] - pP[1]; PL[2] = pL[2] - pP[2];
    LM[0] = pM[0] - pL[0]; LM[1] = pM[1] - pL[1]; LM[2] = pM[2] - pL[2];
  }
  else
    foldL = 0;

  /* segments connected to R */
  if (uvr != NULL) {
    MR[0] = pR[0] - pM[0]; MR[1] = pR[1] - pM[1]; MR[2] = pR[2] - pM[2];
    RP[0] = pP[0] - pR[0]; RP[1] = pP[1] - pR[1]; RP[2] = pP[2] - pR[2];
  }
  else
    foldR = 0;


#ifdef FACE_NORMAL_FOLD_CHECK
  /* get mtype=SFORWARD or mtype=SREVERSE for the face to get topology normal */
  mtype = face->mtype;
#else
  /* use the original triangle normal as the reference */
  if (uvl != NULL) {
    CRXSS(PL, LM, crITL);
  }
  if (uvr != NULL) {
    CRXSS(MR, RP, crITR);
  }
#endif

  /* allow the initial solution to have folds. Hopefully they get fixed. */
  fold = 1;

  ctt     = (1.0 - fact) * (1.0 - fact) / (fact * fact);
  for (it = 0; it < nT; it++) {

    /* perform a line search such that residual decreases from the previous iteration */
    s = 1.0;
    for (i = 0; i < nS; i++) {
      uv[0] = uvIT[0] + s*delta[0];
      uv[1] = uvIT[1] + s*delta[1];
      uv[2] = uvIT[2] + s*delta[2];

      stat = EG_evaluate(face, uv, pIT);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_getSidepoint: EG_evaluate = %d\n", stat);
        s /= 2.0;
        continue;
      }

#ifdef FACE_NORMAL_FOLD_CHECK
      if (uvl != NULL || uvr != NULL) {
        /* get the normal vector at the proposed point */
        CRXSS(pIT+3, pIT+6, crITL);
        crITR[0] = (crITL[0] *= mtype);
        crITR[1] = (crITL[1] *= mtype);
        crITR[2] = (crITL[2] *= mtype);
      }
#endif

      if (uvl != NULL) {
        /* Get vectors */
        LO[0] = pIT[0] - pL[0]; LO[1] = pIT[1] - pL[1]; LO[2] = pIT[2] - pL[2];
        MO[0] = pIT[0] - pM[0]; MO[1] = pIT[1] - pM[1]; MO[2] = pIT[2] - pM[2];

        /* get areas */
        CRXSS(PL, LO, crPL);
        CRXSS(LM, MO, crLM);

        /* check for fold */
        foldL = (DOT(crPL, crITL) < 0.0) || (DOT(crLM, crITL) < 0.0);
      }

      if (uvr != NULL) {
        /* Get vectors */
        RO[0] = pIT[0] - pR[0]; RO[1] = pIT[1] - pR[1]; RO[2] = pIT[2] - pR[2];
        PO[0] = pIT[0] - pP[0]; PO[1] = pIT[1] - pP[1]; PO[2] = pIT[2] - pP[2];

        /* get areas */
        CRXSS(MR, RO, crMR);
        CRXSS(RP, PO, crRP);

        /* check for fold */
        foldR = (DOT(crRP, crITR) < 0.0) || (DOT(crMR, crITR) < 0.0);
      }

      /* Found a solution without a fold! Don't allow a new fold. */
      if (!(foldL || foldR)) fold = 0;

      if ((it > 0) &&
          (uv[0] < range[0] || uv[0] > range[1] ||
           uv[1] < range[2] || uv[1] > range[3] ||
          (fold == 0 && (foldL || foldR)) )) {
        s /= 2.0;
        continue;
      } else if (it == 0 && fold == 1) {
#ifdef DEBUG
        printf(" UV %lf %lf OUT OF RANGE: %lf %lf %lf %lf\n", uv[0], uv[1],
               range[0], range[1], range[2], range[3]);
#endif
        if (uv[0] < range[0] || uv[0] > range[1] ||
            uv[1] < range[2] || uv[1] > range[3] ||
            foldL || foldR ) {

          if (fold == 2) {
            /* already tried better guess with inverse evaluate... */
            printf(" EG_getSidepoint: Initial out-of-range!\n");
            return;
          }

          /* try a better initial guess with an inverse evaluate */
          pIT[0] = fact*pP[0] + (1.0-fact)*pM[0];
          pIT[1] = fact*pP[1] + (1.0-fact)*pM[1];
          pIT[2] = fact*pP[2] + (1.0-fact)*pM[2];
          stat   = EG_invEvaluate(face, pIT, uvOUT, xyz);
          if (stat != EGADS_SUCCESS ||
              uvOUT[0] < range[0] || uvOUT[0] > range[1] ||
              uvOUT[1] < range[2] || uvOUT[1] > range[3]) {
            printf(" EG_getSidepoint: EG_invEvaluate Out-of-Range %d!\n", stat);
            /* we are hopelessly out of range... */
            return;
          }
          fold = 2;
          continue; /* try again */
        }
      }

      r0[0] = pIT[0] - pM[0];
      r0[1] = pIT[1] - pM[1];
      r0[2] = pIT[2] - pM[2];

      r1[0] = pIT[0] - pP[0];
      r1[1] = pIT[1] - pP[1];
      r1[2] = pIT[2] - pP[2];

      l[0]  = DOT(r0, r0);
      l[1]  = DOT(r1, r1);
      dlu0  = 2.0 * (r0[0] * pIT[3] + r0[1] * pIT[4] + r0[2] * pIT[5]);
      dlu1  = 2.0 * (r1[0] * pIT[3] + r1[1] * pIT[4] + r1[2] * pIT[5]);
      dlv0  = 2.0 * (r0[0] * pIT[6] + r0[1] * pIT[7] + r0[2] * pIT[8]);
      dlv1  = 2.0 * (r1[0] * pIT[6] + r1[1] * pIT[7] + r1[2] * pIT[8]);

      /* UPDATE FUNCTION L */
      L[0]  = dlu0 * (l[0] - ctt * uv[2]) + dlu1 * (l[1] + uv[2]);
      L[1]  = dlv0 * (l[0] - ctt * uv[2]) + dlv1 * (l[1] + uv[2]);
      L[2]  = l[1] - ctt * l[0];

      /* check for convergence on the residual */
      rsdnrm  = sqrt(DOT(L,L));

      if (it == 0 || rsdnrm < x3) {
        x3 = rsdnrm; /* save off the last residual and exit */
        break;
      }  else {
        s /= 2.0; /* try again if the residual grew */
        continue;
      }
    }
    if (i == nS) {
      printf(" EG_getSidepoint: LineSearch Failed (%lf %lf)!\n",
             uvIT[0], uvIT[1]);
      pIT[0] = fact*pP[0] + (1.0-fact)*pM[0];
      pIT[1] = fact*pP[1] + (1.0-fact)*pM[1];
      pIT[2] = fact*pP[2] + (1.0-fact)*pM[2];
      stat   = EG_invEvaluate(face, pIT, uvOUT, xyz);
      if (stat != EGADS_SUCCESS ||
          uvOUT[0] < range[0] || uvOUT[0] > range[1] ||
          uvOUT[1] < range[2] || uvOUT[1] > range[3]) {
        printf(" EG_getSidepoint: EG_invEvaluate %d out-of-range!\n", stat);
      }
      return;
    }

#ifdef DEBUG
    if      (it == 0) x0 = rsdnrm;
    else if (it == 1) x1 = rsdnrm;
    else {
      e1 = fabs(x1 / x0);
      e2 = fabs(rsdnrm / x1);
      x0 = x1;
      x1 = rsdnrm;
    }
#endif

    /* update the solution */
    uvIT[0] = uv[0];
    uvIT[1] = uv[1];
    uvIT[2] = uv[2];

    if (i      ==   nS) break;  /* line search failed */
    if (rsdnrm < EPS10) break;  /* converged! */

    /* duu */
    ddl0    = 2.0 * (pIT[9 ] * r0[0] + pIT[3] * pIT[3] +
                     pIT[10] * r0[1] + pIT[4] * pIT[4] +
                     pIT[11] * r0[2] + pIT[5] * pIT[5]);
    ddl1    = 2.0 * (pIT[9 ] * r1[0] + pIT[3] * pIT[3] +
                     pIT[10] * r1[1] + pIT[4] * pIT[4] +
                     pIT[11] * r1[2] + pIT[5] * pIT[5]);
    b       = ddl0 * (l[0] - ctt * uvIT[2]) + ddl1* (l[1] + uvIT[2]);
    J[0][0] = dlu0 * dlu0 + dlu1 * dlu1 + b;

    /* duv */
    ddl0    = 2.0 * (pIT[12] * r0[0] + pIT[3] * pIT[6] +
                     pIT[13] * r0[1] + pIT[4] * pIT[7] +
                     pIT[14] * r0[2] + pIT[5] * pIT[8]);

    ddl1    = 2.0 * (pIT[12] * r1[0] + pIT[3] * pIT[6] +
                     pIT[13] * r1[1] + pIT[4] * pIT[7] +
                     pIT[14] * r1[2] + pIT[5] * pIT[8]);
    b       = ddl0 * (l[0] - ctt * uvIT[2]) + ddl1 * (l[1] + uvIT[2]);
    J[0][1] = dlu0 * dlv0 + dlu1 * dlv1 + b;

    /* dvv */
    ddl0    = 2.0 * (pIT[15] * r0[0] + pIT[6] * pIT[6] +
                     pIT[16] * r0[1] + pIT[7] * pIT[7] +
                     pIT[17] * r0[2] + pIT[8] * pIT[8]);

    ddl1    = 2.0 * (pIT[15] * r1[0] + pIT[6] * pIT[6] +
                     pIT[16] * r1[1] + pIT[7] * pIT[7] +
                     pIT[17] * r1[2] + pIT[8] * pIT[8]);
    b       = ddl0 * (l[0] - ctt * uvIT[2]) + ddl1 * (l[1] + uvIT[2]);
    J[1][1] = dlv0 * dlv0 + dlv1 * dlv1 + b;

    J[0][2] = dlu1 - ctt * dlu0;
    J[1][2] = dlv1 - ctt * dlv0;
    J[1][0] = J[0][1];
    J[2][0] = J[0][2];
    J[2][1] = J[1][2];
    J[2][2] = 0.0;

    /* Solve Linear System: J * delta = - L
       For now: Invert Jacobian directly --> delta = J^-1 L  */
    detJ    = J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] +
    J[0][2] * J[2][1] * J[1][0] - J[0][2] * J[1][1] * J[2][0] -
    J[1][2] * J[2][1] * J[0][0] - J[2][2] * J[1][0] * J[0][1];
    if (detJ == 0.0) break;

    ATJ[0][0] =   J[1][1] * J[2][2] - J[2][1] * J[1][2];
    ATJ[0][1] = -(J[0][1] * J[2][2] - J[2][1] * J[0][2]);
    ATJ[0][2] =   J[0][1] * J[1][2] - J[1][1] * J[0][2];
    ATJ[1][0] = -(J[1][0] * J[2][2] - J[2][0] * J[1][2]);
    ATJ[1][1] =   J[0][0] * J[2][2] - J[2][0] * J[0][2];
    ATJ[1][2] = -(J[0][0] * J[1][2] - J[1][0] * J[0][2]);
    ATJ[2][0] =   J[1][0] * J[2][1] - J[2][0] * J[1][1];
    ATJ[2][1] = -(J[0][0] * J[2][1] - J[2][0] * J[0][1]);
    ATJ[2][2] =   J[0][0] * J[1][1] - J[1][0] * J[0][1];

    detJ      = 1.0 / detJ;
    delta[0]  = -detJ * (ATJ[0][0] * L[0] + ATJ[0][1] * L[1] + ATJ[0][2] * L[2]);
    delta[1]  = -detJ * (ATJ[1][0] * L[0] + ATJ[1][1] * L[1] + ATJ[1][2] * L[2]);
    delta[2]  = -detJ * (ATJ[2][0] * L[0] + ATJ[2][1] * L[1] + ATJ[2][2] * L[2]);

    /* check for convergence on the parameter update */
    delnrm    = sqrt(DOT(delta, delta));
    if (delnrm < EPS10) break;  /* converged! */
  }

#ifdef DEBUG
  if (i != EGADS_SUCCESS) printf("EG_evaluate %d !!\n", i);
  printf("IT %d DELTA SIZE %1.2e < %1.2e L [%lf  %lf  %lf]\n",
         it, rsdnrm, EPS10, L[0], L[1], L[2]);
  if (e1 > EPS10 && e2 > EPS10) printf("CONVERGENCE RATE %lf \n", log(e2)/log(e1));
  b    = sqrt(l[0]) + sqrt(l[1]);
  printf(" --------------------------------------------------- \n");
  printf(" NEW TOTAL ARC %lf with l0 = ||pIT - pM||^2 = %lf l1 = ||PIT - pP||^2 = %lf\n",
         b, l[0], l[1]);
  printf(" ratios %lf %lf (EXPECTED %lf %lf)\n",
         l[0]/(b*b), l[1]/(b*b), fact*fact, (1.0-fact)*(1.0-fact));
  l[0] = sqrt(l[0]);
  l[1] = sqrt(l[1]);
  printf(" As Actual chords ||pIT - pM|| = %lf ||pIT - pP|| = %lf\n", l[0], l[1]);
  printf(" ratios %lf %lf (EXPECTED %lf %lf)\n",
         l[0]/b, l[1]/b, fact, (1.0 - fact));
  printf(" --------------------------------------------------- \n");
#endif

  if (rsdnrm >= EPS10 && delnrm >= EPS10)
    printf(" EG_getSidepoint: not converged -- residual %1.2le delta %1.2le (%1.2le)!\n",
           rsdnrm, delnrm, EPS10);
  
  /* found a solution out of range -- report it (for now) */
  if (uvIT[0] < range[0] || uvIT[0] > range[1] ||
      uvIT[1] < range[2] || uvIT[1] > range[3])
    printf(" EG_getSidepoint: solution out-of-range!\n");

  if (foldL || foldR) {
    printf(" EG_getSidepoint: solution with fold!\n");
#ifdef DEBUGM
    printf("VARIABLES=X, Y, Z, U, V\n");
    printf("ZONE N=5, E=4, ET=TRIANGLE F=FEPOINT\n");
    printf("%f, %f, %f, %f, %f\n", pM[0], pM[1], pM[2], uvm[0], uvm[1]);
    printf("%f, %f, %f, %f, %f\n", pR[0], pR[1], pR[2], uvr[0], uvr[1]);
    printf("%f, %f, %f, %f, %f\n", pP[0], pP[1], pP[2], uvp[0], uvp[1]);
    printf("%f, %f, %f, %f, %f\n", pL[0], pL[1], pL[2], uvl[0], uvl[1]);
    printf("%f, %f, %f, %f, %f\n", pIT[0], pIT[1], pIT[2], uvIT[0], uvIT[1]);
    printf("1 2 5\n");
    printf("2 3 5\n");
    printf("3 4 5\n");
    printf("4 1 5\n");
    printf("ZONE N=4, E=2, ET=TRIANGLE F=FEPOINT\n");
    printf("%f, %f, %f, %f, %f\n", pM[0], pM[1], pM[2], uvm[0], uvm[1]);
    printf("%f, %f, %f, %f, %f\n", pR[0], pR[1], pR[2], uvr[0], uvr[1]);
    printf("%f, %f, %f, %f, %f\n", pP[0], pP[1], pP[2], uvp[0], uvp[1]);
    printf("%f, %f, %f, %f, %f\n", pL[0], pL[1], pL[2], uvl[0], uvl[1]);
    printf("1 2 3\n");
    printf("3 4 1\n");
#endif
    pIT[0] = fact*pP[0] + (1.0-fact)*pM[0];
    pIT[1] = fact*pP[1] + (1.0-fact)*pM[1];
    pIT[2] = fact*pP[2] + (1.0-fact)*pM[2];
    stat   = EG_invEvaluate(face, pIT, uvIT, xyz);
    if (stat != EGADS_SUCCESS ||
        uvIT[0] < range[0] || uvIT[0] > range[1] ||
        uvIT[1] < range[2] || uvIT[1] > range[3]) {
      printf(" EG_getSidepoint: EG_invEvaluate out-of-range %d!\n", stat);
    }
  }

  uvOUT[0] = uvIT[0];
  uvOUT[1] = uvIT[1];
}


__HOST_AND_DEVICE__ static void
EG_quadThread(void *struc)
{
  int     index, stat;
  long    ID;
  EMPquad *qthread;

  qthread = (EMPquad *) struc;

  /* get our identifier */
  ID = EMP_ThreadID();

  /* look for work */
  for (;;) {

    /* only one thread at a time here -- controlled by a mutex! */
    if (qthread->mutex != NULL) EMP_LockSet(qthread->mutex);
    for (index = qthread->index; index < qthread->end; index++) {
      if (qthread->ntess->tess2d[index].tfi ==    1) continue;
      if (qthread->bodydata->qm[index]      == NULL) continue;
      if (qthread->bodydata->qm[index]->fID ==    0) continue;
      break;
    }
    qthread->index = index+1;
    if (qthread->mutex != NULL) EMP_LockRelease(qthread->mutex);
    if (index >= qthread->end) break;

    /* do the work */
    stat = EG_meshRegularization(qthread->bodydata->qm[index]);
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Warning: EG_fullMeshRegularization %d = %d (EG_quadTess)!\n",
             index+1, stat);
  }

  /* exhausted all work -- exit */
  if (ID != qthread->master) EMP_ThreadExit();
}


__HOST_AND_DEVICE__ int
EG_quadTess(const egObject *tess, egObject **quadTess)
{
  int          i, j, k, m, n, nedges, nfaces, stat, npts, ntris, nside, alen;
  int          outLevel, iv, np, nt, is, ie, ien, oclass, mtype, atype;
  int          sum[2], side[4], degens[2], iuv[2], *table, *tris, *senses;
  int          i0, i1, i2, i3, otri, flip;
  const int    *ptype, *pindex, *trs, *trc, *ints;
  double       result[18], xyz[3], uv[2], uvm[2], uvp[2], trange[2], t;
  double       *coords, *parms;
  long         start;
  const double *xyzs, *ts, *uvs, *reals;
  const char   *str;
  egTessel     *btess, *ntess;
  egObject     *obj, *geom, *newTess, **nodes, **edges, **faces, **objs;
  midside      *mid;
  bodyQuad     bodydata;
  EMPquad      qthread;
  void         **threads = NULL;
#ifdef TRIOUT
  FILE         *fp;
  char         filename[100];
#endif
  static int   sides[3][2] = {{1,2}, {2,0}, {0,1}       };
  static int   sideq[4][2] = {{1,2}, {2,5}, {5,0}, {0,1}};
  static int   neigq[4]    = { 0,     3,     4,     2   };

  *quadTess = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (EG_sameThread(tess))          return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(tess);

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_quadTess)!\n");
    return EGADS_NOTFOUND;
  }
  if (btess->done == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Tessellation is open (EG_quadTess)!\n");
    return EGADS_TESSTATE;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_quadTess)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_quadTess)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_quadTess)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_quadTess)!\n");
    return EGADS_NODATA;
  }
  if ((btess->tess2d == NULL) && (btess->nFace != 0)) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_quadTess)!\n");
    return EGADS_NODATA;
  }

  /* initialize the new tessellation object */
  stat = EG_initTessBody(obj, &newTess);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_initTessBody = %d (EG_quadTess)!\n", stat);
    return stat;
  }

  stat = EG_getBodyTopos(obj, NULL, EDGE, &nedges, &edges);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getBodyTopos E = %d (EG_quadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }

  /* rebuild the Edges */

  for (j = i = 0; i < nedges; i++) {
    if (edges[i]->mtype == DEGENERATE) continue;
    stat = EG_getTessEdge(tess, i+1, &npts, &xyzs, &ts);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTessEdge %d = %d (EG_quadTess)!\n",
               i+1, stat);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
    if (npts == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTessEdge %d -- no points (EG_quadTess)!\n",
               i+1);
      EG_free(edges);
      EG_deleteObject(newTess);
      return EGADS_INDEXERR;
    }
    if (npts > j) j = npts;
  }
  /* allocate to the maximum length */
  coords = (double *) EG_alloc(4*(2*j-1)*sizeof(double));
  if (coords == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d points (EG_quadTess)!\n", j);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }
  parms = &coords[3*(2*j-1)];

  for (i = 0; i < nedges; i++) {
    if (edges[i]->mtype == DEGENERATE) continue;
    stat = EG_getTessEdge(tess, i+1, &npts, &xyzs, &ts);
    if (stat != EGADS_SUCCESS) continue;

    for (j = 0; j < npts-1; j++) {
      parms[2*j  ] = ts[j];
/*    parms[2*j+1] = 0.5*(ts[j] + ts[j+1]);  */
      EG_getEdgepoint(edges[i], 0.5, ts[j], ts[j+1], &parms[2*j+1]);
      stat = EG_evaluate(edges[i], &parms[2*j+1], result);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_evaluate Edge %d/%d = %d (EG_quadTess)!\n",
                 i+1, j+1, stat);
        EG_free(coords);
        EG_free(edges);
        EG_deleteObject(newTess);
        return stat;
      }
      coords[6*j  ] = xyzs[3*j  ];
      coords[6*j+1] = xyzs[3*j+1];
      coords[6*j+2] = xyzs[3*j+2];
      coords[6*j+3] = result[0];
      coords[6*j+4] = result[1];
      coords[6*j+5] = result[2];
    }
    j = npts-1;
    parms[2*j  ]  = ts[j];
    coords[6*j  ] = xyzs[3*j  ];
    coords[6*j+1] = xyzs[3*j+1];
    coords[6*j+2] = xyzs[3*j+2];

    stat = EG_setTessEdge(newTess, i+1, 2*npts-1, coords, parms);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessEdge %d = %d (EG_quadTess)!\n",
               i+1, stat);
      EG_free(coords);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
  }
  EG_free(coords);
  ntess = (egTessel *) newTess->blind;

  /* size and allocate temporary arrays */

  stat = EG_getBodyTopos(obj, NULL, FACE, &nfaces, &faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getBodyTopos F = %d (EG_quadTess)!\n", stat);
    EG_free(edges);
    EG_deleteObject(newTess);
    return stat;
  }
  for (nside = npts = ntris = i = 0; i < nfaces; i++) {
    stat = EG_getTessFace(tess, i+1, &np, &xyzs, &uvs, &ptype, &pindex,
                          &nt, &trs, &trc);
    if ((stat != EGADS_SUCCESS) || (nt == 0)) {
      if (outLevel > 0)
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: EG_getTessFace %d = %d (EG_quadTess)!\n",
                 i+1, stat);
        } else {
          printf(" EGADS Error: Face %d has no tessellation (EG_quadTess)!\n",
                 i+1);
        }
      EG_free(faces);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
    for (sum[0] = sum[1] = j = 0; j < nt; j++)
      for (k = 0; k < 3; k++)
        if (trc[3*j+k] > 0) {
          sum[0]++;
        } else {
          sum[1]++;
        }
    k = sum[0]/2 + sum[1];
    if (k  > nside) nside = k;
    if (np > npts)  npts  = np;
    if (nt > ntris) ntris = nt;
  }
  k      = npts + ntris + nside;
  coords = (double *) EG_alloc(5*k*sizeof(double));
  if (coords == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face points (EG_quadTess)!\n", k);
    EG_free(faces);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }
  parms = &coords[3*k];
  tris  = (int *) EG_alloc(6*3*ntris*sizeof(int));
  if (tris == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face tris (EG_quadTess)!\n", 3*ntris);
    EG_free(coords);
    EG_free(faces);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }
  mid = (midside *) EG_alloc(nside*sizeof(midside));
  if (mid == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face mids (EG_quadTess)!\n", nside);
    EG_free(tris);
    EG_free(coords);
    EG_free(faces);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }
  table = (int *) EG_alloc(npts*sizeof(int));
  if (table == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face table (EG_quadTess)!\n", npts);
    EG_free(mid);
    EG_free(tris);
    EG_free(coords);
    EG_free(faces);
    EG_free(edges);
    EG_deleteObject(newTess);
    return EGADS_MALLOC;
  }

  /* fill in the tris, a Face at a time */

  for (i = 0; i < nfaces; i++) {
    if (btess->tess2d[i].tfi == 1) continue;
    stat = EG_getTessFace(tess, i+1, &np, &xyzs, &uvs, &ptype, &pindex,
                          &nt, &trs, &trc);
    if (stat != EGADS_SUCCESS) continue;

    /* find degenerate nodes (if any) */
    degens[0] = degens[1] = 0;
    iuv[0]    = iuv[1]    = 0;
    stat      = EG_getBodyTopos(obj, faces[i], EDGE, &k, &objs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Internal: EG_getBodyTopos on Face %d = %d\n", i+1, stat);
    } else {
      for (j = 0; j < k; j++) {
        stat = EG_getTopology(objs[j], &geom, &oclass, &mtype,
                              trange, &n, &nodes, &senses);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Internal: EG_getTopology on Edge = %d\n", stat);
          continue;
        }
        if (mtype != DEGENERATE) continue;
        stat = EG_getEdgeUVeval(faces[i], objs[j], 0, trange[0], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Internal: EG_getEdgeUVeval = %d\n", stat);
          continue;
        }
        n = EG_indexBodyTopo(obj, nodes[0]);
        if (n > 0) {
          if (degens[0] == 0) {
            degens[0] = n;
            if (result[3] != 0.0) iuv[0] = 1;
          } else if (degens[1] == 0) {
            degens[1] = n;
            if (result[3] != 0.0) iuv[1] = 1;
          } else {
            printf(" EGADS Info: More than 2 Degen Nodes in Face %d!\n", i+1);
          }
        }
      }
      EG_free(objs);
    }

    /* make the vertices */
    for (j = 0; j < np; j++) {
      coords[3*j  ] = xyzs[3*j  ];
      coords[3*j+1] = xyzs[3*j+1];
      coords[3*j+2] = xyzs[3*j+2];
      parms[2*j  ]  = uvs[2*j  ];
      parms[2*j+1]  = uvs[2*j+1];
      table[j]      = NOTFILLED;
    }
    iv = np;

    /* get triangle side midpoint insertions */
    for (is = j = 0; j < nt; j++)
      for (k = 0; k < 3; k++) {
        i1   = trs[3*j+sides[k][0]] - 1;
        i2   = trs[3*j+sides[k][1]] - 1;
        flip = 0;
        if (i2 < i1) {
          stat = i1;
          i1   = i2;
          i2   = stat;
          flip = 1;
        }
        m = EG_findMidSide(i1, i2, table, mid);
        if (m > 0) continue;
        uvm[0] = uvs[2*i1  ];
        uvm[1] = uvs[2*i1+1];
        uvp[0] = uvs[2*i2  ];
        uvp[1] = uvs[2*i2+1];
        if (trc[3*j+k] > 0) {
          /* Interior side */
          otri = trc[3*j+k] - 1;
          i0   = trs[3*j+k] - 1;
          i3   = trs[3*otri] + trs[3*otri+1] + trs[3*otri+2] - i1 - i2 - 3;
          if (flip == 1) {
            stat = i0;
            i0   = i3;
            i3   = stat;
          }
          if ((ptype[i1] == 0) && (pindex[i1] == degens[0])) {
            if (iuv[0] == 0) {
              uvm[0] = uvs[2*i2  ];
            } else {
              uvm[1] = uvs[2*i2+1];
            }
          } else if ((ptype[i1] == 0) && (pindex[i1] == degens[1])) {
            if (iuv[1] == 0) {
              uvm[0] = uvs[2*i2  ];
            } else {
              uvm[1] = uvs[2*i2+1];
            }
          } else if ((ptype[i2] == 0) && (pindex[i2] == degens[0])) {
            if (iuv[0] == 0) {
              uvp[0] = uvs[2*i1  ];
            } else {
              uvp[1] = uvs[2*i1+1];
            }
          } else if ((ptype[i2] == 0) && (pindex[i2] == degens[1])) {
            if (iuv[1] == 0) {
              uvp[0] = uvs[2*i1  ];
            } else {
              uvp[1] = uvs[2*i1+1];
            }
          }
          EG_getSidepoint(faces[i], 0.5, uvm, uvp, &uvs[2*i0], &uvs[2*i3], uv);
          stat = EG_evaluate(faces[i], uv, result);
        } else {
          /* Edge side */
          ie   = -trc[3*j+k] - 1;
          stat = EG_getTopology(edges[ie], &geom, &oclass, &mtype,
                                trange, &n, &objs, &senses);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: EG_getTopology %d/%d = %d (EG_quadTess)!\n",
                     i+1, ie+1, stat);
            EG_free(table);
            EG_free(mid);
            EG_free(tris);
            EG_free(coords);
            EG_free(faces);
            EG_free(edges);
            EG_deleteObject(newTess);
            return stat;
          }
          n = EG_indexBodyTopo(obj, objs[0]);
          if ((ptype[i1] == 0) && (ptype[i2] == 0)) {
            ien = 1;
          } else if (ptype[i1] == 0) {
            if (mtype == ONENODE) {
              if ((xyzs[3*i2  ] == ntess->tess1d[ie].xyz[6]) &&
                  (xyzs[3*i2+1] == ntess->tess1d[ie].xyz[7]) &&
                  (xyzs[3*i2+2] == ntess->tess1d[ie].xyz[8])) {
                ien = 1;
              } else {
                ien = ntess->tess1d[ie].npts - 2;
              }
            } else {
              if (pindex[i1] == n) {
                ien = 1;
              } else {
                ien = 2*ptype[i2] - 1;
              }
            }
            if (ie+1 != pindex[i2])
              printf(" EGADS Info: Edge mismatch %d %d\n", ie+1, pindex[i2]);
          } else if (ptype[i2] == 0) {
            if (mtype == ONENODE) {
              if ((xyzs[3*i1  ] == ntess->tess1d[ie].xyz[6]) &&
                  (xyzs[3*i1+1] == ntess->tess1d[ie].xyz[7]) &&
                  (xyzs[3*i1+2] == ntess->tess1d[ie].xyz[8])) {
                ien = 1;
              } else {
                ien = ntess->tess1d[ie].npts - 2;
              }
            } else {
              if (pindex[i2] == n) {
                ien = 1;
              } else {
                ien = 2*ptype[i1] - 1;
              }
            }
            if (ie+1 != pindex[i1])
              printf(" EGADS Info: Edge mismatch %d %d\n", ie+1, pindex[i1]);
          } else {
            ien = 2*ptype[i1] - 1;
            if (ien > 2*ptype[i2]-1) ien = 2*ptype[i2] - 1;
            if (ien < 0) {
              printf(" EGADS Internal: EG_quadTess ien = %d!\n", ien);
              continue;
            }
            if ((ie+1 != pindex[i1]) || (ie+1 != pindex[i2]))
              printf(" EGADS Info: Edge mismatch %d %d %d\n", ie+1, pindex[i1],
                     pindex[i2]);
          }
          t    = ntess->tess1d[ie].t[ien];
          stat = EG_getEdgeUV(faces[i], edges[ie], 0, t, uv);
          if (stat == EGADS_TOPOERR) {
            /* sense in Face twice! */
            stat = EG_getEdgeUV(faces[i], edges[ie], -1, t, uvm);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Info:  EG_getEdgeUV -1\n");
            } else {
              stat = EG_getEdgeUV(faces[i], edges[ie],  1, t, uvp);
              if (stat != EGADS_SUCCESS) {
                printf(" EGADS Info:  EG_getEdgeUV +1\n");
              } else {
                result[0] = sqrt((uvm[0]-uv[0])*(uvm[0]-uv[0]) +
                                 (uvm[1]-uv[1])*(uvm[1]-uv[1]));
                result[1] = sqrt((uvp[0]-uv[0])*(uvp[0]-uv[0]) +
                                 (uvp[1]-uv[1])*(uvp[1]-uv[1]));
                if (result[0] < result[1]) {
                  uv[0] = uvm[0];
                  uv[1] = uvm[1];
                } else {
                  uv[0] = uvp[0];
                  uv[1] = uvp[1];
                }
              }
            }
          }
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: EG_getEdgeUV %d/%d = %d (EG_quadTess)!\n",
                     i+1, ie+1, stat);
            EG_free(table);
            EG_free(mid);
            EG_free(tris);
            EG_free(coords);
            EG_free(faces);
            EG_free(edges);
            EG_deleteObject(newTess);
            return stat;
          }
          result[0] = ntess->tess1d[ie].xyz[3*ien  ];
          result[1] = ntess->tess1d[ie].xyz[3*ien+1];
          result[2] = ntess->tess1d[ie].xyz[3*ien+2];
        }
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_evaluate %d %d/%d = %d (EG_quadTess)!\n",
                   i+1, j+1, k+1, stat);
          EG_free(table);
          EG_free(mid);
          EG_free(tris);
          EG_free(coords);
          EG_free(faces);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
        }
        if (m == 0) {
          table[i1]      = is;
        } else {
          mid[-m-1].next = is;
        }
        mid[is].vert2 = i2;
        mid[is].nvert = iv + 1;
        mid[is].next  = NOTFILLED;
        is++;
        coords[3*iv  ] = result[0];
        coords[3*iv+1] = result[1];
        coords[3*iv+2] = result[2];
        parms[2*iv  ]  = uv[0];
        parms[2*iv+1]  = uv[1];
        iv++;
      }

    /* fill in the new triangulation */
    for (j = 0; j < nt; j++) {
      for (k = 0; k < 3; k++) {
        i1 = trs[3*j+sides[k][0]] - 1;
        i2 = trs[3*j+sides[k][1]] - 1;
        if (i2 < i1) {
          stat = i1;
          i1   = i2;
          i2   = stat;
        }
        side[k] = EG_findMidSide(i1, i2, table, mid);
        if (side[k] <= 0) {
          if (outLevel > 0)
            printf(" EGADS Error: findMidSide %d = %d (EG_quadTess)!\n",
                   i+1, side[k]);
          EG_free(table);
          EG_free(mid);
          EG_free(tris);
          EG_free(coords);
          EG_free(faces);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
        }
      }
      /* get middle vertex based on mid-sides */
      t    = 1.0/3.0;
      stat = EG_baryInsert(faces[i], t, t, t, &parms[2*(side[0]-1)],
                           &parms[2*(side[1]-1)], &parms[2*(side[2]-1)], uv);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Info: EG_baryInsert = %d (EG_quadTess)!\n", stat);
        xyz[0] = (coords[3*(side[0]-1)  ] + coords[3*(side[1]-1)  ] +
                  coords[3*(side[2]-1)  ])/3.0;
        xyz[1] = (coords[3*(side[0]-1)+1] + coords[3*(side[1]-1)+1] +
                  coords[3*(side[2]-1)+1])/3.0;
        xyz[2] = (coords[3*(side[0]-1)+2] + coords[3*(side[1]-1)+2] +
                  coords[3*(side[2]-1)+2])/3.0;
        uv[0]  = (parms[2*(side[0]-1)  ]  + parms[2*(side[1]-1)  ] +
                  parms[2*(side[2]-1)  ])/3.0;
        uv[1]  = (parms[2*(side[0]-1)+1]  + parms[2*(side[1]-1)+1] +
                  parms[2*(side[2]-1)+1])/3.0;
        EG_getInterior(faces[i], xyz, uv);
      }
      stat = EG_evaluate(faces[i], uv, result);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_evaluate %d %d = %d (EG_quadTess)!\n",
                 i+1, j+1, stat);
        EG_free(table);
        EG_free(mid);
        EG_free(tris);
        EG_free(coords);
        EG_free(faces);
        EG_free(edges);
        EG_deleteObject(newTess);
        return stat;
      }
      coords[3*iv  ] = result[0];
      coords[3*iv+1] = result[1];
      coords[3*iv+2] = result[2];
      parms[2*iv  ]  = uv[0];
      parms[2*iv+1]  = uv[1];
      iv++;

      tris[18*j   ] = trs[3*j  ];
      tris[18*j+ 1] = side[2];
      tris[18*j+ 2] = iv;
      tris[18*j+ 3] = trs[3*j  ];
      tris[18*j+ 4] = iv;
      tris[18*j+ 5] = side[1];

      tris[18*j+ 6] = trs[3*j+1];
      tris[18*j+ 7] = side[0];
      tris[18*j+ 8] = iv;
      tris[18*j+ 9] = trs[3*j+1];
      tris[18*j+10] = iv;
      tris[18*j+11] = side[2];

      tris[18*j+12] = trs[3*j+2];
      tris[18*j+13] = side[1];
      tris[18*j+14] = iv;
      tris[18*j+15] = trs[3*j+2];
      tris[18*j+16] = iv;
      tris[18*j+17] = side[0];
    }

    stat = EG_setTessFace(newTess, i+1, iv, coords, parms, 6*nt, tris);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessFace %d = %d (EG_quadTess)!\n",
               i+1, stat);
      EG_free(table);
      EG_free(mid);
      EG_free(tris);
      EG_free(coords);
      EG_free(faces);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
#ifdef TRIOUT
    sprintf(filename, "Components.%d.i.tri", i+1);
    fp = fopen(filename, "w");
    fprintf(fp," %d %d\n", iv, 6*nt);
    for (j = 0; j < iv; j++)
      fprintf(fp, " %lf %lf %lf\n",coords[3*j  ],coords[3*j+1],coords[3*j+2]);
    for (j = 0; j < 6*nt; j++)
      fprintf(fp, " %d %d %d\n", tris[3*j  ], tris[3*j+1], tris[3*j+2]);
    for (j = 0; j < 6*nt; j++)
      fprintf(fp, " 1\n");
    fclose(fp);
#endif
  }

  /* fill in the quads from TFI Faces */

  for (i = 0; i < nfaces; i++) {
    if (btess->tess2d[i].tfi != 1) continue;
    stat = EG_getTessFace(tess, i+1, &np, &xyzs, &uvs, &ptype, &pindex,
                          &nt, &trs, &trc);
    if (stat != EGADS_SUCCESS) continue;

    /* find degenerate nodes (if any) */
    degens[0] = degens[1] = 0;
    iuv[0]    = iuv[1]    = 0;
    stat = EG_getBodyTopos(obj, faces[i], EDGE, &k, &objs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Internal: EG_getBodyTopos on Face %d = %d\n", i+1, stat);
    } else {
      for (j = 0; j < k; j++) {
        stat = EG_getTopology(objs[j], &geom, &oclass, &mtype,
                              trange, &n, &nodes, &senses);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Internal: EG_getTopology on Edge = %d\n", stat);
          continue;
        }
        if (mtype != DEGENERATE) continue;
        stat = EG_getEdgeUVeval(faces[i], objs[j], 0, trange[0], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Internal: EG_getEdgeUVeval = %d\n", stat);
          continue;
        }
        n = EG_indexBodyTopo(obj, nodes[0]);
        if (n > 0) {
          if (degens[0] == 0) {
            degens[0] = n;
            if (result[3] != 0.0) iuv[0] = 1;
          } else if (degens[1] == 0) {
            degens[1] = n;
            if (result[3] != 0.0) iuv[1] = 1;
          } else {
            printf(" EGADS Info: More than 2 Degen Nodes in Face %d!\n",
                   i+1);
          }
        }
      }
      EG_free(objs);
    }

    /* make the vertices */
    for (j = 0; j < np; j++) {
      coords[3*j  ] = xyzs[3*j  ];
      coords[3*j+1] = xyzs[3*j+1];
      coords[3*j+2] = xyzs[3*j+2];
      parms[2*j  ]  = uvs[2*j  ];
      parms[2*j+1]  = uvs[2*j+1];
      table[j]      = NOTFILLED;
    }
    iv = np;

    /* get quad side midpoint insertions */
    for (is = j = 0; j < nt/2; j++)
      for (k = 0; k < 4; k++) {
        i1 = trs[6*j+sideq[k][0]] - 1;
        i2 = trs[6*j+sideq[k][1]] - 1;
        if (i2 < i1) {
          stat = i1;
          i1   = i2;
          i2   = stat;
        }
        m = EG_findMidSide(i1, i2, table, mid);
        if (m > 0) continue;
        ie = trc[6*j+neigq[k]];
        if (ie > 0) {
          /* Interior side */
          uvm[0] = uvs[2*i1  ];
          uvm[1] = uvs[2*i1+1];
          uvp[0] = uvs[2*i2  ];
          uvp[1] = uvs[2*i2+1];
          if ((ptype[i1] == 0) && (pindex[i1] == degens[0])) {
            if (iuv[0] == 0) {
              uvm[0] = uvs[2*i2  ];
            } else {
              uvm[1] = uvs[2*i2+1];
            }
          } else if ((ptype[i1] == 0) && (pindex[i1] == degens[1])) {
            if (iuv[1] == 0) {
              uvm[0] = uvs[2*i2  ];
            } else {
              uvm[1] = uvs[2*i2+1];
            }
          } else if ((ptype[i2] == 0) && (pindex[i2] == degens[0])) {
            if (iuv[0] == 0) {
              uvp[0] = uvs[2*i1  ];
            } else {
              uvp[1] = uvs[2*i1+1];
            }
          } else if ((ptype[i2] == 0) && (pindex[i2] == degens[1])) {
            if (iuv[1] == 0) {
              uvp[0] = uvs[2*i1  ];
            } else {
              uvp[1] = uvs[2*i1+1];
            }
          }
          /* could add the opposing triangles */
          EG_getSidepoint(faces[i], 0.5, uvm, uvp, NULL, NULL, uv);
          stat = EG_evaluate(faces[i], uv, result);
        } else {
          /* Edge side */
          ie   = -ie - 1;
          stat = EG_getTopology(edges[ie], &geom, &oclass, &mtype,
                                trange, &n, &objs, &senses);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: EG_getTopology %d/%d = %d (EG_quadTess)!\n",
                     i+1, ie+1, stat);
            EG_free(table);
            EG_free(mid);
            EG_free(tris);
            EG_free(coords);
            EG_free(faces);
            EG_free(edges);
            EG_deleteObject(newTess);
            return stat;
          }
          n = EG_indexBodyTopo(obj, objs[0]);
          if ((ptype[i1] == 0) && (ptype[i2] == 0)) {
            ien = 1;
          } else if (ptype[i1] == 0) {
            if (pindex[i1] == n) {
              ien = 1;
            } else {
              ien = 2*ptype[i2] - 1;
            }
            if (ie+1 != pindex[i2])
              printf(" EGADS Info: Edge mismatch %d %d\n", ie+1, pindex[i2]);
          } else if (ptype[i2] == 0) {
            if (pindex[i2] == n) {
              ien = 1;
            } else {
              ien = 2*ptype[i1] - 1;
            }
            if (ie+1 != pindex[i1])
              printf(" EGADS Info: Edge mismatch %d %d\n", ie+1, pindex[i1]);
          } else {
            ien = 2*ptype[i1] - 1;
            if (ien > 2*ptype[i2]-1) ien = 2*ptype[i2] - 1;
            if (ien < 0) {
              printf(" EGADS Internal: EG_quadTess ien = %d!\n", ien);
              continue;
            }
            if ((ie+1 != pindex[i1]) || (ie+1 != pindex[i2]))
              printf(" EGADS Info: Edge mismatch %d %d %d\n", ie+1, pindex[i1],
                     pindex[i2]);
          }
          t    = ntess->tess1d[ie].t[ien];
          stat = EG_getEdgeUV(faces[i], edges[ie], 0, t, uv);
          if (stat == EGADS_TOPOERR) {
            /* sense in Face twice! */
            stat = EG_getEdgeUV(faces[i], edges[ie], -1, t, uvm);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Info:  EG_getEdgeUV -1\n");
            } else {
              stat = EG_getEdgeUV(faces[i], edges[ie],  1, t, uvp);
              if (stat != EGADS_SUCCESS) {
                printf(" EGADS Info:  EG_getEdgeUV +1\n");
              } else {
                result[0] = sqrt((uvm[0]-uv[0])*(uvm[0]-uv[0]) +
                                 (uvm[1]-uv[1])*(uvm[1]-uv[1]));
                result[1] = sqrt((uvp[0]-uv[0])*(uvp[0]-uv[0]) +
                                 (uvp[1]-uv[1])*(uvp[1]-uv[1]));
                if (result[0] < result[1]) {
                  uv[0] = uvm[0];
                  uv[1] = uvm[1];
                } else {
                  uv[0] = uvp[0];
                  uv[1] = uvp[1];
                }
              }
            }
          }
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: EG_getEdgeUV %d/%d = %d (EG_quadTess)!\n",
                     i+1, ie+1, stat);
            EG_free(table);
            EG_free(mid);
            EG_free(tris);
            EG_free(coords);
            EG_free(faces);
            EG_free(edges);
            EG_deleteObject(newTess);
            return stat;
          }
          result[0] = ntess->tess1d[ie].xyz[3*ien  ];
          result[1] = ntess->tess1d[ie].xyz[3*ien+1];
          result[2] = ntess->tess1d[ie].xyz[3*ien+2];
        }
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_evaluate %d %d/%d = %d (EG_quadTess)!\n",
                   i+1, j+1, k+1, stat);
          EG_free(table);
          EG_free(mid);
          EG_free(tris);
          EG_free(coords);
          EG_free(faces);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
        }
        if (m == 0) {
          table[i1]      = is;
        } else {
          mid[-m-1].next = is;
        }
        mid[is].vert2 = i2;
        mid[is].nvert = iv + 1;
        mid[is].next  = NOTFILLED;
        is++;
        coords[3*iv  ] = result[0];
        coords[3*iv+1] = result[1];
        coords[3*iv+2] = result[2];
        parms[2*iv  ]  = uv[0];
        parms[2*iv+1]  = uv[1];
        iv++;
      }

    /* fill in the new triangulation based on Quads */
    for (j = 0; j < nt/2; j++) {
      for (k = 0; k < 4; k++) {
        i1 = trs[6*j+sideq[k][0]] - 1;
        i2 = trs[6*j+sideq[k][1]] - 1;
        if (i2 < i1) {
          stat = i1;
          i1   = i2;
          i2   = stat;
        }
        side[k] = EG_findMidSide(i1, i2, table, mid);
        if (side[k] <= 0) {
          if (outLevel > 0)
            printf(" EGADS Error: findMidSide %d = %d (EG_quadTess)!\n",
                   i+1, side[k]);
          EG_free(table);
          EG_free(mid);
          EG_free(tris);
          EG_free(coords);
          EG_free(faces);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
        }
      }
      /* get middle vertex based on mid-sides */
      EG_minArc4(faces[i], 0.5, 0.5,
                 &parms[2*(side[3]-1)], &parms[2*(side[0]-1)],
                 &parms[2*(side[1]-1)], &parms[2*(side[2]-1)], uv);
      stat  = EG_evaluate(faces[i], uv, result);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_evaluate %d %d = %d (EG_quadTess)!\n",
                 i+1, j+1, stat);
        EG_free(table);
        EG_free(mid);
        EG_free(tris);
        EG_free(coords);
        EG_free(faces);
        EG_free(edges);
        EG_deleteObject(newTess);
        return stat;
      }
      coords[3*iv  ] = result[0];
      coords[3*iv+1] = result[1];
      coords[3*iv+2] = result[2];
      parms[2*iv  ]  = uv[0];
      parms[2*iv+1]  = uv[1];
      iv++;

      tris[24*j   ] = trs[6*j  ];
      tris[24*j+ 1] = side[3];
      tris[24*j+ 2] = iv;
      tris[24*j+ 3] = trs[6*j  ];
      tris[24*j+ 4] = iv;
      tris[24*j+ 5] = side[2];

      tris[24*j+ 6] = trs[6*j+1];
      tris[24*j+ 7] = side[0];
      tris[24*j+ 8] = iv;
      tris[24*j+ 9] = trs[6*j+1];
      tris[24*j+10] = iv;
      tris[24*j+11] = side[3];

      tris[24*j+12] = trs[6*j+2];
      tris[24*j+13] = side[1];
      tris[24*j+14] = iv;
      tris[24*j+15] = trs[6*j+2];
      tris[24*j+16] = iv;
      tris[24*j+17] = side[0];

      tris[24*j+18] = trs[6*j+5];
      tris[24*j+19] = side[2];
      tris[24*j+20] = iv;
      tris[24*j+21] = trs[6*j+5];
      tris[24*j+22] = iv;
      tris[24*j+23] = side[1];
    }

    stat = EG_setTessFace(newTess, i+1, iv, coords, parms, 4*nt, tris);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessFace %d = %d (EG_quadTess)!\n",
               i+1, stat);
      EG_free(table);
      EG_free(mid);
      EG_free(tris);
      EG_free(coords);
      EG_free(faces);
      EG_free(edges);
      EG_deleteObject(newTess);
      return stat;
    }
    /* set tfi flag in new tessellation */
    ntess->tess2d[i].tfi = 1;
#ifdef TRIOUT
    sprintf(filename, "Components.%d.i.tri", i+1);
    fp = fopen(filename, "w");
    fprintf(fp," %d %d\n", iv, 4*nt);
    for (j = 0; j < iv; j++)
      fprintf(fp, " %lf %lf %lf\n",coords[3*j  ],coords[3*j+1],coords[3*j+2]);
    for (j = 0; j < 4*nt; j++)
      fprintf(fp, " %d %d %d\n", tris[3*j  ], tris[3*j+1], tris[3*j+2]);
    for (j = 0; j < 4*nt; j++)
      fprintf(fp, " 1\n");
    fclose(fp);
#endif
  }

  /* clean up temp storage */
  EG_free(table);
  EG_free(mid);
  EG_free(tris);
  EG_free(coords);
  EG_free(faces);
  EG_free(edges);

  /* close up the open tessellation */
  stat = EG_statusTessBody(newTess, &obj, &i, &npts);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_statusTessBody = %d (EG_quadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }
  if (i != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: New Tessellation Object is Open (EG_quadTess)!\n");
    EG_deleteObject(newTess);
    return EGADS_TESSTATE;
  }
#ifndef LITE
  table = (int *) EG_alloc(nfaces*sizeof(int));
  if (table != NULL) {
    for (i = 0; i < nfaces; i++) table[i] = ntess->tess2d[i].ntris/2;
    stat = EG_attributeAdd(newTess, ".mixed", ATTRINT, nfaces,
                           table, NULL, NULL);
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Warning: EG_attributeAdd m = %d !\n", stat);
    EG_free(table);
  }
  stat = EG_attributeAdd(newTess, ".tessType", ATTRSTRING, 4,
                         NULL, NULL, "Quad");
  if (stat != EGADS_SUCCESS)
    if (outLevel > 0)
      printf(" EGADS Warning: EG_attributeAdd Q = %d (EG_quadTess)!\n", stat);
#endif

  /* do we regularize? */
  stat = EG_attributeRet(tess, ".qRegular", &atype, &alen, &ints, &reals, &str);
  if (stat == EGADS_SUCCESS)
    if ((atype == ATTRSTRING) && (str != NULL))
#ifdef __CUDACC__
      if ((EG_strncmp(str,"OFF",3) == 0) || (EG_strncmp(str,"Off",3) == 0) ||
          (EG_strncmp(str,"off",3) == 0)) {
#else
      if ((strcmp(str,"OFF") == 0) || (strcmp(str,"Off") == 0) ||
          (strcmp(str,"off") == 0)) {
#endif
        *quadTess = newTess;
        return EGADS_SUCCESS;
      }
#ifndef REGULAR
  if (stat == EGADS_NOTFOUND) {
    *quadTess = newTess;
    return EGADS_SUCCESS;
  }
#endif

  /* yes! */
  bodydata.tess = newTess;
  stat = EG_getBodyTopos(btess->src, NULL, EDGE, &bodydata.nedges, NULL);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_getBodyTopos E = %d (EG_quadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }
  stat = EG_getBodyTopos(btess->src, NULL, FACE, &bodydata.nfaces,
                         &bodydata.faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_getBodyTopos F = %d (EG_quadTess)!\n", stat);
    EG_deleteObject(newTess);
    return stat;
  }
  stat = EG_createMeshMap(&bodydata);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_createMeshMap = %d (EG_quadTess)!\n", stat);
    EG_free(bodydata.faces);
    *quadTess = newTess;
    return EGADS_SUCCESS;
  }

  /* set the thread storage */
  qthread.mutex    = NULL;
  qthread.master   = EMP_ThreadID();
  qthread.end      = bodydata.nfaces;
  qthread.index    = 0;
  qthread.ntess    = ntess;
  qthread.bodydata = &bodydata;

  np = EMP_Init(&start);
  if (outLevel > 1) printf(" EMP NumProcs = %d!\n", np);

  if (np > 1) {
    /* create the mutex to handle list synchronization */
    qthread.mutex = EMP_LockCreate();
    if (qthread.mutex == NULL) {
      printf(" EMP Error: mutex creation = NULL!\n");
      np = 1;
    } else {
      /* get storage for our extra threads */
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(qthread.mutex);
        np = 1;
      }
    }
  }

  /* create the threads and get going! */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(EG_quadThread, &qthread);
      if (threads[i] == NULL)
        printf(" EMP Error Creating Thread #%d!\n", i+1);
    }
  /* now run the thread block from the original thread */
  EG_quadThread(&qthread);

  /* wait for all others to return */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);

  /* thread cleanup */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);
  if (qthread.mutex != NULL) EMP_LockDestroy(qthread.mutex);
  if (threads != NULL) free(threads);
  if (outLevel > 1)
    printf(" EMP Number of Seconds on Quad Thread Block = %ld\n",
           EMP_Done(&start));

  /* collect all of the Face quads and make the tessellation object */
  stat = EG_makeQuadTess(bodydata, &obj);
  EG_deleteObject(newTess);
  EG_destroyMeshMap(&bodydata);
  EG_free(bodydata.faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: EG_makeQuadTess = %d (EG_quadTess)!\n", stat);
    return stat;
  }

  *quadTess = obj;
  return EGADS_SUCCESS;
}
