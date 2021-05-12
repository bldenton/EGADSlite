/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             produce a limitted body tessellation
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef STANDALONE
#include "egads.h"
/* undocumented internal functions from egadsInternals.h */
extern /*@kept@*/ /*@null@*/ egObject *EG_context( const egObject *object );
extern int  EG_outLevel( const egObject *object );
extern int  EG_makeObject( /*@null@*/ egObject *context, egObject **obj );
extern int  EG_referenceObject( egObject *object,
                                /*@null@*/ const egObject *ref );
extern int  EG_referenceTopObj( egObject *object,
                                /*@null@*/ const egObject *ref );
#else
#include "egadsTypes.h"
#include "egadsInternals.h"
extern int  EG_getBodyTopos( const egObject *body, /*@null@*/ egObject *src,
                             int oclass, int *ntopo, egObject ***topos );
extern int  EG_getTopology( const egObject *topo, egObject **geom, int *oclass,
                            int *type, /*@null@*/ double *limits,
                            int *nChildren, egObject ***children, int **sense );
extern int  EG_indexBodyTopo( const egObject *body, const egObject *src );
extern int  EG_evaluate( const egObject *geom, const double *param,
                         double *results );
#endif
#include "egadsTris.h"
#include "emp.h"


/* undocumented functions from egadsTess.c */
extern int  EG_fillTris(egObject *body, int iFace, egObject *face,
                        egObject *tess, triStruct *ts, fillArea *fa, long tID);
extern void EG_cleanupTess(egTessel *btess);



static int
EG_limitEdges(egTessel *btess, int *elimits)
{
  int      i, j, k, n, npts, stat, outLevel, nedge, oclass, mtype, nnode;
  int      nf, nface, nloop, ntype, ndum, *senses, *finds;
  double   xyz[MAXELEN][3], t[MAXELEN], aux[MAXELEN][3], mindist, dist;
  double   d, dotnrm, dot, limits[4], result[18];
  egObject *body, *geom, *ref, **faces, **loops, **edges, **nodes, **dum;

  body     = btess->src;
  outLevel = EG_outLevel(body);
  dotnrm   = -0.000001;               /* make sure we are atleast 90 degrees */
  
  stat = EG_getBodyTopos(body, NULL, EDGE, &nedge, &edges);
  if (stat != EGADS_SUCCESS) return stat;
  stat = EG_getBodyTopos(body, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) return stat;
  
  btess->tess1d = (egTess1D *) EG_alloc(nedge*sizeof(egTess1D));
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Alloc %d Edges (EG_limitEdges)!\n", nedge);
    EG_free(faces);
    EG_free(edges);
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
  btess->nEdge = nedge;
  
  /* get the face indices (if any) */
  for (i = 0; i < nface; i++) {
    stat = EG_getTopology(faces[i], &geom, &oclass, &mtype, limits,
                          &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (j = 0; j < nloop; j++) {
      stat = EG_getTopology(loops[j], &geom, &oclass, &mtype, limits,
                            &ndum, &dum, &senses);
      if (stat != EGADS_SUCCESS) continue;
      for (k = 0; k < ndum; k++) {
        n = EG_indexBodyTopo(body, dum[k]);
        if (n <= EGADS_SUCCESS) continue;
        if (senses[k] < 0) {
          if (btess->tess1d[n-1].faces[0].nface != 0) {
            if (btess->tess1d[n-1].faces[0].nface == 1) {
              btess->tess1d[n-1].faces[0].faces = (int *) EG_alloc(2*sizeof(int));
              if (btess->tess1d[n-1].faces[0].faces == NULL) {
                if (outLevel > 0)
                  printf(" EGADS Error: Alloc (-) Edge %d (EG_limitEdges)!\n", n);
                EG_free(faces);
                EG_free(edges);
                return EGADS_MALLOC;
              }
              btess->tess1d[n-1].faces[0].faces[0] = btess->tess1d[n-1].faces[0].index;
              btess->tess1d[n-1].faces[0].faces[1] = i+1;
            } else {
              finds = (int *) EG_reall( btess->tess1d[n-1].faces[0].faces,
                                       (btess->tess1d[n-1].faces[0].nface+1)*sizeof(int));
              if (finds == NULL) {
                if (outLevel > 0)
                  printf(" EGADS Error: ReAlloc (-) Edge %d (EG_limitEdges)!\n", n);
                EG_free(faces);
                EG_free(edges);
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
                  printf(" EGADS Error: Alloc (+) Edge %d (EG_limitEdges)!\n", n);
                EG_free(faces);
                EG_free(edges);
                return EGADS_MALLOC;
              }
              btess->tess1d[n-1].faces[1].faces[0] = btess->tess1d[n-1].faces[1].index;
              btess->tess1d[n-1].faces[1].faces[1] = i+1;
            } else {
              finds = (int *) EG_reall( btess->tess1d[n-1].faces[1].faces,
                                       (btess->tess1d[n-1].faces[1].nface+1)*sizeof(int));
              if (finds == NULL) {
                if (outLevel > 0)
                  printf(" EGADS Error: ReAlloc (+) Edge %d (EG_limitEdges)!\n", n);
                EG_free(faces);
                EG_free(edges);
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
  
  /* report any non-manifold Edges */
  if (outLevel > 1)
    for (j = 0; j < nedge; j++) {
      if (btess->tess1d[j].faces[0].nface > 1) {
        printf(" EGADS Internal: Non-manifold Edge %d (-) with Faces", j+1);
        if (btess->tess1d[j].faces[0].faces != NULL)
          for (k = 0; k < btess->tess1d[j].faces[0].nface; k++)
            printf(" %d", btess->tess1d[j].faces[0].faces[k]);
        printf("!\n");
      }
      if (btess->tess1d[j].faces[1].nface > 1) {
        printf(" EGADS Internal: Non-manifold Edge %d (+) with Faces", j+1);
        if (btess->tess1d[j].faces[1].faces != NULL)
          for (k = 0; k < btess->tess1d[j].faces[1].nface; k++)
            printf(" %d", btess->tess1d[j].faces[1].faces[k]);
        printf("!\n");
      }
    }

  /* do the Edges -- one at a time */
  
  for (j = 0; j < nedge; j++) {
    stat = EG_getTopology(edges[j], &geom, &oclass, &mtype, limits,
                          &nnode, &nodes, &senses);
    if (stat != EGADS_SUCCESS) {
      EG_free(faces);
      EG_free(edges);
      return stat;
    }
          
    /* set end points */
    stat = EG_getTopology(nodes[0], &ref, &oclass, &ntype, xyz[0],
                          &ndum, &dum, &senses);
    if (stat != EGADS_SUCCESS) {
      EG_free(faces);
      EG_free(edges);
      return stat;
    }
    npts      = 2;
    t[0]      = limits[0];
    xyz[1][0] = xyz[0][0];
    xyz[1][1] = xyz[0][1];
    xyz[1][2] = xyz[0][2];
    t[1]      = limits[1];
    btess->tess1d[j].nodes[0] = EG_indexBodyTopo(body, nodes[0]);
    btess->tess1d[j].nodes[1] = btess->tess1d[j].nodes[0];
    if (mtype == TWONODE) {
      stat = EG_getTopology(nodes[1], &ref, &oclass, &ntype, xyz[1],
                            &ndum, &dum, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(faces);
        EG_free(edges);
        return stat;
      }
      btess->tess1d[j].nodes[1] = EG_indexBodyTopo(body, nodes[1]);
    }

    /* degenerate -- finish up */
    if (mtype == DEGENERATE) {
      btess->tess1d[j].xyz = (double *) EG_alloc(3*npts*sizeof(double));
      if (btess->tess1d[j].xyz == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Alloc %d Pts Edge %d (EG_limitEdges)!\n",
                 npts, j+1);
        EG_free(faces);
        EG_free(edges);
        return EGADS_MALLOC;  
      }
      btess->tess1d[j].t = (double *) EG_alloc(npts*sizeof(double));
      if (btess->tess1d[j].t == NULL) {
        EG_free(btess->tess1d[j].xyz);
        btess->tess1d[j].xyz = NULL;
        if (outLevel > 0)
          printf(" EGADS Error: Alloc %d Ts Edge %d (EG_limitEdges)!\n",
                 npts, j+1);
        EG_free(faces);
        EG_free(edges);
        return EGADS_MALLOC;  
      }
      for (i = 0; i < npts; i++) {
        btess->tess1d[j].xyz[3*i  ] = xyz[i][0];
        btess->tess1d[j].xyz[3*i+1] = xyz[i][1];
        btess->tess1d[j].xyz[3*i+2] = xyz[i][2];
        btess->tess1d[j].t[i]       = t[i];
      }
      btess->tess1d[j].npts = npts;
      continue;
    }
    
    /* get minimum distance */
    stat = EG_evaluate(edges[j], &t[0], result);
    if (stat != EGADS_SUCCESS) {
      EG_free(faces);
      EG_free(edges);
      return stat;
    }
    mindist = (xyz[0][0]-result[0])*(xyz[0][0]-result[0]) +
              (xyz[0][1]-result[1])*(xyz[0][1]-result[1]) +
              (xyz[0][2]-result[2])*(xyz[0][2]-result[2]);
    stat = EG_evaluate(edges[j], &t[1], result);
    if (stat != EGADS_SUCCESS) {
      EG_free(faces);
      EG_free(edges);
      return stat;
    }
    dist = (xyz[1][0]-result[0])*(xyz[1][0]-result[0]) +
           (xyz[1][1]-result[1])*(xyz[1][1]-result[1]) +
           (xyz[1][2]-result[2])*(xyz[1][2]-result[2]);
    if (dist > mindist) mindist = dist;
    mindist = sqrt(mindist);
    
    /* periodic -- add a vertex */
    if (mtype == ONENODE) {
      xyz[2][0] = xyz[1][0];
      xyz[2][1] = xyz[1][1];
      xyz[2][2] = xyz[1][2];
      aux[2][0] = aux[1][0];
      aux[2][1] = aux[1][1];
      aux[2][2] = aux[1][2];
      t[2]      = t[1];
      t[1]      = 0.5*(t[0]+t[2]);
      stat      = EG_evaluate(edges[j], &t[1], result);
      if (stat != EGADS_SUCCESS) {
        EG_free(faces);
        EG_free(edges);
        return stat;
      }
      dist      = sqrt(result[3]*result[3] + result[4]*result[4] +
                       result[5]*result[5]);
      if (dist == 0) dist = 1.0;
      xyz[1][0] = result[0];
      xyz[1][1] = result[1];
      xyz[1][2] = result[2];
      aux[1][0] = result[3]/dist;
      aux[1][1] = result[4]/dist;
      aux[1][2] = result[5]/dist;
      npts      = 3;
    }

    /* non-linear curve types */    
    if (geom->mtype != LINE) {
      
      /* angle criteria - aux is normalized tangent */

      stat = EG_evaluate(edges[j], &t[0], result);
      if (stat != EGADS_SUCCESS) {
        EG_free(faces);
        EG_free(edges);
        return stat;
      }
      dist = sqrt(result[3]*result[3] + result[4]*result[4] +
                  result[5]*result[5]);
      if (dist == 0) dist = 1.0;
      aux[0][0] = result[3]/dist;
      aux[0][1] = result[4]/dist;
      aux[0][2] = result[5]/dist;
      stat = EG_evaluate(edges[j], &t[npts-1], result);
      if (stat != EGADS_SUCCESS) {
        EG_free(faces);
        EG_free(edges);
        return stat;
      }
      dist = sqrt(result[3]*result[3] + result[4]*result[4] +
                  result[5]*result[5]);
      if (dist == 0) dist = 1.0;
      aux[npts-1][0] = result[3]/dist;
      aux[npts-1][1] = result[4]/dist;
      aux[npts-1][2] = result[5]/dist;

      while (npts < MAXELEN) {
        /* find the segment with the largest angle */
        k   = -1;
        dot =  1.0;
        for (i = 0; i < npts-1; i++) {
          dist = (xyz[i][0]-xyz[i+1][0])*(xyz[i][0]-xyz[i+1][0]) +
                 (xyz[i][1]-xyz[i+1][1])*(xyz[i][1]-xyz[i+1][1]) +
                 (xyz[i][2]-xyz[i+1][2])*(xyz[i][2]-xyz[i+1][2]);
          if (dist < mindist*mindist) continue;
          d = aux[i][0]*aux[i+1][0] + aux[i][1]*aux[i+1][1] +
              aux[i][2]*aux[i+1][2];
          if (d < dot) {
            dot = d;
            k   = i;
          }
        }
        if ((dot > dotnrm) || (k == -1)) break;
        /* insert */
        for (i = npts-1; i > k; i--) {
          xyz[i+1][0] = xyz[i][0];
          xyz[i+1][1] = xyz[i][1];
          xyz[i+1][2] = xyz[i][2];
          aux[i+1][0] = aux[i][0];
          aux[i+1][1] = aux[i][1];
          aux[i+1][2] = aux[i][2];
          t[i+1]      = t[i];
        }
        t[k+1] = 0.5*(t[k]+t[k+2]);
        stat   = EG_evaluate(edges[j], &t[k+1], result);
        if (stat != EGADS_SUCCESS) {
          EG_free(faces);
          EG_free(edges);
          return stat;
        }
        dist   = sqrt(result[3]*result[3] + result[4]*result[4] +
                      result[5]*result[5]);
        if (dist == 0.0) dist = 1.0;
        xyz[k+1][0] = result[0];
        xyz[k+1][1] = result[1];
        xyz[k+1][2] = result[2];
        aux[k+1][0] = result[3]/dist;
        aux[k+1][1] = result[4]/dist;
        aux[k+1][2] = result[5]/dist;
        npts++;
      }
    }

    /* insert until we get to the point count */
    for (i = 0; i < npts-1; i++)
      aux[i][0] = (xyz[i][0]-xyz[i+1][0])*(xyz[i][0]-xyz[i+1][0]) +
                  (xyz[i][1]-xyz[i+1][1])*(xyz[i][1]-xyz[i+1][1]) +
                  (xyz[i][2]-xyz[i+1][2])*(xyz[i][2]-xyz[i+1][2]);
    aux[npts-1][0] = 0.0;
    while ((npts < MAXELEN) && (npts < elimits[j]+2)) {
      /* find the biggest segment */
      k    = 0;
      dist = aux[0][0];
      for (i = 1; i < npts-1; i++) {
        d = aux[i][0];
        if (d > dist) {
          dist = d;
          k    = i;
        }
      }
      /* insert */
      for (i = npts-1; i > k; i--) {
        xyz[i+1][0] = xyz[i][0];
        xyz[i+1][1] = xyz[i][1];
        xyz[i+1][2] = xyz[i][2];
        aux[i+1][0] = aux[i][0];
        t[i+1]      = t[i];
      }
      t[k+1] = 0.5*(t[k]+t[k+2]);
      stat   = EG_evaluate(edges[j], &t[k+1], result);
      if (stat != EGADS_SUCCESS) {
        EG_free(faces);
        EG_free(edges);
        return stat;
      }
      xyz[k+1][0] = result[0];
      xyz[k+1][1] = result[1];
      xyz[k+1][2] = result[2];
      npts++;
      d = (xyz[k][0]-xyz[k+1][0])*(xyz[k][0]-xyz[k+1][0]) +
          (xyz[k][1]-xyz[k+1][1])*(xyz[k][1]-xyz[k+1][1]) +
          (xyz[k][2]-xyz[k+1][2])*(xyz[k][2]-xyz[k+1][2]);
      aux[k][0] = d;
      d = (xyz[k+2][0]-xyz[k+1][0])*(xyz[k+2][0]-xyz[k+1][0]) +
          (xyz[k+2][1]-xyz[k+1][1])*(xyz[k+2][1]-xyz[k+1][1]) +
          (xyz[k+2][2]-xyz[k+1][2])*(xyz[k+2][2]-xyz[k+1][2]);
      aux[k+1][0] = d;
    }
    
    /* fill in the 1D structure */
    btess->tess1d[j].xyz = (double *) EG_alloc(3*npts*sizeof(double));
    if (btess->tess1d[j].xyz == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Alloc %d Pts Edge %d (EG_limitEdges)!\n",
               npts, j+1);
      EG_free(faces);
      EG_free(edges);
      return EGADS_MALLOC;  
    }
    btess->tess1d[j].t = (double *) EG_alloc(npts*sizeof(double));
    if (btess->tess1d[j].t == NULL) {
      EG_free(btess->tess1d[j].xyz);
      btess->tess1d[j].xyz = NULL;
      if (outLevel > 0)
        printf(" EGADS Error: Alloc %d Ts Edge %d (EG_limitEdges)!\n",
               npts, j+1);
      EG_free(faces);
      EG_free(edges);
      return EGADS_MALLOC;  
    }
    nf = btess->tess1d[j].faces[0].nface;
    if (nf > 0) {
      btess->tess1d[j].faces[0].tric = (int *) EG_alloc((nf*(npts-1))*sizeof(int));
      if (btess->tess1d[j].faces[0].tric == NULL) {
        EG_free(btess->tess1d[j].t);
        btess->tess1d[j].t   = NULL;
        EG_free(btess->tess1d[j].xyz);
        btess->tess1d[j].xyz = NULL;
        if (outLevel > 0)
          printf(" EGADS Error: Alloc %d Tric- Edge %d (EG_limitEdges)!\n",
                 npts, j+1);
        EG_free(faces);
        EG_free(edges);
        return EGADS_MALLOC;
      }
    }
    nf = btess->tess1d[j].faces[1].nface;
    if (nf > 0) {
      btess->tess1d[j].faces[1].tric = (int *) EG_alloc((nf*(npts-1))*sizeof(int));
      if (btess->tess1d[j].faces[1].tric == NULL) {
        if (btess->tess1d[j].faces[0].tric != NULL)
          EG_free(btess->tess1d[j].faces[0].tric);
        btess->tess1d[j].faces[0].tric = NULL;
        EG_free(btess->tess1d[j].t);
        btess->tess1d[j].t   = NULL;
        EG_free(btess->tess1d[j].xyz);
        btess->tess1d[j].xyz = NULL;
        if (outLevel > 0)
          printf(" EGADS Error: Alloc %d Tric+ Edge %d (EG_limitEdges)!\n",
                 npts, j+1);
        EG_free(faces);
        EG_free(edges);
        return EGADS_MALLOC;
      }
    }
    for (i = 0; i < npts; i++) {
      btess->tess1d[j].xyz[3*i  ] = xyz[i][0];
      btess->tess1d[j].xyz[3*i+1] = xyz[i][1];
      btess->tess1d[j].xyz[3*i+2] = xyz[i][2];
      btess->tess1d[j].t[i]       = t[i];
    }
    for (i = 0; i < npts-1; i++) {
      nf = btess->tess1d[j].faces[0].nface;
      for (k = 0; k < nf; k++)
        btess->tess1d[j].faces[0].tric[i*nf+k] = 0;
      nf = btess->tess1d[j].faces[1].nface;
      for (k = 0; k < nf; k++)
        btess->tess1d[j].faces[1].tric[i*nf+k] = 0;
    }
    btess->tess1d[j].npts = npts;
  }

  EG_free(faces);
  EG_free(edges);
  return EGADS_SUCCESS;
}


static void
EG_limitThread(void *struc)
{
  int       index, stat, *flimits;
  long      ID;
  triStruct tst;
  fillArea  fast;
  EMPtess   *tthread;
  
  tthread = (EMPtess *) struc;
  flimits = (int *) tthread->ptr;

  /* get our identifier */
  ID = EMP_ThreadID();

  tst.maxlen   = tthread->params[0];
  tst.chord    =  0.0;
  tst.dotnrm   = -1.0;
  tst.minlen   = tthread->tparam[0];
  tst.maxPts   = tthread->tparam[1];
  tst.qparm[0] = tthread->qparam[0];
  tst.qparm[1] = tthread->qparam[1];
  tst.qparm[2] = tthread->qparam[2];
  tst.mverts   = tst.nverts = 0;
  tst.verts    = NULL;
  tst.mtris    = tst.ntris  = 0;
  tst.tris     = NULL;
  tst.msegs    = tst.nsegs  = 0;
  tst.segs     = NULL;
  tst.mframe   = tst.nframe = 0;
  tst.frame    = NULL;
  tst.mloop    = tst.nloop  = 0;
  tst.loop     = NULL;
  tst.numElem  = -1;
  tst.hashTab  = NULL;
  
  fast.pts     = NULL;
  fast.segs    = NULL;
  fast.front   = NULL;
 
  /* look for work */
  for (;;) {
    
    /* only one thread at a time here -- controlled by a mutex! */
    if (tthread->mutex != NULL) EMP_LockSet(tthread->mutex);
    if (tthread->mark == NULL) {
      index = tthread->index;
    } else {
      for (index = tthread->index; index < tthread->end; index++) {
        if (tthread->mark[index] == 0) continue;
        break;
      }
    }
    tthread->index = index+1;
    if (tthread->mutex != NULL) EMP_LockRelease(tthread->mutex);
    if (index >= tthread->end) break;

    /* do the work */
    tst.maxPts = -flimits[index] - 1;
    stat = EG_fillTris(tthread->body, index+1, tthread->faces[index],
                       tthread->tess, &tst, &fast, ID);
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Warning: Face %d -> EG_fillTris = %d (EG_limitTessBody)!\n",
             index+1, stat);
  }
  
  /* exhausted all work -- cleanup & exit */
  if (tst.verts  != NULL) EG_free(tst.verts);
  if (tst.tris   != NULL) EG_free(tst.tris);
  if (tst.segs   != NULL) EG_free(tst.segs);
  if (tst.frame  != NULL) EG_free(tst.frame);
  if (tst.loop   != NULL) EG_free(tst.loop);
  
  if (fast.segs  != NULL) EG_free(fast.segs);
  if (fast.pts   != NULL) EG_free(fast.pts);
  if (fast.front != NULL) EG_free(fast.front);
  
  if (ID != tthread->master) EMP_ThreadExit();
}


/* generates a simple, limited tessellation of the Body (not for WireBodies)
 *
 * where: object  - the Body object to tessellate
 *        elimits - an integer array of sizes (one per Edge in the Body)
 *                  this number indicates the number of vertex points in the
 *                  interior of the Edge (the bounding Nodes are not counted).
 *        flimits - an integer array of sizes (one per Face in the Body)
 *                  this number is the target for the interior vertices placed
 *                  in the triangulation.
 *         tess   - the returned Tessellation object
 *
 * note: the actual number of vertices generated on an object may be larger than
 *       the request if the scheme finds the request is too small to coarsely
 *       represent the object.
 */

int
EG_limitTessBody(egObject *object, int *elimits, int *flimits, egObject **tess)
{
  int      i, j, stat, outLevel, nface, np;
  double   params[3] = {1.e-6, 0.0, 0.0};
  void     **threads = NULL;
  long     start;
  egTessel *btess;
  egObject *ttess, *context, **faces;
  egCntxt  *cntx;
  EMPtess  tthread;

  *tess = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass != BODY)       return EGADS_NOTBODY;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);
  if (context == NULL)              return EGADS_NULLOBJ;
  cntx     = (egCntxt *) context->blind;
  
  btess = (egTessel *) EG_alloc(sizeof(egTessel));
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Blind Malloc (EG_limitTessBody)!\n");
    return EGADS_MALLOC;
  }
  btess->src       = object;
  btess->xyzs      = NULL;
  btess->tess1d    = NULL;
  btess->tess2d    = NULL;
  btess->globals   = NULL;
  btess->nGlobal   = 0;
  btess->nEdge     = 0;
  btess->nFace     = 0;
  btess->nu        = 0;
  btess->nv        = 0;
  btess->done      = 1;
  btess->params[0] = params[0];
  btess->params[1] = params[1];
  btess->params[2] = params[2];
  for (i = 0; i < MTESSPARAM; i++) btess->tparam[i] = 0.0;
  if (cntx != NULL)
    for (i = 0; i < MTESSPARAM; i++) btess->tparam[i] = cntx->tess[i];
  
  /* do the Edges & make the Tessellation Object */
  
  stat = EG_limitEdges(btess, elimits);
  if (stat != EGADS_SUCCESS) {
    EG_cleanupTess(btess);
    EG_free(btess);
    return stat;  
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
  
  /* Wire Body */
  if (object->mtype == WIREBODY) return EGADS_SUCCESS;

  /* not a WireBody */

  stat = EG_getBodyTopos(object, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getBodyTopos = %d (EG_limitTessBody)!\n",
           stat);
    EG_deleteObject(ttess);
    *tess = NULL;
    return stat;
  }
  btess->tess2d = (egTess2D *) EG_alloc(2*nface*sizeof(egTess2D));
 if (btess->tess2d == NULL) {
    printf(" EGADS Error: Alloc %d Faces (EG_limitTessBody)!\n", nface);  
    EG_deleteObject(ttess);
    *tess = NULL;
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
  }
  btess->nFace = nface;
  
  /* set up for explicit multithreading */
  tthread.mutex     = NULL;
  tthread.master    = EMP_ThreadID();
  tthread.index     = 0;
  tthread.end       = nface;
  tthread.mark      = NULL;
  tthread.tess      = ttess;
  tthread.body      = object;
  tthread.faces     = faces;
  tthread.params    = params;
  tthread.tparam    = btess->tparam;
  tthread.qparam[0] = -1.0;                /* no quadding */
  tthread.qparam[1] = tthread.qparam[2] = 0.0;
  tthread.ptr       = flimits;
  
  np = EMP_Init(&start);
  if (outLevel > 1) printf("EMP NumProcs = %d!\n", np);
  
  if (np > 1) {
    /* create the mutex to handle list synchronization */
    tthread.mutex = EMP_LockCreate();
    if (tthread.mutex == NULL) {
      printf(" EMP Error: mutex creation = NULL!\n");
      np = 1;
    } else {
      /* get storage for our extra threads */
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(tthread.mutex);
        np = 1;
      }
    }
  }

  /* create the threads and get going! */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(EG_limitThread, &tthread);
      if (threads[i] == NULL)
        printf(" EMP Error Creating Thread #%d!\n", i+1);
    }
  /* now run the thread block from the original thread */
  EG_limitThread(&tthread);
  
  /* wait for all others to return */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);

  /* cleanup */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);
  if (tthread.mutex != NULL) EMP_LockDestroy(tthread.mutex);
  if (threads != NULL) free(threads);
  EG_free(faces);
  if (outLevel > 0)
    printf("EMP Number of Seconds on Thread Block = %ld\n", EMP_Done(&start));

  return EGADS_SUCCESS;
}


#ifdef STANDALONE

int main(int argc, char *argv[])
{
  int          i, j, n, stat, oclass, mtype, nbody, ibody, nfaces, nedges, per;
  int          *elimit, *flimit, *senses, ntri;
  double       box[6], trange[2], size, alen, data[14];
  const int    *tris, *tric, *ptype, *pindex;
  const char   *OCCrev;
  const double *xyzs, *parms;
  ego          context, model, geom, tess, *bodies, *faces, *edges, *dum;
  
  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %1d.%02d with %s\n\n", i, j, OCCrev);

  if (argc != 2) {
    printf("\n Usage: limit filename\n\n");
    return 1;
  }
  
  /* initialize */
  printf(" EG_open           = %d\n", EG_open(&context));
  printf(" EG_loadModel      = %d\n", EG_loadModel(context, 0, argv[1], &model));
  printf(" EG_getBoundingBox = %d\n", EG_getBoundingBox(model, box));
  printf("       BoundingBox = %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("                     %lf %lf %lf\n", box[3], box[4], box[5]);
  printf(" \n");

                            size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];
  size *= 0.05;
  
  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                        &bodies, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n", stat);
    return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n", nbody);
  
  for (ibody = 0; ibody < nbody; ibody++) {
    stat = EG_getTopology(bodies[ibody], &geom, &oclass,
                          &mtype, NULL, &j, &dum, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" Body %d: EG_getTopology = %d!\n", ibody+1, stat);
      continue;
    }
    if (mtype == WIREBODY) {
      printf(" Body %d: Type = WireBody\n", ibody+1);
      continue;
    } else if (mtype == FACEBODY) {
      printf(" Body %d: Type = FaceBody\n", ibody+1);
    } else if (mtype == SHEETBODY) {
      printf(" Body %d: Type = SheetBody\n", ibody+1);
    } else {
      printf(" Body %d: Type = SolidBody\n", ibody+1);
    }
    
    stat = EG_getBodyTopos(bodies[ibody], NULL, FACE, &nfaces, &faces);
    i    = EG_getBodyTopos(bodies[ibody], NULL, EDGE, &nedges, &edges);
    if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
      printf(" EG_getBodyTopos Face = %d\n", stat);
      printf(" EG_getBodyTopos Edge = %d\n", i);
      if (faces != NULL) EG_free(faces);
      if (edges != NULL) EG_free(edges);
      continue;
    }
    
    elimit = (int *) malloc(nedges*sizeof(int));
    flimit = (int *) malloc(nfaces*sizeof(int));
    if ((elimit == NULL) || (flimit == NULL)) {
      printf(" Malloc error on %d Edges & %d Faces!\n", nedges, nfaces);
      if (flimit != NULL) free(flimit);
      if (elimit != NULL) free(elimit);
      EG_free(faces);
      EG_free(edges);
      continue;
    }
    
    /* set the edge counts */
    for (i = 0; i < nedges; i++) {
      elimit[i] = 0;
      stat = EG_getRange(edges[i], trange, &per);
      if (stat != EGADS_SUCCESS) {
        printf(" Edge %d: EG_getRange = %d!\n", i+1, stat);
        continue;
      }
      stat = EG_arcLength(edges[i], trange[0], trange[1], &alen);
      if (stat != EGADS_SUCCESS) {
        printf(" Edge %d: EG_arcLength = %d!\n", i+1, stat);
        continue;
      }
      elimit[i] = alen/size;
    }
    
    /* set the face counts */
    for (i = 0; i < nfaces; i++) {
      flimit[i] = 0;
      stat = EG_getMassProperties(faces[i], data);
      if (stat != EGADS_SUCCESS) {
        printf(" Face %d: EG_getMassProperties = %d!\n", i+1, stat);
        continue;
      }
      flimit[i] = data[1]/(size*size);
    }
    
    /* tessellate and return counts */
    stat = EG_limitTessBody(bodies[ibody], elimit, flimit, &tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_limitTessBody = %d!\n", stat);
    } else {
      
      printf("\n");
      for (i = 0; i < nedges; i++) {
        stat = EG_getTessEdge(tess, i+1, &n, &xyzs, &parms);
        if (stat != EGADS_SUCCESS) {
          printf(" Edge %d: EG_getTessEdge = %d!\n", i+1, stat);
          continue;
        }
        printf("   Edge %d:  %d %d\n", i+1, elimit[i], n-2);
      }
      
      printf("\n");
      for (i = 0; i < nfaces; i++) {
        stat = EG_getTessFace(tess, i+1, &n, &xyzs, &parms, &ptype, &pindex,
                              &ntri, &tris, &tric);
        if (stat != EGADS_SUCCESS) {
          printf(" Edge %d: EG_getTessFace = %d!\n", i+1, stat);
          continue;
        }
        for (per = j = 0; j < n; j++)
          if (ptype[j] < 0) per++;
        printf("   Face %d:  %d %d\n", i+1, flimit[i], per);
      }
      
      EG_deleteObject(tess);
    }
    printf("\n");

    free(flimit);
    free(elimit);
    EG_free(faces);
    EG_free(edges);
  }

  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}
#endif
