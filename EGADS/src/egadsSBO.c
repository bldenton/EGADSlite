/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Replacement Solid Boolean Operator Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsStack.h"


//#define DEBUG
//#define WRITERESULT

#define PARAMACC 1.0e-4         // parameter accuracy


  typedef struct {
    int      cnt;
    int      index;
    egObject *node;
  } nodeCnt;

  typedef struct {
    egObject *node;
    egObject *entity;
    /*@null@*/
    egObject *degen;
    int      type;
    int      match;
    double   param[2];
  } nodeInfo;

  typedef struct {
    egObject *node;
    egObject *edge;
    egObject *pcurve;
    int      sense;
    int      cnt;
    double   ts[16];
  } degenEdge;

  typedef struct {
    egObject *face;
    egObject *ent;
    egObject *pcurve;
    nodeInfo ni[2];
    int      next;
    int      sense;
    int      shell;
  } edgeInfo;

  typedef struct {
    egObject *edge;
    egObject *pcurve;
    int      sense;
    int      connect[2];
    egObject *nodes[2];
    double   range[2];
    double   dir[2][2];
    double   uvs[2][2];
  } edgeLoop;

  typedef struct {
    int      nSplit;
    egObject **ents;
  } splitEnts;


  extern int  EG_getBody( const egObject *obj, egObject **body );
  extern int  EG_outLevel( const egObject *object );
  extern int  EG_setOutLevel( const egObject *object, int level );
  extern int  EG_getGeometry( const egObject *geom, int *oclass, int *mtype,
                              egObject **refGeom, int **ivec, double **rvec );
  extern int  EG_getTopology( const egObject *topo, egObject **geom,
                              int *oclass, int *type,
                              /*@null@*/ double *limits, int *nChildren,
                              egObject ***children, int **senses );
  extern int  EG_makeTopology( egObject *context, /*@null@*/ egObject *geom,
                               int oclass, int mtype, /*@null@*/ double *limits,
                               int nChildren, /*@null@*/ egObject **children,
                               /*@null@*/ int *senses, egObject **topo );
  extern int  EG_getBodyTopos( const egObject *body, /*@null@*/ egObject *src,
                               int oclass, int *ntopo,
                               /*@null@*/ egObject ***topos );
  extern int  EG_indexBodyTopo( const egObject *body, const egObject *src );
  extern int  EG_getTolerance( const egObject *topo, double *tol );
  extern int  EG_tolerance( const egObject *topo, double *tol );
  extern int  EG_evaluate( const egObject *geom, const double *param,
                           double *result );
  extern int  EG_invEvaluate( const egObject *geom, double *xyz, double *param,
                              double *result );
  extern int  EG_otherCurve( const egObject *surface, const egObject *curve,
                             double tol, egObject **newcurve );
  extern int  EG_getRange( const egObject *geom, double *range, int *pflg );
  extern int  EG_getUVinfo( egObject *face, egObject *loop, double *box,
                            double *area );
#ifdef WRITERESULT
  extern int  EG_saveModel( const egObject *model, const char *name );
  extern int  EG_copyObject( const egObject *object, /*@null@*/ void *ptr,
                             egObject **copy );
#endif


static void
EG_splitInit(int num, splitEnts **splits)
{
  int       i;
  splitEnts *splitVec;
  
  *splits  = NULL;
  splitVec = (splitEnts *) EG_alloc(num*sizeof(splitEnts));
  if (splitVec == NULL) return;
  for (i = 0; i < num; i++) {
    splitVec[i].nSplit = 0;
    splitVec[i].ents   = NULL;
  }
  
  *splits = splitVec;
}


static void
EG_splitFree(int num, splitEnts **splits)
{
  int       i;
  splitEnts *splitVec;

  if (*splits == NULL) return;
  splitVec = *splits;
  *splits  = NULL;
  for (i = 0; i < num; i++)
    if (splitVec[i].ents != NULL) EG_free(splitVec[i].ents);
  EG_free(splitVec);
}


static int
EG_splitAlloc(int index, splitEnts *splits, int num)
{
  int i;
  
  if (splits[index].ents != NULL) {
    printf(" EGADS Internal: Data for index %d already exists!\n", index);
    return EGADS_EXISTS;
  }
  if (num <= 0) {
    printf(" EGADS Internal: numEntities for index %d is %d!\n", index, num);
    return EGADS_RANGERR;
  }
  
  splits[index].ents = (egObject **) EG_alloc(num*sizeof(egObject *));
  if (splits[index].ents == NULL) return EGADS_MALLOC;
  for (i = 0; i < num; i++) splits[index].ents[i] = NULL;
  splits[index].nSplit = num;
  
  return EGADS_SUCCESS;
}


static void
EG_fillNodeCnt(int i, int *eface, edgeInfo *fe, int *n, nodeCnt *nodec)
{
  int j, k, index, last, nnode = 0;
  
  *n    = 0;
  index = eface[i];
  if (index == -1) return;
  
  last = -1;
  while (index != -1) {
    if (fe[index].ni[0].type == FACE) {
      for (j = 0; j < nnode; j++)
        if (fe[index].ni[0].node == nodec[j].node) {
          nodec[j].cnt++;
          break;
        }
      if (j == nnode) {
        nodec[nnode].node  = fe[index].ni[0].node;
        nodec[nnode].cnt   = 1;
        nodec[nnode].index = last;
        nnode++;
      }
    }
    
    if (fe[index].ni[1].type == FACE) {
      for (j = 0; j < nnode; j++)
        if (fe[index].ni[1].node == nodec[j].node) {
          nodec[j].cnt++;
          break;
        }
      if (j == nnode) {
        nodec[nnode].node  = fe[index].ni[1].node;
        nodec[nnode].cnt   = 1;
        nodec[nnode].index = last;
        nnode++;
      }
    }
    
    last  = index;
    index = fe[index].next;
  }
  
  /* check for lone segment */
  for (j = 0; j < nnode-1; j++) {
    if (nodec[j].cnt != 1) continue;
    for (k = j+1; k < nnode; k++) {
      if (nodec[k].cnt != 1) continue;
      if (nodec[j].index != nodec[k].index) continue;
      nodec[k].cnt = 0;
      break;
    }
  }
  
  *n = nnode;
}


static int
EG_delOpenEdges(int nseg, int nfaces, int *eface, edgeInfo *fe)
{
  int     i, j, k, n, index, last;
  nodeCnt *nodec;
  
  nodec = (nodeCnt *) EG_alloc(2*nseg*sizeof(nodeCnt));
  if (nodec == NULL) return EGADS_MALLOC;
  
  for (i = 0; i < nfaces; i++) {
    EG_fillNodeCnt(i, eface, fe, &n, nodec);
/*  if (eface[i] != -1) printf("  Face %d:  n = %d\n", i+1, n);  */
    if (n == 0) continue;
    for (k = j = 0; j < n; j++) {
/*    printf("     node %d: count = %d\n", j+1, nodec[j].cnt);  */
      if (nodec[j].cnt != 1) continue;
      /* remove the segment -- patch up the linked list */
      last = nodec[j].index;
      if (last == -1) {
        index    = eface[i];
        eface[i] = fe[index].next;
      } else {
        index         = fe[last].next;
        fe[last].next = fe[index].next;
      }
      k++;
    }
    if (k        ==  0) continue;
    if (eface[i] == -1) continue;
    
    /* again! */
    i--;
  }

  EG_free(nodec);
  
  return EGADS_SUCCESS;
}


static double
EG_angleEdges(int fsense, double *dir1, double *dir2)
{
  double angle;
  
  angle  = atan2(-dir1[1], -dir1[0]) - atan2(dir2[1], dir2[0]);
  angle *= fsense;
  
  if  (angle < 0.0) angle += 2.0*PI;
  if ((angle < 0.0) || (angle > 2.0*PI))
    printf(" EGADS Internal: EG_angleEdges = %lf\n", angle);
  
  return angle;
}


static int
EG_matchEndpts(int nEdges, edgeLoop *el, int n, double *uv, double tol)
{
  int    i, m;
  double dist;
  
  for (m = i = 0; i < nEdges; i++) {
    if (el[i].connect[0] == -1) {
      dist = sqrt((uv[0]-el[i].uvs[0][0])*(uv[0]-el[i].uvs[0][0]) +
                  (uv[1]-el[i].uvs[0][1])*(uv[1]-el[i].uvs[0][1]));
      if (dist < tol) {
        el[i].connect[0] = n;
        m++;
/*      printf(" distance = %le (%le)\n", dist, tol);  */
      }
    }
    if (el[i].connect[1] == -1) {
      dist = sqrt((uv[0]-el[i].uvs[1][0])*(uv[0]-el[i].uvs[1][0]) +
                  (uv[1]-el[i].uvs[1][1])*(uv[1]-el[i].uvs[1][1]));
      if (dist < tol) {
        el[i].connect[1] = n;
        m++;
/*      printf(" distance = %le (%le)\n", dist, tol);  */
      }
    }
  }

  if (m%2 != 1) return EGADS_INDEXERR;
  return EGADS_SUCCESS;
}


static void
EG_adjPeriodic(double *uv, double *uvs, double *srange, int sper)
{
  
  if ((sper&1) != 0)
    if (fabs(uv[0]-uvs[0]) > 0.5*(srange[1]-srange[0])) {
#ifdef DEBUG
      printf(" adjPeriodic U: %lf %lf, period = %lf\n",
             uv[0], uvs[0], srange[1]-srange[0]);
#endif
      if (uvs[0] > uv[0]) {
        uvs[0] -= srange[1] - srange[0];
      } else {
        uvs[0] += srange[1] - srange[0];
      }
    }

  if ((sper&2) != 0)
    if (fabs(uv[1]-uvs[1]) > 0.5*(srange[3]-srange[2])) {
#ifdef DEBUG
      printf(" adjPeriodic V: %lf %lf, period = %lf\n",
             uv[1], uvs[1], srange[3]-srange[2]);
#endif
      if (uvs[1] > uv[1]) {
        uvs[1] -= srange[3] - srange[2];
      } else {
        uvs[1] += srange[3] - srange[2];
      }
    }

}


static void
EG_adjustPeriod(double *param, double *srange, int sper)
{
  double period;
  
  if ((sper&1) != 0) {
    period = srange[1] - srange[0];
    if ((param[0]+PARAMACC < srange[0]) || (param[0]-PARAMACC > srange[1]))
      if (param[0]+PARAMACC < srange[0]) {
        if (param[0]+period-PARAMACC < srange[1]) param[0] += period;
      } else {
        if (param[0]-period+PARAMACC > srange[0]) param[0] -= period;
      }
  }
  
  if ((sper&2) != 0) {
    period = srange[3] - srange[2];
    if ((param[1]+PARAMACC < srange[2]) || (param[1]-PARAMACC > srange[3]))
      if (param[1]+PARAMACC < srange[2]) {
        if (param[1]+period-PARAMACC < srange[3]) param[1] += period;
      } else {
        if (param[1]-period+PARAMACC > srange[2]) param[1] -= period;
      }
  }
  
}

static int
EG_traceLoops(egObject *context, const egObject *face, int nEdges,
              egObject **edges, int *senses, /*@null@*/ egObject **pcurves,
              int *nLoops, egObject ***loops)
{
  int      i, n, cnt, oclass, mtype, stat, fsense, outLevel, first, next, nNode;
  int      j0, j1, sper, *sens, *sen, *ints;
  double   ang, angle, tol, range[4], srange[4], result[9], u[3], v[3], dir[2];
  double   duv[4], uvs[2], uv[2], *pdata, *reals;
  egObject *geom, *surf, *ref, *sref, **nodes, **objs, **list;
  edgeLoop *el;
  nodeCnt  *nl;
  
  *nLoops  = 0;
  *loops   = NULL;
  el       = NULL;
  objs     = NULL;
  sens     = NULL;
  outLevel = EG_outLevel(context);
  
  /* look at Face */
  stat = EG_getTopology(face, &surf, &oclass, &fsense, range, &n, &nodes,
                        &sen);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: Cannot get Face = %d (EG_traceLoops)!\n", stat);
    return stat;
  }
  u[0] = u[1] = u[2] = 0.0;
  v[0] = v[1] = v[2] = 0.0;
  if (surf->mtype == PLANE) {
    if (pcurves != NULL) {
      printf(" EGADS Error: Planar Surface with PCurves (EG_traceLoops)!\n");
      return EGADS_TOPOERR;
    }
    stat = EG_getGeometry(surf, &oclass, &mtype, &ref, &sen, &pdata);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot get Plane = %d (EG_traceLoops)!\n", stat);
      return stat;
    }
    u[0] = pdata[3];
    u[1] = pdata[4];
    u[2] = pdata[5];
    v[0] = pdata[6];
    v[1] = pdata[7];
    v[2] = pdata[8];
    EG_free(pdata);
  } else {
    if (pcurves == NULL) {
      printf(" EGADS Error: nonPlanar Surf without PCurves (EG_traceLoops)!\n");
      return EGADS_TOPOERR;
    }
  }
  sref = surf;
  if (surf->mtype == TRIMMED) {
    stat = EG_getGeometry(surf, &oclass, &mtype, &sref, &ints, &reals);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot get Geometry = %d (EG_traceLoops)!\n", stat);
      return stat;
    }
    if (ints != NULL) EG_free(ints);
    EG_free(reals);
  }
  stat = EG_getRange(sref, srange, &sper);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: Cannot get Surf Range = %d (EG_traceLoops)!\n", stat);
    return stat;
  }
  
  /* save away the Edge information for processing */
  el = (edgeLoop *) EG_alloc(nEdges*sizeof(edgeLoop));
  if (el == NULL) return EGADS_MALLOC;

  for (i = 0; i < nEdges; i++) {
    stat = EG_getTopology(edges[i], &geom, &oclass, &mtype, range, &n,
                          &nodes, &sen);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot get Edge %d = %d (EG_traceLoops)!\n",
             i+1, stat);
      goto bail;
    }
    el[i].edge       = edges[i];
    el[i].sense      = senses[i];
    el[i].range[0]   = range[0];
    el[i].range[1]   = range[1];
    el[i].connect[0] = el[i].connect[1] = -1;
    el[i].pcurve     = NULL;
    if (pcurves != NULL) el[i].pcurve = pcurves[i];
    if (senses[i] == -1) {
      
      el[i].nodes[0] = el[i].nodes[1] = nodes[0];
      if (n == 2) el[i].nodes[0] = nodes[1];
      if (pcurves != NULL) {
        stat = EG_evaluate(pcurves[i], &range[0], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Cannot evaluate PCurve %d = %d (EG_traceLoops)!\n",
                 i+1, stat);
          goto bail;
        }
        EG_adjustPeriod(result, srange, sper);
        el[i].uvs[1][0] = result[0];
        el[i].uvs[1][1] = result[1];
        result[0]       = result[2];
        result[1]       = result[3];
      } else {
        stat = EG_evaluate(edges[i], &range[0], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Cannot evaluate Edge %d = %d (EG_traceLoops)!\n",
                 i+1, stat);
          goto bail;
        }
        el[i].uvs[1][0] = result[0]*u[0] + result[1]*u[1] + result[2]*u[2];
        el[i].uvs[1][1] = result[0]*v[0] + result[1]*v[1] + result[2]*v[2];
        result[0]       = result[3]*u[0] + result[4]*u[1] + result[5]*u[2];
        result[1]       = result[3]*v[0] + result[4]*v[1] + result[5]*v[2];
      }
      el[i].dir[1][0] = -result[0];
      el[i].dir[1][1] = -result[1];
      
      if (pcurves != NULL) {
        stat = EG_evaluate(pcurves[i], &range[1], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Cannot evaluate PCurve %d = %d (EG_traceLoops)!\n",
                 i+1, stat);
          goto bail;
        }
        EG_adjustPeriod(result, srange, sper);
        el[i].uvs[0][0] = result[0];
        el[i].uvs[0][1] = result[1];
        result[0]       = result[2];
        result[1]       = result[3];
      } else {
        stat = EG_evaluate(edges[i], &range[1], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Cannot evaluate Edge %d = %d (EG_traceLoops)!\n",
                 i+1, stat);
          goto bail;
        }
        el[i].uvs[0][0] = result[0]*u[0] + result[1]*u[1] + result[2]*u[2];
        el[i].uvs[0][1] = result[0]*v[0] + result[1]*v[1] + result[2]*v[2];
        result[0]       = result[3]*u[0] + result[4]*u[1] + result[5]*u[2];
        result[1]       = result[3]*v[0] + result[4]*v[1] + result[5]*v[2];
      }
      el[i].dir[0][0] = -result[0];
      el[i].dir[0][1] = -result[1];
      
    } else {
      
      el[i].nodes[0] = el[i].nodes[1] = nodes[0];
      if (n == 2) el[i].nodes[1] = nodes[1];
      if (pcurves != NULL) {
        stat = EG_evaluate(pcurves[i], &range[0], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Cannot evaluate PCurve %d = %d (EG_traceLoops)!\n",
                 i+1, stat);
          goto bail;
        }
        EG_adjustPeriod(result, srange, sper);
        el[i].uvs[0][0] = result[0];
        el[i].uvs[0][1] = result[1];
        result[0]       = result[2];
        result[1]       = result[3];
      } else {
        stat = EG_evaluate(edges[i], &range[0], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Cannot evaluate Edge %d = %d (EG_traceLoops)!\n",
                 i+1, stat);
          goto bail;
        }
        el[i].uvs[0][0] = result[0]*u[0] + result[1]*u[1] + result[2]*u[2];
        el[i].uvs[0][1] = result[0]*v[0] + result[1]*v[1] + result[2]*v[2];
        result[0]       = result[3]*u[0] + result[4]*u[1] + result[5]*u[2];
        result[1]       = result[3]*v[0] + result[4]*v[1] + result[5]*v[2];
      }
      el[i].dir[0][0] = result[0];
      el[i].dir[0][1] = result[1];
      
      if (pcurves != NULL) {
        stat = EG_evaluate(pcurves[i], &range[1], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Cannot evaluate PCurve %d = %d (EG_traceLoops)!\n",
                 i+1, stat);
          goto bail;
        }
        EG_adjustPeriod(result, srange, sper);
        el[i].uvs[1][0] = result[0];
        el[i].uvs[1][1] = result[1];
        result[0]       = result[2];
        result[1]       = result[3];
      } else {
        stat = EG_evaluate(edges[i], &range[1], result);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Cannot evaluate Edge %d = %d (EG_traceLoops)!\n",
                 i+1, stat);
          goto bail;
        }
        el[i].uvs[1][0] = result[0]*u[0] + result[1]*u[1] + result[2]*u[2];
        el[i].uvs[1][1] = result[0]*v[0] + result[1]*v[1] + result[2]*v[2];
        result[0]       = result[3]*u[0] + result[4]*u[1] + result[5]*u[2];
        result[1]       = result[3]*v[0] + result[4]*v[1] + result[5]*v[2];
      }
      el[i].dir[1][0] = result[0];
      el[i].dir[1][1] = result[1];
    }
  }
  
  /* figure out the uv tolerance for the existing Nodes */
  nl = (nodeCnt *) EG_alloc(2*nEdges*sizeof(nodeCnt));
  if (nl == NULL) {
    printf(" EGADS Error: Cannot Malloc %d Nodes (EG_traceLoops)!\n", 2*nEdges);
    stat = EGADS_MALLOC;
    goto bail;
  }
  tol = 1.e-7;
  for (nNode = i = 0; i < nEdges; i++) {
    for (j0 = 0; j0 < nNode; j0++)
      if (nl[j0].node == el[i].nodes[0]) break;
    if (j0 == nNode) {
      nl[nNode].node  = el[i].nodes[0];
      nl[nNode].index = nNode;
      nl[nNode].cnt   = 0;
      nNode++;
    }
    for (j1 = 0; j1 < nNode; j1++)
      if (nl[j1].node == el[i].nodes[1]) break;
    if (j1 == nNode) {
      nl[nNode].node  = el[i].nodes[1];
      nl[nNode].index = nNode;
      nl[nNode].cnt   = 0;
      nNode++;
    }
    if (el[i].edge->mtype == DEGENERATE) {
      stat = EG_tolerance(edges[i], &ang);
      if ((stat == EGADS_SUCCESS) && (ang > tol)) tol = ang;
      nl[j0].cnt = 1;
      nl[j1].cnt = 1;
    }
  }
#ifdef DEBUG
  printf(" start tol = %le\n", tol);
#endif
  duv[0] = duv[1] = duv[2] = duv[3] = 0.0;
  for (i = 0; i < nNode; i++) {
    if (nl[i].cnt == 1) continue;
    for (j0 = 0; j0 < nEdges; j0++)
      if (nl[i].node == el[j0].nodes[0]) {
        if (nl[i].cnt == 0) {
          duv[0] = duv[1] = el[j0].uvs[0][0];
          duv[2] = duv[3] = el[j0].uvs[0][1];
          nl[i].cnt--;
        } else {
          uvs[0] = el[j0].uvs[0][0];
          uvs[1] = el[j0].uvs[0][1];
          if (sper != 0) {
            uv[0] = 0.5*(duv[0] + duv[1]);
            uv[1] = 0.5*(duv[2] + duv[3]);
            EG_adjPeriodic(uv, uvs, srange, sper);
          }
          if (uvs[0] < duv[0]) duv[0] = uvs[0];
          if (uvs[0] > duv[1]) duv[1] = uvs[0];
          if (uvs[1] < duv[2]) duv[2] = uvs[1];
          if (uvs[1] > duv[3]) duv[3] = uvs[1];
        }
      }
    for (j0 = 0; j0 < nEdges; j0++)
      if (nl[i].node == el[j0].nodes[1]) {
        if (nl[i].cnt == 0) {
          duv[0] = duv[1] = el[j0].uvs[1][0];
          duv[2] = duv[3] = el[j0].uvs[1][1];
          nl[i].cnt--;
        } else {
          uvs[0] = el[j0].uvs[1][0];
          uvs[1] = el[j0].uvs[1][1];
          if (sper != 0) {
            uv[0] = 0.5*(duv[0] + duv[1]);
            uv[1] = 0.5*(duv[2] + duv[3]);
            EG_adjPeriodic(uv, uvs, srange, sper);
          }
          if (uvs[0] < duv[0]) duv[0] = uvs[0];
          if (uvs[0] > duv[1]) duv[1] = uvs[0];
          if (uvs[1] < duv[2]) duv[2] = uvs[1];
          if (uvs[1] > duv[3]) duv[3] = uvs[1];
        }
      }
    ang = sqrt((duv[1]-duv[0])*(duv[1]-duv[0]) +
               (duv[3]-duv[2])*(duv[3]-duv[2]));
    if (ang > tol) tol = 1.01*ang;
  }
#ifdef DEBUG
  printf(" adjst tol = %le\n", tol);
#endif
  EG_free(nl);
  
  tol /= 10.0;
  cnt  = 0;
  do {
    if (cnt == 3) {
      printf(" EGADS Error: Cannot match endPoints (EG_traceLoops)!\n");
      for (i = 0; i < nEdges; i++) {
        printf(" Edge %d: start uv = [%lf,%lf],  end uv = [%lf,%lf]\n",
               i+1, el[i].uvs[0][0], el[i].uvs[0][1],
                    el[i].uvs[1][0], el[i].uvs[1][1]);
      }
      stat = EGADS_INDEXERR;
      goto bail;
    }
    for (i = 0; i < nEdges; i++) el[i].connect[0] = el[i].connect[1] = -1;
    tol *= 10.0;
    
    for (n = i = 0; i < nEdges; i++) {
      if (el[i].connect[0] == -1) {
        el[i].connect[0] = n;
        stat = EG_matchEndpts(nEdges, el, n, el[i].uvs[0], tol);
        if (stat != EGADS_SUCCESS) break;
        n++;
      }
      if (el[i].connect[1] == -1) {
        el[i].connect[1] = n;
        stat = EG_matchEndpts(nEdges, el, n, el[i].uvs[1], tol);
        if (stat != EGADS_SUCCESS) break;
        n++;
      }
    }
    cnt++;
  } while (i != nEdges);
#ifdef DEBUG
  for (i = 0; i < nEdges; i++) {
    printf(" Edge %d: %d start uv = [%lf,%lf],  %d end uv = [%lf,%lf]\n",
           i, el[i].connect[0], el[i].uvs[0][0], el[i].uvs[0][1],
              el[i].connect[1], el[i].uvs[1][0], el[i].uvs[1][1]);
  }
#endif
  
  objs = (egObject **) EG_alloc(2*nEdges*sizeof(egObject *));
  if (objs == NULL) {
    stat = EGADS_MALLOC;
    goto bail;
  }
  sens = (int *) EG_alloc(nEdges*sizeof(int));
  if (sens == NULL) {
    stat = EGADS_MALLOC;
    goto bail;
  }
  
  /* connect the Edges to make Loops */
  n = nEdges;
  while (n != 0) {
    for (i = 0; i < nEdges; i++)
      if (el[i].edge != NULL) break;
    first        = el[i].connect[0];
    next         = el[i].connect[1];
    dir[0]       = el[i].dir[1][0];
    dir[1]       = el[i].dir[1][1];
    sens[0]      = el[i].sense;
    objs[0]      = el[i].edge;
    objs[nEdges] = el[i].pcurve;
    el[i].edge   = NULL;
    n            = 1;
#ifdef DEBUG
    printf("  first   = %d  next = %d  n = %d  type = %d\n",
           i, next, 0, objs[0]->mtype);
#endif
    if (first != next) {
      do {
        for (cnt = i = 0; i < nEdges; i++)
          if ((el[i].connect[0] == next) && (el[i].edge != NULL)) cnt++;
        if (cnt == 0) {
          printf(" EGADS Error: Open Loop (EG_traceLoops)!\n");
          stat = EGADS_DEGEN;
          goto bail;
        }
        if (cnt == 1) {
          for (i = 0; i < nEdges; i++)
            if ((el[i].connect[0] == next) && (el[i].edge != NULL)) break;
        } else {
          angle = 3.0*PI;
          for (cnt = i = 0; i < nEdges; i++)
            if ((el[i].connect[0] == next) && (el[i].edge != NULL)) {
              ang = EG_angleEdges(fsense, dir, el[i].dir[0]);
#ifdef DEBUG
              printf(" %d/%d: ang = %lf (%lf)  type = %d\n",
                     i, nEdges, ang, angle, el[i].edge->mtype);
              printf("       %lf %lf   %lf %lf\n", dir[0], dir[1],
                     el[i].dir[0][0], el[i].dir[0][1]);
#endif
              if ((ang < angle) && (fabs(ang) > 1.e-7)) {
                angle = ang;
                cnt   = i;
              }
            }
          if (angle > 2.0*PI) {
            printf(" EGADS Error: Cannot find Connection (EG_traceLoops)!\n");
            stat = EGADS_DEGEN;
            goto bail;
          }
          i = cnt;
        }
        next           = el[i].connect[1];
        dir[0]         = el[i].dir[1][0];
        dir[1]         = el[i].dir[1][1];
        sens[n]        = el[i].sense;
        objs[n]        = el[i].edge;
        objs[n+nEdges] = el[i].pcurve;
        el[i].edge     = NULL;
#ifdef DEBUG
        stat = -1;
        if (objs[n] != NULL) stat = objs[n]->mtype;
        printf("  current = %d  next = %d  n = %d  type = %d\n",
               i, next, n, stat);
#endif
        n++;
      } while (next != first);
    }
    
    /* make the loop */
    if (*nLoops == 0) {
      list = (egObject **) EG_alloc(sizeof(egObject *));
    } else {
      list = (egObject **) EG_reall(*loops, (*nLoops+1)*sizeof(egObject *));
    }
    if (list == NULL) {
      stat = EGADS_MALLOC;
      goto bail;
    }
    
#ifdef DEBUG
    printf(" EG_traceLoops: making %d Loop w/ %d of %d Edges!\n",
           *nLoops+1, n, nEdges);
#endif
    if (pcurves != NULL) {
      for (i = 0; i < n; i++) objs[i+n] = objs[i+nEdges];
      stat = EG_makeTopology(context, surf, LOOP, CLOSED, NULL, n, objs, sens,
                             &list[*nLoops]);
    } else {
      stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, n, objs, sens,
                             &list[*nLoops]);
    }
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: makeTopology %d = %d (EG_traceLoops)!\n",
             *nLoops+1, stat);
      goto bail;
    }
    if (list[*nLoops]->mtype == OPEN) {
      if (outLevel > 0)
        printf("             Loop %d is Open (EG_traceLoops)!\n", *nLoops+1);
      *nLoops += 1;
      *loops   = list;
      stat     = EGADS_TOPOERR;
      goto bail;
    }
    
    *nLoops += 1;
    *loops   = list;
    
    for (n = i = 0; i < nEdges; i++)
      if (el[i].edge != NULL) n++;
  }
  stat = EGADS_SUCCESS;
  
bail:
  if (stat != EGADS_SUCCESS) {
    list = *loops;
    if (list != NULL)
      for (i = 0; i < *nLoops; i++)
        EG_deleteObject(list[i]);
    if (*loops != NULL) EG_free(*loops);
    *loops  = NULL;
    *nLoops = 0;
  }
  if (sens != NULL) EG_free(sens);
  if (objs != NULL) EG_free(objs);
  if (el   != NULL) EG_free(el);
  return stat;
}


static void
EG_insertT(degenEdge *degenCnt, double t)
{
  int i;
  
  for (i = degenCnt->cnt-1; i > 0; i--) {
    degenCnt->ts[i+1] = degenCnt->ts[i];
    if ((t > degenCnt->ts[i-1]) && (t < degenCnt->ts[i])) {
      degenCnt->ts[i] = t;
      degenCnt->cnt  += 1;
      return;
    }
  }
  printf(" EGADS Internal: %lf NOT in range [%lf,%lf] (EG_insertT)!\n",
         t, degenCnt->ts[0], degenCnt->ts[degenCnt->cnt-1]);
}


static /*@null@*/ egObject *
nodeDegenerate(int nedges, egObject **edges, const egObject *node)
{
  int      i, n, status, oclass, mtype, *senses;
  double   data[4];
  egObject *ref, **nodes;
  
  for (i = 0; i < nedges; i++) {
    if (edges[i]->mtype != DEGENERATE) continue;
    status = EG_getTopology(edges[i], &ref, &oclass, &mtype, data, &n,
                            &nodes, &senses);
    if (status != EGADS_SUCCESS) {
      printf(" EGADS Warning: EG_getTopology = %d (nodeDegenerate)!\n", status);
      return NULL;
    }
    if (node == nodes[0]) return edges[i];
  }
  
  return NULL;
}


static int
EG_cutDegenerate(int iSeg, edgeInfo *fe, egObject *edge)
{
  int j;
  
  j = iSeg;
  while (j != -1) {
    if (fe[j].ni[0].degen == edge) return EGADS_SUCCESS;
    if (fe[j].ni[1].degen == edge) return EGADS_SUCCESS;
    j = fe[j].next;
  }
  
  return EGADS_OUTSIDE;
}


static int
EG_splitFace(egObject *context, int iFace, int iSeg, edgeInfo *fe,
             /*@unused@*/ splitEnts *sFaces, objStack *stack)
{
  int       i, j, k, m, n, stat, oclass, mtype, nLoop, nEdge, fsense, nFace, ol;
  int       nDegen, iper, *sens, *senses, *mark;
  double    t, tol, toler, area, oarea;
  double    uvbox[4], box[4], obox[4], uvs[6], d[2], trang[2], trange[2];
  egObject  *ref, *surf, *pcurve, *face, **objs, **edges, **pcurves, **loops;
  degenEdge *degenCnt;
  
  stat = EG_getTolerance(fe[iSeg].face, &tol);
  if (stat != EGADS_SUCCESS) return stat;
  stat = EG_getTopology(fe[iSeg].face, &surf, &oclass, &fsense, uvbox, &nLoop,
                        &loops, &sens);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Internal: Face %d -- EG_getTopology = %d\n",
           iFace+1, stat);
    return stat;
  }
#ifdef DEBUG
  for (i = 0; i < nLoop; i++) {
    stat = EG_getUVinfo(surf, loops[i], box, &area);
    if (stat == EGADS_SUCCESS)
      printf(" %d: mtype = %d  or = %d  UVinfo = %lf  %lf %lf  %lf %lf\n",
             i+1, surf->mtype, fsense, area, box[0], box[1], box[2], box[3]);
  }
#endif
  
  /* mark the degenerate NODE intersections on Degenerate Edges */
  for (nEdge = i = 0; i < nLoop; i++) {
    stat = EG_getTopology(loops[i], &ref, &oclass, &mtype, uvbox, &n,
                          &objs, &sens);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Internal: Face %d/Loop %d -- EG_getTopology = %d\n",
             iFace+1, i+1, stat);
      return stat;
    }
    nEdge += n;
  }
  edges = (egObject **) EG_alloc(nEdge*sizeof(egObject *));
  if (edges == NULL) {
    printf(" EGADS Internal: Face %d -- Malloc on %d Edge Objects\n",
           iFace+1, nEdge);
    return EGADS_MALLOC;
  }
  for (nEdge = i = 0; i < nLoop; i++) {
    EG_getTopology(loops[i], &ref, &oclass, &mtype, uvbox, &n, &objs, &sens);
    for (j = 0; j < n; j++, nEdge++) edges[nEdge] = objs[j];
  }
  j = iSeg;
  while (j != -1) {
    fe[j].ni[0].degen = fe[j].ni[1].degen = NULL;
    if (fe[j].ni[0].type == NODE)
      fe[j].ni[0].degen = nodeDegenerate(nEdge, edges, fe[j].ni[0].node);
#ifdef DEBUG
    if (fe[j].ni[0].degen != NULL)
      printf(" Input Edge %d:0 is Degenerate!\n", j+1);
#endif
    if (fe[j].ni[1].type == NODE)
      fe[j].ni[1].degen = nodeDegenerate(nEdge, edges, fe[j].ni[1].node);
#ifdef DEBUG
    if (fe[j].ni[1].degen != NULL)
      printf(" Input Edge %d:1 is Degenerate!\n", j+1);
#endif
    j = fe[j].next;
  }
  EG_free(edges);
  
  /* collect all of the Edges */
  for (nEdge = nDegen = i = 0; i < nLoop; i++) {
    stat = EG_getTopology(loops[i], &ref, &oclass, &mtype, uvbox, &n,
                          &objs, &sens);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Internal: Face %d/Loop %d -- EG_getTopology = %d\n",
             iFace+1, i+1, stat);
      return stat;
    }
    for (j = 0; j < n; j++)
      if (objs[j]->mtype == DEGENERATE) nDegen++;
    nEdge += n;
  }
#ifdef DEBUG
  printf(" EG_splitFace: # degenerate Edges = %d\n", nDegen);
#endif
  degenCnt = NULL;
  if (nDegen != 0) {
    degenCnt = (degenEdge *) EG_alloc(nDegen*sizeof(degenEdge));
    if (degenCnt == NULL) {
      printf(" EGADS Internal: Face %d -- Malloc Error on %d degen Edges\n",
             iFace+1, nDegen);
      return EGADS_MALLOC;
    }
    for (nDegen = i = 0; i < nLoop; i++) {
      EG_getTopology(loops[i], &ref, &oclass, &mtype, uvbox, &n, &objs, &sens);
      for (j = 0; j < n; j++) {
        if (objs[j]->mtype != DEGENERATE) continue;
        degenCnt[nDegen].edge   = objs[j];
        degenCnt[nDegen].pcurve = objs[n+j];
        degenCnt[nDegen].sense  = sens[j];
        degenCnt[nDegen].cnt    = 2;
        stat = EG_getRange(objs[j], degenCnt[nDegen].ts, &k);
        if (stat != EGADS_SUCCESS) {
          EG_free(degenCnt);
          printf(" EGADS Internal: Face %d EG_getRange Edge = %d\n",
                 iFace+1, stat);
          return stat;
        }
        nDegen++;
      }
    }
    for (j = 0; j < nDegen; j++) {
      stat = EG_getTopology(degenCnt[j].edge, &ref, &oclass, &mtype, d, &n,
                            &objs, &sens);
      if (stat != EGADS_SUCCESS) {
        EG_free(degenCnt);
        printf(" EGADS Internal: Face %d getTopology Edge = %d\n",
               iFace+1, stat);
        return stat;
      }
      degenCnt[j].node = objs[0];
    }
  }
  j = iSeg;
  m = n = 0;
  while (j != -1) {
    n++;
    if (fe[j].ni[0].degen != NULL) m++;
    if (fe[j].ni[1].degen != NULL) m++;
    j = fe[j].next;
  }
 
#ifdef DEBUG
  printf(" EG_splitFace: nEdges = %d from Face, %d from Splits, %d from Degen\n",
         nEdge, n, m);
#endif
  i = nEdge + 2*n + m;
  pcurves = NULL;
  edges   = (egObject **) EG_alloc(i*sizeof(egObject *));
  senses  = (int *)       EG_alloc(i*sizeof(int));
  if ((edges == NULL) || (senses == NULL)) {
    if (edges  != NULL) EG_free(edges);
    if (senses != NULL) EG_free(senses);
    EG_free(degenCnt);
    printf(" EGADS Internal: Face %d Cannot allocate %d Edges\n", iFace+1, i);
    return EGADS_MALLOC;
  }
  if (surf->mtype != PLANE) {
    pcurves = (egObject **) EG_alloc(i*sizeof(egObject *));
    if (pcurves == NULL) {
      EG_free(senses);
      EG_free(edges);
      EG_free(degenCnt);
      printf(" EGADS Internal: Face %d Cannot alloc %d PCurves\n", iFace+1, i);
      return EGADS_MALLOC;
    }
  }
  
  /* fill up the list of Edges */
  for (nEdge = i = 0; i < nLoop; i++) {
    stat = EG_getTopology(loops[i], &ref, &oclass, &mtype, uvbox, &n,
                          &objs, &sens);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Internal: Face %d/Loop %d -- EG_getTopology = %d\n",
             iFace+1, i+1, stat);
      if (pcurves != NULL) EG_free(pcurves);
      EG_free(senses);
      EG_free(edges);
      EG_free(degenCnt);
      return stat;
    }
    for (j = 0; j < n; j++) {
      if ((degenCnt != NULL) && (objs[j]->mtype == DEGENERATE)) {
        stat = EG_cutDegenerate(iSeg, fe, objs[j]);
        if (stat == EGADS_SUCCESS) {
#ifdef DEBUG
          printf(" Edge -: Loop %d  sense %2d  type %d -- %lx\n",
                 i+1, sens[j], objs[j]->mtype, (long) objs[j]);
#endif
          continue;
        }
      }
      edges[nEdge]  = objs[j];
      senses[nEdge] = sens[j];
      if (pcurves != NULL) pcurves[nEdge] = objs[j+n];
#ifdef DEBUG
      printf(" Edge %d: Loop %d  sense %2d  type %d\n",
             nEdge, i+1, sens[j], objs[j]->mtype);
#endif
      nEdge++;
    }
  }
  j = iSeg;
  while (j != -1) {
    pcurve          = fe[j].pcurve;
    edges[nEdge  ]  = fe[j].ent;
    senses[nEdge  ] = SFORWARD;
    edges[nEdge+1]  = fe[j].ent;
    senses[nEdge+1] = SREVERSE;
    if (pcurves != NULL) {
      pcurves[nEdge  ] = fe[j].pcurve;
      pcurves[nEdge+1] = fe[j].pcurve;
    }
    if ((pcurves != NULL) && (pcurve == NULL)) {
      toler = tol;
      ol    = EG_outLevel(context);
      EG_setOutLevel(context, 0);
      for (i = 0; i < 5; i++) {
        stat = EG_otherCurve(surf, fe[j].ent, toler, &pcurve);
        if (stat == EGADS_SUCCESS) break;
        if (stat != EGADS_CONSTERR) {
          printf(" EGADS Internal: Face %d -- EG_otherCurve = %d with tol = %le %d!\n",
                 iFace+1, stat, toler, i);
          if (pcurves != NULL) EG_free(pcurves);
          EG_free(senses);
          EG_free(edges);
          EG_free(degenCnt);
          return stat;
        }
#ifdef DEBUG
        printf(" EG_splitFace: %d otherCurve with tol = %le\n", i, tol);
#endif
        toler *= 10.0;
      }
      if (i == 5) {
        printf(" EGADS Internal: Face %d -- Cannot create PCurve!\n", iFace+1);
        if (pcurves != NULL) EG_free(pcurves);
        EG_free(senses);
        EG_free(edges);
        EG_free(degenCnt);
        return stat;
      }
      EG_setOutLevel(context, ol);
      if (pcurve != NULL) {
        stat = EG_stackPush(stack, pcurve);
        if (stat != EGADS_SUCCESS) {
          if (pcurves != NULL) EG_free(pcurves);
          EG_free(senses);
          EG_free(edges);
          EG_free(degenCnt);
          return stat;
        }
      }
      pcurves[nEdge  ] = pcurve;
      pcurves[nEdge+1] = pcurve;
      fe[j].pcurve     = pcurve;
    }
    nEdge += 2;
    j = fe[j].next;
  }

  /* patch up our split degenerate Edges */
  if (degenCnt != NULL) {
    j = iSeg;
    while (j != -1) {
      stat = EG_getRange(fe[j].ent, trang, &m);
      if (stat != EGADS_SUCCESS)
        printf(" EGADS Internal: Face %d -- %d getRange = %d\n",
               iFace+1, j, stat);
      if (fe[j].ni[0].degen != NULL) {
        for (k = 0; k < nDegen; k++)
          if (degenCnt[k].edge == fe[j].ni[0].degen) {
            if (degenCnt[k].cnt == 16) {
              printf(" EGADS Internal: Face %d -- > 15 Splits hit Degen Edge!\n",
                     iFace+1);
              if (pcurves != NULL) EG_free(pcurves);
              EG_free(senses);
              EG_free(edges);
              EG_free(degenCnt);
              return EGADS_DEGEN;
            }
            stat = EG_evaluate(fe[j].pcurve, &trang[0], uvs);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Internal: Face %d -- degen evaluate = %d!\n",
                     iFace+1, stat);
              if (pcurves != NULL) EG_free(pcurves);
              EG_free(senses);
              EG_free(edges);
              EG_free(degenCnt);
              return stat;
            }
            stat = EG_invEvaluate(degenCnt[k].pcurve, uvs, &t, d);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Internal: Face %d -- degen invEvaluate = %d!\n",
                     iFace+1, stat);
              if (pcurves != NULL) EG_free(pcurves);
              EG_free(senses);
              EG_free(edges);
              EG_free(degenCnt);
              return stat;
            }
            stat = EG_getRange(degenCnt[k].edge, trange, &iper);
            if (stat == EGADS_SUCCESS)
              if ((t <= trange[0]) || (t >= trange[1])) {
                printf(" EGADS Internal: Face %d -- t range = %lf [%lf, %lf]!\n",
                       iFace+1, t, trange[0], trange[1]);
                t = 0.5*(trange[0] + trange[1]);
              }
#ifdef DEBUG
            printf(" t = %lf (%lf), uvs = [%lf,%lf],  sense = %d\n",
                   t, trang[0], uvs[0], uvs[1], fe[j].sense);
#endif
            EG_insertT(&degenCnt[k], t);
            break;
          }
      }
      if (fe[j].ni[1].degen != NULL) {
        for (k = 0; k < nDegen; k++)
          if (degenCnt[k].edge == fe[j].ni[1].degen) {
            if (degenCnt[k].cnt == 16) {
              printf(" EGADS Internal: Face %d -- > 15 Splits hit Degen Edge!\n",
                     iFace+1);
              if (pcurves != NULL) EG_free(pcurves);
              EG_free(senses);
              EG_free(edges);
              EG_free(degenCnt);
              return EGADS_DEGEN;
            }
            stat = EG_evaluate(fe[j].pcurve, &trang[1], uvs);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Internal: Face %d -- degen evaluate = %d!\n",
                     iFace+1, stat);
              if (pcurves != NULL) EG_free(pcurves);
              EG_free(senses);
              EG_free(edges);
              EG_free(degenCnt);
              return stat;
            }
            stat = EG_invEvaluate(degenCnt[k].pcurve, uvs, &t, d);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Internal: Face %d -- degen invEvaluate = %d!\n",
                     iFace+1, stat);
              if (pcurves != NULL) EG_free(pcurves);
              EG_free(senses);
              EG_free(edges);
              EG_free(degenCnt);
              return stat;
            }
            stat = EG_getRange(degenCnt[k].edge, trange, &iper);
            if (stat == EGADS_SUCCESS)
              if ((t <= trange[0]) || (t >= trange[1])) {
                printf(" EGADS Internal: Face %d -- t range = %lf [%lf, %lf]!\n",
                       iFace+1, t, trange[0], trange[1]);
                t = 0.5*(trange[0] + trange[1]);
              }
#ifdef DEBUG
            printf(" t = %lf (%lf), uvs = [%lf,%lf],  sense = %d\n",
                   t, trang[1], uvs[0], uvs[1], fe[j].sense);
#endif
            EG_insertT(&degenCnt[k], t);
            break;
          }
      }
      j = fe[j].next;
    }
    
    /* make the new degenerate Edges based on t intervals of original */
    for (k = 0; k < nDegen; k++) {
      if (degenCnt[k].cnt == 2) continue;
/*    printf("    %lx  %d\n        ", (long) degenCnt[k].edge, degenCnt[k].cnt);
      for (j = 0; j < degenCnt[k].cnt; j++)
        printf(" %lf", degenCnt[k].ts[j]);
      printf("\n");  */
      for (j = 0; j < degenCnt[k].cnt-1; j++) {
        d[0] = degenCnt[k].ts[j];
        d[1] = degenCnt[k].ts[j+1];
        stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, d, 1,
                               &degenCnt[k].node, NULL, &ref);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Face %d -- EG_makeTopology = %d Degen Edge!\n",
                 iFace+1, stat);
          if (pcurves != NULL) EG_free(pcurves);
          EG_free(senses);
          EG_free(edges);
          EG_free(degenCnt);
          return stat;
        }
        stat = EG_stackPush(stack, ref);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Face %d -- EG_stackPush = %d Degen Edge!\n",
                 iFace+1, stat);
          if (pcurves != NULL) EG_free(pcurves);
          EG_free(senses);
          EG_free(edges);
          EG_free(degenCnt);
          return stat;
        }
        edges[nEdge]  = ref;
        senses[nEdge] = degenCnt[k].sense;
        if (pcurves != NULL) pcurves[nEdge] = degenCnt[k].pcurve;
        nEdge++;
      }
    }
  }
  if (degenCnt != NULL) EG_free(degenCnt);
  
  stat = EG_traceLoops(context, fe[iSeg].face, nEdge, edges, senses, pcurves,
                       &nLoop, &loops);
#ifdef DEBUG
  printf(" EGADS Info: traceLoops = %d, returns %d loops!\n", stat, nLoop);
#endif
  
  /* need to find inners and make the Faces */
  if (loops != NULL) {
    mark = (int *) EG_alloc(nLoop*sizeof(int));
    if (mark == NULL) {
      printf(" EGADS Error: Face %d -- Malloc of %d Loop Table!\n",
             iFace+1, nLoop);
      for (i = 0; i < nLoop; i++)
        EG_deleteObject(loops[i]);
      EG_free(loops);
      if (pcurves != NULL) EG_free(pcurves);
      EG_free(senses);
      EG_free(edges);
      return EGADS_MALLOC;
    }
    for (i = 0; i < nLoop; i++) mark[i] = 0;
    for (nFace = i = 0; i < nLoop; i++) {
      
      /* put the loops on the stack for cleanup later */
      stat = EG_stackPush(stack, loops[i]);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d -- EG_stackPush from Loop %d = %d!\n",
               iFace+1, stat, i+1);
        for (j = i; j < nLoop; j++) EG_deleteObject(loops[j]);
        EG_free(mark);
        EG_free(loops);
        if (pcurves != NULL) EG_free(pcurves);
        EG_free(senses);
        EG_free(edges);
        return stat;
      }
  
      stat = EG_getUVinfo(surf, loops[i], box, &area);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d -- EG_getUVinfo from Loop %d = %d!\n",
               iFace+1, i+1, stat);
        for (j = i+1; j < nLoop; j++) EG_deleteObject(loops[j]);
        EG_free(mark);
        EG_free(loops);
        if (pcurves != NULL) EG_free(pcurves);
        EG_free(senses);
        EG_free(edges);
        return stat;
      }
#ifdef DEBUG
      printf(" %d: mtype = %d  or = %d  UVinfo = %lf  %lf %lf  %lf %lf\n",
             i+1, surf->mtype, fsense, area, box[0], box[1], box[2], box[3]);
#endif
      if (area*fsense > 0.0) {
        nFace++;
        mark[i] = nFace;
      } else {
        for (j = 0; j < nLoop; j++) {
          if (i == j) continue;
          stat = EG_getUVinfo(surf, loops[j], obox, &oarea);
          if (stat != EGADS_SUCCESS) continue;
          if (oarea*fsense < 0.0) continue;
          if (fabs(area+oarea) <= tol) continue;
          if ((box[0] >= obox[0]) && (box[1] <= obox[1]) &&
              (box[2] >= obox[2]) && (box[3] <= obox[3])) {
            if (mark[i] != 0)
              printf(" EGADS Warning: Loop %d marked with %d\n", i+1, mark[i]);
            mark[i] = -j-1;
          }
        }
      }
    }
    for (j = i = 0; i < nLoop; i++) {
      if (mark[i] == 0) j++;
      if (mark[i] >  0) continue;
      k = -mark[i]-1;
      mark[i] = -mark[k];
    }
    if ((j != 0) || (nFace == 0)) {
      printf(" EGADS Internal: Face %d -- UnAssigned = %d  for nFace = %d\n",
             iFace+1, j, nFace);
      EG_free(mark);
      EG_free(loops);
      if (pcurves != NULL) EG_free(pcurves);
      EG_free(senses);
      EG_free(edges);
      return stat;
    }
    
    /* get the storage for the subFaces */
    stat = EG_splitAlloc(iFace, sFaces, nFace);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Internal: Face %d -- EG_splitAlloc = %d\n", iFace+1, stat);
      EG_free(mark);
      EG_free(loops);
      if (pcurves != NULL) EG_free(pcurves);
      EG_free(senses);
      EG_free(edges);
      return stat;
    }
    
    /* get the maximum number of Loops */
    for (m = i = 0; i < nFace; i++) {
      for (j = 0; j < nLoop; j++)
        if (mark[j] == i+1) break;
      k = 1;
      for (j = 0; j < nLoop; j++)
        if (mark[j] == -i-1) k++;
      if (k > m) m = k;
    }
/*  printf(" maxLoop = %d\n", m);  */
    objs = (egObject **) EG_alloc(m*sizeof(egObject *));
    sens = (int *)       EG_alloc(m*sizeof(int));
    if ((sens == NULL) || (objs == NULL)) {
      printf(" EGADS Error: Cannot allocate %d Loop storage\n", m);
      if (sens != NULL) EG_free(sens);
      if (objs != NULL) EG_free(objs);
      EG_free(mark);
      EG_free(loops);
      if (pcurves != NULL) EG_free(pcurves);
      EG_free(senses);
      EG_free(edges);
      return EGADS_MALLOC;
    }
    
    for (i = 0; i < nFace; i++) {
      for (j = 0; j < nLoop; j++)
        if (mark[j] == i+1) break;
      objs[0] = loops[j];
      sens[0] = SFORWARD;
      k = 1;
      for (j = 0; j < nLoop; j++)
        if (mark[j] == -i-1) {
          objs[k] = loops[j];
          sens[k] = SREVERSE;
          k++;
        }
      /* make the Face -- first try without paying attention to PCurves */
      ol   = EG_outLevel(context);
      EG_setOutLevel(context, 0);
      stat = EG_makeTopology(context, surf, FACE, PCURVE*fsense, NULL,
                             k, objs, sens, &face);
      EG_setOutLevel(context, ol);
      if (stat == EGADS_CONSTERR)
        stat = EG_makeTopology(context, surf, FACE, fsense, NULL,
                               k, objs, sens, &face);
#ifdef DEBUG
      EG_tolerance(face, &t);
      printf(" %d: newFace (%d Loops) status = %d  toler = %le\n",
             i+1, k, stat, t);
#endif
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d -- subFace %d makeTopology = %d\n",
               iFace+1, i+1, stat);
        EG_free(sens);
        EG_free(objs);
        EG_free(mark);
        EG_free(loops);
        if (pcurves != NULL) EG_free(pcurves);
        EG_free(senses);
        EG_free(edges);
        return stat;
      }
      stat = EG_stackPush(stack, face);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d -- subFace %d pushStack = %d\n",
               iFace+1, i+1, stat);
        EG_free(sens);
        EG_free(objs);
        EG_free(mark);
        EG_free(loops);
        if (pcurves != NULL) EG_free(pcurves);
        EG_free(senses);
        EG_free(edges);
        return stat;
      }
      EG_attributeDup(fe[iSeg].face, face);
      sFaces[iFace].ents[i] = face;
    }
    
    EG_free(sens);
    EG_free(objs);
    EG_free(mark);
    EG_free(loops);
  }

  if (pcurves != NULL) EG_free(pcurves);
  EG_free(senses);
  EG_free(edges);
  return stat;
}


int
EG_splitBody(const egObject *body, int nseg, const egObject **facEdg,
             egObject **result)
{
  int       i, j, jj, k, kk, l, m, n, *senses, *sen, *newSen, *eface = NULL;
  int       oclass, mtype, oc, mt, len, per, outLevel, status = EGADS_SUCCESS;
  int        nnodes,  nedges,  nfaces,  nshells;
  egObject  **nodes, **edges, **faces, **shells, **newObjs, **newFaces;
  egObject  *context, *ref, *lnodes[2], **children, **dum, **objs, *geom, *obj;
  egObject  *newBody;
  double    d, t, tol, toln, toler, period, xyz[3], data[4], limits[4], uv[2];
  double    *ts, trange[2];
  edgeInfo  *fe;
  objStack  stack;
  splitEnts *sEdges = NULL, *sFaces = NULL;
#ifdef WRITERESULT
  egObject    *model;
  static char fname[13] = "body_a.egads";
#endif
  
  *result = NULL;
  if (body->oclass != BODY)     return EGADS_NOTBODY;
  if (body->mtype  == WIREBODY) return EGADS_TOPOERR;
  
  context = EG_context(body);
  if (context == NULL) return EGADS_NOTCNTX;
  outLevel = EG_outLevel(body);
  
#ifdef WRITERESULT
  printf(" EG_splitBody: nSeg = %d\n", nseg);
  objs = (egObject  **) EG_alloc((nseg+1)*sizeof(egObject  *));
  if (objs == NULL) return EGADS_MALLOC;
  status = EG_copyObject(body, NULL, &objs[0]);
  if (status == EGADS_SUCCESS) {
    for (i = 0; i < nseg; i++) {
      if (facEdg[2*i+1]->oclass == EDGE) {
        jj     = 1;
        status = EG_makeTopology(context, NULL, LOOP, OPEN, NULL, 1,
                                 (egObject **) &facEdg[2*i+1], &jj, &obj);
        if (status != EGADS_SUCCESS) {
          for (j = 0; j <= i; j++) EG_deleteObject(objs[j]);
          EG_free(objs);
          printf(" EGADS Error: Cannot make Loop = %d!\n", status);
          return status;
        }
      } else {
        status = EG_copyObject(facEdg[2*i+1], NULL, &obj);
        if (status != EGADS_SUCCESS) {
          for (j = 0; j <= i; j++) EG_deleteObject(objs[j]);
          EG_free(objs);
          printf(" EGADS Error: Cannot copy Loop = %d!\n", status);
          return status;
        }
      }
      status = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &obj,
                               senses, &objs[i+1]);
      if (status != EGADS_SUCCESS) {
        for (j = 0; j <= i; j++) EG_deleteObject(objs[j]);
        EG_free(objs);
        printf(" EGADS Error: Cannot make Loop = %d!\n", status);
        return status;
      }
    }
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL, nseg+1, objs, NULL,
                             &model);
    if (status == EGADS_SUCCESS) {
      printf(" EGADS Info: writing start.egads = %d!\n",
             EG_saveModel(model, "start.egads"));
      EG_deleteObject(model);
    }
    EG_free(objs);
  }
#endif

  nnodes  = nedges = nfaces = nshells = 0;
  nodes   =  edges =  faces =  shells = NULL;
  fe      = (edgeInfo *) EG_alloc(nseg*sizeof(edgeInfo));
  if (fe == NULL) return EGADS_MALLOC;
  for (i = 0; i < nseg; i++) {
    fe[i].face           = (egObject *) facEdg[2*i  ];
    fe[i].ent            = (egObject *) facEdg[2*i+1];
    fe[i].pcurve         = NULL;
    fe[i].next           = -1;
    fe[i].sense          = fe[i].shell          = 1;
    fe[i].ni[0].entity   = fe[i].ni[0].node     = NULL;
    fe[i].ni[0].degen    = NULL;
    fe[i].ni[0].type     = EMPTY;
    fe[i].ni[0].match    = 0;
    fe[i].ni[0].param[0] = fe[i].ni[0].param[1] = 0.0;
    fe[i].ni[1]          = fe[i].ni[0];
  }
  status  = EG_stackInit(&stack);
  if (status != EGADS_SUCCESS) goto cleanup;
  
  /* get our PCurves if required */
  for (i = 0; i < nseg; i++) {
    status = EG_getTopology(fe[i].face, &ref, &oclass, &mtype, data, &n,
                            &children, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (ref->mtype == PLANE) continue;
    status = EG_otherCurve(fe[i].face, fe[i].ent, 0.0, &fe[i].pcurve);
    if ((status != EGADS_SUCCESS) || (fe[i].pcurve == NULL)) continue;
    status = EG_stackPush(&stack, fe[i].pcurve);
    if (status != EGADS_SUCCESS) goto cleanup;
  }

  /* get the data from the input body */

  status = EG_getBodyTopos(body, NULL, NODE,  &nnodes,  &nodes);
  if (status != EGADS_SUCCESS) goto cleanup;
  if (nodes == NULL) {
    status = EGADS_NULLOBJ;
    goto cleanup;
  }
  status = EG_getBodyTopos(body, NULL, EDGE,  &nedges,  &edges);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_getBodyTopos(body, NULL, FACE,  &nfaces,  &faces);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_getBodyTopos(body, NULL, SHELL, &nshells, &shells);
  if (status != EGADS_SUCCESS) goto cleanup;
  
  if (nshells > 1) {
    printf(" EGADS Internal: # Shells = %d\n", nshells);
    for (i = 0; i < nseg; i++) {
      status = EG_getBodyTopos(body, fe[i].face, SHELL, &n, &objs);
      if (status != EGADS_SUCCESS) goto cleanup;
      if (n != 1) printf(" EGADS Internal: Face is in %d Shells!\n", n);
      fe[i].shell = EG_indexBodyTopo(body, objs[0]);
      EG_free(objs);
    }
  }
  
  /* assign the input Edge/Loop endpoints to their respective entities */

  for (i = 0; i < nseg; i++) {
    status = EG_getTopology(fe[i].ent, &ref, &oclass, &mtype, data, &n,
                            &children, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (oclass == LOOP) {
      /* need to get endpoints of the loop */
      status = EG_getTopology(children[0], &ref, &oc, &mtype, data, &m, &dum,
                              &sen);
      if (status != EGADS_SUCCESS) goto cleanup;
      if (senses[0] == SFORWARD) {
        lnodes[0] = dum[0];
      } else {
        lnodes[0] = dum[1];
      }
      status = EG_getTopology(children[n-1], &ref, &oc, &mtype, data, &m, &dum,
                              &sen);
      if (status != EGADS_SUCCESS) goto cleanup;
      if (senses[0] == SFORWARD) {
        lnodes[1] = dum[1];
      } else {
        lnodes[1] = dum[0];
      }
      children = (egObject **) &lnodes;
      n = 1;
      if (lnodes[0] != lnodes[1]) n = 2;
    } else {
      if (mtype == DEGENERATE) {
        printf(" EGADS Info: Degenerate Edge found in Intersection!\n");
        continue;
      }
    }
    
    status = EG_getTolerance(fe[i].ent, &tol);
    if (status != EGADS_SUCCESS) goto cleanup;
/*  toler = tol;  */
    for (j = 0; j < n; j++) {
      fe[i].ni[j].node = children[j];
      status = EG_getTolerance(children[j], &toler);
      if (status != EGADS_SUCCESS) goto cleanup;
      if (toler < tol) toler = tol;
      status = EG_getTopology(children[j], &ref, &oc, &mtype, xyz, &m, &dum,
                              &sen);
      if (status != EGADS_SUCCESS) goto cleanup;
      
      /* check to see if we have been here */
      for (k = 0; k < i; k++) {
        status = EG_getTopology(fe[k].ni[0].node, &ref, &oc, &mtype, data, &m,
                                &dum, &sen);
        if (status != EGADS_SUCCESS) goto cleanup;
        status = EG_getTolerance(fe[k].ni[0].node, &toln);
        if (status != EGADS_SUCCESS) goto cleanup;
        if (toln < toler) toln = toler;
        d = sqrt((xyz[0]-data[0])*(xyz[0]-data[0]) +
                 (xyz[1]-data[1])*(xyz[1]-data[1]) +
                 (xyz[2]-data[2])*(xyz[2]-data[2]));
        if (d <= toln) {
#ifdef DEBUG
          printf(" Input Edge %d:%d is the same as Input Edge %d:0 (%le)!\n",
                 i+1, j, k+1, d);
#endif
          fe[i].ni[j].match = k+1;
          break;
        }
        status = EG_getTopology(fe[k].ni[1].node, &ref, &oc, &mtype, data, &m,
                                &dum, &sen);
        if (status != EGADS_SUCCESS) goto cleanup;
        status = EG_getTolerance(fe[k].ni[1].node, &toln);
        if (status != EGADS_SUCCESS) goto cleanup;
        if (toln < toler) toln = toler;
        d = sqrt((xyz[0]-data[0])*(xyz[0]-data[0]) +
                 (xyz[1]-data[1])*(xyz[1]-data[1]) +
                 (xyz[2]-data[2])*(xyz[2]-data[2]));
        if (d <= toln) {
#ifdef DEBUG
          printf(" Input Edge %d:%d is the same as Input Edge %d:1 (%le)!\n",
                 i+1, j, k+1, d);
#endif
          fe[i].ni[j].match = -(k+1);
          break;
        }
      }
      
      /* compare to nodes in the body */
      for (k = 0; k < nnodes; k++) {
        status = EG_getTopology(nodes[k], &ref, &oc, &mtype, data, &m, &dum,
                                &sen);
        if (status != EGADS_SUCCESS) goto cleanup;
        status = EG_getTolerance(nodes[k], &toln);
        if (status != EGADS_SUCCESS) goto cleanup;
        if (toln < toler) toln = toler;
        d = sqrt((xyz[0]-data[0])*(xyz[0]-data[0]) +
                 (xyz[1]-data[1])*(xyz[1]-data[1]) +
                 (xyz[2]-data[2])*(xyz[2]-data[2]));
        if (d <= toln) {
#ifdef DEBUG
          printf(" Input Edge %d:%d is the same as Node %d (%le -- %le)!\n",
                 i+1, j, k+1, d, toln);
#endif
          fe[i].ni[j].type = NODE;
          fe[i].ni[j].node = nodes[k];
          break;
        }
      }
      if (fe[i].ni[j].type == NODE) continue;
      
      /* look for splitting an Edge */
      status = EG_getBodyTopos(body, fe[i].face, EDGE, &m,  &dum);
      if (status != EGADS_SUCCESS) goto cleanup;
      for (k = 0; k < m; k++) {
        if (dum[k]->mtype == DEGENERATE) continue;
        status = EG_getTolerance(dum[k], &toler);
        if (status != EGADS_SUCCESS) {
          EG_free(dum);
          goto cleanup;
        }
        if (toler < tol) toler = tol;
        status = EG_invEvaluate(dum[k], xyz, &t, data);
        if (status != EGADS_SUCCESS) {
          EG_free(dum);
          goto cleanup;
        }
        d = sqrt((xyz[0]-data[0])*(xyz[0]-data[0]) +
                 (xyz[1]-data[1])*(xyz[1]-data[1]) +
                 (xyz[2]-data[2])*(xyz[2]-data[2]));
        if (d <= toler) {
          status = EG_getRange(dum[k], trange, &per);
          if (status != EGADS_SUCCESS) {
            EG_free(dum);
            goto cleanup;
          }
          if (per != 0) {
            status = EG_getTopology(dum[k], &ref, &oc, &mt, data, &l, &objs,
                                    &sen);
            if (status != EGADS_SUCCESS) {
              EG_free(dum);
              goto cleanup;
            }
            status = EG_getRange(ref, data, &per);
            if (status != EGADS_SUCCESS) {
              EG_free(dum);
              goto cleanup;
            }
            period = data[1] - data[0];
            while (t < trange[0]) t += period;
            while (t > trange[1]) t -= period;
          }
#ifdef DEBUG
          printf(" Input Edge %d:%d is on Edge %d @ %lf [%lf %lf] (%le)!\n",
                 i+1, j, EG_indexBodyTopo(body, dum[k]), t, trange[0], trange[1],
                 d);
          printf("            node = %lx   %lf %lf %lf\n",
                 (long) fe[i].ni[j].node, xyz[0], xyz[1], xyz[2]);
#endif
          fe[i].ni[j].type     = EDGE;
          fe[i].ni[j].entity   = dum[k];
          fe[i].ni[j].param[0] = t;
          break;
        }
      }
      EG_free(dum);
      if (fe[i].ni[j].type == EDGE) continue;
      
      /* must be interior in the Face */
      status = EG_invEvaluate(fe[i].face, xyz, uv, data);
      if (status != EGADS_SUCCESS) goto cleanup;
#ifdef DEBUG
      d = sqrt((xyz[0]-data[0])*(xyz[0]-data[0]) +
               (xyz[1]-data[1])*(xyz[1]-data[1]) +
               (xyz[2]-data[2])*(xyz[2]-data[2]));
      printf(" Input Edge %d:%d is in Face %d @ [%lf,%lf] (%le)!\n",
             i+1, j, EG_indexBodyTopo(body, fe[i].face), uv[0], uv[1], d);
      printf("            node = %lx   %lf %lf %lf\n",
             (long) fe[i].ni[j].node, xyz[0], xyz[1], xyz[2]);
#endif
      /* do we need to adjust for periodics? */
      fe[i].ni[j].type     = FACE;
      fe[i].ni[j].entity   = fe[i].face;
      fe[i].ni[j].param[0] = uv[0];
      fe[i].ni[j].param[1] = uv[1];
    }
    /* copy node for closed entity */
    if (n == 1) fe[i].ni[1] = fe[i].ni[0];
  }
  
  /* check for consistency */
  
  for (i = 0; i < nseg; i++)
    for (j = 0; j < 2; j++) {
      if (fe[i].ni[j].match == 0) continue;
      k = abs(fe[i].ni[j].match)-1;
      l = 0;
      if  (fe[i].ni[j].match < 0) l = 1;
      if ((fe[i].ni[j].type   == fe[k].ni[l].type) &&
          (fe[i].ni[j].entity == fe[k].ni[l].entity)) {
        fe[i].ni[j] = fe[k].ni[l];    /* make sure only 1 set of parameters */
        continue;
      }
      if (fe[i].ni[j].type > fe[k].ni[l].type) {
        if (fe[i].shell != 0) {
#ifdef WIN32
          printf(" EGADS Internal: Bad Topo matching -- %llx (%d)  %llx (%d)!\n",
                 (long long) fe[i].ni[j].entity, fe[i].ni[j].type,
                 (long long) fe[k].ni[l].entity, fe[k].ni[l].type);
          printf("                 Disabling Segment %d!\n", i+1);
#else
          printf(" EGADS Internal: Bad Topo matching -- %lx (%d)  %lx (%d)!\n",
                 (long) fe[i].ni[j].entity, fe[i].ni[j].type,
                 (long) fe[k].ni[l].entity, fe[k].ni[l].type);
          printf("                 Disabling Segment %d!\n", i+1);
#endif
        }
        fe[i].shell = 0;
        continue;
      }
      if (fe[i].ni[j].type < fe[k].ni[l].type) {
        if (fe[k].shell != 0) {
#ifdef WIN32
          printf(" EGADS Internal: Bad Topo matching -- %llx (%d)  %llx (%d)!\n",
                 (long long) fe[i].ni[j].entity, fe[i].ni[j].type,
                 (long long) fe[k].ni[l].entity, fe[k].ni[l].type);
#else
          printf(" EGADS Internal: Bad Topo matching -- %lx (%d)  %lx (%d)!\n",
                 (long) fe[i].ni[j].entity, fe[i].ni[j].type,
                 (long) fe[k].ni[l].entity, fe[k].ni[l].type);
#endif
          printf("                 Disabling segment %d!\n", k+1);
        }
        fe[k].shell = 0;
        fe[k].ni[l] = fe[i].ni[j];
        continue;
      }
      if (fe[i].ni[j].entity != fe[k].ni[l].entity) {
#ifdef WIN32
        printf(" EGADS Warning: Bad Topo matching -- %llx (%d)  %llx (%d)!\n",
               (long long) fe[i].ni[j].entity, fe[i].ni[j].type,
               (long long) fe[k].ni[l].entity, fe[k].ni[l].type);
#else
        printf(" EGADS Warning: Bad Topo matching -- %lx (%d)  %lx (%d)!\n",
               (long) fe[i].ni[j].entity, fe[i].ni[j].type,
               (long) fe[k].ni[l].entity, fe[k].ni[l].type);
#endif
      }
    }
  
  /* get Edges/Loops that cut each Face */
  
  eface = (int *) EG_alloc(nfaces*sizeof(int));
  if (eface == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }
  k = -1;
  for (i = 0; i < nfaces; i++) eface[i] = -1;
  for (i = 0; i < nfaces; i++)
    for (m = j = 0; j < nseg; j++) {
      if (fe[j].shell == 0) continue;
      if (fe[j].face == faces[i]) {
        if (m == 0) {
          eface[i] = j;
          k        = j;
        } else {
          fe[k].next = j;
          k          = j;
        }
        m++;
      }
    }

  /* check internal Edges for hanging segments */
  status = EG_delOpenEdges(nseg, nfaces, eface, fe);
  if (status != EGADS_SUCCESS) goto cleanup;
  
  /* remake all of the input Edges -- so we know they have the correct Nodes! */
  for (i = 0; i < nseg; i++) {
    if (fe[i].shell == 0) continue;
    status = EG_getTopology(fe[i].ent, &ref, &oclass, &mtype, data, &n,
                            &children, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    if (oclass == LOOP) {
      if (n == 1) {
        status = EG_getTopology(children[0], &geom, &oc, &mt, limits, &l,
                                &dum, &sen);
        if (status != EGADS_SUCCESS) goto cleanup;
        lnodes[0] = dum[0];
        if (l != 1) lnodes[1] = dum[1];
        if (senses[0] == SFORWARD) {
          lnodes[0] = fe[i].ni[0].node;
          lnodes[1] = fe[i].ni[1].node;
        } else {
          lnodes[1] = fe[i].ni[0].node;
          lnodes[0] = fe[i].ni[1].node;
        }
        status = EG_makeTopology(context, geom, oc, mt, limits, l, lnodes,
                                 sen, &obj);
        if (status != EGADS_SUCCESS) goto cleanup;
        status = EG_stackPush(&stack, obj);
        if (status != EGADS_SUCCESS) goto cleanup;
        
        status = EG_makeTopology(context, ref, oclass, mtype, data, n, &obj,
                                 senses, &fe[i].ent);
        if (status != EGADS_SUCCESS) goto cleanup;
        status = EG_stackPush(&stack, fe[i].ent);
        if (status != EGADS_SUCCESS) goto cleanup;
      } else {
        objs = (egObject **) EG_alloc(n*sizeof(egObject *));
        if (objs == NULL) {
          status = EGADS_MALLOC;
          goto cleanup;
        }
        for (j = 1; j < n-1; j++) objs[j] = children[j];
        
        status = EG_getTopology(children[0], &geom, &oc, &mt, limits, &l,
                                &dum, &sen);
        if (status != EGADS_SUCCESS) {
          EG_free(objs);
          goto cleanup;
        }
        lnodes[0] = dum[0];
        if (l != 1) lnodes[1] = dum[1];
        if (senses[0] == SFORWARD) {
          lnodes[0] = fe[i].ni[0].node;
        } else {
          lnodes[1] = fe[i].ni[0].node;
        }
        status = EG_makeTopology(context, geom, oc, mt, limits, l, lnodes,
                                 sen, &objs[0]);
        if (status != EGADS_SUCCESS) {
          EG_free(objs);
          goto cleanup;
        }
        status = EG_stackPush(&stack, objs[0]);
        if (status != EGADS_SUCCESS) {
          EG_free(objs);
          goto cleanup;
        }
        
        status = EG_getTopology(children[n-1], &geom, &oc, &mt, limits, &l,
                                &dum, &sen);
        if (status != EGADS_SUCCESS) {
          EG_free(objs);
          goto cleanup;
        }
        lnodes[0] = dum[0];
        if (l != 1) lnodes[1] = dum[1];
        if (senses[0] == SFORWARD) {
          lnodes[1] = fe[i].ni[1].node;
        } else {
          lnodes[0] = fe[i].ni[1].node;
        }
        status = EG_makeTopology(context, geom, oc, mt, limits, l, lnodes,
                                 sen, &objs[n-1]);
        if (status != EGADS_SUCCESS) {
          EG_free(objs);
          goto cleanup;
        }
        status = EG_stackPush(&stack, objs[n-1]);
        if (status != EGADS_SUCCESS) {
          EG_free(objs);
          goto cleanup;
        }
        
        status = EG_makeTopology(context, ref, oclass, mtype, data, n, objs,
                                 senses, &fe[i].ent);
        EG_free(objs);
        if (status != EGADS_SUCCESS) goto cleanup;
        status = EG_stackPush(&stack, fe[i].ent);
        if (status != EGADS_SUCCESS) goto cleanup;
      }
    } else {
      lnodes[0] = fe[i].ni[0].node;
      lnodes[1] = fe[i].ni[1].node;
      status = EG_makeTopology(context, ref, oclass, mtype, data, n, lnodes,
                               senses, &fe[i].ent);
      if (status != EGADS_SUCCESS) goto cleanup;
      status = EG_stackPush(&stack, fe[i].ent);
      if (status != EGADS_SUCCESS) goto cleanup;
    }
  }
  
  /* make the split edges -- ordered by t! */
  EG_splitInit(nedges, &sEdges);
  if (sEdges == NULL) goto cleanup;
  dum = (egObject **) EG_alloc((2*nseg+1)*sizeof(egObject *));
  if (dum == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }
  ts = (double *) EG_alloc((2*nseg+1)*sizeof(double));
  if (ts == NULL) {
    EG_free(dum);
    status = EGADS_MALLOC;
    goto cleanup;
  }
  for (i = 0; i < nedges; i++) {
    for (n = j = 0; j < nseg; j++) {
      if (fe[j].shell == 0) continue;
      if (fe[j].ni[0].type == EDGE)
        if (edges[i] == fe[j].ni[0].entity) {
          for (m = l = 0; l < n; l++)
            if (dum[l] == fe[j].ni[0].node) m++;
          if (m == 0) {
            dum[n] = fe[j].ni[0].node;
            ts[n]  = fe[j].ni[0].param[0];
            n++;
          }
        }
      if (fe[j].ni[1].type == EDGE)
        if (edges[i] == fe[j].ni[1].entity) {
          for (m = l = 0; l < n; l++)
            if (dum[l] == fe[j].ni[1].node) m++;
          if (m == 0) {
            dum[n] = fe[j].ni[1].node;
            ts[n]  = fe[j].ni[1].param[0];
            n++;
          }
        }
    }
    if (n > nseg) {
      printf(" EGADS Internal: n = %d nseg = %d!\n", n, nseg);
      EG_free(ts);
      EG_free(dum);
      status = EGADS_INDEXERR;
      goto cleanup;
    }
    if (n == 0) continue;
    do {
      for (m = j = 0; j < n-1; j++)
        if (ts[j+1] < ts[j]) {
          t        = ts[j+1];
          ts[j+1]  = ts[j];
          ts[j]    = t;
          ref      = dum[j+1];
          dum[j+1] = dum[j];
          dum[j]   = ref;
          m++;
        }
    } while (m != 0);
    status = EG_splitAlloc(i, sEdges, n+1);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_getTopology(edges[i], &ref, &oclass, &mtype, data, &l,
                            &children, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
    limits[1] = data[0];
    lnodes[1] = children[0];
    for (j = 0; j < n; j++) {
#ifdef DEBUG
      printf(" Edge %3d: split @ %lf, Node = %lx\n", i+1, ts[j], (long) dum[j]);
#endif
      limits[0] = limits[1];
      lnodes[0] = lnodes[1];
      limits[1] = ts[j];
      lnodes[1] = dum[j];
      status    = EG_makeTopology(context, ref, oclass, mtype, limits, 2,
                                  lnodes, senses, &sEdges[i].ents[j]);
#ifdef DEBUG
      printf("           newEdge status = %d\n", status);
#endif
      if (status != EGADS_SUCCESS) {
        EG_free(ts);
        EG_free(dum);
        goto cleanup;
      }
      status = EG_stackPush(&stack, sEdges[i].ents[j]);
      if (status != EGADS_SUCCESS) {
        EG_free(ts);
        EG_free(dum);
        goto cleanup;
      }
    }
    limits[0] = limits[1];
    lnodes[0] = lnodes[1];
    limits[1] = data[1];
    lnodes[1] = children[l-1];
    status    = EG_makeTopology(context, ref, oclass, mtype, limits, 2,
                                lnodes, senses, &sEdges[i].ents[n]);
#ifdef DEBUG
    printf("           newEdge status = %d\n", status);
#endif
    if (status != EGADS_SUCCESS) {
      EG_free(ts);
      EG_free(dum);
      goto cleanup;
    }
    status = EG_stackPush(&stack, sEdges[i].ents[n]);
    if (status != EGADS_SUCCESS) {
      EG_free(ts);
      EG_free(dum);
      goto cleanup;
    }
  }
  EG_free(ts);
  EG_free(dum);
  
  /* rebuild the Faces with our split Edges */
  
  newFaces = (egObject **) EG_alloc(nfaces*sizeof(egObject *));
  if (newFaces == NULL) {
    status = EGADS_MALLOC;
    goto cleanup;
  }
  for (i = 0; i < nfaces; i++) {
    newFaces[i] = NULL;
    status = EG_getBodyTopos(body, faces[i], EDGE, &m,  &dum);
    if (status != EGADS_SUCCESS) {
      EG_free(newFaces);
      goto cleanup;
    }
    for (l = j = 0; j < m; j++) {
      k = EG_indexBodyTopo(body, dum[j])-1;
      if (k < 0) continue;
      if (sEdges[k].nSplit != 0) l++;
    }
    EG_free(dum);
    if (l == 0) continue;
    
    status = EG_getTopology(faces[i], &geom, &oclass, &mtype, data, &n,
                            &children, &senses);
    if (status != EGADS_SUCCESS) goto cleanup;
#ifdef DEBUG
    printf(" Face %d: %d Loop(s) -- number of Edge Splits is %d\n", i+1, n, l);
#endif
    dum = (egObject **) EG_alloc(n*sizeof(egObject *));
    if (dum == NULL) {
      EG_free(newFaces);
      status = EGADS_MALLOC;
      goto cleanup;
    }
    /* adjust each loop */
    for (j = 0; j < n; j++) {
      status = EG_getTopology(children[j], &ref, &oc, &mt, limits, &m,
                              &objs, &sen);
      if (status != EGADS_SUCCESS) goto cleanup;
      for (l = k = 0; k < m; k++) {
        kk = EG_indexBodyTopo(body, objs[k])-1;
        if (kk < 0) continue;
        if (sEdges[kk].nSplit != 0) l += sEdges[kk].nSplit-1;
      }
#ifdef DEBUG
      printf("       Loop %d: ref = %lx  nEdge = %d (%d)\n",
             j+1, (long) ref, m, l);
#endif
      if (l == 0) {
        dum[j] = children[j];
        continue;
      }
      /* make new loop */
      len    = m+l;
      newSen = (int *) EG_alloc(len*sizeof(int));
      if (newSen == NULL) {
        EG_free(dum);
        EG_free(newFaces);
        status = EGADS_MALLOC;
        goto cleanup;
      }
      if (ref == NULL) {
        newObjs = (egObject **) EG_alloc(  len*sizeof(egObject *));
      } else {
        newObjs = (egObject **) EG_alloc(2*len*sizeof(egObject *));
      }
      if (newObjs == NULL) {
        EG_free(newSen);
        EG_free(dum);
        EG_free(newFaces);
        status = EGADS_MALLOC;
        goto cleanup;
      }
      for (l = k = 0; k < m; k++) {
        newObjs[l] = objs[k];
        newSen[l]  = sen[k];
        if (ref != NULL) newObjs[l+len] = objs[m+k];
        l++;
        kk = EG_indexBodyTopo(body, objs[k])-1;
        if (kk < 0) continue;
        if (sEdges[kk].nSplit == 0) continue;
        l--;
        for (jj = 0; jj < sEdges[kk].nSplit; jj++, l++) {
          newSen[l]  = sen[k];
          if (ref != NULL) newObjs[l+len] = objs[m+k];
          if (sen[k] == SFORWARD) {
            newObjs[l] = sEdges[kk].ents[jj];
          } else {
            newObjs[l] = sEdges[kk].ents[sEdges[kk].nSplit-jj-1];
          }
        }
      }
      status = EG_makeTopology(context, ref, oc, mt, NULL, len, newObjs, newSen,
                               &dum[j]);
#ifdef DEBUG
      printf("       newLoop status = %d for %d Edges\n", status, len);
#endif
      EG_free(newObjs);
      EG_free(newSen);
      if (status != EGADS_SUCCESS) {
        EG_free(dum);
        EG_free(newFaces);
        goto cleanup;
      }
      status = EG_stackPush(&stack, dum[j]);
      if (status != EGADS_SUCCESS) {
        EG_free(dum);
        EG_free(newFaces);
        goto cleanup;
      }
    }
    /* make new Faces */
    status = EG_makeTopology(context, geom, oclass, mtype, NULL, n, dum, senses,
                             &obj);
#ifdef DEBUG
    printf("       newFace status = %d for %d Loops\n", status, n);
#endif
    EG_free(dum);
    if (status != EGADS_SUCCESS) {
      EG_free(newFaces);
      goto cleanup;
    }
    EG_attributeDup(faces[i], obj);
    status = EG_stackPush(&stack, obj);
    if (status != EGADS_SUCCESS) {
      EG_free(newFaces);
      goto cleanup;
    }
    /* patch up the segment info with the new Face */
    j = eface[i];
    while (j != -1) {
      fe[j].face = obj;
      j          = fe[j].next;
    }
    newFaces[i]  = obj;
  }

  /* adjust shells */
  for (j = 0; j < nshells; j++) {
    status = EG_getTopology(shells[j], &ref, &oc, &mt, limits, &m, &objs,
                            &sen);
    if (status != EGADS_SUCCESS) goto cleanup;
    for (l = k = 0; k < m; k++) {
      i = EG_indexBodyTopo(body, objs[k])-1;
      if (i < 0) continue;
      if (newFaces[i] != NULL) l++;
    }
    if (l == 0) continue;
    
    dum = (egObject **) EG_alloc(m*sizeof(egObject *));
    if (dum == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    for (k = 0; k < m; k++) {
      dum[k] = objs[k];
      i = EG_indexBodyTopo(body, objs[k])-1;
      if (i < 0) continue;
      if (newFaces[i] != NULL) dum[k] = newFaces[i];
    }
    /* overwrite shell with the new one */
    status = EG_makeTopology(context, ref, oc, mt, NULL, m, dum, sen,
                             &shells[j]);
    EG_free(dum);
#ifdef DEBUG
    printf(" newShell %d: status = %d for %d Faces -- body = %d\n",
           j, status, m, body->mtype);
#endif
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(&stack, shells[j]);
    if (status != EGADS_SUCCESS) goto cleanup;
  }
  /* update our face storage */
  for (i = 0; i < nfaces; i++) {
    if (newFaces[i] != NULL) faces[i] = newFaces[i];
  }
  EG_free(newFaces);

  /* adjust the Loops & split Faces in each Face cut */
  
  EG_splitInit(nfaces, &sFaces);
  if (sFaces == NULL) goto cleanup;

  for (i = 0; i < nfaces; i++) {
    j = eface[i];
    if (j == -1) continue;
    status = EG_splitFace(context, i, j, fe, sFaces, &stack);
#ifdef DEBUG
    printf(" EG_splitFace = %d\n", status);
#endif
    if (status != EGADS_SUCCESS) goto cleanup;
  }
  
  /* make the new body by rebuilding Shells */
  
  newBody = (egObject *) body;
  status  = EG_getTopology(newBody, &geom, &oclass, &mtype, data, &n,
                           &children, &senses);
  if (status != EGADS_SUCCESS) goto cleanup;

  /* promote FaceBody with more than a single result to SheetBody */
  if (mtype == FACEBODY) {
    if (nfaces != 1) {
      printf(" EGADS Internal: FaceBody with nFace = %d!\n", nfaces);
      status = EGADS_TOPOERR;
      goto cleanup;
    }
    if (sFaces[0].nSplit == 0) {
      /* still a single Face -- finish up */
      status = EG_makeTopology(context, NULL, BODY, mtype, NULL, 1, faces,
                               senses, &obj);
#ifdef DEBUG
      printf(" newBody: status = %d \n", status);
#endif
      if (status != EGADS_SUCCESS) goto cleanup;
      
#ifdef WRITERESULT
      status = EG_copyObject(obj, NULL, &ref);
      if (status == EGADS_SUCCESS) {
        status = EG_makeTopology(context, NULL, MODEL, 0, NULL, 1, &ref, senses,
                                 &model);
        if (status == EGADS_SUCCESS) {
          printf(" EGADS Info: writing %s = %d!\n",
                 fname, EG_saveModel(model, fname));
          fname[5]++;
          EG_deleteObject(model);
        }
      }
#endif
      *result = obj;
      goto cleanup;
    }
    
    nshells = 1;
    shells  = (egObject **) EG_alloc(nshells*sizeof(egObject *));
    if (shells == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    status = EG_makeTopology(context, NULL, SHELL, 0, NULL, 1, children, senses,
                             &shells[0]);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(&stack, shells[0]);
    if (status != EGADS_SUCCESS) goto cleanup;
    mtype  = SHEETBODY;
    status = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL, 1, shells,
                             senses, &newBody);
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(&stack, newBody);
    if (status != EGADS_SUCCESS) goto cleanup;
    lnodes[0] = shells[0];
    children  = (egObject **) &lnodes;
  }  else {
    if (n != nshells)
      printf(" EGADS Internal: nShells = %d %d (EG_splitBody)!\n", nshells, n);
  }
  
  /* rebuild shells */
  
  for (j = 0; j < nshells; j++) {
    status = EG_getTopology(children[j], &ref, &oc, &mt, limits, &m, &objs,
                            &sen);
    if (status != EGADS_SUCCESS) goto cleanup;
    for (l = kk = k = 0; k < m; k++) {
      i = EG_indexBodyTopo(body, objs[k])-1;
      if (i < 0) continue;
      l += sFaces[i].nSplit;
      if (sFaces[i].nSplit != 0) kk++;
    }
    if (l == 0) continue;
    
    /* new Faces in the Shell */
    l += m - kk;
    dum = (egObject **) EG_alloc(l*sizeof(egObject *));
    if (dum == NULL) {
      status = EGADS_MALLOC;
      goto cleanup;
    }
    for (l = k = 0; k < m; k++) {
      i = EG_indexBodyTopo(body, objs[k])-1;
      if (i < 0) continue;
      dum[l] = faces[i];
      l++;
      if (sFaces[i].nSplit == 0) continue;
      l--;
      for (kk = 0; kk < sFaces[i].nSplit; kk++, l++)
        dum[l] = sFaces[i].ents[kk];
    }
    /* overwrite shell with the new one */
    status = EG_makeTopology(context, ref, oc, mt, NULL, l, dum, NULL,
                             &shells[j]);
    EG_free(dum);
#ifdef DEBUG
    printf(" newShell %d: status = %d for %d->%d Faces -- body = %d\n",
           j, status, m, l, mtype);
#endif
    if (status != EGADS_SUCCESS) goto cleanup;
    status = EG_stackPush(&stack, shells[j]);
    if (status != EGADS_SUCCESS) goto cleanup;
  }
  
  /* make the body from our new Shells */

  status = EG_makeTopology(context, NULL, BODY, mtype, NULL, nshells,
                           shells, senses, &obj);
#ifdef DEBUG
  printf(" newBody: status = %d \n", status);
#endif
  if (status != EGADS_SUCCESS) goto cleanup;
  
#ifdef WRITERESULT
  status = EG_copyObject(obj, NULL, &ref);
  if (status == EGADS_SUCCESS) {
    status = EG_makeTopology(context, NULL, MODEL, 0, NULL, 1, &ref, senses,
                             &model);
    if (status == EGADS_SUCCESS) {
      printf(" EGADS Info: writing %s = %d!\n",
             fname, EG_saveModel(model, fname));
      fname[5]++;
      EG_deleteObject(model);
    }
  }
#endif
  *result = obj;
  
cleanup:
  EG_splitFree(nedges, &sEdges);
  EG_splitFree(nfaces, &sFaces);
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (EG_splitBody)!\n", i);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);
  if (eface  != NULL) EG_free(eface);
  if (shells != NULL) EG_free(shells);
  if (faces  != NULL) EG_free(faces);
  if (edges  != NULL) EG_free(edges);
  if (nodes  != NULL) EG_free(nodes);
  if (fe     != NULL) EG_free(fe);
  
  return status;
}
