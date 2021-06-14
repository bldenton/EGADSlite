/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Effective (Virtual) Topology Functions
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
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
#include "egadsTris.h"
#include "emp.h"

//#define DEBUG


#define NOTFILLED       -1      /* Not yet filled flag */
#define PI              3.1415926535897931159979635
#define ABS(x)           ((x) < 0 ? -(x) : (x))
#define CROSS(a,b,c)       a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                           a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                           a[2] = (b[0]*c[1]) - (b[1]*c[0])

typedef struct {
  egObject *eedge;              /* the EEdge */
  int      sense;               /* the sense */
  int      iloop;               /* the ELoop index */
  int      iedge;               /* EEdge index in the ELoop */
  int      newel;               /* new ELoop index */
  int      prev;                /* the previous fragment index */
  int      bnode;               /* open node index at the beginning */
  int      enode;               /* open node index at the end */
  int      next;                /* the next fragment index */
} eloopFrag;


extern int  EG_uvmapLocate( void *uvmap, int *trmap, double *uv, int *fID,
                            int *itri, int *verts, double *ws );
extern int  EG_uv2UVmap( void *uvmap, int *trmap, double *fuv, double *fuvs,
                         int *tris, int tbeg, int tend, double *uv );
extern int  EG_uvmapMake( int ntri, int *tris, int *itris, int nvrt,
                          double *vrts, double *range, int **trmap,
                          void **uvmap );
extern void EG_getUVmap( void *uvmap, int index, double *uv );
extern int  EG_uvmapWrite( void *uvmap, int *trmap, FILE *fp );
extern int  EG_uvmapRead( FILE *fp, double *range, void **uvmap, int **trmap );
extern int  EG_uvmapCopy( void *uvsrc, int *trsrc, void **uvmap, int **trmap );
extern int  uvmap_struct_free( void *uvmap );

extern int  EG_fillArea( int ncontours, const int *cntr, const double *vertices,
                         int *triangles, int *n_fig8, int pass, fillArea *fa );
extern void EG_makeConnect(int k1, int k2, int *tri, int *kedge, int *ntable,
                           connect *etable, int face);
extern void EG_readAttrs( egObject *obj, int nattr, FILE *fp );
extern void EG_writeAttr( egAttrs *attrs, FILE *fp );
extern int  EG_writeNumAttr( egAttrs *attrs );
extern int  EG_inTriExact( double *t1, double *t2, double *t3, double *p,
                           double *w );
extern int  EG_computeTessMap( egTessel *btess, int outLevel );
#ifdef DEBUG
extern int  EG_evaluate( const egObject *geom, /*@null@*/ const double *param,
                         double *result );
#endif
extern int  EG_evaluatX( const egObject *geom, /*@null@*/ const double *param,
                         double *result );
extern int  EG_invEvaluatX( const egObject *geom, double *xyz, double *param,
                            double *result );
extern int  EG_getEdgeUVX( const egObject *face, const egObject *edge,
                           int sense, double t, double *result );
extern int  EG_getTopology( const egObject *topo, egObject **geom, int *oclass,
                            int *mtype, /*@null@*/ double *limits, int *nChild,
                            egObject ***children, int **senses );
extern int  EG_indexBodyTopo( const egObject *body, const egObject *src );
extern int  EG_objectBodyTopo( const egObject *body, int oclass, int index,
                               egObject **obj );
extern int  EG_getBodyTopos( const egObject *body, /*@null@*/ egObject *src,
                             int oclass, int *ntopo,
                             /*@null@*/ egObject ***topos );
extern int  EG_attributeRet( const egObject *obj, const char *name, int *atype,
                             int *len, /*@null@*/ const int    **ints,
                                       /*@null@*/ const double **reals,
                                       /*@null@*/ const char   **str );
extern int  EG_attributeAdd( egObject *obj, const char *name, int type, int len,
                             /*@null@*/ const int    *ints,
                             /*@null@*/ const double *reals,
                             /*@null@*/ const char   *str );
extern int  EG_attributeDup( const egObject *src, egObject *dst );
extern int  EG_attributeCommon( const egObject *src, egObject *dst );
extern int  EG_attributeDel( egObject *obj, /*@null@*/ const char *name );
extern int  EG_fullAttrs( const egObject *obj );
extern int  EG_getWindingAngle( egObject *edge, double t, double *angle );
extern int  EG_getRangX( const egObject *geom, double *range, int *periodic );
extern int  EG_getAreX( egObject *object, /*@null@*/ const double *limits,
                        double *area );
extern int  EG_arcLenX( const egObject *geom, double t1, double t2,
                        double *alen );
extern int  EG_inFaceX( const egObject *face, const double *uv,
                        /*@null@*/ double *pt, /*@null@*/ double *uvx );
extern int  EG_getBoundingBX( const egObject *topo, double *bbox );
extern int  EG_curvaturX( const egObject *geom, const double *param,
                          double *result );
extern int  EG_massProperties( int nTopo, egObject **topos, double *data );
extern int  EG_getBody( const egObject *object, egObject **body );


static int
EG_effectNeighbor(egEFace *effect)
{
  int     i, j, n, i0, i1, i2, npts, ntris, *tris, *tric, *vtab, ipat = 0;
  connect *etab;
  
  if (effect->npatch != 1) return EGADS_INDEXERR;

  npts  = effect->patches[ipat].nuvs;
  ntris = effect->patches[ipat].ntris;
  tris  = effect->patches[ipat].uvtris;
  tric  = (int *) EG_alloc(3*ntris*sizeof(int));
  if (tric == NULL) return EGADS_MALLOC;
  vtab  = (int *) EG_alloc(npts*sizeof(int));
  if (vtab == NULL) {
    EG_free(tric);
    return EGADS_MALLOC;
  }
  etab  = (connect *) EG_alloc(ntris*3*sizeof(connect));
  if (etab == NULL) {
    EG_free(vtab);
    EG_free(tric);
    return EGADS_MALLOC;
  }
  
  n = NOTFILLED;
  for (j = 0; j < npts;  j++) vtab[j] = NOTFILLED;
  for (i = 0; i < ntris; i++)
    tric[3*i  ] = tric[3*i+1] = tric[3*i+2] = i+1;
  for (i = 0; i < ntris; i++) {
    i0 = tris[3*i  ];
    i1 = tris[3*i+1];
    i2 = tris[3*i+2];
    EG_makeConnect(i1, i2, &tric[3*i  ], &n, vtab, etab, 0);
    EG_makeConnect(i0, i2, &tric[3*i+1], &n, vtab, etab, 0);
    EG_makeConnect(i0, i1, &tric[3*i+2], &n, vtab, etab, 0);
  }
  /* find any unconnected triangle sides -- mark with zero */
  for (j = 0; j <= n; j++)
    if (etab[j].tri != NULL) *etab[j].tri = 0;

  EG_free(etab);
  EG_free(vtab);
  effect->patches[ipat].uvtric = tric;
  
  return EGADS_SUCCESS;
}


static int
EG_effectWalk(egEFace *effect, double *uv, int *tri, double *w)
{
  int     itri, i0, i1, i2, ntris, *tris, *tric, stat, cnt, last, ipat = 0;
  double  *tuv, tol = -1.e-6;
 
  ntris = effect->patches[ipat].ntris;
  if ((*tri <= 0) || (*tri > ntris)) return EGADS_INDEXERR;
  
  /* start from input triangle */
  itri  = *tri;
  tuv   = effect->patches[ipat].uvs;
  tris  = effect->patches[ipat].uvtris;
  tric  = effect->patches[ipat].uvtric;
  cnt   = 0;
  last  = itri;
  while (itri > 0) {
    cnt++;
    *tri = itri;
    i0   = tris[3*itri-3] - 1;
    i1   = tris[3*itri-2] - 1;
    i2   = tris[3*itri-1] - 1;
    stat = EG_inTriExact(&tuv[2*i0], &tuv[2*i1], &tuv[2*i2], uv, w);
    if ((stat == EGADS_SUCCESS) ||
        ((w[0] > tol) && (w[1] > tol) && (w[2] > tol))) {
 /*   printf(" %d: itri = %d  %d %d %d   %le %le %le  %d\n", cnt, itri,
             tric[3*itri-3], tric[3*itri-2], tric[3*itri-1], w[0], w[1], w[2],
             stat); */
      return EGADS_SUCCESS;
    }
    /* are we hopelessly lost? */
    if (cnt >= ntris) {
/*    printf(" %d/%d: itri = %d  %d %d %d   %le %le %le  %d\n", cnt, ntris, itri,
             tric[3*itri-3], tric[3*itri-2], tric[3*itri-1], w[0], w[1], w[2],
             stat);  */
      return EGADS_EXTRAPOL;
    }
    if ((w[0] < w[1]) && (w[0] < w[2])) {
      itri = tric[3*itri-3];
    } else if ((w[1] < w[0]) && (w[1] < w[2])) {
      itri = tric[3*itri-2];
    } else {
      itri = tric[3*itri-1];
    }
    /* are we flip/flopping? */
    if (itri == last) {
/*    printf(" %d: itri = %d  %d %d %d   %le %le %le -- %d\n", cnt, itri,
             tric[3*itri-3], tric[3*itri-2], tric[3*itri-1], w[0], w[1], w[2],
             last);  */
      return EGADS_EXTRAPOL;
    }
    last = *tri;
  }

  return EGADS_EXTRAPOL;
}


static int
EG_effectInTri(egEFace *effect, const double *uvx, int *itrix, double *w)
{
  int    stat, itri, i1, i2, i3, cls, ipat = 0;
  double uv[2], neg;
  
  uv[0] = uvx[0];
  uv[1] = uvx[1];
  
  /* lets start from last triangle */
  if ((effect->last != 0) && (effect->last <= effect->patches[ipat].ntris)) {
    stat   = EG_effectWalk(effect, uv, &effect->last, w);
    *itrix = effect->last;
    if (stat != EGADS_EXTRAPOL) return stat;
  }
  
  /* no hit -- exhaustive search */
  cls = 0;
  for (itri = 1; itri <= effect->patches[ipat].ntris; itri++) {
    i1   = effect->patches[ipat].uvtris[3*itri-3] - 1;
    i2   = effect->patches[ipat].uvtris[3*itri-2] - 1;
    i3   = effect->patches[ipat].uvtris[3*itri-1] - 1;
    stat = EG_inTriExact(&effect->patches[ipat].uvs[2*i1],
                         &effect->patches[ipat].uvs[2*i2],
                         &effect->patches[ipat].uvs[2*i3], uv, w);
    if (stat == EGADS_SUCCESS) {
      *itrix = effect->last = itri;
      return EGADS_SUCCESS;
    }
    if (w[1] < w[0]) w[0] = w[1];
    if (w[2] < w[0]) w[0] = w[2];
    if (cls == 0) {
      cls = itri;
      neg = w[0];
    } else {
      if (w[0] > neg) {
        cls = itri;
        neg = w[0];
      }
    }
  }

  /* extrapolate */
  i1 = effect->patches[ipat].uvtris[3*cls-3] - 1;
  i2 = effect->patches[ipat].uvtris[3*cls-2] - 1;
  i3 = effect->patches[ipat].uvtris[3*cls-1] - 1;
  EG_inTriExact(&effect->patches[ipat].uvs[2*i1],
                &effect->patches[ipat].uvs[2*i2],
                &effect->patches[ipat].uvs[2*i3], uv, w);
  *itrix = effect->last = cls;

  return EGADS_SUCCESS;
}


static int
EG_eFaceInterior(egEFace *effect, double *uv, egObject **face)
{
  int    stat, ipat, itri, i1, i2, i3, verts[3];
  double w[3];
  
  if (effect->npatch == 1) {
    ipat = 0;
    stat = EG_effectInTri(effect, uv, &itri, w);
    if (stat != EGADS_SUCCESS) return stat;
    i1    = effect->patches[ipat].uvtris[3*itri-3] - 1;
    i2    = effect->patches[ipat].uvtris[3*itri-2] - 1;
    i3    = effect->patches[ipat].uvtris[3*itri-1] - 1;
  } else {
    stat = EG_uvmapLocate(effect->uvmap, effect->trmap, uv, &ipat, &i1, verts,
                          w);
    if (stat != EGADS_SUCCESS) return stat;
    ipat--;
    itri  = i1 - effect->patches[ipat].start;
    i1    = effect->patches[ipat].uvtris[3*itri-3] - 1;
    i2    = effect->patches[ipat].uvtris[3*itri-2] - 1;
    i3    = effect->patches[ipat].uvtris[3*itri-1] - 1;
    uv[0] = w[0]*effect->patches[ipat].uvs[2*i1  ] +
            w[1]*effect->patches[ipat].uvs[2*i2  ] +
            w[2]*effect->patches[ipat].uvs[2*i3  ];
    uv[1] = w[0]*effect->patches[ipat].uvs[2*i1+1] +
            w[1]*effect->patches[ipat].uvs[2*i2+1] +
            w[2]*effect->patches[ipat].uvs[2*i3+1];
  }
  
  *face = effect->patches[ipat].face;
  if (i1 >= effect->patches[ipat].ndeflect) i1 = 0;
  if (i2 >= effect->patches[ipat].ndeflect) i2 = 0;
  if (i3 >= effect->patches[ipat].ndeflect) i3 = 0;
  if (i1+i2+i3 == 0) return EGADS_SUCCESS;
  
  return EGADS_OUTSIDE;
}


static int
EG_effect1EdgeEval(egEEdge *effect, double t, double *result, int *hit)
{
  int    i, ie, stat;
  double w1, dxyz1[3], w2, dxyz2[3];
  
  *hit = 0;
  ie   = effect->segs[0].iedge - 1;
  w1   = dxyz1[0] = dxyz1[1] = dxyz1[2] = 0.0;
  w2   = dxyz2[0] = dxyz2[1] = dxyz2[2] = 0.0;

  if (t < effect->sedges[ie].ts[1]) {
    dxyz1[0] = effect->sedges[ie].dstart[0];
    dxyz1[1] = effect->sedges[ie].dstart[1];
    dxyz1[2] = effect->sedges[ie].dstart[2];
    w1       = 1.0 - (t                        - effect->sedges[ie].ts[0])/
                     (effect->sedges[ie].ts[1] - effect->sedges[ie].ts[0]);
    if (w1 <  0.0) w1   = 0.0;
    if (w1 >  1.0) w1   = 1.0;
    if (w1 != 0.0) *hit = 1;
  }
  i = effect->sedges[ie].npts - 2;
  if (t > effect->sedges[ie].ts[i]) {
    dxyz2[0] = effect->sedges[ie].dend[0];
    dxyz2[1] = effect->sedges[ie].dend[1];
    dxyz2[2] = effect->sedges[ie].dend[2];
    w2       = (t                          - effect->sedges[ie].ts[i])/
               (effect->sedges[ie].ts[i+1] - effect->sedges[ie].ts[i]);
    if (w2 <  0.0) w2   = 0.0;
    if (w2 >  1.0) w2   = 1.0;
    if (w2 != 0.0) *hit = 1;
  }
  stat = EG_evaluatX(effect->sedges[ie].edge, &t, result);
  if (stat != EGADS_SUCCESS) return stat;
  if (*hit == 0) return EGADS_SUCCESS;
  
  result[0] += w1*dxyz1[0] + w2*dxyz2[0];
  result[1] += w1*dxyz1[1] + w2*dxyz2[1];
  result[2] += w1*dxyz1[2] + w2*dxyz2[2];
  
  return EGADS_SUCCESS;
}


static int
EG_effectEdgeEval(egEEdge *effect, double t, double *result, int *hit)
{
  int    stat, i, ie, iseg;
  double tx, w1, dxyz1[3], w2, dxyz2[3];
  
  *hit = 0;
  w1   = dxyz1[0] = dxyz1[1] = dxyz1[2] = 0.0;
  w2   = dxyz2[0] = dxyz2[1] = dxyz2[2] = 0.0;

  for (iseg = 0; iseg < effect->nsegs; iseg++)
    if (t <= effect->segs[iseg].tend) break;
  if (iseg == effect->nsegs) iseg--;
  ie = effect->segs[iseg].iedge - 1;
  
  /* get t in segment */
  tx = t - effect->segs[iseg].tstart;
  if (effect->segs[iseg].sense == SREVERSE) {
    tx  = effect->sedges[ie].ts[effect->sedges[ie].npts-1] - tx;
  } else {
    tx += effect->sedges[ie].ts[0];
  }

  if (tx < effect->sedges[ie].ts[1]) {
    dxyz1[0] = effect->sedges[ie].dstart[0];
    dxyz1[1] = effect->sedges[ie].dstart[1];
    dxyz1[2] = effect->sedges[ie].dstart[2];
    w1       = 1.0 - (tx                       - effect->sedges[ie].ts[0])/
                     (effect->sedges[ie].ts[1] - effect->sedges[ie].ts[0]);
    if (w1 <  0.0) w1   = 0.0;
    if (w1 >  1.0) w1   = 1.0;
    if (w1 != 0.0) *hit = 1;
  }
  i = effect->sedges[ie].npts - 2;
  if (tx > effect->sedges[ie].ts[i]) {
    dxyz2[0] = effect->sedges[ie].dend[0];
    dxyz2[1] = effect->sedges[ie].dend[1];
    dxyz2[2] = effect->sedges[ie].dend[2];
    w2       = (tx                         - effect->sedges[ie].ts[i])/
               (effect->sedges[ie].ts[i+1] - effect->sedges[ie].ts[i]);
    if (w2 <  0.0) w2   = 0.0;
    if (w2 >  1.0) w2   = 1.0;
    if (w2 != 0.0) *hit = 1;
  }
  stat = EG_evaluatX(effect->sedges[ie].edge, &tx, result);
  if (stat != EGADS_SUCCESS) return stat;
  
  /* adjust first derivatives if sense is reversed */
  if (effect->segs[iseg].sense == SREVERSE) {
    result[3] = -result[3];
    result[4] = -result[4];
    result[5] = -result[5];
  }
  if (*hit == 0) return EGADS_SUCCESS;
  
  result[0] += w1*dxyz1[0] + w2*dxyz2[0];
  result[1] += w1*dxyz1[1] + w2*dxyz2[1];
  result[2] += w1*dxyz1[2] + w2*dxyz2[2];
  
  return EGADS_SUCCESS;
}


static int
EG_effectEvalEdge(egEEdge *effect, double tx, double *result)
{
  int    stat, hit;
  double t, tminus[9], tplus[9], uvstep = 1.e-8;
/*
  double dt1[3], dt2[3], dtt[3]
*/
  
  if (effect->nsegs == 1) {

    stat = EG_effect1EdgeEval(effect, tx, result, &hit);
    if (stat != EGADS_SUCCESS) return stat;
    if (hit == 0) return EGADS_SUCCESS;
    
    t = tx - uvstep;
    stat = EG_effect1EdgeEval(effect, t, tminus,  &hit);
    if (stat != EGADS_SUCCESS) return stat;
    
    t = tx + uvstep;
    stat = EG_effect1EdgeEval(effect, t, tplus,   &hit);
    if (stat != EGADS_SUCCESS) return stat;
    
  } else {
    
    stat = EG_effectEdgeEval(effect, tx, result,  &hit);
    if (stat != EGADS_SUCCESS) return stat;
    if (hit == 0) return EGADS_SUCCESS;
    
    t = tx - uvstep;
    stat = EG_effectEdgeEval(effect, t, tminus,   &hit);
    if (stat != EGADS_SUCCESS) return stat;
    
    t = tx + uvstep;
    stat = EG_effectEdgeEval(effect, t, tplus,    &hit);
    if (stat != EGADS_SUCCESS) return stat;
    
  }
  
  result[3] = ( tplus[0]-tminus[0])/(2.0*uvstep);
  result[4] = ( tplus[1]-tminus[1])/(2.0*uvstep);
  result[5] = ( tplus[2]-tminus[2])/(2.0*uvstep);
/*
  dt1[0]    = (result[0]-tminus[0])/uvstep;
  dt2[0]    = ( tplus[0]-result[0])/uvstep;
  dt1[1]    = (result[1]-tminus[1])/uvstep;
  dt2[1]    = ( tplus[1]-result[1])/uvstep;
  dt1[2]    = (result[2]-tminus[2])/uvstep;
  dt2[2]    = ( tplus[2]-result[2])/uvstep;
  result[6] = (   dt2[0]-   dt1[0])/uvstep;
  result[7] = (   dt2[1]-   dt1[1])/uvstep;
  result[8] = (   dt2[2]-   dt1[2])/uvstep;
*/
  
  return EGADS_SUCCESS;
}


static int
EG_effect1FaceEval(egEFace *effect, const double *uv, double *data, int *flag)
{
  int    stat, itri, i1, i2, i3, ipat;
  double w[3], dxyz1[3], dxyz2[3], dxyz3[3];

  *flag = 0;
  if (effect->npatch != 1) return EGADS_TOPOCNT;
  ipat  = 0;                   /* only for single patch EFaces */
  stat  = EG_effectInTri(effect, uv, &itri, w);
  if (stat != EGADS_SUCCESS) return stat;
  stat = EG_evaluatX(effect->patches[ipat].face, uv, data);
  if (stat != EGADS_SUCCESS) return stat;
  i1 = effect->patches[ipat].uvtris[3*itri-3];
  i2 = effect->patches[ipat].uvtris[3*itri-2];
  i3 = effect->patches[ipat].uvtris[3*itri-1];
  if (i1 > effect->patches[ipat].ndeflect) i1 = 0;
  if (i2 > effect->patches[ipat].ndeflect) i2 = 0;
  if (i3 > effect->patches[ipat].ndeflect) i3 = 0;
  if (i1+i2+i3 == 0) return EGADS_SUCCESS;

  *flag    = 1;
  dxyz1[0] = dxyz1[1] = dxyz1[2] = 0.0;
  dxyz2[0] = dxyz2[1] = dxyz2[2] = 0.0;
  dxyz3[0] = dxyz3[1] = dxyz3[2] = 0.0;
  if (i1 != 0) {
    dxyz1[0] = effect->patches[ipat].deflect[3*i1-3];
    dxyz1[1] = effect->patches[ipat].deflect[3*i1-2];
    dxyz1[2] = effect->patches[ipat].deflect[3*i1-1];
  }
  if (i2 != 0) {
    dxyz2[0] = effect->patches[ipat].deflect[3*i2-3];
    dxyz2[1] = effect->patches[ipat].deflect[3*i2-2];
    dxyz2[2] = effect->patches[ipat].deflect[3*i2-1];
  }
  if (i3 != 0) {
    dxyz3[0] = effect->patches[ipat].deflect[3*i3-3];
    dxyz3[1] = effect->patches[ipat].deflect[3*i3-2];
    dxyz3[2] = effect->patches[ipat].deflect[3*i3-1];
  }
/*
  if ((w[0] < 0.0) || (w[1] < 0.0) || (w[2] < 0.0)) {
    double sum;
    if ((w[0] < -1.e-4) || (w[1] < -1.e-4) || (w[2] < -1.e-4)) *flag = 2;
    if (w[0] < 0.0) w[0] = 0.0;
    if (w[1] < 0.0) w[1] = 0.0;
    if (w[2] < 0.0) w[2] = 0.0;
    sum   = w[0] + w[1] + w[2];
    w[0] /= sum;
    w[1] /= sum;
    w[2] /= sum;
  }
*/
  data[0] += dxyz1[0]*w[0] + dxyz2[0]*w[1] + dxyz3[0]*w[2];
  data[1] += dxyz1[1]*w[0] + dxyz2[1]*w[1] + dxyz3[1]*w[2];
  data[2] += dxyz1[2]*w[0] + dxyz2[2]*w[1] + dxyz3[2]*w[2];
  
  return EGADS_SUCCESS;
}


static int
EG_effectFaceEval(egEFace *effect, const double *uvx, double *data, int *flag)
{
  int    i, i1, i2, i3, ix, itri, stat, verts[3];
  double w[3], uv[2], uvf[2], dxyz1[3], dxyz2[3], dxyz3[3];

  *flag = 0;
  uv[0] = uvx[0];
  uv[1] = uvx[1];
  stat  = EG_uvmapLocate(effect->uvmap, effect->trmap, uv, &ix, &itri, verts,
                         w);
  if (stat != EGADS_SUCCESS) {
    if (stat != EGADS_EXTRAPOL)
      printf(" EGADS Error: EG_uvmapLocate = %d\n", stat);
    return stat;
  }
  i      = itri - 1 - effect->patches[ix-1].start;
  i1     = effect->patches[ix-1].uvtris[3*i  ] - 1;
  i2     = effect->patches[ix-1].uvtris[3*i+1] - 1;
  i3     = effect->patches[ix-1].uvtris[3*i+2] - 1;
  uvf[0] = w[0]*effect->patches[ix-1].uvs[2*i1  ] +
           w[1]*effect->patches[ix-1].uvs[2*i2  ] +
           w[2]*effect->patches[ix-1].uvs[2*i3  ];
  uvf[1] = w[0]*effect->patches[ix-1].uvs[2*i1+1] +
           w[1]*effect->patches[ix-1].uvs[2*i2+1] +
           w[2]*effect->patches[ix-1].uvs[2*i3+1];
  stat = EG_evaluatX(effect->patches[ix-1].face, uvf, data);
  if (stat != EGADS_SUCCESS) return stat;

  if (i1 >= effect->patches[ix-1].ndeflect) i1 = 0;
  if (i2 >= effect->patches[ix-1].ndeflect) i2 = 0;
  if (i3 >= effect->patches[ix-1].ndeflect) i3 = 0;
  if (i1+i2+i3 == 0) return EGADS_SUCCESS;

  dxyz1[0] = dxyz1[1] = dxyz1[2] = 0.0;
  dxyz2[0] = dxyz2[1] = dxyz2[2] = 0.0;
  dxyz3[0] = dxyz3[1] = dxyz3[2] = 0.0;
  if (i1 != 0) {
    dxyz1[0] = effect->patches[ix-1].deflect[3*i1-3];
    dxyz1[1] = effect->patches[ix-1].deflect[3*i1-2];
    dxyz1[2] = effect->patches[ix-1].deflect[3*i1-1];
  }
  if (i2 != 0) {
    dxyz2[0] = effect->patches[ix-1].deflect[3*i2-3];
    dxyz2[1] = effect->patches[ix-1].deflect[3*i2-2];
    dxyz2[2] = effect->patches[ix-1].deflect[3*i2-1];
  }
  if (i3 != 0) {
    dxyz3[0] = effect->patches[ix-1].deflect[3*i3-3];
    dxyz3[1] = effect->patches[ix-1].deflect[3*i3-2];
    dxyz3[2] = effect->patches[ix-1].deflect[3*i3-1];
  }
/*
  if ((w[0] < 0.0) || (w[1] < 0.0) || (w[2] < 0.0)) {
    double sum;
    if ((w[0] < -1.e-4) || (w[1] < -1.e-4) || (w[2] < -1.e-4)) *flag = 2;
    if (w[0] < 0.0) w[0] = 0.0;
    if (w[1] < 0.0) w[1] = 0.0;
    if (w[2] < 0.0) w[2] = 0.0;
    sum   = w[0] + w[1] + w[2];
    w[0] /= sum;
    w[1] /= sum;
    w[2] /= sum;
  }
*/
  data[0] += dxyz1[0]*w[0] + dxyz2[0]*w[1] + dxyz3[0]*w[2];
  data[1] += dxyz1[1]*w[0] + dxyz2[1]*w[1] + dxyz3[1]*w[2];
  data[2] += dxyz1[2]*w[0] + dxyz2[2]*w[1] + dxyz3[2]*w[2];
  
  return EGADS_SUCCESS;
}


static int
EG_effectEvalFace(egEFace *effect, const double *uvx, double *result)
{
  int    flag, stat, status = EGADS_SUCCESS;
  double uvstep, uv[2], uminus[18], uplus[18], vminus[18], vplus[18];
  double du1[3], du2[3], dv1[3], dv2[3], uvp[18], uvm[18];
  
  if (effect->npatch == 1) {
    uvstep = 1.e-8;
    
    /* single Face in the object */
    stat   = EG_effect1FaceEval(effect, uvx, result, &flag);
    if (stat != EGADS_SUCCESS) return stat;
    if (flag == 0) return EGADS_SUCCESS;
    if (flag == 2) status = EGADS_EXTRAPOL;
    
    uv[0]  = uvx[0] - uvstep;
    uv[1]  = uvx[1];
    stat   = EG_effect1FaceEval(effect, uv, uminus,  &flag);
    if (stat != EGADS_SUCCESS) return stat;
    
    uv[0] += 2.0*uvstep;
    stat   = EG_effect1FaceEval(effect, uv, uplus,   &flag);
    if (stat != EGADS_SUCCESS) return stat;
    
    uv[0] -= uvstep;
    uv[1] -= uvstep;
    stat   = EG_effect1FaceEval(effect, uv, vminus,  &flag);
    if (stat != EGADS_SUCCESS) return stat;
    
    uv[1] += 2.0*uvstep;
    stat   = EG_effect1FaceEval(effect, uv, vplus,   &flag);
    if (stat != EGADS_SUCCESS) return stat;
    
    result[3] = (uplus[0]-uminus[0])/(2.0*uvstep);
    result[4] = (uplus[1]-uminus[1])/(2.0*uvstep);
    result[5] = (uplus[2]-uminus[2])/(2.0*uvstep);
    result[6] = (vplus[0]-vminus[0])/(2.0*uvstep);
    result[7] = (vplus[1]-vminus[1])/(2.0*uvstep);
    result[8] = (vplus[2]-vminus[2])/(2.0*uvstep);
    
    return status;
  }
    
  /* a composite */

  uvstep = 1.e-5;
  stat   = EG_effectFaceEval(effect, uvx, result, &flag);
  if (stat != EGADS_SUCCESS) return stat;
  if (flag == 2) status = EGADS_EXTRAPOL;
  
  uv[0]  = uvx[0] - uvstep;
  uv[1]  = uvx[1];
  stat   = EG_effectFaceEval(effect, uv,  uminus, &flag);
  if (stat != EGADS_SUCCESS) return stat;
  
  uv[0] += 2.0*uvstep;
  stat   = EG_effectFaceEval(effect, uv,  uplus,  &flag);
  if (stat != EGADS_SUCCESS) return stat;
  
  uv[0] -= uvstep;
  uv[1] -= uvstep;
  stat   = EG_effectFaceEval(effect, uv,  vminus, &flag);
  if (stat != EGADS_SUCCESS) return stat;
  
  uv[1] += 2.0*uvstep;
  stat   = EG_effectFaceEval(effect, uv,  vplus,  &flag);
  if (stat != EGADS_SUCCESS) return stat;

  uv[0] += uvstep;
  stat   = EG_effectFaceEval(effect, uv,  uvp,    &flag);
  if (stat != EGADS_SUCCESS) return stat;
  
  uv[0] -= 2.0*uvstep;
  stat   = EG_effectFaceEval(effect, uv,  uvm,    &flag);
  if (stat != EGADS_SUCCESS) return stat;
  
  result[ 3] = (uplus[0]-uminus[0])/(2.0*uvstep);
  result[ 4] = (uplus[1]-uminus[1])/(2.0*uvstep);
  result[ 5] = (uplus[2]-uminus[2])/(2.0*uvstep);
  result[ 6] = (vplus[0]-vminus[0])/(2.0*uvstep);
  result[ 7] = (vplus[1]-vminus[1])/(2.0*uvstep);
  result[ 8] = (vplus[2]-vminus[2])/(2.0*uvstep);

  du1[0]     = (result[0]-uminus[0])/uvstep;
  du2[0]     = ( uplus[0]-result[0])/uvstep;
  du1[1]     = (result[1]-uminus[1])/uvstep;
  du2[1]     = ( uplus[1]-result[1])/uvstep;
  du1[2]     = (result[2]-uminus[2])/uvstep;
  du2[2]     = ( uplus[2]-result[2])/uvstep;
  result[ 9] = (   du2[0]-   du1[0])/uvstep;
  result[10] = (   du2[1]-   du1[1])/uvstep;
  result[11] = (   du2[2]-   du1[2])/uvstep;
  dv1[0]     = (result[0]-vminus[0])/uvstep;
  dv2[0]     = ( vplus[0]-result[0])/uvstep;
  dv1[1]     = (result[1]-vminus[1])/uvstep;
  dv2[1]     = ( vplus[1]-result[1])/uvstep;
  dv1[2]     = (result[2]-vminus[2])/uvstep;
  dv2[2]     = ( vplus[2]-result[2])/uvstep;
  result[15] = (   dv2[0]-   dv1[0])/uvstep;
  result[16] = (   dv2[1]-   dv1[1])/uvstep;
  result[17] = (   dv2[2]-   dv1[2])/uvstep;
  du1[0]     = (   uvm[0]-uminus[0])/uvstep;
  du2[0]     = (   uvp[0]- uplus[0])/uvstep;
  du1[1]     = (   uvm[1]-uminus[1])/uvstep;
  du2[1]     = (   uvp[1]- uplus[1])/uvstep;
  du1[2]     = (   uvm[2]-uminus[2])/uvstep;
  du2[2]     = (   uvp[2]- uplus[2])/uvstep;
  result[12] = (   du2[0]-   du1[0])/(2.0*uvstep);
  result[13] = (   du2[1]-   du1[1])/(2.0*uvstep);
  result[14] = (   du2[2]-   du1[2])/(2.0*uvstep);
  
  return status;
}


static void
EG_dereference(egObject *eobj, int flag)
{
  int      stat, locked = 0;
  egObject *context;
  egCntxt  *cntx;
  
  context = EG_context(eobj);
  if (context == NULL) return;
  cntx    = (egCntxt *) context->blind;
  if ((cntx->mutex != NULL) && (flag == 1)) {
    if (!EMP_LockTest(cntx->mutex)) {
      locked = 1;
      EMP_LockSet(cntx->mutex);
    }
  }
  stat = EG_dereferenceObject(eobj, context);
  if ((cntx->mutex != NULL) && (locked == 1)) EMP_LockRelease(cntx->mutex);
  if (stat != EGADS_SUCCESS)
    printf(" EGADS Internal: EG_dereferenceObject = %d\n!", stat);
}


static int
EG_effeContained(egObject *obj, egObject *src)
{
  int      i, stat;
  egEShell *eshell;
  egEFace  *eface;
  egELoop  *eloop;
  
  if (src->oclass == ELOOPX) {
    eloop = (egELoop *) src->blind;
    if (obj->oclass == EEDGE) {
      for (i = 0; i < eloop->eedges.nobjs; i++)
        if (eloop->eedges.objs[i] == obj) return EGADS_SUCCESS;
    } else {
      for (i = 0; i < eloop->eedges.nobjs; i++) {
        stat = EG_effeContained(obj, eloop->eedges.objs[i]);
        if (stat == EGADS_SUCCESS) return stat;
      }
    }
  } else if (src->oclass == EFACE) {
    eface = (egEFace *) src->blind;
    if (obj->oclass == ELOOPX) {
      for (i = 0; i < eface->eloops.nobjs; i++)
        if (eface->eloops.objs[i] == obj) return EGADS_SUCCESS;
    } else {
      for (i = 0; i < eface->eloops.nobjs; i++) {
        stat = EG_effeContained(obj, eface->eloops.objs[i]);
        if (stat == EGADS_SUCCESS) return stat;
      }
    }
  } else if (src->oclass == ESHELL) {
    eshell = (egEShell *) src->blind;
    if (obj->oclass == EFACE) {
      for (i = 0; i < eshell->efaces.nobjs; i++)
        if (eshell->efaces.objs[i] == obj) return EGADS_SUCCESS;
    } else {
      for (i = 0; i < eshell->efaces.nobjs; i++) {
        stat = EG_effeContained(obj, eshell->efaces.objs[i]);
        if (stat == EGADS_SUCCESS) return stat;
      }
    }
  }
  
  return EGADS_OUTSIDE;
}


static void
EG_swapEBodies(egObject *ebo1, egObject *ebo2)
{
  int      i;
  egObject *obj;
  egEBody  *ebody1, *ebody2;
  
  /* only works for copies -- but not checked! */
  ebody1 = (egEBody *) ebo1->blind;
  ebody2 = (egEBody *) ebo2->blind;

/*@-kepttrans@*/
  for (i = 0; i < ebody1->eedges.nobjs; i++) {
    obj  = ebody1->eedges.objs[i];
    obj->topObj = ebo2;
    obj  = ebody2->eedges.objs[i];
    obj->topObj = ebo1;
  }
  for (i = 0; i < ebody1->eloops.nobjs; i++) {
    obj  = ebody1->eloops.objs[i];
    obj->topObj = ebo2;
    obj  = ebody2->eloops.objs[i];
    obj->topObj = ebo1;
  }
  for (i = 0; i < ebody1->efaces.nobjs; i++) {
    obj  = ebody1->efaces.objs[i];
    obj->topObj = ebo2;
    obj  = ebody2->efaces.objs[i];
    obj->topObj = ebo1;
  }
  for (i = 0; i < ebody1->eshells.nobjs; i++) {
    obj  = ebody1->eshells.objs[i];
    obj->topObj = ebo2;
    obj  = ebody2->eshells.objs[i];
    obj->topObj = ebo1;
  }
/*@+kepttrans@*/
  
  ebo1->blind   = ebody2;
  ebo2->blind   = ebody1;
}


int
EG_getETopology(const egObject *topo, egObject **geom, int *oclass,
                int *type, /*@null@*/ double *range, int *nChildren,
                egObject ***children, int **senses)
{
  egEBody  *ebody;
  egEShell *eshell;
  egEFace  *eface;
  egELoop  *eloop;
  egEEdge  *eedge;
  egTessel *btess;

  *geom      = NULL;
  *oclass    = *type = 0;
  *nChildren = 0;
  *children  = NULL;
  *senses    = NULL;
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->blind == NULL)        return EGADS_NODATA;
  if (topo->oclass < EEDGE)       return EGADS_NOTTOPO;
  *oclass = topo->oclass;
  *type   = topo->mtype;
  
  if (topo->oclass == EEDGE) {
    eedge = (egEEdge *) topo->blind;
    *children  = eedge->nodes;
    *nChildren = 1;
    if (topo->mtype == TWONODE) *nChildren = 2;
    if (range != NULL) {
      range[0] = eedge->trange[0];
      range[1] = eedge->trange[1];
    }
  } else if (topo->oclass == ELOOPX) {
    eloop      = (egELoop *) topo->blind;
    *nChildren = eloop->eedges.nobjs;
    *children  = eloop->eedges.objs;
    *senses    = eloop->senses;
  } else if (topo->oclass == EFACE) {
    eface      = (egEFace *) topo->blind;
    *nChildren = eface->eloops.nobjs;
    *children  = eface->eloops.objs;
    *senses    = eface->senses;
    if (range != NULL) {
      range[0] = eface->range[0];
      range[1] = eface->range[1];
      range[2] = eface->range[2];
      range[3] = eface->range[3];
    }
  } else if (topo->oclass == ESHELL) {
    eshell     = (egEShell *) topo->blind;
    *nChildren = eshell->efaces.nobjs;
    *children  = eshell->efaces.objs;
  } else {
    ebody = (egEBody *) topo->blind;
    if (ebody->done == 0) {
      btess = (egTessel *) ebody->ref->blind;
      if (btess != NULL) *geom = btess->src;
    } else {
      *geom = ebody->ref;
    }
    if (topo->mtype == FACEBODY) {
      *nChildren = ebody->efaces.nobjs;
      *children  = ebody->efaces.objs;
    } else {
      *nChildren = ebody->eshells.nobjs;
      *children  = ebody->eshells.objs;
      if (topo->mtype == SOLIDBODY) *senses = ebody->senses;
    }
    if (ebody->done == 0) return EGADS_OUTSIDE;
  }
  
  return EGADS_SUCCESS;
}


int
EG_getEBodyTopos(const egObject *body, /*@null@*/ egObject *src,
                 int oclass, int *ntopo, /*@null@*/ egObject ***topos)
{
  int      i, j, k, kk, n, ie;
  egObject **objs, *eobj, *eeobj;
  egEBody  *ebody;
  egEFace  *eface;
  egELoop  *eloop;
  egEEdge  *eedge;
  egEMap   map;

  *ntopo = 0;
  if (topos != NULL) *topos = NULL;
  if (body == NULL)                 return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (body->oclass != EBODY)        return EGADS_NOTBODY;
  if (body->blind == NULL)          return EGADS_NODATA;
  if (src != NULL) {
    /* selective body objects */
    if  (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
    if ((src->oclass < EEDGE) ||
        (src->oclass > ESHELL))     return EGADS_NOTTOPO;
    if  (src->oclass == oclass)     return EGADS_TOPOERR;
    if  (src->blind == NULL)        return EGADS_NODATA;
  }
  ebody = (egEBody *) body->blind;
  
  /* asking for real Topology */

  if ((oclass == NODE) || (oclass == EDGE) || (oclass == FACE)) {
    if (src == NULL) {

      /* all effective body objects */
      if  (oclass == NODE) {
        
        for (n = i = 0; i < ebody->eedges.nobjs; i++) {
          eobj = ebody->eedges.objs[i];
          if (eobj == NULL) continue;
          if (eobj->blind == NULL) continue;
          eedge = (egEEdge *) eobj->blind;
          n++;
          if (eedge->nodes[0] != eedge->nodes[1]) n++;
        }
        objs = (egObject **) EG_alloc(n*sizeof(egObject *));
        if (objs == NULL) return EGADS_MALLOC;
        for (n = i = 0; i < ebody->eedges.nobjs; i++) {
          eobj = ebody->eedges.objs[i];
          if (eobj == NULL) continue;
          if (eobj->blind == NULL) continue;
          eedge = (egEEdge *) eobj->blind;
          for (k = 0; k < n; k++)
             if (objs[k] == eedge->nodes[0]) break;
           if (k == n) {
             objs[n] = eedge->nodes[0];
             n++;
           }
           for (k = 0; k < n; k++)
             if (objs[k] == eedge->nodes[1]) break;
           if (k == n) {
             objs[n] = eedge->nodes[1];
             n++;
           }
        }
        if (topos != NULL) {
          *topos = objs;
        } else {
          EG_free(objs);
        }
        
      } else if (oclass == EDGE) {

        for (n = i = 0; i < ebody->eedges.nobjs; i++) {
          eobj = ebody->eedges.objs[i];
          if (eobj == NULL) continue;
          if (eobj->blind == NULL) continue;
          eedge = (egEEdge *) eobj->blind;
          n += eedge->nsegs;
        }
        if (topos != NULL) {
          objs = (egObject **) EG_alloc(n*sizeof(egObject *));
          if (objs == NULL) return EGADS_MALLOC;
          for (n = i = 0; i < ebody->eedges.nobjs; i++) {
            eobj = ebody->eedges.objs[i];
            if (eobj == NULL) continue;
            if (eobj->blind == NULL) continue;
            eedge = (egEEdge *) eobj->blind;
            for (j = 0; j < eedge->nsegs; j++, n++) {
              ie      = eedge->segs[j].iedge - 1;
              objs[n] = eedge->sedges[ie].edge;
            }
          }
          *topos = objs;
        }
        
      } else {

        for (n = i = 0; i < ebody->efaces.nobjs; i++) {
          eobj = ebody->efaces.objs[i];
          if (eobj == NULL) continue;
          if (eobj->blind == NULL) continue;
          eface = (egEFace *) eobj->blind;
          n += eface->npatch;
        }
        if (topos != NULL) {
          objs = (egObject **) EG_alloc(n*sizeof(egObject *));
          if (objs == NULL) return EGADS_MALLOC;
          for (n = i = 0; i < ebody->efaces.nobjs; i++) {
            eobj = ebody->efaces.objs[i];
            if (eobj == NULL) continue;
            if (eobj->blind == NULL) continue;
            eface = (egEFace *) eobj->blind;
            for (j = 0; j < eface->npatch; j++, n++)
              objs[n] = eface->patches[j].face;
          }
          *topos = objs;
        }
      }

      *ntopo = n;
      return EGADS_SUCCESS;
    }

    /* src contains effective Topology -- only look down */
    if ((oclass == FACE) && (src->oclass < EFACE)) return EGADS_TOPOERR;
    
    if  (oclass == NODE) {
      objs = NULL;
      if (src->oclass == EEDGE) {
        eedge = (egEEdge *) src->blind;
        n     = 2;
        if (eedge->nodes[0] == eedge->nodes[1]) n = 1;
        if (topos != NULL) {
          objs = (egObject **) EG_alloc(n*sizeof(egObject *));
          if (objs == NULL) return EGADS_MALLOC;
          objs[0] = eedge->nodes[0];
          if (n == 2) objs[1] = eedge->nodes[1];
          *topos  = objs;
        }
        *ntopo = n;
      } else {
        eface = (egEFace *) src->blind;
        for (n = i = 0; i < eface->eloops.nobjs; i++) {
          eobj = eface->eloops.objs[i];
          if (eobj == NULL) continue;
          if (eobj->blind == NULL) continue;
          eloop = (egELoop *) eobj->blind;
          for (j = 0; j < eloop->eedges.nobjs; j++) {
            eeobj = eface->eloops.objs[i];
            if (eeobj == NULL) continue;
            if (eeobj->blind == NULL) continue;
            n += 2;
          }
        }
        objs = (egObject **) EG_alloc(n*sizeof(egObject *));
        if (objs == NULL) return EGADS_MALLOC;
        for (i = 0; i < n; i++) objs[i] = NULL;
        for (n = i = 0; i < eface->eloops.nobjs; i++) {
          eobj = eface->eloops.objs[i];
          if (eobj == NULL) continue;
          if (eobj->blind == NULL) continue;
          eloop = (egELoop *) eobj->blind;
          for (j = 0; j < eloop->eedges.nobjs; j++) {
            eeobj = eface->eloops.objs[i];
            if (eeobj == NULL) continue;
            if (eeobj->blind == NULL) continue;
            eedge = (egEEdge *) eeobj->blind;
            for (k = 0; k < n; k++)
              if (objs[k] == eedge->nodes[0]) break;
            if (k == n) {
              objs[n] = eedge->nodes[0];
              n++;
            }
            for (k = 0; k < n; k++)
              if (objs[k] == eedge->nodes[1]) break;
            if (k == n) {
              objs[n] = eedge->nodes[1];
              n++;
            }
          }
        }
      }
      if (topos != NULL) {
        *topos = objs;
      } else {
        if (objs != NULL) EG_free(objs);
      }
      *ntopo = n;
      
    } else if (oclass == EDGE) {
      
      if (src->oclass == EEDGE) {
        eedge = (egEEdge *) src->blind;
        if (topos != NULL) {
          objs = (egObject **) EG_alloc(eedge->nsegs*sizeof(egObject *));
          if (objs == NULL) return EGADS_MALLOC;
          for (i = 0; i < eedge->nsegs; i++) {
            ie      = eedge->segs[i].iedge - 1;
            objs[i] = eedge->sedges[ie].edge;
          }
          *topos = objs;
        }
        *ntopo = eedge->nsegs;
      } else {
        eface = (egEFace *) src->blind;
        for (n = i = 0; i < eface->eloops.nobjs; i++) {
          eobj = eface->eloops.objs[i];
          if (eobj == NULL) continue;
          if (eobj->blind == NULL) continue;
          eloop = (egELoop *) eobj->blind;
          for (j = 0; j < eloop->eedges.nobjs; j++) {
            eeobj = eface->eloops.objs[i];
            if (eeobj == NULL) continue;
            if (eeobj->blind == NULL) continue;
            eedge = (egEEdge *) eeobj->blind;
            n += eedge->nsegs;
          }
        }
        objs = (egObject **) EG_alloc(n*sizeof(egObject *));
        if (objs == NULL) return EGADS_MALLOC;
        for (i = 0; i < n; i++) objs[i] = NULL;
        for (n = i = 0; i < eface->eloops.nobjs; i++) {
          eobj = eface->eloops.objs[i];
          if (eobj == NULL) continue;
          if (eobj->blind == NULL) continue;
          eloop = (egELoop *) eobj->blind;
          for (j = 0; j < eloop->eedges.nobjs; j++) {
            eeobj = eface->eloops.objs[i];
            if (eeobj == NULL) continue;
            if (eeobj->blind == NULL) continue;
            eedge = (egEEdge *) eeobj->blind;
            for (kk = 0; kk < eedge->nsegs; kk++) {
              ie = eedge->segs[kk].iedge - 1;
              for (k = 0; k < n; k++)
                if (eedge->sedges[ie].edge == objs[k]) break;
              if (k == n) {
                objs[n] = eedge->sedges[ie].edge;
                n++;
              }
            }
          }
        }
        if (topos != NULL) {
          *topos = objs;
        } else {
          EG_free(objs);
        }
        *ntopo = n;
      }
      
    } else {
      
      eface = (egEFace *) src->blind;
      if (topos != NULL) {
        objs = (egObject **) EG_alloc(eface->npatch*sizeof(egObject *));
        if (objs == NULL) return EGADS_MALLOC;
        for (i = 0; i < eface->npatch; i++) objs[i] = eface->patches[i].face;
        *topos = objs;
      }
      *ntopo = eface->npatch;
    }
    
    return EGADS_SUCCESS;
  }
  
  /* asking for effective Topology */

  if (oclass == EEDGE) {
    map = ebody->eedges;
  } else if (oclass == ELOOPX) {
    map = ebody->eloops;
  } else if (oclass == EFACE) {
    map = ebody->efaces;
  } else if (oclass == ESHELL) {
    map = ebody->eshells;
  } else {
    return EGADS_NOTTOPO;
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
  
  /* look down the tree */
  if (src->oclass > oclass) {
    for (n = i = 0; i < map.nobjs; i++)
      if (EG_effeContained(map.objs[i], src) == EGADS_SUCCESS) n++;
    if (n == 0) return EGADS_SUCCESS;
    if (topos == NULL) {
      *ntopo = n;
      return EGADS_SUCCESS;
    }
    objs = (egObject **) EG_alloc(n*sizeof(egObject *));
    if (objs == NULL) return EGADS_MALLOC;
    for (n = i = 0; i < map.nobjs; i++)
      if (EG_effeContained(map.objs[i], src) == EGADS_SUCCESS) {
        objs[n] = map.objs[i];
        n++;
      }
    *ntopo = n;
    *topos = objs;
    return EGADS_SUCCESS;
  }
  
  /* look up the tree */
  for (n = i = 0; i < map.nobjs; i++)
    if (EG_effeContained(src, map.objs[i]) == EGADS_SUCCESS) n++;
  if (n == 0) return EGADS_SUCCESS;
  if (topos == NULL) {
    *ntopo = n;
    return EGADS_SUCCESS;
  }
  objs = (egObject **) EG_alloc(n*sizeof(egObject *));
  if (objs == NULL) return EGADS_MALLOC;
  for (n = i = 0; i < map.nobjs; i++)
    if (EG_effeContained(src, map.objs[i]) == EGADS_SUCCESS) {
      objs[n] = map.objs[i];
      n++;
    }
  *ntopo = n;
  *topos = objs;
  return EGADS_SUCCESS;
}


int
EG_indexEBodyTopo(const egObject *body, const egObject *src)
{
  int      i, stat, n, index = 0;
  egEBody  *ebody;
  egObject **objs;
  
  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->oclass != EBODY)      return EGADS_NOTBODY;
  if (body->blind == NULL)        return EGADS_NODATA;
  if (src->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if (src->blind == NULL)         return EGADS_NODATA;
  ebody = (egEBody *) body->blind;
  
  if ((src->oclass == NODE) || (src->oclass == EDGE) || (src->oclass == FACE)) {
    stat = EG_getEBodyTopos(body, NULL, src->oclass, &n, &objs);
    if (stat != EGADS_SUCCESS) return stat;
    for (i = 0; i < n; i++)
      if (src == objs[i]) {
        index = i+1;
        break;
      }
    EG_free(objs);
  } else if ((src->oclass >= EEDGE) && (src->oclass <= ESHELL)) {
    if (src->oclass == EEDGE) {
      for (i = 0; i < ebody->eedges.nobjs; i++)
        if (src == ebody->eedges.objs[i]) {
          index = i+1;
          break;
        }
    } else if (src->oclass == ELOOPX) {
      for (i = 0; i < ebody->eloops.nobjs; i++)
        if (src == ebody->eloops.objs[i]) {
          index = i+1;
          break;
        }
    } else if (src->oclass == EFACE) {
      for (i = 0; i < ebody->efaces.nobjs; i++)
        if (src == ebody->efaces.objs[i]) {
          index = i+1;
          break;
        }
    } else {
      for (i = 0; i < ebody->eshells.nobjs; i++)
        if (src == ebody->eshells.objs[i]) {
          index = i+1;
          break;
        }
    }
  } else {
    return EGADS_NOTTOPO;
  }

  if (index == 0) index = EGADS_NOTFOUND;
  return index;
}


int
EG_objectEBodyTopo(const egObject *body, int oclass, int index, egObject **obj)
{
  int      n, stat;
  egEBody  *ebody;
  egObject **objs;
  egEMap   map;
  
  *obj = NULL;
  if  (body == NULL)                       return EGADS_NULLOBJ;
  if  (body->magicnumber != MAGIC)         return EGADS_NOTOBJ;
  if  (body->oclass != EBODY)              return EGADS_NOTBODY;
  if  (body->blind == NULL)                return EGADS_NODATA;
  if  (index <= 0)                         return EGADS_INDEXERR;
  ebody = (egEBody *) body->blind;
  
  if ((oclass == NODE) || (oclass == EDGE) || (oclass == FACE)) {
    stat = EG_getEBodyTopos(body, NULL, oclass, &n, &objs);
    if (stat != EGADS_SUCCESS) return stat;
    if (index > n) {
      EG_free(objs);
      return EGADS_INDEXERR;
    }
    *obj = objs[index-1];
    EG_free(objs);
  } else if ((oclass >= EEDGE) && (oclass <= ESHELL)) {
    if (oclass == EEDGE) {
      map = ebody->eedges;
    } else if (oclass == ELOOPX) {
      map = ebody->eloops;
    } else if (oclass == EFACE) {
      map = ebody->efaces;
    } else {
      map = ebody->eshells;
    }
    if (index > map.nobjs) return EGADS_INDEXERR;
    *obj = map.objs[index-1];
  } else {
    return EGADS_NOTTOPO;
  }
  
  return EGADS_SUCCESS;
}


int
EG_getERange(const egObject *obj, double *range, int *periodic)
{
  egEFace *eface;
  egEEdge *eedge;
  
  *periodic = 0;
  if  (obj == NULL)               return EGADS_NULLOBJ;
  if  (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((obj->oclass != EEDGE) &&
      (obj->oclass != EFACE))     return EGADS_NOTGEOM;
  if (obj->blind == NULL)         return EGADS_NODATA;
  
  if (obj->oclass == EFACE) {
    eface = (egEFace *) obj->blind;
    if (eface->npatch != 1) {
      range[0] = eface->range[0];
      range[1] = eface->range[1];
      range[2] = eface->range[2];
      range[3] = eface->range[3];
    } else {
      return EG_getRangX(eface->patches[0].face, range, periodic);
    }
  } else {
    eedge    = (egEEdge *) obj->blind;
    range[0] = eedge->trange[0];
    range[1] = eedge->trange[1];
    if (obj->mtype == ONENODE) *periodic = 1;
  }
  
  return EGADS_SUCCESS;
}


int
EG_getEArea(egObject *object, /*@null@*/ const double *limits, double *area)
{
  int     i, stat;
  double  a;
  egEFace *eface;
  
  *area = 0.0;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->blind  == NULL)       return EGADS_NODATA;
  if (object->oclass != EFACE)      return EGADS_TOPOERR;
  
  eface = (egEFace *) object->blind;
  for (i = 0; i < eface->npatch; i++) {
    stat = EG_getAreX(eface->patches[i].face, limits, &a);
    if (stat != EGADS_SUCCESS) return stat;
    *area += fabs(a);
  }
  
  return EGADS_SUCCESS;
}


int
EG_inEFace(const egObject *object, const double *uv)
{
  int     i, stat;
  egEFace *eface;

  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->blind  == NULL)       return EGADS_NODATA;
  if (object->oclass != EFACE)      return EGADS_TOPOERR;
  
  eface = (egEFace *) object->blind;
  for (i = 0; i < eface->npatch; i++) {
    stat = EG_inFaceX(eface->patches[i].face, uv, NULL, NULL);
    if (stat <= EGADS_SUCCESS) return stat;
    if (stat == EGADS_SUCCESS) return stat;
  }
  
  return EGADS_OUTSIDE;
}


int
EG_arcELength(const egObject *object, double t1, double t2, double *alen)
{
  int     i, ie, stat;
  double  len, ts, tf;
  egEEdge *eedge;
  
  *alen = 0.0;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->blind  == NULL)       return EGADS_NODATA;
  if (object->oclass != EEDGE)      return EGADS_TOPOERR;
  
  eedge = (egEEdge *) object->blind;
  for (i = 0; i < eedge->nsegs; i++) {
    ie = eedge->segs[i].iedge - 1;
    if  (t1 >= eedge->segs[i].tend) continue;
    if ((t1 <= eedge->segs[i].tend) && (t1 >= eedge->segs[i].tstart)) {
      ts = t1 - eedge->segs[i].tstart;
      if (eedge->segs[i].sense == SREVERSE) {
        ts  = eedge->sedges[ie].ts[eedge->sedges[ie].npts-1] - ts;
      } else {
        ts += eedge->sedges[ie].ts[0];
      }
    } else {
      if (eedge->segs[i].sense == SREVERSE) {
        ts = eedge->sedges[ie].ts[eedge->sedges[ie].npts-1];
      } else {
        ts = eedge->sedges[ie].ts[0];
      }
    }
    
    if ((t2 <= eedge->segs[i].tend) && (t2 >= eedge->segs[i].tstart)) {
      tf = t2 - eedge->segs[i].tstart;
      if (eedge->segs[i].sense == SREVERSE) {
        tf  = eedge->sedges[ie].ts[eedge->sedges[ie].npts-1] - tf;
      } else {
        tf += eedge->sedges[ie].ts[0];
      }
    } else {
      if (eedge->segs[i].sense != SREVERSE) {
        tf = eedge->sedges[ie].ts[eedge->sedges[ie].npts-1];
      } else {
        tf = eedge->sedges[ie].ts[0];
      }
    }
    if (eedge->segs[i].sense == SREVERSE) {
      len = ts;
      ts  = tf;
      tf  = len;
    }
    stat = EG_arcLenX(eedge->sedges[ie].edge, ts, tf, &len);
    if (stat != EGADS_SUCCESS) return stat;
    *alen += len;
    if (t2 <= eedge->segs[i].tend) break;
  }
  
  return EGADS_SUCCESS;
}


int
EG_eEvaluate(const egObject *object, const double *param, double *result)
{
  egEEdge *eedge;
  egEFace *eface;

  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (object->blind  == NULL)       return EGADS_NODATA;
  if ((object->oclass != EEDGE) &&
      (object->oclass != EFACE))     return EGADS_TOPOERR;
  
  if (object->oclass == EEDGE) {
    if (object->mtype == DEGENERATE) return EGADS_DEGEN;
    eedge = (egEEdge *) object->blind;
    return EG_effectEvalEdge(eedge, param[0], result);
  }
    
  eface = (egEFace *) object->blind;
  return EG_effectEvalFace(eface, param, result);
}


int
EG_invEEvaluate(const egObject *object, double *xyz, double *param,
                double *result)
{
  int     i, stat, tbeg, tend, ie, index = -1;
  double  data[18], uv[2], d, dist = 1.e200;
  egEEdge *eedge;
  egEFace *eface;
  
  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (object->blind  == NULL)       return EGADS_NODATA;
  if ((object->oclass != EEDGE) &&
      (object->oclass != EFACE))     return EGADS_TOPOERR;
  
  if (object->oclass == EEDGE) {
    if (object->mtype == DEGENERATE) return EGADS_DEGEN;
    param[0] = 0.0;
    eedge    = (egEEdge *) object->blind;
    for (i = 0; i < eedge->nsegs; i++) {
      ie   = eedge->segs[i].iedge - 1;
      stat = EG_invEvaluatX(eedge->sedges[ie].edge, xyz, uv, data);
      if (stat != EGADS_SUCCESS) return stat;
      d = sqrt((xyz[0]-data[0])*(xyz[0]-data[0]) +
               (xyz[1]-data[1])*(xyz[1]-data[1]) +
               (xyz[2]-data[2])*(xyz[2]-data[2]));
      if (index == -1) {
        index     = i;
        dist      = d;
        param[0]  = uv[0];
        result[0] = data[0];
        result[1] = data[1];
        result[2] = data[2];
      } else {
        if (d < dist) {
          index     = i;
          dist      = d;
          param[0]  = uv[0];
          result[0] = data[0];
          result[1] = data[1];
          result[2] = data[2];
        }
      }
    }
    ie = eedge->segs[index].iedge - 1;
    if (eedge->segs[index].sense == SREVERSE) {
      param[0]  = eedge->segs[index].tstart +
                  eedge->sedges[ie].ts[eedge->sedges[ie].npts-1] - param[0];
    } else {
      param[0] += eedge->segs[index].tstart - eedge->sedges[ie].ts[0];
    }
  } else {
    eface = (egEFace *) object->blind;
    for (i = 0; i < eface->npatch; i++) {
      stat = EG_invEvaluatX(eface->patches[i].face, xyz, uv, data);
      if (stat != EGADS_SUCCESS) return stat;
      d = sqrt((xyz[0]-data[0])*(xyz[0]-data[0]) +
               (xyz[1]-data[1])*(xyz[1]-data[1]) +
               (xyz[2]-data[2])*(xyz[2]-data[2]));
      if (index == -1) {
        index     = i;
        dist      = d;
        param[0]  = uv[0];
        param[1]  = uv[1];
        result[0] = data[0];
        result[1] = data[1];
        result[2] = data[2];
      } else {
        if (d < dist) {
          index     = i;
          dist      = d;
          param[0]  = uv[0];
          param[1]  = uv[1];
          result[0] = data[0];
          result[1] = data[1];
          result[2] = data[2];
        }
      }
    }
    /* transform back to uvmap */
    if (eface->npatch != 1) {
      uv[0] = param[0];
      uv[1] = param[1];
      tbeg  = eface->patches[index].start;
      tend  = tbeg + eface->patches[index].ntris;
      stat  = EG_uv2UVmap(eface->uvmap, eface->trmap, uv,
                          eface->patches[index].uvs,
                          eface->patches[index].uvtris, tbeg+1, tend, param);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }
  
  return EGADS_SUCCESS;
}


int
EG_getEEdgeUV(const egObject *face, const egObject *topo, int sensx, double t,
              double *uv)
{
  int      i, n, stat, oclass, mtype, iloop, iedge, ie, iseg, *sens, sense;
  int      senses[2], *iuv;
  double   tx, w, range[2], uv1[2], uv2[2];
  egObject *elobj, *geom, **nodes;
  egEFace  *eFace;
  egELoop  *eLoop;
  egEEdge  *eEdge;

  uv[0] = uv[1] = 0.0;
  sense = sensx;
  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->blind == NULL)        return EGADS_NODATA;
  if (face->oclass != EFACE)      return EGADS_NOTTOPO;
  eFace = (egEFace *) face->blind;
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->blind == NULL)        return EGADS_NODATA;
  if (topo->oclass == NODE) {
    if (eFace->npatch == 1)
      return EG_getEdgeUVX(eFace->patches[0].face, topo, sense, t, uv);
    for (iloop = 0; iloop < eFace->eloops.nobjs; iloop++) {
      elobj = eFace->eloops.objs[iloop];
      if (elobj == NULL) continue;
      eLoop = (egELoop *) elobj->blind;
      if (eLoop == NULL) continue;
      for (iedge = 0; iedge < eLoop->nedge; iedge++) {
        stat = EG_getTopology(eLoop->edgeUVs[iedge].edge, &geom, &oclass,
                              &mtype, range, &n, &nodes, &sens);
        if (stat != EGADS_SUCCESS) continue;
        if (nodes[0] == topo) {
          EG_getUVmap(eFace->uvmap, eLoop->edgeUVs[iedge].iuv[0], uv);
          return EGADS_SUCCESS;
        }
        if (n == 2)
          if (nodes[1] == topo) {
            i = eLoop->edgeUVs[iedge].npts-1;
            EG_getUVmap(eFace->uvmap, eLoop->edgeUVs[iedge].iuv[i], uv);
            return EGADS_SUCCESS;
          }
      }
    }
    return EGADS_NOTFOUND;
  }
  if (topo->oclass != EEDGE)      return EGADS_NOTTOPO;
  eEdge = (egEEdge *) topo->blind;
  
  /* find the EEdge in the ELoop */
  senses[0] = senses[1] = 0;
  for (iloop = 0; iloop < eFace->eloops.nobjs; iloop++) {
    elobj = eFace->eloops.objs[iloop];
    if (elobj == NULL) continue;
    eLoop = (egELoop *) elobj->blind;
    if (eLoop == NULL) continue;
    for (iedge = 0; iedge < eLoop->eedges.nobjs; iedge++)
      if (eLoop->eedges.objs[iedge] == topo) {
        if (eLoop->senses[iedge] < 0) senses[0]++;
        if (eLoop->senses[iedge] > 0) senses[1]++;
      }
  }
  if (senses[0]+senses[1] == 0) return EGADS_NOTFOUND;
  if (sense == 0) {
    if ((senses[0] != 0) && (senses[1] != 0)) return EGADS_TOPOERR;
    sense = SREVERSE;
    if (senses[1] != 0) sense = SFORWARD;
  }
  
  for (iloop = 0; iloop < eFace->eloops.nobjs; iloop++) {
    elobj = eFace->eloops.objs[iloop];
    if (elobj == NULL) continue;
    eLoop = (egELoop *) elobj->blind;
    if (eLoop == NULL) continue;
    for (iedge = 0; iedge < eLoop->eedges.nobjs; iedge++) {
      if ((eLoop->eedges.objs[iedge] == topo) &&
          (eLoop->senses[iedge]      == sense)) {
        /* found our EEdge -- get the segment */
        for (iseg = 0; iseg < eEdge->nsegs; iseg++)
          if (t <= eEdge->segs[iseg].tend) break;
        if (iseg == eEdge->nsegs) iseg--;
        ie = eEdge->segs[iseg].iedge - 1;
        /* get t in segment */
        tx = t - eEdge->segs[iseg].tstart;
        if (eEdge->segs[iseg].sense == SREVERSE) {
          tx  = eEdge->sedges[ie].ts[eEdge->sedges[ie].npts-1] - tx;
        } else {
          tx += eEdge->sedges[ie].ts[0];
        }
        /* treat degerenate EEdges */
        if (topo->mtype == DEGENERATE)
          if (eFace->npatch == 1)
            return EG_getEdgeUVX(eFace->patches[0].face,
                                 eEdge->sedges[ie].edge, sense, t, uv);
        if (topo->mtype == DEGENERATE) return EGADS_DEGEN;

        /* interpolate UV */
        iuv = NULL;
        for (i = 0; i < eLoop->nedge; i++)
          if (eEdge->sedges[ie].edge == eLoop->edgeUVs[i].edge)
            if (eEdge->segs[iseg].sense*sense == eLoop->edgeUVs[i].sense) {
              iuv = eLoop->edgeUVs[i].iuv;
              break;
            }
        if (iuv == NULL) return EGADS_NOTFOUND;
        for (i = 0; i < eEdge->sedges[ie].npts-1; i++)
          if ((tx <= eEdge->sedges[ie].ts[i+1]) ||
              (i == eEdge->sedges[ie].npts-2)) {
            w = (tx                        - eEdge->sedges[ie].ts[i]) /
                (eEdge->sedges[ie].ts[i+1] - eEdge->sedges[ie].ts[i]);
            if (w < 0.0) w = 0.0;
            if (w > 1.0) w = 1.0;
            if (eFace->npatch == 1) {
              uv1[0] = eFace->patches[0].uvs[2*iuv[i  ]-2];
              uv1[1] = eFace->patches[0].uvs[2*iuv[i  ]-1];
              uv2[0] = eFace->patches[0].uvs[2*iuv[i+1]-2];
              uv2[1] = eFace->patches[0].uvs[2*iuv[i+1]-1];
            } else {
              EG_getUVmap(eFace->uvmap, iuv[i  ], uv1);
              EG_getUVmap(eFace->uvmap, iuv[i+1], uv2);
            }
            uv[0] = (1.0-w)*uv1[0] + w*uv2[0];
            uv[1] = (1.0-w)*uv1[1] + w*uv2[1];
            return EGADS_SUCCESS;
          }
      }
    }
  }
  
  return EGADS_RANGERR;
}


int
EG_eCurvature(const egObject *geom, const double *param, double *result)
{
  int      i, stat;
  double   data[18], dir[3], d, s, *d1, *d2, norm[3], *der1[2], *der2[3];
  double   a, b, c, d11, d12, d21, d22, g11, g12, g21, g22, ud, vd, len, uv[2];
  egObject *face;
  egEFace  *eface;
  
  if (geom->oclass == EEDGE) {
  
    stat = EG_eEvaluate(geom, param, data);
    if (stat != EGADS_SUCCESS) return stat;
    
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
    
    if (geom->blind == NULL) return EGADS_NODATA;
    eface = (egEFace *) geom->blind;
    uv[0] = param[0];
    uv[1] = param[1];
    stat  = EG_eFaceInterior(eface, uv, &face);
    if (stat == EGADS_SUCCESS) return EG_curvaturX(face, uv, result);
    
    stat = EG_eEvaluate(geom, param, data);
    if (stat != EGADS_SUCCESS) return stat;
    
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
EG_eBoundingBox(const egObject *topo, double *bbox)
{
  int            i, n, stat;
  double         box[6];
  egObject       **objs, *body;
  egEBody        *effect;
  egTessel       *btess;
  const egObject *ebody;
  
  bbox[0] = bbox[1] = bbox[2] =  1.e200;
  bbox[3] = bbox[4] = bbox[5] = -1.e200;
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->blind == NULL)        return EGADS_NODATA;
  if (topo->oclass < EEDGE)       return EGADS_NOTTOPO;
  if (topo->oclass > EBODY)       return EGADS_NOTTOPO;
  ebody = topo;
  if (topo->oclass  != EBODY) ebody = topo->topObj;
  if (ebody->oclass != EBODY)     return EGADS_NOTBODY;
  
  if (topo->oclass == EBODY) {
    effect = (egEBody *) ebody->blind;
    if (effect == NULL) return EGADS_NODATA;
    if (effect->done == 0) {
      btess = (egTessel *) effect->ref->blind;
      if (btess  == NULL) return EGADS_NOTTESS;
      body = btess->src;
    } else {
      body = effect->ref;
    }
    return EG_getBoundingBX(body, bbox);
  } else if (topo->oclass >= EFACE) {
    stat = EG_getEBodyTopos(ebody, (egObject *) topo, FACE, &n, &objs);
  } else {
    stat = EG_getEBodyTopos(ebody, (egObject *) topo, EDGE, &n, &objs);
  }
  if (stat != EGADS_SUCCESS) return stat;
  
  for (i = 0; i < n; i++) {
    stat = EG_getBoundingBX(objs[i], box);
    if (stat != EGADS_SUCCESS) {
      EG_free(objs);
      return stat;
    }
    if (box[0] < bbox[0]) bbox[0] = box[0];
    if (box[1] < bbox[1]) bbox[1] = box[1];
    if (box[2] < bbox[2]) bbox[2] = box[2];
    if (box[3] > bbox[3]) bbox[3] = box[3];
    if (box[4] > bbox[4]) bbox[4] = box[4];
    if (box[5] > bbox[5]) bbox[5] = box[5];
  }
  
  EG_free(objs);
  return EGADS_SUCCESS;
}


int
EG_massEProps(const egObject *topo, double *data)
{
  int      stat, i, j, ie, n = 0;
  egObject **topos, *obj;
  egTessel *btess;
  egEEdge  *eedge;
  egELoop  *eloop;
  egEFace  *eface;
  egEShell *eshell;
  egEBody  *ebody;
  
  if (topo->blind == NULL) return EGADS_NODATA;
  
  if (topo->oclass == EEDGE) {
    eedge = (egEEdge *) topo->blind;
    if (eedge->nsegs == 1) {
      ie = eedge->segs[0].iedge - 1;
      return EG_massProperties(1, &eedge->sedges[ie].edge, data);
    }
    n = eedge->nsegs;
  } else if (topo->oclass == ELOOPX) {
    eloop = (egELoop *) topo->blind;
    for (i = 0; i < eloop->eedges.nobjs; i++) {
      obj = eloop->eedges.objs[i];
      if (obj == NULL) continue;
      eedge = (egEEdge *) obj->blind;
      if (eedge == NULL) continue;
      n += eedge->nsegs;
    }
  } else if (topo->oclass == EFACE) {
    eface = (egEFace *) topo->blind;
    if (eface->npatch == 1)
      return EG_massProperties(1, &eface->patches[0].face, data);
    n = eface->npatch;
  } else if (topo->oclass == ESHELL) {
    eshell = (egEShell *) topo->blind;
    for (i = 0; i < eshell->efaces.nobjs; i++) {
      obj = eshell->efaces.objs[i];
      if (obj == NULL) continue;
      eface = (egEFace *) obj->blind;
      if (eface == NULL) continue;
      n += eface->npatch;
    }
  } else {
    ebody = (egEBody *) topo->blind;
    if (ebody->ref == NULL) return EGADS_NODATA;
    if (ebody->done == 0) {
      btess = (egTessel *) ebody->ref->blind;
      if (btess == NULL) return EGADS_NODATA;
      obj = btess->src;
    } else {
      obj = ebody->ref;
    }
    return EG_massProperties(1, &obj, data);
  }
  
  if (n == 0) return EGADS_TOPOCNT;
  topos = (egObject **) EG_alloc(n*sizeof(egObject *));
  if (topos == NULL) return EGADS_MALLOC;
  
  if (topo->oclass == EEDGE) {
    eedge = (egEEdge *) topo->blind;
    for (i = 0; i < eedge->nsegs; i++) {
      ie       = eedge->segs[i].iedge - 1;
      topos[i] = eedge->sedges[ie].edge;
    }
  } else if (topo->oclass == ELOOPX) {
    eloop = (egELoop *) topo->blind;
    for (n = j = 0; j < eloop->eedges.nobjs; j++) {
      obj = eloop->eedges.objs[j];
      if (obj == NULL) continue;
      eedge = (egEEdge *) obj->blind;
      if (eedge == NULL) continue;
      for (i = 0; i < eedge->nsegs; i++, n++) {
        ie       = eedge->segs[i].iedge - 1;
        topos[n] = eedge->sedges[ie].edge;
      }
    }
  } else if (topo->oclass == EFACE) {
    eface = (egEFace *) topo->blind;
    for (i = 0; i < eface->npatch; i++) topos[i] = eface->patches[i].face;
  } else {
    eshell = (egEShell *) topo->blind;
    for (n = j = 0; j < eshell->efaces.nobjs; j++) {
      obj = eshell->efaces.objs[j];
      if (obj == NULL) continue;
      eface = (egEFace *) obj->blind;
      if (eface == NULL) continue;
      for (i = 0; i < eface->npatch; i++, n++) topos[n] = eface->patches[i].face;
    }
  }
  
  stat = EG_massProperties(n, topos, data);
  EG_free(topos);
  return stat;
}


int
EG_tolEObject(const egObject *topo, double *tol)
{
  int      i, j, k, ie;
  double   toler, tolf;
  egObject *obj;
  egEEdge  *eedge;
  egELoop  *eloop;
  egEFace  *eface;
  egEShell *eshell;
  egEBody  *ebody;
  
  if (topo->blind == NULL) return EGADS_NODATA;
  
  if (topo->oclass == EEDGE) {
    if (topo->mtype == DEGENERATE) return EGADS_SUCCESS;
    eedge = (egEEdge *) topo->blind;
    if (eedge == NULL) return EGADS_SUCCESS;
    for (i = 0; i < eedge->nsegs; i++) {
      ie    = eedge->segs[i].iedge - 1;
      toler = sqrt(eedge->sedges[ie].dstart[0]*eedge->sedges[ie].dstart[0] +
                   eedge->sedges[ie].dstart[1]*eedge->sedges[ie].dstart[1] +
                   eedge->sedges[ie].dstart[2]*eedge->sedges[ie].dstart[2]);
      if (toler > *tol) *tol = toler;
      toler = sqrt(eedge->sedges[ie].dend[0]*eedge->sedges[ie].dend[0] +
                   eedge->sedges[ie].dend[1]*eedge->sedges[ie].dend[1] +
                   eedge->sedges[ie].dend[2]*eedge->sedges[ie].dend[2]);
      if (toler > *tol) *tol = toler;
    }
  } else if (topo->oclass == ELOOPX) {
    eloop = (egELoop *) topo->blind;
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      obj = eloop->eedges.objs[j];
      if (obj == NULL) continue;
      if (obj->mtype == DEGENERATE) continue;
      eedge = (egEEdge *) obj->blind;
      if (eedge == NULL) continue;
      for (i = 0; i < eedge->nsegs; i++) {
        ie    = eedge->segs[i].iedge - 1;
        toler = sqrt(eedge->sedges[ie].dstart[0]*eedge->sedges[ie].dstart[0] +
                     eedge->sedges[ie].dstart[1]*eedge->sedges[ie].dstart[1] +
                     eedge->sedges[ie].dstart[2]*eedge->sedges[ie].dstart[2]);
        if (toler > *tol) *tol = toler;
        toler = sqrt(eedge->sedges[ie].dend[0]*eedge->sedges[ie].dend[0] +
                     eedge->sedges[ie].dend[1]*eedge->sedges[ie].dend[1] +
                     eedge->sedges[ie].dend[2]*eedge->sedges[ie].dend[2]);
        if (toler > *tol) *tol = toler;
      }
    }
  } else if (topo->oclass == EFACE) {
    eface = (egEFace *) topo->blind;
    if (eface == NULL) return EGADS_SUCCESS;
    for (j = 0; j < eface->npatch; j++) {
      tolf = eface->patches[j].tol;
      if (tolf < 0.0) {
        tolf = 0.0;
        for (i = 0; i < eface->patches[j].ndeflect; i++) {
          toler = eface->patches[j].deflect[3*i  ]*
                  eface->patches[j].deflect[3*i  ] +
                  eface->patches[j].deflect[3*i+1]*
                  eface->patches[j].deflect[3*i+1] +
                  eface->patches[j].deflect[3*i+2]*
                  eface->patches[j].deflect[3*i+2];
          if (toler > tolf) tolf = toler;
        }
        tolf = sqrt(tolf);
        eface->patches[j].tol = tolf;
      }
      if (tolf > *tol) *tol = tolf;
    }
  } else if (topo->oclass == ESHELL) {
    eshell = (egEShell *) topo->blind;
    for (k = 0; k < eshell->efaces.nobjs; k++) {
      obj = eshell->efaces.objs[k];
      if (obj == NULL) continue;
      eface = (egEFace *) obj->blind;
      if (eface == NULL) continue;
      for (j = 0; j < eface->npatch; j++) {
        tolf = eface->patches[j].tol;
        if (tolf < 0.0) {
          tolf = 0.0;
          for (i = 0; i < eface->patches[j].ndeflect; i++) {
            toler = eface->patches[j].deflect[3*i  ]*
                    eface->patches[j].deflect[3*i  ] +
                    eface->patches[j].deflect[3*i+1]*
                    eface->patches[j].deflect[3*i+1] +
                    eface->patches[j].deflect[3*i+2]*
                    eface->patches[j].deflect[3*i+2];
            if (toler > tolf) tolf = toler;
          }
          tolf = sqrt(tolf);
          eface->patches[j].tol = tolf;
        }
        if (tolf > *tol) *tol = tolf;
      }
    }
  } else {
    ebody = (egEBody *) topo->blind;
    if (topo->mtype == WIREBODY) {
      if (ebody->eloops.objs[0] == NULL) return EGADS_SUCCESS;
      eloop = (egELoop *) ebody->eloops.objs[0]->blind;
      for (j = 0; j < eloop->eedges.nobjs; j++) {
        obj = eloop->eedges.objs[j];
        if (obj == NULL) continue;
        if (obj->mtype == DEGENERATE) continue;
        eedge = (egEEdge *) obj->blind;
        if (eedge == NULL) continue;
        for (i = 0; i < eedge->nsegs; i++) {
          ie    = eedge->segs[i].iedge - 1;
          toler = sqrt(eedge->sedges[ie].dstart[0]*eedge->sedges[ie].dstart[0] +
                       eedge->sedges[ie].dstart[1]*eedge->sedges[ie].dstart[1] +
                       eedge->sedges[ie].dstart[2]*eedge->sedges[ie].dstart[2]);
          if (toler > *tol) *tol = toler;
          toler = sqrt(eedge->sedges[ie].dend[0]*eedge->sedges[ie].dend[0] +
                       eedge->sedges[ie].dend[1]*eedge->sedges[ie].dend[1] +
                       eedge->sedges[ie].dend[2]*eedge->sedges[ie].dend[2]);
          if (toler > *tol) *tol = toler;
        }
      }
    } else {
      for (k = 0; k < ebody->efaces.nobjs; k++) {
        obj = ebody->efaces.objs[k];
        if (obj == NULL) continue;
        eface = (egEFace *) obj->blind;
        if (eface == NULL) continue;
        for (j = 0; j < eface->npatch; j++) {
          tolf = eface->patches[j].tol;
          if (tolf < 0.0) {
            tolf = 0.0;
            for (i = 0; i < eface->patches[j].ndeflect; i++) {
              toler = eface->patches[j].deflect[3*i  ]*
                      eface->patches[j].deflect[3*i  ] +
                      eface->patches[j].deflect[3*i+1]*
                      eface->patches[j].deflect[3*i+1] +
                      eface->patches[j].deflect[3*i+2]*
                      eface->patches[j].deflect[3*i+2];
              if (toler > tolf) tolf = toler;
            }
            tolf = sqrt(tolf);
            eface->patches[j].tol = tolf;
          }
          if (tolf > *tol) *tol = tolf;
        }
      }
    }
  }
  
  return EGADS_SUCCESS;
}


void
EG_destroyEBody(egObject *ebobj, int flag)
{
  int      i, j, k;
  egObject *eobj;
  egEBody  *ebody;
  egEShell *eshell;
  egEFace  *eface;
  egELoop  *eloop;
  egEEdge  *eedge;
  
  if (ebobj == NULL) {
    printf(" EGADS Internal: NULL Object (EG_destroyEBody)!\n");
    return;
  }
  if (ebobj->magicnumber != MAGIC) {
    printf(" EGADS Internal: Not an Object (EG_destroyEBody)!\n");
    return;
  }
  if (ebobj->oclass != EBODY) {
    printf(" EGADS Internal: Object not an EBody (EG_destroyEBody)!\n");
    return;
  }
  ebody = (egEBody *) ebobj->blind;
  if (ebody == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_destroyEBody)!\n");
    return;
  }
  
  if (ebody->edges != NULL) {
    for (i = 0; i < ebody->nedge; i++)
      EG_free(ebody->edges[i].ts);
    EG_free(ebody->edges);
  }
  
  if (ebody->eshells.objs != NULL) {
    for (i = 0; i < ebody->eshells.nobjs; i++) {
      eobj = ebody->eshells.objs[i];
      if (eobj == NULL) continue;
      if (eobj->blind == NULL) continue;
      eshell = (egEShell *) eobj->blind;
      EG_free(eshell->efaces.objs);
      EG_free(eshell);
      eobj->blind = NULL;
      EG_dereference(eobj, flag);
    }
    EG_free(ebody->eshells.objs);
  }
  
  if (ebody->efaces.objs != NULL) {
    for (i = 0; i < ebody->efaces.nobjs; i++) {
      eobj = ebody->efaces.objs[i];
      if (eobj == NULL) continue;
      if (eobj->blind == NULL) continue;
      eface = (egEFace *) eobj->blind;
      if (eface->trmap != NULL) EG_free(eface->trmap);
      if (eface->uvmap != NULL) uvmap_struct_free(eface->uvmap);
      if (eface->patches != NULL) {
        for (j = 0; j < eface->npatch; j++) {
          EG_free(eface->patches[j].uvtric);
          EG_free(eface->patches[j].uvtris);
          EG_free(eface->patches[j].uvs);
          EG_free(eface->patches[j].deflect);
        }
        EG_free(eface->patches);
      }
      EG_free(eface->eloops.objs);
      if (eface->senses != NULL) EG_free(eface->senses);
      EG_free(eface);
      eobj->blind = NULL;
      EG_dereference(eobj, flag);
    }
    EG_free(ebody->efaces.objs);
  }
  
  if (ebody->eloops.objs != NULL) {
    for (i = 0; i < ebody->eloops.nobjs; i++) {
      eobj = ebody->eloops.objs[i];
      if (eobj == NULL) continue;
      if (eobj->blind == NULL) continue;
      eloop = (egELoop *) eobj->blind;
      if (eloop->edgeUVs != NULL) {
        for (k = 0; k < eloop->nedge; k++) EG_free(eloop->edgeUVs[k].iuv);
        EG_free(eloop->edgeUVs);
      }
      EG_free(eloop->eedges.objs);
      EG_free(eloop->senses);
      EG_free(eloop);
      eobj->blind = NULL;
      EG_dereference(eobj, flag);
    }
    EG_free(ebody->eloops.objs);
  }
  
  if (ebody->eedges.objs != NULL) {
    for (i = 0; i < ebody->eedges.nobjs; i++) {
      eobj = ebody->eedges.objs[i];
      if (eobj == NULL) continue;
      if (eobj->blind == NULL) continue;
      eedge = (egEEdge *) eobj->blind;
      if (eedge->segs != NULL) EG_free(eedge->segs);
      EG_free(eedge);
      eobj->blind = NULL;
      EG_dereference(eobj, flag);
    }
    EG_free(ebody->eedges.objs);
  }

  EG_dereferenceTopObj(ebody->ref, ebobj);
  if (ebody->senses != NULL) EG_free(ebody->senses);
  EG_free(ebody);
}


static void
EG_trimELoop(egObject *loop, egObject *e1, int sense, egObject *e2)
{
  int     j, k;
  egELoop *eloop;
  
  eloop = (egELoop *) loop->blind;
  for (k = j = 0; j < eloop->eedges.nobjs; j++)
    if (eloop->eedges.objs[j] != e2) {
      eloop->eedges.objs[k] = eloop->eedges.objs[j];
      eloop->senses[k]      = eloop->senses[j];
      if (eloop->eedges.objs[j] == e1) eloop->senses[k] *= sense;
      k++;
    }
  eloop->eedges.nobjs = k;
}


static int
EG_eRemoveNode(egObject *EBody)
{
  int          i, j, k, ie, index, stat, nnode, sen1, sen2, last, fullAttr;
  int          atype, alen, *ncnt, *ecnt, ei[2];
  double       result[9], maxdot, tang[3], range[2], dist, dot;
  egEBody      *ebody;
  egELoop      *eloop;
  egEEdge      *eedge, *ee1, *ee2, *eesav;
  egTessel     *btess;
  egObject     *body, *e1, *e2, *esav, *loops[2], **nodes;
  const int    *ints;
  const double *reals;
  const char   *str;

  if (EBody == NULL)                return EGADS_NULLOBJ;
  if (EBody->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if (EBody->oclass != EBODY)       return EGADS_NOTBODY;
  if (EG_sameThread(EBody))         return EGADS_CNTXTHRD;
  ebody = (egEBody *) EBody->blind;
  if (ebody == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_eRemoveNode)!\n");
    return EGADS_NODATA;
  }
  if (ebody->done == 1)             return EGADS_EXISTS;
  if (ebody->ref  == NULL)          return EGADS_NOTTESS;
  btess = (egTessel *) ebody->ref->blind;
  if (btess == NULL)                return EGADS_NODATA;
  body     = btess->src;
  maxdot   = cos(ebody->angle*PI/180.0);
  fullAttr = EG_fullAttrs(body);
  
  stat  = EG_getBodyTopos(body, NULL, NODE, &nnode, &nodes);
  if (stat != EGADS_SUCCESS) return stat;
  ncnt = (int *) EG_alloc(nnode*sizeof(int));
  if (ncnt == NULL) {
    printf(" EGADS Error: Malloc of %d Nodes (EG_eRemoveNode)!\n", nnode);
    return EGADS_MALLOC;
  }
  /* mark the Nodes to be kept */
  for (i = 0; i < nnode; i++) {
    ncnt[i] =  0;
    stat    = EG_attributeRet(nodes[i], ".Keep", &atype, &alen,
                              &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) continue;
    ncnt[i] = -1;
  }
  
  /* count the Number of Nodes w/ 2 EEdges */
  for (i = 0; i < ebody->eedges.nobjs; i++) {
    if (ebody->eedges.objs[i] == NULL) continue;
    eedge = (egEEdge *) ebody->eedges.objs[i]->blind;
    /* Edges with .Keep attribute have Nodes that cannot be removed! */
    if (eedge->nsegs == 1) {
      ie   = eedge->segs[0].iedge - 1;
      stat = EG_attributeRet(eedge->sedges[ie].edge, ".Keep", &atype, &alen,
                             &ints, &reals, &str);
      if (stat == EGADS_SUCCESS) {
        index = EG_indexBodyTopo(body, eedge->nodes[0]);
        if (index <= EGADS_SUCCESS) {
          printf(" EGADS Internal: Node0 indexBodyTopo = %d (EG_eRemoveNode)!\n",
                 index);
        } else {
          ncnt[index-1] = -1;
        }
        index = EG_indexBodyTopo(body, eedge->nodes[1]);
        if (index <= EGADS_SUCCESS) {
          printf(" EGADS Internal: Node1 indexBodyTopo = %d (EG_eRemoveNode)!\n",
                 index);
        } else {
          ncnt[index-1] = -1;
        }
        continue;
      }
    }
    /* Only 2 Node EEdges are candidates */
    if (ebody->eedges.objs[i]->mtype != TWONODE) {
      index = EG_indexBodyTopo(body, eedge->nodes[0]);
      if (index <= EGADS_SUCCESS) {
        printf(" EGADS Internal: Node indexBodyTopo = %d (EG_eRemoveNode)!\n",
               index);
      } else {
        ncnt[index-1] = -1;
      }
      continue;
    }
    /* Mark the 2 Nodes */
    index = EG_indexBodyTopo(body, eedge->nodes[0]);
    if (index <= EGADS_SUCCESS) {
      printf(" EGADS Internal: Node- indexBodyTopo = %d (EG_eRemoveNode)!\n",
             index);
    } else {
      if (ncnt[index-1] >= 0) ncnt[index-1]++;
    }
    index = EG_indexBodyTopo(body, eedge->nodes[1]);
    if (index <= EGADS_SUCCESS) {
      printf(" EGADS Internal: Node+ indexBodyTopo = %d (EG_eRemoveNode)!\n",
             index);
    } else {
      if (ncnt[index-1] >= 0) ncnt[index-1]++;
    }
  }
  
  /* look at exposed EEdges */
  if (body->mtype == SHEETBODY) {
    ecnt = (int *) EG_alloc(ebody->eedges.nobjs*sizeof(int));
    if (ecnt == NULL) {
      EG_free(ncnt);
      EG_free(nodes);
      printf(" EGADS Error: Malloc of %d EEdges (EG_eRemoveNode)!\n",
             ebody->eedges.nobjs);
      return EGADS_MALLOC;
    }
    for (i = 0; i < ebody->eedges.nobjs; i++) ecnt[i] = 0;
    for (j = 0; j < ebody->eloops.nobjs; j++) {
      if (ebody->eloops.objs[j] == NULL) continue;
      eloop = (egELoop *) ebody->eloops.objs[j]->blind;
      if (eloop == NULL) continue;
      for (i = 0; i < eloop->eedges.nobjs; i++) {
        index = EG_indexEBodyTopo(EBody, eloop->eedges.objs[i]);
        if (index <= EGADS_SUCCESS) {
          printf(" EGADS Internal: EEdge indexBodyTopo = %d (EG_eRemoveNode)!\n",
                 index);
        } else {
          ecnt[index-1]++;
        }
      }
    }
    
    /* look at the Node tangents for single exposed EEdges */
    for (i = 0; i < nnode; i++)
      if (ncnt[i] == 2) {
        tang[0] = tang[1] = tang[2] = 0.0;
        for (j = 0; j < ebody->eedges.nobjs; j++) {
          if (ebody->eedges.objs[j] == NULL) continue;
          eedge = (egEEdge *) ebody->eedges.objs[j]->blind;
          if (ebody->eedges.objs[j]->mtype != TWONODE) continue;
          if (ecnt[j] != 1) continue;
          if (i+1 == EG_indexBodyTopo(body, eedge->nodes[0])) {
            stat = EG_eEvaluate(ebody->eedges.objs[j], &eedge->trange[0],
                                result);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Internal: EG_eEvaluate = %d (EG_eRemoveNode)!\n",
                     stat);
              continue;
            }
            dist = sqrt(result[3]*result[3] + result[4]*result[4] +
                        result[5]*result[5]);
            if (dist != 0.0) {
              result[3] /= dist;
              result[4] /= dist;
              result[5] /= dist;
            }
            if ((tang[0] == 0.0) && (tang[1] == 0.0) && (tang[2] == 0.0)) {
              tang[0] = result[3];
              tang[1] = result[4];
              tang[2] = result[5];
            } else {
              dot = result[3]*tang[0] + result[4]*tang[1] + result[5]*tang[2];
/*            printf(" Node %2d: dot = %lf (%lf)\n", i+1, dot, maxdot);  */
              if ((dot < maxdot) && (dot > -maxdot)) ncnt[i] = -2;
            }
          }
          if (i+1 == EG_indexBodyTopo(body, eedge->nodes[1])) {
            stat = EG_eEvaluate(ebody->eedges.objs[j], &eedge->trange[1],
                                result);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Internal: EG_eEvaluate = %d (EG_eRemoveNode)!\n",
                     stat);
              continue;
            }
            dist = sqrt(result[3]*result[3] + result[4]*result[4] +
                        result[5]*result[5]);
            if (dist != 0.0) {
              result[3] /= dist;
              result[4] /= dist;
              result[5] /= dist;
            }
            if ((tang[0] == 0.0) && (tang[1] == 0.0) && (tang[2] == 0.0)) {
              tang[0] = result[3];
              tang[1] = result[4];
              tang[2] = result[5];
            } else {
              dot = result[3]*tang[0] + result[4]*tang[1] + result[5]*tang[2];
/*            printf(" Node %2d: dot = %lf (%lf)\n", i+1, dot, maxdot);  */
              if ((dot < maxdot) && (dot > -maxdot)) ncnt[i] = -2;
            }
          }
        }
      }

    EG_free(ecnt);
  }
  
  /* adjust the Effective Topology */

  for (i = 0; i < nnode; i++)
    if (ncnt[i] == 2) {
#ifdef DEBUG
      printf(" Node %2d is being removed!\n", i+1);
#endif
      ei[0] = ei[1] = 0;
      for (j = 0; j < ebody->eedges.nobjs; j++) {
        if (ebody->eedges.objs[j] == NULL) continue;
        eedge = (egEEdge *) ebody->eedges.objs[j]->blind;
        if (ebody->eedges.objs[j]->mtype != TWONODE) continue;
        if (i+1 == EG_indexBodyTopo(body, eedge->nodes[0]))
          if (ei[0] == 0) {
            ei[0] = j+1;
          } else {
            ei[1] = j+1;
            break;
          }
        if (i+1 == EG_indexBodyTopo(body, eedge->nodes[1]))
          if (ei[0] == 0) {
            ei[0] = j+1;
          } else {
            ei[1] = j+1;
            break;
          }
      }
      if (ei[1] == 0) {
/*      printf(" EGADS Internal: EEdges = %d %d for Node %d (EG_eRemoveNode)!\n",
               ei[0], ei[1], i+1);  */
        continue;
      }
      loops[0] = loops[1] = NULL;
      for (j = 0; j < ebody->eloops.nobjs; j++) {
        if (ebody->eloops.objs[j] == NULL) continue;
        eloop = (egELoop *) ebody->eloops.objs[j]->blind;
        if (eloop == NULL) continue;
        for (k = 0; k < eloop->eedges.nobjs; k++) {
          index = EG_indexEBodyTopo(EBody, eloop->eedges.objs[k]);
          if (index == ei[0]) {
            if (loops[0] == NULL) {
              loops[0] = ebody->eloops.objs[j];
            } else {
              if (loops[0] == ebody->eloops.objs[j]) continue;
              loops[1] = ebody->eloops.objs[j];
              break;
            }
          }
          if (index == ei[1]) {
            if (loops[0] == NULL) {
              loops[0] = ebody->eloops.objs[j];
            } else {
              if (loops[0] == ebody->eloops.objs[j]) continue;
              loops[1] = ebody->eloops.objs[j];
              break;
            }
          }
        }
      }
      if (loops[0] == NULL) {
        printf(" EGADS Internal: Found No Loops for Node %d (EG_eRemoveNode)!\n",
               i+1);
        continue;
      }
      
      /* grab our 2 EEdges */
      e1    = ebody->eedges.objs[abs(ei[0])-1];
      ee1   = (egEEdge *) e1->blind;
      e2    = ebody->eedges.objs[abs(ei[1])-1];
      ee2   = (egEEdge *) e2->blind;
#ifdef DEBUG
      eloop = (egELoop *) loops[0]->blind;
      printf("     ELoop0 %d  %d\n", ei[0], ei[1]);
      printf("            %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[0]));
      for (k = 1; k < eloop->eedges.nobjs; k++)
        printf(" %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[k]));
      printf("\n");
      if (loops[1] != NULL) {
        eloop = (egELoop *) loops[1]->blind;
        printf("     ELoop1 %d  %d\n", ei[0], ei[1]);
        printf("            %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[0]));
        for (k = 1; k < eloop->eedges.nobjs; k++)
          printf(" %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[k]));
        printf("\n");
      }
      printf("     E %3d Nodes = %d %d\n", ei[0],
             EG_indexBodyTopo(body, ee1->nodes[0]),
             EG_indexBodyTopo(body, ee1->nodes[1]));
      printf("     E %3d Nodes = %d %d\n", ei[1],
             EG_indexBodyTopo(body, ee2->nodes[0]),
             EG_indexBodyTopo(body, ee2->nodes[1]));
#endif
      
      /* concatinate the segments */
      eedge = (egEEdge *) EG_alloc(sizeof(egEEdge));
      if (eedge == NULL) {
        printf(" EGADS Error: Malloc on %d EEdge blind (EG_eRemoveNode)!\n",
               i+1);
        EG_free(ncnt);
        EG_free(nodes);
        return EGADS_MALLOC;
      }
      eedge->sedges   = ebody->edges;
      eedge->nodes[0] = eedge->nodes[1] = NULL;
      if (EG_indexBodyTopo(body, ee1->nodes[0]) != i+1)
        eedge->nodes[0] = ee1->nodes[0];
      if (EG_indexBodyTopo(body, ee1->nodes[1]) != i+1)
        eedge->nodes[1] = ee1->nodes[1];
      last = 1;
      if (eedge->nodes[0] == NULL) {
        if (EG_indexBodyTopo(body, ee2->nodes[0]) != i+1)
          eedge->nodes[0] = ee2->nodes[0];
        if (EG_indexBodyTopo(body, ee2->nodes[1]) != i+1)
          eedge->nodes[0] = ee2->nodes[1];
      }
      if (eedge->nodes[1] == NULL) {
        if (EG_indexBodyTopo(body, ee2->nodes[0]) != i+1)
          eedge->nodes[1] = ee2->nodes[0];
        if (EG_indexBodyTopo(body, ee2->nodes[1]) != i+1)
          eedge->nodes[1] = ee2->nodes[1];
        last = 2;
      }
      if ((eedge->nodes[0] == NULL) || (eedge->nodes[1] == NULL)) {
        printf(" EGADS Internal: Cant find Node %d in Edges (EG_eRemoveNode)!\n",
               i+1);
        EG_free(eedge);
        EG_free(ncnt);
        EG_free(nodes);
        return EGADS_TOPOERR;
      }
      if (last == 1) {
        eesav = ee1;
        ee1   = ee2;
        ee2   = eesav;
        esav  = e1;
        e1    = e2;
        e2    = esav;
      }
      eedge->trange[0] = ee1->trange[0];
      eedge->trange[1] = ee1->trange[1] + ee2->trange[1] - ee2->trange[0];
      eedge->nsegs = ee1->nsegs + ee2->nsegs;
      eedge->segs  = (egEEseg *) EG_alloc(eedge->nsegs*sizeof(egEEseg));
      if (eedge->segs == NULL) {
        printf(" EGADS Error: Malloc on %d EEdge %d segs (EG_eRemoveNode)!\n",
               i+1, eedge->nsegs);
        EG_free(eedge);
        EG_free(ncnt);
        EG_free(nodes);
        return EGADS_MALLOC;
      }
      sen1 = sen2 = SFORWARD;
      if (EG_indexBodyTopo(body, ee1->nodes[0]) == i+1) sen1 = SREVERSE;
      if (EG_indexBodyTopo(body, ee2->nodes[1]) == i+1) sen2 = SREVERSE;
#ifdef DEBUG
      printf("     E new Nodes = %d %d   last = %d  fill = %d %d\n",
             EG_indexBodyTopo(body, eedge->nodes[0]),
             EG_indexBodyTopo(body, eedge->nodes[1]), last, sen1, sen2);
#endif

      /* make it go from Node to Node */
      for (j = 0; j < ee1->nsegs; j++) {
        if (sen1 < 0) {
          eedge->segs[j].iedge  =  ee1->segs[ee1->nsegs-j-1].iedge;
          eedge->segs[j].sense  = -ee1->segs[ee1->nsegs-j-1].sense;
          eedge->segs[j].nstart =  ee1->segs[ee1->nsegs-j-1].nstart;
        } else {
          eedge->segs[j].iedge  = ee1->segs[j].iedge;
          eedge->segs[j].sense  = ee1->segs[j].sense;
          eedge->segs[j].nstart = ee1->segs[j].nstart;
        }
        ie       = eedge->segs[j].iedge - 1;
        range[0] = eedge->sedges[ie].ts[0];
        range[1] = eedge->sedges[ie].ts[eedge->sedges[ie].npts-1];
        eedge->segs[j].tstart = range[0];
        eedge->segs[j].tend   = range[1];
        if (j > 0) {
          eedge->segs[j].tstart = eedge->segs[j-1].tend;
          eedge->segs[j].tend   = eedge->segs[j].tstart + range[1] - range[0];
        }
      }
      k = j;
      for (j = 0; j < ee2->nsegs; j++, k++) {
        if (sen2 < 0) {
          eedge->segs[k].iedge  =  ee2->segs[ee2->nsegs-j-1].iedge;
          eedge->segs[k].sense  = -ee2->segs[ee2->nsegs-j-1].sense;
          eedge->segs[k].nstart =  ee2->segs[ee2->nsegs-j-1].nstart;
        } else {
          eedge->segs[k].iedge  = ee2->segs[j].iedge;
          eedge->segs[k].sense  = ee2->segs[j].sense;
          eedge->segs[k].nstart = ee2->segs[j].nstart;
        }
        if (j == 0) eedge->segs[k].nstart = nodes[i];
        ie = eedge->segs[k].iedge - 1;
        range[0]              = eedge->sedges[ie].ts[0];
        range[1]              = eedge->sedges[ie].ts[eedge->sedges[ie].npts-1];
        eedge->segs[k].tstart = eedge->segs[k-1].tend;
        eedge->segs[k].tend   = eedge->segs[k].tstart + range[1] - range[0];
      }
#ifdef DEBUG
      printf(" New EEdge nseg = %d\n", k);
      printf("           tRange = %lf %lf  nsegs = %d   nodes = %d %d\n",
             eedge->trange[0], eedge->trange[1], eedge->nsegs,
             EG_indexBodyTopo(body, eedge->nodes[0]),
             EG_indexBodyTopo(body, eedge->nodes[1]));
      for (j = 0; j < eedge->nsegs; j++) {
        printf("            Edge = %2d: sense = %d  tstart = %lf  tend = %lf\n",
               eedge->segs[j].iedge,  eedge->segs[j].sense,
               eedge->segs[j].tstart, eedge->segs[j].tend);
        if (eedge->segs[j].nstart != NULL)
          printf("                       Internal Node = %d\n",
                 EG_indexBodyTopo(body, eedge->segs[j].nstart));
        ie = eedge->segs[j].iedge - 1;
        printf("                       Ts = %lf %lf\n", eedge->sedges[j].ts[0],
               eedge->sedges[ie].ts[eedge->sedges[ie].npts-1]);
      }
      stat = EG_evaluate(eedge->nodes[0], NULL, result);
      printf("            %lf:  %lf %lf %lf   %d\n",
             eedge->trange[0], result[0], result[1], result[2], stat);
      for (j = 0; j < 11; j++) {
        double t;
        
        t = eedge->trange[0] + j*(eedge->trange[1] - eedge->trange[0])/10.0;
        stat = EG_effectEdgeEval(eedge, t, result, &k);
        if (stat != EGADS_SUCCESS) {
          printf("            %lf:  Error = %d\n", t, stat);
          continue;
        }
        printf("            %lf:  %lf %lf %lf  %d\n",
               t, result[0], result[1], result[2], k);
      }
      stat = EG_evaluate(eedge->nodes[1], NULL, result);
      printf("            %lf:  %lf %lf %lf   %d\n",
             eedge->trange[1], result[0], result[1], result[2], stat);
      
      printf(" Removing EEdge %d\n", EG_indexBodyTopo(EBody, e2));
#endif
      
      /* adjust the loops */
      EG_trimELoop(loops[0], e1, sen1, e2);
      if (loops[1] != NULL) EG_trimELoop(loops[1], e1, sen1, e2);
      
      /* adjust the Edges */
      EG_free(ee1->segs);
      EG_free(ee1);
      e1->blind = eedge;
      if (eedge->nodes[0] == eedge->nodes[1]) e1->mtype = ONENODE;
      for (j = k = 0; j < ebody->eedges.nobjs; j++)
        if (ebody->eedges.objs[j] != e2) {
          ebody->eedges.objs[k] = ebody->eedges.objs[j];
          k++;
        }
      if (fullAttr == 1) {
        EG_attributeDup(e2, e1);
      } else {
        EG_attributeCommon(e2, e1);
      }
      ebody->eedges.nobjs = k;
      EG_free(ee2->segs);
      EG_free(ee2);
      e2->blind = NULL;
      EG_dereference(e2, 1);
#ifdef DEBUG
      eloop = (egELoop *) loops[0]->blind;
      printf("     ELoop0\n");
      printf("            %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[0]));
      for (k = 1; k < eloop->eedges.nobjs; k++)
        printf(" %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[k]));
      printf("\n");
      printf("            ");
      for (k = 0; k < eloop->eedges.nobjs; k++) {
        ee2 = (egEEdge *) eloop->eedges.objs[k]->blind;
        if (eloop->senses[k] == SFORWARD) {
          printf("%d,%d ", EG_indexBodyTopo(body, ee2->nodes[0]),
                           EG_indexBodyTopo(body, ee2->nodes[1]));
        } else {
          printf("%d,%d ", EG_indexBodyTopo(body, ee2->nodes[1]),
                           EG_indexBodyTopo(body, ee2->nodes[0]));
        }
      }
      printf("\n");
      if (loops[1] != NULL) {
        eloop = (egELoop *) loops[1]->blind;
        printf("     ELoop1\n");
        printf("            %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[0]));
        for (k = 1; k < eloop->eedges.nobjs; k++)
          printf(" %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[k]));
        printf("\n");
        printf("            ");
        for (k = 0; k < eloop->eedges.nobjs; k++) {
          ee2 = (egEEdge *) eloop->eedges.objs[k]->blind;
          if (eloop->senses[k] == SFORWARD) {
            printf("%d,%d ", EG_indexBodyTopo(body, ee2->nodes[0]),
                             EG_indexBodyTopo(body, ee2->nodes[1]));
          } else {
            printf("%d,%d ", EG_indexBodyTopo(body, ee2->nodes[1]),
                             EG_indexBodyTopo(body, ee2->nodes[0]));
          }
        }
        printf("\n");
      }
#endif
    }
  
  EG_free(ncnt);
  EG_free(nodes);
  
  return EGADS_SUCCESS;
}


static void
EG_getEdgeIndices(egTessel *btess, int iface, int iedge, int i1, int i2,
                  int off, egEdgeUV **edgsen)
{
  int m1, pt1, pi1, m2, pt2, pi2, isen;
  
  m1  = pt1 = btess->tess2d[iface].ptype[i1];
  pi1 = btess->tess2d[iface].pindex[i1];
  m2  = pt2 = btess->tess2d[iface].ptype[i2];
  pi2 = btess->tess2d[iface].pindex[i2];
  if ((pt1 == 0) && (pt2 == 0)) {
    if (btess->tess1d[iedge].nodes[0] == pi1) {
      m1 = 1;
      m2 = 2;
    } else {
      m2 = 1;
      m1 = 2;
    }
  } else if (pt1 == 0) {
    if ((btess->tess1d[iedge].nodes[0] == pi1) &&
        (btess->tess1d[iedge].nodes[1] == pi1)) {
      if (pt2 == 2) {
        m1 = 1;
      } else {
        m1 = btess->tess1d[iedge].npts;
      }
    } else if (btess->tess1d[iedge].npts > 3) {
      if (pt2 == 2) m1 = 1;
      if (pt2 == btess->tess1d[iedge].npts-1) m1 = btess->tess1d[iedge].npts;
    } else if (btess->tess1d[iedge].nodes[0] == pi1) {
      m1 = 1;
    } else {
      m1 = btess->tess1d[iedge].npts;
    }
  } else if (pt2 == 0) {
    if ((btess->tess1d[iedge].nodes[0] == pi2) &&
        (btess->tess1d[iedge].nodes[1] == pi2)) {
      if (pt1 == 2) {
        m2 = 1;
      } else {
        m2 = btess->tess1d[iedge].npts;
      }
    } else if (btess->tess1d[iedge].npts > 3) {
      if (pt1 == 2) m2 = 1;
      if (pt1 == btess->tess1d[iedge].npts-1) m2 = btess->tess1d[iedge].npts;
    } else if (btess->tess1d[iedge].nodes[0] == pi2) {
      m2 = 1;
    } else {
      m2 = btess->tess1d[iedge].npts;
    }
  }
  isen = 0;
  if (m1 < m2) isen = 1;
  if (edgsen[3*iedge+isen] == NULL) {
    printf(" EGADS Internal: Index %d sense = %d is NULL!\n", iedge+1, isen);
/*  printf(" pt/pi1 = %d/%d pt/pi2 = %d/%d  ms = %d %d\n", pt1, pi2, pt2, pi2,
           m1, m2);
    printf("      nodes = %d %d\n", btess->tess1d[iedge].nodes[0],
           btess->tess1d[iedge].nodes[1]); */
    return;
  }
  
  /* fill the indices */
  if (edgsen[3*iedge+isen]->iuv[m1-1] == 0) {
    edgsen[3*iedge+isen]->iuv[m1-1] = i1 + off + 1;
  } else {
    if (edgsen[3*iedge+isen]->iuv[m1-1] != i1 + off + 1)
      printf(" EGADS Internal: Index %d/%d Mismatch %d %d  sen = %d\n",
             iedge+1, m1, edgsen[3*iedge+isen]->iuv[m1-1], i1 + off + 1, isen);
  }
  if (edgsen[3*iedge+isen]->iuv[m2-1] == 0) {
    edgsen[3*iedge+isen]->iuv[m2-1] = i2 + off + 1;
  } else {
    if (edgsen[3*iedge+isen]->iuv[m2-1] != i2 + off + 1)
      printf(" EGADS Internal: Index %d/%d Mismatch %d %d  sen = %d\n",
             iedge+1, m2, edgsen[3*iedge+isen]->iuv[m2-1], i2 + off + 1, isen);
  }
}


static void
EG_setLoopArea(void *uvmap, egELoop *eloop)
{
  int      i, j, k, m, ie, eesen, first = 0;
  double   ustart[2], uvlast[2], uv[2], area = 0.0;
  egObject *eobj;
  egEEdge  *eedge;

  ustart[0] = ustart[1] = 0.0;
  uvlast[0] = uvlast[1] = 0.0;
  for (i = 0; i < eloop->eedges.nobjs; i++) {
    eesen = eloop->senses[i];
    eobj  = eloop->eedges.objs[i];
    if (eobj == NULL) continue;
    if (eobj->mtype == DEGENERATE) continue;
    eedge = (egEEdge *) eobj->blind;
    if (eedge == NULL) continue;
    for (j = 0; j < eedge->nsegs; j++) {
      for (k = 0; k < eloop->nedge; k++) {
        ie = eedge->segs[j].iedge - 1;
        if (eedge->sedges[ie].edge == eloop->edgeUVs[k].edge)
          if (eesen*eedge->segs[j].sense == eloop->edgeUVs[k].sense) break;
      }
      if (k == eloop->nedge) {
        printf(" EGADS Internal: Cant find Edge/sense (EG_setLoopArea)!");
        continue;
      }
/*    {
        int      stat, oclass, mtype, n, *senses;
        double   range[2];
        egObject *geom, **nodes;
        stat = EG_getTopology(eloop->edgeUVs[k].edge, &geom, &oclass, &mtype,
                              range, &n, &nodes, &senses);
        if (stat != EGADS_SUCCESS) {
          printf(" EG_getTopology = %d\n", stat);
        } else {
          if (n == 1) {
            printf(" Node  = %lx\n", (long) nodes[0]);
          } else if (eloop->edgeUVs[k].sense == SFORWARD) {
            printf(" Nodes = %lx %lx\n", (long) nodes[0], (long) nodes[1]);
          } else {
            printf(" Nodes = %lx %lx\n", (long) nodes[1], (long) nodes[0]);
          }
        }
      }  */
      if (eloop->edgeUVs[k].sense == SFORWARD) {
        for (m = 0; m < eloop->edgeUVs[k].npts-1; m++) {
          EG_getUVmap(uvmap, eloop->edgeUVs[k].iuv[m], uv);
          if (first == 0) {
            ustart[0] = uv[0];
            ustart[1] = uv[1];
            first++;
          } else {
            area += (uv[0]+uvlast[0])*(uv[1]-uvlast[1]);
          }
          uvlast[0] = uv[0];
          uvlast[1] = uv[1];
        }
      } else {
        for (m = eloop->edgeUVs[k].npts-1; m > 0; m--) {
          EG_getUVmap(uvmap, eloop->edgeUVs[k].iuv[m], uv);
          if (first == 0) {
            ustart[0] = uv[0];
            ustart[1] = uv[1];
            first++;
          } else {
            area += (uv[0]+uvlast[0])*(uv[1]-uvlast[1]);
          }
          uvlast[0] = uv[0];
          uvlast[1] = uv[1];
        }
      }
    }
  }
  area += (ustart[0]+uvlast[0])*(ustart[1]-uvlast[1]);
  eloop->area = area/2.0;
}


static int
EG_coarseTris(int neloop, egELoop **eloops, void *uvmap)
{
  int      i, j, k, m, n, jj, nf8, ie, iloop, ntri, ntot, eesen, last;
  int      *nc, *frame;
  double   *uvs;
  egObject *eobj;
  egEEdge  *eedge;
  egELoop  *eloop;
  fillArea fast;
  
  nc = (int *) EG_alloc(neloop*sizeof(int));
  if (nc == NULL) {
    printf(" EGADS Internal: Cannot Malloc %d ELoops (EG_coarseTris)!", neloop);
    EG_free(nc);
    return EGADS_TOPOERR;
  }
  
  /* count the vertices */
  for (last = ntot = iloop = 0; iloop < neloop; iloop++) {
    eloop = eloops[iloop];
    for (i = 0; i < eloop->eedges.nobjs; i++) {
      eesen = eloop->senses[i];
      eobj  = eloop->eedges.objs[i];
      if (eobj == NULL) continue;
      if (eobj->mtype == DEGENERATE) continue;
      eedge = (egEEdge *) eobj->blind;
      if (eedge == NULL) continue;
      for (j = 0; j < eedge->nsegs; j++) {
        for (k = 0; k < eloop->nedge; k++) {
          ie = eedge->segs[j].iedge - 1;
          if (eedge->sedges[ie].edge == eloop->edgeUVs[k].edge)
            if (eesen*eedge->segs[j].sense == eloop->edgeUVs[k].sense) break;
        }
        if (k == eloop->nedge) {
          printf(" EGADS Internal: Cant find Edge/sense (EG_coarseTris)!");
          EG_free(nc);
          return EGADS_TOPOERR;
        }
        ntot += eloop->edgeUVs[k].npts-1;
      }
    }
    if (iloop == 0) {
      nc[0] = ntot;
    } else {
      nc[iloop] = ntot - last;
    }
    last = ntot;
  }
  ntri = ntot-2 + 2*(neloop-1);
#ifdef DEBUG
  printf(" neloop = %d  ntri = %d   ntot = %d\n", neloop, ntri, ntot);
  for (iloop = 0; iloop < neloop; iloop++)
    printf("     eloop %d: len = %d\n", iloop+1, nc[iloop]);
#endif
  uvs = (double *) EG_alloc((2*ntot+2)*sizeof(double));
  if (uvs == NULL) {
    printf(" EGADS Error: Allocating %d uvs (EG_coarseTris)!\n", ntot);
    EG_free(nc);
    return EGADS_MALLOC;
  }
  uvs[0] = uvs[1] = 0.0;
  
  for (ntot = iloop = 0; iloop < neloop; iloop++) {
    eloop = eloops[iloop];
    for (i = 0; i < eloop->eedges.nobjs; i++) {
      eesen = eloop->senses[i];
      eobj  = eloop->eedges.objs[i];
      if (eobj == NULL) continue;
      if (eobj->mtype == DEGENERATE) continue;
      eedge = (egEEdge *) eobj->blind;
      if (eedge == NULL) continue;
      for (jj = 0; jj < eedge->nsegs; jj++) {
        j = jj;
        if (eesen == SREVERSE) j = eedge->nsegs - jj - 1;
        for (k = 0; k < eloop->nedge; k++) {
          ie = eedge->segs[j].iedge - 1;
          if (eedge->sedges[ie].edge == eloop->edgeUVs[k].edge)
            if (eesen*eedge->segs[j].sense == eloop->edgeUVs[k].sense) break;
        }
        if (k == eloop->nedge) {
          printf(" EGADS Internal: Cant find Edge/Sense (EG_coarseTris)!");
          EG_free(uvs);
          EG_free(nc);
          return EGADS_TOPOERR;
        }
        if (eloop->edgeUVs[k].sense == SFORWARD) {
          for (m = 0; m < eloop->edgeUVs[k].npts-1; m++) {
            ntot++;
            EG_getUVmap(uvmap, eloop->edgeUVs[k].iuv[m], &uvs[2*ntot]);
          }
        } else {
          for (m = eloop->edgeUVs[k].npts-1; m > 0; m--) {
            ntot++;
            EG_getUVmap(uvmap, eloop->edgeUVs[k].iuv[m], &uvs[2*ntot]);
          }
        }
      }
    }
  }

#ifdef DEBUG
  {
    char       filename[25];
    static int filecnt = 0;
    FILE       *fp;
    
    snprintf(filename, 25, "composit%d.frame", filecnt);
    filecnt++;
    fp = fopen(filename, "w");
    if (fp == NULL) {
      printf(" CANNOT OPEN frame.dat!\n");
    } else {
      fprintf(fp, "%d\n", neloop);
      for (i = 0; i < neloop; i++)
        fprintf(fp, "%d\n", nc[i]);
      for (i = 1; i <= ntot; i++)
        fprintf(fp, "%20.13le %20.13le\n", uvs[2*i  ], uvs[2*i+1]);
      for (iloop = 0; iloop < neloop; iloop++) {
        eloop = eloops[iloop];
        for (i = 0; i < eloop->eedges.nobjs; i++) {
          eesen = eloop->senses[i];
          eobj  = eloop->eedges.objs[i];
          if (eobj == NULL) continue;
          if (eobj->mtype == DEGENERATE) continue;
          eedge = (egEEdge *) eobj->blind;
          if (eedge == NULL) continue;
          for (jj = 0; jj < eedge->nsegs; jj++) {
            j = jj;
            if (eesen == SREVERSE) j = eedge->nsegs - jj - 1;
            for (k = 0; k < eloop->nedge; k++) {
              ie = eedge->segs[j].iedge - 1;
              if (eedge->sedges[ie].edge == eloop->edgeUVs[k].edge)
                if (eesen*eedge->segs[j].sense == eloop->edgeUVs[k].sense) break;
            }
            if (k == eloop->nedge) {
              printf(" EGADS Internal: Cant Find Edge/Sense (EG_coarseTris)!");
              k--;
            }
            for (m = 0; m < eloop->edgeUVs[k].npts-1; m++)
              fprintf(fp, "%d\n", k+1);
          }
        }
      }
      fclose(fp);
    }
  }
#endif

  /* triangulate */
  fast.pts   = NULL;
  fast.segs  = NULL;
  fast.front = NULL;
  frame      = (int *) EG_alloc(3*ntri*sizeof(int));
  if (frame == NULL) {
    printf(" EGADS Error: Allocating %d frame (EG_coarseTris)!\n", ntri);
    EG_free(uvs);
    EG_free(nc);
    return EGADS_MALLOC;
  }
  n = EG_fillArea(-neloop, nc, uvs, frame, &nf8, 0, &fast);
  if (fast.segs  != NULL) EG_free(fast.segs);
  if (fast.pts   != NULL) EG_free(fast.pts);
  if (fast.front != NULL) EG_free(fast.front);
  EG_free(frame);
  EG_free(uvs);
  EG_free(nc);
  
  if ((nf8 != 0) || (n != ntri)) {
    printf(" EGADS Error: Cannot Triangulate Region %d %d (EG_coarseTris)!\n",
           n, nf8);
    return EGADS_TOPOERR;
  }
  
  return EGADS_SUCCESS;
}


int
EG_writeEBody(egObject *EBody, FILE *fp)
{
  int      i, j, k, stat, nds[2], nattr;
  egObject *body, *obj;
  egAttrs  *attrs;
  egEBody  *ebody;
  egEShell *eshell;
  egEFace  *eface;
  egELoop  *eloop;
  egEEdge  *eedge;
  
  if (EBody == NULL)               return EGADS_NULLOBJ;
  if (EBody->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EBody->oclass != EBODY)      return EGADS_NOTTESS;
  if (EG_sameThread(EBody))        return EGADS_CNTXTHRD;
  
  ebody = (egEBody *) EBody->blind;
  if (ebody == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_writeEBody)!\n");
    return EGADS_NOTFOUND;
  }
  if (ebody->done == 0) {
    printf(" EGADS Error: EBody not finialized (EG_writeEBody)!\n");
    return EGADS_EFFCTOBJ;
  }
  body  = ebody->ref;
  nattr = 0;
  attrs = NULL;
  if (EBody->attrs != NULL) {
    attrs = (egAttrs *) EBody->attrs;
    nattr = EG_writeNumAttr(attrs);
  }
  
  fprintf(fp, "%d %d %d %d %d %d %lf\n", ebody->eedges.nobjs,
          ebody->eloops.nobjs, ebody->efaces.nobjs, ebody->eshells.nobjs, nattr,
          ebody->nedge, ebody->angle);
  if (body->mtype == SOLIDBODY) {
    for (i = 0; i < ebody->eshells.nobjs; i++) {
      fprintf(fp, "%d ", ebody->senses[i]);
      if (((i+1)%20 == 0) && (i != ebody->eshells.nobjs-1)) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  if ((nattr != 0) && (attrs != NULL)) EG_writeAttr(attrs, fp);
  
  /* source Edges */
  for (j = 0; j < ebody->nedge; j++) {
    nds[0] = EG_indexBodyTopo(body, ebody->edges[j].edge);
    if (nds[0] <= EGADS_SUCCESS) {
      printf(" EGADS Error: Source Edge = %d (EG_writeEBody)!\n", nds[0]);
      return EGADS_TOPOERR;
    }
    fprintf(fp, "%d %d %d\n", nds[0], ebody->edges[j].curve,
            ebody->edges[j].npts);
    fprintf(fp, "%19.12le %19.12le %19.12le\n", ebody->edges[j].dstart[0],
            ebody->edges[j].dstart[1], ebody->edges[j].dstart[2]);
    fprintf(fp, "%19.12le %19.12le %19.12le\n", ebody->edges[j].dend[0],
            ebody->edges[j].dend[1], ebody->edges[j].dend[2]);
    for (k = 0; k < ebody->edges[j].npts; k++) {
      fprintf(fp, "%19.12le ", ebody->edges[j].ts[k]);
      if (((k+1)%5 == 0) && (k != ebody->edges[j].npts-1)) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  
  /* EEdges */
  for (i = 0; i < ebody->eedges.nobjs; i++) {
    obj     = ebody->eedges.objs[i];
    nattr   = 0;
    attrs   = NULL;
    if (obj->attrs != NULL) {
      attrs = (egAttrs *) obj->attrs;
      nattr = EG_writeNumAttr(attrs);
    }
    eedge   = (egEEdge *) obj->blind;
    nds[0]  = EG_indexBodyTopo(body, eedge->nodes[0]);
    nds[1]  = EG_indexBodyTopo(body, eedge->nodes[1]);
    if ((nds[0] <= EGADS_SUCCESS) || (nds[1] <= EGADS_SUCCESS)) {
      printf(" EGADS Error: Node indices = %d %d (EG_writeEBody)!\n",
             nds[0], nds[1]);
      return EGADS_TOPOERR;
    }
    fprintf(fp, "%hd %d %d %d %d %19.12le %19.12le\n", obj->mtype, eedge->nsegs,
            nds[0], nds[1], nattr, eedge->trange[0], eedge->trange[1]);
    for (j = 0; j < eedge->nsegs; j++) {
      nds[0] = 0;
      if (eedge->segs[j].nstart != NULL)
        nds[0] = EG_indexBodyTopo(body, eedge->segs[j].nstart);
      if (nds[0] < EGADS_SUCCESS) {
        printf(" EGADS Error: Node start = %d (EG_writeEBody)!\n", nds[0]);
        return EGADS_TOPOERR;
      }
      fprintf(fp, "%d %d %d\n", eedge->segs[j].iedge, eedge->segs[j].sense,
              nds[0]);
      fprintf(fp, "%19.12le %19.12le\n", eedge->segs[j].tstart,
              eedge->segs[j].tend);
    }
    if ((nattr != 0) && (attrs != NULL)) EG_writeAttr(attrs, fp);
  }
  
  /* ELoops */
  for (i = 0; i < ebody->eloops.nobjs; i++) {
    obj     = ebody->eloops.objs[i];
    nattr   = 0;
    attrs   = NULL;
    if (obj->attrs != NULL) {
      attrs = (egAttrs *) obj->attrs;
      nattr = EG_writeNumAttr(attrs);
    }
    eloop   = (egELoop *) obj->blind;
    if (eloop == NULL) {
      fprintf(fp, "%hd %d %d %d %19.12le\n", obj->mtype, 0, 0, 0, 0.0);
      continue;
    }
    fprintf(fp, "%hd %d %d %d %19.12le\n", obj->mtype, eloop->eedges.nobjs,
            eloop->nedge, nattr, eloop->area);
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      k = EG_indexBodyTopo(EBody, eloop->eedges.objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: EEdge %d index = %d in ELoop %d (EG_writeEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      fprintf(fp, "%d ", k);
      if (((j+1)%20 == 0) && (j != eloop->eedges.nobjs-1)) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      fprintf(fp, "%d ", eloop->senses[j]);
      if (((j+1)%20 == 0) && (j != eloop->eedges.nobjs-1)) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j = 0; j < eloop->nedge; j++) {
      k = EG_indexBodyTopo(body, eloop->edgeUVs[j].edge);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Edge %d index = %d in ELoop %d (EG_writeEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      fprintf(fp, "%d %d %d\n", k, eloop->edgeUVs[j].sense,
              eloop->edgeUVs[j].npts);
      for (k = 0; k < eloop->edgeUVs[j].npts; k++) {
        fprintf(fp, "%d ", eloop->edgeUVs[j].iuv[k]);
        if (((k+1)%20 == 0) && (k != eloop->edgeUVs[j].npts-1))
          fprintf(fp, "\n");
      }
      fprintf(fp, "\n");
    }
    if ((nattr != 0) && (attrs != NULL)) EG_writeAttr(attrs, fp);
  }
  
  /* EFaces */
  for (i = 0; i < ebody->efaces.nobjs; i++) {
    obj     = ebody->efaces.objs[i];
    nattr   = 0;
    attrs   = NULL;
    if (obj->attrs != NULL) {
      attrs = (egAttrs *) obj->attrs;
      nattr = EG_writeNumAttr(attrs);
    }
    eface   = (egEFace *) obj->blind;
    fprintf(fp, "%hd %d %d %d %d\n", obj->mtype, eface->npatch,
            eface->eloops.nobjs, eface->last, nattr);
    if (eface->npatch != 1) {
      stat = EG_uvmapWrite(eface->uvmap, eface->trmap, fp);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EFace %d  uvmapWrite = %d (EG_writeEBody)!\n",
               i+1, stat);
        return stat;
      }
    } else {
      fprintf(fp, "%19.12le %19.12le %19.12le %19.12le\n", eface->range[0],
              eface->range[1], eface->range[2], eface->range[3]);
    }
    for (j = 0; j < eface->eloops.nobjs; j++) {
      k = EG_indexBodyTopo(EBody, eface->eloops.objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Loop %d index = %d in EFace %d (EG_writeEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      fprintf(fp, "%d ", k);
      if (((j+1)%20 == 0) && (j != eface->eloops.nobjs-1)) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j = 0; j < eface->eloops.nobjs; j++) {
      fprintf(fp, "%d ", eface->senses[j]);
      if (((j+1)%20 == 0) && (j != eface->eloops.nobjs-1)) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j = 0; j < eface->npatch; j++) {
      k = EG_indexBodyTopo(body, eface->patches[j].face);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d index = %d in EFace %d (EG_writeEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      fprintf(fp, "%d %d %d %d %d\n", k, eface->patches[j].start,
              eface->patches[j].nuvs,  eface->patches[j].ndeflect,
              eface->patches[j].ntris);
      for (k = 0; k < eface->patches[j].ntris; k++)
        fprintf(fp, "%d %d %d\n", eface->patches[j].uvtris[3*k  ],
                                  eface->patches[j].uvtris[3*k+1],
                                  eface->patches[j].uvtris[3*k+2]);
      for (k = 0; k < eface->patches[j].nuvs; k++)
        fprintf(fp, "%19.12le %19.12le\n", eface->patches[j].uvs[2*k  ],
                                           eface->patches[j].uvs[2*k+1]);
      for (k = 0; k < eface->patches[j].ndeflect; k++)
        fprintf(fp, "%19.12le %19.12le %19.12le\n",
                eface->patches[j].deflect[3*k  ],
                eface->patches[j].deflect[3*k+1],
                eface->patches[j].deflect[3*k+2]);
    }
    if ((nattr != 0) && (attrs != NULL)) EG_writeAttr(attrs, fp);
  }
  
  /* EShells */
  for (i = 0; i < ebody->eshells.nobjs; i++) {
    obj     = ebody->eshells.objs[i];
    nattr   = 0;
    attrs   = NULL;
    if (obj->attrs != NULL) {
      attrs = (egAttrs *) obj->attrs;
      nattr = EG_writeNumAttr(attrs);
    }
    eshell  = (egEShell *) obj->blind;
    fprintf(fp, "%hd %d %d\n", obj->mtype, eshell->efaces.nobjs, nattr);
    for (j = 0; j < eshell->efaces.nobjs; j++) {
      k = EG_indexBodyTopo(EBody, eshell->efaces.objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Face %d index = %d in EFace %d (EG_writeEBody)!\n",
               k, j+1, i+1);
        return EGADS_TOPOERR;
      }
      fprintf(fp, "%d ", k);
      if (((j+1)%20 == 0) && (j != eshell->efaces.nobjs-1)) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    if ((nattr != 0) && (attrs != NULL)) EG_writeAttr(attrs, fp);
  }
  
  return EGADS_SUCCESS;
}


int
EG_readEBody(FILE *fp, egObject *body, egObject **EBody)
{
  int      i, j, k, n, stat, nattr, nds[2];
  egObject *context, *eobj, *tobj;
  egEBody  *ebody;
  egEShell *eshell;
  egEFace  *eface;
  egELoop  *eloop;
  egEEdge  *eedge;
  
  *EBody = NULL;
  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->oclass != BODY)       return EGADS_NOTTESS;
  if (body->blind == NULL)        return EGADS_NODATA;
  if (EG_sameThread(body))        return EGADS_CNTXTHRD;
  context = EG_context(body);
  
  ebody = (egEBody *) EG_alloc(sizeof(egEBody));
  if (ebody == NULL) {
    printf(" EGADS Error: Malloc of EBody  (EG_readEBody)!\n");
    return EGADS_MALLOC;
  }
  ebody->ref           = body;
  ebody->eedges.objs   = NULL;
  ebody->eloops.objs   = NULL;
  ebody->efaces.objs   = NULL;
  ebody->eshells.objs  = NULL;
  ebody->senses        = NULL;
  ebody->eedges.nobjs  = 0;
  ebody->eloops.nobjs  = 0;
  ebody->efaces.nobjs  = 0;
  ebody->eshells.nobjs = 0;
  ebody->angle         = 0.0;
  ebody->done          = 1;
  ebody->nedge         = 0;
  ebody->edges         = NULL;
  stat = EG_makeObject(context, &eobj);
  if (stat != EGADS_SUCCESS) {
    EG_free(ebody);
    printf(" EGADS Error: Cannot make EBody Object (EG_readEBody)!\n");
    return stat;
  }
  eobj->oclass = EBODY;
  eobj->mtype  = body->mtype;
  eobj->blind  = ebody;
  EG_referenceObject(eobj, context);
  EG_referenceTopObj(body, eobj);
  
  n = fscanf(fp, "%d %d %d %d %d %d %lf", &ebody->eedges.nobjs,
             &ebody->eloops.nobjs, &ebody->efaces.nobjs, &ebody->eshells.nobjs,
             &nattr, &ebody->nedge, &ebody->angle);
  if (n != 7) goto readerr;
  if (body->mtype == SOLIDBODY) {
    ebody->senses = (int *) EG_alloc(ebody->eshells.nobjs*sizeof(int));
    if (ebody->senses == NULL) {
      printf(" EGADS Error: Malloc on %d Shell senses (EG_readEBody)!\n",
             ebody->eshells.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (i = 0; i < ebody->eshells.nobjs; i++) {
      n = fscanf(fp, "%d", &ebody->senses[i]);
      if (n != 1) goto readerr;
    }
  }
  if (nattr != 0) EG_readAttrs(eobj, nattr, fp);
  
  /* generate source Edge discretization */
  ebody->edges = (egEdVert *) EG_alloc(ebody->nedge*sizeof(egEdVert));
  if (ebody->edges == NULL) {
    printf(" EGADS Error: Malloc on %d Edges (EG_readEBody)!\n", ebody->nedge);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->nedge; i++) {
    ebody->edges[i].edge = NULL;
    ebody->edges[i].npts = 0;
    ebody->edges[i].ts   = NULL;
  }
  for (j = 0; j < ebody->nedge; j++) {
    n = fscanf(fp, "%d %d %d", &nds[0], &ebody->edges[j].curve,
               &ebody->edges[j].npts);
    if (n != 3) goto readerr;
    stat = EG_objectBodyTopo(body, EDGE, nds[0], &ebody->edges[j].edge);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Object on Source Edge %d (EG_readEBody)!\n", i+1);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    n = fscanf(fp, "%le %le %le", &ebody->edges[j].dstart[0],
               &ebody->edges[j].dstart[1], &ebody->edges[j].dstart[2]);
    if (n != 3) goto readerr;
    n = fscanf(fp, "%le %le %le", &ebody->edges[j].dend[0],
               &ebody->edges[j].dend[1], &ebody->edges[j].dend[2]);
    if (n != 3) goto readerr;
    ebody->edges[j].ts = (double *) EG_alloc(ebody->edges[j].npts*
                                             sizeof(double));
    if (ebody->edges[j].ts == NULL) {
      printf(" EGADS Error: Malloc on %d ts source Edge %d (EG_readEBody)!\n",
             ebody->edges[j].npts, j+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (k = 0; k < ebody->edges[j].npts; k++) {
      n = fscanf(fp, "%le", &ebody->edges[j].ts[k]);
      if (n != 1) goto readerr;
    }
  }
  
  /* populate the EEdges */
  ebody->eedges.objs  = (egObject **) EG_alloc(ebody->eedges.nobjs*
                                               sizeof(egObject *));
  if (ebody->eedges.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EEdge Objects (EG_readEBody)!\n",
           ebody->eedges.nobjs);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eedges.nobjs; i++) ebody->eedges.objs[i] = NULL;
  for (i = 0; i < ebody->eedges.nobjs; i++) {
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EEdge %d/%d Object (EG_readEBody)!\n",
             i+1, ebody->eedges.nobjs);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    tobj->oclass          = EEDGE;
    tobj->mtype           = TWONODE;
    tobj->topObj          = eobj;
    ebody->eedges.objs[i] = tobj;
    eedge = (egEEdge *) EG_alloc(sizeof(egEEdge));
    if (eedge == NULL) {
      printf(" EGADS Error: Malloc on %d EEdge blind (EG_readEBody)!\n", i+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    tobj->blind = eedge;
    n = fscanf(fp, "%hd %d %d %d %d %le %le\n", &tobj->mtype, &eedge->nsegs,
               &nds[0], &nds[1], &nattr, &eedge->trange[0], &eedge->trange[1]);
    if (n != 7) goto readerr;
    stat = EG_objectBodyTopo(body, NODE, nds[0], &eedge->nodes[0]);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Object on EEdge %d Node 0 (EG_readEBody)!\n",
             i+1);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    stat = EG_objectBodyTopo(body, NODE, nds[1], &eedge->nodes[1]);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Object on EEdge %d Node 1 (EG_readEBody)!\n",
             i+1);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    eedge->sedges = ebody->edges;
    eedge->segs   = (egEEseg *) EG_alloc(eedge->nsegs*sizeof(egEEseg));
    if (eedge->segs == NULL) {
      printf(" EGADS Error: Malloc EEdge %d nsegs = %d (EG_readEBody)!\n",
             i+1, eedge->nsegs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eedge->nsegs; j++)
      eedge->segs[j].nstart  = NULL;
    for (j = 0; j < eedge->nsegs; j++) {
      n = fscanf(fp, "%d %d %d", &eedge->segs[j].iedge, &eedge->segs[j].sense,
                 &nds[0]);
      if (n != 3) goto readerr;
      if (nds[0] != 0) {
        stat = EG_objectBodyTopo(body, NODE, nds[0], &eedge->segs[j].nstart);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Object on EEdge %d Node 0 (EG_readEBody)!\n",
                 i+1);
          EG_destroyEBody(eobj, 1);
          return stat;
        }
      }
      n = fscanf(fp, "%le %le", &eedge->segs[j].tstart, &eedge->segs[j].tend);
      if (n != 2) goto readerr;
    }
    if (nattr != 0) EG_readAttrs(tobj, nattr, fp);
  }
  
  /* populate the ELoops */
  ebody->eloops.objs  = (egObject **) EG_alloc(ebody->eloops.nobjs*
                                               sizeof(egObject *));
  if (ebody->eloops.objs == NULL) {
    printf(" EGADS Error: Malloc on %d ELoop Objects (EG_readEBody)!\n",
           ebody->eloops.nobjs);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eloops.nobjs; i++) ebody->eloops.objs[i] = NULL;
  for (i = 0; i < ebody->eloops.nobjs; i++) {
    eloop = (egELoop *) EG_alloc(sizeof(egELoop));
    if (eloop == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop blind (EG_readEBody)!\n", i+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make ELoop %d Object (EG_readEBody)!\n",
             i+1);
      EG_free(eloop);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    ebody->eloops.objs[i] = tobj;
    tobj->oclass          = ELOOPX;
/*@-kepttrans@*/
    tobj->topObj          = eobj;
/*@+kepttrans@*/
    tobj->blind           = eloop;
    eloop->eedges.nobjs   = 0;
    eloop->eedges.objs    = NULL;
    eloop->edgeUVs        = NULL;
    eloop->senses         = NULL;
    n = fscanf(fp, "%hd %d %d %d %le\n", &tobj->mtype, &eloop->eedges.nobjs,
               &eloop->nedge, &nattr, &eloop->area);
    if (n != 5) goto readerr;
    if (eloop->eedges.nobjs == 0) {
      EG_free(tobj->blind);
      tobj->blind = NULL;
      continue;
    }
    eloop->eedges.objs = (egObject **) EG_alloc(eloop->eedges.nobjs*
                                                sizeof(egObject *));
    if (eloop->eedges.objs == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop %d Objects (EG_readEBody)!\n",
             i+1, eloop->eedges.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eloop->eedges.nobjs; j++) eloop->eedges.objs[j] = NULL;
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      n = fscanf(fp, "%d", &k);
      if (n != 1) goto readerr;
      eloop->eedges.objs[j] = ebody->eedges.objs[k-1];
    }
    eloop->senses = (int *) EG_alloc(eloop->eedges.nobjs*sizeof(int));
    if (eloop->senses == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop %d senses (EG_readEBody)!\n",
             i+1, eloop->nedge);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      n = fscanf(fp, "%d", &eloop->senses[j]);
      if (n != 1) goto readerr;
    }
    eloop->edgeUVs = (egEdgeUV *) EG_alloc(eloop->nedge*sizeof(egEdgeUV));
    if (eloop->edgeUVs == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop %d edgeUVs (EG_readEBody)!\n",
             i+1, eloop->nedge);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eloop->nedge; j++) {
      eloop->edgeUVs[j].npts  = 0;
      eloop->edgeUVs[j].sense = 0;
      eloop->edgeUVs[j].edge  = NULL;
      eloop->edgeUVs[j].iuv   = NULL;
    }
    for (j = 0; j < eloop->nedge; j++) {
      n = fscanf(fp, "%d %d %d", &k, &eloop->edgeUVs[j].sense,
                 &eloop->edgeUVs[j].npts);
      if (n != 3) goto readerr;
      stat = EG_objectBodyTopo(body, EDGE, k, &eloop->edgeUVs[j].edge);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Object on ELoop %d Edge (EG_readEBody)!\n",
               i+1);
        EG_destroyEBody(eobj, 1);
        return stat;
      }
      eloop->edgeUVs[j].iuv = (int *) EG_alloc(eloop->edgeUVs[j].npts*
                                               sizeof(int));
      if (eloop->edgeUVs == NULL) {
        printf(" EGADS Error: Malloc on %d ELoop %d iUVs %d (EG_readEBody)!\n",
               i+1, j+1, eloop->edgeUVs[j].npts);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < eloop->edgeUVs[j].npts; k++) {
        n = fscanf(fp, "%d", &eloop->edgeUVs[j].iuv[k]);
        if (n != 1) goto readerr;
      }
    }
    if (nattr != 0) EG_readAttrs(tobj, nattr, fp);
  }
  
  /* populate the EFaces */
  ebody->efaces.objs  = (egObject **) EG_alloc(ebody->efaces.nobjs*
                                               sizeof(egObject *));
  if (ebody->efaces.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EFace Objects (EG_readEBody)!\n",
           ebody->efaces.nobjs);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->efaces.nobjs; i++) ebody->efaces.objs[i] = NULL;
  for (i = 0; i < ebody->efaces.nobjs; i++) {
    eface = (egEFace *) EG_alloc(sizeof(egEFace));
    if (eface == NULL) {
      printf(" EGADS Error: Malloc on %d EFace blind (EG_readEBody)!\n", i+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EFace %d/%d Object (EG_readEBody)!\n",
             i+1, ebody->efaces.nobjs);
      EG_free(eface);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    ebody->efaces.objs[i] = tobj;
    tobj->oclass          = EFACE;
/*@-kepttrans@*/
    tobj->topObj          = eobj;
/*@+kepttrans@*/
    tobj->blind           = eface;
    eface->trmap          = NULL;
    eface->uvmap          = NULL;
    eface->npatch         = 0;
    eface->patches        = NULL;
    eface->eloops.nobjs   = 0;
    eface->eloops.objs    = NULL;
    eface->senses         = NULL;
    eface->patches        = NULL;
    eface->sedges         = ebody->edges;
    n = fscanf(fp, "%hd %d %d %d %d", &tobj->mtype, &eface->npatch,
               &eface->eloops.nobjs, &eface->last, &nattr);
    if (n != 5) goto readerr;
    if (eface->npatch != 1) {
      stat = EG_uvmapRead(fp, eface->range, &eface->uvmap, &eface->trmap);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EFace %d  uvmapRead = %d (EG_readEBody)!\n",
               i+1, stat);
        EG_destroyEBody(eobj, 1);
        return stat;
      }
    } else {
      n = fscanf(fp, "%le %le %le %le\n", &eface->range[0],
                 &eface->range[1], &eface->range[2], &eface->range[3]);
      if (n != 4) goto readerr;
    }
    eface->eloops.objs = (egObject **) EG_alloc(eface->eloops.nobjs*
                                                sizeof(egObject *));
    if (eface->eloops.objs == NULL) {
      printf(" EGADS Error: Malloc on %d EFace %d Objects (EG_readEBody)!\n",
             i+1, eface->eloops.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eface->eloops.nobjs; j++) eface->eloops.objs[j] = NULL;
    for (j = 0; j < eface->eloops.nobjs; j++) {
      n = fscanf(fp, "%d", &k);
      if (n != 1) goto readerr;
      eface->eloops.objs[j] = ebody->eloops.objs[k-1];
    }
    eface->senses = (int *) EG_alloc(eface->eloops.nobjs*sizeof(int));
    if (eface->senses == NULL) {
      printf(" EGADS Error: Malloc on %d EFace senses %d (EG_readEBody)!\n",
             i+1, eface->eloops.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eface->eloops.nobjs; j++) {
      n = fscanf(fp, "%d", &eface->senses[j]);
      if (n != 1) goto readerr;
    }
    eface->patches = (egEPatch *) EG_alloc(eface->npatch*sizeof(egEPatch));
    if (eface->patches == NULL) {
      printf(" EGADS Error: Malloc on %d Patch Object %d (EG_readEBody)!\n",
             i+1, eface->npatch);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eface->npatch; j++) {
      eface->patches[j].start    = 0;
      eface->patches[j].nuvs     = 0;
      eface->patches[j].ndeflect = 0;
      eface->patches[j].ntris    = 0;
      eface->patches[j].tol      = -1.0;
      eface->patches[j].uvtris   = NULL;
      eface->patches[j].uvtric   = NULL;
      eface->patches[j].uvs      = NULL;
      eface->patches[j].deflect  = NULL;
      eface->patches[j].face     = NULL;
    }
    for (j = 0; j < eface->npatch; j++) {
      n = fscanf(fp, "%d %d %d %d %d", &nds[0], &eface->patches[j].start,
                 &eface->patches[j].nuvs,  &eface->patches[j].ndeflect,
                 &eface->patches[j].ntris);
      if (n != 5) goto readerr;
      stat = EG_objectBodyTopo(body, FACE, nds[0], &eface->patches[j].face);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: Object on EFace %d Face (EG_readEBody)!\n",
               i+1);
        EG_destroyEBody(eobj, 1);
        return stat;
      }
      eface->patches[j].uvtris = (int *) EG_alloc(3*eface->patches[j].ntris*
                                                  sizeof(int));
      if (eface->patches[j].uvtris == NULL) {
        printf(" EGADS Error: Malloc on %d Patch uvtris %d (EG_readEBody)!\n",
               i+1, eface->patches[j].ntris);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < eface->patches[j].ntris; k++) {
        n = fscanf(fp, "%d %d %d", &eface->patches[j].uvtris[3*k  ],
                                   &eface->patches[j].uvtris[3*k+1],
                                   &eface->patches[j].uvtris[3*k+2]);
        if (n != 3) goto readerr;
      }
      eface->patches[j].uvs = (double *) EG_alloc(2*eface->patches[j].nuvs*
                                                  sizeof(double));
      if (eface->patches[j].uvs == NULL) {
        printf(" EGADS Error: Malloc on %d Patch uvs %d (EG_readEBody)!\n",
               i+1, eface->patches[j].nuvs);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < eface->patches[j].nuvs; k++) {
        n = fscanf(fp, "%le %le", &eface->patches[j].uvs[2*k  ],
                                  &eface->patches[j].uvs[2*k+1]);
        if (n != 2) goto readerr;
      }
      eface->patches[j].deflect = (double *)
                          EG_alloc(3*eface->patches[j].ndeflect*sizeof(double));
      if (eface->patches[j].deflect == NULL) {
        printf(" EGADS Error: Malloc on %d Patch deflect %d (EG_readEBody)!\n",
               i+1, eface->patches[j].ndeflect);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < eface->patches[j].ndeflect; k++) {
        n = fscanf(fp, "%le %le %le\n", &eface->patches[j].deflect[3*k  ],
                                        &eface->patches[j].deflect[3*k+1],
                                        &eface->patches[j].deflect[3*k+2]);
        if (n != 3) goto readerr;
      }
    }
    if (nattr != 0) EG_readAttrs(tobj, nattr, fp);
    if (eface->npatch == 1) {
      stat = EG_effectNeighbor(eface);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: %d Patch EG_effectNeighbor = %d (EG_readEBody)!\n",
               i+1, stat);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
    }
  }

  /* populate the EShells */
  ebody->eshells.objs  = (egObject **) EG_alloc(ebody->eshells.nobjs*
                                                sizeof(egObject *));
  if (ebody->eshells.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EShell Objects (EG_readEBody)!\n",
           ebody->eshells.nobjs);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eshells.nobjs; i++) ebody->eshells.objs[i] = NULL;
  for (i = 0; i < ebody->eshells.nobjs; i++) {
    eshell = (egEShell *) EG_alloc(sizeof(egEShell));
    if (eshell == NULL) {
      printf(" EGADS Error: Malloc on %d EShell blind (EG_readEBody)!\n", i+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EFace %d/%d Object (EG_readEBody)!\n",
             i+1, ebody->efaces.nobjs);
      EG_free(eshell);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    tobj->oclass           = ESHELL;
/*@-kepttrans@*/
    tobj->topObj           = eobj;
/*@+kepttrans@*/
    ebody->eshells.objs[i] = tobj;
    tobj->blind            = eshell;
    eshell->efaces.nobjs   = 0;
    eshell->efaces.objs    = NULL;
    n = fscanf(fp, "%hd %d %d\n", &tobj->mtype, &eshell->efaces.nobjs, &nattr);
    if (n != 3) goto readerr;
    eshell->efaces.objs    = (egObject **) EG_alloc(eshell->efaces.nobjs*
                                                    sizeof(egObject *));
    if (eshell->efaces.objs == NULL) {
      printf(" EGADS Error: Malloc on %dEShell %d Objects (EG_readEBody)!\n",
             i+1, eshell->efaces.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eshell->efaces.nobjs; j++) eshell->efaces.objs[j] = NULL;
    for (j = 0; j < eshell->efaces.nobjs; j++) {
      n = fscanf(fp, "%d", &k);
      if (n != 1) goto readerr;
      eshell->efaces.objs[j] = ebody->efaces.objs[k-1];
    }
    if (nattr != 0) EG_readAttrs(tobj, nattr, fp);
  }

  *EBody = eobj;
  return EGADS_SUCCESS;
  
readerr:
  EG_destroyEBody(eobj, 1);
  return EGADS_READERR;
}


int
EG_copyEBody(const egObject *body, egObject **EBody)
{
  int      i, j, k, stat;
  egObject *context, *eobj, *tobj, *sobj;
  egEBody  *ebody,  *src;
  egEShell *eshell, *sshell;
  egEFace  *eface,  *sface;
  egELoop  *eloop,  *sloop;
  egEEdge  *eedge,  *sedge;
  
  *EBody = NULL;
  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->oclass != EBODY)      return EGADS_NOTTESS;
  if (body->blind == NULL)        return EGADS_NODATA;
  if (EG_sameThread(body))        return EGADS_CNTXTHRD;
  context = EG_context(body);
  
  src   = (egEBody *) body->blind;
  ebody = (egEBody *) EG_alloc(sizeof(egEBody));
  if (ebody == NULL) {
    printf(" EGADS Error: Malloc of EBody  (EG_copyEBody)!\n");
    return EGADS_MALLOC;
  }
  ebody->ref           = src->ref;
  ebody->eedges.objs   = NULL;
  ebody->eloops.objs   = NULL;
  ebody->efaces.objs   = NULL;
  ebody->eshells.objs  = NULL;
  ebody->senses        = NULL;
  ebody->eedges.nobjs  = src->eedges.nobjs;
  ebody->eloops.nobjs  = src->eloops.nobjs;
  ebody->efaces.nobjs  = src->efaces.nobjs;
  ebody->eshells.nobjs = src->eshells.nobjs;
  ebody->angle         = src->angle;
  ebody->done          = src->done;
  ebody->nedge         = src->nedge;
  ebody->edges         = NULL;
  stat = EG_makeObject(context, &eobj);
  if (stat != EGADS_SUCCESS) {
    EG_free(ebody);
    printf(" EGADS Error: Cannot make EBody Object (EG_copyEBody)!\n");
    return stat;
  }
  eobj->oclass = EBODY;
  eobj->mtype  = body->mtype;
  eobj->blind  = ebody;
  EG_referenceObject(eobj, context);
  EG_referenceTopObj(src->ref, eobj);

  if (body->mtype == SOLIDBODY) {
    ebody->senses = (int *) EG_alloc(ebody->eshells.nobjs*sizeof(int));
    if (ebody->senses == NULL) {
      printf(" EGADS Error: Malloc on %d Shell senses (EG_copyEBody)!\n",
             ebody->eshells.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (i = 0; i < ebody->eshells.nobjs; i++)
      ebody->senses[i] = src->senses[i];
  }
  EG_attributeDup(body, eobj);
  
  /* populate the source Edges */
  ebody->edges = (egEdVert *) EG_alloc(ebody->nedge*sizeof(egEdVert));
  if (ebody->edges == NULL) {
    printf(" EGADS Error: Malloc on %d Edges (EG_copyEBody)!\n", ebody->nedge);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->nedge; i++) {
    ebody->edges[i].edge = NULL;
    ebody->edges[i].npts = 0;
    ebody->edges[i].ts   = NULL;
  }
  for (j = 0; j < ebody->nedge; j++) {
    ebody->edges[j].edge      = src->edges[j].edge;
    ebody->edges[j].curve     = src->edges[j].curve;
    ebody->edges[j].npts      = src->edges[j].npts;
    ebody->edges[j].dstart[0] = src->edges[j].dstart[0];
    ebody->edges[j].dstart[1] = src->edges[j].dstart[1];
    ebody->edges[j].dstart[2] = src->edges[j].dstart[2];
    ebody->edges[j].dend[0]   = src->edges[j].dend[0];
    ebody->edges[j].dend[1]   = src->edges[j].dend[1];
    ebody->edges[j].dend[2]   = src->edges[j].dend[2];
    ebody->edges[j].ts        = (double *) EG_alloc(ebody->edges[j].npts*
                                                    sizeof(double));
    if (ebody->edges[j].ts == NULL) {
      printf(" EGADS Error: Malloc on %d ts source Edge %d (EG_copyEBody)!\n",
             ebody->edges[j].npts, j+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (k = 0; k < ebody->edges[j].npts; k++)
      ebody->edges[j].ts[k] = src->edges[j].ts[k];
  }
  
  /* populate the EEdges */
  ebody->eedges.objs = (egObject **) EG_alloc(ebody->eedges.nobjs*
                                              sizeof(egObject *));
  if (ebody->eedges.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EEdge Objects (EG_copyEBody)!\n",
           ebody->eedges.nobjs);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eedges.nobjs; i++) ebody->eedges.objs[i] = NULL;
  for (i = 0; i < ebody->eedges.nobjs; i++) {
    sobj  = src->eedges.objs[i];
    sedge = (egEEdge *) sobj->blind;
    stat  = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EEdge %d/%d Object (EG_copyEBody)!\n",
             i+1, ebody->eedges.nobjs);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    tobj->oclass          = sobj->oclass;
    tobj->mtype           = sobj->mtype;
    tobj->topObj          = eobj;
    ebody->eedges.objs[i] = tobj;
    eedge = (egEEdge *) EG_alloc(sizeof(egEEdge));
    if (eedge == NULL) {
      printf(" EGADS Error: Malloc on %d EEdge blind (EG_copyEBody)!\n", i+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    tobj->blind      = eedge;
    eedge->sedges    = ebody->edges;
    eedge->trange[0] = sedge->trange[0];
    eedge->trange[1] = sedge->trange[1];
    eedge->nodes[0]  = sedge->nodes[0];
    eedge->nodes[1]  = sedge->nodes[1];
    eedge->nsegs     = sedge->nsegs;
    eedge->segs      = (egEEseg *) EG_alloc(eedge->nsegs*sizeof(egEEseg));
    if (eedge->segs == NULL) {
      printf(" EGADS Error: Malloc EEdge %d nsegs = %d (EG_copyEBody)!\n",
             i+1, eedge->nsegs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eedge->nsegs; j++)
      eedge->segs[j].nstart  = NULL;
    for (j = 0; j < eedge->nsegs; j++) {
      eedge->segs[j].iedge  = sedge->segs[j].iedge;
      eedge->segs[j].sense  = sedge->segs[j].sense;
      eedge->segs[j].nstart = sedge->segs[j].nstart;
      eedge->segs[j].tstart = sedge->segs[j].tstart;
      eedge->segs[j].tend   = sedge->segs[j].tend;
    }
    EG_attributeDup(sobj, tobj);
  }
  
  /* populate the ELoops */
  ebody->eloops.objs = (egObject **) EG_alloc(ebody->eloops.nobjs*
                                              sizeof(egObject *));
  if (ebody->eloops.objs == NULL) {
    printf(" EGADS Error: Malloc on %d ELoop Objects (EG_copyEBody)!\n",
           ebody->eloops.nobjs);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eloops.nobjs; i++) ebody->eloops.objs[i] = NULL;
  for (i = 0; i < ebody->eloops.nobjs; i++) {
    sobj  = src->eloops.objs[i];
    sloop = (egELoop *) sobj->blind;
    eloop = (egELoop *) EG_alloc(sizeof(egELoop));
    if (eloop == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop blind (EG_copyEBody)!\n", i+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make ELoop %d Object (EG_copyEBody)!\n",
             i+1);
      EG_free(eloop);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    ebody->eloops.objs[i] = tobj;
    tobj->oclass          = sobj->oclass;
    tobj->mtype           = sobj->mtype;
/*@-kepttrans@*/
    tobj->topObj          = eobj;
/*@+kepttrans@*/
    tobj->blind           = eloop;
    eloop->eedges.nobjs   = sloop->eedges.nobjs;
    eloop->eedges.objs    = NULL;
    eloop->edgeUVs        = NULL;
    eloop->senses         = NULL;
    eloop->nedge          = sloop->nedge;
    eloop->area           = sloop->area;
    if (eloop->eedges.nobjs == 0) {
      EG_free(tobj->blind);
      tobj->blind = NULL;
      continue;
    }
    eloop->eedges.objs = (egObject **) EG_alloc(eloop->eedges.nobjs*
                                                sizeof(egObject *));
    if (eloop->eedges.objs == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop %d Objects (EG_copyEBody)!\n",
             i+1, eloop->eedges.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eloop->eedges.nobjs; j++) eloop->eedges.objs[j] = NULL;
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      k = EG_indexBodyTopo(body, sloop->eedges.objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: EEdge %d index = %d in ELoop %d (EG_copyEBody)!\n",
               k, j+1, i+1);
        EG_destroyEBody(eobj, 1);
        return EGADS_TOPOERR;
      }
      eloop->eedges.objs[j] = ebody->eedges.objs[k-1];
    }
    eloop->senses = (int *) EG_alloc(eloop->eedges.nobjs*sizeof(int));
    if (eloop->senses == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop %d senses (EG_copyEBody)!\n",
             i+1, eloop->nedge);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eloop->eedges.nobjs; j++)
      eloop->senses[j] = sloop->senses[j];
    eloop->edgeUVs = (egEdgeUV *) EG_alloc(eloop->nedge*sizeof(egEdgeUV));
    if (eloop->edgeUVs == NULL) {
      printf(" EGADS Error: Malloc on %d ELoop %d edgeUVs (EG_copyEBody)!\n",
             i+1, eloop->nedge);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eloop->nedge; j++) {
      eloop->edgeUVs[j].npts  = sloop->edgeUVs[j].npts;
      eloop->edgeUVs[j].sense = sloop->edgeUVs[j].sense;
      eloop->edgeUVs[j].edge  = sloop->edgeUVs[j].edge;
      eloop->edgeUVs[j].iuv   = NULL;
    }
    for (j = 0; j < eloop->nedge; j++) {
      eloop->edgeUVs[j].iuv = (int *) EG_alloc(eloop->edgeUVs[j].npts*
                                               sizeof(int));
      if (eloop->edgeUVs == NULL) {
        printf(" EGADS Error: Malloc on %d ELoop %d iUVs %d (EG_copyEBody)!\n",
               i+1, j+1, eloop->edgeUVs[j].npts);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < eloop->edgeUVs[j].npts; k++)
        eloop->edgeUVs[j].iuv[k] = sloop->edgeUVs[j].iuv[k];
    }
    EG_attributeDup(sobj, tobj);
  }
  
  /* populate the EFaces */
  ebody->efaces.objs = (egObject **) EG_alloc(ebody->efaces.nobjs*
                                              sizeof(egObject *));
  if (ebody->efaces.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EFace Objects (EG_readEBody)!\n",
           ebody->efaces.nobjs);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->efaces.nobjs; i++) ebody->efaces.objs[i] = NULL;
  for (i = 0; i < ebody->efaces.nobjs; i++) {
    sobj  = src->efaces.objs[i];
    sface = (egEFace *) sobj->blind;
    eface = (egEFace *) EG_alloc(sizeof(egEFace));
    if (eface == NULL) {
      printf(" EGADS Error: Malloc on %d EFace blind (EG_readEBody)!\n", i+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EFace %d/%d Object (EG_readEBody)!\n",
             i+1, ebody->efaces.nobjs);
      EG_free(eface);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    ebody->efaces.objs[i] = tobj;
    tobj->oclass          = sobj->oclass;
    tobj->mtype           = sobj->mtype;
/*@-kepttrans@*/
    tobj->topObj          = eobj;
/*@+kepttrans@*/
    tobj->blind           = eface;
    eface->trmap          = NULL;
    eface->uvmap          = NULL;
    eface->npatch         = sface->npatch;
    eface->patches        = NULL;
    eface->eloops.nobjs   = sface->eloops.nobjs;
    eface->eloops.objs    = NULL;
    eface->senses         = NULL;
    eface->patches        = NULL;
    eface->range[0]       = sface->range[0];
    eface->range[1]       = sface->range[1];
    eface->range[2]       = sface->range[2];
    eface->range[3]       = sface->range[3];
    eface->last           = sface->last;
    eface->sedges         = ebody->edges;
    if (eface->npatch != 1) {
      stat = EG_uvmapCopy( sface->uvmap,  sface->trmap,
                          &eface->uvmap, &eface->trmap);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: EFace %d  uvmapCopy = %d (EG_copyEBody)!\n",
               i+1, stat);
        EG_destroyEBody(eobj, 1);
        return stat;
      }
    }
    eface->eloops.objs = (egObject **) EG_alloc(eface->eloops.nobjs*
                                                sizeof(egObject *));
    if (eface->eloops.objs == NULL) {
      printf(" EGADS Error: Malloc on %d EFace %d Objects (EG_copyEBody)!\n",
             i+1, eface->eloops.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eface->eloops.nobjs; j++) eface->eloops.objs[j] = NULL;
    for (j = 0; j < eface->eloops.nobjs; j++) {
      k = EG_indexBodyTopo(body, sface->eloops.objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: ELoop %d index = %d in EFace %d (EG_copyEBody)!\n",
               k, j+1, i+1);
        EG_destroyEBody(eobj, 1);
        return EGADS_TOPOERR;
      }
      eface->eloops.objs[j] = ebody->eloops.objs[k-1];
    }
    eface->senses = (int *) EG_alloc(eface->eloops.nobjs*sizeof(int));
    if (eface->senses == NULL) {
      printf(" EGADS Error: Malloc on %d EFace senses %d (EG_copyEBody)!\n",
             i+1, eface->eloops.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eface->eloops.nobjs; j++)
      eface->senses[j] = sface->senses[j];
    eface->patches = (egEPatch *) EG_alloc(eface->npatch*sizeof(egEPatch));
    if (eface->patches == NULL) {
      printf(" EGADS Error: Malloc on %d Patch Object %d (EG_copyEBody)!\n",
             i+1, eface->npatch);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eface->npatch; j++) {
      eface->patches[j].face     = sface->patches[j].face;
      eface->patches[j].start    = sface->patches[j].start;
      eface->patches[j].nuvs     = sface->patches[j].nuvs;
      eface->patches[j].ndeflect = sface->patches[j].ndeflect;
      eface->patches[j].ntris    = sface->patches[j].ntris;
      eface->patches[j].tol      = sface->patches[j].tol;
      eface->patches[j].uvtris   = NULL;
      eface->patches[j].uvtric   = NULL;
      eface->patches[j].uvs      = NULL;
      eface->patches[j].deflect  = NULL;
    }
    for (j = 0; j < eface->npatch; j++) {
      eface->patches[j].uvtris = (int *) EG_alloc(3*eface->patches[j].ntris*
                                                  sizeof(int));
      if (eface->patches[j].uvtris == NULL) {
        printf(" EGADS Error: Malloc on %d Patch uvtris %d (EG_copyEBody)!\n",
               i+1, eface->patches[j].ntris);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < eface->patches[j].ntris; k++) {
        eface->patches[j].uvtris[3*k  ] = sface->patches[j].uvtris[3*k  ];
        eface->patches[j].uvtris[3*k+1] = sface->patches[j].uvtris[3*k+1];
        eface->patches[j].uvtris[3*k+2] = sface->patches[j].uvtris[3*k+2];
      }
      eface->patches[j].uvs = (double *) EG_alloc(2*eface->patches[j].nuvs*
                                                  sizeof(double));
      if (eface->patches[j].uvs == NULL) {
        printf(" EGADS Error: Malloc on %d Patch uvs %d (EG_copyEBody)!\n",
               i+1, eface->patches[j].nuvs);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < eface->patches[j].nuvs; k++) {
        eface->patches[j].uvs[2*k  ] = sface->patches[j].uvs[2*k  ];
        eface->patches[j].uvs[2*k+1] = sface->patches[j].uvs[2*k+1];
      }
      eface->patches[j].deflect = (double *)
                          EG_alloc(3*eface->patches[j].ndeflect*sizeof(double));
      if (eface->patches[j].deflect == NULL) {
        printf(" EGADS Error: Malloc on %d Patch deflect %d (EG_copyEBody)!\n",
               i+1, eface->patches[j].ndeflect);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < eface->patches[j].ndeflect; k++) {
        eface->patches[j].deflect[3*k  ] = sface->patches[j].deflect[3*k  ];
        eface->patches[j].deflect[3*k+1] = sface->patches[j].deflect[3*k+1];
        eface->patches[j].deflect[3*k+2] = sface->patches[j].deflect[3*k+2];
      }
    }
    EG_attributeDup(sobj, tobj);
    if (eface->npatch == 1) {
      eface->patches[0].uvtric = (int *) EG_alloc(3*eface->patches[0].ntris*
                                                  sizeof(int));
      if (eface->patches[0].uvtric == NULL) {
        printf(" EGADS Error: Malloc on %d Patch uvtric %d (EG_copyEBody)!\n",
               i+1, eface->patches[0].ntris);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < eface->patches[0].ntris; k++) {
        eface->patches[0].uvtric[3*k  ] = sface->patches[0].uvtric[3*k  ];
        eface->patches[0].uvtric[3*k+1] = sface->patches[0].uvtric[3*k+1];
        eface->patches[0].uvtric[3*k+2] = sface->patches[0].uvtric[3*k+2];
      }
    }
  }

  /* populate the EShells */
  ebody->eshells.objs = (egObject **) EG_alloc(ebody->eshells.nobjs*
                                               sizeof(egObject *));
  if (ebody->eshells.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EShell Objects (EG_copyEBody)!\n",
           ebody->eshells.nobjs);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eshells.nobjs; i++) ebody->eshells.objs[i] = NULL;
  for (i = 0; i < ebody->eshells.nobjs; i++) {
    sobj   = src->eshells.objs[i];
    sshell = (egEShell *) sobj->blind;
    eshell = (egEShell *) EG_alloc(sizeof(egEShell));
    if (eshell == NULL) {
      printf(" EGADS Error: Malloc on %d EShell blind (EG_copyEBody)!\n", i+1);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EShell %d/%d Object (EG_copyEBody)!\n",
             i+1, ebody->eshells.nobjs);
      EG_free(eshell);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    tobj->oclass           = sobj->oclass;
    tobj->mtype            = sobj->mtype;
/*@-kepttrans@*/
    tobj->topObj           = eobj;
/*@+kepttrans@*/
    ebody->eshells.objs[i] = tobj;
    tobj->blind            = eshell;
    eshell->efaces.nobjs   = sshell->efaces.nobjs;
    eshell->efaces.objs    = (egObject **) EG_alloc(eshell->efaces.nobjs*
                                                    sizeof(egObject *));
    if (eshell->efaces.objs == NULL) {
      printf(" EGADS Error: Malloc on %d EShell %d Objects (EG_copyEBody)!\n",
             i+1, eshell->efaces.nobjs);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < eshell->efaces.nobjs; j++) eshell->efaces.objs[j] = NULL;
    for (j = 0; j < eshell->efaces.nobjs; j++) {
      k = EG_indexBodyTopo(body, sshell->efaces.objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: EFace %d index = %d in EShell %d (EG_copyEBody)!\n",
               k, j+1, i+1);
        EG_destroyEBody(eobj, 1);
        return EGADS_TOPOERR;
      }
      eshell->efaces.objs[j] = ebody->efaces.objs[k-1];
    }
    EG_attributeDup(sobj, tobj);
  }

  *EBody = eobj;
  return EGADS_SUCCESS;
}


/*  ************************* Exposed Entry Points ************************* */

int
EG_initEBody(egObject *tess, double angle, egObject **EBody)
{
  int      i, j, k, m, nchild, nedge, nloop, nface, nshell, index, i0, i1, i2;
  int      stat, oclass, mtype, nbound, *senses;
  double   range[4], result[18];
#ifdef INVEVADJUST
  double   uv[2], xyz[3], dtess, dnode;
#endif
  egTessel *btess;
  egObject *context, *obj, *geom, *eobj, *tobj, **childs, **edges, **loops;
  egObject **faces, **shells;
  egEBody  *ebody;
  egEShell *eshell;
  egEFace  *eface;
  egELoop  *eloop;
  egEEdge  *eedge;
  egEdgeUV **edgsen;

  *EBody = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (tess->blind == NULL)          return EGADS_NODATA;
  if (EG_sameThread(tess))          return EGADS_CNTXTHRD;
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_initEBody)!\n");
    return EGADS_NOTFOUND;
  }
  if (btess->done != 1) {
    printf(" EGADS Error: Tessellation Open (EG_initEBody)!\n");
    return EGADS_TESSTATE;
  }
  if (btess->tess1d == NULL) {
    printf(" EGADS Error: No Edge Tessellations (EG_initEBody)!\n");
    return EGADS_NODATA;
  }
  if (btess->globals == NULL) {
    stat = EG_computeTessMap(btess, EG_outLevel(tess));
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_computeTessMap = %d (EG_initEBody)!\n", stat);
      return stat;
    }
  }
  obj = btess->src;
  if (obj == NULL) {
    printf(" EGADS Error: NULL Source Object (EG_initEBody)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    printf(" EGADS Error: Source Not an Object (EG_initEBody)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    printf(" EGADS Error: Source Not Body (EG_initEBody)!\n");
    return EGADS_NOTBODY;
  }
  if ((obj->mtype != SHEETBODY) && (obj->mtype != SOLIDBODY)) {
    printf(" EGADS Error: Source Body not Solid or Sheet (EG_initEBody)!\n");
    return EGADS_TOPOERR;
  }
  if ((angle < 0.0) || (angle > 90.0)) {
    printf(" EGADS Error: Angle out of range = %lf [0.-90.] (EG_initEBody)!\n",
           angle);
    return EGADS_TOPOERR;
  }
  context = EG_context(obj);
  
  /* fill the EBody with a copy of the translated Body data */
  
  ebody = (egEBody *) EG_alloc(sizeof(egEBody));
  if (ebody == NULL) {
    printf(" EGADS Error: Malloc of EBody  (EG_initEBody)!\n");
    return EGADS_MALLOC;
  }
  ebody->ref           = tess;
  ebody->eedges.objs   = NULL;
  ebody->eloops.objs   = NULL;
  ebody->efaces.objs   = NULL;
  ebody->eshells.objs  = NULL;
  ebody->senses        = NULL;
  ebody->eedges.nobjs  = 0;
  ebody->eloops.nobjs  = 0;
  ebody->efaces.nobjs  = 0;
  ebody->eshells.nobjs = 0;
  ebody->angle         = angle;
  ebody->done          = 0;
  ebody->nedge         = 0;
  ebody->edges         = NULL;
  stat = EG_makeObject(context, &eobj);
  if (stat != EGADS_SUCCESS) {
    EG_free(ebody);
    printf(" EGADS Error: Cannot make EBody Object (EG_initEBody)!\n");
    return stat;
  }
  eobj->oclass = EBODY;
  eobj->mtype  = obj->mtype;
  eobj->blind  = ebody;
  EG_referenceObject(eobj, context);
  EG_referenceTopObj(tess, eobj);
  EG_attributeDup(obj, eobj);
  
  if (obj->mtype == SOLIDBODY) {
    stat = EG_getTopology(obj, &geom, &oclass, &mtype, range, &nchild,
                          &childs, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: getTopology = %d EBody (EG_initEBody)!\n", stat);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    ebody->senses = (int *) EG_alloc(nchild*sizeof(int));
    if (ebody->senses == NULL) {
      printf(" EGADS Error: Malloc on %d Shell senses (EG_initEBody)!\n",
             nchild);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (i = 0; i < nchild; i++) ebody->senses[i] = senses[i];
  }
  
  /* populate the source Edges */
  stat = EG_getBodyTopos(obj, NULL, EDGE, &nedge, &edges);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getBodyTopos Edges = %d (EG_initEBody)!\n", stat);
    EG_destroyEBody(eobj, 1);
    return stat;
  }
  ebody->nedge = nedge;
  ebody->edges = (egEdVert *) EG_alloc(nedge*sizeof(egEdVert));
  if (ebody->edges == NULL) {
    printf(" EGADS Error: Malloc on %d Source Edges (EG_initEBody)!\n", nedge);
    EG_free(edges);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nedge; i++) {
    stat = EG_getTopology(edges[i], &geom, &oclass, &mtype, range, &nchild,
                          &childs, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: getTopology = %d src Edge %d/%d (EG_initEBody)!\n",
             stat, i+1, nedge);
      EG_free(edges);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    ebody->edges[i].edge    = edges[i];
    ebody->edges[i].curve   = 1;
    ebody->edges[i].npts    = 0;
    ebody->edges[i].ts      = NULL;
    ebody->edges[i].dend[0] = ebody->edges[i].dstart[0] = 0.0;
    ebody->edges[i].dend[1] = ebody->edges[i].dstart[1] = 0.0;
    ebody->edges[i].dend[2] = ebody->edges[i].dstart[2] = 0.0;
    ebody->edges[i].ts      = (double *) EG_alloc(btess->tess1d[i].npts*
                                                  sizeof(double));
    if (ebody->edges[i].ts == NULL) {
      printf(" EGADS Error: Malloc on %d ts %d/%d source Edge (EG_initEBody)!\n",
             btess->tess1d[i].npts, i+1, nedge);
      EG_free(edges);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    ebody->edges[i].npts = btess->tess1d[i].npts;
    for (j = 0; j < ebody->edges[i].npts; j++)
      ebody->edges[i].ts[j] = btess->tess1d[i].t[j];
    if (geom->mtype == LINE) ebody->edges[i].curve = 0;
    if (mtype == DEGENERATE) {
      ebody->edges[i].curve = -1;
      continue;
    }
    j    = 0;
    stat = EG_evaluatX(edges[i], &ebody->edges[i].ts[j], result);
    if (stat == EGADS_SUCCESS) {
      ebody->edges[i].dstart[0] = btess->tess1d[i].xyz[3*j  ] - result[0];
      ebody->edges[i].dstart[1] = btess->tess1d[i].xyz[3*j+1] - result[1];
      ebody->edges[i].dstart[2] = btess->tess1d[i].xyz[3*j+2] - result[2];
    }
    j    = ebody->edges[i].npts - 1;
    stat = EG_evaluatX(edges[i], &ebody->edges[i].ts[j], result);
    if (stat == EGADS_SUCCESS) {
      ebody->edges[i].dend[0] = btess->tess1d[i].xyz[3*j  ] - result[0];
      ebody->edges[i].dend[1] = btess->tess1d[i].xyz[3*j+1] - result[1];
      ebody->edges[i].dend[2] = btess->tess1d[i].xyz[3*j+2] - result[2];
    }
  }
  
  /* populate the EEdges */
  ebody->eedges.nobjs = nedge;
  ebody->eedges.objs  = (egObject **) EG_alloc(nedge*sizeof(egObject *));
  if (ebody->eedges.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EEdge Objects (EG_initEBody)!\n", nedge);
    EG_free(edges);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nedge; i++) ebody->eedges.objs[i] = NULL;
  for (i = 0; i < nedge; i++) {
    stat = EG_getTopology(edges[i], &geom, &oclass, &mtype, range, &nchild,
                          &childs, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: getTopology = %d EEdge %d/%d (EG_initEBody)!\n",
             stat, i+1, nedge);
      EG_free(edges);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    if (nchild == 0) continue;
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EEdge %d/%d Object (EG_initEBody)!\n",
             i+1, nedge);
      EG_free(edges);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    tobj->oclass          = EEDGE;
    tobj->mtype           = mtype;
    tobj->topObj          = eobj;
    ebody->eedges.objs[i] = tobj;
    eedge = (egEEdge *) EG_alloc(sizeof(egEEdge));
    if (eedge == NULL) {
      printf(" EGADS Error: Malloc on %d/%d EEdge blind (EG_initEBody)!\n",
             i+1, nedge);
      EG_free(edges);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    tobj->blind      = eedge;
    eedge->sedges    = ebody->edges;
    eedge->trange[0] = range[0];
    eedge->trange[1] = range[1];
    eedge->nodes[0]  = eedge->nodes[1] = childs[0];
    if (mtype == TWONODE) eedge->nodes[1] = childs[1];
    eedge->nsegs     = 0;
    eedge->segs      = (egEEseg *) EG_alloc(sizeof(egEEseg));
    if (eedge->segs == NULL) {
      printf(" EGADS Error: Malloc on %d/%d EEdge segs (EG_initEBody)!\n",
             i+1, nedge);
      EG_free(edges);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    EG_attributeDup(edges[i], tobj);
    eedge->nsegs        = 1;
    eedge->segs->iedge  = i+1;
    eedge->segs->sense  = SFORWARD;
    eedge->segs->nstart = NULL;
    eedge->segs->tstart = range[0];
    eedge->segs->tend   = range[1];
  }
  EG_free(edges);
  
  /* populate the ELoops */
  stat = EG_getBodyTopos(obj, NULL, LOOP, &nloop, &loops);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getBodyTopos Loops = %d (EG_initEBody)!\n", stat);
    EG_destroyEBody(eobj, 1);
    return stat;
  }
  ebody->eloops.nobjs = nloop;
  ebody->eloops.objs  = (egObject **) EG_alloc(nloop*sizeof(egObject *));
  if (ebody->eloops.objs == NULL) {
    printf(" EGADS Error: Malloc on %d ELoop Objects (EG_initEBody)!\n", nloop);
    EG_free(loops);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nloop; i++) ebody->eloops.objs[i] = NULL;
  for (i = 0; i < nloop; i++) {
    stat = EG_getTopology(loops[i], &geom, &oclass, &mtype, range, &nchild,
                          &childs, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: getTopology = %d ELoop %d/%d (EG_initEBody)!\n",
             stat, i+1, nloop);
      EG_free(loops);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    if (nchild == 0) continue;
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make ELoop %d/%d Object (EG_initEBody)!\n",
             i+1, nloop);
      EG_free(loops);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    tobj->oclass          = ELOOPX;
    tobj->mtype           = mtype;
/*@-kepttrans@*/
    tobj->topObj          = eobj;
/*@+kepttrans@*/
    tobj->blind           = NULL;
    ebody->eloops.objs[i] = tobj;
    eloop = (egELoop *) EG_alloc(sizeof(egELoop));
    if (eloop == NULL) {
      printf(" EGADS Error: Malloc on %d/%d ELoop blind (EG_initEBody)!\n",
             i+1, nloop);
      EG_free(loops);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    EG_attributeDup(loops[i], tobj);
    tobj->blind         = eloop;
    eloop->eedges.nobjs = 0;
    eloop->eedges.objs  = NULL;
    eloop->nedge        = nchild;
    eloop->edgeUVs      = NULL;
    eloop->area         = 0.0;
    eloop->senses       = (int *) EG_alloc(nchild*sizeof(int));
    if (eloop->senses == NULL) {
      printf(" EGADS Error: Malloc on %d/%d ELoop %d senses (EG_initEBody)!\n",
             i+1, nloop, nchild);
      EG_free(loops);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    eloop->edgeUVs      = (egEdgeUV *) EG_alloc(nchild*sizeof(egEdgeUV));
    if (eloop->edgeUVs == NULL) {
      printf(" EGADS Error: Malloc on %d/%d ELoop %d edgeUVs (EG_initEBody)!\n",
             i+1, nloop, nchild);
      EG_free(loops);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nchild; j++) {
      eloop->senses[j]        = senses[j];
      eloop->edgeUVs[j].npts  = 0;
      eloop->edgeUVs[j].sense = senses[j];
      eloop->edgeUVs[j].edge  = childs[j];
      eloop->edgeUVs[j].iuv   = NULL;
    }
    eloop->eedges.objs = (egObject **) EG_alloc(nchild*sizeof(egObject *));
    if (eloop->eedges.objs == NULL) {
      printf(" EGADS Error: Malloc on %d/%d ELoop %d Objects (EG_initEBody)!\n",
             i+1, nloop, nchild);
      EG_free(loops);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nchild; j++) eloop->eedges.objs[j] = NULL;
    for (j = 0; j < nchild; j++) {
      k = EG_indexBodyTopo(obj, childs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Loop index = %d %d/%d on %d/%d (EG_initEBody)!\n",
               k, j+1, nchild, i+1, nloop);
        EG_free(loops);
        EG_destroyEBody(eobj, 1);
        return k;
      }
      eloop->eedges.objs[j] = ebody->eedges.objs[k-1];
    }
    eloop->eedges.nobjs = nchild;
    for (m = 0; m < nchild; m++) {
      if (childs[m]->mtype == DEGENERATE) continue;
      j = EG_indexBodyTopo(obj, childs[m]) - 1;
      if (j < EGADS_SUCCESS) continue;
      eloop->edgeUVs[m].npts = btess->tess1d[j].npts;
      if (eloop->edgeUVs[m].npts <= 0) continue;
      eloop->edgeUVs[m].iuv = (int *) EG_alloc(btess->tess1d[j].npts*
                                               sizeof(int));
      if (eloop->edgeUVs[m].iuv == NULL) {
        printf(" EGADS Error: MALLOC on Edge %d/%d for %d UVs (EG_initEBody)!\n",
               i+1, nchild, eloop->edgeUVs[j].npts);
        EG_free(loops);
        EG_destroyEBody(eobj, 1);
        return EGADS_MALLOC;
      }
      for (k = 0; k < btess->tess1d[j].npts; k++) eloop->edgeUVs[m].iuv[k] = 0;
    }
  }
  EG_free(loops);
  
  /* populate the EFaces */
  stat = EG_getBodyTopos(obj, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getBodyTopos Faces = %d (EG_initEBody)!\n", stat);
    EG_destroyEBody(eobj, 1);
    return stat;
  }
  for (j = i = 0; i < nface; i++)
    if (btess->tess2d[i].npts == 0) j++;
  if (j != 0) {
    printf(" EGADS Error: %d Faces w/out Triangles (EG_initEBody)!\n", j);
    EG_free(faces);
    EG_destroyEBody(eobj, 1);
    return stat;
  }
  ebody->efaces.nobjs = nface;
  ebody->efaces.objs  = (egObject **) EG_alloc(nface*sizeof(egObject *));
  if (ebody->efaces.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EFace Objects (EG_initEBody)!\n", nface);
    EG_free(faces);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nface; i++) ebody->efaces.objs[i] = NULL;
  for (i = 0; i < nface; i++) {
    stat = EG_getTopology(faces[i], &geom, &oclass, &mtype, range, &nchild,
                          &childs, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: getTopology = %d EFace %d/%d (EG_initEBody)!\n",
             stat, i+1, nface);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make EFace %d/%d Object (EG_initEBody)!\n",
             i+1, nface);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    tobj->oclass          = EFACE;
    tobj->mtype           = mtype;
/*@-kepttrans@*/
    tobj->topObj          = eobj;
/*@+kepttrans@*/
    ebody->efaces.objs[i] = tobj;
    eface = (egEFace *) EG_alloc(sizeof(egEFace));
    if (eface == NULL) {
      printf(" EGADS Error: Malloc on %d/%d EFace blind (EG_initEBody)!\n",
             i+1, nface);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    EG_attributeDup(faces[i], tobj);
    tobj->blind         = eface;
    eface->sedges       = ebody->edges;
    eface->npatch       = 0;
    eface->patches      = NULL;
    eface->eloops.nobjs = 0;
    eface->eloops.objs  = NULL;
    eface->senses       = NULL;
    eface->uvmap        = NULL;
    eface->trmap        = NULL;
    eface->range[0]     = range[0];
    eface->range[1]     = range[1];
    eface->range[2]     = range[2];
    eface->range[3]     = range[3];
    eface->last         = 0;
/*@-kepttrans@*/
    eface->sedges       = ebody->edges;
/*@+kepttrans@*/
    eface->patches      = (egEPatch *) EG_alloc(sizeof(egEPatch));
    if (eface->patches == NULL) {
      printf(" EGADS Error: Malloc on %d/%d Face Object (EG_initEBody)!\n",
             i+1, nface);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    eface->patches[0].start    = 0;
    eface->patches[0].nuvs     = 0;
    eface->patches[0].ndeflect = 0;
    eface->patches[0].ntris    = 0;
    eface->patches[0].tol      = -1.0;
    eface->patches[0].uvtris   = NULL;
    eface->patches[0].uvtric   = NULL;
    eface->patches[0].uvs      = NULL;
    eface->patches[0].deflect  = NULL;
    eface->patches[0].face     = faces[i];
    eface->npatch              = 1;
    
    eface->eloops.objs   = (egObject **) EG_alloc(nchild*sizeof(egObject *));
    if (eface->eloops.objs == NULL) {
      printf(" EGADS Error: Malloc on %d/%d ELoop %d Objects (EG_initEBody)!\n",
             i+1, nface, nchild);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nchild; j++) eface->eloops.objs[j] = NULL;
    for (j = 0; j < nchild; j++) {
      k = EG_indexBodyTopo(obj, childs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Face index = %d %d/%d on %d/%d (EG_initEBody)!\n",
               k, j+1, nchild, i+1, nface);
        EG_free(faces);
        EG_destroyEBody(eobj, 1);
        return k;
      }
      eface->eloops.objs[j] = ebody->eloops.objs[k-1];
    }
    eface->eloops.nobjs = nchild;
    eface->senses       = (int *) EG_alloc(nchild*sizeof(int));
    if (eface->senses == NULL) {
      printf(" EGADS Error: Malloc on %d/%d senses %d (EG_initEBody)!\n",
             i+1, nface, nchild);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nchild; j++) eface->senses[j] = senses[j];
    eface->patches[0].uvs = (double *)
                            EG_alloc(2*btess->tess2d[i].npts*sizeof(double));
    if (eface->patches[0].uvs == NULL) {
      printf(" EGADS Error: Malloc on Face %d/%d UVs %d (EG_initEBody)!\n",
             i+1, nface, btess->tess2d[i].npts);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < 2*btess->tess2d[i].npts; j++)
      eface->patches[0].uvs[j] = btess->tess2d[i].uv[j];
    eface->patches[0].nuvs   = btess->tess2d[i].npts;
    eface->patches[0].uvtris = (int *)
                               EG_alloc(3*btess->tess2d[i].ntris*sizeof(int));
    if (eface->patches[0].uvtris == NULL) {
      printf(" EGADS Error: Malloc on Face %d/%d UVtris %d (EG_initEBody)!\n",
             i+1, nface, btess->tess2d[i].ntris);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    eface->patches[0].uvtric = (int *)
                               EG_alloc(3*btess->tess2d[i].ntris*sizeof(int));
    if (eface->patches[0].uvtric == NULL) {
      printf(" EGADS Error: Malloc on Face %d/%d UVtric %d (EG_initEBody)!\n",
             i+1, nface, btess->tess2d[i].ntris);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < 3*btess->tess2d[i].ntris; j++) {
      eface->patches[0].uvtris[j] = btess->tess2d[i].tris[j];
      eface->patches[0].uvtric[j] = btess->tess2d[i].tric[j];
    }
    eface->patches[0].ntris = btess->tess2d[i].ntris;
    
    nbound = btess->tess2d[i].frlps[btess->tess2d[i].nfrlps-1];
    eface->patches[0].deflect = (double *) EG_alloc(3*nbound*sizeof(double));
    if (eface->patches[0].deflect == NULL) {
      printf(" EGADS Error: Malloc on Face %d/%d deflect %d (EG_initEBody)!\n",
             i+1, nface, nbound);
      EG_free(faces);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nbound; j++) {
      eface->patches[0].deflect[3*j  ] = 0.0;
      eface->patches[0].deflect[3*j+1] = 0.0;
      eface->patches[0].deflect[3*j+2] = 0.0;
      stat = EG_evaluatX(faces[i], &eface->patches[0].uvs[2*j], result);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Warning: Face %d/%d  evaluate = %d (EG_initEBody)!\n",
               i+1, nface, stat);
      } else {
        eface->patches[0].deflect[3*j  ] = btess->tess2d[i].xyz[3*j  ]-result[0];
        eface->patches[0].deflect[3*j+1] = btess->tess2d[i].xyz[3*j+1]-result[1];
        eface->patches[0].deflect[3*j+2] = btess->tess2d[i].xyz[3*j+2]-result[2];
      }
    }
    eface->patches[0].ndeflect = nbound;
  }
  EG_free(faces);
  
  /* fill in iuv */
  
  edgsen = (egEdgeUV **) EG_alloc(3*btess->nEdge*sizeof(egEdgeUV *));
  if (edgsen == NULL) {
    printf(" EGADS Error: Allocating %d egdeUVs (EG_initEBody)!\n",
           btess->nEdge);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (j = 0; j < btess->nFace; j++) {
    eface = (egEFace *) ebody->efaces.objs[j]->blind;
    for (i = 0; i < 3*btess->nEdge; i++) edgsen[i] = NULL;
    for (i = 0; i < eface->eloops.nobjs; i++) {
      if (eface->eloops.objs[i] == NULL) continue;
      eloop = (egELoop *) eface->eloops.objs[i]->blind;
      for (k = 0; k < eloop->nedge; k++) {
        index = EG_indexBodyTopo(obj, eloop->edgeUVs[k].edge);
        if (index <= EGADS_SUCCESS) continue;
        if (eloop->edgeUVs[k].edge->mtype == DEGENERATE) continue;
        if (eloop->edgeUVs[k].sense == SREVERSE) {
          edgsen[3*index-3] = &eloop->edgeUVs[k];
        } else {
          edgsen[3*index-2] = &eloop->edgeUVs[k];
        }
      }
    }
    for (k = 0; k < btess->tess2d[j].ntris; k++) {
      i0 = btess->tess2d[j].tris[3*k  ] - 1;
      i1 = btess->tess2d[j].tris[3*k+1] - 1;
      i2 = btess->tess2d[j].tris[3*k+2] - 1;
      if (btess->tess2d[j].tric[3*k  ] < 0) {
        index = -btess->tess2d[j].tric[3*k  ] - 1;
        EG_getEdgeIndices(btess, j, index, i1, i2, 0, edgsen);
      }
      if (btess->tess2d[j].tric[3*k+1] < 0) {
        index = -btess->tess2d[j].tric[3*k+1] - 1;
        EG_getEdgeIndices(btess, j, index, i2, i0, 0, edgsen);
      }
      if (btess->tess2d[j].tric[3*k+2] < 0) {
        index = -btess->tess2d[j].tric[3*k+2] - 1;
        EG_getEdgeIndices(btess, j, index, i0, i1, 0, edgsen);
      }
    }
  }
  EG_free(edgsen);
  
  /* populate the EShells */
  stat = EG_getBodyTopos(obj, NULL, SHELL, &nshell, &shells);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getBodyTopos Shell = %d (EG_initEBody)!\n", stat);
    EG_destroyEBody(eobj, 1);
    return stat;
  }
  ebody->eshells.nobjs = nshell;
  ebody->eshells.objs  = (egObject **) EG_alloc(nshell*sizeof(egObject *));
  if (ebody->eshells.objs == NULL) {
    printf(" EGADS Error: Malloc on %d EShell Objects (EG_initEBody)!\n",
           nshell);
    EG_free(shells);
    EG_destroyEBody(eobj, 1);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nshell; i++) ebody->eshells.objs[i] = NULL;
  for (i = 0; i < nshell; i++) {
    stat = EG_getTopology(shells[i], &geom, &oclass, &mtype, range, &nchild,
                          &childs, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: getTopology = %d EShell %d/%d (EG_initEBody)!\n",
             stat, i+1, nshell);
      EG_free(shells);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    stat = EG_makeObject(context, &tobj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Cannot make ELoop %d/%d Object (EG_initEBody)!\n",
             i+1, nloop);
      EG_free(shells);
      EG_destroyEBody(eobj, 1);
      return stat;
    }
    tobj->oclass           = ESHELL;
    tobj->mtype            = mtype;
/*@-kepttrans@*/
    tobj->topObj           = eobj;
/*@+kepttrans@*/
    ebody->eshells.objs[i] = tobj;
    eshell = (egEShell *) EG_alloc(sizeof(egEShell));
    if (eshell == NULL) {
      printf(" EGADS Error: Malloc on %d/%d EShell blind (EG_initEBody)!\n",
             i+1, nshell);
      EG_free(shells);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    EG_attributeDup(shells[i], tobj);
    tobj->blind          = eshell;
    eshell->efaces.nobjs = 0;
    eshell->efaces.objs  = (egObject **) EG_alloc(nchild*sizeof(egObject *));
    if (eshell->efaces.objs == NULL) {
      printf(" EGADS Error: Malloc on %d/%d EShell %d Objects (EG_initEBody)!\n",
             i+1, nshell, nchild);
      EG_free(shells);
      EG_destroyEBody(eobj, 1);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nchild; j++) eshell->efaces.objs[j] = NULL;
    for (j = 0; j < nchild; j++) {
      k = EG_indexBodyTopo(obj, childs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: Shell index = %d %d/%d on %d/%d (EG_initEBody)!\n",
               k, j+1, nchild, i+1, nshell);
        EG_free(shells);
        EG_destroyEBody(eobj, 1);
        return k;
      }
      eshell->efaces.objs[j] = ebody->efaces.objs[k-1];
    }
    eshell->efaces.nobjs = nchild;
  }
  EG_free(shells);
  
  stat = EG_eRemoveNode(eobj);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_eRemoveNode = %d (EG_initEBody)!\n", stat);
    EG_destroyEBody(eobj, 1);
    return stat;
  }
  
  *EBody = eobj;
  return EGADS_SUCCESS;
}


int
EG_finishEBody(egObject *EBody)
{
  int      i, j;
  egObject *eobj;
  egEBody  *ebody;
  egEFace  *eface;
  egTessel *btess;
  
  if (EBody == NULL)               return EGADS_NULLOBJ;
  if (EBody->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EBody->oclass != EBODY)      return EGADS_NOTTESS;
  if (EG_sameThread(EBody))        return EGADS_CNTXTHRD;
  
  ebody = (egEBody *) EBody->blind;
  if (ebody == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_finishEBody)!\n");
    return EGADS_NOTFOUND;
  }
  if (ebody->done == 1) {
    printf(" EGADS Error: EBody already closed (EG_finishEBody)!\n");
    return EGADS_EXISTS;
  }
  
  /* cleanup uvtric storage for multi-patch EFaces */
  for (i = 0; i < ebody->efaces.nobjs; i++) {
    eobj = ebody->efaces.objs[i];
    if (eobj == NULL) continue;
    if (eobj->blind == NULL) continue;
    eface = (egEFace *) eobj->blind;
    if (eface->npatch == 1) continue;
    for (j = 0; j < eface->npatch; j++) {
      EG_free(eface->patches[j].uvtric);
      eface->patches[j].uvtric = NULL;
    }
  }
  
  /* dereference tess & reference body */
  btess = (egTessel *) ebody->ref->blind;
  EG_dereferenceTopObj(ebody->ref, EBody);
  EG_referenceTopObj(btess->src, EBody);
  ebody->ref  = btess->src;
  ebody->done = 1;

  return EGADS_SUCCESS;
}


int
EG_makeEFace(egObject *EBody, int nFace, egObject **Faces, egObject **EFace)
{
  int          i, j, k, m, n, ie, stat, index, neloop, nvrt, ntri, j0n, newel;
  int          atype, alen, i0, i1, i2, *iedge, *map, *tris, *itri, *trmap;
  int          outLevel, fullAttr, iface = -1;
  double       *xyzs, range[4], angle, maxang, big, biga;
  void         *uvmap;
  egObject     *body, *obj, **objs, *eshell, *saved, **efaces, **eloops;
  egEBody      *ebody;
  egEShell     *esh;
  egEFace      *eface, *ef;
  egEdgeUV     **edgsen;
  egELoop      *eloop, **newels;
  egEEdge      *eedge;
  egTessel     *btess;
  eloopFrag    *frags;
  const int    *ints;
  const double *reals;
  const char   *str;
  
  *EFace = NULL;
  if (EBody == NULL)               return EGADS_NULLOBJ;
  if (EBody->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EBody->oclass != EBODY)      return EGADS_NOTBODY;
  if (EG_sameThread(EBody))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(EBody);
  
  ebody = (egEBody *) EBody->blind;
  if (ebody == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_makeEFace)!\n");
    return EGADS_NOTFOUND;
  }
  if (ebody->done == 1) {
    printf(" EGADS Error: EBody finialized (EG_makeEFace)!\n");
    return EGADS_EXISTS;
  }
  btess    = (egTessel *) ebody->ref->blind;
  body     = btess->src;
  fullAttr = EG_fullAttrs(body);
  
  for (i = 0; i < nFace-1; i++)
    for (j = i+1; j < nFace; j++)
      if (Faces[i] == Faces[j]) {
        printf(" EGADS Error: Face %d is same as %d (EG_makeEFace)!\n",
               i+1, j+1);
        return EGADS_NOTTOPO;
      }
  
  /* are the Faces OK and not in another EFace? */
  for (neloop = i = 0; i < nFace; i++) {
    if (Faces[i] == NULL) {
      printf(" EGADS Error: Face %d is NULL (EG_makeEFace)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (Faces[i]->magicnumber != MAGIC) {
      printf(" EGADS Error: Face %d is no Object (EG_makeEFace)!\n", i+1);
      return EGADS_NOTOBJ;
    }
    if (Faces[i]->oclass != FACE) {
      printf(" EGADS Error: Face %d is NOT a Face (EG_makeEFace)!\n", i+1);
      return EGADS_NOTTOPO;
    }
    j = EG_indexBodyTopo(body, Faces[i]);
    if (j <= EGADS_SUCCESS) {
      printf(" EGADS Error: Face %d is NOT in Body = %d (EG_makeEFace)!\n",
             i+1, j);
      return EGADS_TOPOERR;
    }
#ifdef DEBUG
    printf(" Face %d  or = %2d\n", j, Faces[i]->mtype);
#endif
    for (j = 0; j < ebody->efaces.nobjs; j++) {
      obj   = ebody->efaces.objs[j];
      if (obj == NULL) continue;
      eface = (egEFace *) obj->blind;
      if (eface == NULL) continue;
      neloop += eface->eloops.nobjs;
      if (eface->npatch == 1) continue;
      for (k = 0; k < eface->npatch; k++)
        if (Faces[i] == eface->patches[k].face) {
          printf(" EGADS Error: Face %d is already taken (EG_makeEFace)!\n",
                 i+1);
          return EGADS_TOPOERR;
        }
    }
  }
  
  /* are we connected? */
  iedge = (int *) EG_alloc(btess->nEdge*sizeof(int));
  if (iedge == NULL) {
    printf(" EGADS Error: Malloc of %d Edge markers (EG_makeEFace)!\n",
           btess->nEdge);
    return EGADS_MALLOC;
  }
  for (i = 0; i < btess->nEdge; i++) iedge[i] = 0;
  for (i = 0; i < nFace; i++) {
    stat = EG_getBodyTopos(body, Faces[i], EDGE, &n, &objs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_getBodyTopos %d = %d (EG_makeEFace)!\n",
             i+1, stat);
      EG_free(iedge);
      return stat;
    }
    for (j = 0; j < n; j++) {
      if (objs[j]->mtype == DEGENERATE) continue;
      stat = EG_attributeRet(objs[j], ".Keep", &atype, &alen, &ints, &reals,
                             &str);
      if (stat == EGADS_SUCCESS) continue;
      k = EG_indexBodyTopo(body, objs[j]);
      if (k <= EGADS_SUCCESS) {
        printf(" EGADS Error: %d/%d EG_indexBodyTopo = %d (EG_makeEFace)!\n",
               i+1, j+1, k);
        EG_free(objs);
        EG_free(iedge);
        return k;
      }
      iedge[k-1]++;
    }
    EG_free(objs);
  }
  for (i = 0; i < nFace; i++) {
    EG_getBodyTopos(body, Faces[i], EDGE, &n, &objs);
    for (j = 0; j < n; j++) {
      if (objs[j]->mtype == DEGENERATE) continue;
      k = EG_indexBodyTopo(body, objs[j]);
      if (iedge[k-1] > 1) break;
    }
    EG_free(objs);
    if (j == n) {
      printf(" EGADS Error: %d Not Connected (EG_makeEFace)!\n", i+1);
      EG_free(iedge);
      return EGADS_TOPOCNT;
    }
  }
  EG_free(iedge);
  
  /* collect the objects in play */
  efaces = (egObject **) EG_alloc((nFace+neloop)*sizeof(egObject *));
  if (efaces == NULL) {
    printf(" EGADS Error: Malloc of %d %d Objects (EG_makeEFace)!\n",
           nFace, neloop);
    return EGADS_MALLOC;
  }
  eloops = &efaces[nFace];
  for (neloop = i = 0; i < nFace; i++) {
    for (j  = 0; j < ebody->efaces.nobjs; j++) {
      obj   = ebody->efaces.objs[j];
      if (obj == NULL) continue;
      eface = (egEFace *) obj->blind;
      if (eface == NULL) continue;
      if (eface->npatch != 1) continue;
      if (Faces[i] != eface->patches[0].face) continue;
      efaces[i] = obj;
      for (k = 0; k < eface->eloops.nobjs; k++, neloop++)
        eloops[neloop] = eface->eloops.objs[k];
      break;
    }
  }
  /* find our EShell */
  eshell = NULL;
  for (i = 0; i < ebody->eshells.nobjs; i++) {
    stat = EG_getBodyTopos(EBody, ebody->eshells.objs[i], EFACE, &m, &objs);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Warning: %d EG_getBodyTopos = %d (EG_makeEFace)!\n",
             i+1, stat);
      continue;
    }
    for (k = j = 0; j < m; j++)
      for (n = 0; n < nFace; n++)
        if (efaces[n] == objs[j]) k++;
    EG_free(objs);
    if (k != 0)
      if (eshell == NULL) {
        if (k != nFace) {
          printf(" EGADS Error: EShell %d has %d EFaces %d (EG_makeEFace)!\n",
                 i+1, k, nFace);
          EG_free(efaces);
          return EGADS_TOPOERR;
        }
        eshell = ebody->eshells.objs[i];
      } else {
        printf(" EGADS Error: EShell %d has %d EFaces %d (EG_makeEFace)!\n",
               i+1, k, nFace);
        EG_free(efaces);
        return EGADS_TOPOERR;
      }
  }
  if (eshell == NULL) {
    printf(" EGADS Error: EShell not found for EFaces (EG_makeEFace)!\n");
    EG_free(efaces);
    return EGADS_TOPOERR;
  }
  
  /* find the EEdges that disappear */
  iedge = (int *) EG_alloc(ebody->eedges.nobjs*sizeof(int));
  if (iedge == NULL) {
    printf(" EGADS Error: Malloc on %d EEdges (EG_makeEFace)!\n",
           ebody->eedges.nobjs);
    EG_free(efaces);
    return EGADS_MALLOC;
  }
  for (i = 0; i < ebody->eedges.nobjs; i++) iedge[i] = 0;
  for (i = 0; i < neloop; i++) {
    if (eloops[i] == NULL) continue;
    eloop = (egELoop *) eloops[i]->blind;
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      index = EG_indexBodyTopo(EBody, eloop->eedges.objs[j]);
      if (index > EGADS_SUCCESS) iedge[index-1]++;
    }
  }
  /* mark non-manifold Edges w/ more than 2 Faces */
  for (i = 0; i < ebody->eedges.nobjs; i++)
    if (iedge[i] > 2) iedge[i] = -1;
  /* mark kept EEdges */
  for (i = 0; i < neloop; i++) {
    if (eloops[i] == NULL) continue;
    eloop = (egELoop *) eloops[i]->blind;
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      if (eloop->eedges.objs[j] == NULL) continue;
      index = EG_indexBodyTopo(EBody, eloop->eedges.objs[j]);
      if (index <= EGADS_SUCCESS) continue;
      if (iedge[index-1] <= 1)    continue;
      /* Only 2 Node EEdges are candidates
      if (eloop->eedges.objs[j]->mtype != TWONODE) {
        iedge[index-1] = -1;
        continue;
      }  */
      eedge = (egEEdge *) eloop->eedges.objs[j]->blind;
      /* Edges with .Keep attribute cannot be removed! */
      if (eedge->nsegs == 1) {
        ie   = eedge->segs[0].iedge - 1;
        stat = EG_attributeRet(eedge->sedges[ie].edge, ".Keep", &atype, &alen,
                               &ints, &reals, &str);
        if (stat == EGADS_SUCCESS) {
          iedge[index-1] = -1;
          continue;
        }
      }
      /* Check winding angle */
      maxang = 0.0;
      for (k = 0; k < eedge->nsegs; k++) {
        ie = eedge->segs[k].iedge - 1;
        for (m = 1; m < eedge->sedges[ie].npts-1; m++) {
          stat = EG_getWindingAngle(eedge->sedges[ie].edge,
                                    eedge->sedges[ie].ts[m], &angle);
          if (stat != EGADS_SUCCESS) {
            printf(" EGADS Warning: getWindingAngle %d = %d (EG_makeEFace)!\n",
                   m+1, stat);
            continue;
          }
          if (maxang < fabs(angle-180.0)) maxang = fabs(angle-180.0);
        }
        if (eedge->sedges[ie].npts == 2) {
          stat = EG_getWindingAngle(eedge->sedges[ie].edge,
                                    0.5*(eedge->sedges[ie].ts[0] +
                                         eedge->sedges[ie].ts[1]), &angle);
          if (stat != EGADS_SUCCESS) {
            printf(" EGADS Warning: getWindingAngle = %d (EG_makeEFace)!\n",
                   stat);
          } else {
            if (maxang < fabs(angle-180.0)) maxang = fabs(angle-180.0);
          }
        }
      }
      if (maxang > ebody->angle) {
        if (outLevel > 0)
          printf(" EGADS Info: Dihedral deviation = %lf (%lf) for EEdge %d\n",
                 maxang, ebody->angle, index);
        iedge[index-1] = -1;
        continue;
      }
      /* ignore if in ELoop twice */
      for (k = j+1; k < eloop->eedges.nobjs; k++)
        if (eloop->eedges.objs[j] == eloop->eedges.objs[k]) break;
      if (k != eloop->eedges.nobjs) iedge[index-1] = -1;
    }
  }
  
  /* n is the number of retained EEdges */
  for (k = n = i = 0; i < ebody->eedges.nobjs; i++) {
    if (ebody->eedges.objs[i] == NULL) continue;
    if (ebody->eedges.objs[i]->mtype == DEGENERATE) {
      iedge[i] = -2;
      continue;
    }
    if (abs(iedge[i]) == 1) n++;
    if (    iedge[i]  >  1) {
      k++;
#ifdef DEBUG
      printf(" Remove EEdge %d -- %d\n", i+1, iedge[i]);
#endif
    }
  }
  if (k == 0) {
    printf(" EGADS Error: No EEdges removed (EG_makeEFace)!\n");
    EG_free(iedge);
    EG_free(efaces);
    return EGADS_TOPOERR;
  }
  /* conservative -- edges can be in the frag list twice! */
  n *= 2;
  
  /* save away our current state */
  stat = EG_copyEBody(EBody, &saved);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: copyEBody = %d (EG_makeEFace)!\n", stat);
    EG_free(iedge);
    EG_free(efaces);
    return stat;
  }
  
  /* collect the Edge/sense edgeUVs in play */
  edgsen = (egEdgeUV **) EG_alloc(3*btess->nEdge*sizeof(egEdgeUV *));
  if (edgsen == NULL) {
    printf(" EGADS Error: Allocating %d egdeUVs (EG_makeEFace)!\n",
           btess->nEdge);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  for (i = 0; i < 3*btess->nEdge; i++) edgsen[i] = NULL;
  for (i = 0; i < neloop; i++) {
    if (eloops[i] == NULL) continue;
    eloop = (egELoop *) eloops[i]->blind;
    for (j = 0; j < eloop->nedge; j++) {
      index = EG_indexBodyTopo(body, eloop->edgeUVs[j].edge);
      if (index <= EGADS_SUCCESS) continue;
      if (eloop->edgeUVs[j].edge->mtype == DEGENERATE) continue;
      for (m = 0; m < eloop->edgeUVs[j].npts; m++)
        eloop->edgeUVs[j].iuv[m] = 0;
      if (eloop->edgeUVs[j].sense == SREVERSE) {
        edgsen[3*index-3] = &eloop->edgeUVs[j];
      } else {
        edgsen[3*index-2] = &eloop->edgeUVs[j];
      }
    }
  }
  /* mark the removed Edges */
  for (i = 0; i < ebody->eedges.nobjs; i++)
    if (iedge[i] > 1) {
      eedge = (egEEdge *) ebody->eedges.objs[i]->blind;
      for (j = 0; j < eedge->nsegs; j++) {
        index             = eedge->segs[j].iedge;
        edgsen[3*index-1] = edgsen[3*index-2];
      }
    }
/*
  for (i = 0; i < btess->nEdge; i++)
    printf(" %2d: %lx %lx %lx\n", i+1, (long) edgsen[3*i  ],
           (long) edgsen[3*i+1], (long) edgsen[3*i+2]);  */
  for (m = i = 0; i < nFace; i++) {
    j  = EG_indexBodyTopo(body, Faces[i]);
    m += btess->tess2d[j-1].npts;
  }
  map = (int *) EG_alloc(m*sizeof(int));
  if (map == NULL) {
    printf(" EGADS Error: Malloc of %d verts (EG_makeEFace)!\n", m);
    EG_free(edgsen);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  for (i = 0; i < m; i++) map[i] = 0;
  for (m = i = 0; i < nFace; i++) {
    j = EG_indexBodyTopo(body, Faces[i]);
    for (k = 0; k < btess->tess2d[j-1].ntris; k++) {
      i0 = btess->tess2d[j-1].tris[3*k  ] - 1;
      i1 = btess->tess2d[j-1].tris[3*k+1] - 1;
      i2 = btess->tess2d[j-1].tris[3*k+2] - 1;
      if (btess->tess2d[j-1].tric[3*k  ] < 0) {
        index = -btess->tess2d[j-1].tric[3*k  ] - 1;
        EG_getEdgeIndices(btess, j-1, index, i1, i2, m, edgsen);
      }
      if (btess->tess2d[j-1].tric[3*k+1] < 0) {
        index = -btess->tess2d[j-1].tric[3*k+1] - 1;
        EG_getEdgeIndices(btess, j-1, index, i2, i0, m, edgsen);
      }
      if (btess->tess2d[j-1].tric[3*k+2] < 0) {
        index = -btess->tess2d[j-1].tric[3*k+2] - 1;
        EG_getEdgeIndices(btess, j-1, index, i0, i1, m, edgsen);
      }
    }
    m += btess->tess2d[j-1].npts;
  }
  /* mark the verts to be removed */
  for (i = 0; i < btess->nEdge; i++) {
    if (edgsen[3*i+2] == NULL) continue;
    for (j = 0; j < edgsen[3*i]->npts; j++) {
      i1 = edgsen[3*i  ]->iuv[j];
      i2 = edgsen[3*i+1]->iuv[j];
      if (i2 < i1) {
        i0 = i2;
        i2 = i1;
        i1 = i0;
      }
      map[i2-1] = i1;
    }
  }
  for (nvrt = i = 0; i < m; i++)
    if (map[i] == 0) {
      nvrt++;
      map[i] = nvrt;
    } else {
      map[i] = -map[i];
    }
  /* update the vertex indices */
  for (i = 0; i < btess->nEdge; i++) {
    if (edgsen[3*i+2] != NULL) continue;
    if (edgsen[3*i  ] != NULL) {
      for (j = 0; j < edgsen[3*i  ]->npts; j++) {
        k = edgsen[3*i  ]->iuv[j];
        if (k == 0) {
          printf(" EGADS Warning: 0 in Edge- %d Mapping %d/%d (makeEFace)!\n",
                 i+1, j+1, edgsen[3*i  ]->npts);
          continue;
        }
        if (map[k-1] < 0) k = -map[k-1];
        edgsen[3*i  ]->iuv[j] = map[k-1];
      }
    }
    if (edgsen[3*i+1] == NULL) continue;
    for (j = 0; j < edgsen[3*i+1]->npts; j++) {
      k = edgsen[3*i+1]->iuv[j];
      if (k == 0) {
        printf(" EGADS Warning: 0 in Edge+ %d Mapping %d/%d (makeEFace)!\n",
               i+1, j+1, edgsen[3*i+1]->npts);
        continue;
      }
      if (map[k-1] < 0) k = -map[k-1];
      edgsen[3*i+1]->iuv[j] = map[k-1];
    }
  }

  for (ntri = i = 0; i < nFace; i++) {
    j = EG_indexBodyTopo(body, Faces[i]);
    ntri += btess->tess2d[j-1].ntris;
  }
  tris = (int *) EG_alloc(4*ntri*sizeof(int));
  if (tris == NULL) {
    printf(" EGADS Error: Malloc of %d triangles (EG_makeEFace)!\n", ntri);
    EG_free(edgsen);
    EG_free(map);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  itri = &tris[3*ntri];
  for (ntri = m = i = 0; i < nFace; i++) {
    j = EG_indexBodyTopo(body, Faces[i]);
    for (k = 0; k < btess->tess2d[j-1].ntris; k++, ntri++) {
      itri[ntri] = i+1;
      index = btess->tess2d[j-1].tris[3*k  ] + m;
      while (map[index-1] < 0) index = -map[index-1];
      tris[3*ntri  ] = map[index-1];
      index = btess->tess2d[j-1].tris[3*k+1] + m;
      while (map[index-1] < 0) index = -map[index-1];
      tris[3*ntri+1] = map[index-1];
      index = btess->tess2d[j-1].tris[3*k+2] + m;
      while (map[index-1] < 0) index = -map[index-1];
      tris[3*ntri+2] = map[index-1];
    }
    m += btess->tess2d[j-1].npts;
  }

  xyzs = (double *) EG_alloc(3*nvrt*sizeof(double));
  if (xyzs == NULL) {
    printf(" EGADS Error: Malloc of %d points (EG_makeEFace)!\n", nvrt);
    EG_free(tris);
    EG_free(edgsen);
    EG_free(map);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  for (nvrt = m = i = 0; i < nFace; i++) {
    j = EG_indexBodyTopo(body, Faces[i]);
    for (k = 0; k < btess->tess2d[j-1].npts; k++) {
      if (map[k+m] < 0) continue;
      xyzs[3*nvrt  ] = btess->tess2d[j-1].xyz[3*k  ];
      xyzs[3*nvrt+1] = btess->tess2d[j-1].xyz[3*k+1];
      xyzs[3*nvrt+2] = btess->tess2d[j-1].xyz[3*k+2];
      nvrt++;
    }
    m += btess->tess2d[j-1].npts;
  }
#ifdef DEBUG
  printf(" New triangulation: npts = %d,  ntris = %d!\n", nvrt, ntri);
#endif
  EG_free(edgsen);
  EG_free(map);
#ifdef DEBUG
  /* write out the triangulation
  {
    FILE *fp;
    fp = fopen("composite.dat", "w");
    if (fp == NULL) {
      printf(" CANNOT OPEN composite.dat!\n");
    } else {
      fprintf(fp, "1\n");
      fprintf(fp, "1\n");
      fprintf(fp, " %d %d\n", nvrt, ntri);
      for (i = 0; i < nvrt; i++)
        fprintf(fp, "%20.13le %20.13le %20.13le\n",
                xyzs[3*i  ], xyzs[3*i+1], xyzs[3*i+2]);
      for (i = 0; i < ntri; i++)
        fprintf(fp, "%d %d %d\n", tris[3*i  ], tris[3*i+1], tris[3*i+2]);
      fclose(fp);
    }
  }  */
#endif
/*
  for (i = 0; i < ntri; i++)
    if ((tris[3*i  ] < 1) || (tris[3*i  ] > nvrt) ||
        (tris[3*i+1] < 1) || (tris[3*i+1] > nvrt) ||
        (tris[3*i+2] < 1) || (tris[3*i+2] > nvrt))
      printf(" Bad tri %d/%d = %d %d %d [%d]\n", i, ntri, tris[3*i  ],
             tris[3*i+1], tris[3*i+2], nvrt);
 */
  stat = EG_uvmapMake(ntri, tris, itri, nvrt, xyzs, range, &trmap, &uvmap);
  EG_free(xyzs);
  EG_free(tris);
  if ((stat != EGADS_SUCCESS) || (range[0]*0.0 != 0.0)) {
    printf(" EGADS Error: EG_uvmapMake = %d (EG_makeEFace)!\n", stat);
    if (stat == EGADS_SUCCESS) {
      printf("              range = %lf %lf  %lf %lf\n", range[0], range[1],
                                                         range[2], range[3]);
      if (trmap != NULL) EG_free(trmap);
      uvmap_struct_free(uvmap);
      stat = EGADS_DEGEN;
    }
    EG_free(iedge);
    EG_free(efaces);
    goto restore;
  }
  
  /* fill the new eface structure */
  eface = (egEFace *) EG_alloc(sizeof(egEFace));
  if (eface == NULL) {
    printf(" EGADS Error: Malloc on EFace blind (EG_makeEFace)!\n");
    if (trmap != NULL) EG_free(trmap);
    uvmap_struct_free(uvmap);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  eface->eloops.nobjs = 0;
  eface->eloops.objs  = NULL;
  eface->senses       = NULL;
  eface->trmap        = trmap;
  eface->uvmap        = uvmap;
  eface->range[0]     = range[0];
  eface->range[1]     = range[1];
  eface->range[2]     = range[2];
  eface->range[3]     = range[3];
  eface->last         = 0;
  eface->sedges       = ebody->edges;
  eface->npatch       = nFace;
  eface->patches      = (egEPatch *) EG_alloc(nFace*sizeof(egEPatch));
  if (eface->patches == NULL) {
    printf(" EGADS Error: Malloc on %d Patches (EG_makeEFace)!\n", nFace);
    if (eface->trmap != NULL) EG_free(eface->trmap);
    EG_free(eface);
    uvmap_struct_free(uvmap);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  for (j = i = 0; i < nFace; i++) {
    ef = (egEFace *) efaces[i]->blind;
    eface->patches[i] = ef->patches[0];
    eface->patches[i].start = j;
    j += eface->patches[i].ntris;
  }

  /* build the EEdge fragment list */
  frags = (eloopFrag *) EG_alloc(n*sizeof(eloopFrag));
  if (frags == NULL) {
    printf(" EGADS Error: Malloc on %d Fragements (EG_makeEFace)!\n", n);
    EG_free(eface->patches);
    if (eface->trmap != NULL) EG_free(eface->trmap);
    EG_free(eface);
    uvmap_struct_free(uvmap);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  for (n = i = 0; i < neloop; i++) {
    if (eloops[i] == NULL) continue;
    eloop = (egELoop *) eloops[i]->blind;
    j0n = -1;
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      if (eloop->eedges.objs[j] == NULL) continue;
      eedge = (egEEdge *) eloop->eedges.objs[j]->blind;
      index = EG_indexBodyTopo(EBody, eloop->eedges.objs[j]);
      if (abs(iedge[index-1]) == 1) {
        frags[n].eedge = eloop->eedges.objs[j];
        frags[n].sense = eloop->senses[j];
        frags[n].iloop = i;
        frags[n].iedge = j;
        frags[n].newel = -1;
        frags[n].prev  = -1;
        if (j == 0) {
          j0n   = n;
        } else {
          index = EG_indexBodyTopo(EBody, eloop->eedges.objs[j-1]);
          if (abs(iedge[index-1]) == 1) frags[n].prev = n - 1;
        }
        if (eloop->senses[j] > 0) {
          frags[n].bnode = EG_indexBodyTopo(body, eedge->nodes[0]);
          frags[n].enode = EG_indexBodyTopo(body, eedge->nodes[1]);
        } else {
          frags[n].bnode = EG_indexBodyTopo(body, eedge->nodes[1]);
          frags[n].enode = EG_indexBodyTopo(body, eedge->nodes[0]);
        }
        frags[n].next = -1;
        if (j == eloop->eedges.nobjs-1) {
          if (j0n != -1) {
            frags[n].next   = j0n;
            frags[j0n].prev = n;
          }
        } else {
          index = EG_indexBodyTopo(EBody, eloop->eedges.objs[j+1]);
          if (abs(iedge[index-1]) == 1) frags[n].next = n + 1;
        }
        n++;
      }
    }
  }
  /* check for proper node connections */
  stat = EG_getBodyTopos(body, NULL, NODE, &nvrt, NULL);
  if (stat != EGADS_SUCCESS) {
    EG_free(frags);
    EG_free(iedge);
    EG_free(eface->patches);
    if (eface->trmap != NULL) EG_free(eface->trmap);
    EG_free(eface);
    uvmap_struct_free(uvmap);
    EG_free(efaces);
    goto restore;
  }
  map = (int *) EG_alloc(nvrt*sizeof(int));
  if (map == NULL) {
    printf(" EGADS Error: Malloc of %d Node mapping (EG_makeEFace)!\n",
           nvrt);
    EG_free(frags);
    EG_free(eface->patches);
    if (eface->trmap != NULL) EG_free(eface->trmap);
    EG_free(eface);
    uvmap_struct_free(uvmap);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  for (i = 0; i < nvrt; i++) map[i] = 0;
  for (i = 0; i < n; i++) {
    if (frags[i].prev == -1) map[frags[i].bnode-1]++;
    if (frags[i].next == -1) map[frags[i].enode-1]++;
  }
  for (j = 0; j < nvrt; j++)
    if ((map[j] != 0) && (map[j] != 2)) {
      printf(" EGADS Error: Incorrect Node count %d = %d (EG_makeEFace)!\n",
             j+1, map[j]);
      for (i = 0; i < n; i++) {
        printf(" %3d: ELoop/EEdge = %d/%d  Nodes = %3d %3d   prev = %3d   next = %3d\n",
               i, frags[i].iloop, frags[i].iedge, frags[i].bnode, frags[i].enode,
               frags[i].prev, frags[i].next);
      }
      EG_free(map);
      EG_free(frags);
      EG_free(eface->patches);
      if (eface->trmap != NULL) EG_free(eface->trmap);
      EG_free(eface);
      uvmap_struct_free(uvmap);
      EG_free(iedge);
      EG_free(efaces);
      stat = EGADS_TOPOERR;
      goto restore;
    }
  EG_free(map);
  
  /* build the new ELoops */
  for (i = 0; i < n; i++)
    if (frags[i].prev == -1) {
      for (j = 0; j < n; j++)
        if (frags[j].next == -1)
          if (frags[i].bnode == frags[j].enode) {
            frags[i].prev = j;
            frags[j].next = i;
            break;
          }
    }
  for (j = 0; j < n; j++)
    if ((frags[j].next == -1) || (frags[j].next == -1)) {
      printf(" EGADS Error: ELoop not closed (EG_makeEFace)!\n");
      for (i = 0; i < n; i++) {
        printf(" %3d: ELoop/EEdge = %d/%d  Nodes = %3d %3d   prev = %3d   next = %3d\n",
               i, frags[i].iloop, frags[i].iedge, frags[i].bnode, frags[i].enode,
               frags[i].prev, frags[i].next);
      }
      EG_free(frags);
      EG_free(eface->patches);
      if (eface->trmap != NULL) EG_free(eface->trmap);
      EG_free(eface);
      uvmap_struct_free(uvmap);
      EG_free(iedge);
      EG_free(efaces);
      stat = EGADS_TOPOERR;
      goto restore;
    }
  newel = 0;
  for (i = 0; i < n; i++) {
    if (frags[i].newel != -1) continue;
    frags[i].newel = newel;
    j = frags[i].next;
    while (frags[j].newel == -1) {
      frags[j].newel = newel;
      j = frags[j].next;
    }
    newel++;
  }
  newels = (egELoop **) EG_alloc(newel*sizeof(egELoop *));
  if (newels == NULL) {
    printf(" EGADS Error: Malloc on %d ELoop list (EG_makeEFace)!\n", newel);
    EG_free(frags);
    EG_free(eface->patches);
    if (eface->trmap != NULL) EG_free(eface->trmap);
    EG_free(eface);
    uvmap_struct_free(uvmap);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  for (i = 0; i < newel; i++) {
    newels[i] = (egELoop *) EG_alloc(sizeof(egELoop));
    if (newels[i] == NULL) {
      printf(" EGADS Error: Malloc on new ELoop %d (EG_makeEFace)!\n", i+1);
      for (j = 0; j < i; j++) EG_free(newels[j]);
      EG_free(newels);
      EG_free(frags);
      EG_free(eface->patches);
      if (eface->trmap != NULL) EG_free(eface->trmap);
      EG_free(eface);
      uvmap_struct_free(uvmap);
      EG_free(iedge);
      EG_free(efaces);
      stat = EGADS_MALLOC;
      goto restore;
    }
    newels[i]->eedges.nobjs = 0;
    newels[i]->eedges.objs  = NULL;
    newels[i]->senses       = NULL;
    newels[i]->area         = 0.0;
    newels[i]->nedge        = 0;
    newels[i]->edgeUVs      = NULL;
    for (k = j = 0; j < n; j++)
      if (frags[j].newel == i) k++;
#ifdef DEBUG
    printf(" New ELoop %d has %d EEdges\n", i+1, k);
#endif
    newels[i]->senses = (int *) EG_alloc(k*sizeof(int));
    if (newels[i]->senses == NULL) {
      printf(" EGADS Error: Malloc on new ELoop %d sense (EG_makeEFace)!\n",
             i+1);
      for (j = 0; j < i; j++) {
        EG_free(newels[j]->edgeUVs);
        EG_free(newels[j]->eedges.objs);
        EG_free(newels[j]->senses);
        EG_free(newels[j]);
      }
      EG_free(newels);
      EG_free(frags);
      EG_free(eface->patches);
      if (eface->trmap != NULL) EG_free(eface->trmap);
      EG_free(eface);
      uvmap_struct_free(uvmap);
      EG_free(iedge);
      EG_free(efaces);
      stat = EGADS_MALLOC;
      goto restore;
    }
    newels[i]->eedges.objs = (egObject **) EG_alloc(k*sizeof(egObject *));
    if (newels[i]->eedges.objs == NULL) {
      printf(" EGADS Error: Malloc on new ELoop %d objs (EG_makeEFace)!\n",
             i+1);
      for (j = 0; j < i; j++) {
        EG_free(newels[j]->edgeUVs);
        EG_free(newels[j]->eedges.objs);
        EG_free(newels[j]->senses);
        EG_free(newels[j]);
      }
      EG_free(newels);
      EG_free(frags);
      EG_free(eface->patches);
      if (eface->trmap != NULL) EG_free(eface->trmap);
      EG_free(eface);
      uvmap_struct_free(uvmap);
      EG_free(iedge);
      EG_free(efaces);
      stat = EGADS_MALLOC;
      goto restore;
    }
    newels[i]->eedges.nobjs = k;
    for (k = j = 0; j < n; j++) {
      if (frags[j].newel == i) {
        m = j;
        while (k != newels[i]->eedges.nobjs) {
          newels[i]->eedges.objs[k] = frags[m].eedge;
          newels[i]->senses[k]      = frags[m].sense;
          m = frags[m].next;
          k++;
        }
        break;
      }
      if (k == newels[i]->eedges.nobjs) break;
    }
    /* find all of the Edges */
    map = (int *) EG_alloc(2*btess->nEdge*sizeof(int));
    if (map == NULL) {
      printf(" EGADS Error: Malloc of %d EDGE mapping (EG_makeEFace)!\n",
             btess->nEdge);
      for (j = 0; j < i; j++) {
        EG_free(newels[j]->edgeUVs);
        EG_free(newels[j]->eedges.objs);
        EG_free(newels[j]->senses);
        EG_free(newels[j]);
      }
      EG_free(newels);
      EG_free(frags);
      EG_free(eface->patches);
      if (eface->trmap != NULL) EG_free(eface->trmap);
      EG_free(eface);
      uvmap_struct_free(uvmap);
      EG_free(iedge);
      EG_free(efaces);
      stat = EGADS_MALLOC;
      goto restore;
    }
    for (j = 0; j < 2*btess->nEdge; j++) map[j] = 0;
    for (j = 0; j < newels[i]->eedges.nobjs; j++) {
      if (newels[i]->eedges.objs[j] == NULL) continue;
      eedge = (egEEdge *) newels[i]->eedges.objs[j]->blind;
      for (k = 0; k < eedge->nsegs; k++) {
        index = eedge->segs[k].iedge;
        if (eedge->segs[k].sense*newels[i]->senses[j] == -1) {
          map[2*index-2]++;
        } else {
          map[2*index-1]++;
        }
      }
    }
    for (k = j = 0; j < 2*btess->nEdge; j++) if (map[j] != 0) k++;
#ifdef DEBUG
    printf("   nEdges/Senses = %d\n", k);
#endif
    newels[i]->nedge   = k;
    newels[i]->edgeUVs = (egEdgeUV *) EG_alloc(k*sizeof(egEdgeUV));
    if (newels[i]->edgeUVs == NULL) {
      printf(" EGADS Error: Malloc of %d edgeUVs (EG_makeEFace)!\n", k);
      EG_free(map);
      for (j = 0; j < i; j++) {
        EG_free(newels[j]->edgeUVs);
        EG_free(newels[j]->eedges.objs);
        EG_free(newels[j]->senses);
        EG_free(newels[j]);
      }
      EG_free(newels);
      EG_free(frags);
      EG_free(eface->patches);
      if (eface->trmap != NULL) EG_free(eface->trmap);
      EG_free(eface);
      uvmap_struct_free(uvmap);
      EG_free(iedge);
      EG_free(efaces);
      stat = EGADS_MALLOC;
      goto restore;
    }
    for (m = j = 0; j < neloop; j++) {
      if (eloops[j] == NULL) continue;
      eloop = (egELoop *) eloops[j]->blind;
      for (k = 0; k < eloop->nedge; k++) {
        index = EG_indexBodyTopo(body, eloop->edgeUVs[k].edge);
        if (index <= EGADS_SUCCESS) {
          printf(" EGADS Warning: indexBodyTopo = %d (EG_makeEFace)!\n",
                 index);
          continue;
        }
        if (eloop->edgeUVs[k].sense == -1) {
          index = 2*index-2;
        } else {
          index = 2*index-1;
        }
        if (map[index] != 0) {
          newels[i]->edgeUVs[m] = eloop->edgeUVs[k];
          m++;
          map[index] = 0;
        }
      }
    }
    if (m != newels[i]->nedge)
      printf(" EGADS Internal: EdgeUVs fill mismatch = %d %d\n",
             m, newels[i]->nedge);
    EG_free(map);
    EG_setLoopArea(uvmap, newels[i]);
  }
/*
  for (i = 0; i < n; i++) {
    printf(" %3d: ELoop/EEdge = %d/%d  Nodes = %3d %3d   prev/next = %3d %3d  %d\n",
           i, frags[i].iloop, frags[i].iedge, frags[i].bnode, frags[i].enode,
           frags[i].prev, frags[i].next, frags[i].newel);
  }  */
  EG_free(frags);
  
  /* get ELoop storage for the new EFace */
  eface->eloops.nobjs = newel;
  eface->senses       = (int *) EG_alloc(newel*sizeof(int));
  eface->eloops.objs  = (egObject **) EG_alloc(newel*sizeof(egObject *));
  if ((eface->senses == NULL) || (eface->eloops.objs == NULL)) {
    printf(" EGADS Error: Malloc of %d ELoop Objs (EG_makeEFace)!\n", newel);
    if (eface->senses != NULL) EG_free(eface->senses);
    if (eface->eloops.objs != NULL) EG_free(eface->eloops.objs);
    for (j = 0; j < newel; j++) {
      EG_free(newels[j]->edgeUVs);
      EG_free(newels[j]->eedges.objs);
      EG_free(newels[j]->senses);
      EG_free(newels[j]);
    }
    EG_free(newels);
    EG_free(eface->patches);
    if (eface->trmap != NULL) EG_free(eface->trmap);
    EG_free(eface);
    uvmap_struct_free(uvmap);
    EG_free(iedge);
    EG_free(efaces);
    stat = EGADS_MALLOC;
    goto restore;
  }
  if (newel == 1) {
    eface->senses[0] = 1;
    if (newels[0]->area <= 0.0)
      printf(" EGADS Internal: single loop area = %le (EG_makeEFace)!\n",
             newels[0]->area);
  } else {
    /* specify inner/outer flag */
    big  = 0.0;
    if (newels[0]->area > 0.0) big = newels[0]->area;
    biga = fabs(newels[0]->area);
    j    = 0;
    for (i = 1; i < newel; i++) {
      eface->senses[i] = 0;
      if (     newels[i]->area  > big)  big = newels[i]->area;
      if (fabs(newels[i]->area) > biga) {
        biga = fabs(newels[i]->area);
        j    = i;
      }
    }
    if (big != biga) {
      printf(" EGADS Internal: multi-loop areas = %le %le(EG_makeEFace)!\n",
             big, biga);
      for (i = 0; i < newel; i++)
        printf("                 area = %le\n", newels[i]->area);
    }
    for (i = 0; i < newel; i++)
      if (i == j) {
        eface->senses[i] = SFORWARD;
      } else {
        eface->senses[i] = SREVERSE;
      }
    stat = EG_coarseTris(newel, newels, uvmap);
    if (stat != EGADS_SUCCESS) {
      EG_free(eface->senses);
      EG_free(eface->eloops.objs);
      for (j = 0; j < newel; j++) {
        EG_free(newels[j]->edgeUVs);
        EG_free(newels[j]->eedges.objs);
        EG_free(newels[j]->senses);
        EG_free(newels[j]);
      }
      EG_free(newels);
      EG_free(eface->patches);
      if (eface->trmap != NULL) EG_free(eface->trmap);
      EG_free(eface);
      uvmap_struct_free(uvmap);
      EG_free(iedge);
      EG_free(efaces);
      goto restore;
    }
  }

  /* remove edgeUV storage from old ELoops w/ removed EEdges */
  for (i = 0; i < neloop; i++) {
    if (eloops[i] == NULL) continue;
    eloop = (egELoop *) eloops[i]->blind;
    for (j = 0; j < eloop->eedges.nobjs; j++) {
      if (eloop->eedges.objs[j] == NULL) continue;
      index = EG_indexBodyTopo(EBody, eloop->eedges.objs[j]);
      if (index <= EGADS_SUCCESS) continue;
      if (iedge[index-1] <= 1)    continue;
      eedge = (egEEdge *) eloop->eedges.objs[j]->blind;
      for (m = 0; m < eedge->nsegs; m++)
        for (k = 0; k < eloop->nedge; k++) {
          ie = eedge->segs[m].iedge - 1;
          if (eloop->edgeUVs[k].edge == eedge->sedges[ie].edge) {
            EG_free(eloop->edgeUVs[k].iuv);
            eloop->edgeUVs[k].iuv = NULL;
            break;
          }
        }
    }
  }
  /* remove EEdges */
  for (i = 0; i < ebody->eedges.nobjs; i++) {
    if (iedge[i] <= 1) continue;
    eedge = (egEEdge *) ebody->eedges.objs[i]->blind;
    EG_free(eedge->segs);
    EG_free(eedge);
    ebody->eedges.objs[i]->blind = NULL;
    EG_dereference(ebody->eedges.objs[i], 1);
    ebody->eedges.objs[i] = NULL;
  }
  for (k = i = 0; i < ebody->eedges.nobjs; i++) {
    if (ebody->eedges.objs[i] == NULL) continue;
    ebody->eedges.objs[k] = ebody->eedges.objs[i];
    k++;
  }
#ifdef DEBUG
  printf(" New Number of EEdges = %d (%d)\n", k, ebody->eedges.nobjs);
#endif
  ebody->eedges.nobjs = k;
  EG_free(iedge);
  /* remove old ELoops & insert new */
  if (neloop < newel)
    printf(" EGADS Internal: Old # ELoops = %d, New = %d (EG_makeEFace)!\n",
           neloop, newel);
  for (j = i = 0; i < neloop; i++) {
    if (eloops[i] == NULL) continue;
    eloop = (egELoop *) eloops[i]->blind;
    EG_free(eloop->senses);
    EG_free(eloop->edgeUVs);
    EG_free(eloop->eedges.objs);
    EG_free(eloop);
    index = EG_indexBodyTopo(EBody, eloops[i]);
    if (index <= EGADS_SUCCESS) {
      printf(" EGADS Internal: Cant find ELoop in EBody %d (EG_makeEFace)!\n",
             index);
      continue;
    }
    if (j < newel) {
      EG_attributeDel(ebody->eloops.objs[index-1], NULL);
      ebody->eloops.objs[index-1]->blind = newels[i];
      eface->eloops.objs[j] = ebody->eloops.objs[index-1];
      j++;
    } else {
      ebody->eloops.objs[index-1]->blind = NULL;
      EG_dereference(ebody->eloops.objs[index-1], 1);
    }
  }
  if (j != newel)
    printf(" EGADS Internal: # ELoop mismatch %d %d (EG_makeEFace)!\n",
           j, newel);
  EG_free(newels);

  /* remove single EFaces & make new EFace */
  for (k = i = 0; i < nFace; i++) {
    index = EG_indexBodyTopo(EBody, efaces[i]);
    if (index <= EGADS_SUCCESS) {
      printf(" EGADS Internal: Cant find EFace in EBody %d (EG_makeEFace)!\n",
             index);
      continue;
    }
    ef = (egEFace *) efaces[i]->blind;
    if (ef->trmap   != NULL) EG_free(ef->trmap);
    if (ef->uvmap   != NULL) uvmap_struct_free(ef->uvmap);
    if (ef->patches != NULL) EG_free(ef->patches);
    EG_free(ef->eloops.objs);
    if (ef->senses  != NULL) EG_free(ef->senses);
    EG_free(ef);
    efaces[i]->blind = NULL;
    if (k == 0) {
      ebody->efaces.objs[index-1]->blind = eface;
      ebody->efaces.objs[index-1]->mtype = SFORWARD;
      *EFace = ebody->efaces.objs[index-1];
      iface  = index-1;
      k++;
    } else {
      if (*EFace != NULL)
        if (fullAttr == 1) {
          EG_attributeDup(ebody->efaces.objs[index-1], *EFace);
        } else {
          EG_attributeCommon(ebody->efaces.objs[index-1], *EFace);
        }
      ebody->efaces.objs[index-1]->blind = NULL;
      EG_dereference(ebody->efaces.objs[index-1], 1);
    }
  }
  for (j = i = 0; i < ebody->efaces.nobjs; i++)
    if (ebody->efaces.objs[i]->blind != NULL) {
      ebody->efaces.objs[j] = ebody->efaces.objs[i];
      j++;
    }
#ifdef DEBUG
  printf(" New Number of EBody  EFaces = %d (%d)\n", j, ebody->efaces.nobjs);
#endif
  ebody->efaces.nobjs = j;
  for (j = i = 0; i < ebody->eloops.nobjs; i++)
    if (ebody->eloops.objs[i]->blind != NULL) {
      ebody->eloops.objs[j] = ebody->eloops.objs[i];
      j++;
    }
  ebody->eloops.nobjs = j;
  
  /* update our EShell */
  esh = (egEShell *) eshell->blind;
  for (j = i = 0; i < esh->efaces.nobjs; i++) {
    for (k = 0; k < nFace; k++)
      if (esh->efaces.objs[i] == efaces[k]) {
        esh->efaces.objs[j] = esh->efaces.objs[i];
        if (k == 0) j++;
        break;
      }
    if (k == nFace) {
      esh->efaces.objs[j] = esh->efaces.objs[i];
      j++;
    }
  }
  
#ifdef DEBUG
  printf(" New Number of EShell EFaces = %d (%d)\n", j, esh->efaces.nobjs);
  for (i = 0; i < eface->eloops.nobjs; i++) {
    if (eface->eloops.objs[i]        == NULL) continue;
    if (eface->eloops.objs[i]->blind == NULL) continue;
    printf(" ELoop %d   Outer/Inner = %d\n", i+1, eface->senses[i]);
    eloop = (egELoop *) eface->eloops.objs[i]->blind;
    printf("         %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[0]));
    for (k = 1; k < eloop->eedges.nobjs; k++)
      printf(" %d", EG_indexEBodyTopo(EBody, eloop->eedges.objs[k]));
    printf("\n");
    printf("         ");
    for (k = 0; k < eloop->eedges.nobjs; k++) {
      eedge = (egEEdge *) eloop->eedges.objs[k]->blind;
      if (eloop->senses[k] == SFORWARD) {
        printf("%d,%d ", EG_indexBodyTopo(body, eedge->nodes[0]),
                         EG_indexBodyTopo(body, eedge->nodes[1]));
      } else {
        printf("%d,%d ", EG_indexBodyTopo(body, eedge->nodes[1]),
                         EG_indexBodyTopo(body, eedge->nodes[0]));
      }
    }
    printf("\n");
  }
#endif
  esh->efaces.nobjs = j;
  EG_free(efaces);

  /* remove extraneous Nodes */
  EG_deleteObject(saved);
  stat = EG_copyEBody(EBody, &saved);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Warning: Second copyEBody = %d (EG_makeEFace)!\n", stat);
    return EGADS_SUCCESS;
  }
  stat = EG_eRemoveNode(EBody);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Info: EG_eRemoveNode = %d (EG_makeEFace)!\n", stat);
    EG_swapEBodies(EBody, saved);
    ebody  = (egEBody *) EBody->blind;
    *EFace = ebody->efaces.objs[iface];
  }
  EG_deleteObject(saved);
  
  return EGADS_SUCCESS;
  
restore:
  *EFace = NULL;
  EG_swapEBodies(EBody, saved);
  EG_deleteObject(saved);
  return stat;
}


int
EG_makeAttrEFaces(egObject *EBody, const char *attrName, int *nEFace,
                  /*@null@*/ egObject ***EFaces)
{
  int          i, j, k, m, n, stat, atype, alen, atx, alx, nFace, index, *mark;
  egEBody      *ebody;
  egEFace      *eface, *efacex;
  egObject     *body, *obj, **faces, *EFace, **efaces, **tmp;
  egTessel     *btess;
  const int    *ints,  *intx;
  const double *reals, *realx;
  const char   *str,   *strx;
  
  *nEFace = 0;
  efaces  = NULL;
  if (EFaces != NULL) *EFaces = NULL;

  if (EBody == NULL)               return EGADS_NULLOBJ;
  if (EBody->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EBody->oclass != EBODY)      return EGADS_NOTBODY;
  if (EG_sameThread(EBody))        return EGADS_CNTXTHRD;
  if (attrName == NULL)            return EGADS_NONAME;
 
  ebody = (egEBody *) EBody->blind;
  if (ebody == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_makeAttrEFaces)!\n");
    return EGADS_NOTFOUND;
  }
  if (ebody->done == 1) {
    printf(" EGADS Error: EBody finialized (EG_makeAttrEFaces)!\n");
    return EGADS_EXISTS;
  }
  btess = (egTessel *) ebody->ref->blind;
  body  = btess->src;
  faces = (egObject **) EG_alloc(btess->nFace*sizeof(egObject *));
  if (faces == NULL) {
    printf(" EGADS Error: Malloc on %d Faces (EG_makeAttrEFaces)!\n",
           btess->nFace);
    return EGADS_MALLOC;
  }
  mark = (int *) EG_alloc(btess->nFace*sizeof(int));
  if (mark == NULL) {
    printf(" EGADS Error: Malloc on %d Marks (EG_makeAttrEFaces)!\n",
           btess->nFace);
    EG_free(faces);
    return EGADS_MALLOC;
  }
  for (i = 0; i < btess->nFace; i++) mark[i] = 0;
  
  /* look for attributes */
  for (n = i = 0; i < ebody->efaces.nobjs-1; i++) {
    obj = ebody->efaces.objs[i];
    if (obj == NULL) continue;
    eface = (egEFace *) obj->blind;
    if (eface == NULL) continue;
    if (eface->npatch != 1) continue;
    index = EG_indexBodyTopo(body, eface->patches[0].face);
    if (index <= EGADS_SUCCESS) continue;
    if (mark[index-1] != 0) continue;
    stat = EG_attributeRet(eface->patches[0].face, attrName, &atype, &alen,
                           &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) {
      if (stat != EGADS_NOTFOUND)
        printf(" EGADS Warning: EG_attributeRet = %d (EG_makeAttrEFaces)!\n",
               stat);
      continue;
    }
    nFace = 0;
    for (j = i+1; j < ebody->efaces.nobjs; j++) {
      obj = ebody->efaces.objs[j];
      if (obj == NULL) continue;
      efacex = (egEFace *) obj->blind;
      if (efacex == NULL) continue;
      if (efacex->npatch != 1) continue;
      index = EG_indexBodyTopo(body, efacex->patches[0].face);
      if (index <= EGADS_SUCCESS) continue;
      if (mark[index-1] != 0) continue;
      stat = EG_attributeRet(efacex->patches[0].face, attrName, &atx, &alx,
                             &intx, &realx, &strx);
      if (stat != EGADS_SUCCESS) {
        if (stat != EGADS_NOTFOUND)
          printf(" EGADS Warning: EG_attributeRet = %d (EG_makeAttrEFaces)!\n",
                 stat);
        continue;
      }
      if (atx != atype) continue;
      if (alx != alen)  continue;
      if (atype == ATTRINT) {
        for (m = k = 0; k < alen; k++)
          if (ints[k] == intx[k]) m++;
        if (m != alen) continue;
      } else if ((atype == ATTRREAL) || (atype == ATTRCSYS)) {
        for (m = k = 0; k < alen; k++)
          if (reals[k] == realx[k]) m++;
        if (m != alen) continue;
      } else if (atype == ATTRPTR) {
        if (str != strx) continue;
      } else {
        if ((str != NULL) || (strx != NULL)) {
          if (str  == NULL) continue;
          if (strx == NULL) continue;
          if (strcmp(str,strx) != 0) continue;
        }
      }
      if (nFace == 0) {
        faces[0] = eface->patches[0].face;
        index    = EG_indexBodyTopo(body, faces[0]);
        if (index > EGADS_SUCCESS) {
          mark[index-1]++;
#ifdef DEBUG
          printf(" ...Using Face %d\n", index);
#endif
        } else {
          printf(" EGADS Warning: EG_indexBodyTopo = %d (EG_makeAttrEFaces)!\n",
                 stat);
        }
        nFace++;
      }
      faces[nFace] = efacex->patches[0].face;
      index        = EG_indexBodyTopo(body, faces[nFace]);
      if (index > EGADS_SUCCESS) {
        mark[index-1]++;
#ifdef DEBUG
        printf(" ...Using Face %d\n", index);
#endif
      } else {
        printf(" EGADS Warning: EG_indexBodyTopo = %d (EG_makeAttrEFaces)!\n",
               stat);
      }
      nFace++;
    }
    if (nFace == 0) continue;
#ifdef DEBUG
    printf(" ...total Faces = %d\n", nFace);
#endif
    stat = EG_makeEFace(EBody, nFace, faces, &EFace);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Warning: EG_makeEFace = %d (EG_makeAttrEFaces)!\n",
             stat);
      ebody = (egEBody *) EBody->blind;
      continue;
    }
    if (EFaces != NULL)
      if (n == 0) {
        efaces = (egObject **) EG_alloc(sizeof(egObject *));
        if (efaces == NULL) {
          printf(" EGADS Error: allocating return (EG_makeAttrEFaces)!\n");
          printf("              Effective Body Updated -- EFace return invalid\n");
        }
      } else {
        tmp = (egObject **) EG_reall(efaces, (n+1)*sizeof(egObject *));
        if (tmp == NULL) {
          printf(" EGADS Error: allocating %d return (EG_makeAttrEFaces)!\n",
                 n+1);
          printf("              Effective Body Updated -- EFace return invalid\n");
          EG_free(efaces);
          efaces = NULL;
          n      = 0;
        } else {
          efaces = tmp;
        }
      }
    *nEFace += 1;
    if (efaces != NULL) {
      efaces[n] = EFace;
      n++;
    }
  }
  EG_free(mark);
  EG_free(faces);
  if (*nEFace == 0) {
    if (efaces != NULL) EG_free(efaces);
    return EGADS_NOTFOUND;
  }
  
  if (EFaces != NULL) {
    *EFaces = efaces;
    *nEFace = n;
  }
  
  return EGADS_SUCCESS;
}


int
EG_effectiveMap(egObject *EObject, double *eparam, egObject **Object,
                double *param)
{
  int     stat, i1, i2, i3, iseg, ie, ix, itri, verts[3];
  double  tx, w[3];
  egEEdge *eedge;
  egEFace *eface;
  
  *Object = NULL;
  if  (EObject == NULL)               return EGADS_NULLOBJ;
  if  (EObject->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((EObject->oclass != EFACE) &&
      (EObject->oclass != EEDGE))     return EGADS_NOTTOPO;
  if  (EObject->blind == NULL)        return EGADS_NOTFOUND;
  
  if (EObject->oclass == EFACE) {
    
    eface = (egEFace *) EObject->blind;
    if (eface->npatch == 1) {
      *Object  = eface->patches[0].face;
      param[0] = eparam[0];
      param[1] = eparam[1];
      return EGADS_SUCCESS;
    }
    stat = EG_uvmapLocate(eface->uvmap, eface->trmap, eparam, &ix, &itri, verts,
                          w);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_uvmapLocate = %d (EG_effectiveMap)!\n", stat);
      return stat;
    }
    itri    -= eface->patches[ix-1].start + 1;
    *Object  = eface->patches[ix-1].face;
    i1       = eface->patches[ix-1].uvtris[3*itri  ] - 1;
    i2       = eface->patches[ix-1].uvtris[3*itri+1] - 1;
    i3       = eface->patches[ix-1].uvtris[3*itri+2] - 1;
    param[0] = w[0]*eface->patches[ix-1].uvs[2*i1  ] +
               w[1]*eface->patches[ix-1].uvs[2*i2  ] +
               w[2]*eface->patches[ix-1].uvs[2*i3  ];
    param[1] = w[0]*eface->patches[ix-1].uvs[2*i1+1] +
               w[1]*eface->patches[ix-1].uvs[2*i2+1] +
               w[2]*eface->patches[ix-1].uvs[2*i3+1];
  
  } else {
  
    eedge = (egEEdge *) EObject->blind;
    /* get the segment */
    for (iseg = 0; iseg < eedge->nsegs; iseg++)
      if (eparam[0] <= eedge->segs[iseg].tend) break;
    if (iseg == eedge->nsegs) iseg--;
    ie = eedge->segs[iseg].iedge - 1;
    /* get t in segment */
    tx = eparam[0] - eedge->segs[iseg].tstart;
    if (eedge->segs[iseg].sense == SREVERSE) {
      tx  = eedge->sedges[ie].ts[eedge->sedges[ie].npts-1] - tx;
    } else {
      tx += eedge->sedges[ie].ts[0];
    }
    *Object  = eedge->sedges[ie].edge;
    param[0] = tx;

  }
  
  return EGADS_SUCCESS;
}


int
EG_effectiveEdgeList(egObject *EEdge, int *nedges, egObject ***edges,
                     int **senses, double **tstart)
{
  int      iseg, ie, *msenses;
  double   *mtstart;
  egObject **medges;
  egEEdge  *eedge;
  
  *nedges = 0;
  *edges  = NULL;
  *senses = NULL;
  *tstart = NULL;
  if (EEdge == NULL)               return EGADS_NULLOBJ;
  if (EEdge->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EEdge->oclass != EEDGE)      return EGADS_NOTTOPO;
  if (EEdge->blind == NULL)        return EGADS_NOTFOUND;
  
  eedge   = (egEEdge *)   EEdge->blind;
  medges  = (egObject **) EG_alloc(eedge->nsegs*sizeof(egObject *));
  msenses = (int *)       EG_alloc(eedge->nsegs*sizeof(int));
  mtstart = (double *)    EG_alloc(eedge->nsegs*sizeof(double));
  if ((medges == NULL) || (msenses == NULL) || (mtstart == NULL)) {
    if (medges  != NULL) EG_free(medges);
    if (msenses != NULL) EG_free(msenses);
    if (mtstart != NULL) EG_free(mtstart);
    return EGADS_MALLOC;
  }
  
  for (iseg = 0; iseg < eedge->nsegs; iseg++) {
    ie            = eedge->segs[iseg].iedge - 1;
    medges[iseg]  = eedge->sedges[ie].edge;
    msenses[iseg] = eedge->segs[iseg].sense;
    mtstart[iseg] = eedge->segs[iseg].tstart;
  }
  *nedges = eedge->nsegs;
  *edges  = medges;
  *senses = msenses;
  *tstart = mtstart;

  return EGADS_SUCCESS;
}


int
EG_effectiveTri(egObject *EObj, double *uv, int *fIndex, int *itri, double *w)
{
  int      stat, ipat, i, verts[3];
  egObject *body;
  egEFace  *effect;
  
  w[0]    = w[1] = w[2] = 0.0;
  *fIndex = 0;
  *itri   = 0;
  if (EObj == NULL)               return EGADS_NULLOBJ;
  if (EObj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EObj->oclass != EFACE)      return EGADS_NOTTOPO;
  if (EObj->blind == NULL)        return EGADS_NOTFOUND;
  effect = (egEFace *) EObj->blind;
  
  if (effect->npatch == 1) {
    stat = EG_effectInTri(effect, uv, itri, w);
    if (stat != EGADS_SUCCESS) return stat;
    ipat = 0;
  } else {
    stat = EG_uvmapLocate(effect->uvmap, effect->trmap, uv, &ipat, &i, verts, w);
    if (stat != EGADS_SUCCESS) return stat;
    ipat--;
    *itri = i - effect->patches[ipat].start;
  }
  stat = EG_getBody(effect->patches[ipat].face, &body);
  if (stat != EGADS_SUCCESS) return stat;

  stat = EG_indexBodyTopo(body, effect->patches[ipat].face);
  if (stat > EGADS_SUCCESS) {
    *fIndex = stat;
    stat    = EGADS_SUCCESS;
  }
  
  return stat;
}


int
EG_getTessEFace(egObject *tess, int index, egObject **faces, double *uvs)
{
  int      i, i1, i2, i3, itri, ix, stat, verts[3];
  double   w[3];
  egTessel *btess;
  egObject *ebobj, *efobj;
  egEFace  *eface;
  
  if (tess == NULL)                          return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)            return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION)          return EGADS_NOTTESS;
  btess = (egTessel *) tess->blind;
  if (btess == NULL)                         return EGADS_NOTFOUND;
  ebobj = btess->src;
  if (ebobj == NULL)                         return EGADS_NULLOBJ;
  if (ebobj->magicnumber != MAGIC)           return EGADS_NOTOBJ;
  if (ebobj->oclass != EBODY)                return EGADS_EFFCTOBJ;
  if (btess->tess2d == NULL)                 return EGADS_NODATA;
  if ((index < 1) || (index > btess->nFace)) return EGADS_INDEXERR;

  for (i = 0; i < btess->tess2d[index-1].npts; i++) {
    faces[i] = NULL;
    uvs[2*i] = uvs[2*i+1] = 0.0;
  }
  stat = EG_objectEBodyTopo(ebobj, EFACE, index, &efobj);
  if (stat != EGADS_SUCCESS)                 return stat;
  eface = (egEFace  *) efobj->blind;
  if (eface == NULL)                         return EGADS_NULLOBJ;
  
  if (eface->npatch == 1) {
    for (i = 0; i < btess->tess2d[index-1].npts; i++) {
      faces[i]   = eface->patches[0].face;
      uvs[2*i  ] = btess->tess2d[index-1].uv[2*i  ];
      uvs[2*i+1] = btess->tess2d[index-1].uv[2*i+1];
    }
    return EGADS_SUCCESS;
  }
  
  for (i = 0; i < btess->tess2d[index-1].npts; i++) {
    stat = EG_uvmapLocate(eface->uvmap, eface->trmap,
                          &btess->tess2d[index-1].uv[2*i],
                          &ix, &itri, verts, w);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: %d/%d EG_uvmapLocate = %d (EG_getTessEFace)!\n",
             i+1, btess->tess2d[index-1].npts, stat);
      return stat;
    }
    itri      -= eface->patches[ix-1].start + 1;
    faces[i]   = eface->patches[ix-1].face;
    i1         = eface->patches[ix-1].uvtris[3*itri  ] - 1;
    i2         = eface->patches[ix-1].uvtris[3*itri+1] - 1;
    i3         = eface->patches[ix-1].uvtris[3*itri+2] - 1;
    uvs[2*i  ] = w[0]*eface->patches[ix-1].uvs[2*i1  ] +
                 w[1]*eface->patches[ix-1].uvs[2*i2  ] +
                 w[2]*eface->patches[ix-1].uvs[2*i3  ];
    uvs[2*i+1] = w[0]*eface->patches[ix-1].uvs[2*i1+1] +
                 w[1]*eface->patches[ix-1].uvs[2*i2+1] +
                 w[2]*eface->patches[ix-1].uvs[2*i3+1];
  }
  
  return EGADS_SUCCESS;
}
