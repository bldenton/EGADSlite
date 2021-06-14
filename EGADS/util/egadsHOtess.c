/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Functions to enhance tessellations for High-Order
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */
 
#include "egads.h"
#include <string.h>
#include <math.h>

//#define DEBUG
#define REPOSITION

#define MXSIDE           20
#define FUZZ             1.e-12  /* allowable negative weight */


  extern int  EG_getEdgeUVeval( const ego face, const ego topo, int sense,
                                double t, double *result );
#ifdef REPOSITION
  extern void EG_getSidepoint( const ego face, double fac, const double *uvm,
                               const double *uvp, /*@null@*/ const double *uvl,
                               /*@null@*/ const double *uvr, double *uv );
  extern void EG_getEdgepoint( const ego edge, double w, double tm, double tp,
                               double *tOUT );
  extern void EG_minArc4( const ego face, double fact1, double fact2,
                          const double *uv0, const double *uv1,
                          const double *uv2, const double *uv3, double *uvOUT );
  extern int  EG_baryInsert( const ego face, double w1, double w2, double w3,
                             const double *uvA, const double *uvB,
                             const double *uvC, double *uvOUT );
  extern void EG_getInterior( const ego face, const double *xyz, double *uv );
#endif


static void
EG_mdTFI(double xi, double et, int dim, const double *xl, const double *xu,
         const double *el,  const double *eu,  const double *uv0,
         const double *uv1, const double *uv2, const double *uv3,
         double *result)
{
  int i;
  
  for (i = 0; i < dim; i++) {
    result[i] = (1.0-xi)            * el[i]  +
                (    xi)            * eu[i]  +
                           (1.0-et) * xl[i]  +
                           (    et) * xu[i]  -
                (1.0-xi) * (1.0-et) * uv0[i] -
                (1.0-xi) * (    et) * uv3[i] -
                (    xi) * (1.0-et) * uv1[i] -
                (    xi) * (    et) * uv2[i];
  }
}


static int
EG_findSideIndex(int side, const double weight, int nst,
                 const int *type, const double *frac)
{
  int i, ret = -1;
  
  for (i = 0; i < nst; i++) {
    if (type[i] >= 0) continue;
    if (-type[i] != side+1) continue;
    if (fabs(frac[i]-weight) > FUZZ) continue;
    ret = i;
    break;
  }
    
  return ret;
}


#ifdef REPOSITION
static void
EG_correctEndPts(int i1, int i2, const int *degens, const int *iuv,
                 const double *uvs, const int *ptype, const int *pindex,
                 double *uvm, double *uvp)
{
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
}

#else

static void
EG_correctUV(double *uv, int i1, int i2, const int *degens, const int *iuv,
             const double *uvs, const int *ptype, const int *pindex)
{
  if (degens[0] == 0) return;
  
  if ((ptype[i1] == 0) && (pindex[i1] == degens[0])) {
    if (iuv[0] == 0) {
      uv[0] = uvs[2*i2  ];
    } else {
      uv[1] = uvs[2*i2+1];
    }
  } else if ((ptype[i1] == 0) && (pindex[i1] == degens[1])) {
    if (iuv[1] == 0) {
      uv[0] = uvs[2*i2  ];
    } else {
      uv[1] = uvs[2*i2+1];
    }
  } else if ((ptype[i2] == 0) && (pindex[i2] == degens[0])) {
    if (iuv[0] == 0) {
      uv[0] = uvs[2*i1  ];
    } else {
      uv[1] = uvs[2*i1+1];
    }
  } else if ((ptype[i2] == 0) && (pindex[i2] == degens[1])) {
    if (iuv[1] == 0) {
      uv[0] = uvs[2*i1  ];
    } else {
      uv[1] = uvs[2*i1+1];
    }
  }
}

#endif


static void
EG_correctUVq(double *u0, double *u1, double *u2, double *u3,
              int i0, int i1, int i2, int i3, const int *degens, const int *iuv,
              const int *ptype, const int *pindex)
{
  int i;
  
  if (degens[0] == 0) return;
  
  for (i = 0; i < 2; i++) {
    if (degens[i] == 0) return;
      if ((ptype[i0] == 0) && (pindex[i0] == degens[i])) {
        if (iuv[i] == 0) {
          u0[0] = 0.5*(u1[0]+u3[0]);
        } else {
          u0[1] = 0.5*(u1[1]+u3[1]);
        }
      } else if((ptype[i1] == 0) && (pindex[i1] == degens[i])) {
        if (iuv[i] == 0) {
          u1[0] = 0.5*(u0[0]+u2[0]);
        } else {
          u1[1] = 0.5*(u0[1]+u2[1]);
        }
      } else if((ptype[i2] == 0) && (pindex[i2] == degens[i])) {
        if (iuv[i] == 0) {
          u2[0] = 0.5*(u1[0]+u3[0]);
        } else {
          u2[1] = 0.5*(u1[1]+u3[1]);
        }
      } else if((ptype[i3] == 0) && (pindex[i3] == degens[i])) {
        if (iuv[i] == 0) {
          u3[0] = 0.5*(u2[0]+u0[0]);
        } else {
          u3[1] = 0.5*(u2[1]+u0[1]);
        }
      }
  }
}


static void
EG_evalEdgeSeg(const ego body, const ego face, const ego *edges,
               const egTessel *tessel, int ie, int i0, int i1, double weight,
               const double *uvs, const int *ptype, const int *pindex,
               const int *degens, double *uv)
{
  int    n, stat, oclass, mtype, pt0, pi0, pt1, pi1, nodes[2], *senses;
  double t, t0, t1, uvm[2], uvp[2], uvx[2], result[18];
  ego    geom, *objs;
  
  pt0    = ptype[i0];
  pi0    = pindex[i0];
  pt1    = ptype[i1];
  pi1    = pindex[i1];
  uvx[0] = uvs[2*i0  ];
  uvx[1] = uvs[2*i0+1];
  if ((pt0 == 0) && ((pi0 == degens[0]) || (pi0 == degens[1]))) {
    uvx[0] = uvs[2*i1  ];
    uvx[1] = uvs[2*i1+1];
  }
  
  stat = EG_getTopology(edges[ie], &geom, &oclass, &mtype, result, &n, &objs,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getTopo %d = %d (EG_evalEdgeSeg)!\n", ie+1, stat);
    return;
  }
  nodes[0] = nodes[1] = EG_indexBodyTopo(body, objs[0]);
  if (mtype == TWONODE) nodes[1] = EG_indexBodyTopo(body, objs[1]);
  
  if ((pt0 == 0) && (pt1 == 0)) {
    if (pi0 == nodes[0]) {
      t0 = tessel->tess1d[ie].t[0];
      t1 = tessel->tess1d[ie].t[1];
    } else {
      t0 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-1];
      t1 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-2];
    }
  } else if (pt0 == 0) {
    if (mtype == TWONODE) {
      if (pi0 == nodes[0]) {
        t0 = tessel->tess1d[ie].t[0];
        t1 = tessel->tess1d[ie].t[1];
      } else {
        t0 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-1];
        t1 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-2];
      }
    } else {
      if (pt1 == 2) {
        t0 = tessel->tess1d[ie].t[0];
        t1 = tessel->tess1d[ie].t[1];
      } else {
        t0 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-1];
        t1 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-2];
      }
    }
  } else if (pt1 == 0) {
    if (mtype == TWONODE) {
      if (pi1 == nodes[0]) {
        t1 = tessel->tess1d[ie].t[0];
        t0 = tessel->tess1d[ie].t[1];
      } else {
        t1 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-1];
        t0 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-2];
      }
    } else {
      if (pt0 == 2) {
        t0 = tessel->tess1d[ie].t[0];
        t1 = tessel->tess1d[ie].t[1];
      } else {
        t0 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-1];
        t1 = tessel->tess1d[ie].t[tessel->tess1d[ie].npts-2];
      }
    }
  } else if ((pi0 == ie+1) && (pi1 == ie+1)) {
    t0 = tessel->tess1d[ie].t[pt0-1];
    t1 = tessel->tess1d[ie].t[pt1-1];
  } else  {
     printf(" EGADS Info:  Iedge = %d   %d/%d   %d/%d (EG_evalEdgeSeg)!\n",
            ie+1, pt0, pi0, pt1, pi1);
    return;
  }

#ifdef REPOSITION
  if (t0 < t1) {
    EG_getEdgepoint(edges[ie],      weight,  t0, t1, &t);
  } else {
    EG_getEdgepoint(edges[ie], (1.0-weight), t1, t0, &t);
  }
#else
  t = (1.0-weight)*t0 + weight*t1;
#endif
  stat = EG_getEdgeUV(face, edges[ie], 0, t, uv);
  if (stat == EGADS_TOPOERR) {
    /* sense in Face twice! */
    uv[0] = uvx[0];
    uv[1] = uvx[1];
    stat = EG_getEdgeUV(face, edges[ie], -1, t, uvm);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Info:  EG_getEdgeUV -1 (EG_evalEdgeSeg)!\n");
    } else {
      stat = EG_getEdgeUV(face, edges[ie], 1, t, uvp);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Info:  EG_getEdgeUV +1\n (EG_evalEdgeSeg)!");
      } else {
        result[0] = sqrt((uvm[0]-uvx[0])*(uvm[0]-uvx[0]) +
                         (uvm[1]-uvx[1])*(uvm[1]-uvx[1]));
        result[1] = sqrt((uvp[0]-uvx[0])*(uvp[0]-uvx[0]) +
                         (uvp[1]-uvx[1])*(uvp[1]-uvx[1]));
        if (result[0] < result[1]) {
          uv[0] = uvm[0];
          uv[1] = uvm[1];
        } else {
          uv[0] = uvp[0];
          uv[1] = uvp[1];
        }
      }
    }
#ifdef DEBUG
    printf(" **** Double Edge %d: %lf %lf (%lf %lf) ***\n",
           ie+1, uv[0], uv[1], uvx[0], uvx[1]);
#endif
  }
  if (stat != EGADS_SUCCESS)
    printf(" EGADS Info:  EG_getEdgeUV = %d (EG_evalEdgeSeg)!", stat);
  
}


static void
EG_interiorTri(const ego body, const ego face, const ego *edges,
               const egTessel *tessel, const int *trs, const int *trc,
               const int *degens, const int *iuv,
               const double *uvs, const int *ptype, const int *pindex,
               double *w, double *uv0, double *uv1, double *uv2)
{
  int    ie;
  double dist, theta, suv[6];
#ifdef REPOSITION
  double uvm[2], uvp[2];
#endif
  
  /* side 0 */
  theta = w[1] + w[2];
  if (theta == 0.0) theta = 1.0;
  dist  = w[2]/theta;
  if (trc[0] > 0) {
#ifdef REPOSITION
    EG_correctEndPts(trs[1]-1, trs[2]-1, degens, iuv, uvs, ptype, pindex,
                     uvm, uvp);
    EG_getSidepoint(face, 1.0-dist, uvm, uvp, NULL, NULL, &suv[0]);
#else
    suv[0] = uv1[0]  + dist*(uv2[0]  - uv1[0]);
    suv[1] = uv1[1]  + dist*(uv2[1]  - uv1[1]);
    EG_correctUV(&suv[0], trs[1]-1, trs[2]-1, degens, iuv, uvs, ptype, pindex);
#endif
  } else {
    ie = -trc[0]-1;
    EG_evalEdgeSeg(body, face, edges, tessel, ie, trs[1]-1, trs[2]-1, 1.0-dist,
                   uvs, ptype, pindex, degens, &suv[0]);
  }
  
  /* side 1 */
  theta = atan2(w[2], 1.0-w[1]);
  dist  = w[2]    + w[1]*tan(theta);
  if (trc[1] > 0) {
#ifdef REPOSITION
    EG_correctEndPts(trs[0]-1, trs[2]-1, degens, iuv, uvs, ptype, pindex,
                     uvm, uvp);
    EG_getSidepoint(face, 1.0-dist, uvm, uvp, NULL, NULL, &suv[2]);
#else
    suv[2] = uv0[0]  + dist*(uv2[0]  - uv0[0]);
    suv[3] = uv0[1]  + dist*(uv2[1]  - uv0[1]);
    EG_correctUV(&suv[2], trs[0]-1, trs[2]-1, degens, iuv, uvs, ptype, pindex);
#endif
  } else {
    ie = -trc[1]-1;
    EG_evalEdgeSeg(body, face, edges, tessel, ie, trs[0]-1, trs[2]-1, 1.0-dist,
                   uvs, ptype, pindex, degens, &suv[2]);
  }
  
  /* side 2 */
  theta = atan2(w[1], 1.0-w[2]);
  dist  = w[1] + w[2]*tan(theta);
  if (trc[2] > 0) {
#ifdef REPOSITION
    EG_correctEndPts(trs[0]-1, trs[1]-1, degens, iuv, uvs, ptype, pindex,
                     uvm, uvp);
    EG_getSidepoint(face, 1.0-dist, uvm, uvp, NULL, NULL, &suv[4]);
#else
    suv[4] = uv0[0]  + dist*(uv1[0]  - uv0[0]);
    suv[5] = uv0[1]  + dist*(uv1[1]  - uv0[1]);
    EG_correctUV(&suv[4], trs[0]-1, trs[1]-1, degens, iuv, uvs, ptype, pindex);
#endif
  } else {
    ie = -trc[2]-1;
    EG_evalEdgeSeg(body, face, edges, tessel, ie, trs[0]-1, trs[1]-1, 1.0-dist,
                   uvs, ptype, pindex, degens, &suv[4]);
  }
  
  /* set up smaller side-based triangle */
  w[0]   = 0.5*(1.0 - w[0]);
  w[1]   = 0.5*(1.0 - w[1]);
  w[2]   = 0.5*(1.0 - w[2]);
  uv0[0] = suv[0];
  uv0[1] = suv[1];
  uv1[0] = suv[2];
  uv1[1] = suv[3];
  uv2[0] = suv[4];
  uv2[1] = suv[5];
}


static int
EG_fillEdgeSeg(const ego face, const ego *edges, const egTessel *tessel,
               int nIns, const int *corner, int ie, const double *xyzs,
               const double *uvs, const int *ptype, const int *pindex,
               int side, int nst, const int *type, const double *frac,
               const double *sinsert, const int *degens, int ielem, int *elems,
               int *np, double *coords, double *parms)
{
  int    i, j, k, stat, pt0, pi0, pt1, pi1, dir, in;
  double xyz0[3], xyz1[3], uv[2], uvm[2], uvp[2], uvx[2], result[2];
#ifdef DEBUG
  double d;
#endif
  
  if (nIns == 0) return EGADS_SUCCESS;
  i       = corner[0] - 1;
  xyz0[0] = xyzs[3*i  ];
  xyz0[1] = xyzs[3*i+1];
  xyz0[2] = xyzs[3*i+2];
  uvx[0]  = uvs[2*i  ];
  uvx[1]  = uvs[2*i+1];
  pt0     = ptype[i];
  pi0     = pindex[i];
  i       = corner[1] - 1;
  xyz1[0] = xyzs[3*i  ];
  xyz1[1] = xyzs[3*i+1];
  xyz1[2] = xyzs[3*i+2];
  pt1     = ptype[i];
  pi1     = pindex[i];
  if ((pt0 == 0) && ((pi0 == degens[0]) || (pi0 == degens[1]))) {
    uvx[0] = uvs[2*i  ];
    uvx[1] = uvs[2*i+1];
  }
  
  j = dir = 0;
  if ((pt0 == 0) && (pt1 == 0)) {
    if (pi0 == corner[2]) {
      dir =  1;
      j   =  0;
    } else {
      dir = -1;
      j   = tessel->tess1d[ie].npts-1;
    }
  } else if (pt0 == 0) {
    if (corner[2] == corner[3]) {
      k = nIns+1;
      if ((xyz1[0] == tessel->tess1d[ie].xyz[3*k  ]) &&
          (xyz1[1] == tessel->tess1d[ie].xyz[3*k+1]) &&
          (xyz1[2] == tessel->tess1d[ie].xyz[3*k+2])) {
        dir =  1;
        j   =  0;
      } else {
        dir = -1;
        j   = tessel->tess1d[ie].npts-1;
      }
    } else {
      if (pi0 == corner[2]) {
        dir =  1;
        j   =  0;
      } else {
        dir = -1;
        j   = tessel->tess1d[ie].npts-1;
      }
    }
  } else if (pt1 == 0) {
    j = (pt0-1)*(nIns+1);
    if (corner[2] == corner[3]) {
      k = nIns+1;
      if ((xyz0[0] == tessel->tess1d[ie].xyz[3*k  ]) &&
          (xyz0[1] == tessel->tess1d[ie].xyz[3*k+1]) &&
          (xyz0[2] == tessel->tess1d[ie].xyz[3*k+2])) {
        dir = -1;
      } else {
        dir =  1;
      }
    } else {
      dir = -1;
      if (pi1 == corner[3]) dir = 1;
    }
  } else if ((pi0 == ie+1) && (pi1 == ie+1)) {
    dir = 1;
    if (pt0 > pt1) dir = -1;
    j = (pt0-1)*(nIns+1);
  }
  if (dir == 0) return EGADS_INDEXERR;

#ifdef DEBUG
  d = sqrt((xyz0[0]-tessel->tess1d[ie].xyz[3*j  ])*
           (xyz0[0]-tessel->tess1d[ie].xyz[3*j  ]) +
           (xyz0[1]-tessel->tess1d[ie].xyz[3*j+1])*
           (xyz0[1]-tessel->tess1d[ie].xyz[3*j+1]) +
           (xyz0[2]-tessel->tess1d[ie].xyz[3*j+2])*
           (xyz0[2]-tessel->tess1d[ie].xyz[3*j+2]));
  if (d > FUZZ) {
    printf("    %lf %lf %lf  %lf %lf %lf   %d\n", xyz0[0], xyz0[1], xyz0[2],
           tessel->tess1d[ie].xyz[3*j  ], tessel->tess1d[ie].xyz[3*j+1],
           tessel->tess1d[ie].xyz[3*j+2], dir);
    return EGADS_INDEXERR;
  }
#endif

  j+=dir;
  for (i = 0; i < nIns; i++, j+=dir) {
    stat = EG_getEdgeUV(face, edges[ie], 0, tessel->tess1d[ie].t[j], uv);
    if (stat == EGADS_TOPOERR) {
      /* sense in Face twice! */
      uv[0] = uvx[0];
      uv[1] = uvx[1];
      stat = EG_getEdgeUV(face, edges[ie], -1, tessel->tess1d[ie].t[j], uvm);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Info:  EG_getEdgeUV -1 (EG_fillEdgeSeg)!\n");
      } else {
        stat = EG_getEdgeUV(face, edges[ie], 1, tessel->tess1d[ie].t[j], uvp);
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Info:  EG_getEdgeUV +1\n (EG_fillEdgeSeg)!");
        } else {
          result[0] = sqrt((uvm[0]-uvx[0])*(uvm[0]-uvx[0]) +
                           (uvm[1]-uvx[1])*(uvm[1]-uvx[1]));
          result[1] = sqrt((uvp[0]-uvx[0])*(uvp[0]-uvx[0]) +
                           (uvp[1]-uvx[1])*(uvp[1]-uvx[1]));
          if (result[0] < result[1]) {
            uv[0] = uvm[0];
            uv[1] = uvm[1];
          } else {
            uv[0] = uvp[0];
            uv[1] = uvp[1];
          }
        }
      }
#ifdef DEBUG
      printf(" **** Double Edge %d: %lf %lf (%lf %lf) %d  %d %d ***\n",
             ie+1, uv[0], uv[1], uvx[0], uvx[1], corner[0], pt0, pi0);
#endif
    }
    if (stat != EGADS_SUCCESS) return stat;
    coords[*np*3  ] = tessel->tess1d[ie].xyz[3*j  ];
    coords[*np*3+1] = tessel->tess1d[ie].xyz[3*j+1];
    coords[*np*3+2] = tessel->tess1d[ie].xyz[3*j+2];
    parms[*np*2  ]  = uv[0];
    parms[*np*2+1]  = uv[1];
    *np += 1;
    in = EG_findSideIndex(side, sinsert[i], nst, type, frac);
    if (in < 0) return EGADS_INDEXERR;
#ifdef DEBUG
    if (elems[ielem+in] != 0) printf(" double Hit!\n");
#endif
    elems[ielem+in] = *np;
  }
  
#ifdef DEBUG
  d = sqrt((xyz1[0]-tessel->tess1d[ie].xyz[3*j  ])*
           (xyz1[0]-tessel->tess1d[ie].xyz[3*j  ]) +
           (xyz1[1]-tessel->tess1d[ie].xyz[3*j+1])*
           (xyz1[1]-tessel->tess1d[ie].xyz[3*j+1]) +
           (xyz1[2]-tessel->tess1d[ie].xyz[3*j+2])*
           (xyz1[2]-tessel->tess1d[ie].xyz[3*j+2]));
  if (d > FUZZ) {
    printf("    BAD end %lf %lf %lf  %lf %lf %lf   %d\n", xyz1[0], xyz1[1], xyz1[2],
           tessel->tess1d[ie].xyz[3*j  ], tessel->tess1d[ie].xyz[3*j+1],
           tessel->tess1d[ie].xyz[3*j+2], dir);
    return EGADS_INDEXERR;
  }
#endif
  
  return EGADS_SUCCESS;
}


/*
 * builds a new tessellation object that inserts the High Order vertices based
 *         the specified internal positions
 *
 * where: tess    - the source tessellation object
 *        nst     - number of positions required for each HO triangle (contains 
 *                  linear corners (- is quads from a quad tessellation)
 *        nItri   - number of internal triangles created per source tri
 *                  (- for quads indicates the new tessellation is triangular)
 *        iTris   - the internal triangle indices (1-bias)
 *                  3*nItri in length
 *        st      - the weights that define the positions within the triangles 
 *                  2*nst in length
 *        newTess - the resultant tessellation object where every triangle is
 *                  replaced with nItri triangles (and in order)
 *
 */

int
EG_tessHOverts(const ego tess, int nstx, int nItrix, const int *iTris,
               const double *st, ego *nTess)
{
  int          i, j, k, n, stat, outLevel, nst, nItri, atype, alen, corner[4];
  int          i0, i1, i2, i3, nIns, nedges, nfaces, sum[2], degens[2], iuv[2];
  int          np, nt, npts, ntris, nside, nei, oclass, mtype, *senses, *type;
  int          nmid = 0, quad = 0, qout = 0, *elems = NULL, *tris = NULL;
  double       area, d, *parms, sinsert[MXSIDE], result[18], trange[2], uv[2];
  double       w[3], u0[2], u1[2], u2[2], u3[2], *frac = NULL, *coords = NULL;
#ifdef REPOSITION
  double       uvm[2], uvp[2], xyz[3];
  const double *uvl, *uvr;
#endif
  ego          body, context, geom, *objs, *nodes;
  ego          *edges = NULL, *faces = NULL, newTess = NULL;
  egTessel     *btess, *tessel;
  const int    *ints, *ptype, *pindex, *trs, *trc;
  const double *reals, *xyzs, *ts, *uvs;
  const char   *str;
  static int   sidet[3][2] = {{1,2}, {2,0}, {0,1}};
  static int   sideq[4][2] = {{1,2}, {2,5}, {5,0}, {0,1}};
  static int   neigq[4]    = { 0,     3,     4,     2   };

  *nTess = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  stat = EG_getContext(tess, &context);
  if (stat != EGADS_SUCCESS) return stat;
  outLevel = EG_setOutLevel(context, 1);
  if (outLevel < EGADS_SUCCESS) return outLevel;
  stat = EG_setOutLevel(context, outLevel);
  if (stat < EGADS_SUCCESS) return stat;

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_tessHOverts)!\n");
    return EGADS_NOTFOUND;
  }
  body = btess->src;
  if (body == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_tessHOverts)!\n");
    return EGADS_NULLOBJ;
  }
  if (body->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_tessHOverts)!\n");
    return EGADS_NOTOBJ;
  }
  if (body->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_tessHOverts)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess2d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_tessHOverts)!\n");
    return EGADS_NODATA;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_tessHOverts)!\n");
    return EGADS_NODATA;
  }
  nst   = nstx;
  nItri = abs(nItrix);
  if (nstx < 0) {
    nst  = -nstx;
    quad = 1;
    stat = EG_attributeRet(tess, ".tessType", &atype, &alen,
                           &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: attributeRet on tessType = %d (EG_tessHOverts)!\n",
               stat);
      return stat;
    }
    if (atype != ATTRSTRING) {
      if (outLevel > 0)
        printf(" EGADS Error: attribute Type on tessType (EG_tessHOverts)!\n");
      return EGADS_ATTRERR;
    }
    if (strcmp(str, "Quad") != 0) {
      if (outLevel > 0)
        printf(" EGADS Error: attribute tessType = %s (EG_tessHOverts)!\n", str);
      return EGADS_ATTRERR;
    }
    if (nItrix > 0) qout = 1;
  }
  
  /* make sure we don't have duplicate verts */
  for (i = 0; i < nst-1; i++)
    for (j = i+1; j < nst; j++) {
      if ((fabs(st[2*i  ]-st[2*j  ]) < FUZZ) &&
          (fabs(st[2*i+1]-st[2*j+1]) < FUZZ)) {
        if (outLevel > 0)
          printf(" EGADS Error: st %d matches st %d (EG_tessHOverts)!\n",
                 i+1, j+1);
        return EGADS_DEGEN;
      }
    }
  
  type = (int *) EG_alloc(nst*sizeof(int));
  if (type == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocating %d ints (EG_tessHOverts)!\n", nst);
    return EGADS_MALLOC;
  }
  frac = (double *) EG_alloc(nst*sizeof(double));
  if (frac == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocating %d doubles (EG_tessHOverts)!\n", nst);
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  
  /* make sure we have all "corner verts" */
  corner[0] = corner[1] = corner[2] = corner[3] = -1;
  for (i = 0; i < nst; i++) {
    type[i] = 0;
    frac[i] = 0.0;
    if ((fabs(    st[2*i  ]) < FUZZ) && (fabs(    st[2*i+1]) < FUZZ)) {
      corner[0] = i;
      type[i]   = 1;
    }
    if ((fabs(1.0-st[2*i  ]) < FUZZ) && (fabs(    st[2*i+1]) < FUZZ)) {
      corner[1] = i;
      type[i]   = 2;
    }
    if ((fabs(    st[2*i  ]) < FUZZ) && (fabs(1.0-st[2*i+1]) < FUZZ)) {
      corner[2] = i;
      type[i]   = 3;
    }
    if ((fabs(1.0-st[2*i  ]) < FUZZ) && (fabs(1.0-st[2*i+1]) < FUZZ)) {
      corner[3] = i;
      type[i]   = 4;
    }
  }
  if (quad == 1) {
    i         = corner[2];
    type[i]   = 4;
    corner[2] = corner[3];
    corner[3] = i;
    i         = corner[2];
    type[i]   = 3;
  }
  for (i = 0; i < 3+quad; i++)
    if (corner[i] == -1) {
      if (outLevel > 0)
        printf(" EGADS Error: Corner %d not found (EG_tessHOverts)!\n", i+1);
      stat = EGADS_NOTFOUND;
      goto cleanup;
    }
  
  /* check triangle area */
  area = 0;
  for (i = 0; i < nItri; i++) {
    i0 = iTris[3*i  ] - 1;
    i1 = iTris[3*i+1] - 1;
    i2 = iTris[3*i+2] - 1;
    d  = ((st[2*i0  ]-st[2*i2  ])*(st[2*i1+1]-st[2*i2+1]) -
          (st[2*i0+1]-st[2*i2+1])*(st[2*i1  ]-st[2*i2  ]));
    if (d <= 0.0) {
      if (outLevel > 0)
        printf(" EGADS Error: Triangle %d has area = %lf (EG_tessHOverts)!\n",
               i+1, d);
      stat = EGADS_INDEXERR;
      goto cleanup;
    }
    area += 0.5*d;
/*  printf(" area = %lf %lf\n", d, area);  */
  }
  if (quad == 0) area *= 2.0;
  if (fabs(area-1.0) > FUZZ) {
    if (outLevel > 0)
      printf(" EGADS Error: Area != 1.0 (%lf) (EG_tessHOverts)!\n", area);
    stat = EGADS_DEGEN;
    goto cleanup;
  }
  
  /* what is our side matching? */
  for (nIns = i = 0; i < nst; i++) {
    if ((corner[0] == i) || (corner[1] == i)) continue;
    if (fabs(st[2*i+1]) >= FUZZ) continue;
    if (nIns >= MXSIDE) {
      if (outLevel > 0)
        printf(" EGADS Error: Too many side insertions (%d) (EG_tessHOverts)!\n",
               nIns);
      stat = EGADS_INDEXERR;
      goto cleanup;
    }
    type[i] = -3;
    frac[i] = st[2*i  ];
    if (nIns == 0) {
      sinsert[nIns] = st[2*i  ];
    } else {
      i0 = nIns;
      for (j = 0; j < nIns; j++)
        if (sinsert[j] > st[2*i  ]) {
          i0 = j;
          break;
        }
      for (j = nIns-1; j >= i0; j--) sinsert[j+1] = sinsert[j];
      sinsert[i0] = st[2*i  ];
    }
    nIns++;
  }
#ifdef DEBUG
  for (i = 0; i < nIns; i++)
    printf(" side insert %d = %lf\n", i+1, sinsert[i]);
#endif
  /* test other sides */
  if (quad == 0) {
    for (i = 0; i < nst; i++) {
      if ((corner[0] == i) || (corner[2] == i)) continue;
      if (fabs(st[2*i  ]) >= FUZZ) continue;
      for (i0 = j = 0; j < nIns; j++)
        if (fabs(st[2*i+1]-sinsert[j]) < FUZZ) i0++;
      if (i0 == 1) {
        type[i] = -2;
        frac[i] = 1.0-st[2*i+1];
      } else {
        if (outLevel > 0)
          printf(" EGADS Error: Side 2 Bad insert (%lf) %d (EG_tessHOverts)!\n",
                 st[2*i+1], i0);
        stat = EGADS_INDEXERR;
        goto cleanup;
      }
    }
    for (i = 0; i < nst; i++) {
      if ((corner[1] == i) || (corner[2] == i)) continue;
      d = 1.0 - st[2*i  ] - st[2*i+1];
      if (fabs(d) >= FUZZ) {
        if ((fabs(st[2*i  ]) >= FUZZ) && (fabs(st[2*i+1]) >= FUZZ)) nmid++;
        continue;
      }
      for (i0 = j = 0; j < nIns; j++)
        if (fabs(st[2*i+1]-sinsert[j]) < FUZZ) i0++;
      if (i0 == 1) {
        type[i] = -1;
        frac[i] = st[2*i+1];
      } else {
        if (outLevel > 0)
          printf(" EGADS Error: Side 1 Bad insert (%lf) %d (EG_tessHOverts)!\n",
                 st[2*i+1], i0);
        stat = EGADS_INDEXERR;
        goto cleanup;
      }
    }
  } else {
    for (i = 0; i < nst; i++) {
      if ((corner[0] == i) || (corner[1] == i) || (corner[2] == i) ||
          (corner[3] == i)) continue;
      d = 0.0;
      if (fabs(st[2*i  ]    ) < FUZZ) {
        d       = 1.0-st[2*i+1];
        type[i] = -3;
      } else if (fabs(st[2*i  ]-1.0) < FUZZ) {
        d       =     st[2*i+1];
        type[i] = -1;
      } else if (fabs(st[2*i+1]    ) < FUZZ) {
        d       =     st[2*i  ];
        type[i] = -4;
      } else if (fabs(st[2*i+1]-1.0) < FUZZ) {
        d       = 1.0-st[2*i  ];
        type[i] = -2;
      } else {
        nmid++;
        continue;
      }
      frac[i] = d;
      for (i0 = j = 0; j < nIns; j++)
        if (fabs(d-sinsert[j]) < FUZZ) i0++;
      if (i0 != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Quad Bad insert (%lf,%lf) %d (EG_tessHOverts)!\n",
                 st[2*i  ], st[2*i+1], i0);
        stat = EGADS_INDEXERR;
        goto cleanup;
      }
    }
  }
  
  /* are we quad-in and quad-out -- if so we need pairs of triangles */
  if ((quad == 1) && (qout == 1)) {
    if (nItri%2 != 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Odd number of Tris = %d (EG_tessHOverts)!\n",
               nItri);
      stat = EGADS_INDEXERR;
      goto cleanup;
    }
    for (i = 0; i < nItri; i+=2) {
      for (i2 = i0 = 0; i0 < 3; i0++)
        for (i1 = 0; i1 < 3; i1++)
          if (iTris[3*i+i0] == iTris[3*i+i1+3]) i2++;
      if (i2 != 2) {
        if (outLevel > 0)
          printf(" EGADS Error: Triangles not in pairs = %d (EG_tessHOverts)!\n",
                 i+1);
        stat = EGADS_INDEXERR;
        goto cleanup;
      }
    }
  }
  
  /* initialize the new tessellation object */
  stat = EG_initTessBody(body, &newTess);
  if ((stat != EGADS_SUCCESS) || (newTess == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_initTessBody = %d (EG_tessHOverts)!\n", stat);
    if (stat == EGADS_SUCCESS) stat = EGADS_NULLOBJ;
    goto cleanup;
  }
  
  stat = EG_getBodyTopos(body, NULL, EDGE, &nedges, &edges);
  if ((stat != EGADS_SUCCESS) || (edges == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getBodyTopos E = %d (EG_tessHOverts)!\n", stat);
    if (stat == EGADS_SUCCESS) stat = EGADS_TOPOERR;
    goto cleanup;
  }
  
  /* rebuild the Edges */
  for (j = i = 0; i < nedges; i++) {
    if (edges[i]->mtype == DEGENERATE) continue;
    stat = EG_getTessEdge(tess, i+1, &npts, &xyzs, &ts);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTessEdge %d = %d (EG_tessHOverts)!\n",
               i+1, stat);
      goto cleanup;
    }
    if (npts == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTessEdge %d -- no points (EG_tessHOverts)!\n",
               i+1);
      stat = EGADS_INDEXERR;
      goto cleanup;
    }
    if (npts > j) j = npts;
  }

  /* allocate to the maximum length */
  coords = (double *) EG_alloc(4*((nIns+1)*j)*sizeof(double));
  if (coords == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d points (EG_tessHOverts)!\n", j);
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  parms = &coords[3*(nIns+1)*j];
  
  for (i = 0; i < nedges; i++) {
    if (edges[i]->mtype == DEGENERATE) continue;
    stat = EG_getTessEdge(tess, i+1, &npts, &xyzs, &ts);
    if (stat != EGADS_SUCCESS) continue;
    
    for (i0 = j = 0; j < npts-1; j++) {
      parms[i0]      = ts[j];
      coords[3*i0  ] = xyzs[3*j  ];
      coords[3*i0+1] = xyzs[3*j+1];
      coords[3*i0+2] = xyzs[3*j+2];
      i0++;
      for (k = 0; k < nIns; k++, i0++) {
        
#ifdef REPOSITION
        EG_getEdgepoint(edges[i], sinsert[k], ts[j], ts[j+1], &parms[i0]);
#else
        parms[i0] = (1.0-sinsert[k])*ts[j] + sinsert[k]*ts[j+1];
#endif
        stat = EG_evaluate(edges[i], &parms[i0], result);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_evaluate Edge %d/%d = %d (EG_tessHOverts)!\n",
                   i+1, j+1, stat);
          goto cleanup;
        }
        coords[3*i0  ] = result[0];
        coords[3*i0+1] = result[1];
        coords[3*i0+2] = result[2];
      }
    }
    j              = npts-1;
    parms[i0]      = ts[j];
    coords[3*i0  ] = xyzs[3*j  ];
    coords[3*i0+1] = xyzs[3*j+1];
    coords[3*i0+2] = xyzs[3*j+2];
    i0++;
    
    stat = EG_setTessEdge(newTess, i+1, i0, coords, parms);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessEdge %d = %d (EG_tessHOverts)!\n",
               i+1, stat);
      goto cleanup;
    }
  }
  EG_free(coords);
  coords = NULL;
#ifdef DEBUG
  printf(" quads = %d %d,  nIns = %d,  nmid = %d\n", quad, qout, nIns, nmid);
  for (i = 0; i < nst; i++)
    if (type[i] >= 0) {
      printf("   %d: type = %2d\n", i, type[i]);
    } else {
      printf("   %d: type = %2d   frac = %lf\n", i, type[i], frac[i]);
    }
#endif
  
  /* size and allocate temporary Face arrays */
  
  tessel = (egTessel *) newTess->blind;
  stat   = EG_getBodyTopos(body, NULL, FACE, &nfaces, &faces);
  if ((stat != EGADS_SUCCESS) || (faces == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getBodyTopos F = %d (EG_tessHOverts)!\n", stat);
    if (stat == EGADS_SUCCESS) stat = EGADS_TOPOERR;
    goto cleanup;
  }
  for (nside = npts = ntris = i = 0; i < nfaces; i++) {
    stat = EG_getTessFace(tess, i+1, &np, &xyzs, &uvs, &ptype, &pindex,
                          &nt, &trs, &trc);
    if ((stat != EGADS_SUCCESS) || (nt == 0)) {
      if (outLevel > 0)
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: EG_getTessFace %d = %d (EG_tessHOverts)!\n",
                 i+1, stat);
        } else {
          printf(" EGADS Error: Face %d has no tessellation (EG_tessHOverts)!\n",
                 i+1);
        }
      goto cleanup;
    }
    for (sum[0] = sum[1] = j = 0; j < nt; j++)
      for (k = 0; k < 3; k++)
        if (trc[3*j+k] > 0) {
          sum[0]++;
        } else {
          sum[1]++;
        }
    k = sum[0]/2 + sum[1];
    if (quad == 1) k -= nt/2;
    if (k  > nside) nside = k;
    if (np > npts)  npts  = np;
    if (nt > ntris) ntris = nt;
  }
#ifdef DEBUG
  printf(" Max npts = %d, nside = %d, ntris = %d\n", npts, nside, ntris);
#endif
  
  if (quad == 0) {
    k = npts + 3*nIns*nside + ntris*nmid;
  } else {
    k = npts + 4*nIns*nside + ntris*nmid/2;
  }
  coords = (double *) EG_alloc(5*k*sizeof(double));
  if (coords == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face points (EG_tessHOverts)!\n", k);
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  parms = &coords[3*k];
  tris  = (int *) EG_alloc(nItri*3*ntris*sizeof(int));
  if (tris == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face tris (EG_tessHOverts)!\n",
             3*ntris);
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  k = nst*ntris;
  if (quad == 1) k /= 2;
  elems  = (int *) EG_alloc(k*sizeof(int));
  if (elems == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc for %d Face elems (EG_tessHOverts)!\n", k);
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  
  /* fill in the Faces */
  
  for (i = 0; i < nfaces; i++) {
    stat = EG_getTessFace(tess, i+1, &np, &xyzs, &uvs, &ptype, &pindex,
                          &nt, &trs, &trc);
    if (stat != EGADS_SUCCESS) continue;
    
    /* find degenerate nodes (if any) */
    degens[0] = degens[1] = 0;
    iuv[0]    = iuv[1]    = 0;
    stat      = EG_getBodyTopos(body, faces[i], EDGE, &k, &objs);
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
        n = EG_indexBodyTopo(body, nodes[0]);
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
#ifdef DEBUG
    if (degens[0] != 0)
      printf(" EGADS Info: Face %d has degenerate Node(s) = %d (%d)  %d (%d)\n",
             i+1, degens[0], iuv[0], degens[1], iuv[1]);
#endif
    
    /* clear element vert positions */
    for (i0 = j = 0; j < nt/(quad+1); j++)
      for (k = 0; k < nst; k++, i0++) elems[i0] = 0;
    
    /* copy source verts */
    for (j = 0; j < np; j++) {
      coords[3*j  ] = xyzs[3*j  ];
      coords[3*j+1] = xyzs[3*j+1];
      coords[3*j+2] = xyzs[3*j+2];
      parms[2*j  ]  = uvs[2*j  ];
      parms[2*j+1]  = uvs[2*j+1];
    }
    
    /* fill in corner verts */
    if (quad == 0) {
      for (i0 = j = 0; j < nt; j++)
        for (k = 0; k < nst; k++, i0++)
          if (type[k] > 0) elems[i0] = trs[3*j+type[k]-1];
    } else {
      for (i0 = j = 0; j < nt; j+=2) {
        corner[0] = trs[3*j  ];
        corner[1] = trs[3*j+1];
        corner[2] = trs[3*j+2];
        corner[3] = trs[3*j+5];
        for (k = 0; k < nst; k++, i0++)
          if (type[k] > 0) elems[i0] = corner[type[k]-1];
      }
    }
    
    /* fill in the side verts */
    if (quad == 0) {
      for (j = 0; j < nt; j++) {
        for (k = 0; k < 3; k++) {
          nei       = abs(trc[3*j+k]) - 1;
          corner[0] = trs[3*j+sidet[k][0]];
          corner[1] = trs[3*j+sidet[k][1]];
          uvl = uvr = NULL;
          if (trc[3*j+k] < 0) {
            stat = EG_getTopology(edges[nei], &geom, &oclass, &mtype,
                                  trange, &n, &objs, &senses);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: EG_getTopo %d = %d (EG_tessHOverts)!\n",
                       nei+1, stat);
              goto cleanup;
            }
            corner[2] = corner[3] = EG_indexBodyTopo(body, objs[0]);
            if (mtype == TWONODE) corner[3] = EG_indexBodyTopo(body, objs[1]);
            stat = EG_fillEdgeSeg(faces[i], edges, tessel, nIns, corner,
                                  nei, xyzs, uvs, ptype, pindex, k, nst,
                                  type, frac, sinsert, degens, j*nst, elems,
                                  &np, coords, parms);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: EG_fillEdgeSeg %d %d/%d = %d (EG_tessHOverts)!\n",
                       i+1, j+1, k+1, stat);
              goto cleanup;
            }
          } else {
            if (nei < j) continue;
            uvl = &uvs[2*trs[3*j+k]-2];
            i0  = trs[3*nei]+trs[3*nei+1]+trs[3*nei+2]-corner[0]-corner[1];
            uvr = &uvs[2*i0-2];
            i1  = -1;
            if (trs[3*nei  ] == i0) i1 = 0;
            if (trs[3*nei+1] == i0) i1 = 1;
            if (trs[3*nei+2] == i0) i1 = 2;
            if (i1 == -1) {
              printf(" FATAL: *** Can't find other side! ***\n");
              stat = EGADS_INDEXERR;
              goto cleanup;
            }
#ifdef REPOSITION
            EG_correctEndPts(corner[0]-1, corner[1]-1, degens, iuv, uvs, ptype,
                             pindex, uvm, uvp);
#endif
            for (i0 = 0; i0 < nIns; i0++) {
#ifdef REPOSITION
              EG_getSidepoint(faces[i], sinsert[i0], uvm, uvp, uvl, uvr, uv);
#else
              uv[0] = (1.0-sinsert[i0])*uvs[2*corner[0]-2] +
                           sinsert[i0] *uvs[2*corner[1]-2];
              uv[1] = (1.0-sinsert[i0])*uvs[2*corner[0]-1] +
                           sinsert[i0] *uvs[2*corner[1]-1];
              EG_correctUV(uv, corner[0]-1, corner[1]-1, degens, iuv,
                           uvs, ptype, pindex);
#endif
              stat = EG_evaluate(faces[i], uv, result);
              if (stat != EGADS_SUCCESS) {
                if (outLevel > 0)
                  printf(" EGADS Error: evaluate %d %d/%d = %d (EG_tessHOverts)!\n",
                         i+1, j+1, k+1, stat);
                goto cleanup;
              }
              coords[3*np  ] = result[0];
              coords[3*np+1] = result[1];
              coords[3*np+2] = result[2];
              parms[2*np  ]  = uv[0];
              parms[2*np+1]  = uv[1];
              np++;
              i2 = EG_findSideIndex(k, sinsert[i0], nst, type, frac);
              if (i2 < 0) {
                if (outLevel > 0)
                  printf(" EGADS Error: Cannot find %d %lf (EG_tessHOverts)!\n",
                         k, sinsert[i0]);
                stat = EGADS_INDEXERR;
                goto cleanup;
              }
              elems[j*nst+i2] = np;
              i2 = EG_findSideIndex(i1, 1.0-sinsert[i0], nst, type, frac);
              if (i2 < 0) {
                if (outLevel > 0)
                  printf(" EGADS Error: Cannot find %d %lf (EG_tessHOverts)!\n",
                         k, sinsert[i0]);
                stat = EGADS_INDEXERR;
                goto cleanup;
              }
#ifdef DEBUG
              if (elems[nei*nst+i2] != 0) printf(" double hit!\n");
#endif
              elems[nei*nst+i2] = np;
            }
          }
        }
      }
    } else {
      for (j = 0; j < nt; j+=2) {
        for (k = 0; k < 4; k++) {
          nei       = abs(trc[3*j+neigq[k]]) - 1;
          corner[0] = trs[3*j+sideq[k][0]];
          corner[1] = trs[3*j+sideq[k][1]];
          if (trc[3*j+neigq[k]] < 0) {
            stat = EG_getTopology(edges[nei], &geom, &oclass, &mtype,
                                  trange, &n, &objs, &senses);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: EG_getTopo %d = %d (EG_tessHOverts)!\n",
                       nei+1, stat);
              goto cleanup;
            }
            corner[2] = corner[3] = EG_indexBodyTopo(body, objs[0]);
            if (mtype == TWONODE) corner[3] = EG_indexBodyTopo(body, objs[1]);
            stat = EG_fillEdgeSeg(faces[i], edges, tessel, nIns, corner,
                                  nei, xyzs, uvs, ptype, pindex, k, nst,
                                  type, frac, sinsert, degens, j*nst/2, elems,
                                  &np, coords, parms);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: EG_fillEdgeSeg %d %d/%d = %d (EG_tessHOverts)!\n",
                       i+1, j+1, k+1, stat);
              goto cleanup;
            }
          } else {
            if (nei < j) continue;
            if (nei%2 == 1) nei--;
            i1 = -1;
            for (i0 = 0; i0 < 4; i0++) {
              corner[2] = trs[3*nei+sideq[i0][0]];
              corner[3] = trs[3*nei+sideq[i0][1]];
              if ((corner[0] == corner[2]) && (corner[1] == corner[3])) {
                i1 = i0;
                break;
              }
              if ((corner[0] == corner[3]) && (corner[1] == corner[2])) {
                i1 = i0;
                break;
              }
            }
            if (i1 == -1) {
              printf(" FATAL: *** Can't find other Q side! ***\n");
              stat = EGADS_INDEXERR;
              goto cleanup;
            }
#ifdef REPOSITION
            EG_correctEndPts(corner[0]-1, corner[1]-1, degens, iuv, uvs, ptype,
                             pindex, uvm, uvp);
#endif
            for (i0 = 0; i0 < nIns; i0++) {
#ifdef REPOSITION
              EG_getSidepoint(faces[i], sinsert[i0], uvm, uvp, NULL, NULL, uv);
#else
              uv[0] = (1.0-sinsert[i0])*uvs[2*corner[0]-2] +
                           sinsert[i0] *uvs[2*corner[1]-2];
              uv[1] = (1.0-sinsert[i0])*uvs[2*corner[0]-1] +
                           sinsert[i0] *uvs[2*corner[1]-1];
              EG_correctUV(uv, corner[0]-1, corner[1]-1, degens, iuv,
                           uvs, ptype, pindex);
#endif
              stat = EG_evaluate(faces[i], uv, result);
              if (stat != EGADS_SUCCESS) {
                if (outLevel > 0)
                  printf(" EGADS Error: Evaluate %d %d/%d = %d (EG_tessHOverts)!\n",
                         i+1, j+1, k+1, stat);
                goto cleanup;
              }
              coords[3*np  ] = result[0];
              coords[3*np+1] = result[1];
              coords[3*np+2] = result[2];
              parms[2*np  ]  = uv[0];
              parms[2*np+1]  = uv[1];
              np++;
              i2 = EG_findSideIndex(k, sinsert[i0], nst, type, frac);
              if (i2 < 0) {
                if (outLevel > 0)
                  printf(" EGADS Error: Cannot find %d %lf (EG_tessHOverts)!\n",
                         k, sinsert[i0]);
                stat = EGADS_INDEXERR;
                goto cleanup;
              }
              elems[j*nst/2+i2] = np;
              i2 = EG_findSideIndex(i1, 1.0-sinsert[i0], nst, type, frac);
              if (i2 < 0) {
                if (outLevel > 0)
                  printf(" EGADS Error: Cannot find %d %lf (EG_tessHOverts)!\n",
                         k, sinsert[i0]);
                stat = EGADS_INDEXERR;
                goto cleanup;
              }
#ifdef DEBUG
              if (elems[nei*nst/2+i2] != 0) printf(" double hit!\n");
#endif
              elems[nei*nst/2+i2] = np;
            }
          }
        }
      }
    }
    
    /* fill in the interior verts */
    if (nmid != 0)
      if (quad == 0) {
        double uv0[2], uv1[2], uv2[2];
        for (k = 0; k < nst; k++) {
          if (type[k] != 0) continue;
          for (j = 0; j < nt; j++) {
            w[1]    = st[2*k  ];
            w[2]    = st[2*k+1];
            w[0]    = 1.0 - w[1] - w[2];
            i0      = trs[3*j  ] - 1;
            i1      = trs[3*j+1] - 1;
            i2      = trs[3*j+2] - 1;
            uv0[0]  = uvs[2*i0  ];
            uv0[1]  = uvs[2*i0+1];
            uv1[0]  = uvs[2*i1  ];
            uv1[1]  = uvs[2*i1+1];
            uv2[0]  = uvs[2*i2  ];
            uv2[1]  = uvs[2*i2+1];
            EG_interiorTri(body, faces[i], edges, btess, &trs[3*j], &trc[3*j],
                           degens, iuv, uvs, ptype, pindex, w, uv0, uv1, uv2);
            uv[0]  = w[0]*uv0[0] + w[1]*uv1[0] + w[2]*uv2[0];
            uv[1]  = w[0]*uv0[1] + w[1]*uv1[1] + w[2]*uv2[1];
#ifdef REPOSITION
            stat = EG_baryInsert(faces[i], w[0], w[1], w[2], uv0, uv1, uv2, uv);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Info: EG_baryInsert = %d (EG_tessHOverts)!\n",
                     stat);
              uv[0]  = w[0]*uv0[0] + w[1]*uv1[0] + w[2]*uv2[0];
              uv[1]  = w[0]*uv0[1] + w[1]*uv1[1] + w[2]*uv2[1];
              xyz[0] = w[0]*xyzs[3*i0  ] + w[1]*xyzs[3*i1  ] + w[2]*xyzs[3*i2  ];
              xyz[1] = w[0]*xyzs[3*i0+1] + w[1]*xyzs[3*i1+1] + w[2]*xyzs[3*i2+1];
              xyz[2] = w[0]*xyzs[3*i0+2] + w[1]*xyzs[3*i1+2] + w[2]*xyzs[3*i2+2];
              EG_getInterior(faces[i], xyz, uv);
            }
#endif
            stat = EG_evaluate(faces[i], uv, result);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: eval %d %d/%d = %d (EG_tessHOverts)!\n",
                       i+1, j+1, k+1, stat);
              goto cleanup;
            }
            coords[3*np  ] = result[0];
            coords[3*np+1] = result[1];
            coords[3*np+2] = result[2];
            parms[2*np  ]  = uv[0];
            parms[2*np+1]  = uv[1];
            np++;
#ifdef DEBUG
            if (elems[j*nst+k] != 0) printf(" double hit!\n");
#endif
            elems[j*nst+k] = np;
          }
        }
      } else {
        double el[3], eu[3], xl[3], xu[3];
        for (k = 0; k < nst; k++) {
          if (type[k] != 0) continue;
          for (j = 0; j < nt; j+=2) {
            i0    = trs[3*j  ] - 1;
            i1    = trs[3*j+1] - 1;
            i2    = trs[3*j+2] - 1;
            i3    = trs[3*j+5] - 1;
            u0[0] = uvs[2*i0  ];
            u0[1] = uvs[2*i0+1];
            u1[0] = uvs[2*i1  ];
            u1[1] = uvs[2*i1+1];
            u2[0] = uvs[2*i2  ];
            u2[1] = uvs[2*i2+1];
            u3[0] = uvs[2*i3  ];
            u3[1] = uvs[2*i3+1];
#ifdef REPOSITION
            EG_correctEndPts(i0, i1, degens, iuv, uvs, ptype, pindex, uvm, uvp);
            EG_getSidepoint(faces[i], st[2*k  ], uvm, uvp, NULL, NULL, xl);
            EG_correctEndPts(i3, i2, degens, iuv, uvs, ptype, pindex, uvm, uvp);
            EG_getSidepoint(faces[i], st[2*k  ], uvm, uvp, NULL, NULL, xu);
            EG_correctEndPts(i0, i3, degens, iuv, uvs, ptype, pindex, uvm, uvp);
            EG_getSidepoint(faces[i], st[2*k+1], uvm, uvp, NULL, NULL, el);
            EG_correctEndPts(i1, i2, degens, iuv, uvs, ptype, pindex, uvm, uvp);
            EG_getSidepoint(faces[i], st[2*k+1], uvm, uvp, NULL, NULL, eu);
#else
            xl[0] = (1.0-st[2*k  ])*u0[0]  + st[2*k  ]*u1[0];
            xl[1] = (1.0-st[2*k  ])*u0[1]  + st[2*k  ]*u1[1];
            EG_correctUV(xl, i0, i1, degens, iuv, uvs, ptype, pindex);
            xu[0] = (1.0-st[2*k  ])*u3[0]  + st[2*k  ]*u2[0];
            xu[1] = (1.0-st[2*k  ])*u3[1]  + st[2*k  ]*u2[1];
            EG_correctUV(xu, i3, i2, degens, iuv, uvs, ptype, pindex);
            el[0] = (1.0-st[2*k+1])*u0[0]  + st[2*k+1]*u3[0];
            el[1] = (1.0-st[2*k+1])*u0[1]  + st[2*k+1]*u3[1];
            EG_correctUV(el, i0, i3, degens, iuv, uvs, ptype, pindex);
            eu[0] = (1.0-st[2*k+1])*u1[0]  + st[2*k+1]*u2[0];
            eu[1] = (1.0-st[2*k+1])*u1[1]  + st[2*k+1]*u2[1];
            EG_correctUV(eu, i1, i2, degens, iuv, uvs, ptype, pindex);
#endif
            EG_correctUVq(u0, u1, u2, u3, i0, i1, i2, i3, degens, iuv,
                          ptype, pindex);
            /* side(s) on Edge(s)? */
            if (trc[3*j  ] < 0)
              EG_evalEdgeSeg(body, faces[i], edges, btess, -trc[3*j  ]-1,
                             i1, i2, st[2*k+1], uvs, ptype, pindex, degens, eu);
            if (trc[3*j+3] < 0)
              EG_evalEdgeSeg(body, faces[i], edges, btess, -trc[3*j+3]-1,
                             i3, i2, st[2*k  ], uvs, ptype, pindex, degens, xu);
            if (trc[3*j+4] < 0)
              EG_evalEdgeSeg(body, faces[i], edges, btess, -trc[3*j+4]-1,
                             i0, i3, st[2*k+1], uvs, ptype, pindex, degens, el);
            if (trc[3*j+2] < 0)
              EG_evalEdgeSeg(body, faces[i], edges, btess, -trc[3*j+2]-1,
                             i0, i1, st[2*k  ], uvs, ptype, pindex, degens, xl);
            EG_mdTFI(st[2*k], st[2*k+1], 2, xl, xu, el, eu, u0, u1, u2, u3, uv);
#ifdef REPOSITION
/*
            double xmid[4][18];
            stat = EG_evaluate(faces[i], xl, xmid[0]);
            stat = EG_evaluate(faces[i], xu, xmid[1]);
            stat = EG_evaluate(faces[i], el, xmid[2]);
            stat = EG_evaluate(faces[i], eu, xmid[3]);
            EG_mdTFI(st[2*k], st[2*k+1], 3, xmid[0], xmid[1], xmid[2], xmid[3],
                     &xyzs[3*i0], &xyzs[3*i1], &xyzs[3*i2], &xyzs[3*i3], xyz);
            EG_getInterior(faces[i], xyz, uv);  */
            EG_minArc4(faces[i], st[2*k], st[2*k+1], xl, eu, xu, el, uv);
#endif
            stat = EG_evaluate(faces[i], uv, result);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: eval %d %d/%d = %d (EG_tessHOverts)!\n",
                       i+1, j+1, k+1, stat);
              goto cleanup;
            }
            coords[3*np  ] = result[0];
            coords[3*np+1] = result[1];
            coords[3*np+2] = result[2];
            parms[2*np  ]  = uv[0];
            parms[2*np+1]  = uv[1];
            np++;
#ifdef DEBUG
            if (elems[j*nst/2+k] != 0) printf(" double hit!\n");
#endif
            elems[j*nst/2+k] = np;
          }
        }
      }
    
    /* fill up the triangles */
    i1 = nt;
    if (quad == 1) i1 /= 2;
    for (ntris = i0 = j = 0; j < i1; j++, i0+=nst)
      for (k = 0; k < nItri; k++, ntris++) {
        tris[3*ntris  ] = elems[i0+iTris[3*k  ]-1];
        tris[3*ntris+1] = elems[i0+iTris[3*k+1]-1];
        tris[3*ntris+2] = elems[i0+iTris[3*k+2]-1];
      }
    
#ifdef DEBUG
    printf(" Face %d: npts = %d, ntris = %d\n", i+1, np, ntris);
#endif
    stat = EG_setTessFace(newTess, i+1, np, coords, parms, ntris, tris);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_setTessFace %d = %d (EG_tessHOverts)!\n",
               i+1, stat);
      goto cleanup;
    }
  }
  
  /* close up the open tessellation */
  stat = EG_statusTessBody(newTess, &geom, &i, &npts);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_statusTessBody = %d (EG_tessHOverts)!\n", stat);
    goto cleanup;
  }
  if (i != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: New Tessellation Object is Open (EG_tessHOverts)!\n");
    stat = EGADS_TESSTATE;
    goto cleanup;
  }

  /* mark tessellation as quadded if appropriate */
  if ((quad == 1) && (qout == 1)) {
    senses = (int *) EG_alloc(nfaces*sizeof(int));
    if (senses != NULL) {
      for (i = 0; i < nfaces; i++) senses[i] = tessel->tess2d[i].ntris/2;
      stat = EG_attributeAdd(newTess, ".mixed", ATTRINT, nfaces,
                             senses, NULL, NULL);
      if (stat != EGADS_SUCCESS)
        if (outLevel > 0)
          printf(" EGADS Warning: EG_attributeAdd m = %d (EG_tessHOverts)!\n",
                 stat);
      EG_free(senses);
    }
    stat = EG_attributeAdd(newTess, ".tessType", ATTRSTRING, 4,
                           NULL, NULL, "Quad");
    if (stat != EGADS_SUCCESS)
      if (outLevel > 0)
        printf(" EGADS Warning: EG_attributeAdd Q = %d (EG_tessHOverts)!\n",
               stat);
  }
  
  *nTess  = newTess;
  newTess = NULL;
  stat    = EGADS_SUCCESS;
  
cleanup:
  EG_free(type);
  if (frac    != NULL) EG_free(frac);
  if (elems   != NULL) EG_free(elems);
  if (tris    != NULL) EG_free(tris);
  if (faces   != NULL) EG_free(faces);
  if (coords  != NULL) EG_free(coords);
  if (edges   != NULL) EG_free(edges);
  if (newTess != NULL) EG_deleteObject(newTess);
  return stat;
}
