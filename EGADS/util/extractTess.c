/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Extract tesselations from matched input body tessellations
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "egads.h"


#define CROSS(a,b,c)       a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                           a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                           a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

  typedef struct {
    ego body;
    int num;
    int *objs;
  } objMatch;

  typedef struct {
    int    hit;
    double xyz[3];
  } nodeMatch;



static int
EG_lookupEdge(int index, objMatch ematch)
{
  int i;
  
  for (i = 0; i < ematch.num; i++)
    if (ematch.objs[2*i+1] == index) return ematch.objs[2*i];
  
  return 0;
}


static int
EG_findEdge(ego body, const double *xyz, int is, const int *tris,
            const int *tric, const int *ptype, const int *pindex,
            objMatch ematch, egTessel *btess, double *t)
{
  int    i, m, n, ie, stat, oclass, mtype, nnds, *senses;
  double trange[2], coord[3], dmin, dmax;
  ego    edge, ref, *nds, *objs;
  
  for (i = 0; i < 3; i++) {
    if (i == is) continue;
    if (tric[i] < 0) break;
  }
  if (i == 3) return EGADS_NOTFOUND;
  
  ie   = EG_lookupEdge(-tric[i], ematch);
  if (ie == 0) return ie;
  stat = EG_objectBodyTopo(body, EDGE, ie, &edge);
  if (stat != EGADS_SUCCESS) {
        printf(" Error: EG_objectBodyTopo = %d for Edge %d\n", stat, ie);
    return stat;
  }
  
  stat = EG_getTopology(edge, &ref, &oclass, &mtype, trange, &nnds, &nds,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" Error: EG_getTopology = %d for Edge %d\n", stat, ie);
    return stat;
  }
  
  if (nnds == 2) {
    
    /* find the closest Node */
    stat = EG_getTopology(nds[0], &ref, &oclass, &mtype, coord, &n, &objs,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" Error: EG_getTopology = %d for start Node\n", stat);
      return stat;
    }
    dmin = sqrt((coord[0]-xyz[0])*(coord[0]-xyz[0]) +
                (coord[1]-xyz[1])*(coord[1]-xyz[1]) +
                (coord[2]-xyz[2])*(coord[2]-xyz[2]));
    
    stat = EG_getTopology(nds[1], &ref, &oclass, &mtype, coord, &n, &objs,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" Error: EG_getTopology = %d for finis Node\n", stat);
      return stat;
    }
    dmax = sqrt((coord[0]-xyz[0])*(coord[0]-xyz[0]) +
                (coord[1]-xyz[1])*(coord[1]-xyz[1]) +
                (coord[2]-xyz[2])*(coord[2]-xyz[2]));
    *t = trange[0];
    if (dmax < dmin) *t = trange[1];
    
  } else {
    
    /* look at the other vert's parameter */
    m = tris[0] + tris[1] + tris[2] - tris[is] - tris[i];
    if (pindex[m-1] != ie) {
      printf(" Error: Edge mismatch - %d %d\n", pindex[m-1], ie);
      return EGADS_TOPOERR;
    }
    if (btess->tess1d[ie-1].npts <= 3)
      printf(" Info: Ambiguous Edge tessellation - %d\n", ie);
    if (ptype[m-1] == 2) {
      *t = trange[0];
    } else {
      *t = trange[1];
    }
    
  }
  
  return ie;
}


/* generates a Body tessellation copying matching Edge/Face tessellation from
 *           input tessellations
 *
 * where: niTess - the number of input Body tessellations
 *        iTess  - an array of input Body tessellations (niTess in length)
 *        body   - the Body Object to tessellate
 *        params - the 3 tessellation parameters (see EG_makeBodyTess)
 *        tess   - the returned Tessellation object
 */

int
EG_extractTess(int niTess, ego *iTess, ego body, double *params, ego *tess)
{
  int          i, j, k, l, m, n, nobj, stat, npts, oclass, mtype, tesstate, len;
  int          *mark, *senses, *trix, *cnt, ntri, ie, is, aType, aLen, *qints;
  double       t, trange[2], xyz[3], coord[3], uv[2], result[18], an[3], ds[3];
  double       x0[3], x1[3], *u, *v, *nprms, *uvs;
  const int    *ptype, *pindex, *tris, *tric, *aInts;
  const double *xyzs, *prms, *aReals;
  const char   *aStr;
  ego          top, prev, next, *objs, *nds;
  egTessel     *btess, *otess;
  objMatch     *matches, *fmatch;
  nodeMatch    *nodes;
  
  /* check inputs */
  if (niTess <= 0) return EGADS_INDEXERR;
  for (i = 0; i < niTess; i++) {
    stat = EG_getInfo(iTess[i], &oclass, &mtype, &top, &prev, &next);
    if (stat   != EGADS_SUCCESS) return stat;
    if (oclass != TESSELLATION)  return EGADS_NOTTESS;
  }
  
  /* allocate space to put the matches */
  stat = EG_initTessBody(body, tess);
  if (stat != EGADS_SUCCESS) return stat;
  
  matches = (objMatch *) EG_alloc(2*niTess*sizeof(objMatch));
  if (matches == NULL) {
    EG_deleteObject(*tess);
    return EGADS_MALLOC;
  }
  for (i = 0; i < 2*niTess; i++) {
    matches[i].num  = 0;
    matches[i].objs = NULL;
  }
  fmatch = &matches[niTess];
  nprms  = NULL;
  mark   = NULL;
  nodes  = NULL;
  uvs    = NULL;
  objs   = NULL;
  
  /* get the Nodes we may expect */
  stat   = EG_getBodyTopos(body, NULL, NODE, &nobj, NULL);
  if (stat != EGADS_SUCCESS) goto cleanup;
  nodes  = (nodeMatch *) EG_alloc(nobj*sizeof(nodeMatch));
  if (nodes == NULL) {
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  for (i = 0; i < nobj; i++) {
    nodes[i].hit    = 0;
    nodes[i].xyz[0] = nodes[i].xyz[1] = nodes[i].xyz[2] = 0.0;
  }
  
  /* get the Edge matches */
  for (n = i = 0; i < niTess; i++) {
    stat = EG_statusTessBody(iTess[i], &matches[i].body, &tesstate, &npts);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_matchBodyEdges(body, matches[i].body, 0.0, &matches[i].num,
                             &matches[i].objs);
    if (stat != EGADS_SUCCESS) goto cleanup;
    n += matches[i].num;
  }
#ifdef REPORT
  printf(" *** Edge Matches = %d ***\n", n);
#endif
  stat = EGADS_INDEXERR;
  if (n == 0) goto cleanup;
  
  /* fill in the matched Edges */
  
  stat = EG_getBodyTopos(body, NULL, EDGE, &nobj, &objs);
  if ((stat != EGADS_SUCCESS) || (objs == NULL)) goto cleanup;
  mark = (int *) EG_alloc(nobj*sizeof(int));
  stat = EGADS_MALLOC;
  if (mark == NULL) goto cleanup;
  for (j = 0; j < nobj; j++) mark[j] = 0;
  for (i = 0; i < niTess; i++)
    for (j = 0; j < matches[i].num; j++) {
      m = matches[i].objs[2*j];
      if (mark[m-1] != 0) continue;
      stat = EG_getTopology(objs[m-1], &top, &oclass, &mtype, trange, &npts,
                            &nds, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      
      /* get the Edge tessellation */
      stat = EG_getTessEdge(iTess[i], matches[i].objs[2*j+1], &len, &xyzs,
                            &prms);
      if (stat != EGADS_SUCCESS) goto cleanup;
      mark[m-1] = 1;
      for (k = 0; k < npts; k++) {
        stat = EG_indexBodyTopo(body, nds[k]);
        if (stat < EGADS_SUCCESS) goto cleanup;
        if (nodes[stat-1].hit != 0) continue;
        if (k == 0) {
          nodes[stat-1].xyz[0] = xyzs[0];
          nodes[stat-1].xyz[1] = xyzs[1];
          nodes[stat-1].xyz[2] = xyzs[2];
        } else {
          nodes[stat-1].xyz[0] = xyzs[3*len-3];
          nodes[stat-1].xyz[1] = xyzs[3*len-2];
          nodes[stat-1].xyz[2] = xyzs[3*len-1];
        }
        nodes[stat-1].hit = 1;
      }
      
      /* are we the same? */
      stat = EG_objectBodyTopo(matches[i].body, EDGE, matches[i].objs[2*j+1],
                               &top);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (EG_isEquivalent(objs[m-1], top) == EGADS_SUCCESS) {
        /* yes! */
#ifdef REPORT
        printf(" *** Tess Edge %d = Body Edge %d ***\n",
               matches[i].objs[2*j+1], m);
#endif
        stat = EG_setTessEdge(*tess, m, len, xyzs, prms);
        if (stat != EGADS_SUCCESS) goto cleanup;
        
      } else {
        /* no -- recompute the parameters */
#ifdef REPORT
        printf(" *** Tess Edge %d ~ Body Edge %d ***\n",
               matches[i].objs[2*j+1], m);
#endif
        nprms = (double *) EG_alloc(len*sizeof(double));
        if (nprms == NULL) {
          stat = EGADS_MALLOC;
          goto cleanup;
        }
        nprms[0]     = trange[0];
        nprms[len-1] = trange[1];
        for (k = 1; k < len-1; k++) {
          coord[0] = xyzs[3*k  ];
          coord[1] = xyzs[3*k+1];
          coord[2] = xyzs[3*k+2];
          stat     = EG_invEvaluate(objs[m-1], coord, &nprms[k], xyz);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }
        stat = EG_setTessEdge(*tess, m, len, xyzs, nprms);
        if (stat != EGADS_SUCCESS) goto cleanup;
        EG_free(nprms);
        nprms = NULL;
      }
    }

  EG_free(mark);
  mark  = NULL;
  EG_free(objs);
  objs  = NULL;
  btess = (egTessel *) (*tess)->blind;
 
  /* get the Face matches */
  
  for (n = i = 0; i < niTess; i++) {
    stat = EG_statusTessBody(iTess[i], &fmatch[i].body, &tesstate, &npts);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_matchBodyFaces(body, fmatch[i].body, 0.0, &fmatch[i].num,
                             &fmatch[i].objs);
    if (stat != EGADS_SUCCESS) goto cleanup;
    n += fmatch[i].num;
  }
#ifdef REPORT
  printf(" *** Face Matches = %d ***\n", n);
#endif

  if (n != 0) {

    /* fill in the matched Faces */
    stat = EG_getBodyTopos(body, NULL, FACE, &nobj, &objs);
    if ((stat != EGADS_SUCCESS) || (objs == NULL)) goto cleanup;
    mark = (int *) EG_alloc(2*nobj*sizeof(int));
    stat = EGADS_MALLOC;
    if (mark == NULL) goto cleanup;
    qints = &mark[nobj];
    for (j = 0; j < nobj; j++) mark[j] = qints[j] = 0;
    for (i = 0; i < niTess; i++) {
      otess = (egTessel *) iTess[i]->blind;
      for (j = 0; j < fmatch[i].num; j++) {
        m = fmatch[i].objs[2*j];
        if (mark[m-1] != 0) continue;
        /* get the Face tessellation */
        stat = EG_getTessFace(iTess[i], fmatch[i].objs[2*j+1], &len, &xyzs,
                              &prms, &ptype, &pindex, &ntri, &tris, &tric);
        if (stat != EGADS_SUCCESS) goto cleanup;
        if (ntri == 0) continue;
        mark[m-1] = 1;
        uvs       = (double *) EG_alloc(2*ntri*sizeof(double));
        stat      = EGADS_MALLOC;
        if (uvs == NULL) goto cleanup;
        /* is the orientation correct? */
        for (npts = k = 0; k < ntri; k++) {
          coord[0] = (xyzs[3*(tris[3*k  ]-1)  ] + xyzs[3*(tris[3*k+1]-1)  ] +
                      xyzs[3*(tris[3*k+2]-1)  ])/3.0;
          coord[1] = (xyzs[3*(tris[3*k  ]-1)+1] + xyzs[3*(tris[3*k+1]-1)+1] +
                      xyzs[3*(tris[3*k+2]-1)+1])/3.0;
          coord[2] = (xyzs[3*(tris[3*k  ]-1)+2] + xyzs[3*(tris[3*k+1]-1)+2] +
                      xyzs[3*(tris[3*k+2]-1)+2])/3.0;
          stat = EG_invEvaluate(objs[m-1], coord, uv, xyz);
          if (stat != EGADS_SUCCESS) goto cleanup;
          uvs[2*k  ] = uv[0];
          uvs[2*k+1] = uv[1];
          stat = EG_evaluate(objs[m-1], uv, result);
          if (stat != EGADS_SUCCESS) goto cleanup;
          u = &result[3];
          v = &result[6];
          CROSS(an, u, v);
          t = sqrt(DOT(an, an))*objs[m-1]->mtype;
          if (t != 0.0) {
            an[0] /= t;
            an[1] /= t;
            an[2] /= t;
          }
          x0[0] = xyzs[3*(tris[3*k  ]-1)  ] - xyzs[3*(tris[3*k+2]-1)  ];
          x1[0] = xyzs[3*(tris[3*k+1]-1)  ] - xyzs[3*(tris[3*k+2]-1)  ];
          x0[1] = xyzs[3*(tris[3*k  ]-1)+1] - xyzs[3*(tris[3*k+2]-1)+1];
          x1[1] = xyzs[3*(tris[3*k+1]-1)+1] - xyzs[3*(tris[3*k+2]-1)+1];
          x0[2] = xyzs[3*(tris[3*k  ]-1)+2] - xyzs[3*(tris[3*k+2]-1)+2];
          x1[2] = xyzs[3*(tris[3*k+1]-1)+2] - xyzs[3*(tris[3*k+2]-1)+2];
          CROSS(ds, x0, x1);
          t = sqrt(DOT(ds, ds));
          if (t != 0.0) {
            ds[0] /= t;
            ds[1] /= t;
            ds[2] /= t;
          }
          if (DOT(ds, an) > 0.0) npts++;
        }
        trix = (int *) tris;
        if (npts < ntri/2) {
#ifdef REPORT
          printf("   reorienting tessellation\n");
#endif
          trix = (int *) EG_alloc(3*ntri*sizeof(int));
          stat = EGADS_MALLOC;
          if (trix == NULL) goto cleanup;
          if (otess->tess2d[fmatch[i].objs[2*j+1]-1].tfi == 1) {
            for (k = 0; k < ntri/2; k++) {
              trix[6*k  ] = tris[6*k  ];
              trix[6*k+1] = tris[6*k+5];
              trix[6*k+2] = tris[6*k+2];
              trix[6*k+3] = tris[6*k+3];
              trix[6*k+4] = tris[6*k+4];
              trix[6*k+5] = tris[6*k+1];
            }
          } else {
            for (k = 0; k < ntri; k++) {
              trix[3*k  ] = tris[3*k  ];
              trix[3*k+1] = tris[3*k+2];
              trix[3*k+2] = tris[3*k+1];
            }
          }
        }
        
        /* are we the same? */
        stat = EG_objectBodyTopo(fmatch[i].body, FACE, fmatch[i].objs[2*j+1],
                                 &top);
        if (stat != EGADS_SUCCESS) {
          if (npts < ntri/2) EG_free(trix);
          goto cleanup;
        }
        if (EG_isEquivalent(objs[m-1], top) == EGADS_SUCCESS) {
          /* yes! */
#ifdef REPORT
          printf(" *** Tess Face %d = Body Face %d ***\n",
                 fmatch[i].objs[2*j+1], m);
#endif
          stat = EG_setTessFace(*tess, m, len, xyzs, prms, ntri, trix);
          if (stat != EGADS_SUCCESS) {
            if (npts < ntri/2) EG_free(trix);
            goto cleanup;
          }
        } else {
          /* no -- recompute the parameters */
#ifdef REPORT
          printf(" *** Tess Face %d ~ Body Face %d ***\n",
                 fmatch[i].objs[2*j+1], m);
#endif
          nprms = (double *) EG_alloc(2*len*sizeof(double)+len*sizeof(int));
          if (nprms == NULL) {
            stat = EGADS_MALLOC;
            if (npts < ntri/2) EG_free(trix);
            goto cleanup;
          }
          cnt = (int *) &nprms[2*len];
          for (k = 0; k < len; k++) cnt[k] = 0;
          for (l = 0; l < ntri; l++)
            for (is = 0; is < 3; is++) {
              k = tris[3*l+is]-1;
              if (cnt[k] != 0) continue;
              if (ptype[k] >= 0) {
                t = 0.0;
                if (ptype[k] == 0) {
                  /* Node -- find Edge */
                  ie = EG_findEdge(body, &xyzs[3*k], is, &tris[3*l], &tric[3*l],
                                   ptype, pindex, matches[i], btess, &t);
                  if (ie == EGADS_NOTFOUND) continue;
                  if (ie <= EGADS_SUCCESS) {
                    printf(" Error: %d Cannot find match for tess %d / Node %d\n",
                           ie, i+1, pindex[k]);
                    if (npts < ntri/2) EG_free(trix);
                    goto cleanup;
                  }
                } else {
                  /* Edge */
                  ie = EG_lookupEdge(pindex[k], matches[i]);
                  if (ie == 0) {
                    printf(" Error: Cannot find match for tess %d / Edge %d\n",
                           i+1, pindex[k]);
                    if (npts < ntri/2) EG_free(trix);
                    goto cleanup;
                  }
                  t = btess->tess1d[ie-1].t[ptype[k]-1];
                }
                top  = btess->tess1d[ie-1].obj;
                stat = EG_getEdgeUV(objs[m-1], top, 0, t, &nprms[2*k]);
                if (stat == EGADS_TOPOERR) {
                  stat = EG_getEdgeUV(objs[m-1], top,  1, t, &nprms[2*k]);
                  if (stat != EGADS_SUCCESS) {
                    if (npts < ntri/2) EG_free(trix);
                    goto cleanup;
                  }
                  stat = EG_getEdgeUV(objs[m-1], top, -1, t, uv);
                  if (stat != EGADS_SUCCESS) {
                    if (npts < ntri/2) EG_free(trix);
                    goto cleanup;
                  }
                  x0[0] = sqrt((nprms[2*k  ]-uvs[2*l  ])*
                               (nprms[2*k  ]-uvs[2*l  ]) +
                               (nprms[2*k+1]-uvs[2*l+1])*
                               (nprms[2*k+1]-uvs[2*l+1]));
                  x0[1] = sqrt((uv[0]-uvs[2*l  ])*(uv[0]-uvs[2*l  ]) +
                               (uv[1]-uvs[2*l+1])*(uv[1]-uvs[2*l+1]));
                  if (x0[1] < x0[0]) {
                    nprms[2*k  ] = uv[0];
                    nprms[2*k+1] = uv[1];
                  }
                } else if (stat != EGADS_SUCCESS) {
                  if (npts < ntri/2) EG_free(trix);
                  goto cleanup;
                }
                cnt[k] = 1;
              } else {
                /* interior -- do inverse evaluation */
                coord[0] = xyzs[3*k  ];
                coord[1] = xyzs[3*k+1];
                coord[2] = xyzs[3*k+2];
                stat     = EG_invEvaluate(objs[m-1], coord, &nprms[2*k], xyz);
                if (stat != EGADS_SUCCESS) {
                  if (npts < ntri/2) EG_free(trix);
                  goto cleanup;
                }
                cnt[k] = 1;
              }
            }
          for (k = 0; k < len; k++)
            if (cnt[k] == 0) {
              printf(" Error: Face %d - Vertex %d/%d not hit!\n", m, k+1, len);
              if (npts < ntri/2) EG_free(trix);
              goto cleanup;
            }

          stat = EG_setTessFace(*tess, m, len, xyzs, nprms, ntri, trix);
          if (stat != EGADS_SUCCESS) {
            if (npts < ntri/2) EG_free(trix);
            goto cleanup;
          }
        }
        /* TFI flag */
        if (otess->tess2d[fmatch[i].objs[2*j+1]-1].tfi == 1) {
          btess->tess2d[m-1].tfi = 1;
          qints[m-1] = btess->tess2d[m-1].ntris/2;
        }
        
        /* do we transfer mixed info? */
        stat = EG_attributeRet(iTess[i], ".tessType", &aType, &aLen, &aInts,
                               &aReals, &aStr);
        if (stat == EGADS_SUCCESS) {
          if (aType == ATTRSTRING)
            if (strcmp(aStr, "Quad") == 0)
              qints[m-1] = btess->tess2d[m-1].ntris/2;
        } else {
          stat = EG_attributeRet(iTess[i], ".mixed", &aType, &aLen, &aInts,
                                 &aReals, &aStr);
          if (stat == EGADS_SUCCESS) {
            if ((aType == ATTRINT) && (aLen == nobj))
              qints[m-1] = aInts[fmatch[i].objs[2*j+1]-1];
          }
        }
        
        if (npts < ntri/2) EG_free(trix);
        EG_free(uvs);
        uvs = NULL;
      }
    }
    
    /* write the tessellation attributes? */
#ifndef LITE
    for (k = j = 0; j < nobj; j++)
      if (qints[j] != 0) k++;
    
    if (k == nobj) {
      stat = EG_attributeAdd(*tess, ".mixed", ATTRINT, nobj, qints, NULL, NULL);
      if (stat != EGADS_SUCCESS)
        printf(" EGADS Warning: EG_attributeAdd m = %d !\n", stat);
      stat = EG_attributeAdd(*tess, ".tessType", ATTRSTRING, 4, NULL, NULL,
                             "Quad");
      if (stat != EGADS_SUCCESS)
        printf(" EGADS Warning: EG_attributeAdd Q = %d !\n", stat);
    } else if (k != 0) {
      stat = EG_attributeAdd(*tess, ".mixed", ATTRINT, nobj, qints, NULL, NULL);
      if (stat != EGADS_SUCCESS) {
         printf(" EGADS Warning: EG_attributeAdd M = %d !\n", stat);
      } else {
        stat = EG_attributeAdd(*tess, ".tessType", ATTRSTRING, 5, NULL, NULL,
                               "Mixed");
        if (stat != EGADS_SUCCESS)
          printf(" EGADS Warning: EG_attributeAdd T = %d!\n", stat);
      }
    }
  }
#endif

  /* complete the tessellation with all shared Faces */
  stat = EG_finishTess(*tess, params);
  if (stat != EGADS_SUCCESS) goto cleanup;
  
  /* fix up Node positions */
  for (i = 0; i < btess->nEdge; i++) {
    for (j = 0; j < 2; j++) {
      if (btess->tess1d[i].nodes[j]              == 0) continue;
      if (nodes[btess->tess1d[i].nodes[j]-1].hit == 0) continue;
      k = j*(btess->tess1d[i].npts - 1);
      btess->tess1d[i].xyz[3*k  ] = nodes[btess->tess1d[i].nodes[j]-1].xyz[0];
      btess->tess1d[i].xyz[3*k+1] = nodes[btess->tess1d[i].nodes[j]-1].xyz[1];
      btess->tess1d[i].xyz[3*k+2] = nodes[btess->tess1d[i].nodes[j]-1].xyz[2];
    }
  }
  for (i = 0; i < btess->nFace; i++) {
    for (j = 0; j < btess->tess2d[i].npts; j++) {
      if (btess->tess2d[i].ptype[j]               != 0) continue;
      if (nodes[btess->tess2d[i].pindex[j]-1].hit == 0) continue;
      btess->tess2d[i].xyz[3*j  ] = nodes[btess->tess2d[i].pindex[j]-1].xyz[0];
      btess->tess2d[i].xyz[3*j+1] = nodes[btess->tess2d[i].pindex[j]-1].xyz[1];
      btess->tess2d[i].xyz[3*j+2] = nodes[btess->tess2d[i].pindex[j]-1].xyz[2];
    }
  }
  
cleanup:
  if (uvs   != NULL) EG_free(uvs);
  if (nprms != NULL) EG_free(nprms);
  if (mark  != NULL) EG_free(mark);
  if (objs  != NULL) EG_free(objs);
  if (nodes != NULL) EG_free(nodes);
  for (i = 0; i < 2*niTess; i++)
    if (matches[i].objs != NULL) EG_free(matches[i].objs);
  EG_free(matches);
  if (stat != EGADS_SUCCESS) {
    EG_deleteObject(*tess);
    *tess = NULL;
  }
  
  return stat;
}


#ifdef STANDALONE
void writeTris(FILE *fp, ego tess, ego body)
{
  int          i, j, n, stat, len, ntri;
  const int    *ptype, *pindex, *tris, *tric;
  const double *xyzs, *prms;
  
  stat = EG_getBodyTopos(body, NULL, FACE, &n, NULL);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getBodyTopos = %d\n", stat);
    return;
  }
  fprintf(fp, " %d\n", n);
  for (i = 1; i <= n; i++) {
    stat = EG_getTessFace(tess, i, &len, &xyzs, &prms, &ptype, &pindex, &ntri,
                          &tris, &tric);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos = %d\n", stat);
      return;
    }
    fprintf(fp, " %d  %d\n", len, ntri);
    for (j = 0; j < len; j++)
      fprintf(fp, " %lf %lf %lf\n", xyzs[3*j  ],xyzs[3*j+1], xyzs[3*j+2]);
    for (j = 0; j < ntri; j++)
      fprintf(fp, "  %d  %d  %d\n", tris[3*j  ],tris[3*j+1], tris[3*j+2]);
  }
}


int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  int    stat;
  double data[7], params[3] = {0.1, 0.01, 12.0};
  ego    context, body1, body2, tess1, tess2;
  FILE   *fp;
  
  /* create an EGADS context */
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_open return = %d\n", stat);
    return 1;
  }
  
  /* create a cylinder */
  data[0] = data[2] = data[3] = data[5] = 0.0;
  data[1] = -2.0;
  data[4] =  2.0;
  data[6] =  0.5;
  stat = EG_makeSolidBody(context, CYLINDER, data, &body1);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody 1 return = %d\n", stat);
    return 1;
  }
  data[1] = 2.0;
  data[4] = 4.0;
  stat = EG_makeSolidBody(context, CYLINDER, data, &body2);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeSolidBody 2 return = %d\n", stat);
    return 1;
  }
  
  stat = EG_makeTessBody(body1, params, &tess1);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeTessBody return = %d\n", stat);
    return 1;
  }
  stat = EG_extractTess(1, &tess1, body2, params, &tess2);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_extractTess return = %d\n", stat);
    return 1;
  }
  
  fp = fopen("extract.tris", "w");
  if (fp != NULL) {
    fprintf(fp, " 2\n");
    writeTris(fp, tess1, body1);
    writeTris(fp, tess2, body2);
    fclose(fp);
  }
  
  EG_deleteObject(tess2);
  EG_deleteObject(tess1);
  EG_deleteObject(body2);
  EG_deleteObject(body1);
  EG_close(context);
  
  return 0;
}
#endif
