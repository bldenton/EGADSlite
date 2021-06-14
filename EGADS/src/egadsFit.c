/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             BSpline Fit Functions
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "egads.h"

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

extern /*@kept@*/ /*@null@*/ egObject *EG_context( const ego object );
extern int EG_outLevel( const ego object );
extern int EG_spline2dAppx( ego context,     int    endc,
                            /*@null@*/ const double *uknot,
                            /*@null@*/ const double *vknot,
                            /*@null@*/ const int    *vdata,
                            /*@null@*/ const double *wesT,
                            /*@null@*/ const double *easT,
                            /*@null@*/ const double *south,
                            /*@null@*/       double *snor,
                            /*@null@*/ const double *north,
                            /*@null@*/       double *nnor, int imax, int jmax,
                            const double *xyz, double tol,
                            ego *esurf );


int
EG_adjustCPs(const ego body, const ego face, double *CPs, ego *newBody,
             ego *newFace)
{
  int    i, j, k, stat, outLevel, oclass, mtype, nLoop, nEdge, iper, faceOr;
  int    *senses, *fsense, *header;
  double t, result[6], uvbox[4], trange[2], *rvec, *cntrl;
  ego    context, surf, ref, newSurf, loop, replace[2], *faces, *loops, *edges;

  *newBody = *newFace = NULL;
  if  (body == NULL)               return EGADS_NULLOBJ;
  if  (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (body->oclass != BODY)       return EGADS_NOTBODY;
  if ((body->mtype != SHEETBODY) &&
      (body->mtype != SOLIDBODY))  return EGADS_NOTTOPO;
  if  (body->blind == NULL)        return EGADS_NODATA;
  if  (face == NULL)               return EGADS_NULLOBJ;
  if  (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (face->oclass != FACE)       return EGADS_NOTTOPO;
  if  (face->blind == NULL)        return EGADS_NODATA;
  if  (CPs == NULL)                return EGADS_NODATA;
  outLevel = EG_outLevel(body);
  context  = EG_context(body);
  if  (context == NULL)            return EGADS_NOTCNTX;

  /* check the Topology */
  stat = EG_indexBodyTopo(body, face);
  if (stat < EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Face not in Body indexBodyTopo %d (EG_adjustCPs)!\n",
             stat);
    return stat;
  }
  stat = EG_getTopology(face, &surf, &oclass, &faceOr, uvbox, &nLoop, &loops,
                        &fsense);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: getTopology on Face = %d (EG_adjustCPs)!\n", stat);
    return stat;
  }
  if (nLoop != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Face has %d Loops (EG_adjustCPs)!\n", nLoop);
    return EGADS_TOPOERR;
  }
  if (surf->mtype != BSPLINE) {
    if (outLevel > 0)
      printf(" EGADS Error: Surf Geometry not a BSpline (EG_adjustCPs)!\n");
    return EGADS_GEOMERR;
  }
  stat = EG_getTopology(loops[0], &ref, &oclass, &mtype, NULL, &nEdge, &edges,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: getTopology on Loop = %d (EG_adjustCPs)!\n", stat);
    return stat;
  }
  if (nEdge != 4) {
    if (outLevel > 0)
      printf(" EGADS Error: Loop has %d Edges (EG_adjustCPs)!\n", nEdge);
    return EGADS_TOPOERR;
  }
  /* is the center of the Edge (via the PCurve) on an IsoCline? */
  for (j = i = 0; i < 4; i++) {
    stat = EG_getRange(edges[i], trange, &iper);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: getRange on Edge %d = %d (EG_adjustCPs)!\n",
               i+1, stat);
      return stat;
    }
    t    = trange[0] + (trange[1] - trange[0])/3.0;
    stat = EG_evaluate(edges[i+4], &t, result);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: evaluate on PCurve 1/3 %d = %d (EG_adjustCPs)!\n",
               i+1, stat);
      return stat;
    }
    if ((fabs(result[2]) > 1.e-10) && (fabs(result[3]) > 1.e-10)) j++;
/*  printf(" EGADS Info:  PCurve %d 1 -- dUV/dt = %le %le (EG_adjustCPs)!\n",
           i+1, result[2], result[3]);  */
    t    = trange[0] + 2.0*(trange[1] - trange[0])/3.0;
    stat = EG_evaluate(edges[i+4], &t, result);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: evaluate on PCurve 2/3 %d = %d (EG_adjustCPs)!\n",
               i+1, stat);
      return stat;
    }
    if ((fabs(result[2]) > 1.e-10) && (fabs(result[3]) > 1.e-10)) j++;
/*  printf(" EGADS Info:  PCurve %d 2 -- dUV/dt = %le %le (EG_adjustCPs)!\n",
           i+1, result[2], result[3]);  */
  }
  if (j != 0) {
    if (outLevel > 0)
      printf(" EGADS Error: PCurve(s) not on IsoCline(s) (EG_adjustCPs)!\n");
    return EGADS_NOTORTHO;
  }

  /* get the Geometry */
  stat = EG_getGeometry(surf, &oclass, &mtype, &ref, &header, &rvec);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: getGeometry on Surf = %d (EG_adjustCPs)!\n", stat);
    return stat;
  }

  /* refill the control points */
  cntrl = &rvec[header[3]+header[6]];
  for (k = j = 0; j < header[5]; j++) {
    for (i = 0; i < header[2]; i++, k++) {
      if ((i == 0) || (i == header[2]-1) || (j == 0) || (j == header[5]-1))
        if ((CPs[3*k  ] != cntrl[3*k  ]) || (CPs[3*k+1] != cntrl[3*k+1]) ||
            (CPs[3*k+2] != cntrl[3*k+2])) {
          printf(" EGADS Error: CP[%d,%d] = %lf %lf %lf -- %lf %lf %lf (EG_adjustCPs)!\n",
                 i+1, j+1, CPs[3*k], cntrl[3*k], CPs[3*k+1], cntrl[3*k+1],
                 CPs[3*k+2], cntrl[3*k+2]);
          EG_free(header);
          EG_free(rvec);
          return EGADS_CONSTERR;
        }
      cntrl[3*k  ] = CPs[3*k  ];
      cntrl[3*k+1] = CPs[3*k+1];
      cntrl[3*k+2] = CPs[3*k+2];
    }
  }
  if ((header[0]&2) != 0) {
    k *= 3;
    for (j = 0; j < header[5]; j++) {
      for (i = 0; i < header[2]; i++, k++) {
        if ((i == 0) || (i == header[2]-1) || (j == 0) || (j == header[5]-1))
          if (CPs[k] != cntrl[k]) {
            printf(" EGADS Error: W[%d,%d] = %lf -- %lf (EG_adjustCPs)!\n",
                   i+1, j+1, CPs[k], cntrl[k]);
            EG_free(header);
            EG_free(rvec);
            return EGADS_CONSTERR;
          }
        cntrl[k] = CPs[k];
      }
    }
  }

  /* make the new surface, Loop & Face */
  stat = EG_makeGeometry(context, oclass, mtype, ref, header, rvec, &newSurf);
  EG_free(header);
  EG_free(rvec);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeGeometry on newSurf = %d (EG_adjustCPs)!\n",
             stat);
    return stat;
  }
  stat = EG_makeTopology(context, newSurf, LOOP, CLOSED, NULL, 4, edges, senses,
                         &loop);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeTopology on newLoop = %d (EG_adjustCPs)!\n",
             stat);
    EG_deleteObject(newSurf);
    return stat;
  }
  stat = EG_makeTopology(context, newSurf, FACE, faceOr, NULL, 1, &loop, fsense,
                         &replace[1]);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeTopology on newFace = %d (EG_adjustCPs)!\n",
             stat);
    EG_deleteObject(loop);
    EG_deleteObject(newSurf);
    return stat;
  }

  /* replace the Face */
  replace[0] = face;
  stat = EG_replaceFaces(body, 1, replace, newBody);
  if ((stat != EGADS_SUCCESS) || (*newBody == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: replaceFaces = %d (EG_adjustCPs)!\n", stat);
    EG_deleteObject(replace[1]);
    EG_deleteObject(loop);
    EG_deleteObject(newSurf);
    return stat;
  }

  /* get newFace from newBody */
  stat = EG_getBodyTopos(*newBody, NULL, FACE, &k, &faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: getBodyTopos on newBody = %d (EG_adjustCPs)!\n",
             stat);
    EG_deleteObject(*newBody);
    EG_deleteObject(replace[1]);
    EG_deleteObject(loop);
    EG_deleteObject(newSurf);
    *newBody = NULL;
    return stat;
  }

  for (j = i = 0; i < k; i++) {
    stat = EG_isEquivalent(faces[i], replace[1]);
    if (stat < EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: isEquivalent on Face %d = %d (EG_adjustCPs)!\n",
               i+1, stat);
      EG_free(faces);
      EG_deleteObject(*newBody);
      EG_deleteObject(replace[1]);
      EG_deleteObject(loop);
      EG_deleteObject(newSurf);
      *newBody = NULL;
      return stat;
    }
    if (stat == EGADS_SUCCESS) {
      if (*newFace != NULL)
        printf(" EGADS Warning: Face multi-match %d %d (EG_adjustCPs)!\n",
               i+1, j);
      *newFace = faces[i];
      j = i+1;
    }
  }
  if (*newFace == NULL) {
    printf(" EGADS Warning: newFace not found in newBody (EG_adjustCPs)!\n");
  } else {
    /* copy the attributes to the newFace */
    stat = EG_attributeDup(face, *newFace);
    if (stat < EGADS_SUCCESS)
      if (outLevel > 0)
        printf(" EGADS Warning: attributeDup = %d (EG_adjustCPs)!\n", stat);
  }

  EG_free(faces);
  EG_deleteObject(replace[1]);
  EG_deleteObject(loop);
  EG_deleteObject(newSurf);

  return EGADS_SUCCESS;
}


/*
 ************************************************************************
 *                                                                      *
 *   EG_spline2d - create 2d cubic spline from input data               *
 *                                                                      *
 ************************************************************************
 */

int EG_spline2d(ego context, int endc, /*@null@*/ const double **drnd,
                int imax, int jmax, const double *xyz, double tol, ego *esurf)
{
  int          i, j, k, stat, outLevel;
  double       dist, dy, *norT, *souT, *easT, *wesT;
  double       x0[3], x1[3], nnor[3], snor[3], enor[3], wnor[3], *rot;
  const double *north, *south, *east, *west;
  
  *esurf = NULL;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if ((imax < 2) || (jmax < 2))      return EGADS_DEGEN;
  if ((endc < 0) || (endc > 2))      return EGADS_RANGERR;
  outLevel = EG_outLevel(context);
  
  /* check for degenerate sides */
  north = south = east = west = NULL;
  norT  = souT  = easT = wesT = NULL;
  i  = 0;
  dy = 0.0;
  for (j = 1; j < jmax; j++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  spline2d: Imin (west) degenerate!\n");
#endif
    if (drnd != NULL)
      if (drnd[0] != NULL) {
        west  = drnd[0];
        x0[0] = drnd[0][1];
        x0[1] = drnd[0][2];
        x0[2] = drnd[0][3];
        x1[0] = drnd[0][5];
        x1[1] = drnd[0][6];
        x1[2] = drnd[0][7];
        CROSS(wnor, x0, x1);
        dist  = DOT(wnor, wnor);
        if ((dist == 0.0) || (DOT(x0, x1) > tol) ||
            (fabs(1.0-DOT(x0, x0)) > tol) ||
            (fabs(1.0-DOT(x1, x1)) > tol)) {
          if (outLevel > 0)
            printf(" EGADS Error: BAD West Axes (EG_spline2d)!\n");
          return EGADS_NOTORTHO;
        }
        dist   = 1.0/sqrt(dist);
        wnor[0] *= dist;
        wnor[1] *= dist;
        wnor[2] *= dist;
        wesT     = wnor;
#ifdef DEBUG
        printf("            with normal = %lf %lf %lf!\n",
               wnor[0], wnor[1], wnor[2]);
#endif
      }
  }
  j  = 0;
  dy = 0.0;
  for (i = 1; i < imax; i++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  spline2d: Jmin (south) degenerate!\n");
#endif
    if (drnd != NULL)
      if (drnd[1] != NULL) {
        south = drnd[1];
        x0[0] = drnd[1][1];
        x0[1] = drnd[1][2];
        x0[2] = drnd[1][3];
        x1[0] = drnd[1][5];
        x1[1] = drnd[1][6];
        x1[2] = drnd[1][7];
        CROSS(snor, x0, x1);
        dist  = DOT(snor, snor);
        if ((dist == 0.0) || (DOT(x0, x1) > tol) ||
            (fabs(1.0-DOT(x0, x0)) > tol) ||
            (fabs(1.0-DOT(x1, x1)) > tol)) {
          if (outLevel > 0)
            printf(" EGADS Error: BAD South Axes (EG_spline2d)!\n");
          return EGADS_NOTORTHO;
        }
        dist   = 1.0/sqrt(dist);
        snor[0] *= dist;
        snor[1] *= dist;
        snor[2] *= dist;
        souT     = snor;
#ifdef DEBUG
        printf("            with normal = %lf %lf %lf!\n",
               snor[0], snor[1], snor[2]);
#endif
      }
  }
  i  = imax-1;
  dy = 0.0;
  for (j = 1; j < jmax; j++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  spline2d: Imax (east) degenerate!\n");
#endif
    if (drnd != NULL)
      if (drnd[2] != NULL) {
        east  = drnd[2];
        x0[0] = drnd[2][1];
        x0[1] = drnd[2][2];
        x0[2] = drnd[2][3];
        x1[0] = drnd[2][5];
        x1[1] = drnd[2][6];
        x1[2] = drnd[2][7];
        CROSS(enor, x0, x1);
        dist  = DOT(enor, enor);
        if ((dist == 0.0) || (DOT(x0, x1) > tol)||
            (fabs(1.0-DOT(x0, x0)) > tol) ||
            (fabs(1.0-DOT(x1, x1)) > tol)) {
          if (outLevel > 0)
            printf(" EGADS Error: BAD East Axes (EG_spline2d)!\n");
          return EGADS_NOTORTHO;
        }
        dist   = 1.0/sqrt(dist);
        enor[0] *= dist;
        enor[1] *= dist;
        enor[2] *= dist;
        easT     = enor;
#ifdef DEBUG
        printf("            with normal = %lf %lf %lf!\n",
               enor[0], enor[1], enor[2]);
#endif
      }
  }
  j  = jmax-1;
  dy = 0.0;
  for (i = 1; i < imax; i++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  spline2d: Jmax (north) degenerate!\n");
#endif
    if (drnd != NULL)
      if (drnd[3] != NULL) {
        north = drnd[3];
        x0[0] = drnd[3][1];
        x0[1] = drnd[3][2];
        x0[2] = drnd[3][3];
        x1[0] = drnd[3][5];
        x1[1] = drnd[3][6];
        x1[2] = drnd[3][7];
        CROSS(nnor, x0, x1);
        dist  = DOT(nnor, nnor);
        if ((dist == 0.0) || (DOT(x0, x1) > tol)||
            (fabs(1.0-DOT(x0, x0)) > tol) ||
            (fabs(1.0-DOT(x1, x1)) > tol)) {
          if (outLevel > 0)
            printf(" EGADS Error: BAD North Axes (EG_spline2d)!\n");
          return EGADS_NOTORTHO;
        }
        dist   = 1.0/sqrt(dist);
        nnor[0] *= dist;
        nnor[1] *= dist;
        nnor[2] *= dist;
        norT     = nnor;
#ifdef DEBUG
        printf("            with normal = %lf %lf %lf!\n",
               nnor[0], nnor[1], nnor[2]);
#endif
      }
  }
  if ((north != NULL) && ((east != NULL) || (west != NULL)))
    return EGADS_DEGEN;
  if ((south != NULL) && ((east != NULL) || (west != NULL)))
    return EGADS_DEGEN;
  
  /* call approx as is */
  if ((east == NULL) && (west == NULL))
    return EG_spline2dAppx(context, endc, NULL, NULL, NULL, NULL, NULL,
                           south, souT, north, norT, imax, jmax, xyz, tol,
                           esurf);
  
  /* rotate to get special treatment as north/south */
  rot = (double *) EG_alloc(3*imax*jmax*sizeof(double));
  if (rot == NULL) return EGADS_MALLOC;
  
  for (k = i = 0; i < imax; i++)
    for (j = jmax-1; j >= 0; j--, k++) {
      rot[3*k  ] = xyz[3*(i+j*imax)  ];
      rot[3*k+1] = xyz[3*(i+j*imax)+1];
      rot[3*k+2] = xyz[3*(i+j*imax)+2];
    }
#ifdef DEBUG
  printf(" spline2d: rotate -- south=west, north=east!\n");
#endif
  stat = EG_spline2dAppx(context, endc, NULL, NULL, NULL, NULL, NULL,
                         west, wesT, east, easT, jmax, imax, rot, tol,
                         esurf);
  
  EG_free(rot);
  return stat;
}
