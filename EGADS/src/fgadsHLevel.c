/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             FORTRAN Bindings for High Level Functions
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>

#define INT8 unsigned long long

#include "egadsTypes.h"
#include "egadsInternals.h"

  extern int EG_fuseSheets(const egObject *src, const egObject *tool,
                           egObject **sheet);
  extern int EG_generalBoolean(egObject *src, egObject *tool, int oper,
                               double tol, egObject **model);
  extern int EG_solidBoolean(const egObject *src, const egObject *tool,
                             int oper, egObject **model);
  extern int EG_intersection(const egObject *src, const egObject *tool,
                             int *nedge, /*@null@*/ egObject ***facEdg,
                             egObject **model);
  extern int EG_imprintBody(const egObject *src, int nedge, 
                            egObject **facEdg, egObject **result);
  extern int EG_filletBody(const egObject *src, int nedge, 
                           const egObject **edges, double radius,
                                 egObject **result, int **faceMap);
  extern int EG_chamferBody(const egObject *src, int nedge, 
                            const egObject **edges, const egObject **faces,
                            double dis1, double dis2, egObject **result,
                            int **faceMap);
  extern int EG_hollowBody(const egObject *src, int nface, 
                           const egObject **faces, double off, int join,
                                 egObject **result, int **faceMap);
  extern int EG_extrude(const egObject *src, double dist, 
                        const double *dir, egObject **result);
  extern int EG_rotate(const egObject *src, double angle, 
                       const double *axis, egObject **result);
  extern int EG_sweep(const egObject *src, const egObject *spine, int mode,
                            egObject **result);
  extern int EG_loft(int nsec, const egObject **secs, int opt, 
                                     egObject **result);
  extern int EG_blend(int nsec, const egObject **secs, /*@null@*/ double *rc1,
                      /*@null@*/ double *rcN, ego *result);
  extern int EG_ruled(int nsec, const egObject **secs, egObject **result);


int
#ifdef WIN32
IG_FUSESHEETS (INT8 *isrc, INT8 *itool, INT8 *isheet)
#else
ig_fusesheets_(INT8 *isrc, INT8 *itool, INT8 *isheet)
#endif
{
  int      stat;
  egObject *src, *tool, *sheet;
  
  *isheet = 0;
  src     = (egObject *) *isrc;
  tool    = (egObject *) *itool;
  stat    = EG_fuseSheets(src, tool, &sheet);
  if (stat == EGADS_SUCCESS) *isheet = (INT8) sheet;
  return stat;
}


int
#ifdef WIN32
IG_GENERALBOOLEAN (INT8 *isrc, INT8 *itool, int *oper, double *tol,
                   INT8 *imodel)
#else
ig_generalboolean_(INT8 *isrc, INT8 *itool, int *oper, double *tol,
                   INT8 *imodel)
#endif
{
  int      stat;
  egObject *src, *tool, *model;
  
  *imodel = 0;
  src     = (egObject *) *isrc;
  tool    = (egObject *) *itool;
  stat    = EG_generalBoolean(src, tool, *oper, *tol, &model);
  if (stat == EGADS_SUCCESS) *imodel = (INT8) model;
  return stat;
}


int
#ifdef WIN32
IG_SOLIDBOOLEAN (INT8 *isrc, INT8 *itool, int *oper, INT8 *imodel)
#else
ig_solidboolean_(INT8 *isrc, INT8 *itool, int *oper, INT8 *imodel)
#endif
{
  int      stat;
  egObject *src, *tool, *model;
  
  *imodel = 0;
  src     = (egObject *) *isrc;
  tool    = (egObject *) *itool;
  stat    = EG_solidBoolean(src, tool, *oper, &model);
  if (stat == EGADS_SUCCESS) *imodel = (INT8) model;
  return stat;
}


int
#ifdef WIN32
IG_INTERSECTION (INT8 *isrc, INT8 *itool, int *nedge, 
                 INT8 **facedg8, INT8 *imodel)
#else
ig_intersection_(INT8 *isrc, INT8 *itool, int *nedge, 
                 INT8 **facedg8, INT8 *imodel)
#endif
{
  int      i, stat;
  INT8     *cobjs;
  egObject *src, *tool, **facEdg, *model;
  
  *nedge  = 0;
  *imodel = 0;
  src     = (egObject *) *isrc;
  tool    = (egObject *) *itool;
  stat    = EG_intersection(src, tool, nedge, &facEdg, &model);
  if (stat == EGADS_SUCCESS) {
    *facedg8 = cobjs = (INT8 *) EG_alloc(*nedge*2*sizeof(INT8));
    if (cobjs == NULL) {
      EG_free(facEdg);
      return EGADS_MALLOC;
    }
    for (i = 0; i < *nedge*2; i++) cobjs[i] = (INT8) facEdg[i];
    EG_free(facEdg);
    *imodel = (INT8) model;
  }
  return stat;
}


int
#ifdef WIN32
IG_IMPRINTBODY (INT8 *isrc, int *nedge, INT8 *facEdg, INT8 *irslt)
#else
ig_imprintbody_(INT8 *isrc, int *nedge, INT8 *facEdg, INT8 *irslt)
#endif
{
  int      i, stat;
  egObject *src, *result, **objs = NULL;
  
  *irslt = 0;
  if (*nedge <= 0) return EGADS_RANGERR;
  src  = (egObject *) *isrc;
  objs = (egObject **) EG_alloc(*nedge*2*sizeof(egObject *));
  if (objs == NULL) return EGADS_MALLOC;
  for (i = 0; i < *nedge*2; i++)
    objs[i] = (egObject *) facEdg[i];

  stat = EG_imprintBody(src, *nedge, objs, &result);
  if (objs != NULL) EG_free(objs);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}


int
#ifdef WIN32
IG_FILLETBODY (INT8 *isrc, int *nedge, INT8 *edges, double *radius,
               INT8 *irslt, int **facemap)
#else
ig_filletbody_(INT8 *isrc, int *nedge, INT8 *edges, double *radius,
               INT8 *irslt, int **facemap)
#endif
{
  int            i, stat;
  egObject       *src, *result;
  const egObject **objs = NULL;
  
  *irslt = 0;
  if (*nedge <= 0) return EGADS_RANGERR;
  src  = (egObject *) *isrc;
  objs = (const egObject **) EG_alloc(*nedge*sizeof(egObject *));
  if (objs == NULL) return EGADS_MALLOC;
  for (i = 0; i < *nedge; i++)
    objs[i] = (egObject *) edges[i];

  stat = EG_filletBody(src, *nedge, objs, *radius, &result, facemap);
  if (objs != NULL) EG_free((void *) objs);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}


int
#ifdef WIN32
IG_CHAMFERBODY (INT8 *isrc, int *nedge, INT8 *edges, INT8 *faces,
                double *dis1, double *dis2, INT8 *irslt, int **facemap)
#else
ig_chamferbody_(INT8 *isrc, int *nedge, INT8 *edges, INT8 *faces,
                double *dis1, double *dis2, INT8 *irslt, int **facemap)
#endif
{
  int            i, stat;
  egObject       *src, *result;
  const egObject **objs = NULL;
  
  *irslt = 0;
  if (*nedge <= 0) return EGADS_RANGERR;
  src  = (egObject *) *isrc;
  objs = (const egObject **) EG_alloc(*nedge*2*sizeof(egObject *));
  if (objs == NULL) return EGADS_MALLOC;
  for (i = 0; i < *nedge; i++) {
    objs[       i] = (egObject *) edges[i];
    objs[*nedge+i] = (egObject *) faces[i];
  }

  stat = EG_chamferBody(src, *nedge, objs, &objs[*nedge], 
                        *dis1, *dis2, &result, facemap);
  if (objs != NULL) EG_free((void *) objs);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}


int
#ifdef WIN32
IG_HOLLOWBODY (INT8 *isrc, int *nface, INT8 *faces, double *offset,
               int *join, INT8 *irslt, int **facemap)
#else
ig_hollowbody_(INT8 *isrc, int *nface, INT8 *faces, double *offset,
               int *join, INT8 *irslt, int **facemap)
#endif
{
  int            i, stat;
  egObject       *src, *result;
  const egObject **objs = NULL;
  
  *irslt = 0;
  if (*nface <= 0) return EGADS_RANGERR;
  src  = (egObject *) *isrc;
  objs = (const egObject **) EG_alloc(*nface*sizeof(egObject *));
  if (objs == NULL) return EGADS_MALLOC;
  for (i = 0; i < *nface; i++)
    objs[i] = (egObject *) faces[i];

  stat = EG_hollowBody(src, *nface, objs, *offset, *join, &result, facemap);
  if (objs != NULL) EG_free((void *) objs);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}


int
#ifdef WIN32
IG_EXTRUDE (INT8 *isrc, double *dist, double *dir, INT8 *irslt)
#else
ig_extrude_(INT8 *isrc, double *dist, double *dir, INT8 *irslt)
#endif
{
  int      stat;
  egObject *src, *result;
  
  *irslt = 0;
  src    = (egObject *) *isrc;
  stat   = EG_extrude(src, *dist, dir, &result);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}


int
#ifdef WIN32
IG_ROTATE (INT8 *isrc, double *angle, double *axis, INT8 *irslt)
#else
ig_rotate_(INT8 *isrc, double *angle, double *axis, INT8 *irslt)
#endif
{
  int      stat;
  egObject *src, *result;
  
  *irslt = 0;
  src    = (egObject *) *isrc;
  stat   = EG_rotate(src, *angle, axis, &result);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}


int
#ifdef WIN32
IG_SWEEP (INT8 *isrc, INT8 *ispine, int *mode, INT8 *irslt)
#else
ig_sweep_(INT8 *isrc, INT8 *ispine, int *mode, INT8 *irslt)
#endif
{
  int      stat;
  egObject *src, *spine, *result;
  
  *irslt = 0;
  src    = (egObject *) *isrc;
  spine  = (egObject *) *ispine;
  stat   = EG_sweep(src, spine, *mode, &result);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}


int
#ifdef WIN32
IG_LOFT (int *nsec, INT8 *secs, int *opt, INT8 *irslt)
#else
ig_loft_(int *nsec, INT8 *secs, int *opt, INT8 *irslt)
#endif
{
  int            i, stat;
  egObject       *result;
  const egObject **objs = NULL;
  
  *irslt = 0;
  if (*nsec <= 1) return EGADS_RANGERR;
  objs = (const egObject **) EG_alloc(*nsec*sizeof(egObject *));
  if (objs == NULL) return EGADS_MALLOC;
  for (i = 0; i < *nsec; i++)
    objs[i] = (egObject *) secs[i];

  stat = EG_loft(*nsec, objs, *opt, &result);
  if (objs != NULL) EG_free((void *) objs);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}


int
#ifdef WIN32
IG_BLEND (int *nsec, INT8 *secs, double *rc1, double *rcN, INT8 *irslt)
#else
ig_blend_(int *nsec, INT8 *secs, double *rc1, double *rcN, INT8 *irslt)
#endif
{
  int            i, stat;
  double         *ptr1, *ptrN;
  egObject       *result;
  const egObject **objs = NULL;
  
  *irslt = 0;
  if (*nsec <= 1) return EGADS_RANGERR;
  objs = (const egObject **) EG_alloc(*nsec*sizeof(egObject *));
  if (objs == NULL) return EGADS_MALLOC;
  for (i = 0; i < *nsec; i++)
    objs[i] = (egObject *) secs[i];
  ptr1 = rc1;
  ptrN = rcN;
  if (rc1[0] <= 0.0) ptr1 = NULL;
  if (rcN[0] <= 0.0) ptrN = NULL;

  stat = EG_blend(*nsec, objs, ptr1, ptrN, &result);
  if (objs != NULL) EG_free((void *) objs);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}


int
#ifdef WIN32
IG_RULED (int *nsec, INT8 *secs, INT8 *irslt)
#else
ig_ruled_(int *nsec, INT8 *secs, INT8 *irslt)
#endif
{
  int            i, stat;
  egObject       *result;
  const egObject **objs = NULL;
  
  *irslt = 0;
  if (*nsec <= 1) return EGADS_RANGERR;
  objs = (const egObject **) EG_alloc(*nsec*sizeof(egObject *));
  if (objs == NULL) return EGADS_MALLOC;
  for (i = 0; i < *nsec; i++)
    objs[i] = (egObject *) secs[i];
  
  stat = EG_ruled(*nsec, objs, &result);
  if (objs != NULL) EG_free((void *) objs);
  if (stat == EGADS_SUCCESS) *irslt = (INT8) result;
  return stat;
}
