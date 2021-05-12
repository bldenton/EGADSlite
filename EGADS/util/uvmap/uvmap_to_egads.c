#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_to_egads.c,v 1.12 2020/06/13 20:43:17 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_to_egads (
  INT_ nbface,
  INT_ nnode,
  INT_ *idibf,
  INT_3D *inibf,
  DOUBLE_2D *u,
  DOUBLE_3D *x,
  int *ntria,
  int *nvert,
  int **idmap,
  int **tria,
  double **uv,
  double **xyz)
{
  // Convert from AFLR to EGADS style arrays.

  INT_ i;
  INT_ err = 0;

  *ntria = (int) nbface;
  *nvert = (int) nnode;

  *idmap = NULL;
  *tria = NULL;
  *uv = NULL;
  *xyz = NULL;

  if (idibf) *idmap = (int *) uvmap_malloc (&err, ((*ntria)+1)*sizeof(int));
  *tria = (int *) uvmap_malloc (&err, 3*((*ntria)+1)*sizeof(int));
  if (u) *uv = (double *) uvmap_malloc (&err, 2*((*nvert)+1)*sizeof(double));
  if (x) *xyz = (double *) uvmap_malloc (&err, 3*((*nvert)+1)*sizeof(double));

  if (err) {
    uvmap_error_message ("*** ERROR 103517 : unable to allocate required memory ***");
    return 103517;
  }

  if (idibf) {
    for (i = 0; i < *ntria; i++) {
      (*idmap)[i] = (int) idibf[i+1];
    }
  }

  for (i = 0; i < *ntria; i++) {
    (*tria)[3*i]   = (int) inibf[i+1][0];
    (*tria)[3*i+1] = (int) inibf[i+1][1];
    (*tria)[3*i+2] = (int) inibf[i+1][2];
  }

  if (u) {
    for (i = 0; i < *nvert; i++) {
      (*uv)[2*i]   = u[i+1][0];
      (*uv)[2*i+1] = u[i+1][1];
    }
  }

  if (x) {
    for (i = 0; i < *nvert; i++) {
      (*xyz)[3*i]   = x[i+1][0];
      (*xyz)[3*i+1] = x[i+1][1];
      (*xyz)[3*i+2] = x[i+1][2];
    }
  }

  return 0;
}
