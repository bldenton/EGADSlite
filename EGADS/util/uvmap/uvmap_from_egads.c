#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_from_egads.c,v 1.12 2020/06/11 23:34:56 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_from_egads (
  int ntria,
  int nvert,
  int *idmap,
  int *tria,
  double *uv,
  double *xyz,
  INT_ *nbface,
  INT_ *nnode,
  INT_ **idibf,
  INT_3D **inibf,
  DOUBLE_2D **u,
  DOUBLE_3D **x)
{
  // Convert from EGADS to AFLR style arrays.

  INT_ i;
  INT_ status = 0;

  *nbface = (INT_) ntria;
  *nnode = (INT_) nvert;

  *idibf = NULL;
  *inibf = NULL;
  *u = NULL;
  *x = NULL;

  if (idmap) *idibf = (INT_ *) uvmap_malloc (&status, ((*nbface)+1)*sizeof(INT_));
  *inibf = (INT_3D *) uvmap_malloc (&status, ((*nbface)+1)*sizeof(INT_3D));
  if (uv) *u = (DOUBLE_2D *) uvmap_malloc (&status, ((*nnode)+1)*sizeof(DOUBLE_2D));
  if (xyz) *x = (DOUBLE_3D *) uvmap_malloc (&status, ((*nnode)+1)*sizeof(DOUBLE_3D));

  if (status) {
    uvmap_error_message ("*** ERROR 103503 : unable to allocate required memory ***");
    return 103503;
  }

  if (idmap) {
    for (i = 0; i < *nbface; i++) {
      (*idibf)[i+1] = (INT_) idmap[i];
    }
  }

  for (i = 0; i < *nbface; i++) {
    (*inibf)[i+1][0] = (INT_) tria[3*i];
    (*inibf)[i+1][1] = (INT_) tria[3*i+1];
    (*inibf)[i+1][2] = (INT_) tria[3*i+2];
  }

  if (uv) {
    for (i = 0; i < *nnode; i++) {
      (*u)[i+1][0] = uv[2*i];
      (*u)[i+1][1] = uv[2*i+1];
    }
  }

  if (xyz) {
    for (i = 0; i < *nnode; i++) {
      (*x)[i+1][0] = xyz[3*i];
      (*x)[i+1][1] = xyz[3*i+1];
      (*x)[i+1][2] = xyz[3*i+2];
    }
  }

  return 0;
}
