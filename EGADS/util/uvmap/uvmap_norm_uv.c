#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_norm_uv.c,v 1.7 2020/06/04 02:32:57 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

void uvmap_norm_uv (
  INT_ nnode,
  INT_ nnodei,
  DOUBLE_2D *u)
{
  // Normalize uv coordinates so that all uv values are between zero and one.

  INT_ inode;

  double rdu, u1max, u1min, u2max, u2min;

  inode = 1;

  u1max = u[inode][0];
  u2max = u[inode][1];
  u1min = u[inode][0];
  u2min = u[inode][1];

  for (inode = 2; inode <= nnodei; inode++) {

    u1max = MAX (u1max, u[inode][0]);
    u2max = MAX (u2max, u[inode][1]);
    u1min = MIN (u1min, u[inode][0]);
    u2min = MIN (u2min, u[inode][1]);
  }

  rdu = 1.0 / MAX (u1max-u1min, u2max-u2min);

  for (inode = 1; inode <= nnode; inode++) {

    u[inode][0] = (u[inode][0] - u1min) * rdu;
    u[inode][1] = (u[inode][1] - u2min) * rdu;
  }

  return;
}
