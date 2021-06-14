#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_inl_uv_bnd.c,v 1.4 2020/06/04 02:32:57 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

void uvmap_inl_uv_bnd (
  INT_ nbedge,
  INT_ nnode,
  INT_2D *ibeibe,
  INT_ *iccibe,
  INT_2D *inibe,
  DOUBLE_2D *u)
{
  // Generate a UV surface mesh given an XYZ surface mesh.

  INT_ i, ibedge, ibedge0, inode, nbedge0;

  double ang, dang;

  // set all uv's to zero

  for (inode = 1; inode <= nnode; inode++) {
    u[inode][0] = 0.5;
    u[inode][1] = 0.5;
  }

  // choose closed curve 1 to be the outer boundary-edge curve
  // count number of boundary-edges and set first edge on outer curve

  ibedge0 = 0;
  nbedge0 = 0;

  for (ibedge = 1; ibedge <= nbedge; ibedge++) {

    if (iccibe[ibedge] == 1) {

      nbedge0++;

      if (ibedge0 == 0)
        ibedge0 = ibedge;
    }
  }

  // set uv coordinates on outer curve to a circle

  dang = 8.0 * atan (1.0) / ((double) nbedge0);

  ibedge = ibedge0;

  i = 1;

  ang = 0.0;

  do {

    inode = inibe[ibedge][0];

    u[inode][0] = 0.5 * (cos (ang) + 1.0);
    u[inode][1] = 0.5 * (sin (ang) + 1.0);

    ibedge = ibeibe[ibedge][0];

    ang += dang;

    i++;
  }
  while (i <= nbedge && ibedge != ibedge0);

  return;
}
