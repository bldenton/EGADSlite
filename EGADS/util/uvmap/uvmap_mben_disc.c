#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_mben_disc.c,v 1.6 2020/06/13 20:43:17 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_mben_disc (
  INT_ nbedge,
  INT_ nnode,
  INT_2D *ibeibe,
  INT_2D *inibe, 
  INT_ **mben_disc,
  double angdbe,
  DOUBLE_3D *x)
{
  // Flag discontinuous (set to 1) boundary-edge nodes.

  INT_ ibedge, inode1, inode2, jbedge, jnode2;
  INT_ status = 0;

  double cosmin, cosmins, dsi, dsj, dx1i, dx2i, dx3i, dx1j, dx2j, dx3j, w;

  // set cosine for discontinuity angle

  cosmin = cos (angdbe * atan (1.0) / 45.0);
  cosmins = cosmin * fabs (cosmin);

  // allocate node flag

  *mben_disc = (INT_ *) uvmap_malloc (&status, (nnode+1)*sizeof(INT_));

  if (status) {
    uvmap_error_message ("*** ERROR 103514 : unable to allocate required memory ***");
    return 103514;
  }

  // initialize node flag to zero

  memset (*mben_disc, 0, (nnode+1) * sizeof (INT_));

  // return if coordinates are not set

  if (x == NULL)
    return 0;

  // loop over boundary-edges

  for (ibedge = 1; ibedge <= nbedge; ++ibedge) {

    jbedge = ibeibe[ibedge][0];

    inode1 = inibe[ibedge][0];
    inode2 = inibe[ibedge][1];
    jnode2 = inibe[jbedge][1];

    // set edge vectors for boundary-edge and neighbor

    dx1i = x[inode2][0] - x[inode1][0];
    dx2i = x[inode2][1] - x[inode1][1];
    dx3i = x[inode2][2] - x[inode1][2];
    dx1j = x[jnode2][0] - x[inode2][0];
    dx2j = x[jnode2][1] - x[inode2][1];
    dx3j = x[jnode2][2] - x[inode2][2];

    dsi = dx1i * dx1i + dx2i * dx2i + dx3i * dx3i;
    dsj = dx1j * dx1j + dx2j * dx2j + dx3j * dx3j;

    // if dot product between edge vectors has an angle greater than the
    // discontinuity angle then set node flag to one

    w = dx1i * dx1j + dx2i * dx2j + dx3i * dx3j;

    if (w * fabs (w) < cosmins * dsi * dsj)
      (*mben_disc)[inode2] = 1;
  }

  return 0;
}

