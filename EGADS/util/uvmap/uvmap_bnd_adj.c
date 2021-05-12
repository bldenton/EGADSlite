#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_bnd_adj.c,v 1.8 2020/06/04 02:32:57 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

void uvmap_bnd_adj (
  INT_ icc,
  INT_ nbedge,
  INT_ nnodei,
  INT_2D *ibeibe,
  INT_ *ibfibe,
  INT_ *ibfin,
  INT_ *iccibe,
  INT_2D *inibe,
  INT_3D *inibf,
  INT_ *libfin,
  INT_ *mben_disc,
  double *dumax,
  double urelaxb,
  double urelaxb_i,
  double w_ortho_i,
  DOUBLE_2D *u)
{
  // Do one pass of closed curve boundary-edge node uv mapping adjustment for
  // orthogonality.

  double cdfh, cdfhm, dusmin, du1, du2,
         du211, du212, du21s, du311, du312, du31s, du321, du322, du32s,
         u1sum, u2sum, w;

  INT_ i, ibedge, ibedge0, ibface, inode, inode1, inode2, inode3, isum,
      jbedge, jnodem, loc, loc1, loc2, urelaxb_it;
  INT_ n_urelaxb_it = 3;

  cdfh  = 0.5 * sqrt (3.0);
  cdfhm = 0.7 * cdfh;

  // set first boundary-edge on closed curve 

  ibedge0 = 0;
  ibedge  = 1;

  do {

    if (iccibe[ibedge] == icc)
      ibedge0 = ibedge;

    ibedge++;
  }
  while (ibedge <= nbedge && ibedge0 == 0);

  // adjust uv mapping on closed curve icc boundary-edge nodes for orthogonality

  ibedge = ibedge0;

  i = 1;

  do {

    inode = inibe[ibedge][0];

    u1sum = 0.0;
    u2sum = 0.0;

    isum = 0;

    loc1 = libfin[inode];
    loc2 = libfin[inode+1];

    // skip nodes that are attached to only one tria-face

    if (loc2-1 > loc1) {

      // add orthogonality contribution of all surrounding tria-faces

      for (loc = loc1; loc < loc2; loc++) {

        ibface = ibfin[loc];

        // if the closed curve is an inner curve and the tria-face
        // contains a hole node then skip contribution 

        if (icc > 1) {
          jnodem = MAX (inibf[ibface][0], inibf[ibface][1]);
          jnodem = MAX (jnodem, inibf[ibface][2]);
        }
        else
          jnodem = 0;

        if (jnodem <= nnodei) {

          if (inode == inibf[ibface][0]) {
            inode1 = inibf[ibface][0];
            inode2 = inibf[ibface][1];
            inode3 = inibf[ibface][2];
          }
          else if (inode == inibf[ibface][1]) {
            inode1 = inibf[ibface][1];
            inode2 = inibf[ibface][2];
            inode3 = inibf[ibface][0];
          }
          else {
            inode1 = inibf[ibface][2];
            inode2 = inibf[ibface][0];
            inode3 = inibf[ibface][1];
          }

          du211 = u[inode2][0] - u[inode1][0];
          du212 = u[inode2][1] - u[inode1][1];
          du311 = u[inode3][0] - u[inode1][0];
          du312 = u[inode3][1] - u[inode1][1];
          du321 = u[inode3][0] - u[inode2][0];
          du322 = u[inode3][1] - u[inode2][1];

          du21s = du211 * du211 + du212 * du212;
          du31s = du311 * du311 + du312 * du312;
          du32s = du321 * du321 + du322 * du322;

          dusmin = MIN (du21s, du31s);
          dusmin = w_ortho_i * w_ortho_i * dusmin;
          dusmin = MIN (dusmin, du32s);

          w = cdfhm * sqrt (dusmin / du32s);

          du1 = -w * du322;
          du2 = w * du321;

          u1sum = u1sum + 0.5 * (u[inode2][0] + u[inode3][0]) + du1;
          u2sum = u2sum + 0.5 * (u[inode2][1] + u[inode3][1]) + du2;

          isum++;
        }
      }

      // add orthogonality contribution at nodes that had contributions
      // from more than one tria-face

      if (isum > 1) {

        du1 = urelaxb * (u1sum / ((double) isum) - u[inode][0]);
        du2 = urelaxb * (u2sum / ((double) isum) - u[inode][1]);

        *dumax = MAX (MAX (fabs (du1), fabs (du2)), *dumax);

        u[inode][0] = u[inode][0] + du1;
        u[inode][1] = u[inode][1] + du2;
      }
    }

    ibedge = ibeibe[ibedge][0];

    i++;
  }
  while (i <= nbedge && ibedge != ibedge0);

  // add orthogonality contribution at nodes that are attached to only
  // one tria-face

  ibedge = ibedge0;

  i = 1;

  do {

    jbedge = ibeibe[ibedge][1];

    if (ibfibe[ibedge] == ibfibe[jbedge]) {

      inode  = inibe[ibedge][0];
      inode2 = inibe[ibedge][1];
      inode3 = inibe[jbedge][0];

      u[inode][0] = 0.5 * (u[inode2][0] + u[inode3][0])
                  + cdfh * (u[inode2][1] - u[inode3][1]);
      u[inode][1] = 0.5 * (u[inode2][1] + u[inode3][1])
                  + cdfh * (u[inode3][0] - u[inode2][0]);
    }

    ibedge = ibeibe[ibedge][0];

    i++;
  }
  while (i <= nbedge && ibedge != ibedge0);

  // add multi-pass boundary under relaxation to orthogonality adjustment

  for (urelaxb_it = 1; urelaxb_it <= n_urelaxb_it; urelaxb_it++) {

    // add secondary under relaxation to closed curve boundary-edge node
    // adjustment

    ibedge = ibedge0;

    i = 1;

    do {

      jbedge = ibeibe[ibedge][1];

      inode  = inibe[ibedge][0];
      inode2 = inibe[ibedge][1];
      inode1 = inibe[jbedge][0];

      // skip discontinuous boundary-edge nodes

      if (mben_disc[inode] == 0) {

        du1 = urelaxb_i * (0.5 * (u[inode1][0] + u[inode2][0]) - u[inode][0]);
        du2 = urelaxb_i * (0.5 * (u[inode1][1] + u[inode2][1]) - u[inode][1]);

        *dumax = MAX (MAX (fabs (du1), fabs (du2)), *dumax);

        u[inode][0] = u[inode][0] + du1;
        u[inode][1] = u[inode][1] + du2;
      }

      ibedge = ibeibe[ibedge][0];

      i++;
    }
    while (i <= nbedge && ibedge != ibedge0);
  }

  return;
}
