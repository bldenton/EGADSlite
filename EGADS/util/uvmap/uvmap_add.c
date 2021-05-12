#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_add.c,v 1.8 2020/06/11 23:34:56 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_add (
  INT_ nbedge,
  INT_ ncc, 
  INT_ *nbface,
  INT_ *nnode,
  INT_ *iccibe,
  INT_2D *inibe,
  INT_3D **inibf,
  DOUBLE_3D **x)
{
  // Add a node at the centroid of all inner closed boundary-edge curves and
  // add tria-faces that connect to the centroid node and the nodes of the
  // associated closed boundary-edge curve.

  INT_ icc, ibedge, ibface, inode, inode1, inode2, isum;
  INT_ status = 0;

  double x1sum, x2sum, x3sum;

  if (ncc == 1)
    return 0;

  // reallocate tria-connectivity and xyz coordinates for addition of nodes

  *inibf = (INT_3D *) uvmap_realloc (&status, *inibf, ((*nbface)+nbedge+1)*sizeof (INT_3D));
  if (*x) *x = (DOUBLE_3D *) uvmap_realloc (&status, *x, ((*nnode)+ncc)*sizeof (DOUBLE_3D));

  if (status) {
    uvmap_error_message ("*** ERROR 103500 : unable to allocate required memory ***");
    return 103500;
  }

  // loop over inner closed boundary-edge curves

  ibface = *nbface;
  inode = *nnode;

  for (icc = 2; icc <= ncc; icc++) {

    inode++;

    // loop over boundary-edges of closed curve and add tria-faces

    for (ibedge = 1; ibedge <= nbedge; ibedge++) {

      if (iccibe[ibedge] == icc) {

        // add new tria-face for each boundary-edge

        ibface++;

        inode1 = inibe[ibedge][0];
        inode2 = inibe[ibedge][1];

        (*inibf)[ibface][0] = inode2;
        (*inibf)[ibface][1] = inode1;
        (*inibf)[ibface][2] = inode;
      }
    }

    // loop over boundary-edges of closed curve and determine new coordinate

    if (*x) {

      isum  = 0;
      x1sum = 0.0;
      x2sum = 0.0;
      x3sum = 0.0;

      for (ibedge = 1; ibedge <= nbedge; ibedge++) {

        if (iccibe[ibedge] == icc) {

          // sum boundary-edge coordinates

          isum++;

          inode1 = inibe[ibedge][0];

          x1sum = x1sum + (*x)[inode1][0];
          x2sum = x2sum + (*x)[inode1][1];
          x3sum = x3sum + (*x)[inode1][2];
        }
      }

      // add new node at centroid of closed curve

      (*x)[inode][0] = x1sum / ((double) isum);
      (*x)[inode][1] = x2sum / ((double) isum);
      (*x)[inode][2] = x3sum / ((double) isum);
    }
  }

  // set new number of tria-faces and nodes

  *nbface = ibface;
  *nnode = inode;

  return 0;
}
