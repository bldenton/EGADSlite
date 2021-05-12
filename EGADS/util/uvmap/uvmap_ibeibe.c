#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_ibeibe.c,v 1.11 2020/10/20 05:45:48 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_ibeibe (
  INT_ nbedge,
  INT_ nnode,
  INT_2D *inibe,
  INT_2D **ibeibe)
{
  // Determine boundary-edge neighbor map.
  // Only simply connected closed curve boundary edges are considered.

  INT_2D *ibein = NULL;

  INT_ ibedge, inode1, inode2, jbedge1, jbedge2;
  INT_ status = 0;

  // allocate boundary-edge neighbor map and temporary boundary-edge node map

  *ibeibe = (INT_2D *) uvmap_malloc (&status, (nbedge+1)*sizeof(INT_2D));
  ibein = (INT_2D *) uvmap_malloc (&status, (nnode+1)*sizeof(INT_2D));

  if (status) {
    uvmap_free (ibein);
    uvmap_error_message ("*** ERROR 103506 : unable to allocate required memory ***");
    return 103506;
  }

  // loop over nodes and initialize boundary-edge node map

  memset (ibein, 0, (nnode+1) * sizeof (INT_2D));

  // loop over boundary-edges and set boundary-edge node map

  for (ibedge = 1; ibedge <= nbedge; ibedge++) {

    inode1 = inibe[ibedge][0];
    inode2 = inibe[ibedge][1];

    ibein[inode1][0] = ibedge;
    ibein[inode2][1] = ibedge;
  }

  // loop over boundary-edges and set boundary-edge neighbor map

  for (ibedge = 1; ibedge <= nbedge; ibedge++) {

    inode1 = inibe[ibedge][0];
    inode2 = inibe[ibedge][1];

    (*ibeibe)[ibedge][0] = ibein[inode2][0];
    (*ibeibe)[ibedge][1] = ibein[inode1][1];
  }

  // free temporary boundary-edge node map

  uvmap_free (ibein);

  // check neighbor map

  for (ibedge = 1; ibedge <= nbedge; ibedge++) {

    jbedge1 = (*ibeibe)[ibedge][0];
    jbedge2 = (*ibeibe)[ibedge][1];

    if (inibe[jbedge1][0] != inibe[ibedge][1] ||inibe[jbedge2][1] != inibe[ibedge][0]) {
      uvmap_error_message ("*** ERROR 3504 : boundary-edge ordering is wrong ***");
      return 3504;
    }
  }

  return 0;
}
