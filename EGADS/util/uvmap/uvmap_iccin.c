#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_iccin.c,v 1.6 2020/06/13 20:43:17 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_iccin (
  INT_ nbedge,
  INT_ nnode,
  INT_ *iccibe,
  INT_ **iccin,
  INT_2D *inibe)
{
  // Set boundary-edge closed curve ID at each node.

  INT_ ibedge, inode1, inode2;
  INT_ status = 0;

  // allocate node flag

  *iccin = (INT_ *) uvmap_malloc (&status, (nnode+1)*sizeof(INT_));

  if (status) {
    uvmap_error_message ("*** ERROR 103511 : unable to allocate required memory ***");
    return 103511;
  }

  // initialize node flag to zero

  memset (*iccin, 0, (nnode+1) * sizeof (INT_));

  // loop over boundary-edges and set boundary-edge closed curve ID

  for (ibedge = 1; ibedge <= nbedge; ++ibedge) {

    inode1 = inibe[ibedge][0];
    inode2 = inibe[ibedge][1];

    (*iccin)[inode1] = iccibe[ibedge];
    (*iccin)[inode2] = iccibe[ibedge];
  }

  return 0;
}
