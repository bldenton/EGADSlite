#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_inibe.c,v 1.11 2020/06/13 20:43:17 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_inibe (
  INT_ nbface,
  INT_ *nbedge,
  INT_ **ibfibe,
  INT_3D *ibfibf,
  INT_3D *inibf,
  INT_2D **inibe)
{
  // Determine boundary-edge connectivity and tria-face boundary-edge map.

  INT_ ibedge, ibface, ibfn1, ibfn2, ibfn3;
  INT_ status = 0;

  // count number of boundary-edges

  *nbedge = 0;

  for (ibface = 1; ibface <= nbface; ibface++) {
    for (ibfn1 = 0; ibfn1 < 3; ibfn1++) {
      if (ibfibf[ibface][ibfn1] <= 0) (*nbedge)++;
    }
  }

  // allocate boundary-edge connectivity and tria-face boundary-edge map

  *ibfibe = (INT_ *) uvmap_malloc (&status, ((*nbedge)+1)*sizeof(INT_));
  *inibe = (INT_2D *) uvmap_malloc (&status, ((*nbedge)+1)*sizeof(INT_2D));

  if (status) {
    uvmap_error_message ("*** ERROR 103513 : unable to allocate required memory ***");
    return 103513;
  }

  // loop over tria-faces and set boundary-edges

  ibedge = 0;

  for (ibface = 1; ibface <= nbface; ibface++) {
    for (ibfn1 = 0; ibfn1 < 3; ibfn1++) {

      // assume that if there is no neighbor tria-face then the neighbor is a
      // boundary-edge

      if (ibfibf[ibface][ibfn1] <= 0) {

        ibfn2 = (ibfn1 < 2) ? ibfn1+1: 0;
        ibfn3 = (ibfn2 < 2) ? ibfn2+1: 0;

        ibedge++;

        ibfibf[ibface][ibfn1] = -ibedge;

        (*ibfibe)[ibedge] = ibface;

        (*inibe)[ibedge][0] = inibf[ibface][ibfn2];
        (*inibe)[ibedge][1] = inibf[ibface][ibfn3];
      }
    }
  }

  return 0;
}
