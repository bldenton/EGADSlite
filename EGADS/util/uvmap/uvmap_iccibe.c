#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_iccibe.c,v 1.11 2020/06/09 19:56:01 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_iccibe (
  INT_ nbedge,
  INT_ *ncc, 
  INT_2D *ibeibe,
  INT_ **iccibe)
{
  // Set boundary-edge closed curve ID.

  INT_ ibedge, icc, jcc, nbedgei;
  INT_ nbedgei_max = 0;
  INT_ status = 0;

  // allocate boundary-edge closed curve ID

  status = uvmap_idibe (nbedge, ncc, ibeibe, iccibe);

  if (status)
    return status;

  // set closed curve ID with the largest number of edges to one

  jcc = 1;

  for (icc = 1; icc <= *ncc; icc++) {

    nbedgei = 0;

    for (ibedge = 1; ibedge <= nbedge; ibedge++) {
      if ((*iccibe)[ibedge] == icc) ++nbedgei;
    }

    if (nbedgei > nbedgei_max) {
      jcc = icc;
      nbedgei_max = nbedgei;
    }
  }

  if (jcc != 1) {

    for (ibedge = 1; ibedge <= nbedge; ibedge++) {
      if ((*iccibe)[ibedge] == 1)
        (*iccibe)[ibedge] = -jcc;
      else if ((*iccibe)[ibedge] == jcc)
        (*iccibe)[ibedge] = -1;
    }

    for (ibedge = 1; ibedge <= nbedge; ibedge++) {
      if ((*iccibe)[ibedge] < 0)
        (*iccibe)[ibedge] = -((*iccibe)[ibedge]);
    }
  }

  return 0;
}

