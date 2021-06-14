#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_idibe.c,v 1.7 2020/10/20 05:45:48 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_idibe (
  INT_ nbedge,
  INT_ *nids, 
  INT_2D *ibeibe,
  INT_ **idibe)
{
  // Set boundary-edge closed curve ID.

  INT_ found, i, ibedge, id, jbedge, nbedgei;
  INT_ status = 0;

  // allocate boundary-edge closed curve ID

  *idibe = (INT_ *) uvmap_malloc (&status, (nbedge+1)*sizeof(INT_));

  if (status) {
    uvmap_error_message ("*** ERROR 103512 : unable to allocate required memory ***");
    return 103512;
  }

  // initialize boundary-edge closed curve ID

  memset (*idibe, 0, (nbedge+1) * sizeof (INT_));

  // loop over closed curves (there may only be one)

  id = 1;
  nbedgei = 0;

  do {

    // find a boundary-edge with an unset closed curve ID

    found = 0;
    ibedge = 1;

    while (ibedge <= nbedge && found == 0) {

      if ((*idibe)[ibedge] == 0)
        found = 1;
      else
        ibedge++;
    }

    // if a boundary-edge was found

    if (found) {

      // set boundary-edge closed curve ID for id and count number of edges
      // in curve

      (*idibe)[ibedge] = id;

      jbedge = ibedge;
      i = 1;

      do {

        nbedgei++;

        (*idibe)[ibedge] = id;

        ibedge = ibeibe[ibedge][0];

        i++;
      }
      while (i <= nbedge && ibedge != jbedge);
    }

    id++;
  }
  while (id <= nbedge && found == 1 && nbedgei < nbedge);

  // set number of closed curves

  *nids = id-1;

  return 0;
}

