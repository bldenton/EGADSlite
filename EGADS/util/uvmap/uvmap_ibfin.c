#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_ibfin.c,v 1.10 2020/06/13 20:43:17 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_ibfin (
  INT_ nbface,
  INT_ nnode,
  INT_ *nbfpnt,
  INT_3D *inibf,
  INT_ **ibfin,
  INT_ **libfin)
{
  // Determine a list of tria-faces attached to a node.

  INT_ ibface, inode, inode1, inode2, inode3, loc, n;
  INT_ status = 0;

  *nbfpnt = 0;

  // allocate location list
  // note that size is nnode+2 as there is one extra ending location stored

  *libfin = (INT_ *) uvmap_realloc (&status, *libfin, (nnode+2)*sizeof(INT_));

  if (status) {
    uvmap_error_message ("*** ERROR 103509 : unable to allocate required memory ***");
    return 103509;
  }

  // initialize location list

  memset (*libfin, 0, (nnode+2) * sizeof (INT_));

  // set location list equal to number of tria-faces attached to given node

  for (ibface = 1; ibface <= nbface; ibface++) {

    inode1 = inibf[ibface][0];
    inode2 = inibf[ibface][1];
    inode3 = inibf[ibface][2];

    (*libfin)[inode1]++;
    (*libfin)[inode2]++;
    (*libfin)[inode3]++;
  }

  // sum tria-faces attached to each node

  for (inode = 1; inode <= nnode; inode++) {
    *nbfpnt = (*nbfpnt) + (*libfin)[inode];
  }

  // allocate tria-face list

  *ibfin = (INT_ *) uvmap_realloc (&status, *ibfin, ((*nbfpnt)+1)*sizeof (INT_));

  if (status) {
    uvmap_error_message ("*** ERROR 103510 : unable to allocate required memory ***");
    return 103510;
  }

  // set location list equal to start of tria-face list for each node

  loc = 1;

  n = 0;

  for (inode = 1; inode <= nnode; inode++) {

    loc = loc + n;

    n = (*libfin)[inode];

    (*libfin)[inode] = loc;
  }

  // loop over tria-faces and set tria-face list and increment location list

  for (ibface = 1; ibface <= nbface; ibface++) {

    inode1 = inibf[ibface][0];
    inode2 = inibf[ibface][1];
    inode3 = inibf[ibface][2];

    loc = (*libfin)[inode1];

    (*ibfin)[loc] = ibface;

    loc++;

    (*libfin)[inode1] = loc;

    loc = (*libfin)[inode2];

    (*ibfin)[loc] = ibface;

    loc++;

    (*libfin)[inode2] = loc;

    loc = (*libfin)[inode3];

    (*ibfin)[loc] = ibface;

    loc++;

    (*libfin)[inode3] = loc;
  }

  // reset location list to start of tria-face list for each node
  // note that there is one extra ending location stored at end of list

  for (inode = nnode+1; inode >= 2; --inode) {
    (*libfin)[inode] = (*libfin)[inode-1];
  }

  (*libfin)[1] = 1;

  return 0;
}
