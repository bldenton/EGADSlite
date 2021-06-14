#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_ibfibf.c,v 1.13 2020/10/20 05:45:48 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_ibfibf (
  INT_ nbface,
  INT_ *ibfin,
  INT_3D *inibf,
  INT_ *libfin,
  INT_3D **ibfibf)
{
  // Determine tria-face neighbor map. Only simply connected surfaces are
  // considered (only two tria-faces may share an edge).

  INT_ *ibfichk = NULL;
  INT_ *mchk = NULL;

  INT_ found, ibface, ibfn1, ibfn2, ibfn3, ichk,
      jbface, jbface2, jbfn1, jbfn2, jbfn3, jnode2, jnode3, loc, loc1, loc2,
      nchk;
  INT_ reorder = 0;
  INT_ status = 0;

  // allocate tria-face neighbor map

  *ibfibf = (INT_3D *) uvmap_realloc (&status, *ibfibf, (nbface+1)*sizeof(INT_3D));

  if (status) {
    uvmap_error_message ("*** ERROR 103507 : unable to allocate required memory ***");
    return 103507;
  }

  // initialize tria-face neighbor map

  memset (*ibfibf, 0, (nbface+1) * sizeof (INT_3D));

  // loop over primary tria-faces

  for (jbface = 1; jbface <= nbface; jbface++) {

    // loop over primary tria-face nodes

    for (jbfn1 = 0; jbfn1 <= 2; jbfn1++) {

      // if tria-face neighbor is not set then find it

      if ((*ibfibf)[jbface][jbfn1] == 0) {

        // set primary tria-face nodes to pair with neighbor tria-face

        jbfn2 = (jbfn1 < 2) ? jbfn1+1: 0;
        jbfn3 = (jbfn2 < 2) ? jbfn2+1: 0;

        jnode2 = inibf[jbface][jbfn2];
        jnode3 = inibf[jbface][jbfn3];

        // loop over all tria-faces attached to primary tria-face node

        loc1 = libfin[jnode3];
        loc2 = libfin[jnode3+1];

        found = 0;

        loc = loc1;

        while (loc < loc2 && found == 0) {

          // tria-face attached to node jnode3 of primary tria-face

          ibface = ibfin[loc];

          // skip if attached tria-face is the primary tria-face

          if (ibface != jbface) {

            // loop over attached tria-face nodes

            ibfn2 = 0;

            while (ibfn2 < 3 && found == 0) {

              // only consider orientation where neighbor tria-face node
              // matches node jnode3 of primary tria-face

              if (inibf[ibface][ibfn2] == jnode3) {

                // search other neighbor tria-face nodes for match to node
                // jnode2 of primary tria-face

                ibfn3 = (ibfn2 < 2) ? ibfn2+1: 0;
                ibfn1 = (ibfn3 < 2) ? ibfn3+1: 0;

                // if found then set tria-face neighbor map

                if (inibf[ibface][ibfn3] == jnode2) {

                  if (inibf[ibface][ibfn1] == inibf[jbface][jbfn1])  {
                    uvmap_error_message ("*** ERROR 3505 : duplicate tria-face found ***");
                    return 3505;
                  }

                  found = 1;

                  (*ibfibf)[jbface][jbfn1] = ibface;
                  (*ibfibf)[ibface][ibfn1] = jbface;
                }

                // if opposite ordering configuration is found then set
                // tria-face neighbor map and set the re-ordering flag

                else if (inibf[ibface][ibfn1] == jnode2) {

                  if (inibf[ibface][ibfn3] == inibf[jbface][jbfn1]) {
                    uvmap_error_message ("*** ERROR 3506 : duplicate tria-face found ***");
                    return 3506;
                  }

                  found = 1;

                  reorder = 1;

                  (*ibfibf)[jbface][jbfn1] = ibface;
                  (*ibfibf)[ibface][ibfn3] = jbface;
                }

                // if a matching pair was not found then skip testing other
                // nodes of this tria-face

                else
                  found = -1;
              }

              ibfn2++;
            }

            if (found < 0) found = 0;
          }

          loc++;
        }
      }
    }
  }

  // if inconsistent ordering was found then reorder tria-faces

  if (reorder) {

    // allocate tria-face neighbor map

    ibfichk = (INT_ *) uvmap_malloc (&status, (nbface+1)*sizeof(INT_));
    mchk = (INT_ *) uvmap_malloc (&status, (nbface+1)*sizeof(INT_));

    if (status) {
      uvmap_free (ibfichk);
      uvmap_free (mchk);
      uvmap_error_message ("*** ERROR 103508 unable to allocate required memory ***");
      return 103508;
    }

    // initialize tria-face flag

    mchk[0] = 1;

    for (ibface = 1; ibface <= nbface; ibface++) {
      mchk[ibface] = 0;
    }

    // set starting tria-face

    ibface = 1;

    mchk[ibface] = 1;

    ibfichk[0] = ibface;

    nchk = 1;

    // loop over check list
    // note that the check list grows within the loop until all tria-faces
    // have been checked

    for (ichk = 0; ichk < nchk; ichk++) {

      ibface = ibfichk[ichk];

      // check tria-face neighbors

      for (ibfn1 = 0; ibfn1 < 3; ibfn1++) {

        jbface = (*ibfibf)[ibface][ibfn1];

        // continue if tria-face neighbors hasn't been checked

        if (mchk[jbface] == 0) {

          // add neighbor to the check list

          ibfichk[nchk] = jbface;

          mchk[jbface] = 1;

          nchk++;

          ibfn2 = (ibfn1 < 2) ? ibfn1+1: 0;

          // check tria-face neighbor ordering

          jbfn1 = 0;

          found = 0;

          while (jbfn1 < 3 && found == 0) {

            if ((*ibfibf)[jbface][jbfn1] == ibface) {

              found = 1;

              jbfn2 = (jbfn1 < 2) ? jbfn1+1: 0;
              jbfn3 = (jbfn2 < 2) ? jbfn2+1: 0;

              jnode2 = inibf[jbface][jbfn2];

              // if ordering is not the same then reorder tria-face neighbor

              if (jnode2 == inibf[ibface][ibfn2]) {

                inibf[jbface][jbfn2] = inibf[jbface][jbfn3];
                inibf[jbface][jbfn3] = jnode2;

                jbface2 = (*ibfibf)[jbface][jbfn2];

                (*ibfibf)[jbface][jbfn2] = (*ibfibf)[jbface][jbfn3];
                (*ibfibf)[jbface][jbfn3] = jbface2;
              }
            }

            jbfn1++;
          }
        }
      }
    }

    for (ibface = 1; ibface <= nbface; ibface++) {

      if (mchk[ibface] == 0) {
        uvmap_error_message ("*** ERROR 3507 : found tria-faces that were not checked for ordering ***");
        status = 3507;
      }
    }

    uvmap_free (ibfichk);
    uvmap_free (mchk);

    if (status)
      return status;
  }

  // check tria-face ordering

  for (ibface = 1; ibface <= nbface; ibface++) {

    for (ibfn1 = 0; ibfn1 < 3; ibfn1++) {

      jbface = (*ibfibf)[ibface][ibfn1];

      if (jbface) {

        ibfn2 = (ibfn1 < 2) ? ibfn1+1: 0;
        ibfn3 = (ibfn2 < 2) ? ibfn2+1: 0;

        found = 0;

        while (jbfn1 < 3 && found == 0) {

          if ((*ibfibf)[jbface][jbfn1] == ibface) {

            found = 1;

            jbfn2 = (jbfn1 < 2) ? jbfn1+1: 0;
            jbfn3 = (jbfn2 < 2) ? jbfn2+1: 0;

            if (inibf[jbface][jbfn2] != inibf[ibface][ibfn3] ||
                inibf[jbface][jbfn3] != inibf[ibface][ibfn2]) {
              uvmap_error_message ("*** ERROR 3508 : tria-face ordering is wrong ***");
              return 3508;
            }
          }

          jbfn1++;
        }
      }
    }
  }

  return 0;
}
