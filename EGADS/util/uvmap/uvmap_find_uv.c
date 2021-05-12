#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_find_uv.c,v 1.18 2020/07/04 22:41:39 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Find location of given UV coordinates.
--------------------------------------------------------------------------------

INT_ uvmap_find_uv (
  INT_ idef,
  double u_[2],
  void *ptr,
  INT_ *local_idef,
  INT_ *ibface,
  INT_ inode_[3],
  double s[3]);


INPUT ARGUMENTS
---------------

idef		Surface ID label.

u_		UV coordinate location to find (2 in length).

ptr		UV mapping data structure.


RETURN VALUE
------------

0		UV coordinate location was found.
-1		UV coordinate location was not found.
>0		An error occurred.


OUTPUT ARGUMENTS
----------------

local_idef	Local surface ID label of the tria-face that contains the given
		UV coordinates.
		If the surface idef is a virtual/composite surface then it is
		composed of one or more local surfaces with differing surface ID
		labels.
		If the surface idef is a standard surface then the local surface
		ID label found is set to idef. 

ibface		Tria-face index on surface idef of the tria-face that contains
		the given UV coordinates.

inode_		Node/Vertex of the tria-face that contains the given UV
		coordinates (3 in length).

s		Linear interpolation shape functions for the tria-face that
		contains the given UV coordinates (3 in length). For example,
		given data in array data stored at nodes/vertices of the surface
		mesh, the interpolated value at the location found can be
		determined from the following expression.

		  data_intp = s[0] * data[inode_[0]]
		            + s[1] * data[inode_[1]]
		            + s[2] * data[inode_[2]]

*/

INT_ uvmap_find_uv (
  INT_ idef,
  double u_[2],
  void *ptr,
  INT_ *local_idef,
  INT_ *ibface,
  INT_ inode_[3],
  double s[3])
{
  // Find location of given UV coordinates.

  uvmap_struct *uvmap_struct_ptr;

  INT_ *idibf = NULL;
  INT_ *msrch = NULL;
  INT_3D *ibfibf = NULL;
  INT_3D *inibf = NULL;

  DOUBLE_2D *u = NULL;

  INT_ found, ibface_save, index, inode, isrch, j, jbface, k, nbface;
  INT_ status = 0;

  double area[3], area_min, area_sum, du[3][2];
  double smin = 1.0e-12;
  double smin2 = 0.1;

  uvmap_struct_ptr = (uvmap_struct *) ptr;

  if (uvmap_struct_ptr == NULL) {
    uvmap_error_message ("*** ERROR 3503 mapping surface structure not set ***");
    return 3503;
  }

  // get data from UV mapping data structure

  status = uvmap_struct_get_entry (idef, &index, &isrch, ibface, &nbface,
                                   &idibf, &msrch, &inibf, &ibfibf, &u, 
                                   uvmap_struct_ptr);

  if (status)
    return status;

  // save starting tria-face index

  ibface_save = *ibface;

  // area coordinate search loop

  jbface = *ibface;

  do
  {
    *ibface = jbface;

    // set search flag

    msrch[*ibface] = isrch;

    // set UV coordinate deltas

    for (j = 0; j < 3; j++) {
      inode = inibf[*ibface][j];
      for (k = 0; k < 2; k++) {
        du[j][k] = u[inode][k] - u_[k];
      }
    }

    // set areas for area coordinates

    area[0] = du[1][0] * du[2][1] - du[1][1] * du[2][0];
    area[1] = du[2][0] * du[0][1] - du[2][1] * du[0][0];
    area[2] = du[0][0] * du[1][1] - du[0][1] * du[1][0];

    area_sum = area[0] + area[1] + area[2];

    // set minimum area

    area_min = MIN (area[0], area[1]);
    area_min = MIN (area[2], area_min);

    // check if tria-face contains the given UV coordinates

    found = (area_min + smin * area_sum >= 0.0) ? 1 : -2;

    // if not found then set next tria-face to search

    j = 0;

    while (j < 3 && found < -1) {

      if (area[j] < 0.0) {

        jbface = ibfibf[*ibface][j];

        found = (jbface > 0) ? ((msrch[jbface] != isrch) ? -1 : -2) : -3;
      }

      j++;
    }
  }
  while (found == -1);

  // if search is stuck at a boundary tria-face then check with a larger
  // tolerance

  if (found == -3 && smin2 > smin && area_min + smin2 * area_sum >= 0.0)
    found = 1;

  // if not found then do a brute force global search for containing tria-face

  if (found < 0)
  {
    found = -1;

    // loop over tria-faces

    *ibface = 0;

    do {

      (*ibface)++;

      if (msrch[*ibface] != isrch) {

        // set search flag

        msrch[*ibface] = isrch;

        // set UV coordinate deltas

        for (j = 0; j < 3; j++) {
          inode = inibf[*ibface][j];
          for (k = 0; k < 2; k++) {
            du[j][k] = u[inode][k] - u_[k];
          }
        }

        // set areas for area coordinates

        area[0] = du[1][0] * du[2][1] - du[1][1] * du[2][0];
        area[1] = du[2][0] * du[0][1] - du[2][1] * du[0][0];
        area[2] = du[0][0] * du[1][1] - du[0][1] * du[1][0];

        area_sum = area[0] + area[1] + area[2];

        // set minimum area

        area_min = MIN (area[0], area[1]);
        area_min = MIN (area[2], area_min);

        // check if tria-face contains the given UV coordinates
        // if search is at a boundary tria-face then also check with a larger
        // tolerance

        if (area_min + smin * area_sum >= 0.0 ||
            (smin2 > smin && area_min + smin2 * area_sum >= 0.0 &&
             ((area[0] == area_min && ibfibf[*ibface][0] <= 0) ||
              (area[1] == area_min && ibfibf[*ibface][1] <= 0) ||
              (area[2] == area_min && ibfibf[*ibface][2] <= 0))))
          found = 1;
      }
    }
    while (*ibface < nbface && found == -1);
  }

  // set search data in UV mapping data structure

  if (found == -1)
    *ibface = ibface_save;

  isrch++;

  uvmap_struct_set_srch_data (index, isrch, *ibface, uvmap_struct_ptr);

  // if found then for the containing tria-face set nodes/vertices, local
  // surface ID label, and shape-functions

  if (found == 1) {

    inode_[0] = inibf[*ibface][0];
    inode_[1] = inibf[*ibface][1];
    inode_[2] = inibf[*ibface][2];

    if (idibf)
      *local_idef = idibf[*ibface];
    else
      *local_idef = idef;

    s[0] = area[0] / (area[0] + area[1] + area[2]);
    s[1] = area[1] / (area[0] + area[1] + area[2]);
    s[2] = area[2] / (area[0] + area[1] + area[2]);
  }

  // if not found then set return value, default output argument values, and
  // reset tria-face starting index

  else {

    status = -1;

    *ibface = -1;

    inode_[0] = -1;
    inode_[1] = -1;
    inode_[2] = -1;

    *local_idef = -1;

    s[0] = -1.0;
    s[1] = -1.0;
    s[2] = -1.0;
  }

  return status;
}
