#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: EG_uvmapTest.c,v 1.3 2021/04/30 18:55:42 marcum Exp $
 * Copyright 1994-2021, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Test uvmap with EGADS APIs.
--------------------------------------------------------------------------------

int EG_uvmapTest (char Case_Name[], int verbosity);


INPUT ARGUMENTS
---------------

Case_Name	Case name for input and output SURF type surface mesh files.

verbosity	Message flag.
		If verbosity=0 then do not output progress information.
		If verbosity=1 then output progress information to stdout.
		If verbosity=2 then output progress and additional CPU usage
		information to stdout.


RETURN VALUE
------------

EGADS_SUCCESS	Normal completion without errors.
EGADS_MALLOC	Unable to allocate required memory.
EGADS_UVMAP	An error occurred.

*/

int EG_uvmapTest (char Case_Name[], int verbosity)
{
  // Test uvmap with EGADS APIs.

  void *ptr = NULL;

  int *local_idef = NULL;
  int *tria = NULL;
  double *uv = NULL;
  double *xyz = NULL;

  int idef = 1;
  int ntria = 0;
  int nvert = 0;
  int set_struct = 1;
  int status = 0;

  // interpolation testing variables

  char Text[512];

  int ivert[3];
  int i, j, k, itriai, local_idefi;
  int location_found = 0;
  int correct_shape_function = 0;
  int correct_tria = 0;
  int correct_vertex = 0;
  int correct_vertices = 0;

  double uvi[2], s[3];
  double tol = 1.0e-6;

  // read xyz surface mesh

  status = EG_uvmap_Read (Case_Name, &ntria, &nvert, &local_idef, &tria, &xyz);

  // generate uv mapping

  if (status == 0)
    status = EG_uvmapGen (idef, ntria, nvert, set_struct, verbosity,
                           local_idef, tria, xyz, &uv, &ptr);

  // check interpolation

  if (status == 0) {

    uvmap_message ("");
    uvmap_message ("UVMAP    : CHECKING INTERPOLATION");
    uvmap_message ("");

    // loop over tria-faces and check interpolation at tria-face centroids

    i = 0;

    while (i < ntria && (status == EGADS_SUCCESS || status == EGADS_NOTFOUND)) {

      // set UV coordinates of interpolation location equal to the centroid of
      // tria-face

      for (j = 0; j < 2; j++) {
        uvi[j] = (uv[2*(tria[3*i  ]-1)+j]
                + uv[2*(tria[3*i+1]-1)+j]
                + uv[2*(tria[3*i+2]-1)+j]) / 3.0;
      }

      // find interpolation location

      status = EG_uvmapFindUV (idef, uvi, ptr, &local_idefi, &itriai, ivert, s); 

      // if location was found then check interpolation results

      if (status == EGADS_SUCCESS) {

        location_found++;

        // check if tria-face index is correct

        if (itriai == i+1)
          correct_tria++;

        // check if nodes/vertices are correct
        // note that all orientations are checked as the tria-face connectivity
        // within the UV mapping structure may be reordered for ordering
        // consistency

        correct_vertex = 0;
        for (k = 0; k < 3; k++) {
          for (j = 0; j < 3; j++) {
            if (ivert[k] == tria[3*i+j]) correct_vertex++;
          }
        }
        if (correct_vertex == 3)
          correct_vertices++;

        // check if shape functions are correct

        if (sqrt ((s[0]-1.0/3.0)*(s[0]-1.0/3.0)
                + (s[1]-1.0/3.0)*(s[1]-1.0/3.0)
                + (s[2]-1.0/3.0)*(s[2]-1.0/3.0)) <= tol)
          correct_shape_function++;
      }

      i++;
    }

    // output interpolation check results

    if (location_found == ntria)
      snprintf (Text, 512, "UVMAP    : All %d locations were found", ntria);
    else
      snprintf (Text, 512, "UVMAP    : Only %d of %d locations were found", location_found, ntria);
    uvmap_message (Text);

    if (correct_tria == ntria)
      snprintf (Text, 512, "UVMAP    : All %d locations have correct tria-face index", ntria);
    else
      snprintf (Text, 512, "UVMAP    : Only %d of %d locations have correct tria-face index", correct_tria, ntria);
    uvmap_message (Text);

    if (correct_vertices == ntria)
      snprintf (Text, 512, "UVMAP    : All %d locations have correct nodes/vertices", ntria);
    else
      snprintf (Text, 512, "UVMAP    : Only %d of %d locations have correct nodes/vertices", correct_vertices, ntria);
    uvmap_message (Text);

    if (correct_shape_function == ntria)
      snprintf (Text, 512, "UVMAP    : All %d locations have correct shape functions", ntria);
    else
      snprintf (Text, 512, "UVMAP    : Only %d of %d locations have correct shape functions", correct_shape_function, ntria);
    uvmap_message (Text);
  }

  // write uv surface mesh

  if (status == 0)
    status = EG_uvmap_Write (Case_Name, ntria, nvert, local_idef, tria, uv, NULL);

  // free UV mapping data structure

  uvmap_struct_free (ptr);

  // free starting data arrays

  uvmap_free (local_idef);
  uvmap_free (tria);
  uvmap_free (uv);
  uvmap_free (xyz);

  return status;
}
