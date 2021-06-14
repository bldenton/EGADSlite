#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_test.c,v 1.11 2021/04/27 18:34:23 marcum Exp $
 * Copyright 1994-2021, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Test uvmap.
--------------------------------------------------------------------------------

INT_ uvmap_test (char Case_Name[], INT_ verbosity);


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

0		Normal completion without errors.
>0		An error occurred.

*/

INT_ uvmap_test (char Case_Name[], INT_ verbosity)
{
  // Test uvmap with EGADS APIs.

  void *ptr = NULL;

  INT_ *idibf = NULL;
  INT_3D *inibf = NULL;
  DOUBLE_2D *u = NULL;
  DOUBLE_3D *x = NULL;

  INT_ idef = 1;
  INT_ nbface = 0;
  INT_ nnode = 0;
  INT_ set_struct = 1;
  INT_ status = 0;

  // interpolation testing variables

  char Text[512];

  INT_ inode_[3];
  INT_ ibface, ibface_, j, local_idefi;
  INT_ location_found = 0;
  INT_ correct_shape_function = 0;
  INT_ correct_tria = 0;
  INT_ correct_vertex = 0;
  INT_ correct_vertices = 0;

  double u_[2], s[3];
  double tol = 1.0e-6;

  // read xyz surface mesh

  status = uvmap_read (Case_Name, &nbface, &nnode, &idibf, &inibf, &x);

  // generate uv mapping

  if (status == 0)
    status = uvmap_gen (idef, nbface, nnode, set_struct, verbosity,
                        idibf, &inibf, &x, &u, &ptr);

  // check interpolation

  if (status == 0) {

    uvmap_message ("");
    uvmap_message ("UVMAP    : CHECKING INTERPOLATION");
    uvmap_message ("");

    // loop over tria-faces and check interpolation at tria-face centroids

    ibface = 1;

    while (ibface <= nbface && status <= 0) {

      // set UV coordinates of interpolation location equal to the centroid of
      // tria-face

      for (j = 0; j < 2; j++) {
        u_[j] = (u[inibf[ibface][0]][j]
               + u[inibf[ibface][1]][j]
               + u[inibf[ibface][2]][j]) / 3.0;
      }

      // find interpolation location

      status = uvmap_find_uv (idef, u_, ptr, &local_idefi, &ibface_, inode_, s); 
      // if location was found then check interpolation results

      if (status == EGADS_SUCCESS) {

        location_found++;

        // check if tria-face index is correct

        if (ibface_ == ibface)
          correct_tria++;

        // check if nodes/vertices are correct
        // note that all orientations are checked as ordering may have changed

        correct_vertex = 0;
        for (j = 0; j < 3; j++) {
          if (inode_[j] == inibf[ibface][j]) correct_vertex++;
        }
        if (correct_vertex == 3)
          correct_vertices++;

        // check if shape functions are correct

        if (sqrt ((s[0]-1.0/3.0)*(s[0]-1.0/3.0)
                + (s[1]-1.0/3.0)*(s[1]-1.0/3.0)
                + (s[2]-1.0/3.0)*(s[2]-1.0/3.0)) <= tol)
          correct_shape_function++;
      }

      ibface++;
    }

    // output interpolation check results

    if (location_found == nbface)
      snprintf (Text, 512, "UVMAP    : All %d locations were found", (int) nbface);
    else
      snprintf (Text, 512, "UVMAP    : Only %d of %d locations were found", (int) location_found, (int) nbface);
    uvmap_message (Text);

    if (correct_tria == nbface)
      snprintf (Text, 512, "UVMAP    : All %d locations have correct tria-face index", (int) nbface);
    else
      snprintf (Text, 512, "UVMAP    : Only %d of %d locations have correct tria-face index", (int) correct_tria, (int) nbface);
    uvmap_message (Text);

    if (correct_vertices == nbface)
      snprintf (Text, 512, "UVMAP    : All %d locations have correct nodes/vertices", (int) nbface);
    else
      snprintf (Text, 512, "UVMAP    : Only %d of %d locations have correct nodes/vertices", (int) correct_vertices, (int) nbface);
    uvmap_message (Text);

    if (correct_shape_function == nbface)
      snprintf (Text, 512, "UVMAP    : All %d locations have correct shape functions", (int) nbface);
    else
      snprintf (Text, 512, "UVMAP    : Only %d of %d locations have correct shape functions", (int) correct_shape_function, (int) nbface);
    uvmap_message (Text);
  }

  // write uv surface mesh

  if (status == 0)
    status = uvmap_write (Case_Name, nbface, nnode, idibf, inibf, u, NULL);

  // free UV mapping data structure

  uvmap_struct_free (ptr);

  // free starting data arrays

  uvmap_free (idibf);
  uvmap_free (inibf);
  uvmap_free (u);
  uvmap_free (x);

  return status;
}
