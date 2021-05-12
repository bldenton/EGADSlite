#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_write.c,v 1.26 2021/02/25 23:27:27 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */
/*

--------------------------------------------------------------------------------
Write data to a SURF type surface mesh file.
--------------------------------------------------------------------------------

INT_ uvmap_write (
  char Case_Name[],
  INT_ nbface,
  INT_ nnode,
  INT_ *idibf,
  INT_3D *inibf,
  DOUBLE_2D *u,
  DOUBLE_3D *x);


INPUT ARGUMENTS
---------------

Case_Name	Case name for SURF type surface mesh file, Case_Name_uv.surf or
		Case_Name_xyz.surf.

nbface		Number of tria-faces.

nnode		Number of nodes/vertices.

idibf		Surface ID label for each tria-face (nbface+1 in length).

inibf		Tria-face connectivity (nbface+1 in length).

u		UV coordinates (nnode+1 in length).
		If u is allocated then a Case_Name_uv.surf file will be
		written with UV coordinates for the XY coordinates within the
		output file and zero for the Z.

x		XYZ coordinates (nnode+1 in length).
		If x is allocated then a Case_Name_xyz.surf file will be
		written with XYZ coordinates.


RETURN VALUE
------------

0		Normal completion without errors.
>0		An error occurred.


*/

INT_ uvmap_write (
  char Case_Name[],
  INT_ nbface,
  INT_ nnode,
  INT_ *idibf,
  INT_3D *inibf,
  DOUBLE_2D *u,
  DOUBLE_3D *x)
{
  // Write surface mesh file.

  FILE *Grid_File;

  char File_Name[512], Text[512];

  INT_ i, n, type;

  for (type = 1; type <= 2; type++) {

    if ((type == 1 && u) || (type == 2 && x)) {

      strcpy (File_Name, Case_Name);
      if (type == 1)
        strcat (File_Name, "_uv.surf");
      else
        strcat (File_Name, "_xyz.surf");

      uvmap_message ("");
      snprintf (Text, 512, "UVMAP    : Writing Surf File = %s", File_Name);
      uvmap_message (Text);
      uvmap_message ("");
      snprintf (Text, 512, "UVMAP    : Tria Surface Faces=%10d", (int) nbface);
      uvmap_message (Text);
      snprintf (Text, 512, "UVMAP    : Nodes             =%10d", (int) nnode);
      uvmap_message (Text);
      uvmap_message ("");

      Grid_File = fopen (File_Name, "w");

      if (Grid_File == NULL)
      {
        uvmap_error_message ("*** ERROR 3514 : error opening surface mesh file ***");
        return 3514;
      }

      n = fprintf (Grid_File, "%d 0 %d\n", (int) nbface, (int) nnode);

      if (type == 1) {
        for (i = 1; i <= nnode; i++) {
          n = fprintf (Grid_File, "%lf %lf 0.0\n", u[i][0], u[i][1]);
        }
      }
      else {
        for (i = 1; i <= nnode; i++) {
          n = fprintf (Grid_File, "%lf %lf %lf\n", x[i][0], x[i][1], x[i][2]);
        }
      }

      if (idibf) {
        for (i = 1; i <= nbface; i++) {
          n = fprintf (Grid_File, "%d %d %d %d 0 0\n", (int) inibf[i][0], (int) inibf[i][1], (int) inibf[i][2], (int) idibf[i]);
        }
      }
      else {
        for (i = 1; i <= nbface; i++) {
          n = fprintf (Grid_File, "%d %d %d 1 0 0\n", (int) inibf[i][0], (int) inibf[i][1], (int) inibf[i][2]);
        }
      }

      if (n < 0) {
        uvmap_error_message ("*** ERROR 3515 : error writing SURF surface mesh file ***");
        fclose (Grid_File);
        return 3515;
      }

      fclose (Grid_File);
    }
  }

  return 0;
}
