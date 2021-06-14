#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_read.c,v 1.29 2021/04/27 18:34:23 marcum Exp $
 * Copyright 1994-2021, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Read and allocate data from a SURF type surface mesh file.
--------------------------------------------------------------------------------

INT_ uvmap_read (
  char Case_Name[],
  INT_ *nbface,
  INT_ *nnode,
  INT_ **idibf,
  INT_3D **inibf,
  DOUBLE_3D **x);


INPUT ARGUMENTS
---------------

Case_Name	Case name for SURF type surface mesh file, Case_Name.surf.


RETURN VALUE
------------

0		Normal completion without errors.
>0		An error occurred.


OUTPUT ARGUMENTS
----------------

nbface		Number of tria-faces.

nnode		Number of nodes/vertices.

idibf		Surface ID label for each tria-face (nbface+1 in length).

inibf		Tria-face connectivity (nbface+1 in length).

x		XYZ coordinates (nnode+1 in length).

*/

INT_ uvmap_read (
  char Case_Name[],
  INT_ *nbface,
  INT_ *nnode,
  INT_ **idibf,
  INT_3D **inibf,
  DOUBLE_3D **x)
{
  // Read and allocate data from a SURF type surface mesh file.

  FILE *Grid_File;

  char File_Name[500], Text_Line[512], Text[512];
  char *Read_Label;

  int i1, i2, i3, i4, i5, i6;
  INT_ i, n;
  INT_ status = 0;

  double d1, d2;

  *nbface = 0;
  *nnode = 0;
  *idibf = NULL;
  *inibf = NULL;
  *x = NULL;

  strcpy (File_Name, Case_Name);
  strcat (File_Name, ".surf");

  uvmap_message ("");
  snprintf (Text, 512, "UVMAP    : Reading Surf File = %s", File_Name);
  uvmap_message (Text);

  Grid_File = fopen (File_Name, "r");

  if (Grid_File == NULL)
  {
    uvmap_error_message ("*** ERROR 3509 : error opening surface mesh file ***");
    return 3509;
  }

  n = fscanf (Grid_File, "%d %d %d", &i1, &i2, &i3);

  *nbface = (INT_) i1;
  *nnode  = (INT_) i3;

  if (n == EOF) {
    uvmap_error_message ("*** ERROR 3510 : error reading SURF surface mesh file ***");
    status = 3510;
  }

  if (status == 0) {
    *idibf = (INT_ *) uvmap_malloc (&status, ((*nbface)+1)*sizeof(INT_));
    *inibf = (INT_3D *) uvmap_malloc (&status, ((*nbface)+1)*sizeof(INT_3D));
    *x = (DOUBLE_3D *) uvmap_malloc (&status, ((*nnode)+1)*sizeof(DOUBLE_3D));

    if (status) {
      uvmap_error_message ("*** ERROR 103515 : unable to allocate required memory ***");
      status = 103515;
    }
  }

  if (status == 0) {
#ifndef __clang_analyzer__
    Read_Label = fgets (Text_Line, 512, Grid_File);
#endif
    Read_Label = fgets (Text_Line, 512, Grid_File);

    if (Read_Label == NULL) {
      uvmap_error_message ("*** ERROR 3511 : error reading SURF surface mesh file ***");
      status = 3511;
    }
  }

  if (status == 0) {
    i = 1;
    n = sscanf (Text_Line, "%lf %lf %lf %lf %lf",
                &((*x)[i][0]), &((*x)[i][1]), &((*x)[i][2]), &d1, &d2);

    if (n == 3) {
      for (i = 2; i <= *nnode; i++) {
        n = fscanf (Grid_File, "%lf %lf %lf",
                    &((*x)[i][0]), &((*x)[i][1]), &((*x)[i][2]));
      }
    }

    else if (n == 4) {
      for (i = 2; i <= *nnode; i++) {
        n = fscanf (Grid_File, "%lf %lf %lf %lf",
                    &((*x)[i][0]), &((*x)[i][1]), &((*x)[i][2]), &d1);
      }
    }

    else if (n == 5) {
      for (i = 2; i <= *nnode; i++) {
        n = fscanf (Grid_File, "%lf %lf %lf %lf %lf",
                    &((*x)[i][0]), &((*x)[i][1]), &((*x)[i][2]), &d1, &d2);
      }
    }

    else
      n = EOF;

    if (n == EOF) {
      uvmap_error_message ("*** ERROR 3512 : error reading SURF surface mesh file ***");
      status = 3512;
    }
  }

  if (status == 0) {
    for (i = 1; i <= *nbface; i++) {
      n = fscanf (Grid_File, "%d %d %d %d %d %d", &i1, &i2, &i3, &i4, &i5, &i6);

      (*inibf)[i][0] = (INT_) i1;
      (*inibf)[i][1] = (INT_) i2;
      (*inibf)[i][2] = (INT_) i3;
      (*idibf)[i]    = (INT_) i4;
    }

    if (n == EOF) {
      uvmap_error_message ("*** ERROR 3513 : error reading SURF surface mesh file ***");
      status = 3513;
    }
  }

  fclose (Grid_File);

  if (status == 0) {
    uvmap_message ("");
    snprintf (Text, 512, "UVMAP    : Tria Surface Faces=%10d", (int) *nbface);
    uvmap_message (Text);
    snprintf (Text, 512, "UVMAP    : Nodes             =%10d", (int) *nnode);
    uvmap_message (Text);
    uvmap_message ("");
  }

  return 0;
}
