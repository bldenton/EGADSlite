#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: EG_uvmap_Write.c,v 1.6 2021/04/30 18:58:36 marcum Exp $
 * Copyright 1994-2021, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Write data to a SURF type surface mesh file.
EGADS style data are used in this API.
--------------------------------------------------------------------------------

int EG_uvmap_Write (
  char Case_Name[],
  int ntria,
  int nvert,
  int *idmap,
  int *tria,
  double *uv,
  double *xyz);


INPUT ARGUMENTS
---------------

Case_Name	Case name for SURF type surface mesh file, Case_Name_uv.surf or
		Case_Name_xyz.surf.

ntria		Number of tria-faces.

nvert		Number of nodes/vertices.

idmap		Surface ID label for each tria-face (ntria in length).

tria		Tria-face connectivity (3*ntria in length).

uv		UV coordinates (2*nvert in length).
		If uv is allocated then a Case_Name_uv.surf file will be
		written with UV coordinates for the XY coordinates within the
		output file and zero for the Z.

xyz		XYZ coordinates (3*nvert in length).
		If xyz is allocated then a Case_Name_xyz.surf file will be
		written with XYZ coordinates.


RETURN VALUE
------------

EGADS_SUCCESS	Normal completion without errors.
EGADS_MALLOC	Unable to allocate required memory.
EGADS_WRITERR	An error occurred.


OUTPUT ARGUMENTS
----------------

*/

int EG_uvmap_Write (
  char Case_Name[],
  int ntria,
  int nvert,
  int *idmap,
  int *tria,
  double *uv,
  double *xyz)
{
  // Write data to a SURF type surface mesh file.

  INT_ *idibf = NULL;
  INT_3D *inibf = NULL;

  DOUBLE_2D *u = NULL;
  DOUBLE_3D *x = NULL;

  INT_ nbface = 0;
  INT_ nnode = 0;
  INT_ status = 0;

  status = uvmap_from_egads (ntria, nvert, idmap, tria, uv, xyz,
                             &nbface, &nnode, &idibf, &inibf, &u, &x);

  if (status == 0)
    status = uvmap_write (Case_Name, nbface, nnode, idibf, inibf, u, x);

  uvmap_free (idibf);
  uvmap_free (inibf);
  uvmap_free (u);
  uvmap_free (x);

  // set return value

  if (status > 100000)
    status = EGADS_MALLOC;
  else if (status)
    status = EGADS_WRITERR;
  else
    status = EGADS_SUCCESS;

  return (int) status;
}
