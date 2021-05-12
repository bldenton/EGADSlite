#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: EG_uvmap_Read.c,v 1.6 2021/04/30 18:58:36 marcum Exp $
 * Copyright 1994-2021, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Read and allocate data from a SURF type surface mesh file.
EGADS style data are used in this API.
--------------------------------------------------------------------------------

int EG_uvmap_Read (
  char Case_Name[],
  int *ntria,
  int *nvert,
  int **idmap,
  int **tria,
  double **xyz);


INPUT ARGUMENTS
---------------

Case_Name	Case name for SURF type surface mesh file, Case_Name.surf.


RETURN VALUE
------------

EGADS_SUCCESS	Normal completion without errors.
EGADS_MALLOC	Unable to allocate required memory.
EGADS_READERR	An error occurred.


OUTPUT ARGUMENTS
----------------

ntria		Number of tria-faces.

nvert		Number of nodes/vertices.

idmap		Surface ID label for each tria-face (ntria in length).

tria		Tria-face connectivity (3*ntria in length).

xyz		XYZ coordinates (3*nvert in length).

*/

int EG_uvmap_Read (
  char Case_Name[],
  int *ntria,
  int *nvert,
  int **idmap,
  int **tria,
  double **xyz)
{
  // Read and allocate data from a SURF type surface mesh file.

  INT_ *idibf = NULL;
  INT_3D *inibf = NULL;

  double *uv = NULL;
  DOUBLE_3D *x = NULL;

  INT_ nbface = 0;
  INT_ nnode = 0;
  INT_ status = 0;

  *ntria = 0;
  *nvert = 0;
  *idmap = NULL;
  *tria = NULL;
  *xyz = NULL;

  status = uvmap_read (Case_Name, &nbface, &nnode, &idibf, &inibf, &x);

  if (status == 0)
    status = uvmap_to_egads (nbface, nnode, idibf, inibf, NULL, x,
                             ntria, nvert, idmap, tria, &uv, xyz);

  uvmap_free (idibf);
  uvmap_free (inibf);
  uvmap_free (x);

  // set return value

  if (status > 100000)
    status = EGADS_MALLOC;
  else if (status)
    status = EGADS_READERR;
  else
    status = EGADS_SUCCESS;

  return (int) status;
}
