#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: EG_uvmapGen.c,v 1.1 2021/02/26 22:46:08 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Generate a UV coordinate map for a given tria-face triangulation, create a UV
mapping data structure, and store a copy of UV mapping data within structure.
EGADS style data are used in this API.
--------------------------------------------------------------------------------

int EG_uvmapGen (
  int idef,
  int ntria,
  int nvert,
  int set_struct,
  int verbosity,
  int *local_idef,
  int *tria,
  double *xyz,
  double **uv,
  void **ptr);


INPUT ARGUMENTS
---------------

idef		Surface ID label.
		Not used if set_struct=0.

ntria		Number of tria-faces.

nvert		Number of nodes/vertices.

set_struct	UV mapping data structure flag.
		If set_struct=1, then create UV mapping data structure.
		If set_struct=0, then do not create UV mapping data structure.

verbosity	Message flag.
		If verbosity=0, then do not output progress information.
		If verbosity=1, then output progress information to stdout.
		If verbosity=2, then output progress and additional CPU usage
		information to stdout.

local_idef	Local surface ID label for each tria-face of surface idef
		(ntria in length).
		If the surface idef is a virtual/composite surface, then it is
		composed of one or more local surface with differing surface ID
		labels.
		If the surface idef is a standard surface, then the local
		surface ID label need not be set as it is the same as that for
		surface idef. 
		Not used if set_struct=0.

tria		Tria-face connectivity (3*ntria in length).

xyz		XYZ coordinates (3*nvert in length).
		The XYZ coordinates are used only to determine discontinuous
		locations on the outer and inner (if any) boundary edges.
		If the XYZ coordinates are NULL on input, then discontinuity on
		the outer and inner (if any) boundary edges is not considered.

ptr		UV mapping data structure.
		Set ptr to NULL if this is the first call to this routine and
		use previously set ptr on subsequent calls.
		Not used if set_struct=0.


RETURN VALUE
------------

EGADS_SUCCESS	Normal completion without errors.
EGADS_MALLOC	Unable to allocate required memory.
EGADS_UVMAP	An error occurred.


OUTPUT ARGUMENTS
----------------

uv		Generated UV coordinates (2*nvert in length).

ptr		UV mapping data structure.
		If the structure ptr is NULL on input, then the structure is
		allocated and an entry for surface idef is added.
		If the structure ptr is not NULL on input (from a previous call
		to this routine), then the structure is reallocated and an entry
		for surface idef is added.
		Not used if set_struct=0.
		Note that a copy of the local surface ID label local_idef,
		connectivity tria, and UV coordinates uv are saved within the
		UV mapping data structure.
		Also, note that the tria-face connectivity that is stored within
		the UV mapping structure may have been reordered for ordering
		consistency and therefore may differ from that of the input
		connectivity, tria.

*/

int EG_uvmapGen (
  int idef,
  int ntria,
  int nvert,
  int set_struct,
  int verbosity,
  int *local_idef,
  int *tria,
  double *xyz,
  double **uv,
  void **ptr)
{
  // Generate a UV coordinate map for a given tria-face triangulation, create a
  // UV mapping data structure, and store a copy of UV mapping data within
  // structure.

  INT_ *idibf = NULL;
  INT_3D *inibf = NULL;

  DOUBLE_2D *u = NULL;
  DOUBLE_3D *x = NULL;

  INT_ i, nbface, nnode;
  INT_ err = 0;
  INT_ status = 0;

  // convert from EGADS to AFLR style arrays

  status = uvmap_from_egads (ntria, nvert, local_idef, tria, *uv, xyz,
                             &nbface, &nnode, &idibf, &inibf, &u, &x);

  // generate a UV coordinate map for a given tria-face triangulation, create a
  // UV mapping data structure, and store a copy of UV mapping data within
  // structure

  if (status == 0)
    status = (int) uvmap_gen ((INT_) idef, nbface, nnode,
                              (INT_) set_struct, (INT_) verbosity,
                              idibf, &inibf, &x, &u, ptr);

  // allocate generated UV coordinates

  if (status == 0) {
    *uv = (double *) uvmap_malloc (&err, 2*nvert*sizeof(double));

    if (err) {
      uvmap_error_message ("*** ERROR 103504 : unable to allocate required memory ***");
      status = 103504;
    }
  }

  // copy generated UV coordinates for the output arguments

  if (status == 0) {
    for (i = 0; i < nnode; i++) {
      (*uv)[2*i]   = u[i+1][0];
      (*uv)[2*i+1] = u[i+1][1];
    }
  }

  // free temporary data arrays

  uvmap_free (idibf);
  uvmap_free (inibf);
  uvmap_free (u);
  uvmap_free (x);

  if (status) {
    uvmap_free (*uv);
    *uv = NULL;
  }

  // set return value

  if (status > 100000)
    status = EGADS_MALLOC;
  else if (status)
    status = EGADS_UVMAP;
  else
    status = EGADS_SUCCESS;

  return (int) status;
}
