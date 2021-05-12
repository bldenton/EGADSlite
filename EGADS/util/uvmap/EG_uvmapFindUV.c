#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: EG_uvmapFindUV.c,v 1.1 2021/02/26 22:46:08 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Find location of given UV coordinates.
EGADS style data are used in this API.
--------------------------------------------------------------------------------

int EG_uvmapFindUV (
  int idef,
  double uv[2],
  void *ptr,
  int *local_idef,
  int *itria,
  int ivertex[3],
  double s[3]);


INPUT ARGUMENTS
---------------

idef		Surface ID label.

uv		UV coordinate location to find (2 in length).

ptr		UV mapping data structure.


RETURN VALUE
------------

EGADS_SUCCESS	UV coordinate location was found.
EGADS_NOTFOUND	UV coordinate location was not found.
EGADS_MALLOC	Unable to allocate required memory.
EGADS_UVMAP	An error occurred.


OUTPUT ARGUMENTS
----------------

local_idef	Local surface ID label of the tria-face that contains the given
		UV coordinates.
		If the surface idef is a virtual/composite surface then it is
		composed of one or more local surfaces with differing surface ID
		labels.
		If the surface idef is a standard surface then the local surface
		ID label found is set to idef. 

itria		Tria-face index on surface idef of the tria-face that contains
		the given UV coordinates.

ivertex		Node/Vertex of the tria-face that contains the given UV
		coordinates (3 in length).

s		Linear interpolation shape functions for the tria-face that
		contains the given UV coordinates (3 in length). For example,
		given data in array data stored at nodes/vertices of the surface
		mesh, the interpolated value at the location found can be
		determined from the following expression.

		  data_intp = s[0] * data[ivertex[0]]
		            + s[1] * data[ivertex[1]]
		            + s[2] * data[ivertex[2]]

*/

int EG_uvmapFindUV (
  int idef,
  double uv[2],
  void *ptr,
  int *local_idef,
  int *itria,
  int ivertex[3],
  double s[3])
{
  // Find location of given UV coordinates.

  INT_ inode[3];
  INT_ ibface, local_idef_;
  int status = 0;

  // find location of given UV coordinates

  status = (int) uvmap_find_uv ((INT_) idef, uv, ptr,
                                &local_idef_, &ibface, inode, s);

  *itria = (int) ibface;
  *local_idef = (int) local_idef_;

  ivertex[0] = (int) inode[0];
  ivertex[1] = (int) inode[1];
  ivertex[2] = (int) inode[2];

  // set return value

  if (status > 100000)
    status = EGADS_MALLOC;
  else if (status > 0)
    status = EGADS_UVMAP;
  else if (status == -1)
    status = EGADS_NOTFOUND;
  else
    status = EGADS_SUCCESS;

  return status;
}
