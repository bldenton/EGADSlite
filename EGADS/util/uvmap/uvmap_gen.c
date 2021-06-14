#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_gen.c,v 1.62 2021/02/28 22:11:00 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Generate a UV coordinate map for a given tria-face triangulation, create a UV
mapping data structure, and store a copy of UV mapping data within structure.
--------------------------------------------------------------------------------

INT_ uvmap_gen (
  INT_ idef,
  INT_ nbface,
  INT_ nnode,
  INT_ set_struct,
  INT_ verbosity,
  INT_ *idibf,
  INT_3D **inibf,
  DOUBLE_3D **x,
  DOUBLE_2D **u,
  void **ptr);


INPUT ARGUMENTS
---------------

idef		Surface ID label.
		Not used if set_struct=0.

nbface		Number of tria-faces.

nnode		Number of nodes/vertices.

set_struct	UV mapping data structure flag.
		If set_struct=1, then create UV mapping data structure.
		If set_struct=0, then do not create UV mapping data structure.

verbosity	Message flag.
		If verbosity=0, then do not output progress information.
		If verbosity=1, then output progress information to stdout.
		If verbosity=2, then output progress and additional CPU usage
		information to stdout.

idibf		Local surface ID label for each tria-face of surface idef
		(nbface+1 in length).
		If the surface idef is a virtual/composite surface, then it is
		composed of one or more local surface with differing surface ID
		labels.
		If the surface idef is a standard surface, then the local
		surface ID label need not be set as it is the same as that for
		surface idef. 
		Not used if set_struct=0.

inibf		Tria-face connectivity (nbface+1 in length).

x		XYZ coordinates (nnode+1 in length).
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

0		Normal completion without errors.
>0		An error occurred.


OUTPUT ARGUMENTS
----------------

inibf		Tria-face connectivity (nbface+1 in length).
		Connectivity may be reordered for ordering consistency.
		The address may change as input arrays are temporarily
		reallocated to fill holes if there are inner closed curves.

x		XYZ coordinates (nnode+1 in length).
		The address may change as input arrays are temporarily
		reallocated to fill holes if there are inner closed curves.

u		Generated UV coordinates (nnode+1 in length).

ptr		UV mapping data structure.
		If the structure ptr is NULL on input, then the structure is
		allocated and an entry for surface idef is added.
		If the structure ptr is not NULL on input (from a previous call
		to this routine), then the structure is reallocated and an entry
		for surface idef is added.
		Not used if set_struct=0.
		Note that a copy of the local surface ID label idibf,
		connectivity inibf, and UV coordinates u are saved within the
		UV mapping data structure.

*/

/*

GLOBAL VARIABLES
----------------

WriteUVMapOut	Output file flag.
		If set to a non zero value then an output file is written for
		the given XYZ surface mesh and generated UV surface mesh.

UVmapCaseName	Output file case name.
		If set to anything other than "_null_" then the output file
		names will be UVmapCaseName.uvmap_#_uv.surf and
		UVmapCaseName.uvmap_#_xyz.surf, where # is the Surface ID label
		idef.
		If set to "_null_" then the output file names will be
		UVmapOut_#_uv.surf and UVmapOut_#_xyz.surf, where # is 1,2,3...
		and # is incremented each time uvmap_gen is called.

*/

int WriteUVmapOut = 0;

char UVmapCaseName[512] = "_null_";

INT_ uvmap_gen (
  INT_ idef,
  INT_ nbface,
  INT_ nnode,
  INT_ set_struct,
  INT_ verbosity,
  INT_ *idibf,
  INT_3D **inibf,
  DOUBLE_3D **x,
  DOUBLE_2D **u,
  void **ptr)
{
  // Generate a UV coordinate map for a given tria-face triangulation, create a
  // UV mapping data structure, and store a copy of UV mapping data within
  // structure.

  char caseName[512];

  INT_   *idibf_ = NULL;
  INT_3D *ibfibf = NULL;
  INT_3D *inibf_ = NULL;
  INT_2D *inibe = NULL;

  DOUBLE_2D *u_ = NULL;

  INT_ ibface, inode, nbedge;
  INT_ err = 0;
  INT_ status = 0;

  static int outFileCnt = 0;

  // generate a UV surface mesh given an XYZ surface mesh

  if (status == 0)
    status = uvmap_gen_uv (verbosity, &nbedge, &nbface, &nnode,
                           &ibfibf, &inibe, inibf, x, &u_);

  // allocate a copy of the surface ID and tria-face connectivity for the UV
  // mapping structure and generated UV coordinates for the output arguments

  if (status == 0) {
    if (idibf) idibf_ = (INT_ *) uvmap_malloc (&err, (nbface+1)*sizeof(INT_));
    inibf_ = (INT_3D *) uvmap_malloc (&err, (nbface+1)*sizeof(INT_3D));
    *u = (DOUBLE_2D *) uvmap_malloc (&err, (nnode+1)*sizeof(DOUBLE_2D));

    if (err) {
      uvmap_error_message ("*** ERROR 103504 : unable to allocate required memory ***");
      status = 103504;
    }
  }

  // copy surface ID and tria-face connectivity for the UV mapping structure

  if (status == 0) {
    if (idibf) {
      for (ibface = 1; ibface <= nbface; ibface++) {
        idibf_[ibface] = idibf[ibface];
      }
    }

    for (ibface = 1; ibface <= nbface; ibface++) {
      inibf_[ibface][0] = (*inibf)[ibface][0];
      inibf_[ibface][1] = (*inibf)[ibface][1];
      inibf_[ibface][2] = (*inibf)[ibface][2];
    }
  }

  // copy generated UV coordinates for the output arguments

  if (status == 0) {
    for (inode = 1; inode <= nnode; inode++) {
      (*u)[inode][0] = u_[inode][0];
      (*u)[inode][1] = u_[inode][1];
    }
  }

  // create and/or add entry for surface idef in UV mapping data structure

  if (status == 0 && set_struct)
    status = uvmap_struct_add_entry (idef, nbface, idibf_, inibf_, ibfibf, u_,
                                     (uvmap_struct **) ptr);

  // free temporary data arrays

  uvmap_free (inibe);

  // write input xyz and output uv coordinates to a SURF file

  if (status == 0 && WriteUVmapOut) {

    if (strcmp (UVmapCaseName, "_null_"))
      snprintf (caseName, 512, "%s.uvmap_%d", UVmapCaseName, idef);
    else {
      outFileCnt++;
      snprintf (caseName, 512, "UVmapOut_%d", outFileCnt);
    }

    status = (int) uvmap_write (caseName, nbface, nnode, idibf_, inibf_, u_, *x);
  }

  return status;
}
