#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: EG_uvmapStructFree.c,v 1.1 2021/02/26 23:09:19 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

/*

--------------------------------------------------------------------------------
Free UV mapping data structure for all surfaces.
--------------------------------------------------------------------------------

void EG_uvmapStructFree (void *ptr);


INPUT ARGUMENTS
---------------

ptr	UV mapping data structure.

*/

void EG_uvmapStructFree (void *ptr)
{
  uvmap_struct_free (ptr);
  return;
}
