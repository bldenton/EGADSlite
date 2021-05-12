#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_struct_tasks.c,v 1.2 2020/06/15 18:08:07 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_struct_add_entry (
  INT_ idef,
  INT_ nbface,
  INT_ *idibf,
  INT_3D *inibf,
  INT_3D *ibfibf,
  DOUBLE_2D *u,
  uvmap_struct **uvmap_struct_ptr)
{
  // Add mapping structure entry and set data for surface idef.

  uvmap_struct *struct_ptr;

  INT_ *msrch = NULL;

  INT_ err = 0;
  INT_ found = 0;
  INT_ index = 0;
  INT_ jndex = 0;
  INT_ ndef = 0;

  msrch = uvmap_malloc (&err, (nbface+1)*sizeof (INT_));

  if (err) {
    uvmap_error_message ("*** ERROR 103501 : unable to allocate required memory ***");
    return 103501;
  }

  if (*uvmap_struct_ptr) {

    struct_ptr = &((*uvmap_struct_ptr)[0]);

    ndef = struct_ptr->ndef;

    do {

      struct_ptr = &((*uvmap_struct_ptr)[index]);

      if (struct_ptr->mdef == 0)
        found = 1;

      index++;
    }
    while (index < ndef && found == 0);

    if (found)
      index--;
    else
      index = ndef;
  }

  if (index == ndef) {

    ndef++;

    *uvmap_struct_ptr = (uvmap_struct *) uvmap_realloc (&err, *uvmap_struct_ptr, ndef*sizeof (uvmap_struct));

    if (err) {
      uvmap_free (msrch);
      uvmap_error_message ("*** ERROR 103502 : unable to allocate required memory ***");
      return 103502;
    }

    for (jndex = 0; jndex < ndef; jndex++) {
      struct_ptr = &((*uvmap_struct_ptr)[jndex]);
      struct_ptr->ndef = ndef;
    }
  }

  struct_ptr = &((*uvmap_struct_ptr)[index]);

  struct_ptr->mdef = 1;
  struct_ptr->idef = idef;
  struct_ptr->isrch = 0;
  struct_ptr->ibface = 1;
  struct_ptr->nbface = nbface;

  struct_ptr->idibf = idibf;
  struct_ptr->msrch = msrch;
  struct_ptr->inibf = inibf;
  struct_ptr->ibfibf = ibfibf;
  struct_ptr->u = u;

  uvmap_struct_set_srch_data (index, 0, 1, *uvmap_struct_ptr);

  return 0;
}

INT_ uvmap_struct_get_entry (
  INT_ idef,
  INT_ *index,
  INT_ *isrch,
  INT_ *ibface,
  INT_ *nbface,
  INT_ **idibf,
  INT_ **msrch,
  INT_3D **inibf,
  INT_3D **ibfibf,
  DOUBLE_2D **u,
  uvmap_struct *uvmap_struct_ptr)
{
  // Get mapping data structure data for surface idef.

  uvmap_struct *struct_ptr;

  INT_ status = 0;

  status = uvmap_struct_find_entry (idef, index, uvmap_struct_ptr);

  if (status == 0 && *index == -1) {
    uvmap_error_message ("*** ERROR 3501 : unable to find mapping surface ***");
    status = 3501;
  }

  if (status)
    return status;

  struct_ptr = &(uvmap_struct_ptr[*index]);

  *isrch = struct_ptr->isrch;
  *ibface = struct_ptr->ibface;
  *nbface = struct_ptr->nbface;

  *idibf = struct_ptr->idibf;
  *msrch = struct_ptr->msrch;
  *inibf = struct_ptr->inibf;
  *ibfibf = struct_ptr->ibfibf;
  *u = struct_ptr->u;

  return 0;
}

void uvmap_struct_set_srch_data (
  INT_ index,
  INT_ isrch,
  INT_ ibface,
  uvmap_struct *uvmap_struct_ptr)
{
  // Set searching data.

  uvmap_struct *struct_ptr;

  INT_ i;
  INT_ nsrch = 1000000;

  if (uvmap_struct_ptr == NULL)
    return;

  struct_ptr = &(uvmap_struct_ptr[index]);

  if (isrch <= 0 || isrch > nsrch) {

    isrch = 1;

    for (i = 1; i <= struct_ptr->nbface; i++) {
      struct_ptr->msrch[i] = 0;
    }
  }

  struct_ptr->isrch = isrch;
  struct_ptr->ibface = ibface;

  return;
}

INT_ uvmap_struct_find_entry (
  INT_ idef,
  INT_ *index,
  uvmap_struct *uvmap_struct_ptr)
{
  // Find mapping structure index for surface idef.

  uvmap_struct *struct_ptr;

  INT_ ndef;

  if (uvmap_struct_ptr == NULL) {
    uvmap_error_message ("*** ERROR 3502 : mapping surface structure not set ***");
    return 3502;
  }

  struct_ptr = &(uvmap_struct_ptr[0]);

  ndef = struct_ptr->ndef;

  *index = 0;

  do {

    struct_ptr = &(uvmap_struct_ptr[*index]);

    (*index)++;
  }
  while (*index < ndef && (struct_ptr->idef != idef || struct_ptr->mdef == 0));

  (*index)--;

  if (struct_ptr->idef != idef || struct_ptr->mdef == 0)
    *index = -1;

  return 0;
}

void uvmap_struct_free (void *ptr)
{
  // Free UV mapping data structure for all surfaces.

  uvmap_struct *struct_ptr;
  uvmap_struct *uvmap_struct_ptr;

  INT_ index, ndef;

  if (ptr == NULL)
    return;

  uvmap_struct_ptr = (uvmap_struct *) ptr;

  struct_ptr = &(uvmap_struct_ptr[0]);

  ndef = struct_ptr->ndef;

  for (index = 0; index < ndef; index++) {
    uvmap_struct_free_index (index, uvmap_struct_ptr);
  }

  uvmap_free (uvmap_struct_ptr);

  return;
}

void uvmap_struct_free_idef (INT_ idef, uvmap_struct *uvmap_struct_ptr)
{
  // Free mapping data structure for surface idef.

  INT_ index = 0;

  if (uvmap_struct_ptr == NULL)
    return;

  uvmap_struct_find_entry (idef, &index, uvmap_struct_ptr);

  if (index >= 0)
    uvmap_struct_free_index (index, uvmap_struct_ptr);

  return;
}

void uvmap_struct_free_index (INT_ index, uvmap_struct *uvmap_struct_ptr)
{
  // Free mapping data structure for surface at location index.

  uvmap_struct *struct_ptr;

  struct_ptr = &(uvmap_struct_ptr[0]);

  if (index < 0 || index >= struct_ptr->ndef)
    return;

  struct_ptr = &(uvmap_struct_ptr[index]);

  if (struct_ptr->mdef == 0)
    return;

  uvmap_free (struct_ptr->idibf);
  uvmap_free (struct_ptr->msrch);
  uvmap_free (struct_ptr->inibf);
  uvmap_free (struct_ptr->ibfibf);
  uvmap_free (struct_ptr->u);

  struct_ptr->idef = 0;
  struct_ptr->mdef = 0;
  struct_ptr->isrch = 0;
  struct_ptr->ibface = 0;
  struct_ptr->nbface = 0;

  struct_ptr->idibf = NULL;
  struct_ptr->msrch = NULL;
  struct_ptr->inibf = NULL;
  struct_ptr->ibfibf = NULL;
  struct_ptr->u = NULL;

  return;
}
