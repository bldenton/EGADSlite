#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_malloc.c,v 1.2 2020/06/04 02:32:57 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

void (*ext_uvmap_free_) (void *ptr) = NULL;
void * (*ext_uvmap_malloc_) (INT_ *err_flag, size_t size) = NULL;
void * (*ext_uvmap_realloc_) (INT_ *err_flag, void *ptr, size_t size) = NULL;

void uvmap_register_ext_free (void (*ext_free_routine) (void *ptr)) {
  ext_uvmap_free_ = ext_free_routine;
  return;
}

void uvmap_register_ext_malloc (void * (*ext_malloc_routine) (INT_ *err_flag, size_t size)) {
  ext_uvmap_malloc_ = ext_malloc_routine;
  return;
}

void uvmap_register_ext_realloc (void * (*ext_realloc_routine) (INT_ *err_flag, void *ptr, size_t size)) {
  ext_uvmap_realloc_ = ext_realloc_routine;
  return;
}

void uvmap_free (void *ptr) {
  if (ptr) {
    if (ext_uvmap_free_)
      ext_uvmap_free_ (ptr);
    else
      free (ptr);
  }
  return;
}

void * uvmap_malloc (INT_ *err_flag, size_t size) {
  void *pointer = NULL;
  if (size) {
    if (ext_uvmap_malloc_)
      pointer = ext_uvmap_malloc_ (err_flag, size);
    else {
      pointer = malloc (size);
      if (pointer == NULL)
        (*err_flag)++;
    }
  }
  return pointer;
}

void * uvmap_realloc (INT_ *err_flag, void *ptr, size_t size) {
  void *pointer = NULL;
  if (size) {
    if (ext_uvmap_realloc_)
      pointer = ext_uvmap_realloc_ (err_flag, ptr, size);
    else {
      pointer = realloc (ptr, size);
      if (pointer == NULL)
        (*err_flag)++;
    }
  }
  return pointer;
}
