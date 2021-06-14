#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_message.c,v 1.5 2021/02/28 22:04:56 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

FILE *UVMAP_Output_File = NULL;

// Write a message to standard output and a file if specified.

void uvmap_message (char *Text) {
  fprintf (stdout, "%s\n", Text);
  fflush (stdout);
  if (UVMAP_Output_File) {
    fprintf (UVMAP_Output_File, "%s\n", Text);
    fflush (UVMAP_Output_File);
  }
  return;
}

// Write a message to standard error and a file if specified.

void uvmap_error_message (char *Text) {
  fprintf (stderr, "%s\n", Text);
  fflush (stderr);
  if (UVMAP_Output_File) {
    fprintf (UVMAP_Output_File, "%s\n", Text);
    fflush (UVMAP_Output_File);
  }
  return;
}
