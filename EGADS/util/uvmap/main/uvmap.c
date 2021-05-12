#include "../UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap.c,v 1.37 2021/04/27 22:15:57 marcum Exp $
 * Copyright 1994-2021, David L. Marcum
 */

int main (int argc, char *argv[]) {

  extern FILE *UVMAP_Output_File;

  extern char UVmapCaseName[];

  extern int WriteUVmapOut;

  char Case_Name[512] = "_null_";
  char File_Name[522],
       Compile_Date[41], Compile_OS[41], Version_Date[41], Version_Number[41];

  int i;
  int ver_info = 0;
  int help_info = 0;

  int api = 0;
  int log_file = 0;
  int status = 0;
  int verbosity = 1;

  // get program arguments

  if (argc < 2)
    help_info = 1;

  // set input parameters

  for (i = 1; i < argc; i++) {

    if (strcmp (argv[i], "-aflr_api") == 0)
      api = 0;
    else if (strcmp (argv[i], "-egads_api") == 0)
      api = 1;
    else if (strcmp (argv[i], "-log") == 0)
      log_file = 1;
    else if (strcmp (argv[i], "-cpu") == 0)
      verbosity = 2;
    else if (strcmp (argv[i], "-uvout") == 0)
      WriteUVmapOut = 1;
    else if (strcmp (argv[i], "-ver") == 0 || strcmp (argv[i], "--ver") == 0 ||
             strcmp (argv[i], "-build") == 0 || strcmp (argv[i], "--build") == 0)
      ver_info = 1;
    else if (strcmp (argv[i], "-version") == 0 || strcmp (argv[i], "--version") == 0)
      ver_info = 2;
    else if (strcmp (argv[i], "-h") == 0 || strcmp (argv[i], "-help") == 0)
      help_info = 1;
    else if (strncmp (argv[i], "-", 1) == 0) {
      help_info = 1;
      printf ("\nUNKNOWN ARGUMENT %s\n", argv[i]);
    }
    else {
      if (strcmp (Case_Name, "_null_")) {
        if (help_info <= 1)
      	  printf ("\nMORE THAN ONE CASE NAME SET %s\n", Case_Name);
        help_info = 2;
      	printf ("\nMORE THAN ONE CASE NAME SET %s\n", argv[i]);
      }
      else
        strcpy (Case_Name, argv[i]);
    }
  }

  if (ver_info == 0 && strcmp (Case_Name, "_null_") == 0) {
    help_info = 1;
    printf ("\nCASE NAME IS NOT SPECIFIED\n");
  }

  if (help_info == 0 && WriteUVmapOut)
    strcpy (UVmapCaseName, Case_Name);

  // output usage information

  if (help_info) {
    printf ("\n\
uvmap case_name [options]\n\
\n\
Parameter Name        Description\n\
___________________   _____________________________________________________\n\
\n\
-aflr_api		: Use API with AFLR style arrays (default).\n\
-egads_api		: Use API with EGADS style arrays.\n\
-log			: Generate a log file with all output messages.\n\
-cpu			: Generate full CPU usage output messages.\n\
-uvout			: Write internal output files after uv generation\n\
			  case_name.uvmap_ID_uv.surf and\n\
			  case_name.uvmap_ID_xyz.surf.\n\
-h, -help		: Output summary of options.\n\
-ver, --ver		: Output version number.\n\
-version,--version	: Output version number information.\n\
-build,--build		: Output build number.\n\
case_name		: Case Name for input file case_name.surf,\n\
			  output file case_name_uv.surf,\n\
			  and log file (if -log is used) case_name.uvmap.log.\n\
");
    exit (0);
  }

  // output version number

  if (ver_info == 1) {
    uvmap_version (Compile_Date, Compile_OS, Version_Date, Version_Number);
    uvmap_message (Version_Number);
    exit (0);
  }

  // output version information

  if (ver_info == 2) {
    uvmap_version (Compile_Date, Compile_OS, Version_Date, Version_Number);
    printf ("\n");
    printf ("UVMAP    : Version Number %s\n", Version_Number);
    printf ("UVMAP    : Version Date   %s\n", Version_Date);
    printf ("\n");
    printf ("Lib Info : Library Name   uvmap\n");
    printf ("Lib Info : Version Number %s\n", Version_Number);
    printf ("Lib Info : Version Date   %s\n", Version_Date);
    printf ("Lib Info : Compile OS     %s\n", Compile_OS);
    printf ("Lib Info : Compile Date   %s\n", Compile_Date);
    exit (0);
  }

  // open output log file

  if (log_file) {
    sprintf (File_Name, "%s.uvmap.log", Case_Name);
    UVMAP_Output_File = fopen (File_Name, "w");
  }

  // test uvmap with EGADS API

  if (api)
    status = EG_uvmapTest (Case_Name, verbosity);

  // test uvmap with AFLR API

  else
    status = (int) uvmap_test (Case_Name, (INT_) verbosity);

  // close output log file

  if (log_file)
    fclose (UVMAP_Output_File);

  exit (status);
}
