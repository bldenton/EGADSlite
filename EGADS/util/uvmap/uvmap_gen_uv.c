#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_gen_uv.c,v 1.13 2021/02/18 06:15:22 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

/*

-----------------------------------------------------
Generate a UV surface mesh given an XYZ surface mesh.
-----------------------------------------------------

INT_ uvmap_gen_uv (
  INT_ verbosity,
  INT_ *nbedge,
  INT_ *nbface,
  INT_ *nnode,
  INT_3D **ibfibf,
  INT_2D **inibe,
  INT_3D **inibf,
  DOUBLE_3D **x,
  DOUBLE_2D **u)


INPUT ARGUMENTS
---------------

verbosity	Message flag.
		If verbosity=0 then do not output progress information.
		If verbosity=1 then output progress information to stdout.
		If verbosity=2 then output progress and additional CPU usage
		information to stdout.

nbface		Number of tria-faces.

nnode		Number of nodes/vertices.

inibf		Tria-face connectivity (nbface+1 in length).

x		XYZ coordinates (nnode+1 in length).
		The XYZ coordinates are used only to determine discontinuous
		locations on the outer and inner (if any) boundary edges.


RETURN VALUE
------------

0		Normal completion without errors.
>0		An error occurred.


OUTPUT ARGUMENTS
----------------

nbedge		Number of outer and inner boundary-edges.

ibfibf		Tria-face neighbor map (nbface+1 in length).

inibe		Boundary-edge connectivity for outer and inner closed curves
		(nbedge+1 in length).

inibf		Tria-face connectivity (nbface+1 in length).
		Connectivity may be reordered for ordering consistency. The
		address may also change as input arrays are temporarily
		reallocated to fill holes if there are inner closed curves.

x		XYZ coordinates (nnode+1 in length).
		The address may change as input arrays are temporarily
		reallocated to fill holes if there are inner closed curves.

u		Generated UV coordinates (nnode+1 in length).

*/

INT_ uvmap_gen_uv (
  INT_ verbosity,
  INT_ *nbedge,
  INT_ *nbface,
  INT_ *nnode,
  INT_3D **ibfibf,
  INT_2D **inibe,
  INT_3D **inibf,
  DOUBLE_3D **x,
  DOUBLE_2D **u)
{
  // Generate a UV surface mesh given an XYZ surface mesh.

  char Compile_Date[41], Compile_OS[41], Version_Date[41], Version_Number[41];
  char Text[512];

  INT_2D *ibeibe = NULL;
  INT_ *ibfibe = NULL;
  INT_ *ibfin = NULL;
  INT_ *iccibe = NULL; 
  INT_ *iccin = NULL;
  INT_ *libfin = NULL;
  INT_ *mben_disc = NULL;

  INT_ bnd_flag, icc, it, nbfacei, nit_min, nit_max, nneg, nnodei, pass, try;
  INT_ c_nit_min = 2;
  INT_ c_nit_max = 3;
  INT_ cpu_timer = 0;
  INT_ nbfpnt = 0;
  INT_ nit = 1000;
  INT_ ncc = 1;
  INT_ npass = 5;
  INT_ ntry = 10;
  INT_ status = 0;
  INT_ xyz_scale = 0;

  double dumax, dumaxb, relax, relax_i, urelaxb_i, w, w_ortho_i;
  double ducnvm = 0.0;
  double angdbe = 30.0;
  double ducnv1 = 0.001;
  double ducnv2 = 0.1;
  double ducnv3 = 0.01;
  double max_ratio_lim = 20.0;
  double tol = 1.0e-14;
  double urelaxb = 0.5;
  double w_ortho = 1000.0;

  w = sin (4.0 * atan (1.0) / (2.0 * sqrt ((double) (*nnode))));
  relax = 0.9 * 2.0 * (1.0 - w) / (1.0 - w * w);
  relax_i = MAX (0.75 * relax, 1.0);

  cpu_timer = (verbosity == 2) ? 1: 0;

  *u = NULL;

  // check maximum edge length ratio and set xyz scaling

  if (*x)
    xyz_scale = uvmap_chk_edge_ratio (*nbface, *inibf, max_ratio_lim, *x);

  // output heading

  if (verbosity) {

    uvmap_cpu_message ("");

    uvmap_version (Compile_Date, Compile_OS, Version_Date, Version_Number);

    uvmap_message ("UVMAP    : ---------------------------------------");
    uvmap_message ("UVMAP    : UVMAP LIBRARY");
    uvmap_message ("UVMAP    : UV MAPPING GENERATOR");
    snprintf (Text, 512, "UVMAP    : Version Number %s", Version_Number);
    uvmap_message (Text);
    snprintf (Text, 512, "UVMAP    : Version Date   %s", Version_Date);
    uvmap_message (Text);
    snprintf (Text, 512, "UVMAP    : Compile OS     %s", Compile_OS);
    uvmap_message (Text);
    snprintf (Text, 512, "UVMAP    : Compile Date   %s", Compile_Date);
    uvmap_message (Text);
    uvmap_message ("UVMAP    : Copyright 1994-2020, D.L. Marcum");
    uvmap_message ("UVMAP    : ---------------------------------------");
    uvmap_message ("");
  }

  if (cpu_timer) {
    uvmap_cpu_timer ("create", "uvmap_bnd_adj");
    uvmap_cpu_timer ("create", "uvmap_solve");
  }

  // determine list of tria-faces attached to a node

  if (status == 0)
    status = uvmap_ibfin (*nbface, *nnode, &nbfpnt, *inibf, &ibfin, &libfin);

  // determine tria-face neighbor map

  if (status == 0)
    status = uvmap_ibfibf (*nbface, ibfin, *inibf, libfin, ibfibf);

  // determine boundary-edge connectivity and tria-face boundary-edge map

  if (status == 0)
    status = uvmap_inibe (*nbface, nbedge, &ibfibe, *ibfibf, *inibf, inibe);

  // determine boundary-edge neighbor map

  if (status == 0)
    status = uvmap_ibeibe (*nbedge, *nnode, *inibe, &ibeibe);

  // set boundary-edge closed curve ID

  if (status == 0)
    status = uvmap_iccibe (*nbedge, &ncc, ibeibe, &iccibe);

  // add a node at the centroid of all inner closed boundary-edge curves and
  // add tria-faces that connect to the centroid node and the nodes of the
  // associated closed boundary-edge curve

  nbfacei = *nbface;
  nnodei = *nnode;

  if (status == 0)
    status = uvmap_add (*nbedge, ncc, nbface, nnode, iccibe, *inibe, inibf, x);

  // if nodes were added then reset list of tria-faces attached to a node

  if (status == 0 && *nnode != nnodei)
    status = uvmap_ibfin (*nbface, *nnode, &nbfpnt, *inibf, &ibfin, &libfin);

  // flag inner (-1) and outer (1) boundary-edge nodes

  if (status == 0)
    status = uvmap_iccin (*nbedge, *nnode, iccibe, &iccin, *inibe);

  // flag discontinuous (1) boundary-edge nodes

  if (status == 0)
    status = uvmap_mben_disc (*nbedge, *nnode, ibeibe, *inibe,  &mben_disc, angdbe, *x);

  // allocate uv coordinates

  *u = (DOUBLE_2D *) uvmap_malloc (&status, ((*nnode)+1)*sizeof(DOUBLE_2D));

  if (status) {
    uvmap_error_message ("*** ERROR 103505 : unable to allocate required memory ***");
    status = 103505;
  }

  // if an error occurred then free temporary arrays and exit

  if (status) {
    uvmap_free (ibeibe);
    uvmap_free (ibfibe);
    uvmap_free (ibfin);
    uvmap_free (iccibe); 
    uvmap_free (iccin);
    uvmap_free (libfin);
    uvmap_free (mben_disc);
    return status;
  }

  // loop over multiple attempts if uv mapping is invalid at completion
  // additional attempts use a progressively reduced boundary under relaxation

  if (verbosity) {
    uvmap_message ("UVMAP    : UV MAPPING GENERATION");
    uvmap_message ("");
    if (xyz_scale) {
      uvmap_message ("UVMAP    : Using xyz scaling");
      uvmap_message ("");
    }
  }

  try = 1;

  do {

    // set uv coordinates on outer curve to a unit circle and all others to the
    // circle center

    uvmap_inl_uv_bnd (*nbedge, *nnode, ibeibe, iccibe, *inibe, *u);

    // iterate solution to convergence

    it = 0;

    do {

      it++;

      dumax = 0.0;

      // include interior boundary curve nodes in solver
      // freeze outer boundary curve nodes in solver

      bnd_flag = 1;

      // solve uncoupled biharmonic Laplacian for initial uv mapping

      if (cpu_timer) uvmap_cpu_timer ("start", "uvmap_solve");

      uvmap_solve (bnd_flag, *nnode, xyz_scale,
                   ibfin, iccin, *inibf, libfin, &dumax, relax, *u, *x);

      if (cpu_timer) uvmap_cpu_timer ("stop", "uvmap_solve");

      if (it == 1) ducnvm = dumax * ducnv1;
    }
    while (it < nit && dumax >= ducnvm);

    if (verbosity) {
      snprintf (Text, 512, "UVMAP    : Iterations, Try   =%10d%6d     interior solution", (int) it, (int) try);
      uvmap_message (Text);
    }

    // set minimum and maximum number of boundary-movement solution iterations

    nit_min = c_nit_min * it;
    nit_max = c_nit_max * it;

    // do at least one overall boundary-movement passes
    // and more if needed to obtain a valid uv mapping

    pass = 1;

    do {

      // use unweighted orthogonality correction if pass = 1
      // use weighted orthogonality correction if pass > 1

      w_ortho_i = (pass == 1) ? 1.0: w_ortho;

      // iterate solution to convergence
      // use boundary orthogonality coupling with outer curve

      it = 0;

      do {

        it++;

        dumax = 0.0;
        dumaxb = 0.0;

        // adjust outer boundary curve nodes for orthogonality

        icc = 1;

        if (cpu_timer) uvmap_cpu_timer ("start", "uvmap_bnd_adj");

        uvmap_bnd_adj (icc, *nbedge, nnodei,
                       ibeibe, ibfibe, ibfin, iccibe, *inibe, *inibf,
                       libfin, mben_disc,
                       &dumaxb, urelaxb, urelaxb, w_ortho_i, *u);

        if (cpu_timer) uvmap_cpu_timer ("stop", "uvmap_bnd_adj");

        // include interior boundary curve nodes in solver if pass = 1
        // freeze interior boundary curve nodes in solver if pass > 1
        // freeze outer boundary curve nodes in solver

        bnd_flag = (pass == 1) ? 1: 0;

        // solve biharmonic Laplacian for uv mapping

        if (cpu_timer) uvmap_cpu_timer ("start", "uvmap_solve");

        uvmap_solve (bnd_flag, *nnode, xyz_scale,
                     ibfin, iccin, *inibf, libfin, &dumax, relax, *u, *x);

        if (cpu_timer) uvmap_cpu_timer ("stop", "uvmap_solve");

        if (it == 1) ducnvm = dumax * ducnv2;
      }
      while (it < nit_max && (it < nit_min || dumax >= ducnvm));

      if (verbosity) {
        snprintf (Text, 512, "UVMAP    : Iterations, Pass  =%10d%6d     outer orthogonal-BC solution", (int) it, (int) pass);
        uvmap_message (Text);
      }

      // if there are inner closed curves

      if (ncc > 1) {

        // use lowered secondary boundary under relaxation after pass 1

        urelaxb_i = (pass == 1) ? urelaxb: 0.5 * urelaxb;

        // if this is the 2nd or more attempts and if this after pass 2 then
        // check validity and generate another uv mapping solution if the
        // existing one is invalid

        if (pass > 2 && try > 1) {

          // check tria-faces for invalid negative area in uv space

          nneg = uvmap_chk_area_uv (0, nbfacei, *inibf, tol, *u);

          // if uv mapping is invalid then generate secondary interior solution

          if (nneg) {

            // iterate solution to convergence

            it = 0;

            do {

              it++;

              dumax = 0.0;

              // include interior boundary curve nodes in solver
              // freeze outer boundary curve nodes in solver
  
              bnd_flag = 1;

              // solve uncoupled biharmonic Laplacian for uv mapping

              if (cpu_timer) uvmap_cpu_timer ("start", "uvmap_solve");

              uvmap_solve (bnd_flag, *nnode, xyz_scale,
                           ibfin, iccin, *inibf, libfin, &dumax, relax, *u, *x);

              if (cpu_timer) uvmap_cpu_timer ("stop", "uvmap_solve");

              if (it == 1) ducnvm = dumax * ducnv1;
            }
            while (it < nit && dumax >= ducnvm);

            if (verbosity) {
              snprintf (Text, 512, "UVMAP    : Iterations, Pass  =%10d%6d     secondary interior solution", (int) it, (int) pass);
              uvmap_message (Text);
            }

            // check tria-faces for invalid negative area in uv space

            nneg = uvmap_chk_area_uv (0, nbfacei, *inibf, tol, *u);

            if (nneg)
              urelaxb_i = urelaxb;
          }
        }

        // determining modified inner boundary & interior node mapping
        // use boundary orthogonality coupling with inner closed curves

        it = 0;

        do {

          it++;

          dumax = 0.0;
          dumaxb = 0.0;

          // loop over interior closed curves

          for (icc = 2; icc <= ncc; icc++) {

            // adjust interior boundary curve icc nodes for orthogonality

            if (cpu_timer) uvmap_cpu_timer ("start", "uvmap_bnd_adj");

            uvmap_bnd_adj (icc, *nbedge, nnodei,
                           ibeibe, ibfibe, ibfin, iccibe, *inibe, *inibf,
                           libfin, mben_disc,
                           &dumaxb, urelaxb, urelaxb_i, w_ortho, *u);

            if (cpu_timer) uvmap_cpu_timer ("stop", "uvmap_bnd_adj");
          }

          // freeze interior boundary curve nodes in solver
          // freeze outer boundary curve nodes in solver

          bnd_flag = 0;

          // solve biharmonic Laplacian for uv mapping

          if (cpu_timer) uvmap_cpu_timer ("start", "uvmap_solve");

          uvmap_solve (bnd_flag, *nnode, xyz_scale,
                       ibfin, iccin, *inibf, libfin, &dumax, relax_i, *u, *x);

          if (cpu_timer) uvmap_cpu_timer ("stop", "uvmap_solve");

          if (it == 1) ducnvm = dumax * ducnv3;

          // exit iteration loop if at least nit_min iterations are complete and
          // the convergence tolerance is met
        }
        while (it < nit_max && (it < nit_min || dumax >= ducnvm));

        if (verbosity) {
          snprintf (Text, 512, "UVMAP    : Iterations, Pass  =%10d%6d     inner orthogonal-BC solution", (int) it, (int) pass);
          uvmap_message (Text);
        }
      }

      // check tria-faces for invalid negative area in uv space

      nneg = uvmap_chk_area_uv (0, nbfacei, *inibf, tol, *u);

      pass++;
    }
    while (pass <= npass && nneg);

    // if uv mapping is invalid then reduce boundary under relaxation factor
    // for another attempt

    if (nneg)
      urelaxb = 0.5 * urelaxb;

    try++;
  }
  while (try <= ntry && nneg);

  // normalize uv coordinates so that all uv values are between zero and one

  uvmap_norm_uv (*nnode, nnodei, *u);

  // reset number of tria-faces and nodes

  *nbface = nbfacei;
  *nnode = nnodei;

  // free temporary array space

  uvmap_free (ibeibe);
  uvmap_free (ibfibe);
  uvmap_free (ibfin);
  uvmap_free (iccibe); 
  uvmap_free (iccin);
  uvmap_free (libfin);
  uvmap_free (mben_disc);

  // check tria-faces for invalid negative area in uv space

  status = uvmap_chk_area_uv (1, *nbface, *inibf, tol, *u);

  if (status)
    return status;

  // output CPU usuage information

  if (cpu_timer) {
    uvmap_message ("");
    uvmap_message ("UVMAP    : CPU USAGE SUMMARY");
    uvmap_message ("");
    uvmap_cpu_timer ("output", "uvmap_bnd_adj");
    uvmap_cpu_timer ("output", "uvmap_solve");
  }

  if (verbosity) {
    uvmap_message ("");
    uvmap_cpu_message ("UVMAP    :");
  }

  return 0;
}
