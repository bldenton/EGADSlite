#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_chk_edge_ratio.c,v 1.1 2021/02/18 06:16:19 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_chk_edge_ratio (
  INT_ nbface,
  INT_3D *inibf,
  double max_ratio_limit,
  DOUBLE_3D *x)
{
  // Check maximum edge length ratio.

  INT_ ibface, inode1, inode2, inode3;

  double dx211, dx212, dx213, dx311, dx312, dx313, dx321, dx322, dx323,
         dxs21, dxs31, dxs32, max_ratio_limit_s, ratio_s;

  max_ratio_limit_s = max_ratio_limit * max_ratio_limit;

  // loop over faces and check maximum edge length ratio

  for (ibface = 1; ibface <= nbface; ibface++) {

    inode1 = inibf[ibface][0];
    inode2 = inibf[ibface][1];
    inode3 = inibf[ibface][2];

    dx211 = x[inode2][0] - x[inode1][0];
    dx212 = x[inode2][1] - x[inode1][1];
    dx213 = x[inode2][2] - x[inode1][2];
    dx311 = x[inode3][0] - x[inode1][0];
    dx312 = x[inode3][1] - x[inode1][1];
    dx313 = x[inode3][2] - x[inode1][2];
    dx321 = x[inode3][0] - x[inode2][0];
    dx322 = x[inode3][1] - x[inode2][1];
    dx323 = x[inode3][2] - x[inode2][2];

    dxs21 = dx211 * dx211 + dx212 * dx212 + dx213 * dx213;
    dxs31 = dx311 * dx311 + dx312 * dx312 + dx313 * dx313;
    dxs32 = dx321 * dx321 + dx322 * dx322 + dx323 * dx323;

    ratio_s = MAX (MAX (dxs21, dxs31), dxs32) / MIN (MIN (dxs21, dxs31), dxs32);

    if (ratio_s >= max_ratio_limit_s)
      return 1;
  }

  return 0;
}
