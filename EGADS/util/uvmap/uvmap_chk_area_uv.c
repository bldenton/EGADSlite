#include "UVMAP_LIB.h"

/*
 * UVMAP : TRIA-FACE SURFACE MESH UV MAPPING GENERATOR
 *         DERIVED FROM AFLR4, UG, UG2, and UG3 LIBRARIES
 * $Id: uvmap_chk_area_uv.c,v 1.4 2020/06/13 20:43:17 marcum Exp $
 * Copyright 1994-2020, David L. Marcum
 */

INT_ uvmap_chk_area_uv (
  INT_ err_flag,
  INT_ nbface,
  INT_3D *inibf,
  double tol,
  DOUBLE_2D *u)
{
  // Check tria-face area and return number of tria-faces with negative area in
  // uv space.

  INT_ ibface, inode1, inode2, inode3;
  INT_ n = 0;

  double umax[2], umin[2];
  double area, du;

  // loop over all tria-faces

  for (ibface = 1; ibface <= nbface; ibface++) {

    inode1 = inibf[ibface][0];
    inode2 = inibf[ibface][1];
    inode3 = inibf[ibface][2];

    // set 2X area in uv space

    area = (u[inode2][0] - u[inode1][0]) * (u[inode3][1] - u[inode1][1])
         - (u[inode2][1] - u[inode1][1]) * (u[inode3][0] - u[inode1][0]);

    // set maximum length scale

    umax[0] = MAX (u[inode1][0], u[inode2][0]);
    umax[0] = MAX (u[inode3][0], umax[0]);
    umax[1] = MAX (u[inode1][1], u[inode2][1]);
    umax[1] = MAX (u[inode3][1], umax[1]);
    umin[0] = MIN (u[inode1][0], u[inode2][0]);
    umin[0] = MIN (u[inode3][0], umin[0]);
    umin[1] = MIN (u[inode1][1], u[inode2][1]);
    umin[1] = MIN (u[inode3][1], umin[1]);

    du = MAX (umax[0] - umin[0], umax[1] - umin[1]);

    // check if area is less than tolerance

    if (area < tol * du * du) {
      if (err_flag) {
        uvmap_error_message ("*** ERROR 3500 : tria-face found with negative area in uv space ***");
        return 3500;
      }
      else
        n++;
    }
  }

  // return number of tria-faces with negative area

  return n;
}
