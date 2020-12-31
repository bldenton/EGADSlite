/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Tessellate Sensitivity Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"
#include <math.h>

#include "Surreal/SurrealS.h"


extern "C" int EG_sameThread( const ego object );



extern "C" int EG_tessMassProps(const ego tess, double *props)
{
  int      iedge, iface, ipnt, itri, ip0, ip1, ip2, *tris;
  double   len1, area1, xa, ya, za, xb, yb, zb;
  double   xbar, ybar, zbar, areax, areay, areaz, *xyz;
  double   len=0.0, area=0.0, vol=0.0, xcg=0.0, ycg=0.0, zcg=0.0;
  double   Ixx=0.0, Ixy=0.0,  Ixz=0.0, Iyy=0.0, Iyz=0.0, Izz=0.0;
  ego      body;
  egTessel *btess;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (EG_sameThread(tess))          return EGADS_CNTXTHRD;
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_tessMassProps)!\n");
    return EGADS_NOTFOUND;
  }
  body = btess->src;
  if (body == NULL) {
    printf(" EGADS Error: NULL Source Object (EG_tessMassProps)!\n");
    return EGADS_NULLOBJ;
  }
  if (body->magicnumber != MAGIC) {
    printf(" EGADS Error: Source Not an Object (EG_tessMassProps)!\n");
    return EGADS_NOTOBJ;
  }
  if (body->oclass != BODY) {
    printf(" EGADS Error: Source Not Body (EG_tessMassProps)!\n");
    return EGADS_NOTBODY;
  }
  
  if (body->mtype == WIREBODY) {
    
    for (iedge = 0; iedge < btess->nEdge; iedge++) {
      xyz = btess->tess1d[iedge].xyz;
      for (ipnt = 1; ipnt < btess->tess1d[iedge].npts; ipnt++) {
        ip0 = 3 * (ipnt - 1);
        ip1 = 3 * (ipnt    );
        
        len1 = sqrt((xyz[ip1  ]-xyz[ip0  ]) * (xyz[ip1  ]-xyz[ip0  ]) +
                    (xyz[ip1+1]-xyz[ip0+1]) * (xyz[ip1+1]-xyz[ip0+1]) +
                    (xyz[ip1+2]-xyz[ip0+2]) * (xyz[ip1+2]-xyz[ip0+2]));
        
        len += len1;
        xcg += (xyz[ip1  ] + xyz[ip0  ]) * len1 / 2;
        ycg += (xyz[ip1+1] + xyz[ip0+1]) * len1 / 2;
        zcg += (xyz[ip1+2] + xyz[ip0+2]) * len1 / 2;
      }
    }
    
    xcg /= len;
    ycg /= len;
    zcg /= len;
    
    area = len;
    
  } else if (body->mtype != SOLIDBODY) {
    
    /* FaceBody & SheetBody */
    for (iface = 0; iface < btess->nFace; iface++) {
      xyz  = btess->tess2d[iface].xyz;
      tris = btess->tess2d[iface].tris;
      
      for (itri = 0; itri < btess->tess2d[iface].ntris; itri++) {
        ip0 = 3 * (tris[3*itri  ] - 1);
        ip1 = 3 * (tris[3*itri+1] - 1);
        ip2 = 3 * (tris[3*itri+2] - 1);
        
        xa = xyz[ip1  ] - xyz[ip0  ];
        ya = xyz[ip1+1] - xyz[ip0+1];
        za = xyz[ip1+2] - xyz[ip0+2];
        
        xb = xyz[ip2  ] - xyz[ip0  ];
        yb = xyz[ip2+1] - xyz[ip0+1];
        zb = xyz[ip2+2] - xyz[ip0+2];
        
        xbar = xyz[ip0  ] + xyz[ip1  ] + xyz[ip2  ];
        ybar = xyz[ip0+1] + xyz[ip1+1] + xyz[ip2+1];
        zbar = xyz[ip0+2] + xyz[ip1+2] + xyz[ip2+2];
        
        areax = ya * zb - za * yb;
        areay = za * xb - xa * zb;
        areaz = xa * yb - ya * xb;
        area1 = sqrt(areax*areax + areay*areay + areaz*areaz) / 2;
        
        area += area1;
        xcg  += (xbar * area1) / 3;
        ycg  += (ybar * area1) / 3;
        zcg  += (zbar * area1) / 3;
      }
    }
    
    xcg /= area;
    ycg /= area;
    zcg /= area;
    
  } else {
    
    /* SolidBody */
    for (iface = 0; iface < btess->nFace; iface++) {
      xyz  = btess->tess2d[iface].xyz;
      tris = btess->tess2d[iface].tris;
      
      for (itri = 0; itri < btess->tess2d[iface].ntris; itri++) {
        ip0 = 3 * (tris[3*itri  ] - 1);
        ip1 = 3 * (tris[3*itri+1] - 1);
        ip2 = 3 * (tris[3*itri+2] - 1);
        
        xa = xyz[ip1  ] - xyz[ip0  ];
        ya = xyz[ip1+1] - xyz[ip0+1];
        za = xyz[ip1+2] - xyz[ip0+2];
        
        xb = xyz[ip2  ] - xyz[ip0  ];
        yb = xyz[ip2+1] - xyz[ip0+1];
        zb = xyz[ip2+2] - xyz[ip0+2];
        
        xbar = xyz[ip0  ] + xyz[ip1  ] + xyz[ip2  ];
        ybar = xyz[ip0+1] + xyz[ip1+1] + xyz[ip2+1];
        zbar = xyz[ip0+2] + xyz[ip1+2] + xyz[ip2+2];
        
        areax = ya * zb - za * yb;
        areay = za * xb - xa * zb;
        areaz = xa * yb - ya * xb;
        
        area += sqrt(areax*areax + areay*areay + areaz*areaz) / 2;
        vol  += (       xbar*areax +        ybar*areay +        zbar*areaz)/18;
        xcg  += (xbar/2*xbar*areax + xbar  *ybar*areay + xbar  *zbar*areaz)/54;
        ycg  += (ybar  *xbar*areax + ybar/2*ybar*areay + ybar  *zbar*areaz)/54;
        zcg  += (zbar  *xbar*areax + zbar  *ybar*areay + zbar/2*zbar*areaz)/54;
        Ixx  += (ybar*ybar*ybar*areay + zbar*zbar*zbar*areaz)/162;
        Iyy  += (xbar*xbar*xbar*areax + zbar*zbar*zbar*areaz)/162;
        Izz  += (xbar*xbar*xbar*areax + ybar*ybar*ybar*areay)/162;
        
        Ixy  -= (xbar*ybar*xbar*areax/2 + ybar*xbar*ybar*areay/2 +
                 xbar*ybar*zbar*areaz  )/162;
        Ixz  -= (xbar*zbar*xbar*areax/2 + xbar*zbar*ybar*areay   +
                 zbar*xbar*zbar*areaz/2)/162;
        Iyz  -= (ybar*zbar*xbar*areax   + ybar*zbar*ybar*areay/2 +
                 zbar*ybar*zbar*areaz/2)/162;
      }
    }
    
    xcg /= vol;
    ycg /= vol;
    zcg /= vol;
    
    /* parallel-axis theorem */
    Ixx -= vol * (            ycg * ycg + zcg * zcg);
    Iyy -= vol * (xcg * xcg             + zcg * zcg);
    Izz -= vol * (xcg * xcg + ycg * ycg            );
    
    Ixy += vol * (xcg * ycg      );
    Ixz += vol * (xcg *       zcg);
    Iyz += vol * (      ycg * zcg);
  }
  
  /* store the properties in the array */
  props[ 0] = vol;
  props[ 1] = area;
  props[ 2] = xcg;
  props[ 3] = ycg;
  props[ 4] = zcg;
  props[ 5] = Ixx;
  props[ 6] = Ixy;
  props[ 7] = Ixz;
  props[ 8] = Ixy;
  props[ 9] = Iyy;
  props[10] = Iyz;
  props[11] = Ixz;
  props[12] = Iyz;
  props[13] = Izz;
  
  return EGADS_SUCCESS;
}


static int
EG_local2Global(const egTessel *btess, int index, int local)
{
  if (index < 0) {
    if (btess->tess1d[-index-1].global == NULL) return EGADS_DEGEN;
    if (local > btess->tess1d[-index-1].npts) return EGADS_RANGERR;
    return btess->tess1d[-index-1].global[local] - 1;
  } else {
    if (btess->tess2d[ index-1].global == NULL) return EGADS_DEGEN;
    if (local > btess->tess2d[ index-1].npts) return EGADS_RANGERR;
    return btess->tess2d[ index-1].global[local] - 1;
  }
}


extern "C" int EG_tessMassProps_dot(const ego tess, double *xyz_dot,
                                    double *props, double *props_dot)
{
  int         iedge, iface, ipnt, itri, ip0, ip1, ip2, global, *tris;
  SurrealS<1> len1, area1, xa, ya, za, xb, yb, zb;
  SurrealS<1> xbar, ybar, zbar, areax, areay, areaz, xyz0[3], xyz1[3], xyz2[3];
  SurrealS<1> len=0.0, area=0.0, vol=0.0, xcg=0.0, ycg=0.0, zcg=0.0;
  SurrealS<1> Ixx=0.0, Ixy=0.0,  Ixz=0.0, Iyy=0.0, Iyz=0.0, Izz=0.0;
  double      *xyz;
  ego         body;
  egTessel    *btess;
  
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (EG_sameThread(tess))          return EGADS_CNTXTHRD;
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    printf(" EGADS Error: NULL Blind Object (EG_tessMassProps)!\n");
    return EGADS_NOTFOUND;
  }
  body = btess->src;
  if (body == NULL) {
    printf(" EGADS Error: NULL Source Object (EG_tessMassProps)!\n");
    return EGADS_NULLOBJ;
  }
  if (body->magicnumber != MAGIC) {
    printf(" EGADS Error: Source Not an Object (EG_tessMassProps)!\n");
    return EGADS_NOTOBJ;
  }
  if (body->oclass != BODY) {
    printf(" EGADS Error: Source Not Body (EG_tessMassProps)!\n");
    return EGADS_NOTBODY;
  }
  
  if (body->mtype == WIREBODY) {
    
    for (iedge = 0; iedge < btess->nEdge; iedge++) {
      xyz = btess->tess1d[iedge].xyz;
      for (ipnt = 1; ipnt < btess->tess1d[iedge].npts; ipnt++) {
        ip0             = ipnt - 1;
        ip1             = ipnt;
        global          = EG_local2Global(btess, -iedge-1, ip0);
        if (global < EGADS_SUCCESS) {
          printf(" EGADS Error: %d %d EG_local2Global = %d (EG_tessMassProps)!\n",
                 -iedge-1, ip0+1, global);
          return global;
        }
        xyz0[0]         = xyz[3*ip0  ];
        xyz0[1]         = xyz[3*ip0+1];
        xyz0[2]         = xyz[3*ip0+2];
        xyz0[0].deriv() = xyz_dot[3*global  ];
        xyz0[1].deriv() = xyz_dot[3*global+1];
        xyz0[2].deriv() = xyz_dot[3*global+2];
        global          = EG_local2Global(btess, -iedge-1, ip1);
        if (global < EGADS_SUCCESS) {
          printf(" EGADS Error: %d %d EG_local2Global = %d (EG_tessMassProps)!\n",
                 -iedge-1, ip1+1, global);
          return global;
        }
        xyz1[0]         = xyz[3*ip1  ];
        xyz1[1]         = xyz[3*ip1+1];
        xyz1[2]         = xyz[3*ip1+2];
        xyz1[0].deriv() = xyz_dot[3*global  ];
        xyz1[1].deriv() = xyz_dot[3*global+1];
        xyz1[2].deriv() = xyz_dot[3*global+2];
        
        len1 = sqrt((xyz1[0]-xyz0[0]) * (xyz1[0]-xyz0[0]) +
                    (xyz1[1]-xyz0[1]) * (xyz1[1]-xyz0[1]) +
                    (xyz1[2]-xyz0[2]) * (xyz1[2]-xyz0[2]));
        
        len += len1;
        xcg += (xyz1[0] + xyz0[0]) * len1 / 2;
        ycg += (xyz1[1] + xyz0[1]) * len1 / 2;
        zcg += (xyz1[2] + xyz0[2]) * len1 / 2;
      }
    }
    
    xcg /= len;
    ycg /= len;
    zcg /= len;
    
    area = len;
    
  } else if (body->mtype != SOLIDBODY) {
    
    /* FaceBody & SheetBody */
    for (iface = 0; iface < btess->nFace; iface++) {
      xyz  = btess->tess2d[iface].xyz;
      tris = btess->tess2d[iface].tris;
      
      for (itri = 0; itri < btess->tess2d[iface].ntris; itri++) {
        ip0             = tris[3*itri  ] - 1;
        ip1             = tris[3*itri+1] - 1;
        ip2             = tris[3*itri+2] - 1;
        global          = EG_local2Global(btess, iface+1, ip0);
        if (global < EGADS_SUCCESS) {
          printf(" EGADS Error: %d %d EG_local2Global = %d (EG_tessMassProps)!\n",
                 iface+1, ip0+1, global);
          return global;
        }
        xyz0[0]         = xyz[3*ip0  ];
        xyz0[1]         = xyz[3*ip0+1];
        xyz0[2]         = xyz[3*ip0+2];
        xyz0[0].deriv() = xyz_dot[3*global  ];
        xyz0[1].deriv() = xyz_dot[3*global+1];
        xyz0[2].deriv() = xyz_dot[3*global+2];
        global          = EG_local2Global(btess, iface+1, ip1);
        if (global < EGADS_SUCCESS) {
          printf(" EGADS Error: %d %d EG_local2Global = %d (EG_tessMassProps)!\n",
                 iface+1, ip1+1, global);
          return global;
        }
        xyz1[0]         = xyz[3*ip1  ];
        xyz1[1]         = xyz[3*ip1+1];
        xyz1[2]         = xyz[3*ip1+2];
        xyz1[0].deriv() = xyz_dot[3*global  ];
        xyz1[1].deriv() = xyz_dot[3*global+1];
        xyz1[2].deriv() = xyz_dot[3*global+2];
        global          = EG_local2Global(btess, iface+1, ip2);
        if (global < EGADS_SUCCESS) {
          printf(" EGADS Error: %d %d EG_local2Global = %d (EG_tessMassProps)!\n",
                 iface+1, ip2+1, global);
          return global;
        }
        xyz2[0]         = xyz[3*ip2  ];
        xyz2[1]         = xyz[3*ip2+1];
        xyz2[2]         = xyz[3*ip2+2];
        xyz2[0].deriv() = xyz_dot[3*global  ];
        xyz2[1].deriv() = xyz_dot[3*global+1];
        xyz2[2].deriv() = xyz_dot[3*global+2];
        
        xa = xyz1[0] - xyz0[0];
        ya = xyz1[1] - xyz0[1];
        za = xyz1[2] - xyz0[2];
        
        xb = xyz2[0] - xyz0[0];
        yb = xyz2[1] - xyz0[1];
        zb = xyz2[2] - xyz0[2];
        
        xbar = xyz0[0] + xyz1[0] + xyz2[0];
        ybar = xyz0[1] + xyz1[1] + xyz2[1];
        zbar = xyz0[2] + xyz1[2] + xyz2[2];
        
        areax = ya * zb - za * yb;
        areay = za * xb - xa * zb;
        areaz = xa * yb - ya * xb;
        area1 = sqrt(areax*areax + areay*areay + areaz*areaz) / 2;
        
        area += area1;
        xcg  += (xbar * area1) / 3;
        ycg  += (ybar * area1) / 3;
        zcg  += (zbar * area1) / 3;
      }
    }
    
    xcg /= area;
    ycg /= area;
    zcg /= area;
    
  } else {
    
    /* SolidBody */
    for (iface = 0; iface < btess->nFace; iface++) {
      xyz  = btess->tess2d[iface].xyz;
      tris = btess->tess2d[iface].tris;
      
      for (itri = 0; itri < btess->tess2d[iface].ntris; itri++) {
        ip0             = tris[3*itri  ] - 1;
        ip1             = tris[3*itri+1] - 1;
        ip2             = tris[3*itri+2] - 1;
        global          = EG_local2Global(btess, iface+1, ip0);
        if (global < EGADS_SUCCESS) {
          printf(" EGADS Error: %d %d EG_local2Global = %d (EG_tessMassProps)!\n",
                 iface+1, ip0+1, global);
          return global;
        }
        xyz0[0]         = xyz[3*ip0  ];
        xyz0[1]         = xyz[3*ip0+1];
        xyz0[2]         = xyz[3*ip0+2];
        xyz0[0].deriv() = xyz_dot[3*global  ];
        xyz0[1].deriv() = xyz_dot[3*global+1];
        xyz0[2].deriv() = xyz_dot[3*global+2];
        global          = EG_local2Global(btess, iface+1, ip1);
        if (global < EGADS_SUCCESS) {
          printf(" EGADS Error: %d %d EG_local2Global = %d (EG_tessMassProps)!\n",
                 iface+1, ip1+1, global);
          return global;
        }
        xyz1[0]         = xyz[3*ip1  ];
        xyz1[1]         = xyz[3*ip1+1];
        xyz1[2]         = xyz[3*ip1+2];
        xyz1[0].deriv() = xyz_dot[3*global  ];
        xyz1[1].deriv() = xyz_dot[3*global+1];
        xyz1[2].deriv() = xyz_dot[3*global+2];
        global          = EG_local2Global(btess, iface+1, ip2);
        if (global < EGADS_SUCCESS) {
          printf(" EGADS Error: %d %d EG_local2Global = %d (EG_tessMassProps)!\n",
                 iface+1, ip2+1, global);
          return global;
        }
        xyz2[0]         = xyz[3*ip2  ];
        xyz2[1]         = xyz[3*ip2+1];
        xyz2[2]         = xyz[3*ip2+2];
        xyz2[0].deriv() = xyz_dot[3*global  ];
        xyz2[1].deriv() = xyz_dot[3*global+1];
        xyz2[2].deriv() = xyz_dot[3*global+2];
        
        xa = xyz1[0] - xyz0[0];
        ya = xyz1[1] - xyz0[1];
        za = xyz1[2] - xyz0[2];
        
        xb = xyz2[0] - xyz0[0];
        yb = xyz2[1] - xyz0[1];
        zb = xyz2[2] - xyz0[2];
        
        xbar = xyz0[0] + xyz1[0] + xyz2[0];
        ybar = xyz0[1] + xyz1[1] + xyz2[1];
        zbar = xyz0[2] + xyz1[2] + xyz2[2];
        
        areax = ya * zb - za * yb;
        areay = za * xb - xa * zb;
        areaz = xa * yb - ya * xb;
        
        area += sqrt(areax*areax + areay*areay + areaz*areaz) / 2;
        vol  += (       xbar*areax +        ybar*areay +        zbar*areaz)/18;
        xcg  += (xbar/2*xbar*areax + xbar  *ybar*areay + xbar  *zbar*areaz)/54;
        ycg  += (ybar  *xbar*areax + ybar/2*ybar*areay + ybar  *zbar*areaz)/54;
        zcg  += (zbar  *xbar*areax + zbar  *ybar*areay + zbar/2*zbar*areaz)/54;
        Ixx  += (ybar*ybar*ybar*areay + zbar*zbar*zbar*areaz)/162;
        Iyy  += (xbar*xbar*xbar*areax + zbar*zbar*zbar*areaz)/162;
        Izz  += (xbar*xbar*xbar*areax + ybar*ybar*ybar*areay)/162;
        
        Ixy  -= (xbar*ybar*xbar*areax/2 + ybar*xbar*ybar*areay/2 +
                 xbar*ybar*zbar*areaz  )/162;
        Ixz  -= (xbar*zbar*xbar*areax/2 + xbar*zbar*ybar*areay   +
                 zbar*xbar*zbar*areaz/2)/162;
        Iyz  -= (ybar*zbar*xbar*areax   + ybar*zbar*ybar*areay/2 +
                 zbar*ybar*zbar*areaz/2)/162;
      }
    }
    
    xcg /= vol;
    ycg /= vol;
    zcg /= vol;
    
    /* parallel-axis theorem */
    Ixx -= vol * (            ycg * ycg + zcg * zcg);
    Iyy -= vol * (xcg * xcg             + zcg * zcg);
    Izz -= vol * (xcg * xcg + ycg * ycg            );
    
    Ixy += vol * (xcg * ycg      );
    Ixz += vol * (xcg *       zcg);
    Iyz += vol * (      ycg * zcg);
  }
  
  /* store the properties in the array */
  props[ 0]     = vol.value();
  props[ 1]     = area.value();
  props[ 2]     = xcg.value();
  props[ 3]     = ycg.value();
  props[ 4]     = zcg.value();
  props[ 5]     = Ixx.value();
  props[ 6]     = Ixy.value();
  props[ 7]     = Ixz.value();
  props[ 8]     = Ixy.value();
  props[ 9]     = Iyy.value();
  props[10]     = Iyz.value();
  props[11]     = Ixz.value();
  props[12]     = Iyz.value();
  props[13]     = Izz.value();
  
  props_dot[ 0] = vol.deriv();
  props_dot[ 1] = area.deriv();
  props_dot[ 2] = xcg.deriv();
  props_dot[ 3] = ycg.deriv();
  props_dot[ 4] = zcg.deriv();
  props_dot[ 5] = Ixx.deriv();
  props_dot[ 6] = Ixy.deriv();
  props_dot[ 7] = Ixz.deriv();
  props_dot[ 8] = Ixy.deriv();
  props_dot[10] = Iyz.deriv();
  props_dot[11] = Ixz.deriv();
  props_dot[12] = Iyz.deriv();
  props_dot[13] = Izz.deriv();
  
  return EGADS_SUCCESS;
}
