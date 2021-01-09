/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADS HO Tessellation using wv (from vTess)
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include <string.h>
#include <unistd.h>		// usleep

#include "egads.h"
#include "wsserver.h"

//#define NOVIZ

#ifdef WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#ifndef snprintf
#define snprintf _snprintf
#endif
#endif
#include <winsock2.h>
#endif


/* structure to hold on to the EGADS triangulation per Body */
typedef struct {
  ego *faces;
  ego *edges;
  ego body;
  ego tess;
  int mtype;
  int nfaces;
  int nedges;
} bodyData;


/* extrenal functions not in headers */
  extern int EG_tessHOverts(const ego tess, int nstx, int nItri,
                            const int *iTris, const double *st, ego *newTess);


/* globals used in these functions */
#ifndef NOVIZ
  static wvContext *cntxt;
#endif
  static bodyData  *bodydata;



int main(int argc, char *argv[])
{
  int          i, j, nbody, ibody, stat, oclass, mtype, quads, *senses;
  float        arg, focus[4];
  double       box[6], size, tol, params[3];
  char         *startapp;
  const char   *OCCrev;
  ego          context, model, geom, *bodies, *dum;
  ego          tess;
#ifndef NOVIZ
  int          k, ngp, atype, quad, alen, nseg, len, ntri, sum, *segs;
  float        color[3];
  char         gpname[34];
  const int    *tris, *tric, *ptype, *pindex, *ints;
  const double *xyzs, *uvs, *ts, *reals;
  const char   *string;
  wvData       items[5];
  float        eye[3]      = {0.0, 0.0, 7.0};
  float        center[3]   = {0.0, 0.0, 0.0};
  float        up[3]       = {0.0, 1.0, 0.0};
  static int   sides[3][2] = {{1,2}, {2,0}, {0,1}};
  static int   sideq[4][2] = {{1,2}, {2,5}, {5,0}, {0,1}};
  static int   neigq[4]    = {  0,     3,     4,     2};
#endif

  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if ((argc != 2) && (argc != 5)) {
    printf("\n Usage: vHOtess filename [angle maxlen sag]\n\n");
    return 1;
  }
  
  printf(" Tris (0) or Quads (1): ");
  scanf("%d", &quads);

  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

  /* initialize */
  printf(" EG_open           = %d\n", EG_open(&context));
  printf(" EG_loadModel      = %d\n", EG_loadModel(context, 0, argv[1], &model));
  printf(" EG_getBoundingBox = %d\n", EG_getBoundingBox(model, box));
  printf("       BoundingBox = %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("                     %lf %lf %lf\n", box[3], box[4], box[5]);
  printf(" \n");

                            size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];

  focus[0] = 0.5*(box[0]+box[3]);
  focus[1] = 0.5*(box[1]+box[4]);
  focus[2] = 0.5*(box[2]+box[5]);
  focus[3] = size;

  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                        &bodies, &senses);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getTopology = %d\n", stat);
    return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n", nbody);
  bodydata = (bodyData *) malloc(nbody*sizeof(bodyData));
  if (bodydata == NULL) {
    printf(" MALLOC Error on Body storage!\n");
    return 1;
  }

  params[0] =  0.10*size;
  params[1] =  0.01*size;
  params[2] = 30.0;
  if (argc == 5) {
    sscanf(argv[2], "%f", &arg);
    params[2] = arg;
    sscanf(argv[3], "%f", &arg);
    params[0] = arg;
    sscanf(argv[4], "%f", &arg);
    params[1] = arg;
    printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf\n",
           params[2], params[0], params[1]);
    params[0] *= size;
    params[1] *= size;
  }
  printf(" Reference size = %le\n", size);

  /* fill our structure a body at at time */
  for (ibody = 0; ibody < nbody; ibody++) {
    mtype = 0;
    EG_getTopology(bodies[ibody], &geom, &oclass,
                   &mtype, NULL, &j, &dum, &senses);
    stat = EG_tolerance(bodies[ibody], &tol);
    bodydata[ibody].body  = bodies[ibody];
    bodydata[ibody].mtype = mtype;
    if (mtype == WIREBODY) {
      printf(" Body %2d:  Type = WireBody   tol = %le\n", ibody+1, tol);
    } else if (mtype == FACEBODY) {
      printf(" Body %2d:  Type = FaceBody   tol = %le\n", ibody+1, tol);
    } else if (mtype == SHEETBODY) {
      printf(" Body %2d:  Type = SheetBody  tol = %le\n", ibody+1, tol);
    } else {
      printf(" Body %2d:  Type = SolidBody  tol = %le\n", ibody+1, tol);
    }

    stat = EG_getBodyTopos(bodies[ibody], NULL, FACE,
                           &bodydata[ibody].nfaces, &bodydata[ibody].faces);
    i    = EG_getBodyTopos(bodies[ibody], NULL, EDGE,
                           &bodydata[ibody].nedges, &bodydata[ibody].edges);
    if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
      printf(" EG_getBodyTopos Face = %d\n", stat);
      printf(" EG_getBodyTopos Edge = %d\n", i);
      return 1;
    }

    stat = EG_makeTessBody(bodies[ibody], params, &bodydata[ibody].tess);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTessBody %d = %d\n", ibody, stat);
      continue;
    }
    if (quads == 1) {

      int  itri[24] = { 1,5,9, 1,9,8, 5,2,6, 5,6,9, 8,9,7, 8,7,4, 9,6,3, 9,3,7 };
      double st[18] = { 0.0,0.0, 1.0,0.0, 1.0,1.0, 0.0,1.0,
                        0.5,0.0, 1.0,0.5, 0.5,1.0, 0.0,0.5, 0.5,0.5 };
      
      tess = bodydata[ibody].tess;
      stat = EG_quadTess(tess, &bodydata[ibody].tess);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_quadTess %d = %d  -- reverting...\n", ibody, stat);
        bodydata[ibody].tess = tess;
        continue;
      }
      EG_deleteObject(tess);
      
      tess = bodydata[ibody].tess;
      stat = EG_tessHOverts(tess, -9, 8, itri, st, &bodydata[ibody].tess);
/*
      int  itri[24] = { 1,5,8, 8,5,7, 8,7,4, 5,2,6, 5,6,7, 6,3,7 };
      double st[16] = { 0.0,0.0, 1.0,0.0, 1.0,1.0, 0.0,1.0,
                        0.5,0.0, 1.0,0.5, 0.5,1.0, 0.0,0.5 };
      
      tess = bodydata[ibody].tess;
      stat = EG_quadTess(tess, &bodydata[ibody].tess);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_quadTess %d = %d  -- reverting...\n", ibody, stat);
        bodydata[ibody].tess = tess;
        continue;
      }
      EG_deleteObject(tess);
      
      tess = bodydata[ibody].tess;
      stat = EG_tessHOverts(tess, -8, -6, itri, st, &bodydata[ibody].tess);
*/
      printf(" EG_tessHOverts (quads) = %d\n", stat);
    } else {
/*
      int  itri[12] = { 1,6,5, 6,2,4, 4,5,6, 5,4,3 };
      double st[12] = { 0.0,0.0, 1.0,0.0, 0.0,1.0, 0.5,0.5, 0.0,0.5, 0.5,0.0 };
      
      tess = bodydata[ibody].tess;
      stat = EG_tessHOverts(tess, 6, 4, itri, st, &bodydata[ibody].tess);
*/
      int  itri[27] = { 1,8,7,  8,10,7, 9,10,8, 9,4,10, 9,2,4,
                        7,10,6, 10,5,6, 4,5,10, 6,5,3 };
      double st[20] = { 0.0,0.0, 1.0,0.0, 0.0,1.0, 2./3.,1./3., 1./3.,2./3.,
                        0.0,2./3., 0.0,1./3., 1./3.,0.0, 2./3.,0.0, 1./3.,1./3. };
      
      tess = bodydata[ibody].tess;
      stat = EG_tessHOverts(tess, 10, 9, itri, st, &bodydata[ibody].tess);

      printf(" EG_tessHOverts (tris) = %d\n", stat);
    }
    EG_deleteObject(tess);
    if (stat != EGADS_SUCCESS) exit(1);
  }
  printf(" \n");

#ifndef NOVIZ
  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) {
    printf(" failed to create wvContext!\n");
    for (ibody = 0; ibody < nbody; ibody++) {
      EG_deleteObject(bodydata[ibody].tess);
      EG_free(bodydata[ibody].edges);
      EG_free(bodydata[ibody].faces);
    }
    free(bodydata);

    printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
    printf(" EG_close          = %d\n", EG_close(context));
    return 1;
  }

  /* make the scene */
  for (ngp = sum = stat = ibody = 0; ibody < nbody; ibody++) {

    quad = 0;
    stat = EG_attributeRet(bodydata[ibody].tess, ".tessType", &atype,
                           &alen, &ints, &reals, &string);
    if (stat == EGADS_SUCCESS)
      if (atype == ATTRSTRING)
        if (strcmp(string, "Quad") == 0) quad = 1;

    /* get faces */
    for (i = 0; i < bodydata[ibody].nfaces; i++) {
      stat = EG_getTessFace(bodydata[ibody].tess, i+1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) continue;
      snprintf(gpname, 34, "Body %d Face %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs,  WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      color[0]  = 1.0;
      color[1]  = ibody;
      color[1] /= nbody;
      color[2]  = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      if (quad == 0) {
        for (nseg = j = 0; j < ntri; j++)
          for (k = 0; k < 3; k++)
            if (tric[3*j+k] < j+1) nseg++;
      } else {
        for (nseg = j = 0; j < ntri/2; j++)
          for (k = 0; k < 4; k++)
            if (tric[6*j+neigq[k]] < 2*j+1) nseg++;
      }
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d Sides!\n", nseg);
        continue;
      }
      if (quad == 0) {
        for (nseg = j = 0; j < ntri; j++)
          for (k = 0; k < 3; k++)
            if (tric[3*j+k] < j+1) {
              segs[2*nseg  ] = tris[3*j+sides[k][0]];
              segs[2*nseg+1] = tris[3*j+sides[k][1]];
              nseg++;
            }
      } else {
        for (nseg = j = 0; j < ntri/2; j++)
          for (k = 0; k < 4; k++)
            if (tric[6*j+neigq[k]] < 2*j+1) {
              segs[2*nseg  ] = tris[6*j+sideq[k][0]];
              segs[2*nseg+1] = tris[6*j+sideq[k][1]];
              nseg++;
            }
      }
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_LINDICES, &items[3]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 3!\n", i, gpname);
      free(segs);
/*    color[0] = color[1] = color[2] = 0.8;  */
      color[0] = color[1] = color[2] = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_LCOLOR, &items[4]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 4!\n", i, gpname);
      stat = wv_addGPrim(cntxt, gpname, WV_TRIANGLE,
                         WV_ON|WV_ORIENTATION, 5, items);
      if (stat < 0)
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      if (stat > 0) ngp = stat+1;
      sum += ntri;
    }

    /* get edges */
    color[0] = color[1] = 0.0;
    color[2] = 1.0;
    for (i = 0; i < bodydata[ibody].nedges; i++) {
      stat  = EG_getTessEdge(bodydata[ibody].tess, i+1, &len, &xyzs, &ts);
      if (stat != EGADS_SUCCESS) continue;
      if (len == 0) continue;
      nseg = len-1;
      segs = (int *) malloc(2*nseg*sizeof(int));
      if (segs == NULL) {
        printf(" Can not allocate %d segments for Body %d Edge %d\n",
               nseg, ibody, i+1);
        continue;
      }
      for (j = 0; j < len-1; j++) {
        segs[2*j  ] = j + 1;
        segs[2*j+1] = j + 2;
      }

      snprintf(gpname, 34, "Body %d Edge %d", ibody+1, i+1);
      stat = wv_setData(WV_REAL64, len, (void *) xyzs, WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", i, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS,   &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", i, gpname);
      stat = wv_setData(WV_INT32, 2*nseg, (void *) segs, WV_INDICES,  &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", i, gpname);
      free(segs);
      stat = wv_addGPrim(cntxt, gpname, WV_LINE, WV_ON, 3, items);
      if (stat < 0) {
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
      } else {
        if (cntxt != NULL)
          if (cntxt->gPrims != NULL) {
            cntxt->gPrims[stat].lWidth = 1.5;
            if (wv_addArrowHeads(cntxt, stat, 0.05, 1, &nseg) != 0)
              printf(" wv_addArrowHeads Error\n");
          }
        ngp = stat+1;
      }
    }
  }
  printf(" ** %d gPrims with %d triangles **\n", ngp, sum);

  /* start the server code */

  stat = 0;
  if (wv_startServer(7681, NULL, NULL, NULL, 0, cntxt) == 0) {

    /* we have a single valid server -- stay alive a long as we have a client */
    while (wv_statusServer(0)) {
      usleep(500000);
      if (stat == 0) {
        if (startapp != NULL) system(startapp);
        stat++;
      }
    }
  }
  wv_cleanupServers();
#endif

  /* finish up */
  for (ibody = 0; ibody < nbody; ibody++) {
    EG_deleteObject(bodydata[ibody].tess);
    EG_free(bodydata[ibody].edges);
    EG_free(bodydata[ibody].faces);
  }
  free(bodydata);

  printf(" EG_deleteObject   = %d\n", EG_deleteObject(model));
  printf(" EG_close          = %d\n", EG_close(context));
  return 0;
}


/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *wsi, char *text, /*@unused@*/ int lena)
{
  int          i, j, iBody, ient, stat, nattr, atype, alen;
  const int    *pints;
  char         tag[5];
  const char   *name, *pstr;
  const double *preals;
  ego          obj;

  if (strncmp(text,"Picked: ", 8) == 0) {
    iBody = 0;
    sscanf(&text[12], "%d %s %d", &iBody, tag, &ient);
    if (iBody != 0) {
      printf(" Picked: iBody = %d, type = %s, index = %d\n", iBody, tag, ient);
      if (strcmp(tag,"Face") == 0) {
        obj = bodydata[iBody-1].faces[ient-1];
      } else {
        obj = bodydata[iBody-1].edges[ient-1];
      }
      nattr = 0;
      stat  = EG_attributeNum(obj, &nattr);
      if ((stat == EGADS_SUCCESS) && (nattr != 0)) {
        for (i = 1; i <= nattr; i++) {
          stat = EG_attributeGet(obj, i, &name, &atype, &alen,
                                 &pints, &preals, &pstr);
          if (stat != EGADS_SUCCESS) continue;
          printf("   %s: ", name);
          if ((atype == ATTRREAL) || (atype == ATTRCSYS)) {
            for (j = 0; j < alen; j++) printf("%lf ", preals[j]);
          } else if (atype == ATTRSTRING) {
            printf("%s", pstr);
          } else {
            for (j = 0; j < alen; j++) printf("%d ", pints[j]);
          }
          printf("\n");
        }
      }
    }
  }

}
