#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>		// usleep

#include "wsserver.h"

#ifdef WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#ifndef snprintf
#define snprintf _snprintf
#endif
#endif
#include <winsock2.h>
#endif


int main(int argc, char *argv[])
{
  int       i, j, k, n, nf, stat, nvert, mvert, ntri, mtri, tri[3], *segs, *tris;
  float     focus[4], color[3];
  double    size, box[6], xyz[3], *xyzs;
  char      gpname[34], *startapp;
  wvContext *cntxt;
  wvData    items[5];
  float     eye[3]      = {0.0, 0.0, 7.0};
  float     center[3]   = {0.0, 0.0, 0.0};
  float     up[3]       = {0.0, 1.0, 0.0};
  FILE      *fp;

  /* get our starting application line
   *
   * for example on a Mac:
   * setenv WV_START "open -a /Applications/Firefox.app ../client/wv.html"
   */
  startapp = getenv("WV_START");

  if (argc != 2) {
    printf("\n Usage: triServer datafile\n\n");
    return 1;
  }

  fp = fopen(argv[1], "r");
  if (fp == NULL) {
    printf("\n ERROR: Opening %s!\n\n", argv[1]);
    return 1;
  }
  
  /* pre-scan the file */
  mtri   = mvert = 0;
  box[0] = box[1] = box[2] =  1.e100;
  box[3] = box[4] = box[5] = -1.e100;
  fscanf(fp, "%d", &n);
  printf("\n reading %d Bodies...\n", n);
  if (n <= 0) {
    printf("\n ERROR: no data\n\n");
    fclose(fp);
    return 1;
  }
  
  for (i = 0; i < n; i++) {
    fscanf(fp, "%d", &nf);
    printf(" Body %d -- nface = %d\n", i+1, nf);
    for (j = 0; j < nf; j++) {
      fscanf(fp, "%d %d", &nvert, &ntri);
      if (nvert > mvert) mvert = nvert;
      if (ntri  > mtri)  mtri  = ntri;
      for (k = 0; k < nvert; k++) {
        fscanf(fp, "%lf %lf %lf", &xyz[0], &xyz[1], &xyz[2]);
        if (xyz[0] < box[0]) box[0] = xyz[0];
        if (xyz[0] > box[3]) box[3] = xyz[0];
        if (xyz[1] < box[1]) box[1] = xyz[1];
        if (xyz[1] > box[4]) box[4] = xyz[1];
        if (xyz[2] < box[2]) box[2] = xyz[2];
        if (xyz[2] > box[5]) box[5] = xyz[2];
      }
      for (k = 0; k < ntri; k++)
        fscanf(fp, "%d %d %d", &tri[0], &tri[1], &tri[2]);
    }
  }

  printf("\n mVert = %d   mTri = %d\n", mvert, mtri);
  xyzs = (double *) malloc(3*mvert*sizeof(double));
  tris = (int    *) malloc(3*mtri *sizeof(int));
  segs = (int    *) malloc(6*mtri *sizeof(int));
  if ((xyzs == NULL) || (tris == NULL) || (segs == NULL)) {
    printf("\n ERROR: malloc!\n\n");
    if (xyzs != NULL) free(xyzs);
    if (tris != NULL) free(tris);
    if (segs != NULL) free(segs);
    fclose(fp);
    return 1;
  }
  rewind(fp);
  fscanf(fp, "%d", &n);

                            size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];

  focus[0] = 0.5*(box[0]+box[3]);
  focus[1] = 0.5*(box[1]+box[4]);
  focus[2] = 0.5*(box[2]+box[5]);
  focus[3] = size;
  
  /* create the WebViewer context */
  cntxt = wv_createContext(1, 30.0, 1.0, 10.0, eye, center, up);
  if (cntxt == NULL) {
    printf(" failed to create wvContext!\n");
    fclose(fp);
    free(xyzs);
    free(tris);
    free(segs);
    return 1;
  }

  /* make the scene */
  
  for (i = 0; i < n; i++) {
    fscanf(fp, "%d", &nf);
    for (j = 0; j < nf; j++) {
      snprintf(gpname, 34, "Body %d Face %d", i+1, j+1);
      fscanf(fp, "%d %d", &nvert, &ntri);
      
      for (k = 0; k < nvert; k++)
        fscanf(fp, "%lf %lf %lf", &xyzs[3*k  ], &xyzs[3*k+1], &xyzs[3*k+2]);
      for (k = 0; k < ntri; k++) {
        fscanf(fp, "%d %d %d", &tris[3*k  ], &tris[3*k+1], &tris[3*k+2]);
        segs[6*k  ] = tris[3*k  ];
        segs[6*k+1] = tris[3*k+1];
        segs[6*k+2] = tris[3*k+1];
        segs[6*k+3] = tris[3*k+2];
        segs[6*k+4] = tris[3*k+2];
        segs[6*k+5] = tris[3*k  ];
      }
        
      stat = wv_setData(WV_REAL64, nvert, (void *) xyzs,  WV_VERTICES, &items[0]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 0!\n", stat, gpname);
      wv_adjustVerts(&items[0], focus);
      stat = wv_setData(WV_INT32, 3*ntri, (void *) tris, WV_INDICES, &items[1]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 1!\n", stat, gpname);
      /* set the foreground color -- red */
      color[0]  = 1.0;
      color[1]  = 0.0;
      color[2]  = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_COLORS, &items[2]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 2!\n", stat, gpname);
      stat = wv_setData(WV_INT32, 6*ntri, (void *) segs, WV_LINDICES, &items[3]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 3!\n", stat, gpname);
      /* line color -- black */
      color[0] = 0.0;
      color[1] = 0.0;
      color[2] = 0.0;
      stat = wv_setData(WV_REAL32, 1, (void *) color,  WV_LCOLOR, &items[4]);
      if (stat < 0) printf(" wv_setData = %d for %s/item 4!\n", stat, gpname);
      stat = wv_addGPrim(cntxt, gpname, WV_TRIANGLE,
                         WV_ON|WV_ORIENTATION, 5, items);
      if (stat < 0)
        printf(" wv_addGPrim = %d for %s!\n", stat, gpname);
    }
  }
  fclose(fp);

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

  /* finish up */
  free(xyzs);
  free(tris);
  free(segs);
  return 0;
}


/* call-back invoked when a message arrives from the browser */

void browserMessage(/*@unused@*/ void *wsi, /*@unused@*/ char *text,
                    /*@unused@*/ int  lena)
{

}
