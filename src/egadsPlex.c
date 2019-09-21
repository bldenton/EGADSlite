/*
*        Program is intended to read the Geometric Topology of an EGADSlite geometry file and store it
*        in a Petsc Plex
*
*/

#include "egads.h"

int main(int argc, char *argv[])
{
  int i, j, k, n, nn, stat, oclass, mtype, nbodies, *senses;
  ego context, model, geom, *bodies, *objs, *nobjs;
  
  if (argc != 2) {
    printf(" Usage: liteTest liteFile\n\n");
    exit(EXIT_FAILURE);
  }
  /* initialize */
  printf(" EG_open          = %d\n", EG_open(&context));
  printf(" EG_loadModel     = %d  %s\n", EG_loadModel(context, 0, argv[1],
                                                      &model), argv[1]);
  
  /* test bodyTopo functions */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies,
                        &bodies, &senses);
                        
  printf(" Number of BODIES (nbodies): %d \n", nbodies);
  
  for (i = 0; i < nbodies; i++)
    {
    stat = EG_getBodyTopos(bodies[i], NULL, SHELL, &n, &objs);
    printf("   Number of SHELLS (n): %d \n", n);
    }

  /* Close EGADSlite file */
  printf(" EG_close         = %d\n", EG_close(context));
  return 0;
}