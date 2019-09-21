/*
*        Program is intended to read the Geometric Topology of an EGADSlite geometry file and store it
*        in a Petsc Plex
*
*/

#include "egads.h"

int main(int argc, char *argv[])
{
  int i, j, k, n, nn, mm, index, stat, oclass, mtype, nbodies, *senses;
  double limits[4];
  ego context, model, geom, *bodies, *objs, *nobjs, *mobjs;
  
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
    stat = EG_getBodyTopos(bodies[i], NULL, SHELL, &n, &objs);  // Get number of SHELLS
    printf("   Number of SHELLS (n): %d \n", n);
    
    stat = EG_getBodyTopos(bodies[i], NULL, FACE, &n, &objs);  // Get number of FACES
    printf("     Number of FACES (n): %d \n", n);
    
    stat = EG_getBodyTopos(bodies[i], NULL, LOOP, &n, &objs);  // Get number of LOOPS
    printf("       Number of LOOPS (n): %d \n", n);
    
    stat = EG_getBodyTopos(bodies[i], NULL, EDGE, &n, &objs);  // Get number of EDGES
    printf("         Number of EDGES (n): %d \n", n);
    
    // Loop through EDGES
    for (j = 0; j < n; j++)
      {
      index = EG_indexBodyTopo(bodies[i], objs[j]);    // Print out EDGE IDs
      printf("          EDGE ID: %d \n", index);
      
      stat = EG_getTopology(objs[j], &geom, &oclass, &mtype, NULL, &nn,
                        &nobjs, &senses);
      
      // Loop through NODES
      for (k = 0; k < nn; nn++)
        {
        stat = EG_getTopology(nobjs[k], &geom, &oclass, &mtype, limits, &mm,
                        &mobjs, &senses);
        
        index = EG_indexBodyTopo(bodies[i], nobjs[k]);    // Print out NODE IDs & coordinates
        printf("          NODE ID: %d \n", index);
        printf("             (x, y, z) = ( %d, %d, %d) \n", limits[0], limits[1], limits[2]);
        } 
      }
    
    stat = EG_getBodyTopos(bodies[i], NULL, NODE, &n, &objs);  // Get number of NODES
    printf("           Number of NODES (n): %d \n", n);
    }

  /* Close EGADSlite file */
  printf(" EG_close         = %d\n", EG_close(context));
  return 0;
}