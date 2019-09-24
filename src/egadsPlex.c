/*
*        Program is intended to read the Geometric Topology of an EGADSlite geometry file and store it
*        in a Petsc Plex
*
*/

#include "egads.h"
#include "petsc.h"
#include "petscdmplex.h"
#include "petscsys.h"

int main(int argc, char *argv[])
{
  // Define Variables
  int i, j, k, l, n, ll, nn, mm, nloops, index, stat, oclass, mtype, nbodies, *senses;
  PetscInt dim, numNodes, numLoops;
  PetscInt *plexCells;
  PetscReal *plexNodeCoord;
  double limits[4];
  ego context, model, geom, *bodies, *objs, *nobjs, *mobjs, *lobjs;
  
  // Check for the right number or arguments
  if (argc != 2) {
    printf(" Usage: liteTest liteFile\n\n");
    exit(EXIT_FAILURE);
  }
  
  // Open EGADs file and load EGADs model data
  printf(" EG_open          = %d\n", EG_open(&context));
  printf(" EG_loadModel     = %d  %s\n", EG_loadModel(context, 0, argv[1],
                                                      &model), argv[1]);
  
  /* test bodyTopo functions */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies,
                        &bodies, &senses);
                        
  printf(" Number of BODIES (nbodies): %d \n", nbodies);
  
  // Loop through BODIES
  for (i = 0; i < nbodies; i++)
    {
    // Output Basic Model Topology
    stat = EG_getBodyTopos(bodies[i], NULL, SHELL, &n, &objs);  // Get number of SHELLS
    printf("   Number of SHELLS (n): %d \n", n);
    
    stat = EG_getBodyTopos(bodies[i], NULL, FACE, &n, &objs);  // Get number of FACES
    printf("     Number of FACES (n): %d \n", n);
    
    stat = EG_getBodyTopos(bodies[i], NULL, LOOP, &nloops, &lobjs);  // Get number of LOOPS
    printf("       Number of LOOPS (n): %d \n", nloops);

    stat = EG_getBodyTopos(bodies[i], NULL, EDGE, &l, &objs);  // Get number of EDGES
    printf("         Number of EDGES (n): %d \n", l);
    
    stat = EG_getBodyTopos(bodies[i], NULL, NODE, &n, &objs);  // Get number of NODES
    printf("           Number of NODES (n): %d \n", n);
    
    // Cycle through LOOPS
    for (ll = 0; ll < nloops; ll++)
      {
      index = EG_indexBodyTopo(bodies[i], lobjs[ll]);    // Print out Loop IDs
      printf("          LOOP ID: %d \n", index);
      
      // Get EDGE info which associated with the current LOOP
      stat = EG_getTopology(lobjs[ll], &geom, &oclass, &mtype, NULL, &n,
                        &objs, &senses);
      
      // Cycle through EDGES
      for (j = 0; j < n; j++)
        {
        index = EG_indexBodyTopo(bodies[i], objs[j]);    // Print out EDGE IDs
        printf("            EDGE ID: %d \n", index);
        
        // Get NODE info which associated with the current EDGE
        stat = EG_getTopology(objs[j], &geom, &oclass, &mtype, NULL, &nn,
                          &nobjs, &senses);
        
        // Cycle through NODES
        for (k = 0; k < nn; k++)
          {
          // Get Current NODE data
          stat = EG_getTopology(nobjs[k], &geom, &oclass, &mtype, limits, &mm,
                          &mobjs, &senses);
          
          index = EG_indexBodyTopo(bodies[i], nobjs[k]);    // Print out NODE IDs & coordinates
          printf("              NODE ID: %d \n", index);
          printf("                 (x, y, z) = ( %lf, %lf, %lf) \n", limits[0], limits[1], limits[2]);
          } 
        }
      }
    }
  
  // ---------------------------------------------------------------------------------------------------  
  // Generate Petsc Plex
  //    Get all Nodes in model, record coordinates in a correctly formated array
  //    Cycle through bodies, cycle through loops, recorde NODE IDs in a correctly formateed array
  
  // Get All NODEs in a model
  stat = EG_getTopology(model, &geom, &oclass, &mtype, limits, &nbodies,
                          &bodies, &senses);
  numLoops = 0;
  numNodes = 0;
  for (i = 0; i < nbodies; i++)
    {
    stat = EG_getBodyTopos(bodies[i], NULL, LOOP, &nloops, &lobjs);  // Get LOOP data of current body
    stat = EG_getBodyTopos(bodies[i], NULL, NODE, &n, &nobjs); // Get NODE data of current Body
    
    // Determine the total number of LOOPs in the model
    for (j = 0; j < nloops; j++)
      {
      index = EG_indexBodyTopo(bodies[i], lobjs[j]);  // Get ID of current LOOP
      
      if (index > numLoops)
        {
        numLoops = index;
        }
      else
        {
        // Do Nothing
        }
      }
    
    // Determine the total number of NODEs in the model
    for (j = 0; j < n; j++)
      {
      index = EG_indexBodyTopo(bodies[i], nobjs[j]);  // Get ID of current NODE
      
      if (index > numNodes) 
        {
        numNodes = index;
        }
      else
        {
        // Do Nothing
        }
      }
    }
    
  printf(" PLEX Input Array Checkouts \n");
  
  // Output the total number of unique LOOPs and NODEs
  printf(" Total Number of Unique Loops = %d \n", numLoops);
  printf(" Total Number of Unique Nodes = %d \n", numNodes);
  
  // Define plexNodeCoord[] Array size
  dim = 3;    // Assumed 3D Models :: Need to update to handle 2D Models in the future
  PetscMalloc1(dim*numNodes, &plexNodeCoord);
  PetscMalloc1(3*numLoops, &plexCells);  // Only Triangle Cells current understood
                                         // Also all cells must have the same number of corner vertices
                                         // based on the current functionality of Petsc Plexes
  
  // Get Current NODE coordinates data by cycling through BODIES
  // and load plexNodeCoord for plex
  for (i = 0; i < nbodies; i++)
    {
    stat = EG_getBodyTopos(bodies[i], NULL, NODE, &n, &nobjs); // Get NODE data of curren Body
    
    for (j = 0; j < n; j++)
      {
      stat = EG_getTopology(nobjs[j], &geom, &oclass, &mtype, limits, &mm,
                      &mobjs, &senses);
                      
      index = EG_indexBodyTopo(bodies[i], nobjs[j]);    // Print out NODE IDs & coordinates
      
      plexNodeCoord[dim*(index-1)+0] = limits[0];  // Node x-coordinate
      plexNodeCoord[dim*(index-1)+1] = limits[1];  // Node y-coordinate
      plexNodeCoord[dim*(index-1)+2] = limits[2];  // Node z-coordinate
      
      // Checkout Statement
      printf("    Node ID = %d \n", index);
      printf("      (x,y,z) = (%lf, %lf, %lf) \n", plexNodeCoord[dim*(index-1)+0],plexNodeCoord[dim*(index-1)+1],plexNodeCoord[dim*(index-1)+2]);   
      }
    }
    
   // Get all LOOPs in the Model
   // We are assuming LOOPs are equivalent to FACEs or SURFACEs   
    
    
    
  /* Close EGADSlite file */
  printf(" EG_close         = %d\n", EG_close(context));
  return 0;
}