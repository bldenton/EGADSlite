static const char help[] = "Test of EGADSLite CAD functionality";

#include <petscdmplex.h>

#include <egads.h>
#include <petsc.h>

typedef struct {
  char filename[PETSC_MAX_PATH_LEN];
} AppCtx;

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  options->filename[0] = '\0';

  ierr = PetscOptionsBegin(comm, "", "EGADSPlex Problem Options", "EGADSLite");CHKERRQ(ierr);
  ierr = PetscOptionsString("-filename", "The EGADSLite file", "ex9.c", options->filename, options->filename, PETSC_MAX_PATH_LEN, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[])
{
  DMLabel        bodyLabel, faceLabel, edgeLabel;
  PetscInt       cStart, cEnd, c;
  /* EGADSLite variables */
  ego            context, model, geom, *bodies, *objs, *nobjs, *mobjs, *lobjs, *fobjs, *eobjs, *shobjs;
  int            oclass, mtype, nbodies, *senses, *bsenses, *shsenses, *fsenses, *lsenses, *esenses;
  int            b;
  /* PETSc variables */
  DM             dm;
  PetscInt       dim = -1, cdim = -1, numCorners = 0, numVertices = 0, numCells = 0, depth = 0;
  PetscInt       numFaces = 0, numEdges = 0;
  PetscInt      *cells  = NULL, *coneOrient = NULL, *cones = NULL, *coneSize = NULL, *numPoints = NULL;
  PetscReal     *coords = NULL;
  MPI_Comm       comm;
  PetscMPIInt    rank;
  AppCtx         ctx;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help); if (ierr) return ierr;
  comm = PETSC_COMM_WORLD;
  ierr = ProcessOptions(comm, &ctx);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) {
    /* Open EGADs file and load EGADs model data */
    ierr = EG_open(&context);CHKERRQ(ierr);
    ierr = EG_loadModel(context, 0, ctx.filename, &model);CHKERRQ(ierr);

    /* test bodyTopo functions */
    ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &bsenses);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", nbodies);CHKERRQ(ierr);

    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nsh, sh, Nf, f, Nl, l, Ne, e, Nv, v;

      /* Output Basic Model Topology */
      ierr = EG_getBodyTopos(body, NULL, SHELL, &Nsh, &shobjs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of SHELLS: %d \n", Nsh);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of FACES: %d \n", Nf);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, LOOP,  &Nl, &lobjs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of LOOPS: %d \n", Nl);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, EDGE,  &Ne, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of EDGES: %d \n", Ne);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, NODE,  &Nv, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of NODES: %d \n", Nv);CHKERRQ(ierr);
      
      /* Get SHELL info which associated with the current BODY */
      ierr = EG_getTopology(body, &geom, &oclass, &mtype, NULL, &Nsh, &shobjs, &shsenses);CHKERRQ(ierr);

      for (sh = 0; sh < Nsh; ++sh){
        ego shell = shobjs[sh];
        int shsense = shsenses[sh];
        
        id   = EG_indexBodyTopo(body, shell);
        ierr = PetscPrintf(PETSC_COMM_SELF, "      SHELL ID: %d :: sense = %d\n", id, shsense);CHKERRQ(ierr);
        
        /* Get FACE info which associated with the current SHELL */
        ierr = EG_getTopology(shell, &geom, &oclass, &mtype, NULL, &Nf, &fobjs, &fsenses);CHKERRQ(ierr);
        
        for (f = 0; f < Nf; ++f){
          ego face = fobjs[f];
          int fsense = fsenses[f];
          
          id   = EG_indexBodyTopo(body, face);
          ierr = PetscPrintf(PETSC_COMM_SELF, "        FACE ID: %d \n", id);CHKERRQ(ierr);
          
          /* Get LOOP info which associated with the current FACE */
          ierr = EG_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &objs, &lsenses);CHKERRQ(ierr);
          
          for (l = 0; l < Nl; ++l) {
            ego loop = objs[l];
            int lsense = lsenses[l];
  \
            id   = EG_indexBodyTopo(body, loop);
            ierr = PetscPrintf(PETSC_COMM_SELF, "          LOOP ID: %d :: sense = %d\n", id, lsense);CHKERRQ(ierr);
    
            /* Get EDGE info which associated with the current LOOP */
            ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &esenses);CHKERRQ(ierr);
    
            for (e = 0; e < Ne; ++e) {
              ego edge = objs[e];
  			      int esense = esenses[e];
    
              id = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_SELF, "            EDGE ID: %d :: sense = %d\n", id, esense);CHKERRQ(ierr);
    
              double range[4] = {0., 0., 0., 0.};
              double point[3] = {0., 0., 0.};
              int    peri;
    
              ierr = EG_getRange(objs[e], range, &peri);
              ierr = PetscPrintf(PETSC_COMM_SELF, "              Range = %lf, %lf, %lf, %lf \n", range[0], range[1], range[2], range[3]);
    
              /* Get NODE info which associated with the current EDGE */
              ierr = EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);CHKERRQ(ierr);
    
              for (v = 0; v < Nv; ++v) {
                ego    vertex = nobjs[v];
                double limits[4];
                int    dummy;
    
                ierr = EG_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses);CHKERRQ(ierr);
                id = EG_indexBodyTopo(body, vertex);
                ierr = PetscPrintf(PETSC_COMM_SELF, "              NODE ID: %d \n", id);CHKERRQ(ierr);
                ierr = PetscPrintf(PETSC_COMM_SELF, "                 (x, y, z) = (%lf, %lf, %lf) \n", limits[0], limits[1], limits[2]);
    
                point[0] = point[0] + limits[0];
                point[1] = point[1] + limits[1];
                point[2] = point[2] + limits[2];
              }
    
#if 0
              point[0] = point[0]/2.;
              point[1] = point[1]/2.;
              point[2] = point[2]/2.;
    
              double trange[2];
    
              trange[0] = 0.;
              trange[1] = 0.;
              double *params[4];
              double *result[3];
              double *xyzresult[9];
              double t=0.;
    
              ierr = EG_nearestOnCurve(objs[e], point, range, t, &xyzresult);
              ierr = PetscPrintf(PETSC_COMM_SELF, " (t1, t2) = (%lf, %lf) \n", params[0], params[1]);
              ierr = PetscPrintf(PETSC_COMM_SELF, " (x, y, z) = (%lf, %lf, %lf) \n", xyzresult[0], xyzresult[1], xyzresult[2]);
#endif
            }
          }
        }
      }
    }
  }
  
  /* Set up DMPlex DAG */
  depth = 2;                                                // define depth of DAG. Here its just FACE --> EGDES --> VERTICES
  ierr = DMCreate(PETSC_COMM_WORLD, &dm); CHKERRQ(ierr);    // Create dm
  ierr = DMSetDimension(dm, depth); CHKERRQ(ierr);          // Set Topological Dimension
  ierr = DMSetCoordinateDim(dm, 3); CHKERRQ(ierr);          // Set Spacial/Coordinate Dimension
  ierr = DMSetType(dm, DMPLEX); CHKERRQ(ierr);            // Set DM Type
  
  /* Determine Number of Total FACES & NODES in model */
  ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &bsenses);CHKERRQ(ierr);
  int Nf, Ne, Nv;
  for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
      ierr = EG_getBodyTopos(body, NULL, EDGE,  &Ne, &eobjs);CHKERRQ(ierr);
      ierr = EG_getBodyTopos(body, NULL, NODE,  &Nv, &nobjs);CHKERRQ(ierr);
      numFaces = numFaces + Nf;
      numEdges = numEdges + Ne;
      numVertices = numVertices + Nv;
  }
  ierr = PetscPrintf(PETSC_COMM_SELF, "numFaces = %d \n", numFaces);CHKERRQ(ierr);    // Print out results
  ierr = PetscPrintf(PETSC_COMM_SELF, "numEdges = %d \n", numEdges);CHKERRQ(ierr);    // Print out results
  ierr = PetscPrintf(PETSC_COMM_SELF, "numVertices = %d \n\n", numVertices);CHKERRQ(ierr);    // Print out results
  
  /* Load numPoints */
  ierr = PetscMalloc1(depth+1, &numPoints); CHKERRQ(ierr);  // Set Size of numPoints Array
  numPoints[0] = numVertices;
  numPoints[1] = numEdges;
  numPoints[2] = numFaces;

  /* Load coneSize */   //Currently this won't work with multiple bodies -- need to rework....
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n coneSize Array \n");CHKERRQ(ierr);
  ierr = PetscMalloc1(numFaces+numEdges+numVertices, &coneSize); CHKERRQ(ierr);     // Set Size of coneSize Array
  int orientSize = 0;    //Defined here to utilize coneSize elements to calculate it coneOrient Array's size
  int vStartIndex = 0;
  for (int ii = 0; ii < numFaces + numEdges + numVertices; ++ii){
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      ego *faceObj, *eobjs;
      int numEnt = 0, numEdge = 0;
      if (ii < numFaces){
        ierr = EG_objectBodyTopo(body, FACE, ii+1, &faceObj); CHKERRQ(ierr);   // Get FACE objects for given ID
        ierr = EG_getBodyTopos(body, faceObj, EDGE, &numEnt, &objs);         // Get Number of NODES and their objects connected to the current FACE    
        coneSize[ii] = numEnt;
        vStartIndex = vStartIndex + numEnt;
      } else if(ii >= numFaces && ii < (numEdges+numFaces)){
        ierr = EG_objectBodyTopo(body, EDGE, ii+1-numFaces, &eobjs); CHKERRQ(ierr);   // Get EDGE objects for given ID
        ierr = EG_getBodyTopos(body, eobjs, NODE, &numEnt, &objs);         // Get Number of NODES and their objects connected to the current FACE    
        coneSize[ii] = numEnt;
      } else {
        coneSize[ii] = 0;    // Store coneSize for Nodes. Note: This is always zero since Nodes are the base (top of the chart) entities
      }
      orientSize = orientSize + numEnt;
    }   
    ierr = PetscPrintf(PETSC_COMM_SELF, "coneSize[%d] = %d \n", ii, coneSize[ii]);CHKERRQ(ierr);    // Print out results
  }
  
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n vStartIndex = %d \n", vStartIndex);CHKERRQ(ierr);    // Print out results
 
   /* Load coneOrient Array */
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n coneOrient Array\n");CHKERRQ(ierr);
  ierr = PetscMalloc1(orientSize, &coneOrient); CHKERRQ(ierr);    // Set Size of coneSize Array
  for (int ii = 0; ii < orientSize; ++ii){
    coneOrient[ii] = 0;
    //ierr = PetscPrintf(PETSC_COMM_SELF, "coneOrient[%d] = %d \n", ii, coneOrient[ii]);CHKERRQ(ierr);    // Print out results
  }
  
  /* Load cones Array */     // Need updating to respect EDGE sense also probably won't work correctly for multiple bodies
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n cones Array\n");CHKERRQ(ierr);
  ierr = PetscMalloc1(orientSize, &cones); CHKERRQ(ierr);    // Set Size of cones Array (same size as coneOrient)
  int ecntr = 0, eID = 0, vcntr = 0;
  for (b = 0; b < nbodies; ++b) {
      //ierr = PetscPrintf(PETSC_COMM_SELF, "b = %d \n", b);CHKERRQ(ierr);
      ego body = bodies[b];
      ego *faceObj = NULL, *nobjs = NULL, *eobjs = NULL, *lobjs = NULL;
      int numNodes = 0, id = 0, Nsh = 0, Nl = 0, Ne = 0;
      
      ierr = EG_getBodyTopos(body, NULL, SHELL, &Nsh, &shobjs);    // Get SHELL info
      //ierr = PetscPrintf(PETSC_COMM_SELF, "Nsh = %d \n", Nsh);CHKERRQ(ierr);
      
      for (int ii = 0; ii < Nsh; ++ii){
        ierr = EG_getBodyTopos(body, shobjs[ii], FACE, &Nf, &faceObj);    // Get FACE info for SHELL
        //ierr = PetscPrintf(PETSC_COMM_SELF, "Nf = %d \n", Nf);CHKERRQ(ierr);
        
        for (int jj = 0; jj < Nf; ++jj){
          ierr = EG_getBodyTopos(body, faceObj[jj], LOOP, &Nl, &lobjs);CHKERRQ(ierr);   // Get LOOP info for FACE
          //ierr = PetscPrintf(PETSC_COMM_SELF, "Nl = %d \n", Nl);CHKERRQ(ierr);
          
          for (int kk = 0; kk < Nl; ++kk){
            ego loop = lobjs[kk];
            
            ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &esenses);CHKERRQ(ierr);
            //ierr = PetscPrintf(PETSC_COMM_SELF, "Ne = %d \n", Ne);CHKERRQ(ierr);
            
            for (int ll = 0; ll < Ne; ++ll){      // one less edge since the DAG assumes closure for what I'm using it for
              ego edge = eobjs[ll];
  			      int esense = esenses[ll];
             
              eID = EG_indexBodyTopo(body, edge);
              cones[ecntr] = eID + (numFaces - 1);
              //ierr = PetscPrintf(PETSC_COMM_SELF, "cones[%d] = %d \n", ecntr, cones[ecntr]);CHKERRQ(ierr);
              ecntr = ecntr + 1;
              //ierr = PetscPrintf(PETSC_COMM_SELF, "ecntr = %d \n", ecntr);CHKERRQ(ierr);
           
               // Added to see if coneOrient is fixed
               if (esense < 0 ){
                 coneOrient[ecntr] = -1;
               }       
               
               ierr = EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);CHKERRQ(ierr);
               
               vcntr = 0;
               // Need to update this loop to correctly assign the cones for edges
               // Have to respect EDGE sense
               //if (esense > 0){
                 for (int mm = 0; mm < Nv; ++mm){
                   ego vertex = nobjs[mm];
                   
                   cones[vStartIndex + (2*(eID-1)) + mm] = EG_indexBodyTopo(body, vertex) + (numFaces + numEdges - 1);
                   //ierr = PetscPrintf(PETSC_COMM_SELF, "cones[%d] = %d \n", vStartIndex + vcntr - 1, cones[vStartIndex + vcntr - 1]);CHKERRQ(ierr);
                   //vcntr = eID;
                   //ierr = PetscPrintf(PETSC_COMM_SELF, "vcntr = %d \n", vcntr);CHKERRQ(ierr);
                 }
               //} else {
               //  //for (int mm = Nv-1; mm >= 0; --mm){
               //  //  ego vertex = nobjs[mm];
               //  //  
               //  //  cones[vStartIndex + (2*(eID-1)) + mm] = EG_indexBodyTopo(body, vertex) + (numFaces + numEdges - 1);
               //  //  //ierr = PetscPrintf(PETSC_COMM_SELF, "cones[%d] = %d \n", vStartIndex + vcntr - 1, cones[vStartIndex + vcntr - 1]);CHKERRQ(ierr);
               //  //  //vcntr = eID;
               //  //  //ierr = PetscPrintf(PETSC_COMM_SELF, "vcntr = %d \n", vcntr);CHKERRQ(ierr);
               //  //}
               //  
               //  cones[vStartIndex + (2*(eID-1))] = EG_indexBodyTopo(body, nobjs[1]) + (numFaces + numEdges - 1);
               //  cones[vStartIndex + (2*(eID-1)) + 1] = EG_indexBodyTopo(body, nobjs[0]) + (numFaces + numEdges - 1);
               //}
               
               
               /*if (ll == 0){
                 //ierr = PetscPrintf(PETSC_COMM_SELF, "counter = %d \n", counter);CHKERRQ(ierr);
                 //if (esense > 0){
                   cones[counter] = EG_indexBodyTopo(body, nobjs[0]) + (numFaces - 1);
                   cones[counter+1] = EG_indexBodyTopo(body, nobjs[1]) + (numFaces - 1);
                 //} else {
                 //  cones[counter] = EG_indexBodyTopo(body, nobjs[1]) + (numFaces - 1);
                 //  cones[counter+1] = EG_indexBodyTopo(body, nobjs[0]) + (numFaces - 1);
                 //}
                 counter = counter + 2;
               } else {
                 //ierr = PetscPrintf(PETSC_COMM_SELF, "counter = %d \n", counter);CHKERRQ(ierr);
                 //if (esense > 0){
                   cones[counter] = EG_indexBodyTopo(body, nobjs[1]) + (numFaces - 1);
                 //} else {
                 //  cones[counter] = EG_indexBodyTopo(body, nobjs[0]) + (numFaces - 1);
                 //}              
                 counter = counter + 1;
               }*/  // To restore previous work only delete /* and */ from the code block above
               //ierr = PetscPrintf(PETSC_COMM_SELF, "counter = %d \n", counter);CHKERRQ(ierr);
            }
          }
        }
      }        
  }
  for (int ii = 0; ii < orientSize; ++ii){      // originally counter = orientSize
    ierr = PetscPrintf(PETSC_COMM_SELF, "cones[%d] = %d \n", ii, cones[ii]);CHKERRQ(ierr);    // Print out results
  }
  
//  /* Load coneOrient Array */
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n coneOrient Array\n");CHKERRQ(ierr);
//  ierr = PetscMalloc1(orientSize, &coneOrient); CHKERRQ(ierr);    // Set Size of coneSize Array
  for (int ii = 0; ii < orientSize; ++ii){
//    coneOrient[ii] = 0;
    ierr = PetscPrintf(PETSC_COMM_SELF, "coneOrient[%d] = %d \n", ii, coneOrient[ii]);CHKERRQ(ierr);    // Print out results
  }
  
  /* Load Vertices */     // Might need to add an IF statment to only save unique vertices
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n coords Array \n");CHKERRQ(ierr);    // Print out results
  ierr = PetscMalloc1(numVertices*3, &coords); CHKERRQ(ierr);    // Set Size of coneSize Array
  for (b = 0; b < nbodies; ++b) {
    ego body = bodies[b];
    int id, Nv, v;

    ierr = EG_getBodyTopos(body, NULL, NODE, &Nv, &nobjs);CHKERRQ(ierr);
    for (v = 0; v < Nv; ++v) {
      ego    vertex = nobjs[v];
      double limits[4];
      int    dummy;

      ierr = EG_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses);CHKERRQ(ierr);
      id   = EG_indexBodyTopo(body, vertex);CHKERRQ(ierr);
      coords[(id-1)*cdim+0] = limits[0];
      coords[(id-1)*cdim+1] = limits[1];
      coords[(id-1)*cdim+2] = limits[2];
      ierr = PetscPrintf(PETSC_COMM_SELF, "    Node ID = %d \n", id);
      ierr = PetscPrintf(PETSC_COMM_SELF, "      (x,y,z) = (%lf, %lf, %lf) \n \n", coords[(id-1)*cdim+0], coords[(id-1)*cdim+1],coords[(id-1)*cdim+2]);
    }
  }
  
  ierr = DMPlexCreateFromDAG(dm, depth, numPoints, coneSize, cones, coneOrient, coords); CHKERRQ(ierr);  
  
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dm \n");CHKERRQ(ierr);
  ierr = DMView(dm, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);  

  ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
  
  ierr = DMViewFromOptions(dm, NULL, "-dm_view");CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);

  /* Close EGADSlite file */
  ierr = EG_close(context);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

  test:
    suffix: sphere_0
    args: -filename ${wPETSC_DIR}/share/petsc/datafiles/meshes/unit_sphere.egadslite -dm_view ::ascii_info_detail

TEST*/
