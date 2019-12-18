static char help[] = "Program is intended to read the Geometric Topology of an EGADSlite geometry file\n\
and store it in a Petsc Plex\n\n\n";

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
  /* EGADSLite variables */
  ego            context, model, geom, *bodies, *objs, *nobjs, *mobjs, *lobjs, *eobjs, *fobjs;
  int            oclass, mtype, nbodies, *senses;
  int            b, maxNumEdgeConnects = 0, maxNumFaceConnects =0;
  /* PETSc variables */
  DM             dm, dmNE, dmNF;
  PetscInt       dim = -1, cdim = -1, numCorners = 0, numVertices = 0, numCells = 0;
  PetscInt      *cells  = NULL, *numPoints = NULL, *cSize = NULL, *cones = NULL, *coneOrient = NULL;
  PetscInt      *numPointsNF = NULL, *cSizeNF = NULL, *conesNF = NULL, *coneOrientNF = NULL;
  PetscInt      *pCone = NULL, *pConeSize = NULL, *pConeNF = NULL, *pConeSizeNF = NULL;
  PetscInt       concatSum = 0, concatSumNF = 0;
  PetscReal     *coords = NULL;
  MPI_Comm       comm;
  PetscMPIInt    rank;
  AppCtx         ctx;
  PetscErrorCode ierr;
  Mat neConnectP, nfConnectP;
  PetscContainer egadsNE, egadsNF, egadsModel;

  ierr = PetscInitialize(&argc, &argv, NULL, help); if (ierr) return ierr;
  comm = PETSC_COMM_WORLD;
  ierr = ProcessOptions(comm, &ctx);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) {
    /* Open EGADs file and load EGADs model data */
    ierr = EG_open(&context);CHKERRQ(ierr);
    ierr = EG_loadModel(context, 0, ctx.filename, &model);CHKERRQ(ierr);

    /* test bodyTopo functions */
    ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", nbodies);CHKERRQ(ierr);
    
    // ---------------------------------------------------------------------------------------------------------------------BD
    /* Create DMPlexes (DAGs) to be used later */
    ierr = DMCreate(PETSC_COMM_WORLD, &dmNE); CHKERRQ(ierr);
    ierr = DMCreate(PETSC_COMM_WORLD, &dmNF); CHKERRQ(ierr);
    
    /* Set Topological Dimension of DMPlexes (DAGs) to be used later */
    ierr = DMSetDimension(dmNE, 1); CHKERRQ(ierr);
    ierr = DMSetDimension(dmNF, 1); CHKERRQ(ierr);
    
    ierr = PetscMalloc1(2, &numPoints); CHKERRQ(ierr);
    ierr = PetscMalloc1(2, &numPointsNF); CHKERRQ(ierr);
    
    /* Set Coordinate/Spacial Dimension of DMPlexes (DAGs) to be used later */
    /* NOTE: EGADS currenlty only recognizes 3D geometry. spacedim ==4 */
    ierr = DMSetCoordinateDim(dmNE, 3); CHKERRQ(ierr);
    ierr = DMSetCoordinateDim(dmNF, 3); CHKERRQ(ierr);    
    
    /* Set Types for DMPlexes (DAGs) to be used later */
    ierr = DMSetType(dmNE, DMPLEX); CHKERRQ(ierr);
    ierr = DMSetType(dmNF, DMPLEX); CHKERRQ(ierr);
    
    // --------------------------------------------------------------------------------------------------------------------BD-end

    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nsh, Nf, Nl, l, Ne, e, Nv, v;    // Ne, Nv removed and moved up
      
      /* Output Basic Model Topology */
      ierr = EG_getBodyTopos(body, NULL, SHELL, &Nsh, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of SHELLS: %d \n", Nsh);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of FACES: %d \n", Nf);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, LOOP,  &Nl, &lobjs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of LOOPS: %d \n", Nl);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, EDGE,  &Ne, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of EDGES: %d \n", Ne);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, NODE,  &Nv, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of NODES: %d \n", Nv);CHKERRQ(ierr);
      
      //-----------------------------------------------------------------------------------------------------------------------------------------------------BD
      /* Get the NODE -> EDGE and NODE -> FACE/SURFACE relations */
      for (v = 0; v < Nv; ++v){
        ego    vertex = objs[v];
                
        int nEdge, nFace;
        
        id = EG_indexBodyTopo(body, vertex);
        ierr = EG_getBodyTopos(body, vertex, EDGE, &nEdge, &eobjs);CHKERRQ(ierr);
        ierr = EG_getBodyTopos(body, vertex, FACE, &nFace, &fobjs);CHKERRQ(ierr);
        
        ierr = PetscPrintf(PETSC_COMM_SELF, " NODE ID: %d \n Is connected to %d EDGES and %d FACES \n", id, nEdge, nFace);CHKERRQ(ierr);
        
        if (nEdge > maxNumEdgeConnects) maxNumEdgeConnects = nEdge;
        if (nFace > maxNumFaceConnects) maxNumFaceConnects = nFace;
      }
      
      ierr = PetscPrintf(PETSC_COMM_SELF, " maxNumEdgeConnects = %d \n maxNumFaceConnects = %d \n", maxNumEdgeConnects, maxNumFaceConnects);CHKERRQ(ierr);
      
      // Set up Petsc matrices -- Not Working
      //Mat neConnectP, nfConnectP;    //Define vectors to store connections
      
      MatCreate(PETSC_COMM_WORLD, &neConnectP);
      MatSetSizes(neConnectP,PETSC_DECIDE, PETSC_DECIDE, Nv, maxNumEdgeConnects+2);
      MatSetType(neConnectP,MATAIJ);
      MatSetUp(neConnectP);
      
      MatCreate(PETSC_COMM_WORLD, &nfConnectP);
      MatSetSizes(nfConnectP,PETSC_DECIDE, PETSC_DECIDE, Nv, maxNumFaceConnects+2);
      MatSetType(nfConnectP,MATAIJ);
      MatSetUp(nfConnectP);
        
      // Define C matrices - may not need anymore    
      int neConnect[Nv][maxNumEdgeConnects+2];
      int nfConnect[Nv][maxNumFaceConnects+2];
        
      // Load matrices
      for (v = 0; v < Nv; ++v){
        ego vertex = objs[v];
        int nEdge, nFace;
        id = EG_indexBodyTopo(body, vertex);
        ierr = EG_getBodyTopos(body, vertex, EDGE, &nEdge, &eobjs); CHKERRQ(ierr);
        ierr = EG_getBodyTopos(body, vertex, FACE, &nFace, &fobjs); CHKERRQ(ierr);
        
        //neConnect[v][0] = id;
        //neConnect[v][1] = nEdge;
        
        //ierr = MatSeqAIJSetPreallocation(neConnectP,maxNumEdgeConnects+2,NULL); CHKERRQ(ierr);
        ierr = MatSetValue(neConnectP,v,0,id,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(neConnectP,v,1,nEdge,INSERT_VALUES); CHKERRQ(ierr);
        
        //nfConnect[v][0] = id;
        //nfConnect[v][1] = nFace;
        
        ierr = MatSetValue(nfConnectP,v,0,id,INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValue(nfConnectP,v,1,nFace,INSERT_VALUES); CHKERRQ(ierr);
        
        int iii;
        for (iii = 0; iii < nEdge; ++iii){
          id = EG_indexBodyTopo(body, eobjs[iii]);
          //neConnect[v][iii+2] = id;
          ierr = MatSetValue(neConnectP, v, iii+2, id, INSERT_VALUES); CHKERRQ(ierr);
        }
        
        for (iii = 0; iii < nFace; ++iii){
          id = EG_indexBodyTopo(body, fobjs[iii]);
          //nfConnect[v][iii+2] = id;
          ierr = MatSetValue(nfConnectP, v, iii+2, id, INSERT_VALUES); CHKERRQ(ierr);
        }
      }
      
      MatAssemblyBegin(neConnectP,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(neConnectP,MAT_FINAL_ASSEMBLY);
      
      MatAssemblyBegin(nfConnectP,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(nfConnectP,MAT_FINAL_ASSEMBLY);
      
      //// output what was stored - DEBUG Statement
      //for (v = 0; v < Nv; ++v){
      //  int iii;
      //  for (iii = 0; iii < maxNumEdgeConnects+2; ++iii){
      //    ierr = PetscPrintf(PETSC_COMM_SELF, "neConnect(%d, %d) = %d \n", v, iii, neConnect[v][iii]);CHKERRQ(ierr);
      //  }
      //}
      //
      //for (v = 0; v < Nv; ++v){
      //  int iii;
      //  for (iii = 0; iii < maxNumFaceConnects+2; ++iii){
      //    ierr = PetscPrintf(PETSC_COMM_SELF, "nfConnect(%d, %d) = %d \n", v, iii, nfConnect[v][iii]);CHKERRQ(ierr);
      //  }
      //}
      
      // View Data Stored in neConnectP & nfConnectP matrices
      ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
      ierr = MatView(neConnectP, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "]\n");CHKERRQ(ierr);
      
      ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
      ierr = MatView(nfConnectP, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "]\n");CHKERRQ(ierr);
      
      //PetscContainer egadsNE, egadsNF;
      ierr = PetscContainerCreate(comm, &egadsNE); CHKERRQ(ierr);
      ierr = PetscContainerCreate(comm, &egadsNF); CHKERRQ(ierr);
      ierr = PetscContainerCreate(comm, &egadsModel); CHKERRQ(ierr);
      ierr = PetscContainerSetPointer(egadsNE, neConnectP); CHKERRQ(ierr);
      ierr = PetscContainerSetPointer(egadsNF, nfConnectP); CHKERRQ(ierr);
      ierr = PetscContainerSetPointer(egadsModel, model); CHKERRQ(ierr);      // model is an EGADS ego object
      //ierr = PetscObjectCompose((PetscObject) neConnectP, "egadsNE", (PetscObject) egadsNE); CHKERRQ(ierr);
      //ierr = PetscObjectCompose((PetscObject) nfConnectP, "egadsNF", (PetscObject) egadsNF); CHKERRQ(ierr);
      //ierr = PetscContainerDestroy(&egadsNE); CHKERRQ(ierr);
      //ierr = PetscContainerDestroy(&egadsNF); CHKERRQ(ierr);
      
      
      /* For DMPlex(DAGs) trial - Single Body*/
      numPoints[0] = Nv;
      numPoints[1] = Ne;
      ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "numPoints = [%d, %d] \n", numPoints[0], numPoints[1]);CHKERRQ(ierr);
      
      ierr = PetscMalloc1(Nv+Ne, &cSize); CHKERRQ(ierr);
      
      int iii, nTrial, idTrial, kkk;
      ego *objTrial, *objTrial2;
      //PetscInt concatSum = 0;
      
      for (iii = 0; iii < Ne+Nv; ++iii) {
        if (iii < Ne){
          ierr = EG_objectBodyTopo(body, EDGE, iii+1, &objTrial); CHKERRQ(ierr);  // Get EDGE object for given ID
          idTrial = EG_indexBodyTopo(body, objTrial);    // Get ID of EDGE
          ierr = PetscPrintf(PETSC_COMM_SELF, "\n EDGE ID = %d = %d \n", idTrial, iii+1);CHKERRQ(ierr);    // Print out results
          ierr = EG_getBodyTopos(body, objTrial, NODE, &nTrial, &objTrial2);    // Get Number of NODES and their objects connected to the current EDGE
          ierr = PetscPrintf(PETSC_COMM_SELF, "%d NODES connected to EDGE %d (%d) \n", nTrial, iii+1, idTrial);CHKERRQ(ierr);
          for (kkk = 0; kkk < nTrial; ++kkk){
            ego nodeTrial = objTrial2[kkk];
            idTrial = EG_indexBodyTopo(body, nodeTrial);
            ierr = PetscPrintf(PETSC_COMM_SELF, "   Node %d = %d \n", kkk+1, idTrial);CHKERRQ(ierr);  
          }
          cSize[iii] = nTrial;  // Store coneSize for Edges   /// Issue with this line of code. Is this protected??
          concatSum = concatSum + nTrial; // Determine required cones array size
        } else {
          cSize[iii] = 0;    // Store coneSize for Nodes. Note: This is always zero sense Nodes are the base (top of the chart) entities
        }
      }
      
      ierr = PetscMalloc2(concatSum, &cones, concatSum, &coneOrient); CHKERRQ(ierr);
      
      // Load cones array with point data
      PetscInt counter = 0; 
      for (iii = 0; iii < Ne; ++iii) {
        ierr = EG_objectBodyTopo(body, EDGE, iii+1, &objTrial); CHKERRQ(ierr);  // Get EDGE object for given ID
        idTrial = EG_indexBodyTopo(body, objTrial);    // Get ID of EDGE
        ierr = EG_getBodyTopos(body, objTrial, NODE, &nTrial, &objTrial2);    // Get Number of NODES and their objects connected to the current EDGE
          for (kkk = 0; kkk < nTrial; ++kkk){
            ego nodeTrial = objTrial2[kkk];
            idTrial = EG_indexBodyTopo(body, nodeTrial);
            cones[counter] = idTrial+(Ne-1);
            coneOrient[counter] = 0;
            counter = counter + 1;  
          }        
      }
          
      // Check what's stored in coneSize[]
      for (iii = 0; iii < Ne+Nv; ++iii) {
        ierr = PetscPrintf(PETSC_COMM_SELF, " coneSize(%d) = %d \n", iii, cSize[iii]); CHKERRQ(ierr);
      }
      
      // Check cones array size
      ierr = PetscPrintf(PETSC_COMM_SELF, " cone array size = %d \n", concatSum); CHKERRQ(ierr);
      
      // Check contents of cones array
      for (iii = 0; iii < concatSum; ++iii) {
        ierr = PetscPrintf(PETSC_COMM_SELF, " cone(%d) = %d \n", iii, cones[iii]); CHKERRQ(ierr);      
      }

      
      numPointsNF[0] = Nv;
      numPointsNF[1] = Nf;
      ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "numPointsNF = [%d, %d] \n", numPointsNF[0], numPointsNF[1]);CHKERRQ(ierr);
      
      ierr = PetscMalloc1(Nv+Nf, &cSizeNF); CHKERRQ(ierr);
      
      int nTrialNF, idTrialNF;
      ego *objTrialNF, *objTrialNF2;
      //PetscInt concatSum = 0;
      
      for (iii = 0; iii < Nf+Nv; ++iii) {
        if (iii < Nf){
          ierr = EG_objectBodyTopo(body, FACE, iii+1, &objTrialNF); CHKERRQ(ierr);  // Get FACE object for given ID
          idTrialNF = EG_indexBodyTopo(body, objTrialNF);    // Get ID of FACE
          ierr = PetscPrintf(PETSC_COMM_SELF, "\n FACE ID = %d = %d \n", idTrialNF, iii+1);CHKERRQ(ierr);    // Print out results
          ierr = EG_getBodyTopos(body, objTrialNF, NODE, &nTrialNF, &objTrialNF2);    // Get Number of NODES and their objects connected to the current FACE
          ierr = PetscPrintf(PETSC_COMM_SELF, "%d NODES connected to FACE %d (%d) \n", nTrialNF, iii+1, idTrialNF);CHKERRQ(ierr);
          for (kkk = 0; kkk < nTrialNF; ++kkk){
            ego nodeTrialNF = objTrialNF2[kkk];
            idTrialNF = EG_indexBodyTopo(body, nodeTrialNF);
            ierr = PetscPrintf(PETSC_COMM_SELF, "   Node %d = %d \n", kkk+1, idTrialNF);CHKERRQ(ierr);  
          }
          cSizeNF[iii] = nTrialNF;  // Store coneSizeNF for FACEs
          concatSumNF = concatSumNF + nTrialNF; // Determine required cones array size
        } else {
          cSizeNF[iii] = 0;    // Store coneSize for Nodes. Note: This is always zero sense Nodes are the base (top of the chart) entities
        }
      }
      
      ierr = PetscMalloc2(concatSumNF, &conesNF, concatSumNF, &coneOrientNF); CHKERRQ(ierr);
      
      // Load cones array with point data
      counter = 0; 
      for (iii = 0; iii < Nf; ++iii) {
        ierr = EG_objectBodyTopo(body, FACE, iii+1, &objTrialNF); CHKERRQ(ierr);  // Get FACE object for given ID
        idTrialNF = EG_indexBodyTopo(body, objTrialNF);    // Get ID of FACE
        ierr = EG_getBodyTopos(body, objTrialNF, NODE, &nTrialNF, &objTrialNF2);    // Get Number of NODES and their objects connected to the current FACE
          for (kkk = 0; kkk < nTrialNF; ++kkk){
            ego nodeTrialNF = objTrialNF2[kkk];
            idTrialNF = EG_indexBodyTopo(body, nodeTrialNF);
            //ierr = PetscPrintf(PETSC_COMM_SELF, " idTrialNF = %d \n", idTrialNF); CHKERRQ(ierr);
            conesNF[counter] = idTrialNF+(Nf-1);
            coneOrientNF[counter] = 0;
            counter = counter + 1;  
          }        
      }
          
      // Check what's stored in coneSize[]
      for (iii = 0; iii < Nf+Nv; ++iii) {
        ierr = PetscPrintf(PETSC_COMM_SELF, " coneSizeNF(%d) = %d \n", iii, cSizeNF[iii]); CHKERRQ(ierr);
      }
      
      // Check cones array size
      ierr = PetscPrintf(PETSC_COMM_SELF, " cone array size = %d \n", concatSumNF); CHKERRQ(ierr);
      
      // Check contents of cones array
      for (iii = 0; iii < concatSumNF; ++iii) {
        ierr = PetscPrintf(PETSC_COMM_SELF, " conesNF(%d) = %d \n", iii, conesNF[iii]); CHKERRQ(ierr);      
      }
      
      
      /* ---------------------- */
       
      //----------------------------------------------------------------------------------------------------------------------------------------BD-end

      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];

        id   = EG_indexBodyTopo(body, loop);
        ierr = PetscPrintf(PETSC_COMM_SELF, "          LOOP ID: %d\n", id);CHKERRQ(ierr);

        /* Get EDGE info which associated with the current LOOP */
        ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses);CHKERRQ(ierr);

        for (e = 0; e < Ne; ++e) {
          ego edge = objs[e];
          
          id = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_SELF, "            EDGE ID: %d\n", id);CHKERRQ(ierr);

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
          }
        }
      }
    }

    /* ---------------------------------------------------------------------------------------------------
    Generate Petsc Plex
      Get all Nodes in model, record coordinates in a correctly formatted array
      Cycle through bodies, cycle through loops, record NODE IDs in a correctly formatted array */

    /* Calculate cell and vertex sizes */
    ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses);CHKERRQ(ierr);
    numCells    = 0;
    numVertices = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nl, l, Nv, v;

      ierr = EG_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs);CHKERRQ(ierr);
      ierr = EG_getBodyTopos(body, NULL, NODE, &Nv, &nobjs);CHKERRQ(ierr);
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];

        id = EG_indexBodyTopo(body, loop);
        /* TODO: Instead of assuming contiguous ids, we could use a hash table */
        numCells = PetscMax(id, numCells);
      }
      for (v = 0; v < Nv; ++v) {
        ego vertex = nobjs[v];

        id = EG_indexBodyTopo(body, vertex);
        /* TODO: Instead of assuming contiguous ids, we could use a hash table */
        numVertices = PetscMax(id, numVertices);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_SELF, "\nPLEX Input Array Checkouts\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " Total Number of Unique Cells    = %d \n", numCells);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " Total Number of Unique Vertices = %d \n", numVertices);CHKERRQ(ierr);

    dim        = 2; /* Assume 3D Models :: Need to update to handle 2D Models in the future */
    cdim       = 3; /* Assume 3D Models :: Need to update to handle 2D Models in the future */
    numCorners = 3; /* TODO Check number of cell corners from EGADSLite */
    ierr = PetscMalloc2(numVertices*cdim, &coords, numCells*numCorners, &cells);CHKERRQ(ierr);

    /* Get vertex coordinates */
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

    /* Get cell vertices by traversing loops */
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nl, l;

      ierr = EG_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs);CHKERRQ(ierr);
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int lid, Ne, e, nc = 0, c;

        lid  = EG_indexBodyTopo(body, loop);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    LOOP ID: %d \n", lid);CHKERRQ(ierr);
        ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses);CHKERRQ(ierr);

        for (e = 0; e < Ne; ++e) {
          ego edge = objs[e];
          int Nv, v;

          id   = EG_indexBodyTopo(body, edge);
          ierr = PetscPrintf(PETSC_COMM_SELF, "      EDGE ID: %d \n", id);CHKERRQ(ierr);
          if (mtype == DEGENERATE) {ierr = PetscPrintf(PETSC_COMM_SELF, "        EGDE %d is DEGENERATE \n", id);CHKERRQ(ierr);}
          ierr = EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);

          /* Add unique vertices to cells, this handles mtype == DEGENERATE fine */
          for (v = 0; v < Nv; ++v) {
            ego vertex = nobjs[v];

            id = EG_indexBodyTopo(body, vertex);
            for (c = 0; c < nc; ++c) if (cells[(lid-1)*numCorners+c] == id-1) break;
            if (c == nc) cells[(lid-1)*numCorners+nc++] = id-1;
          }
        }
        if (nc != numCorners) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid number of cell corners %D, should be %D", nc, numCorners);
        ierr = PetscPrintf(PETSC_COMM_SELF, "      LOOP Corner NODEs (");
        for (c = 0; c < numCorners; ++c) {
          if (c > 0) {ierr = PetscPrintf(PETSC_COMM_SELF, ", ");}
          ierr = PetscPrintf(PETSC_COMM_SELF, "%D", cells[(lid-1)*numCorners+c]);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, ")\n");
      }
    }
  }
  ierr = DMPlexCreateFromCellList(PETSC_COMM_WORLD, dim, numCells, numVertices, numCorners, PETSC_TRUE, cells, cdim, coords, &dm);CHKERRQ(ierr);
// ------------------------------------------------------------------------------------------------------------------------BD 12-6-19
  ierr = DMPlexCreateFromDAG(dmNE, 1, numPoints, cSize, cones, coneOrient, coords); CHKERRQ(ierr);
  ierr = DMPlexCreateFromDAG(dmNF, 1, numPointsNF, cSizeNF, conesNF, coneOrientNF, coords); CHKERRQ(ierr);
// ----------------------------------------------------------------------------------------------------------------------- BD-end 12-6-19
  ierr = PetscFree2(coords, cells);CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dm \n");CHKERRQ(ierr);
  ierr = DMView(dm, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);

// -----------------------------------------------------------------------------------------------------------------------BD 12-6-19 
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmNE \n");CHKERRQ(ierr);
  ierr = DMView(dmNE, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  ierr = DMPlexGetConeSize(dmNE, 0, &pConeSize); CHKERRQ(ierr);
  ierr = DMPlexGetCone(dmNE, 0, &pCone);CHKERRQ(ierr);
  
  PetscInt kk;
  for (kk = 0; kk<pConeSize; ++kk){
    ierr =PetscPrintf(PETSC_COMM_SELF, "\n Plex ID %d connection(%d) = %d", 1, kk, pCone[kk]); CHKERRQ(ierr);
  }
  
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmNF \n");CHKERRQ(ierr);
  ierr = DMView(dmNF, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  ierr = DMPlexGetConeSize(dmNF, 0, &pConeSizeNF); CHKERRQ(ierr);
  ierr = DMPlexGetCone(dmNF, 0, &pConeNF);CHKERRQ(ierr);
  
  //PetscInt kk;
  for (kk = 0; kk<pConeSizeNF; ++kk){
    ierr =PetscPrintf(PETSC_COMM_SELF, "\n Plex ID %d connection(%d) = %d", 1, kk, pConeNF[kk]); CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
// -----------------------------------------------------------------------------------------------------------------------BD 12-6-19 
  
  ierr = PetscObjectCompose((PetscObject) dm, "egadsNE", (PetscObject) neConnectP); CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject) dm, "egadsNF", (PetscObject) nfConnectP); CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject) dm, "egadsModel", (PetscObject) egadsModel); CHKERRQ(ierr);
  ierr = PetscContainerDestroy(&egadsNE); CHKERRQ(ierr);
  ierr = PetscContainerDestroy(&egadsNF); CHKERRQ(ierr);
  ierr = PetscContainerDestroy(&egadsModel); CHKERRQ(ierr);

  //ierr = DMPlexSetRefinementUniform(dm, PETSC_TRUE);CHKERRQ(ierr);
  ierr = DMSetFromOptions(dm);CHKERRQ(ierr);	// refinement
  
  // -----------------------------------BD
  // Get Coordinate Values in the Refined DM
  Vec dmCoords;
  ierr = DMGetCoordinates(dm, &dmCoords); CHKERRQ(ierr);
  
  PetscInt vSize;
  ierr = VecGetSize(dmCoords, &vSize); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n Refined dmCoord size = %d \n \n", vSize);
  
  ierr = VecAssemblyBegin(dmCoords);
  ierr = VecAssemblyEnd(dmCoords);
  
  PetscInt ix[vSize];
  
  PetscInt ii;
  for (ii = 0; ii < vSize; ++ii){
    ix[ii] = ii;
  }
  
  PetscScalar values[vSize];
  ierr = VecGetValues(dmCoords, vSize, ix, &values);
  
  for (ii = 0; ii < vSize; ++ii){
    ierr = PetscPrintf(PETSC_COMM_SELF, "values[%d] = %lf \n", ii, values[ii]);
  }
  
  // Now back out the original nodes from which the new nodes were derived (from here -> below, we ultimately want to incorporate this into plexrefine.c)
  // Right now this is a proof of concept exercise until we can get petsc to link with the egadslite library
  PetscInt refinedVertexStart;
  refinedVertexStart = 3. * numVertices;    //calcute the index of the first refinded node coordinate in values[]
  
  ierr = PetscPrintf(PETSC_COMM_SELF, " refinedVertexStart = %d \n", refinedVertexStart);
  
  /* load coordinates of 1st refined vertex */
  double xCoord, yCoord, zCoord;
  double xCheck, yCheck, zCheck;
  int nodeAID;
  xCoord = values[refinedVertexStart];
  yCoord = values[refinedVertexStart + 1];
  zCoord = values[refinedVertexStart + 2];
  
  ierr = PetscPrintf(PETSC_COMM_SELF, " (x, y, z) = (%lf, %lf, %lf) \n", xCoord, yCoord, zCoord);
  
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  ierr = DMView(dm, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "]\n");CHKERRQ(ierr);
  
  
  //------------------------------------BD-end
  
  ierr = DMViewFromOptions(dm, NULL, "-dm_view");CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);

  /* Close EGADSlite file */
  ierr = EG_close(context);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
