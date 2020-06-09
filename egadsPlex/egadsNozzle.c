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
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel; // markerLabel;
  PetscInt       cStart, cEnd;
  /* EGADSLite variables */
  ego            context, model, geom, *bodies, *objs, *nobjs, *mobjs, *lobjs, *fobjs, *eobjs, *shobjs;
  ego            *tess;
  int            oclass, mtype, nbodies, *senses, *bsenses, *shsenses, *fsenses, *lsenses, *esenses;
  int            b;
  /* PETSc variables */
  DM             dm, dmNozzle, dmMesh;
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
              
              /* Debug Trial */
              //double xyzCheck[18];
              //double tm[1];
              //tm[0] = range[0];
              //ierr = EG_evaluate(objs[e], tm, &xyzCheck); CHKERRQ(ierr);
              //ierr = PetscPrintf(PETSC_COMM_SELF, "              tmin = %lf :: (x, y, z) = (%lf, %lf, %lf) \n", range[0], xyzCheck[0], xyzCheck[1], xyzCheck[2]); 
              //tm[0] = range[1];
              //ierr = EG_evaluate(objs[e], tm, &xyzCheck); CHKERRQ(ierr);
              //ierr = PetscPrintf(PETSC_COMM_SELF, "              tmax = %lf :: (x, y, z) = (%lf, %lf, %lf) \n", range[1], xyzCheck[0], xyzCheck[1], xyzCheck[2]);
    
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
  
  /* Matrix to hold face mid-points*/
  
    ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &bsenses);CHKERRQ(ierr);
    int Nf=0, Ne=0, Nv=0;
    for (b = 0; b < nbodies; ++b) {
        ego body = bodies[b];
        ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
        ierr = EG_getBodyTopos(body, NULL, EDGE,  &Ne, &eobjs);CHKERRQ(ierr);
        ierr = EG_getBodyTopos(body, NULL, NODE,  &Nv, &nobjs);CHKERRQ(ierr);
        numFaces = numFaces + Nf;
        numEdges = numEdges + Ne;
        numVertices = numVertices + Nv;
    }
    
    int Nftotal = numFaces;
    int Netotal = numEdges;
    int Nvtotal = numVertices;
    
    Mat    faceMidP;
    MatCreate(PETSC_COMM_WORLD, &faceMidP);
    MatSetSizes(faceMidP,PETSC_DECIDE, PETSC_DECIDE, Nf, 4);    // x, y, z, of face center
    MatSetType(faceMidP,MATAIJ);
    MatSetUp(faceMidP);
    
    Mat    edgeMidP;
    MatCreate(PETSC_COMM_WORLD, &edgeMidP);
    MatSetSizes(edgeMidP,PETSC_DECIDE, PETSC_DECIDE, Ne, 4);    // x, y, z, of edge center
    MatSetType(edgeMidP,MATAIJ);
    MatSetUp(edgeMidP); 
    
    // load the face array
    for (b = 0; b < nbodies; ++b){
        ego body = bodies[b];
        
        /* MidPoint on Faces*/
        ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
        
        for (int ii = 0; ii < Nf; ++ii){
          ego face = fobjs[ii];
          double range[4];
          int    *periodic;
          
          ierr = EG_getRange(face, &range, &periodic); CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_SELF, "FACE %d range = %lf, %lf, %lf, %lf \n", ii+1, range[0], range[1], range[2], range[3]);CHKERRQ(ierr);
          
          double avgUV[2];
          avgUV[0] = (range[0] + range[1]) / 2.;
          avgUV[1] = (range[2] + range[3]) / 2.;
          
          double cntrPnt[18];
          ierr = EG_evaluate(face, avgUV, &cntrPnt);
          //ierr = PetscPrintf(PETSC_COMM_SELF, "FACE %d cntrPnt = %lf, %lf, %lf \n", ii+1, cntrPnt[0], cntrPnt[1], cntrPnt[2]);CHKERRQ(ierr);
          
          ierr = MatSetValue(faceMidP,ii,0,Nv+Ne+ii,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValue(faceMidP,ii,1,cntrPnt[0],INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValue(faceMidP,ii,2,cntrPnt[1],INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValue(faceMidP,ii,3,cntrPnt[2],INSERT_VALUES); CHKERRQ(ierr);
          //faceMidP[ii][0] = cntrPnt[0];
          //faceMidP[ii][1] = cntrPnt[0];
          //faceMidP[ii][2] = cntrPnt[0];
        }
        
        /* MidPoint on Edges */
        ierr = EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
        
        for (int ii = 0; ii < Ne; ++ii){
        ego edge = eobjs[ii];
          double range[2];
          int    *periodic;
          
          ierr = EG_getRange(edge, &range, &periodic); CHKERRQ(ierr);
          //ierr = PetscPrintf(PETSC_COMM_SELF, "EDGE %d range = %lf, %lf \n", ii+1, range[0], range[1]);CHKERRQ(ierr);
          
          double avgt[1];
          avgt[0] = (range[0] + range[1]) / 2.;
          
          double cntrPnt[9];
          ierr = EG_evaluate(edge, avgt, &cntrPnt);
          //ierr = PetscPrintf(PETSC_COMM_SELF, "EDGE %d cntrPnt = %lf, %lf, %lf \n", ii+1, cntrPnt[0], cntrPnt[1], cntrPnt[2]);CHKERRQ(ierr);
          
          ierr = MatSetValue(edgeMidP,ii,0,Nv+ii,INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValue(edgeMidP,ii,1,cntrPnt[0],INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValue(edgeMidP,ii,2,cntrPnt[1],INSERT_VALUES); CHKERRQ(ierr);
          ierr = MatSetValue(edgeMidP,ii,3,cntrPnt[2],INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    
    MatAssemblyBegin(faceMidP,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(faceMidP,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(edgeMidP,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(edgeMidP,MAT_FINAL_ASSEMBLY);
    
    //ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
    //ierr = MatView(faceMidP, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_SELF, "]\n");CHKERRQ(ierr);
    //
    //ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
    //ierr = MatView(edgeMidP, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_SELF, "]\n");CHKERRQ(ierr);
    
    
    // -------------------------------------------------------------------------------
    // Start to setup DMPlex
    // -------------------------------------------------------------------------------
    dim = 2;
    cdim = 3;
    numCorners = 3;        // Assumes Triangle cells
    numVertices = Nv + Ne + Nf;
    
    // Determine Number of Cells assuming Triangular Cells (3 nodes)
    numCells = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      ierr = EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs); CHKERRQ(ierr);
      
      for(int f = 0; f < Nf; ++f){
        ego face = fobjs[f];
        ierr = EG_getBodyTopos(body, face, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
        
        numCells = numCells + (2 * Ne); 
      }
    }
    
    //ierr = PetscPrintf(PETSC_COMM_SELF, "numCells = %d \n", numCells);CHKERRQ(ierr);
    PetscReal *uParams = NULL, *vParams = NULL, *tParams = NULL;
    ierr = PetscMalloc2(numVertices*cdim, &coords, numCells*numCorners, &cells);CHKERRQ(ierr);
    ierr = PetscMalloc2(numVertices, &uParams, numVertices, &vParams); CHKERRQ(ierr);
    ierr = PetscMalloc1(numVertices, &tParams); CHKERRQ(ierr);
    
    /* Load coordinate array */
    // First vertices from CAD model
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
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    Node ID = %d \n", id-1);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "      (x,y,z) = (%lf, %lf, %lf) \n \n", coords[(id-1)*cdim+0], coords[(id-1)*cdim+1],coords[(id-1)*cdim+2]);
      }
    }
    
    // Second EDGE midpoints
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      ierr = EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
          
      for (int ii = 0; ii < Ne; ++ii){
        ego edge = eobjs[ii];
        double range[2];
        int    *periodic;
        int     id;
        
        id   = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);
        ierr = EG_getRange(edge, &range, &periodic); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "EDGE %d range = %lf, %lf \n", ii+1, range[0], range[1]);CHKERRQ(ierr);
        
        double avgt[1];
        avgt[0] = (range[0] + range[1]) / 2.;
        
        double cntrPnt[9];
        ierr = EG_evaluate(edge, avgt, &cntrPnt);
        
        // changed ii to id-1
        coords[(Nv+id-1)*cdim+0] = cntrPnt[0];
        coords[(Nv+id-1)*cdim+1] = cntrPnt[1];
        coords[(Nv+id-1)*cdim+2] = cntrPnt[2];
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    Node ID = %d \n", ii+Nv);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "      (x,y,z) = (%lf, %lf, %lf) \n \n", coords[(ii+Nv)*cdim+0], coords[(ii+Nv)*cdim+1],coords[(ii+Nv)*cdim+2]);
      }
    }
    
    // Third FACE midpoints
    for (b = 0; b < nbodies; ++b){
      ego body = bodies[b];
      
      /* MidPoint on Faces*/
      ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
      
      for (int ii = 0; ii < Nf; ++ii){
        ego face = fobjs[ii];
        double range[4];
        int    *periodic;
        int     id;
        
        id   = EG_indexBodyTopo(body, face);CHKERRQ(ierr);
        ierr = EG_getRange(face, &range, &periodic); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "FACE %d range = %lf, %lf, %lf, %lf \n", ii+1, range[0], range[1], range[2], range[3]);CHKERRQ(ierr);
        
        double avgUV[2];
        avgUV[0] = (range[0] + range[1]) / 2.;
        avgUV[1] = (range[2] + range[3]) / 2.;
        
        double cntrPnt[18];
        ierr = EG_evaluate(face, avgUV, &cntrPnt);
        
        // changed ii to id-1
        coords[(Nv+Ne+id-1)*cdim+0] = cntrPnt[0];
        coords[(Nv+Ne+id-1)*cdim+1] = cntrPnt[1];
        coords[(Nv+Ne+id-1)*cdim+2] = cntrPnt[2];
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    Node ID = %d \n", ii+Nv+Ne);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "      (x,y,z) = (%lf, %lf, %lf) \n \n", coords[(ii+Nv+Ne)*cdim+0], coords[(ii+Nv+Ne)*cdim+1],coords[(ii+Nv+Ne)*cdim+2]);
      }
    }
    
    // Define Cells  -- This is most likely not optimized   
    //ierr = PetscPrintf(PETSC_COMM_SELF, " -------------------------------\n");
    //ierr = PetscPrintf(PETSC_COMM_SELF, "        DEFINING CELLS          \n");
    //ierr = PetscPrintf(PETSC_COMM_SELF, " -------------------------------\n");
    int cellCntr = 0;
    for (b = 0; b < nbodies; ++b){
      ego body = bodies[b];
      ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
      
      for (int f = 0; f < Nf; ++f){
        ego face = fobjs[f];
        int Nl;
        //ierr = EG_getBodyTopos(body, NULL, LOOP,  &Nl, &lobjs);CHKERRQ(ierr);      // Either this line or the next line. Using the next line works better (don't know why)
        ierr = EG_getBodyTopos(body, face, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
        //ierr = EG_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &lobjs, &lsenses);CHKERRQ(ierr);
        ego loop = lobjs[0];
         
          int fid = EG_indexBodyTopo(body, face);CHKERRQ(ierr);
          
          int midFaceID = Nvtotal + Netotal + fid-1;    // fid-1 was fid
                    
          for (int e = 0; e < Ne; ++e){
            ego edge = eobjs[e];
            int id   = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);  // ID of current edge
            
            int midPntID = Nvtotal + id - 1;
            
            ierr = EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);CHKERRQ(ierr);
            
            int startID = 0, endID = 0;
            
            if (esenses[e] > 0){
              startID = EG_indexBodyTopo(body, nobjs[0]);CHKERRQ(ierr);  // ID of EDGE start NODE
              endID = EG_indexBodyTopo(body, nobjs[1]);CHKERRQ(ierr);  // ID of EDGE end NODE
            } else {
              startID = EG_indexBodyTopo(body, nobjs[1]);CHKERRQ(ierr);  // ID of EDGE start NODE
              endID = EG_indexBodyTopo(body, nobjs[0]);CHKERRQ(ierr);  // ID of EDGE end NODE
            }
            
            cells[cellCntr*numCorners + 0] = midFaceID;
            cells[cellCntr*numCorners + 1] = startID-1;
            cells[cellCntr*numCorners + 2] = midPntID;
            
            cells[cellCntr*numCorners + 3] = midFaceID;
            cells[cellCntr*numCorners + 4] = midPntID;
            cells[cellCntr*numCorners + 5] = endID-1;   
            
            cellCntr = cellCntr + 2;     
          } 
        //}  
      }
    }
       
    //Build DMPlex  
    ierr = DMPlexCreateFromCellList(PETSC_COMM_WORLD, dim, numCells, numVertices, numCorners, PETSC_TRUE, cells, cdim, coords, &dmNozzle);CHKERRQ(ierr); 
    ierr = PetscFree2(coords, cells);CHKERRQ(ierr);   
    
    // -------------------------------------------------------------------------
    // Attempt to get Coordinates back out
    //   - Trial for plexegads.c -
    // -------------------------------------------------------------------------
    Vec coordsVec;
    ierr = VecCreate(PETSC_COMM_WORLD, &coordsVec); CHKERRQ(ierr);
    ierr = VecSetType(coordsVec, VECSTANDARD); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dmNozzle, &coordsVec); CHKERRQ(ierr);
    PetscInt vecSize;
    ierr = VecGetSize(coordsVec, &vecSize); CHKERRQ(ierr);
    
    int vecCntr = 1;
    for ( int ii = 0; ii < vecSize; ii+=3){
      PetscInt ix[3], ni;
      ni = 3;
      PetscScalar xx[3]; //, yy, zz;
      ix[0] = ii;
      ix[1] = ii+1;
      ix[2] = ii+2;
      ierr = VecGetValues(coordsVec, ni, ix, &xx);    // was &zz

      vecCntr = vecCntr + 1;
    }
    
    // Check out finding support information
    PetscInt sStart, sEnd;
    PetscInt sConeSize = 0, sConeSizeN = 0, cConeSize = 0;
    PetscInt *sCone = NULL, *sConeN = NULL, *sConeOrient = NULL, *cCone = NULL, *cConeOrient = NULL;
    
    ierr = DMPlexGetHeightStratum(dmNozzle, 1, &sStart, &sEnd);CHKERRQ(ierr);
    for (int ii = sStart; ii < sEnd; ++ii){
      ierr = DMPlexGetSupportSize(dmNozzle, ii, &sConeSize); CHKERRQ(ierr);
      ierr = DMPlexGetSupport(dmNozzle, ii, &sCone); CHKERRQ(ierr);      
      ierr = DMPlexGetConeSize(dmNozzle, sCone[0], &cConeSize); CHKERRQ(ierr);
      ierr = DMPlexGetCone(dmNozzle, sCone[0], &cCone); CHKERRQ(ierr);
      ierr = DMPlexGetConeOrientation(dmNozzle, sCone[0], &cConeOrient); CHKERRQ(ierr);
    }
    
  
  {
    PetscContainer modelObj;

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &modelObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(modelObj, model);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject) dmNozzle, "EGADS Model", (PetscObject) modelObj);CHKERRQ(ierr);
    ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);
  }
  
  ierr = DMCreateLabel(dmNozzle, "EGADS Body ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmNozzle, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dmNozzle, "EGADS Face ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmNozzle, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dmNozzle, "EGADS Edge ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmNozzle, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dmNozzle, "EGADS Vertex ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmNozzle, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
  
  /* Set Label Values - EGADS body*/
  for (int jj = 0; jj < 3; ++jj){
    ierr = DMPlexGetHeightStratum(dmNozzle, jj, &cStart, &cEnd);CHKERRQ(ierr);
    for (int ii = cStart; ii < cEnd; ++ii){
      ierr = DMLabelSetValue(bodyLabel, ii, 0); CHKERRQ(ierr);    // Need to change this to be more flexible for multi-body parts
      ierr = DMLabelSetValue(edgeLabel, ii, -1); CHKERRQ(ierr);
      ierr = DMLabelSetValue(faceLabel, ii, -1); CHKERRQ(ierr);
      ierr = DMLabelSetValue(vertexLabel, ii, -1); CHKERRQ(ierr);
    }
  }
  
  /* Get Number of DAG Nodes at each level */
  int fDAGlevel, eDAGlevel, nDAGlevel;
  int fStart, fEnd, eStart, eEnd, nStart, nEnd;
  
  ierr = DMPlexGetHeightStratum(dmNozzle, 0, &fStart, &fEnd);CHKERRQ(ierr);
  fDAGlevel = fEnd - fStart;
  //ierr = PetscPrintf(PETSC_COMM_SELF, "    fStart = %d \n", fStart);
  //ierr = PetscPrintf(PETSC_COMM_SELF, "    fEnd = %d \n", fEnd);
  
  ierr = DMPlexGetHeightStratum(dmNozzle, 1, &eStart, &eEnd);CHKERRQ(ierr);
  eDAGlevel = eEnd - eStart;
  //ierr = PetscPrintf(PETSC_COMM_SELF, "    eStart = %d \n", eStart);
  //ierr = PetscPrintf(PETSC_COMM_SELF, "    eEnd = %d \n", eEnd);
  
  ierr = DMPlexGetHeightStratum(dmNozzle, 2, &nStart, &nEnd);CHKERRQ(ierr);
  nDAGlevel = nEnd - nStart;
  //ierr = PetscPrintf(PETSC_COMM_SELF, "    nStart = %d \n", nStart);
  //ierr = PetscPrintf(PETSC_COMM_SELF, "    nEnd = %d \n", nEnd);
  
  
  /* Set Label Values to EGADS faces & edges */
  int eCntr = 0;
  for (b = 0; b < nbodies; ++b){
    ego body = bodies[b];
    ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
    
    for (int f = 0; f < Nf; ++f){
      ego face = fobjs[f];
      int fID;
      
      fID = EG_indexBodyTopo(body, face);CHKERRQ(ierr);    // face ID
      
      ierr = EG_getBodyTopos(body, face, EDGE, &Ne, &eobjs);CHKERRQ(ierr);
      
      for (int e = 0; e < Ne; ++e){
        ego edge = eobjs[e];
        int eID;
        
        eID = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);    // edge ID
      
        ierr = EG_getBodyTopos(body, edge, NODE, &Nv, &nobjs);CHKERRQ(ierr);
        
        for (int v = 0; v < Nv; ++v){
          ego vertex = nobjs[v];
          int vID;
          
          vID = EG_indexBodyTopo(body, vertex);CHKERRQ(ierr);    // vertex ID
          
          // Edge Endnodes
          ierr = DMLabelSetValue(vertexLabel, nStart + vID - 1, vID); CHKERRQ(ierr);
        }
        // Edge MidPoint
        ierr = DMLabelSetValue(edgeLabel, nStart + Nvtotal + eID - 1, eID);CHKERRQ(ierr);
      }
      // Face Center Node
      ierr = DMLabelSetValue(faceLabel, nStart + Nvtotal + Netotal + fID - 1, fID);CHKERRQ(ierr); 
    }  
  }
  
  /* Define Cells faceLabels */
  cellCntr = 0;
  for (b = 0; b < nbodies; ++b){
    ego body = bodies[b];
    ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
    
    for (int f = 0; f < Nf; ++f){
      ego face = fobjs[f];
      int fID;
      
      fID = EG_indexBodyTopo(body, face);CHKERRQ(ierr);    // face ID
      
      ierr = EG_getBodyTopos(body, face, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
      
      for (int e = 0; e < Ne; ++e){
        ego edge = eobjs[e];
        int eID = EG_indexBodyTopo(body, edge); CHKERRQ(ierr);  // edge ID
      
        for (int jj = 0; jj < 2; ++jj){
          PetscInt coneSize = 0, coneSizeN = 0;
          PetscInt *cone = NULL, *coneN = NULL, *coneOrient = NULL;
          
          ierr = DMLabelSetValue(faceLabel, cellCntr, fID);CHKERRQ(ierr);
          ierr = DMPlexGetConeSize(dmNozzle, cellCntr, &coneSize); CHKERRQ(ierr);
          ierr = DMPlexGetCone(dmNozzle, cellCntr, &cone); CHKERRQ(ierr);
          ierr = DMPlexGetConeOrientation(dmNozzle, cellCntr, &coneOrient); CHKERRQ(ierr);
          
          ierr = DMLabelSetValue(faceLabel, cone[0], fID);CHKERRQ(ierr);
          ierr = DMLabelSetValue(edgeLabel, cone[1], eID);CHKERRQ(ierr);
          ierr = DMLabelSetValue(faceLabel, cone[2], fID);CHKERRQ(ierr);

          cellCntr = cellCntr + 1;
        }
      }
    }
  }
    
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmNozzle \n");CHKERRQ(ierr);
  ierr = DMView(dmNozzle, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  ierr = DMLabelView(bodyLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMLabelView(faceLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMLabelView(edgeLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMLabelView(vertexLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  //ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
  ierr = DMSetFromOptions(dmNozzle); CHKERRQ(ierr);
  
  // Remove when not using Tetgen
  ierr = DMPlexGenerate(dmNozzle, "tetgen", PETSC_TRUE, &dmMesh); CHKERRQ(ierr);
  
  /* Attached EGADS model to Volumetric Mesh DMPlex */
{
  PetscContainer modelObj;
  ierr = PetscContainerCreate(PETSC_COMM_SELF, &modelObj);CHKERRQ(ierr);
  ierr = PetscContainerSetPointer(modelObj, model);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject) dmMesh, "EGADS Model", (PetscObject) modelObj);CHKERRQ(ierr);
  ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n Attached EGADS Model \n");CHKERRQ(ierr);
}
  
  /* Inflate Mesh to EGADS Geometry */
  //ierr = DMPlexInflateToGeomModel(dmMesh); CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n Inflated dmMesh \n");CHKERRQ(ierr);
  
  /* Output State of DMLabels for dmMesh after Volumetric Mesh generated */
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh \n");CHKERRQ(ierr);
  ierr = DMView(dmMesh, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  ierr = DMCreateLabel(dmMesh, "EGADS Body ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmMesh, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dmMesh, "EGADS Face ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmMesh, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dmMesh, "EGADS Edge ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmMesh, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dmMesh, "EGADS Vertex ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmMesh, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
  
  ierr = DMLabelView(bodyLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMLabelView(faceLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMLabelView(edgeLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMLabelView(vertexLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view3");CHKERRQ(ierr);
  
  /* Refine Volumetric Mesh (dmMesh) */
  // Petsc Refinement
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh Created Trying Refinement \n");CHKERRQ(ierr);
  ierr = DMSetFromOptions(dmMesh);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  // Tetgen Refinement
  //DMRefine(dmMesh,PETSC_COMM_WORLD,&dm);
  
  
  
  /* Print Out Start Location of Mesh Entities */
  ierr = PetscPrintf(PETSC_COMM_SELF, "    -- dmNozzle Entity Start IDs  -- \n");
  ierr = DMPlexGetHeightStratum(dmNozzle, 0, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "    fStart = %d \n", fStart);
  
  ierr = DMPlexGetHeightStratum(dmNozzle, 1, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "    eStart = %d \n", eStart);
  
  ierr = DMPlexGetHeightStratum(dmNozzle, 2, &nStart, &nEnd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "    nStart = %d \n", nStart);

  
  ierr = PetscPrintf(PETSC_COMM_SELF, "    -- dmMesh Entity Start IDs  -- \n");
  ierr = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "    cStart = %d \n", cStart);
  
  ierr = DMPlexGetHeightStratum(dmMesh, 1, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "    fStart = %d \n", fStart);
  
  ierr = DMPlexGetHeightStratum(dmMesh, 2, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "    eStart = %d \n", eStart);
  
  ierr = DMPlexGetHeightStratum(dmMesh, 3, &nStart, &nEnd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "    nStart = %d \n", nStart);
  
  
  
  /* Output Refined dmMesh Information */
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh \n");CHKERRQ(ierr);
  //ierr = DMView(dmMesh, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  //ierr = DMCreateLabel(dmMesh, "EGADS Body ID");CHKERRQ(ierr);
  //ierr = DMGetLabel(dmMesh, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  //ierr = DMCreateLabel(dmMesh, "EGADS Face ID");CHKERRQ(ierr);
  //ierr = DMGetLabel(dmMesh, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  //ierr = DMCreateLabel(dmMesh, "EGADS Edge ID");CHKERRQ(ierr);
  //ierr = DMGetLabel(dmMesh, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  //ierr = DMCreateLabel(dmMesh, "EGADS Vertex ID");CHKERRQ(ierr);
  //ierr = DMGetLabel(dmMesh, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
  
  //ierr = DMLabelView(bodyLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //ierr = DMLabelView(faceLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //ierr = DMLabelView(edgeLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //ierr = DMLabelView(vertexLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  //ierr = DMViewFromOptions(dm, NULL, "-dm_view2");CHKERRQ(ierr);        // Use when Revine dmMesh 1st
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view2");CHKERRQ(ierr);    // Use when Refine dmNozzle 1st
  ierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view");CHKERRQ(ierr);
  
  //ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = DMDestroy(&dmMesh);CHKERRQ(ierr);
  ierr = DMDestroy(&dmNozzle);CHKERRQ(ierr);

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
