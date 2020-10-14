static char help[] = "Heat Equation in 2d and 3d with finite elements.\n\
We solve the heat equation in a rectangular\n\
domain, using a parallel unstructured mesh (DMPLEX) to discretize it.\n\
Contributed by: Julian Andrej <juan@tf.uni-kiel.de>\n\n\n";

#include <petscdmplex.h>
#include <petscds.h>
#include <petscts.h>
#include <egads.h>

/*
  Heat equation:

    du/dt - \Delta u = -1 * dim

  Exact 2D solution:

    u = 2t + x^2 + y^2

    2 - (2 + 2) + 2 = 0

  Exact 3D solution:

    u = 3t + x^2 + y^2 + z^2

    3 - (2 + 2 + 2) + 3 = 0
*/

typedef struct {
  PetscInt          dim;
  PetscBool         simplex;
  char filename[PETSC_MAX_PATH_LEN];
  PetscErrorCode (**exactFuncs)(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx);
} AppCtx;

static PetscErrorCode analytic_temp(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar *u, void *ctx)
{
  PetscInt d;

  *u = dim*time;
  for (d = 0; d < dim; ++d) *u += x[d]*x[d];
  return 0;
}

static void f0_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] = u_t[0] + (PetscScalar) dim;
}

static void f1_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) {
    f1[d] = u_x[d];
  }
}

static void g3_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) {
    g3[d*dim+d] = 1.0;
  }
}

static void g0_temp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[])
{
  g0[0] = u_tShift*1.0;
}

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  options->dim     = 2;
  options->simplex = PETSC_TRUE;
  options->filename[0] = '\0';

  ierr = PetscOptionsBegin(comm, "", "Heat Equation Options w/EGADS", "DMPLEX");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-dim", "The topological mesh dimension", "ex45.c", options->dim, &options->dim, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-simplex", "Simplicial (true) or tensor (false) mesh", "ex45.c", options->simplex, &options->simplex, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-filename", "The EGADSLite file", "ex9.c", options->filename, options->filename, PETSC_MAX_PATH_LEN, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateBCLabel(DM dm, const char name[])
{
  DM             plex;
  DMLabel        label;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMCreateLabel(dm, name);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, name, &label);CHKERRQ(ierr);
  ierr = DMConvert(dm, DMPLEX, &plex);CHKERRQ(ierr);
  ierr = DMPlexMarkBoundaryFaces(plex, 1, label);CHKERRQ(ierr);
  ierr = DMDestroy(&plex);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* -----------------------------------------------------------------------*/
//int main(int argc, char *argv[])
static PetscErrorCode ReadCADGenerateDMPlex(MPI_Comm comm, DM* dmMesh, AppCtx* ctx)
{
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel; // markerLabel;
  PetscInt       cStart, cEnd;
  /* EGADSLite variables */
  ego            context, model, geom, *bodies, *objs, *nobjs, *mobjs, *lobjs, *fobjs, *eobjs, *shobjs;
  ego            *tess;
  int            oclass, mtype, nbodies, *senses, *bsenses, *shsenses, *fsenses, *lsenses, *esenses;
  int            b, maxNumEdges;
  /* PETSc variables */
  DM             dmNozzle; //dm, dmMesh;
  PetscInt       dim = -1, cdim = -1, numCorners = 0, numVertices = 0, numCells = 0, depth = 0;
  PetscInt       numFaces = 0, numEdges = 0, numPoints = 0;
  PetscInt      *cells  = NULL, *coneOrient = NULL, *cones = NULL, *coneSize = NULL;
  PetscReal     *coords = NULL, PARAMACC = 1.0e-4;
  //E MPI_Comm       comm;
  PetscMPIInt    rank;
  //E AppCtx         ctx;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;	
  //E ierr = PetscInitialize(&argc, &argv, NULL, help); if (ierr) return ierr;
  //E comm = PETSC_COMM_WORLD;
  //E ierr = ProcessOptions(comm, &ctx);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) {
    /* Open EGADs file and load EGADs model data */
    ierr = EG_open(&context);CHKERRQ(ierr);
    ierr = EG_loadModel(context, 0, ctx->filename, &model);CHKERRQ(ierr);

    /* test bodyTopo functions */
    ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &bsenses);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", nbodies);CHKERRQ(ierr);
	
	/* Initialize Variable to determine the maximum number of EDGES in any BODY */
	maxNumEdges = 0;

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
	  
	  /* Check to see if current BODY has more EDGEs than the previous BODY */
	  maxNumEdges = PetscMax(maxNumEdges, Ne);

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
              int    peri;
    
              ierr = EG_getRange(objs[e], range, &peri);
              ierr = PetscPrintf(PETSC_COMM_SELF, "              Peri = %d :: Range = %lf, %lf, %lf, %lf \n", peri, range[0], range[1], range[2], range[3]);
    
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
      }   
    }
  }
	
	/* Debug */
	//Bierr = PetscPrintf(PETSC_COMM_SELF, "  maxNumEdges = %d \n", maxNumEdges);CHKERRQ(ierr);
	
	/* Define Matrix edgeIDrelate[nbodies][edgeID] to hold postion in coords[] for cells[] declaration for DMPlex */
	int edgeIDrelate[(const) nbodies][(const) maxNumEdges];
	//Bierr = PetscPrintf(PETSC_COMM_SELF, "  edgeIDrelate[%d][%d] = %d \n", nbodies-1, maxNumEdges-1, edgeIDrelate[nbodies-1][maxNumEdges-1]);CHKERRQ(ierr);
  
	/* Caculate Total Number of Model Entities in the EGADS Model */
	ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &bsenses);CHKERRQ(ierr);
	int Nf=0, Ne=0, Nv=0;
	int bodyIndexStart[nbodies], bodyVertexIndexStart[nbodies], bodyEdgeIndexStart[nbodies], bodyFaceIndexStart[nbodies];
	for (b = 0; b < nbodies; ++b) {
		ego body = bodies[b];
		bodyIndexStart[b] = numFaces + numEdges + numVertices;		// May need to offset by 1 to get right index
		bodyVertexIndexStart[b] = numVertices;
		bodyEdgeIndexStart[b] = numEdges;
		bodyFaceIndexStart[b] = numFaces;
		
		ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
		ierr = EG_getBodyTopos(body, NULL, EDGE,  &Ne, &eobjs);CHKERRQ(ierr);
		ierr = EG_getBodyTopos(body, NULL, NODE,  &Nv, &nobjs);CHKERRQ(ierr);
		
		int Netemp = 0, id;
		for (PetscInt e = 0; e < Ne; ++e) {
			ego edge = eobjs[e];
			ego *topRef, *prev, *next;
			ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
			id   = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);
			if (mtype == DEGENERATE) {
				edgeIDrelate[b][id-1] = -1;
			} else {
				++Netemp;
				edgeIDrelate[b][id-1] = Netemp;	
			}
	//B		ierr = PetscPrintf(PETSC_COMM_SELF, "  edgeIDrelate[%d][%d] = %d \n", b, id-1, edgeIDrelate[b][id-1]);CHKERRQ(ierr);
		}
		numFaces = numFaces + Nf;
		numEdges = numEdges + Netemp;
		numVertices = numVertices + Nv;
		
	//B	ierr = PetscPrintf(PETSC_COMM_SELF, "  bodyIndexStart[%d] = %d \n", b, bodyIndexStart[b]);CHKERRQ(ierr);
	//B	ierr = PetscPrintf(PETSC_COMM_SELF, "  bodyVertexIndexStart[%d] = %d \n", b, bodyVertexIndexStart[b]);CHKERRQ(ierr);
	//B	ierr = PetscPrintf(PETSC_COMM_SELF, "  bodyEdgeIndexStart[%d] = %d \n", b, bodyEdgeIndexStart[b]);CHKERRQ(ierr);
	//B	ierr = PetscPrintf(PETSC_COMM_SELF, "  bodyFaceIndexStart[%d] = %d \n", b, bodyFaceIndexStart[b]);CHKERRQ(ierr);
	}
    
    int Nftotal = numFaces;
    int Netotal = numEdges;
    int Nvtotal = numVertices;
	
	//Bierr = PetscPrintf(PETSC_COMM_SELF, "(Nftotal, Netotal, Nvtotal) = (%d, %d, %d) \n", Nftotal, Netotal, Nvtotal);CHKERRQ(ierr);

    // -------------------------------------------------------------------------------
    // Start to setup DMPlex
    // -------------------------------------------------------------------------------
    dim = 2;
    cdim = 3;
    numCorners = 3;        // Assumes Triangle cells
    //numVertices = Nv + Ne + Nf;		// Original Code - Calculated above
	numPoints = Nftotal + Netotal + Nvtotal;	// New Code
    
    // Determine Number of Cells assuming Triangular Cells (3 nodes)
    numCells = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      ierr = EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs); CHKERRQ(ierr);		// Doesn't handle FACEs with Holes yet
      
      for(int f = 0; f < Nf; ++f){
        ego face = fobjs[f];
        ierr = EG_getBodyTopos(body, face, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
		
		int Netemp = 0;
		for (PetscInt e = 0; e < Ne; ++e) {
			ego edge = eobjs[e];
			ego *topRef, *prev, *next;
			ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
			if (mtype != DEGENERATE) {++Netemp;}
		}
        
        //numCells = numCells + (2 * Ne); 		// Original Code
		numCells = numCells + (2 * Netemp);
      }
    }
	
	/* Debug */
    //Bierr = PetscPrintf(PETSC_COMM_SELF, "(numCells) = (%d) \n", numCells); CHKERRQ(ierr);		// This is the right number
	
	/* Allocate Memory for coords[] and cells[] */
    ierr = PetscMalloc2(numPoints*cdim, &coords, numCells*numCorners, &cells);CHKERRQ(ierr);
    
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
        coords[(bodyVertexIndexStart[b] + id-1)*cdim+0] = limits[0];
        coords[(bodyVertexIndexStart[b] + id-1)*cdim+1] = limits[1];
        coords[(bodyVertexIndexStart[b] + id-1)*cdim+2] = limits[2];
    //B    ierr = PetscPrintf(PETSC_COMM_SELF, "    Node ID = %d \n", bodyVertexIndexStart[b] + (id-1));
    //B    ierr = PetscPrintf(PETSC_COMM_SELF, "      (x,y,z) = (%lf, %lf, %lf) \n \n", coords[(bodyVertexIndexStart[b] + id-1)*cdim+0], coords[(bodyVertexIndexStart[b] + id-1)*cdim+1],coords[(bodyVertexIndexStart[b] + id-1)*cdim+2]);
      }
    }
    
    // Second EDGE midpoints
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      ierr = EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
      
	  int Netemp = 0;
      for (int ii = 0; ii < Ne; ++ii){			// may want to change ii to e???
        ego edge = eobjs[ii];
        double range[2];
        int    *periodic;
        int     id, locate;
        
		/* if EDGE is DEGENERATE, Skip to next EDGE evaluation */
		ego *topRef, *prev, *next;
		ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
		if (mtype == DEGENERATE) {
			continue;
		} /*else {
			++Netemp;
		}*/
		
        id   = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);
		locate = edgeIDrelate[b][id-1];
        ierr = EG_getRange(edge, &range, &periodic); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "EDGE %d range = %lf, %lf \n", ii+1, range[0], range[1]); CHKERRQ(ierr);
        
        double avgt[1];
        avgt[0] = (range[0] + range[1]) / 2.;
        
        double cntrPnt[9];
        ierr = EG_evaluate(edge, avgt, &cntrPnt);
        
        // changed ii to id-1
        //coords[(Nvtotal+bodyEdgeIndexStart[b]+id-1)*cdim+0] = cntrPnt[0];				// Original Code Nv new code Nvtotal
        //coords[(Nvtotal+bodyEdgeIndexStart[b]+id-1)*cdim+1] = cntrPnt[1];
        //coords[(Nvtotal+bodyEdgeIndexStart[b]+id-1)*cdim+2] = cntrPnt[2];
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    Node ID = %d \n", (Nvtotal+bodyEdgeIndexStart[b]+id-1));
        //ierr = PetscPrintf(PETSC_COMM_SELF, "      (x,y,z) = (%lf, %lf, %lf) \n \n", coords[(Nvtotal+bodyEdgeIndexStart[b]+id-1)*cdim+0], coords[(Nvtotal+bodyEdgeIndexStart[b]+id-1)*cdim+1],coords[(Nvtotal+bodyEdgeIndexStart[b]+id-1)*cdim+2]);
      
        //coords[(Nvtotal+bodyEdgeIndexStart[b]+Netemp-1)*cdim+0] = cntrPnt[0];			// 2nd Attempt - limited success
        //coords[(Nvtotal+bodyEdgeIndexStart[b]+Netemp-1)*cdim+1] = cntrPnt[1];
        //coords[(Nvtotal+bodyEdgeIndexStart[b]+Netemp-1)*cdim+2] = cntrPnt[2];
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    Node ID = %d \n", (Nvtotal+bodyEdgeIndexStart[b]+Netemp-1));
        //ierr = PetscPrintf(PETSC_COMM_SELF, "      (x,y,z) = (%lf, %lf, %lf) \n \n", coords[(Nvtotal+bodyEdgeIndexStart[b]+Netemp-1)*cdim+0], 
		//																			 coords[(Nvtotal+bodyEdgeIndexStart[b]+Netemp-1)*cdim+1],
		//																			 coords[(Nvtotal+bodyEdgeIndexStart[b]+Netemp-1)*cdim+2]);
		
		coords[(Nvtotal+bodyEdgeIndexStart[b]+locate-1)*cdim+0] = cntrPnt[0];			// 3rd Attempt - better?
        coords[(Nvtotal+bodyEdgeIndexStart[b]+locate-1)*cdim+1] = cntrPnt[1];
        coords[(Nvtotal+bodyEdgeIndexStart[b]+locate-1)*cdim+2] = cntrPnt[2];
    //B    ierr = PetscPrintf(PETSC_COMM_SELF, "    Node ID = %d \n", (Nvtotal+bodyEdgeIndexStart[b]+locate-1));
    //B    ierr = PetscPrintf(PETSC_COMM_SELF, "      (x,y,z) = (%lf, %lf, %lf) \n \n", coords[(Nvtotal+bodyEdgeIndexStart[b]+locate-1)*cdim+0], 
	//B																				 coords[(Nvtotal+bodyEdgeIndexStart[b]+locate-1)*cdim+1],
	//B																				 coords[(Nvtotal+bodyEdgeIndexStart[b]+locate-1)*cdim+2]);   
	
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
        int    peri;
        int    id;
        
        id   = EG_indexBodyTopo(body, face);CHKERRQ(ierr);
        ierr = EG_getRange(face, &range, &peri); CHKERRQ(ierr);
    //B    ierr = PetscPrintf(PETSC_COMM_SELF, "FACE %d peri = %d range = %lf, %lf, %lf, %lf \n", ii+1, peri, range[0], range[1], range[2], range[3]);CHKERRQ(ierr);
		
        double avgUV[2];
        avgUV[0] = (range[0] + range[1]) / 2.;
        avgUV[1] = (range[2] + range[3]) / 2.;
        
        double cntrPnt[18];
        ierr = EG_evaluate(face, avgUV, &cntrPnt);
        
        // changed ii to id-1
        coords[(Nvtotal+Netotal+bodyFaceIndexStart[b]+id-1)*cdim+0] = cntrPnt[0];
        coords[(Nvtotal+Netotal+bodyFaceIndexStart[b]+id-1)*cdim+1] = cntrPnt[1];
        coords[(Nvtotal+Netotal+bodyFaceIndexStart[b]+id-1)*cdim+2] = cntrPnt[2];
    //B    ierr = PetscPrintf(PETSC_COMM_SELF, "    Node ID = %d \n", (Nvtotal+Netotal+bodyFaceIndexStart[b]+id-1));
    //B    ierr = PetscPrintf(PETSC_COMM_SELF, "      (x,y,z) = (%lf, %lf, %lf) \n \n", coords[(Nvtotal+Netotal+bodyFaceIndexStart[b]+id-1)*cdim+0],
	//B																				 coords[(Nvtotal+Netotal+bodyFaceIndexStart[b]+id-1)*cdim+1],
	//B																				 coords[(Nvtotal+Netotal+bodyFaceIndexStart[b]+id-1)*cdim+2]);
	//B	ierr = PetscPrintf(PETSC_COMM_SELF, "      (u, v) = (%lf, %lf) \n \n", avgUV[0], avgUV[1]);
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
		  
		int midFaceID = Nvtotal + Netotal + bodyFaceIndexStart[b] + fid-1;    // fid-1 was fid	//was bodyVertexIndexStart[b]
		
		int Netemp = 0;		
		for (int e = 0; e < Ne; ++e){
			int locate;
			ego edge = eobjs[e];

			/* if EDGE is DEGENERATE, Skip to next EDGE evaluation */
			ego *topRef, *prev, *next;
			ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
			if (mtype == DEGENERATE) {
				continue;
			} /*else {
				++Netemp;
			}*/
			
			int id   = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);  // ID of current edge
			locate = edgeIDrelate[b][id-1];
			
			//int midPntID = Nvtotal + bodyEdgeIndexStart[b] + Netemp - 1;		// Netemp was id 
			int midPntID = Nvtotal + bodyEdgeIndexStart[b] + locate - 1;		// locate was Netemp
			
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
			cells[cellCntr*numCorners + 1] = bodyVertexIndexStart[b] + startID-1;		//Original Attempt
			cells[cellCntr*numCorners + 2] = midPntID;
			
			cells[cellCntr*numCorners + 3] = midFaceID;
			cells[cellCntr*numCorners + 4] = midPntID;
			cells[cellCntr*numCorners + 5] = bodyVertexIndexStart[b] + endID-1;			// Original Attempt
			
			cellCntr = cellCntr + 2;     
		} 
        //}  
      }
    }
       
    //Build DMPlex  
    ierr = DMPlexCreateFromCellList(PETSC_COMM_WORLD, dim, numCells, numPoints, numCorners, PETSC_TRUE, cells, cdim, coords, &dmNozzle);CHKERRQ(ierr); 
    ierr = PetscFree2(coords, cells);CHKERRQ(ierr);   
    
  /* Attached the EGADS Model to the Surface Plex */
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
  
  /* NOT Necessary - All DMLabels are initialized to -1 */
 // /* Set Label Values - EGADS body*/
 // for (int jj = 0; jj < 3; ++jj){
 //   ierr = DMPlexGetHeightStratum(dmNozzle, jj, &cStart, &cEnd);CHKERRQ(ierr);
 //   for (int ii = cStart; ii < cEnd; ++ii){
 //     ierr = DMLabelSetValue(bodyLabel, ii, 0); CHKERRQ(ierr);    // Need to change this to be more flexible for multi-body parts
 //     ierr = DMLabelSetValue(edgeLabel, ii, -1); CHKERRQ(ierr);
 //     ierr = DMLabelSetValue(faceLabel, ii, -1); CHKERRQ(ierr);
 //     ierr = DMLabelSetValue(vertexLabel, ii, -1); CHKERRQ(ierr);
 //   }
 // }
  
  /* MAY NOT NEED THIS SECTION */
  /* Get Number of DAG Nodes at each level */ 
  int fDAGlevel, eDAGlevel, nDAGlevel;
  int fStart, fEnd, eStart, eEnd, nStart, nEnd;
  
  ierr = DMPlexGetHeightStratum(dmNozzle, 0, &fStart, &fEnd);CHKERRQ(ierr);
  fDAGlevel = fEnd - fStart;
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "    fStart = %d \n", fStart);
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "    fEnd = %d \n", fEnd);
  
  ierr = DMPlexGetHeightStratum(dmNozzle, 1, &eStart, &eEnd);CHKERRQ(ierr);
  eDAGlevel = eEnd - eStart;
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "    eStart = %d \n", eStart);
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "    eEnd = %d \n", eEnd);
  
  ierr = DMPlexGetHeightStratum(dmNozzle, 2, &nStart, &nEnd);CHKERRQ(ierr);
  nDAGlevel = nEnd - nStart;
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "    nStart = %d \n", nStart);
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "    nEnd = %d \n", nEnd);
  
  
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
		
		/* If Edges are DEGENERATE */
		ego *topRef, *prev, *next;
		ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
		if (mtype == DEGENERATE) {
			continue;
		} /*else {
			++Netemp;
		}*/
        
        eID = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);    // edge ID
		int locate = edgeIDrelate[b][eID-1];
      
        ierr = EG_getBodyTopos(body, edge, NODE, &Nv, &nobjs);CHKERRQ(ierr);
        
        for (int v = 0; v < Nv; ++v){
          ego vertex = nobjs[v];
          int vID;
          
          vID = EG_indexBodyTopo(body, vertex);CHKERRQ(ierr);    // vertex ID
          
          // Edge Endnodes
          //ierr = DMLabelSetValue(vertexLabel, nStart + vID - 1, vID); CHKERRQ(ierr);		// Original Code
		  ierr = DMLabelSetValue(bodyLabel, nStart + bodyVertexIndexStart[b] + vID - 1, b); CHKERRQ(ierr);
		  ierr = DMLabelSetValue(vertexLabel, nStart + bodyVertexIndexStart[b] + vID - 1, vID); CHKERRQ(ierr);
        }
        // Edge MidPoint
        //ierr = DMLabelSetValue(edgeLabel, nStart + Nvtotal + eID - 1, eID);CHKERRQ(ierr);		// Original Code
		ierr = DMLabelSetValue(bodyLabel, nStart + Nvtotal + bodyEdgeIndexStart[b] + locate - 1, b);CHKERRQ(ierr);
		ierr = DMLabelSetValue(edgeLabel, nStart + Nvtotal + bodyEdgeIndexStart[b] + locate - 1, eID);CHKERRQ(ierr);
      }
      // Face Center Node
      //ierr = DMLabelSetValue(faceLabel, nStart + Nvtotal + Netotal + fID - 1, fID);CHKERRQ(ierr); 	// Original Code
	  ierr = DMLabelSetValue(bodyLabel, nStart + Nvtotal + Netotal + bodyFaceIndexStart[b] + fID - 1, b);CHKERRQ(ierr);
	  ierr = DMLabelSetValue(faceLabel, nStart + Nvtotal + Netotal + bodyFaceIndexStart[b] + fID - 1, fID);CHKERRQ(ierr);
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
		
		/* If Edges are DEGENERATE */
		ego *topRef, *prev, *next;
		ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
		if (mtype == DEGENERATE) {
			continue;
		} /*else {
			++Netemp;
		}*/
        
        int eID = EG_indexBodyTopo(body, edge); CHKERRQ(ierr);  // edge ID
		//locate = edgeIDrelate[b][eID-1];
		
        for (int jj = 0; jj < 2; ++jj){
          //PetscInt coneSize = 0, coneSizeN = 0;
          PetscInt *cone = NULL, *coneN = NULL, *coneOrient = NULL;
          
		  ierr = DMLabelSetValue(bodyLabel, cellCntr, b); CHKERRQ(ierr);
          ierr = DMLabelSetValue(faceLabel, cellCntr, fID);CHKERRQ(ierr);
          //ierr = DMPlexGetConeSize(dmNozzle, cellCntr, &coneSize); CHKERRQ(ierr);
          ierr = DMPlexGetCone(dmNozzle, cellCntr, &cone); CHKERRQ(ierr);
          //ierr = DMPlexGetConeOrientation(dmNozzle, cellCntr, &coneOrient); CHKERRQ(ierr);
          
		  ierr = DMLabelSetValue(bodyLabel, cone[0], b); CHKERRQ(ierr);
          ierr = DMLabelSetValue(faceLabel, cone[0], fID);CHKERRQ(ierr);
		  
		  ierr = DMLabelSetValue(bodyLabel, cone[1], b); CHKERRQ(ierr);
          ierr = DMLabelSetValue(edgeLabel, cone[1], eID);CHKERRQ(ierr);
		  
		  ierr = DMLabelSetValue(bodyLabel, cone[2], b);CHKERRQ(ierr);
          ierr = DMLabelSetValue(faceLabel, cone[2], fID);CHKERRQ(ierr);

          cellCntr = cellCntr + 1;
        }
      }
    }
  }
    
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmNozzle \n");CHKERRQ(ierr);
  ierr = DMView(dmNozzle, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  ///////
  //E ierr = DMLabelView(bodyLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //E ierr = DMLabelView(faceLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //E ierr = DMLabelView(edgeLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //E ierr = DMLabelView(vertexLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  //ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
  //ierr = DMSetFromOptions(dmNozzle); CHKERRQ(ierr);
  ierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view");CHKERRQ(ierr);
  
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n PRE-TETGEN \n");CHKERRQ(ierr);
  
  // Generate Volumetric Mesh using TETGEN - Remove when not using Tetgen
  ierr = DMPlexGenerate(dmNozzle, "tetgen", PETSC_TRUE, dmMesh); CHKERRQ(ierr);		// was &dmMesh
  
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n POST-TETGEN \n");CHKERRQ(ierr);
  
  /* Attached EGADS model to Volumetric Mesh DMPlex - Don't need. Moved to TETGEN code */
//{
//  PetscContainer modelObj;
//  ierr = PetscContainerCreate(PETSC_COMM_SELF, &modelObj);CHKERRQ(ierr);
//  ierr = PetscContainerSetPointer(modelObj, model);CHKERRQ(ierr);
//  ierr = PetscObjectCompose((PetscObject) dmMesh, "EGADS Model", (PetscObject) modelObj);CHKERRQ(ierr);
//  ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);
//  ierr = PetscPrintf(PETSC_COMM_SELF, "\n Attached EGADS Model \n");CHKERRQ(ierr);
//}
  
  /* Inflate Mesh to EGADS Geometry */
  //ierr = DMPlexInflateToGeomModel(dmMesh); CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n Inflated dmMesh \n");CHKERRQ(ierr);
  
  /* Output State of DMLabels for dmMesh after Volumetric Mesh generated */
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh \n");CHKERRQ(ierr);
  ierr = DMView(*dmMesh, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  //ierr = DMCreateLabel(dmMesh, "EGADS Body ID");CHKERRQ(ierr);
  //E ierr = DMGetLabel(dmMesh, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  //ierr = DMCreateLabel(dmMesh, "EGADS Face ID");CHKERRQ(ierr);
  //E ierr = DMGetLabel(dmMesh, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  //ierr = DMCreateLabel(dmMesh, "EGADS Edge ID");CHKERRQ(ierr);
  //E ierr = DMGetLabel(dmMesh, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  //ierr = DMCreateLabel(dmMesh, "EGADS Vertex ID");CHKERRQ(ierr);
  //E ierr = DMGetLabel(dmMesh, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
  //
  //E ierr = DMLabelView(bodyLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //E ierr = DMLabelView(faceLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //E ierr = DMLabelView(edgeLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //E ierr = DMLabelView(vertexLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  
  ierr = DMViewFromOptions(*dmMesh, NULL, "-dm_view3");CHKERRQ(ierr);
  //Bierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view3");CHKERRQ(ierr);
  
  // Calculate 2D :: Surface Area;  3D :: Volume & Surface Area
  //BsurfArea(dmNozzle);
  //E surfArea(dmMesh);
  
  /* Refine Volumetric Mesh (dmMesh) */
  // Petsc Refinement
  // 1st time
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh Created Trying 1st Refinement \n");CHKERRQ(ierr);
  //E ierr = DMSetFromOptions(dmMesh);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  //Bierr = DMSetFromOptions(dmNozzle);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  //E ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view4");CHKERRQ(ierr);
  //Bierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view4");CHKERRQ(ierr);
  
  // Calculate Surface Area
  //BsurfArea(dmNozzle);
  //E surfArea(dmMesh);
  
  // 2nd Time
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh Created Trying 2nd Refinement \n");CHKERRQ(ierr);
  //E ierr = DMSetFromOptions(dmMesh);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric 
  //Bierr = DMSetFromOptions(dmNozzle);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  //E ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view5");CHKERRQ(ierr);
  //Bierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view5");CHKERRQ(ierr);
  
  // Calculate Surface Area
  //BsurfArea(dmNozzle);
  //E surfArea(dmMesh);
  
  // 3rd Time
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh Created Trying 3rd Refinement \n");CHKERRQ(ierr);
  //E ierr = DMSetFromOptions(dmMesh);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  //Bierr = DMSetFromOptions(dmNozzle);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  //E ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view6");CHKERRQ(ierr);
  //Bierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view6");CHKERRQ(ierr);
  
  // Calculate Surface Area
  //BsurfArea(dmNozzle);
  //E surfArea(dmMesh);

  // 4th Time
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh Created Trying 4th Refinement \n");CHKERRQ(ierr);
  //E ierr = DMSetFromOptions(dmMesh);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  //Bierr = DMSetFromOptions(dmNozzle);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  //E ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view7");CHKERRQ(ierr);
  //Bierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view7");CHKERRQ(ierr);
  
  // Calculate Surface Area
  //BsurfArea(dmNozzle);
  //E surfArea(dmMesh);
  
  // 5th Time
  //E ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh Created Trying 5th Refinement \n");CHKERRQ(ierr);
  //E ierr = DMSetFromOptions(dmMesh);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  //Bierr = DMSetFromOptions(dmNozzle);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  //E ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view8");CHKERRQ(ierr);
  //Bierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view8");CHKERRQ(ierr);
  
  // Calculate Surface Area
  //BsurfArea(dmNozzle);
  //E surfArea(dmMesh);
  
  // Tetgen Refinement - Do Not Use
  //DMRefine(dmMesh,PETSC_COMM_WORLD,&dm);
  
  /* Inflate Mesh to EGADS Geometry */
  //ierr = DMPlexInflateToGeomModel(dmMesh); CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n Inflated Refined dmMesh \n");CHKERRQ(ierr);
   
  
  /* Print Out Start Location of Mesh Entities */
  //Bierr = PetscPrintf(PETSC_COMM_SELF, "\n    -- dmNozzle Entity Start IDs  -- \n");
  //Bierr = DMPlexGetHeightStratum(dmNozzle, 0, &fStart, &fEnd);CHKERRQ(ierr);
  //Bierr = PetscPrintf(PETSC_COMM_SELF, "    fStart = %d \n", fStart);
  
  //Bierr = DMPlexGetHeightStratum(dmNozzle, 1, &eStart, &eEnd);CHKERRQ(ierr);
  //Bierr = PetscPrintf(PETSC_COMM_SELF, "    eStart = %d \n", eStart);
  
  //Bierr = DMPlexGetHeightStratum(dmNozzle, 2, &nStart, &nEnd);CHKERRQ(ierr);
  //Bierr = PetscPrintf(PETSC_COMM_SELF, "    nStart = %d \n", nStart);

  
  //Aierr = PetscPrintf(PETSC_COMM_SELF, "    -- dmMesh Entity Start IDs  -- \n");
  //Aierr = DMPlexGetHeightStratum(dmMesh, 0, &cStart, &cEnd);CHKERRQ(ierr);
  //Aierr = PetscPrintf(PETSC_COMM_SELF, "    cStart = %d \n", cStart);
  
  //Aierr = DMPlexGetHeightStratum(dmMesh, 1, &fStart, &fEnd);CHKERRQ(ierr);
  //Aierr = PetscPrintf(PETSC_COMM_SELF, "    fStart = %d \n", fStart);
  
  //Aierr = DMPlexGetHeightStratum(dmMesh, 2, &eStart, &eEnd);CHKERRQ(ierr);
  //Aierr = PetscPrintf(PETSC_COMM_SELF, "    eStart = %d \n", eStart);
  
  //Aierr = DMPlexGetHeightStratum(dmMesh, 3, &nStart, &nEnd);CHKERRQ(ierr);
  //Aierr = PetscPrintf(PETSC_COMM_SELF, "    nStart = %d \n", nStart);
  
  
  
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
  //ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view2");CHKERRQ(ierr);    // Use when Refine dmNozzle 1st
  //ierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view");CHKERRQ(ierr);
  
  //ierr = DMDestroy(&dm);CHKERRQ(ierr);
  //E ierr = DMDestroy(&dmMesh);CHKERRQ(ierr);
  ierr = DMDestroy(&dmNozzle);CHKERRQ(ierr);		// Strange error

  /* Close EGADSlite file */
  ierr = EG_close(context);CHKERRQ(ierr);
  //Eierr = PetscFinalize();
  //Ereturn ierr;
  ierr = PetscPrintf(PETSC_COMM_SELF, " -- End of DMPlex Generation ----\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/* -----------------------------------------------------------------------*/

static PetscErrorCode CreateMesh(MPI_Comm comm, DM *dm, AppCtx *ctx)
{
  DM             pdm = NULL;
  //const PetscInt dim = ctx->dim;
  PetscInt       dim;
  PetscBool      hasLabel;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  //E ierr = DMPlexCreateBoxMesh(comm, dim, ctx->simplex, NULL, NULL, NULL, NULL, PETSC_TRUE, dm);CHKERRQ(ierr);
  ierr = ReadCADGenerateDMPlex(comm, dm, ctx); CHKERRQ(ierr);  
  ierr = DMGetCoordinateDim(*dm, &dim);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, " ---- Have DM Dim = %d ----\n", dim);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) dm, "Mesh");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, " ---- After PetscObjectSetName() ----\n");CHKERRQ(ierr);
  /* If no boundary marker exists, mark the whole boundary */
  ierr = DMHasLabel(*dm, "marker", &hasLabel);CHKERRQ(ierr);
  if (!hasLabel) {ierr = CreateBCLabel(*dm, "marker");CHKERRQ(ierr);}
  /* Distribute mesh over processes */
  ierr = DMPlexDistribute(*dm, 0, NULL, &pdm);CHKERRQ(ierr);
  if (pdm) {
    ierr = DMDestroy(dm);CHKERRQ(ierr);
    *dm  = pdm;
  }
  //ierr = DMSetFromOptions(*dm);CHKERRQ(ierr);
  ierr = DMViewFromOptions(*dm, NULL, "-dm_view");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, " ---- End of CreateMesh() ----\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupProblem(DM dm, AppCtx *ctx)
{
  PetscDS        prob;
  const PetscInt id = 1;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);
  ierr = PetscDSSetResidual(prob, 0, f0_temp, f1_temp);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob, 0, 0, g0_temp, NULL, NULL, g3_temp);CHKERRQ(ierr);
  ctx->exactFuncs[0] = analytic_temp;
  ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "wall", "marker", 0, 0, NULL, (void (*)(void)) ctx->exactFuncs[0], 1, &id, ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupDiscretization(DM dm, AppCtx* ctx)
{
  DM             cdm = dm;
  const PetscInt dim = ctx->dim;
  PetscFE        fe;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* Create finite element */
  ierr = PetscFECreateDefault(PetscObjectComm((PetscObject) dm), dim, 1, ctx->simplex, "temp_", -1, &fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe, "temperature");CHKERRQ(ierr);
  /* Set discretization and boundary conditions for each mesh */
  ierr = DMSetField(dm, 0, NULL, (PetscObject) fe);CHKERRQ(ierr);
  ierr = DMCreateDS(dm);CHKERRQ(ierr);
  ierr = SetupProblem(dm, ctx);CHKERRQ(ierr);
  while (cdm) {
    PetscBool hasLabel;

    ierr = DMHasLabel(cdm, "marker", &hasLabel);CHKERRQ(ierr);
    if (!hasLabel) {ierr = CreateBCLabel(cdm, "marker");CHKERRQ(ierr);}
    ierr = DMCopyDisc(dm, cdm);CHKERRQ(ierr);
    ierr = DMGetCoarseDM(cdm, &cdm);CHKERRQ(ierr);
  }
  ierr = PetscFEDestroy(&fe);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
  AppCtx         ctx;
  DM             dm;
  TS             ts;
  Vec            u, r;
  PetscReal      t       = 0.0;
  PetscReal      L2error = 0.0;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help);if (ierr) return ierr;
  ierr = ProcessOptions(PETSC_COMM_WORLD, &ctx);CHKERRQ(ierr);
  ierr = CreateMesh(PETSC_COMM_WORLD, &dm, &ctx);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(dm, &ctx);CHKERRQ(ierr);
  ierr = PetscMalloc1(1, &ctx.exactFuncs);CHKERRQ(ierr);
  ierr = SetupDiscretization(dm, &ctx);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm, &u);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) u, "temperature");CHKERRQ(ierr);
  ierr = VecDuplicate(u, &r);CHKERRQ(ierr);

  ierr = TSCreate(PETSC_COMM_WORLD, &ts);CHKERRQ(ierr);
  ierr = TSSetDM(ts, dm);CHKERRQ(ierr);
  ierr = DMTSSetBoundaryLocal(dm, DMPlexTSComputeBoundary, &ctx);CHKERRQ(ierr);
  ierr = DMTSSetIFunctionLocal(dm, DMPlexTSComputeIFunctionFEM, &ctx);CHKERRQ(ierr);
  ierr = DMTSSetIJacobianLocal(dm, DMPlexTSComputeIJacobianFEM, &ctx);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  ierr = DMProjectFunction(dm, t, ctx.exactFuncs, NULL, INSERT_ALL_VALUES, u);CHKERRQ(ierr);
  ierr = TSSolve(ts, u);CHKERRQ(ierr);

  ierr = TSGetTime(ts, &t);CHKERRQ(ierr);
  ierr = DMComputeL2Diff(dm, t, ctx.exactFuncs, NULL, u, &L2error);CHKERRQ(ierr);
  if (L2error < 1.0e-11) {ierr = PetscPrintf(PETSC_COMM_WORLD, "L_2 Error: < 1.0e-11\n");CHKERRQ(ierr);}
  else                   {ierr = PetscPrintf(PETSC_COMM_WORLD, "L_2 Error: %g\n", (double)L2error);CHKERRQ(ierr);}
  ierr = VecViewFromOptions(u, NULL, "-sol_vec_view");CHKERRQ(ierr);

  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  //ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFree(ctx.exactFuncs);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

  # Full solves
  test:
    suffix: 2d_p1_r1
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 1 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p1_r3
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p2_r1
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 1 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_p2_r3
    requires: triangle
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dm_refine 3 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q1_r1
    requires: !single
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 1 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q1_r3
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q2_r1
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 1 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 2d_q2_r3
    requires: !single
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -simplex 0 -dm_refine 3 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p1_r1
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 1 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p1_r2
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 2 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p2_r1
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 1 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_p2_r2
    requires: ctetgen
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -dm_refine 2 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q1_r1
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 1 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q1_r2
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 2 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q2_r1
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 1 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor
  test:
    suffix: 3d_q2_r2
    filter: sed -e "s~ATOL~RTOL~g" -e "s~ABS~RELATIVE~g"
    args: -dim 3 -simplex 0 -dm_refine 2 -temp_petscspace_degree 2 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor

TEST*/
