#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/

#include <tetgen.h>

#ifdef PETSC_HAVE_EGADS
#include <egads.h>
#endif

/* This is to fix the tetrahedron orientation from TetGen */
static PetscErrorCode DMPlexInvertCells_Internal(PetscInt dim, PetscInt numCells, PetscInt numCorners, int cells[])
{
  PetscInt       bound = numCells*numCorners, coff;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (coff = 0; coff < bound; coff += numCorners) {
    ierr = DMPlexInvertCell(dim, numCorners, &cells[coff]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PETSC_EXTERN PetscErrorCode DMPlexGenerate_Tetgen(DM boundary, PetscBool interpolate, DM *dm)
{
  MPI_Comm       comm;
  DM_Plex       *mesh      = (DM_Plex *) boundary->data;
  const PetscInt dim       = 3;
  ::tetgenio     in;
  ::tetgenio     out;
  DMLabel        label;
  PetscInt       vStart, vEnd, v, fStart, fEnd, f, eStart, eEnd, e, eConeSize, c, cStart, cEnd;
  const PetscInt *eCone = NULL;
  PetscMPIInt    rank;
  PetscErrorCode ierr;
  
#ifdef PETSC_HAVE_EGADS
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
  const char    *labelName;
  const PetscInt vertexOffset = 100000000, edgeOffset = 200000000, faceOffset = 300000000, bodyOffset = 1000000;
  PetscInt       numLabels;
#else
  const char    *labelName = "marker";
#endif


  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)boundary,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  
#ifdef PETSC_HAVE_EGADS
  ierr = DMGetNumLabels(boundary, &numLabels); CHKERRQ(ierr);
#else
  ierr = DMGetLabel(boundary, labelName, &label);CHKERRQ(ierr);
#endif
  
  /* Work with Points */
  ierr = DMPlexGetDepthStratum(boundary, 0, &vStart, &vEnd);CHKERRQ(ierr);
  
  in.numberofpoints = vEnd - vStart;
  if (in.numberofpoints > 0) {
    PetscSection coordSection;
    Vec          coordinates;
    PetscScalar *array;

    in.pointlist       = new double[in.numberofpoints*dim];
    in.pointmarkerlist = new int[in.numberofpoints];

    ierr = DMGetCoordinatesLocal(boundary, &coordinates);CHKERRQ(ierr);
    ierr = DMGetCoordinateSection(boundary, &coordSection);CHKERRQ(ierr);
    ierr = VecGetArray(coordinates, &array);CHKERRQ(ierr);
    for (v = vStart; v < vEnd; ++v) {
      const PetscInt idx = v - vStart;
      PetscInt       off, d;
	  
      ierr = PetscSectionGetOffset(coordSection, v, &off);CHKERRQ(ierr);
      for (d = 0; d < dim; ++d) in.pointlist[idx*dim + d] = PetscRealPart(array[off+d]);
	  
#ifdef PETSC_HAVE_EGADS
	  PetscInt	val, bodyOffsetVal;
	  ierr = DMGetLabel(boundary, "EGADS Body ID", &label); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(label, v, &val); CHKERRQ(ierr);
	  bodyOffsetVal = val * bodyOffset;
	  
      for (int ii = 0; ii < numLabels; ++ii){
        ierr = DMGetLabelName(boundary, ii, &labelName); CHKERRQ(ierr);
		
		//PetscInt val;
        ierr = DMGetLabel(boundary, labelName, &label); CHKERRQ(ierr);
		ierr = DMLabelGetValue(label, v, &val); CHKERRQ(ierr);
		
		if (val >= 0) {
			if (strcmp(labelName, "EGADS Vertex ID") == 0) { in.pointmarkerlist[idx] = val + vertexOffset + bodyOffsetVal; }		// Store EGADS Vertex date
			if (strcmp(labelName, "EGADS Edge ID") == 0) { in.pointmarkerlist[idx] = val + edgeOffset + bodyOffsetVal; }			// Store EGADS Edge data 
			if (strcmp(labelName, "EGADS Face ID") == 0) { in.pointmarkerlist[idx] = val + faceOffset + bodyOffsetVal; }			// Store EGADS Face data 
		}
	  }
#else
      if (label) {
        PetscInt val;

        ierr = DMLabelGetValue(label, v, &val);CHKERRQ(ierr);
        in.pointmarkerlist[idx] = (int) val;
      }
#endif
    }
    ierr = VecRestoreArray(coordinates, &array);CHKERRQ(ierr);
  }
  
  /* Work with Edges */
  ierr = DMPlexGetHeightStratum(boundary, 1, &eStart, &eEnd);CHKERRQ(ierr);
  
  in.numberofedges = eEnd - eStart;
  if (in.numberofedges > 0) {
    in.edgelist = new int[in.numberofedges * 2];
	in.edgemarkerlist = new int[in.numberofedges];
	for (e = eStart; e < eEnd; ++e) {
	  const PetscInt idx     = e - eStart;
	  
	  //* Get ID's of Vertices at each end of the Edge
      ierr = DMPlexGetConeSize(boundary, e, &eConeSize); CHKERRQ(ierr);
      ierr = DMPlexGetCone(boundary, e, &eCone); CHKERRQ(ierr);
	  in.edgelist[idx*2] = eCone[0] - vStart;
	  in.edgelist[idx*2 + 1] = eCone[1] - vStart;

#ifdef PETSC_HAVE_EGADS
	  PetscInt	val, bodyOffsetVal;
	  ierr = DMGetLabel(boundary, "EGADS Body ID", &label); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(label, e, &val); CHKERRQ(ierr);
	  bodyOffsetVal = val * bodyOffset;
	  
      for (int ii = 0; ii < numLabels; ++ii){
        ierr = DMGetLabelName(boundary, ii, &labelName); CHKERRQ(ierr);
		
        ierr = DMGetLabel(boundary, labelName, &label); CHKERRQ(ierr);
		ierr = DMLabelGetValue(label, e, &val); CHKERRQ(ierr);
		
		if (val >= 0) {
			if (strcmp(labelName, "EGADS Vertex ID")== 0) { in.edgemarkerlist[idx] = val + vertexOffset + bodyOffsetVal; }		// Store EGADS Vertex date
			if (strcmp(labelName, "EGADS Edge ID") == 0) { in.edgemarkerlist[idx] = val + edgeOffset + bodyOffsetVal; }			// Store EGADS Edge data 
			if (strcmp(labelName, "EGADS Face ID") == 0) { in.edgemarkerlist[idx] = val + faceOffset + bodyOffsetVal; }			// Store EGADS Face data 
		}
	  }
#else
      if (label) {
        PetscInt val;

        ierr = DMLabelGetValue(label, e, &val);CHKERRQ(ierr);
        in.edgemarkerlist[idx] = (int) val;
      }
#endif
	}
  }
 
  /* Work with Faces */
  ierr  = DMPlexGetHeightStratum(boundary, 0, &fStart, &fEnd);CHKERRQ(ierr);

  in.numberoffacets = fEnd - fStart;
  if (in.numberoffacets > 0) {
    in.facetlist       = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (f = fStart; f < fEnd; ++f) {
      const PetscInt idx     = f - fStart;
      PetscInt      *points = NULL, numPoints, p, numVertices = 0, v;

      in.facetlist[idx].numberofpolygons = 1;
      in.facetlist[idx].polygonlist      = new tetgenio::polygon[in.facetlist[idx].numberofpolygons];
      in.facetlist[idx].numberofholes    = 0;
      in.facetlist[idx].holelist         = NULL;

      ierr = DMPlexGetTransitiveClosure(boundary, f, PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
      for (p = 0; p < numPoints*2; p += 2) {
        const PetscInt point = points[p];
        if ((point >= vStart) && (point < vEnd)) points[numVertices++] = point;
      }

      tetgenio::polygon *poly = in.facetlist[idx].polygonlist;
      poly->numberofvertices = numVertices;
      poly->vertexlist       = new int[poly->numberofvertices];
      for (v = 0; v < numVertices; ++v) {
        const PetscInt vIdx = points[v] - vStart;
        poly->vertexlist[v] = vIdx;
      }
	  
#ifdef PETSC_HAVE_EGADS
	  PetscInt	val, bodyOffsetVal;
	  ierr = DMGetLabel(boundary, "EGADS Body ID", &label); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(label, f, &val); CHKERRQ(ierr);
	  bodyOffsetVal = val * bodyOffset;
	  
      for (int ii = 0; ii < numLabels; ++ii){
        ierr = DMGetLabelName(boundary, ii, &labelName); CHKERRQ(ierr);
		
        ierr = DMGetLabel(boundary, labelName, &label); CHKERRQ(ierr);
		ierr = DMLabelGetValue(label, f, &val); CHKERRQ(ierr);
		
		if (val >= 0) {
			if (strcmp(labelName, "EGADS Vertex ID")== 0) { in.facetmarkerlist[idx] = val + vertexOffset + bodyOffsetVal; }		// Store EGADS Vertex date
			if (strcmp(labelName, "EGADS Edge ID") == 0) { in.facetmarkerlist[idx] = val + edgeOffset + bodyOffsetVal; }		// Store EGADS Edge data 
			if (strcmp(labelName, "EGADS Face ID") == 0) { in.facetmarkerlist[idx] = val + faceOffset + bodyOffsetVal; }		// Store EGADS Face data 
		}
	  }
#else
      if (label) {
        PetscInt val;

        ierr = DMLabelGetValue(label, f, &val);CHKERRQ(ierr);
        in.facetmarkerlist[idx] = (int) val;
      }
#endif

      ierr = DMPlexRestoreTransitiveClosure(boundary, f, PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
    }
  }
  
  /* Generate Mesh using Tetgen */
  if (!rank) {
    char args[32];
	
    /* Take away 'Q' for verbose output */ 
#ifdef PETSC_HAVE_EGADS
    /* Add Y to preserve Surface Mesh for EGADS */
    ierr = PetscStrcpy(args, "pYqezQ");CHKERRQ(ierr);
#else
	ierr = PetscStrcpy(args, "pqezQ");CHKERRQ(ierr);
#endif
    if (mesh->tetgenOpts) {::tetrahedralize(mesh->tetgenOpts, &in, &out);}
    else                  {::tetrahedralize(args, &in, &out);}
  }
  {
    const PetscInt numCorners  = 4;
    const PetscInt numCells    = out.numberoftetrahedra;
    const PetscInt numVertices = out.numberofpoints;
    const double   *meshCoords = out.pointlist;
    int            *cells      = out.tetrahedronlist;

#ifdef PETSC_HAVE_EGADS
	PetscContainer modelObj;
	ego           *bodies;
	ego            model, geom, body;
	int            Nb, oclass, mtype, *senses;
#else
	DMLabel         glabel     = NULL;
#endif

    ierr = DMPlexInvertCells_Internal(dim, numCells, numCorners, cells);CHKERRQ(ierr);
    ierr = DMPlexCreateFromCellList(comm, dim, numCells, numVertices, numCorners, interpolate, cells, dim, meshCoords, dm);CHKERRQ(ierr);

#ifdef PETSC_HAVE_EGADS
	ierr = DMCreateLabel(*dm, "EGADS Body ID");CHKERRQ(ierr);
	ierr = DMGetLabel(*dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
	ierr = DMCreateLabel(*dm, "EGADS Face ID");CHKERRQ(ierr);
	ierr = DMGetLabel(*dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
	ierr = DMCreateLabel(*dm, "EGADS Edge ID");CHKERRQ(ierr);
	ierr = DMGetLabel(*dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
	ierr = DMCreateLabel(*dm, "EGADS Vertex ID");CHKERRQ(ierr);
	ierr = DMGetLabel(*dm, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
	
	/* Get Attached EGADS Model from Original DMPlex*/
	ierr = PetscObjectQuery((PetscObject) boundary, "EGADS Model", (PetscObject *) &modelObj);CHKERRQ(ierr);
	ierr = PetscContainerGetPointer(modelObj, (void **) &model);CHKERRQ(ierr);
	ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
	
	/* Transfer EGADS Model to Volumetric Mesh */
	ierr = PetscObjectCompose((PetscObject) *dm, "EGADS Model", (PetscObject) modelObj);CHKERRQ(ierr);
	ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);
#else
    if (label) {ierr = DMCreateLabel(*dm, labelName);
	ierr = DMGetLabel(*dm, labelName, &glabel);}
#endif
	
    /* -- Set labels -- */
	/* Vertices Labels */		//rewrite to utilize newly created DMPlex from Tetgen. I think it will stabilize the mesh generation labeling
    for (v = 0; v < numVertices; ++v) {
#ifdef PETSC_HAVE_EGADS
	  PetscInt val, bodyVal, bodyValOffset;
	  
	  val = out.pointmarkerlist[v];
	  if (val >= vertexOffset && val < edgeOffset) {
		val = val - vertexOffset;
		bodyVal = val / bodyOffset;		// Integer Division
		bodyValOffset = bodyVal * bodyOffset;
		ierr = DMLabelSetValue(vertexLabel, v+numCells, out.pointmarkerlist[v] - vertexOffset - bodyValOffset); CHKERRQ(ierr);
		ierr = DMLabelSetValue(bodyLabel, v+numCells, bodyVal); CHKERRQ(ierr);
	  }
	  if (val >= edgeOffset && val < faceOffset) {
		val = val - edgeOffset;
		bodyVal = val / bodyOffset;		// Integer Division
		bodyValOffset = bodyVal * bodyOffset;
		ierr = DMLabelSetValue(edgeLabel, v+numCells, out.pointmarkerlist[v] - edgeOffset - bodyValOffset); CHKERRQ(ierr);
		ierr = DMLabelSetValue(bodyLabel, v+numCells, bodyVal); CHKERRQ(ierr);
	  }
	  if (val >= faceOffset) {
		val = val - faceOffset;
		bodyVal = val / bodyOffset;		// Integer Division
		bodyValOffset = bodyVal * bodyOffset;
	    ierr = DMLabelSetValue(faceLabel, v+numCells, out.pointmarkerlist[v] - faceOffset - bodyValOffset); CHKERRQ(ierr);
		ierr = DMLabelSetValue(bodyLabel, v+numCells, bodyVal); CHKERRQ(ierr);
	  }
#else
      if (out.pointmarkerlist[v]) {
        if (glabel) {ierr = DMLabelSetValue(glabel, v+numCells, out.pointmarkerlist[v]);CHKERRQ(ierr);}
      }
#endif
    }
	
	
    if (interpolate) {
      /* This check is never actually executed for ctetgen (which never returns edgemarkers) and seems to be broken for
       * tetgen */ /* Original Code didn't load edgelist[] or edgemarkerlist[] and is probably why this didn't work previously */
	  /* Set Edge Labels */
	  ierr = DMPlexGetHeightStratum(*dm, 2, &eStart, &eEnd);CHKERRQ(ierr);
      for (e = 0; e < out.numberofedges; e++) {
        if (out.edgemarkerlist[e]) {
          const PetscInt  vertices[2] = {out.edgelist[e*2+0]+numCells, out.edgelist[e*2+1]+numCells};
          const PetscInt *edges;
          PetscInt        numEdges;

          ierr = DMPlexGetJoin(*dm, 2, vertices, &numEdges, &edges);CHKERRQ(ierr);
          if (numEdges != 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Two vertices must cover only one edge, not %D", numEdges);
		  
#ifdef PETSC_HAVE_EGADS
		  PetscInt val, bodyVal, bodyValOffset;
	  
		  val = out.edgemarkerlist[e];
		  if (val >= vertexOffset && val < edgeOffset) {
			val = val - vertexOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
		    ierr = DMLabelSetValue(vertexLabel, edges[0], out.edgemarkerlist[e] - vertexOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, edges[0], bodyVal); CHKERRQ(ierr);
		  }
		  if (val >= edgeOffset && val < faceOffset) {
			val = val - edgeOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
	        ierr = DMLabelSetValue(edgeLabel, edges[0], out.edgemarkerlist[e] - edgeOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, edges[0], bodyVal); CHKERRQ(ierr);
		  }
		  if (val >= faceOffset) {
			val = val - faceOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
		    ierr = DMLabelSetValue(faceLabel, edges[0], out.edgemarkerlist[e] - faceOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, edges[0], bodyVal); CHKERRQ(ierr);
		  }
#else
          if (glabel) {ierr = DMLabelSetValue(glabel, edges[0], out.edgemarkerlist[e]);CHKERRQ(ierr);}
#endif

          ierr = DMPlexRestoreJoin(*dm, 2, vertices, &numEdges, &edges);CHKERRQ(ierr);
        }
      }
	  
	  /* Set Face Labels */
	  ierr = DMPlexGetHeightStratum(*dm, 1, &fStart, &fEnd);CHKERRQ(ierr);
      for (f = 0; f < out.numberoftrifaces; f++) {
        if (out.trifacemarkerlist[f]) {
          const PetscInt  vertices[3] = {out.trifacelist[f*3+0]+numCells, out.trifacelist[f*3+1]+numCells, out.trifacelist[f*3+2]+numCells};
          const PetscInt *faces;
          PetscInt        numFaces;

          ierr = DMPlexGetFullJoin(*dm, 3, vertices, &numFaces, &faces);CHKERRQ(ierr);
          if (numFaces != 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Three vertices must cover only one face, not %D", numFaces);
		  
#ifdef PETSC_HAVE_EGADS
		  PetscInt val, bodyVal, bodyValOffset;
	  
		  val = out.trifacemarkerlist[f];
		  if (val >= vertexOffset && val < edgeOffset) {
			val = val - vertexOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
			ierr = DMLabelSetValue(vertexLabel, faces[0], out.trifacemarkerlist[f] - vertexOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, faces[0], bodyVal); CHKERRQ(ierr);
		  }
		  if (val >= edgeOffset && val < faceOffset) {
			val = val - edgeOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
			ierr = DMLabelSetValue(edgeLabel, faces[0], out.trifacemarkerlist[f] - edgeOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, faces[0], bodyVal); CHKERRQ(ierr);
		  }
		  if (val >= faceOffset) {
			val = val - faceOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
		    ierr = DMLabelSetValue(faceLabel, faces[0], out.trifacemarkerlist[f] - faceOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, faces[0], bodyVal); CHKERRQ(ierr);
		  }
#else
          if (glabel) {ierr = DMLabelSetValue(glabel, faces[0], out.trifacemarkerlist[f]);CHKERRQ(ierr);}
#endif
          ierr = DMPlexRestoreJoin(*dm, 3, vertices, &numFaces, &faces);CHKERRQ(ierr);
        }
      }
    }
	
	/* Set Cell Labels */
#ifdef PETSC_HAVE_EGADS
	/* NEED TO UPDATE TO DEAL WITH WHAT BODY THE CELL IS LOCATED IN */
	/*   1 - Cycle Through Cells 									*/
	/*   2 - Calculate the centroid (x,y,z) of tetrahedrals			*/
	/*   3 - Perfrom EG_inTopology() function to see which body		*/
	/*   4 - Perform transitive closure								*/
	/*   5 - if attached enitites don't have a bodyID assigned		*/
	/*       assign the current bodyID								*/
	
	ierr = DMPlexGetHeightStratum(*dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
	ierr = DMPlexGetHeightStratum(*dm, 1, &fStart, &fEnd);CHKERRQ(ierr);
	ierr = DMPlexGetHeightStratum(*dm, 2, &eStart, &eEnd);CHKERRQ(ierr);

	for ( c = cStart; c < cEnd; c++) {
		PetscInt 	   val=-1, eVal, fVal;
		PetscReal	   vol;
		PetscReal	   centroid[dim], normal[dim];
		PetscInt      *points = NULL, numPoints;
		
		ierr = DMPlexGetTransitiveClosure(*dm, c, PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);		// Get all DMPlex points associated with the current cell
		ierr = DMPlexComputeCellGeometryFVM(*dm, c, &vol, centroid, normal); CHKERRQ(ierr);		
		
		/* Deterimine what body the cell's centroid is located in */
		for (int ii = 0; ii < Nb; ++ii){
			body = bodies[ii];
			int cgCheck = EG_inTopology(body, centroid); 
			if (cgCheck == EGADS_SUCCESS) {
				val = ii;
				break;
			}
		}

		if (val >= 0) {
			ierr = DMLabelSetValue(bodyLabel, c, val); CHKERRQ(ierr);
			
			
			for (int ii = 0; ii < numPoints; ++ii){
				if (points[ii] >= eStart && points[ii] < eEnd) {
					ierr = DMLabelGetValue(bodyLabel, points[ii], &eVal); CHKERRQ(ierr);		// We don't want to overwrite the work done above
					if (eVal < 0) {ierr = DMLabelSetValue(bodyLabel, points[ii], val); CHKERRQ(ierr);}
				} else {
					// Do Nothing - Skip
				}
				
				if (points[ii] >= fStart && points[ii] < fEnd) {
					ierr = DMLabelGetValue(bodyLabel, points[ii], &fVal); CHKERRQ(ierr);		// We don't want to overwrite the work done above
					if (fVal < 0) {ierr = DMLabelSetValue(bodyLabel, points[ii], val); CHKERRQ(ierr);}
				} else {
					// Do Nothing - Skip
				}
			}
			
			/*for (int ii = 0; ii < numPoints; ++ii){
				ierr = DMLabelSetValue(bodyLabel, points[ii], val); CHKERRQ(ierr);
			}*/
		} else {
			// Do Nothing
		}		
		ierr = DMPlexRestoreTransitiveClosure(*dm, c, PETSC_TRUE, &numPoints, &points); CHKERRQ(ierr);
	}
#else
	// Do Nothing
#endif
	
    ierr = DMPlexSetRefinementUniform(*dm, PETSC_FALSE);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}












PETSC_EXTERN PetscErrorCode DMPlexRefine_Tetgen(DM dm, double *maxVolumes, DM *dmRefined)
{  
  MPI_Comm       comm;
  const PetscInt dim       = 3;
  ::tetgenio     in;
  ::tetgenio     out;
  DMLabel        label;
  PetscInt       vStart, vEnd, v, fStart, fEnd, f, eStart, eEnd, e, eConeSize, cStart, cEnd, c;
  PetscInt       depth, depthGlobal;
  const PetscInt *eCone = NULL;
  PetscMPIInt    rank;
  PetscErrorCode ierr;
  
#ifdef PETSC_HAVE_EGADS
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
  const char    *labelName;
  const PetscInt vertexOffset = 100000000, edgeOffset = 200000000, faceOffset = 300000000, bodyOffset = 1000000;
  PetscInt       numLabels;
#else
  const char    *labelName = "marker";
#endif

  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject)dm,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = DMPlexGetDepth(dm, &depth);CHKERRQ(ierr);
  ierr = MPIU_Allreduce(&depth, &depthGlobal, 1, MPIU_INT, MPI_MAX, comm);CHKERRQ(ierr);
  
#ifdef PETSC_HAVE_EGADS
  ierr = DMGetNumLabels(dm, &numLabels); CHKERRQ(ierr);
#else
  ierr = DMGetLabel(dm, labelName, &label);CHKERRQ(ierr);
#endif  
  
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);

  /* Work with Point */
  in.numberofpoints = vEnd - vStart;
  if (in.numberofpoints > 0) {
    PetscSection coordSection;
    Vec          coordinates;
    PetscScalar *array;

    in.pointlist       = new double[in.numberofpoints*dim];
    in.pointmarkerlist = new int[in.numberofpoints];

    ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
    ierr = DMGetCoordinateSection(dm, &coordSection);CHKERRQ(ierr);
    ierr = VecGetArray(coordinates, &array);CHKERRQ(ierr);
    for (v = vStart; v < vEnd; ++v) {
      const PetscInt idx = v - vStart;
      PetscInt       off, d;

      ierr = PetscSectionGetOffset(coordSection, v, &off);CHKERRQ(ierr);
      for (d = 0; d < dim; ++d) in.pointlist[idx*dim + d] = PetscRealPart(array[off+d]);
	  
#ifdef PETSC_HAVE_EGADS
	  PetscInt	val, bodyOffsetVal;
	  ierr = DMGetLabel(dm, "EGADS Body ID", &label); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(label, v, &val); CHKERRQ(ierr);
	  bodyOffsetVal = val * bodyOffset;
	  
      for (int ii = 0; ii < numLabels; ++ii){
        ierr = DMGetLabelName(dm, ii, &labelName); CHKERRQ(ierr);
		
        ierr = DMGetLabel(dm, labelName, &label); CHKERRQ(ierr);
		ierr = DMLabelGetValue(label, v, &val); CHKERRQ(ierr);
		
		if (val >= 0) {
			if (strcmp(labelName, "EGADS Vertex ID") == 0) { in.pointmarkerlist[idx] = val + vertexOffset + bodyOffsetVal; }		// Store EGADS Vertex date
			if (strcmp(labelName, "EGADS Edge ID") == 0) { in.pointmarkerlist[idx] = val + edgeOffset + bodyOffsetVal; }			// Store EGADS Edge data 
			if (strcmp(labelName, "EGADS Face ID") == 0) { in.pointmarkerlist[idx] = val + faceOffset + bodyOffsetVal; }			// Store EGADS Face data 
		}
	  }	  
	  
#else
      if (label) {
        PetscInt val;

        ierr = DMLabelGetValue(label, v, &val);CHKERRQ(ierr);
        in.pointmarkerlist[idx] = (int) val;
      }
#endif
    }
    ierr = VecRestoreArray(coordinates, &array);CHKERRQ(ierr);
  }
  
  /* Work with Edges */
  ierr = DMPlexGetHeightStratum(dm, 2, &eStart, &eEnd);CHKERRQ(ierr);
  
  in.numberofedges = eEnd - eStart;
  if (in.numberofedges > 0) {
    in.edgelist = new int[in.numberofedges * 2];
	in.edgemarkerlist = new int[in.numberofedges];
	for (e = eStart; e < eEnd; ++e) {
	  const PetscInt idx     = e - eStart;
	  
	  //* Get ID's of Vertices at each end of the Edge
      ierr = DMPlexGetConeSize(dm, e, &eConeSize); CHKERRQ(ierr);
      ierr = DMPlexGetCone(dm, e, &eCone); CHKERRQ(ierr);
	  in.edgelist[idx*2] = eCone[0] - vStart;
	  in.edgelist[idx*2 + 1] = eCone[1] - vStart;

#ifdef PETSC_HAVE_EGADS
	  PetscInt	val, bodyOffsetVal;
	  ierr = DMGetLabel(dm, "EGADS Body ID", &label); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(label, e, &val); CHKERRQ(ierr);
	  bodyOffsetVal = val * bodyOffset;
	  
      for (int ii = 0; ii < numLabels; ++ii){
        ierr = DMGetLabelName(dm, ii, &labelName); CHKERRQ(ierr);
		
        ierr = DMGetLabel(dm, labelName, &label); CHKERRQ(ierr);
		ierr = DMLabelGetValue(label, e, &val); CHKERRQ(ierr);
		
		if (val >= 0) {
			if (strcmp(labelName, "EGADS Vertex ID")== 0) { in.edgemarkerlist[idx] = val + vertexOffset + bodyOffsetVal; }		// Store EGADS Vertex date
			if (strcmp(labelName, "EGADS Edge ID") == 0) { in.edgemarkerlist[idx] = val + edgeOffset + bodyOffsetVal; }			// Store EGADS Edge data 
			if (strcmp(labelName, "EGADS Face ID") == 0) { in.edgemarkerlist[idx] = val + faceOffset + bodyOffsetVal; }			// Store EGADS Face data 
		}
	  }
#else
      if (label) {
        PetscInt val;

        ierr = DMLabelGetValue(label, e, &val);CHKERRQ(ierr);
        in.edgemarkerlist[idx] = (int) val;
      }
#endif
	}
  }
 
  /* Work with Faces */
  ierr  = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);CHKERRQ(ierr);

  in.numberoffacets = fEnd - fStart;
  if (in.numberoffacets > 0) {
    in.facetlist       = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (f = fStart; f < fEnd; ++f) {
      const PetscInt idx     = f - fStart;
      PetscInt      *points = NULL, numPoints, p, numVertices = 0, v;

      in.facetlist[idx].numberofpolygons = 1;
      in.facetlist[idx].polygonlist      = new tetgenio::polygon[in.facetlist[idx].numberofpolygons];
      in.facetlist[idx].numberofholes    = 0;
      in.facetlist[idx].holelist         = NULL;

      ierr = DMPlexGetTransitiveClosure(dm, f, PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
      for (p = 0; p < numPoints*2; p += 2) {
        const PetscInt point = points[p];
        if ((point >= vStart) && (point < vEnd)) points[numVertices++] = point;
      }

      tetgenio::polygon *poly = in.facetlist[idx].polygonlist;
      poly->numberofvertices = numVertices;
      poly->vertexlist       = new int[poly->numberofvertices];
      for (v = 0; v < numVertices; ++v) {
        const PetscInt vIdx = points[v] - vStart;
        poly->vertexlist[v] = vIdx;
      }
	  
#ifdef PETSC_HAVE_EGADS
	  PetscInt	val, bodyOffsetVal;
	  ierr = DMGetLabel(dm, "EGADS Body ID", &label); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(label, f, &val); CHKERRQ(ierr);
	  bodyOffsetVal = val * bodyOffset;
	  
      for (int ii = 0; ii < numLabels; ++ii){
        ierr = DMGetLabelName(dm, ii, &labelName); CHKERRQ(ierr);
		
        ierr = DMGetLabel(dm, labelName, &label); CHKERRQ(ierr);
		ierr = DMLabelGetValue(label, f, &val); CHKERRQ(ierr);
		
		if (val >= 0) {
			if (strcmp(labelName, "EGADS Vertex ID")== 0) { in.facetmarkerlist[idx] = val + vertexOffset + bodyOffsetVal; }		// Store EGADS Vertex date
			if (strcmp(labelName, "EGADS Edge ID") == 0) { in.facetmarkerlist[idx] = val + edgeOffset + bodyOffsetVal; }		// Store EGADS Edge data 
			if (strcmp(labelName, "EGADS Face ID") == 0) { in.facetmarkerlist[idx] = val + faceOffset + bodyOffsetVal; }		// Store EGADS Face data 
		}
	  }
#else
      if (label) {
        PetscInt val;

        ierr = DMLabelGetValue(label, f, &val);CHKERRQ(ierr);
        in.facetmarkerlist[idx] = (int) val;
      }
#endif

      ierr = DMPlexRestoreTransitiveClosure(dm, f, PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);
    }
  }
  
  /* Work with Tetrahedrals of the Current Mesh */
  ierr  = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);

  in.numberofcorners       = 4;
  in.numberoftetrahedra    = cEnd - cStart;
  in.tetrahedronvolumelist = (double*) maxVolumes;
  if (in.numberoftetrahedra > 0) {
    in.tetrahedronlist = new int[in.numberoftetrahedra*in.numberofcorners];
    for (c = cStart; c < cEnd; ++c) {
      const PetscInt idx      = c - cStart;
      PetscInt      *closure = NULL;
      PetscInt       closureSize;

      ierr = DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
      if ((closureSize != 5) && (closureSize != 15)) SETERRQ1(comm, PETSC_ERR_ARG_WRONG, "Mesh has cell which is not a tetrahedron, %D vertices in closure", closureSize);
      for (v = 0; v < 4; ++v) {
        in.tetrahedronlist[idx*in.numberofcorners + v] = closure[(v+closureSize-4)*2] - vStart;
      }
      ierr = DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
    }
  }


  /* Geenerate Refined Mesh */
  if (!rank) {
    char args[32];

#if 1
    /* Take away 'Q' for verbose output */
    ierr = PetscStrcpy(args, "qezQra");CHKERRQ(ierr);
#else
    ierr = PetscStrcpy(args, "qezraVVVV");CHKERRQ(ierr);
#endif
    ::tetrahedralize(args, &in, &out);
  }
  
  in.tetrahedronvolumelist = NULL;

  {
    const PetscInt numCorners  = 4;
    const PetscInt numCells    = out.numberoftetrahedra;
    const PetscInt numVertices = out.numberofpoints;
    const double   *meshCoords = out.pointlist;
    int            *cells      = out.tetrahedronlist;

    PetscBool      interpolate = depthGlobal > 1 ? PETSC_TRUE : PETSC_FALSE;
	
#ifdef PETSC_HAVE_EGADS
	PetscContainer modelObj;
	ego           *bodies;
	ego            model, geom, body;
	int            Nb, oclass, mtype, *senses;
#else
	DMLabel         rlabel     = NULL;
#endif

    ierr = DMPlexInvertCells_Internal(dim, numCells, numCorners, cells);CHKERRQ(ierr);
    ierr = DMPlexCreateFromCellList(comm, dim, numCells, numVertices, numCorners, interpolate, cells, dim, meshCoords, dmRefined);CHKERRQ(ierr);

#ifdef PETSC_HAVE_EGADS
	ierr = DMCreateLabel(*dmRefined, "EGADS Body ID");CHKERRQ(ierr);
	ierr = DMGetLabel(*dmRefined, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
	ierr = DMCreateLabel(*dmRefined, "EGADS Face ID");CHKERRQ(ierr);
	ierr = DMGetLabel(*dmRefined, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
	ierr = DMCreateLabel(*dmRefined, "EGADS Edge ID");CHKERRQ(ierr);
	ierr = DMGetLabel(*dmRefined, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
	ierr = DMCreateLabel(*dmRefined, "EGADS Vertex ID");CHKERRQ(ierr);
	ierr = DMGetLabel(*dmRefined, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
	
	/* Get Attached EGADS Model from Original DMPlex*/
	ierr = PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj);CHKERRQ(ierr);
	ierr = PetscContainerGetPointer(modelObj, (void **) &model);CHKERRQ(ierr);
	ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
	
	/* Transfer EGADS Model to Volumetric Mesh */
	ierr = PetscObjectCompose((PetscObject) *dmRefined, "EGADS Model", (PetscObject) modelObj);CHKERRQ(ierr);
	ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);
#else
    if (label) {ierr = DMCreateLabel(*dmRefined, labelName);
	ierr = DMGetLabel(*dmRefined, labelName, &rlabel);}
#endif

    /* -- Set labels -- */
	/* Vertices Labels */
    for (v = 0; v < numVertices; ++v) {
#ifdef PETSC_HAVE_EGADS
	  PetscInt val, bodyVal, bodyValOffset;
	  
	  val = out.pointmarkerlist[v];
	  if (val >= vertexOffset && val < edgeOffset) {
		val = val - vertexOffset;
		bodyVal = val / bodyOffset;		// Integer Division
		bodyValOffset = bodyVal * bodyOffset;
		ierr = DMLabelSetValue(vertexLabel, v+numCells, out.pointmarkerlist[v] - vertexOffset - bodyValOffset); CHKERRQ(ierr);
		ierr = DMLabelSetValue(bodyLabel, v+numCells, bodyVal); CHKERRQ(ierr);
	  }
	  if (val >= edgeOffset && val < faceOffset) {
		val = val - edgeOffset;
		bodyVal = val / bodyOffset;		// Integer Division
		bodyValOffset = bodyVal * bodyOffset;
		ierr = DMLabelSetValue(edgeLabel, v+numCells, out.pointmarkerlist[v] - edgeOffset - bodyValOffset); CHKERRQ(ierr);
		ierr = DMLabelSetValue(bodyLabel, v+numCells, bodyVal); CHKERRQ(ierr);
	  }
	  if (val >= faceOffset) {
		val = val - faceOffset;
		bodyVal = val / bodyOffset;		// Integer Division
		bodyValOffset = bodyVal * bodyOffset;
	    ierr = DMLabelSetValue(faceLabel, v+numCells, out.pointmarkerlist[v] - faceOffset - bodyValOffset); CHKERRQ(ierr);
		ierr = DMLabelSetValue(bodyLabel, v+numCells, bodyVal); CHKERRQ(ierr);
	  }
#else
      if (out.pointmarkerlist[v]) {
        if (rlabel) {ierr = DMLabelSetValue(rlabel, v+numCells, out.pointmarkerlist[v]);CHKERRQ(ierr);}
      }
#endif
    }
	
    if (interpolate) {
      /* This check is never actually executed for ctetgen (which never returns edgemarkers) and seems to be broken for
       * tetgen */ /* Original Code didn't load edgelist[] or edgemarkerlist[] and is probably why this didn't work previously */
	  /* Set Edge Labels */
	  ierr = DMPlexGetHeightStratum(*dmRefined, 2, &eStart, &eEnd);CHKERRQ(ierr);
      for (e = 0; e < out.numberofedges; e++) {
		ierr = DMLabelSetValue(bodyLabel, e + eStart, 0); CHKERRQ(ierr);
        if (out.edgemarkerlist[e]) {
          const PetscInt  vertices[2] = {out.edgelist[e*2+0]+numCells, out.edgelist[e*2+1]+numCells};
          const PetscInt *edges;
          PetscInt        numEdges;

          ierr = DMPlexGetJoin(*dmRefined, 2, vertices, &numEdges, &edges);CHKERRQ(ierr);
          if (numEdges != 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Two vertices must cover only one edge, not %D", numEdges);
		  
#ifdef PETSC_HAVE_EGADS
		  PetscInt val, bodyVal, bodyValOffset;
	  
		  val = out.edgemarkerlist[e];
		  if (val >= vertexOffset && val < edgeOffset) {
			val = val - vertexOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
		    ierr = DMLabelSetValue(vertexLabel, edges[0], out.edgemarkerlist[e] - vertexOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, edges[0], bodyVal); CHKERRQ(ierr);
		  }
		  if (val >= edgeOffset && val < faceOffset) {
			val = val - edgeOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
	        ierr = DMLabelSetValue(edgeLabel, edges[0], out.edgemarkerlist[e] - edgeOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, edges[0], bodyVal); CHKERRQ(ierr);
		  }
		  if (val >= faceOffset) {
			val = val - faceOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
		    ierr = DMLabelSetValue(faceLabel, edges[0], out.edgemarkerlist[e] - faceOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, edges[0], bodyVal); CHKERRQ(ierr);
		  }
#else
          if (rlabel) {ierr = DMLabelSetValue(rlabel, edges[0], out.edgemarkerlist[e]);CHKERRQ(ierr);}
#endif
          ierr = DMPlexRestoreJoin(*dmRefined, 2, vertices, &numEdges, &edges);CHKERRQ(ierr);
        }
      }
 
	  /* Set Face Labels */
	  ierr = DMPlexGetHeightStratum(*dmRefined, 1, &fStart, &fEnd);CHKERRQ(ierr);
      for (f = 0; f < out.numberoftrifaces; f++) {
		ierr = DMLabelSetValue(bodyLabel, f + fStart, 0); CHKERRQ(ierr);
        if (out.trifacemarkerlist[f]) {
          const PetscInt  vertices[3] = {out.trifacelist[f*3+0]+numCells, out.trifacelist[f*3+1]+numCells, out.trifacelist[f*3+2]+numCells};
          const PetscInt *faces;
          PetscInt        numFaces;

          ierr = DMPlexGetFullJoin(*dmRefined, 3, vertices, &numFaces, &faces);CHKERRQ(ierr);
          if (numFaces != 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Three vertices must cover only one face, not %D", numFaces);
		  
#ifdef PETSC_HAVE_EGADS
		  PetscInt val, bodyVal, bodyValOffset;
	  
		  val = out.trifacemarkerlist[f];
		  if (val >= vertexOffset && val < edgeOffset) {
			val = val - vertexOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
			ierr = DMLabelSetValue(vertexLabel, faces[0], out.trifacemarkerlist[f] - vertexOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, faces[0], bodyVal); CHKERRQ(ierr);
		  }
		  if (val >= edgeOffset && val < faceOffset) {
			val = val - edgeOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
			ierr = DMLabelSetValue(edgeLabel, faces[0], out.trifacemarkerlist[f] - edgeOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, faces[0], bodyVal); CHKERRQ(ierr);
		  }
		  if (val >= faceOffset) {
			val = val - faceOffset;
			bodyVal = val / bodyOffset;		// Integer Division
			bodyValOffset = bodyVal * bodyOffset;
		    ierr = DMLabelSetValue(faceLabel, faces[0], out.trifacemarkerlist[f] - faceOffset - bodyValOffset); CHKERRQ(ierr);
			ierr = DMLabelSetValue(bodyLabel, faces[0], bodyVal); CHKERRQ(ierr);
		  }
#else
          if (rlabel) {ierr = DMLabelSetValue(rlabel, faces[0], out.trifacemarkerlist[f]);CHKERRQ(ierr);}
#endif
          ierr = DMPlexRestoreJoin(*dmRefined, 3, vertices, &numFaces, &faces);CHKERRQ(ierr);
        }
      }
    }
	
	/* Set Cell Labels */
#ifdef PETSC_HAVE_EGADS
	/* NEED TO UPDATE TO DEAL WITH WHAT BODY THE CELL IS LOCATED IN */
	/*   1 - Cycle Through Cells 									*/
	/*   2 - Calculate the centroid (x,y,z) of tetrahedrals			*/
	/*   3 - Perfrom EG_inTopology() function to see which body		*/
	/*   4 - Perform transitive closure								*/
	/*   5 - if attached enitites don't have a bodyID assigned		*/
	/*       assign the current bodyID								*/
	
	ierr = DMPlexGetHeightStratum(*dmRefined, 0, &cStart, &cEnd);CHKERRQ(ierr);
	ierr = DMPlexGetHeightStratum(*dmRefined, 1, &fStart, &fEnd);CHKERRQ(ierr);
	ierr = DMPlexGetHeightStratum(*dmRefined, 2, &eStart, &eEnd);CHKERRQ(ierr);

	for ( c = cStart; c < cEnd; c++) {
		PetscInt 	   val = -1, eVal, fVal;
		PetscReal	   vol;
		PetscReal	   centroid[dim], normal[dim];
		PetscInt      *points = NULL, numPoints;
		
		ierr = DMPlexGetTransitiveClosure(*dmRefined, c, PETSC_TRUE, &numPoints, &points);CHKERRQ(ierr);		// Get all DMPlex points associated with the current cell
		ierr = DMPlexComputeCellGeometryFVM(*dmRefined, c, &vol, centroid, normal); CHKERRQ(ierr);		
		
		/* Deterimine what body the cell's centroid is located in */
		for (int ii = 0; ii < Nb; ++ii){
			body = bodies[ii];
			int cgCheck = EG_inTopology(body, centroid); 
			if (cgCheck == EGADS_SUCCESS) {
				val = ii;
				break;
			}
		}

		if (val >= 0) {
			ierr = DMLabelSetValue(bodyLabel, c, val); CHKERRQ(ierr);
			
			for (int ii = 0; ii < numPoints; ++ii){
				if (points[ii] >= eStart && points[ii] < eEnd) {
					ierr = DMLabelGetValue(bodyLabel, points[ii], &eVal); CHKERRQ(ierr);		// We don't want to overwrite the work done above
					if (eVal < 0) {ierr = DMLabelSetValue(bodyLabel, points[ii], val); CHKERRQ(ierr);}
				} else {
					// Do Nothing - Skip
				}
				
				if (points[ii] >= fStart && points[ii] < fEnd) {
					ierr = DMLabelGetValue(bodyLabel, points[ii], &fVal); CHKERRQ(ierr);		// We don't want to overwrite the work done above
					if (fVal < 0) {ierr = DMLabelSetValue(bodyLabel, points[ii], val); CHKERRQ(ierr);}
				} else {
					// Do Nothing - Skip
				}
			}
		} else {
			// Do Nothing
		}		
		ierr = DMPlexRestoreTransitiveClosure(*dmRefined, c, PETSC_TRUE, &numPoints, &points); CHKERRQ(ierr);
	}
#else
	// Do Nothing
#endif	
    ierr = DMPlexSetRefinementUniform(*dmRefined, PETSC_FALSE);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
