#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/

#ifdef PETSC_HAVE_EGADS
#include <egads.h>
#endif

PetscErrorCode DMPlexSnapToGeomModel(DM dm, PetscInt p, const PetscScalar initCoords[], PetscScalar coords[])
{
  DMLabel         bodyLabel, faceLabel, edgeLabel, vertexLabel;
  PetscContainer  modelObj;
  PetscInt        cdim, d, bodyID, faceID, edgeID, vertexID;
  PetscInt        faceA_ID, faceB_ID;
  PetscInt        edgeA_ID, edgeB_ID;
  PetscInt        vertexA_ID, vertexB_ID;
  PetscInt        vecSize, vStart, vEnd, eConeSize, ni;
  PetscInt        ix[3];
  const PetscInt *eCone = NULL;
  PetscScalar     xyzA[3], xyzB[3];
  Vec             vCoords;
  ego            *bodies;
  ego             model, geom, body, face, edge;
  double          point[3], params[3], result[3];
  double          resultA[3], resultB[3];
  double          eval[18], evalC[18];
  double          errTolr = 1.0e-4;
  double          delta;
  int             Nb, oclass, mtype, *senses, pntCntr;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDim(dm, &cdim);CHKERRQ(ierr);

#ifdef PETSC_HAVE_EGADS
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
  ierr = DMLabelGetValue(bodyLabel, p, &bodyID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel, p, &faceID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, p, &edgeID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(vertexLabel, p, &vertexID);CHKERRQ(ierr);
  ierr = PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj);CHKERRQ(ierr);
  ierr = PetscContainerGetPointer(modelObj, (void **) &model);CHKERRQ(ierr);
  ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
  if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
  body = bodies[bodyID];
  
  //for (d = 0; d < cdim; ++d) point[d] = initCoords[d];  // Don't need
  
  /* Get Vector with all DMPlex Node Coordinates */
  ierr = VecCreate(PETSC_COMM_WORLD, &vCoords); CHKERRQ(ierr);
  ierr = VecSetType(vCoords, VECSTANDARD); CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &vCoords); CHKERRQ(ierr);    
  ierr = VecGetSize(vCoords, &vecSize); CHKERRQ(ierr);
    
  ierr = VecAssemblyBegin(vCoords); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vCoords); CHKERRQ(ierr);
    
  /* Get DMPlex Start ID for Vertices */
  ierr = DMPlexGetHeightStratum(dm, 2, &vStart, &vEnd); CHKERRQ(ierr); CHKERRQ(ierr);
 
  /* Get ID's of Vertices at each end of the Edge */
  ierr = DMPlexGetConeSize(dm, p, &eConeSize); CHKERRQ(ierr);
  ierr = DMPlexGetCone(dm, p, &eCone); CHKERRQ(ierr);
  
  
  /* Get & Report vertexID for the End Points */
  ierr = PetscPrintf(PETSC_COMM_SELF, "  p = %d \n", p); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "   (faceID, edgeID, vertexID) = (%d, %d, %d) \n", faceID, edgeID, vertexID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel, eCone[0], &faceA_ID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, eCone[0], &edgeA_ID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(vertexLabel, eCone[0], &vertexA_ID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "   eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      (faceID, edgeID, vertexID) = (%d, %d, %d) \n", faceA_ID, edgeA_ID, vertexA_ID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel, eCone[1], &faceB_ID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, eCone[1], &edgeB_ID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(vertexLabel, eCone[1], &vertexB_ID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "   eCone[1] = %d \n", eCone[1]); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      (faceID, edgeID, vertexID) = (%d, %d, %d) \n", faceB_ID, edgeB_ID, vertexB_ID); CHKERRQ(ierr);

 
  /* Get Coordinates of Vertices defining the 1st Edge on Face 1 */
  ni = cdim;     // Assume 3D Space ni = 3, need to update to reference cdim
  pntCntr = 0;
  for (int ii = (eCone[0]-vStart)*cdim; ii < (eCone[0]-vStart)*cdim+cdim; ++ii){
    ix[pntCntr] = ii;    //was ix[ii]
    pntCntr = pntCntr + 1;
  }
  ierr = VecGetValues(vCoords, ni, ix, xyzA); CHKERRQ(ierr);    // Vertex A coordinates (x,y,z)

  pntCntr = 0;
  for (int ii = (eCone[1]-vStart)*cdim; ii < (eCone[1]-vStart)*cdim+cdim; ++ii){
    ix[pntCntr] = ii;
    pntCntr = pntCntr + 1;
  }
  ierr = VecGetValues(vCoords, ni, ix, xyzB); CHKERRQ(ierr);  // Vertex B coordinates (x,y,z)
 
    
  /* Snap new Vertex to Geometry */
  if (faceA_ID > 0 || faceB_ID > 0) {
    /* We are on a Geometry Face */
    if (faceA_ID > 0) faceID = faceA_ID;
    if (faceB_ID > 0) faceID = faceB_ID;

    double   paramsA[2], paramsB[2], paramsC[2];
    double   paramsNew[2];
    ierr = EG_objectBodyTopo(body, FACE, faceID, &face);CHKERRQ(ierr);
    
    /* Get Corresponding (u,v) parameters for both endpoints */
    ierr = EG_invEvaluate(face, xyzA, &paramsA, &resultA); CHKERRQ(ierr);    // Get (u,v) parameters for Vertex A
    ierr = EG_invEvaluate(face, xyzB, &paramsB, &resultB); CHKERRQ(ierr);    // Get (u,v) parameters for Vertex B
    
    /* Check if we get the same Point A (x, y, z) Location back -- if not we have to correct */
    ierr = EG_evaluate(face, paramsA, &eval); CHKERRQ(ierr);
    double dx, dy, dz, du, dv;
    double dudx, dudy, dudz, dvdx, dvdy, dvdz;
    
    double range[4];    // [umin, umax, vmin, vmax]
    int    peri;
    ierr = EG_getRange(face, range, &peri); CHKERRQ(ierr);
    
    delta = 0.;   
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzA[ii]) * (eval[ii]-xyzA[ii]))/xyzA[ii];
    }
    delta = sqrt(delta);
    
    if (delta > errTolr){      
      for (int jj=0; jj<6; ++jj){
        ierr = PetscPrintf(PETSC_COMM_SELF, "Point A doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "  jj = %d \n", jj); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    delta = %.6e \n", delta); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsA[0], paramsA[1]); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    range (umin, umax) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    range (vmin, vmax) = (%lf, %lf) \n", range[2], range[3]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    peri = %d \n", peri);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    faceId = %d \n", faceID);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    edgeId = %d \n", edgeA_ID);CHKERRQ(ierr);
        
        dx = eval[0] - xyzA[0];
        dy = eval[1] - xyzA[1];
        dz = eval[2] - xyzA[2];
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dx = %lf \n", dx); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dy = %lf \n", dy); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dz = %lf \n", dz); CHKERRQ(ierr);
        
        if (eval[3] != 0.){
          dudx = 1. / eval[3];
        } else {
          dudx = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dudx = %lf \n", dudx); CHKERRQ(ierr);
        
        if (eval[4] != 0.){
          dudy = 1. / eval[4];
        } else {
          dudy = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dudy = %lf \n", dudy); CHKERRQ(ierr);
        
        if (eval[5] != 0.){
          dudz = 1. / eval[5];
        } else {
          dudz = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dudz = %lf \n", dudz); CHKERRQ(ierr);
        
        if (eval[6] != 0.){
          dvdx = 1. / eval[6];
        } else {
          dvdx = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dvdx = %lf \n", dvdx); CHKERRQ(ierr);
        
        if (eval[7] != 0.){
          dvdy = 1. / eval[7];
        } else {
          dvdy = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dvdy = %lf \n", dvdy); CHKERRQ(ierr);
        
        if (eval[8] != 0.){
          dvdz = 1. / eval[8];
        } else {
          dvdz = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dvdz = %lf \n", dvdz); CHKERRQ(ierr);
        
        du = dudx*dx + dudy*dy + dudz*dz;
        dv = dvdx*dx + dvdy*dy + dvdz*dz;
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    du = %.6e \n", du); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dv = %.6e \n", dv); CHKERRQ(ierr);
        
        if (paramsA[0] - du < range[0]){
          //paramsA[0] = (paramsA[0] + range[0]) / 2.;
          //paramsA[0] = paramsA[0]; // - (du/abs(du))*1.e-4;
          paramsA[0] = range[0];
        } else if (paramsA[0] - du > range[1]) {
          //paramsA[0] = (paramsA[0] + range[1]) / 2.;
          //paramsA[0] = paramsA[0]; // - (du/abs(du))*1.e-4;
          paramsA[0] = range[1];
        } else {
          paramsA[0] = paramsA[0] - du;
        }
        
        if (paramsA[1] - dv < range[2]){
          //paramsA[1] = (paramsA[1] + range[2]) / 2.;
          //paramsA[1] = paramsA[1]; // - (dv/abs(dv))*1.e-4;
          paramsA[1] = range[2];
        } else if (paramsA[1] - dv > range[3]) {
          //paramsA[1] = (paramsA[1] + range[3]) / 2.;
          //paramsA[1] = paramsA[1]; // - (dv/abs(dv))*1.e-4;
          paramsA[1] = range[3];
        } else {
          paramsA[1] = paramsA[1] - dv;
        }
        
        ierr = EG_evaluate(face, paramsA, &eval); CHKERRQ(ierr);
        
        delta = 0.;
        for (int ii=0; ii<3; ++ii){
          delta = delta + ((eval[ii]-xyzA[ii]) * (eval[ii]-xyzA[ii]));
        }
        delta = sqrt(delta);  
      }
    }
    
    /* Check if we get the same Point B (x, y, z) Location back -- if not we have to correct */
    ierr = EG_evaluate(face, paramsB, eval); CHKERRQ(ierr);
    delta = 0.;
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzB[ii]) * (eval[ii]-xyzB[ii]))/xyzB[ii];
    }
    delta = sqrt(delta);
    
    if (delta > errTolr){
      for (int jj=0; jj<6; ++jj){
        ierr = PetscPrintf(PETSC_COMM_SELF, "Point B doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "  jj = %d \n", jj); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    delta = %.6e \n", delta); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[1] = %d \n", eCone[1]); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsB[0], paramsB[1]); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    range (umin, umax) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    range (vmin, vmax) = (%lf, %lf) \n", range[2], range[3]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    peri = %d \n", peri);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    faceId = %d \n", faceID);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    edgeId = %d \n", edgeB_ID);CHKERRQ(ierr);
        
        dx = eval[0] - xyzB[0];
        dy = eval[1] - xyzB[1];
        dz = eval[2] - xyzB[2];
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dx = %lf \n", dx); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dy = %lf \n", dy); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dz = %lf \n", dz); CHKERRQ(ierr);
        
        if (eval[3] != 0.){
          dudx = 1. / eval[3];
        } else {
          dudx = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dudx = %lf \n", dudx); CHKERRQ(ierr);
        
        if (eval[4] != 0.){
          dudy = 1. / eval[4];
        } else {
          dudy = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dudy = %lf \n", dudy); CHKERRQ(ierr);
        
        if (eval[5] != 0.){
          dudz = 1. / eval[5];
        } else {
          dudz = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dudz = %lf \n", dudz); CHKERRQ(ierr);
        
        if (eval[6] != 0.){
          dvdx = 1. / eval[6];
        } else {
          dvdx = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dvdx = %lf \n", dvdx); CHKERRQ(ierr);
        
        if (eval[7] != 0.){
          dvdy = 1. / eval[7];
        } else {
          dvdy = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dvdy = %lf \n", dvdy); CHKERRQ(ierr);
        
        if (eval[8] != 0.){
          dvdz = 1. / eval[8];
        } else {
          dvdz = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dvdz = %lf \n", dvdz); CHKERRQ(ierr);
        
        du = dudx*dx + dudy*dy + dudz*dz;
        dv = dvdx*dx + dvdy*dy + dvdz*dz;
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    du = %.6e \n", du); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dv = %.6e \n", dv); CHKERRQ(ierr);
        
        if (paramsB[0] - du < range[0]){
          //paramsB[0] = (paramsB[0] + range[0]) / 2.;
          //paramsB[0] = paramsB[0]; // - (du/abs(du))*1.e-4;
          paramsB[0] = range[0];
        } else if (paramsB[0] - du > range[1]) {
          //paramsB[0] = (paramsB[0] + range[1]) / 2.;
          //paramsB[0] = paramsB[0]; // - (du/abs(du))*1.e-4;
          paramsB[0] = range[1];
        } else {
          paramsB[0] = paramsB[0] - du;
        }
        
        if (paramsB[1] - dv < range[2]){
          //paramsB[1] = (paramsB[1] + range[2]) / 2.;
          //paramsB[1] = paramsB[1]; // - (dv/abs(dv))*1.e-4;
          paramsB[1] = range[2];
        } else if (paramsB[1] - dv > range[3]) {
          //paramsB[1] = (paramsB[1] + range[3]) / 2.;
          //paramsB[1] = paramsB[1]; // - (dv/abs(dv))*1.e-4;
          paramsB[1] = range[3];
        } else {
          paramsB[1] = paramsB[1] - dv;
        }
        
        ierr = EG_evaluate(face, paramsB, &eval); CHKERRQ(ierr);
        
        delta = 0.;
        for (int ii=0; ii<3; ++ii){
          delta = delta + ((eval[ii]-xyzB[ii]) * (eval[ii]-xyzB[ii]));
        }
        delta = sqrt(delta);
      }  
    }   
    
    /* Calculate (u, v) parameters for refined point */
    paramsNew[0] = (paramsA[0] + paramsB[0]) / 2.;
    paramsNew[1] = (paramsA[1] + paramsB[1]) / 2.;
    
    //int paramTest = 0;
    double xMax, xMin, yMax, yMin, zMax, zMin;
    //if ((paramsNew[0] >= range[0] && paramsNew[0] <= range[1]) && (paramsNew[1] >= range[2] && paramsNew[1] <= range[3])){
      ierr = EG_evaluate(face, paramsNew, &eval); CHKERRQ(ierr);     
    //} else {
    //  ierr = EG_invEvaluate(face, initCoords, &paramsC, &eval); CHKERRQ(ierr);
    //}
    
    //if ((eCone[0] == vStart + 53 || eCone[1] == vStart + 53)) { // && (eCone[0] == vStart - 1 + 519 || eCone[1] == vStart - 1 + 519)){
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "  --------------------- \n"); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "  -- FACE Refinement -- \n"); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "  --------------------- \n"); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    point p = %d \n", p); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    Face ID = %d \n", faceID); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    Edge ID = %d \n", edgeID); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    vStart = %d \n", vStart); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    peri = %d \n", peri); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsA[0], paramsA[1]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[1] = %d \n", eCone[1]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsB[0], paramsB[1]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsNew[0], paramsNew[1]); CHKERRQ(ierr);
    //  //ierr = PetscPrintf(PETSC_COMM_SELF, "    initCoords (x, y, z) = (%lf, %lf, %lf) \n", initCoords[0], initCoords[1], initCoords[2]);CHKERRQ(ierr);
    //  //ierr = PetscPrintf(PETSC_COMM_SELF, "    evalC (x, y, z) = (%lf, %lf, %lf) \n", evalC[0], evalC[1], evalC[2]);CHKERRQ(ierr);
    //  //ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsC[0], paramsC[1]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    range(ustart, uend) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "         (vstart, vend) = (%lf, %lf) \n", range[2], range[3]);CHKERRQ(ierr);
    //}
    
    /* Set Correct BodyID, FaceID & edgeID for refinement point */
    ierr = DMLabelSetValue(bodyLabel, p, bodyID);CHKERRQ(ierr);
    ierr = DMLabelSetValue(faceLabel, p, faceID);CHKERRQ(ierr);
    //ierr = DMLabelSetValue(edgeLabel, p, edgeID);CHKERRQ(ierr);      


  } else {
  //if ((vertexA_ID > 0 || vertexB_ID > 0) && edgeA_ID > 0 && edgeB_ID > 0){
    /* We are on an Geometry Edge */
    //if (vertexA_ID > 0) {edgeID = edgeB_ID;}
    //if (vertexB_ID > 0) {edgeID = edgeA_ID;}
    if (edgeA_ID > 0) {edgeID = edgeA_ID;}
    if (edgeB_ID > 0) {edgeID = edgeB_ID;}
    
    //ierr = PetscPrintf(PETSC_COMM_SELF, "      p = %d \n", p); CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_SELF, "      edgeID = %d \n", edgeID); CHKERRQ(ierr);
    
    /* Get (t) parameter for each endpoint */
    double   paramsA[1], paramsB[1], paramsC[1];
    double   paramsNew[1];
    ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge);CHKERRQ(ierr);
    ierr = EG_invEvaluate(edge, xyzA, paramsA, resultA); CHKERRQ(ierr);   // t- parameter for Vertex A
    ierr = EG_invEvaluate(edge, xyzB, paramsB, resultB); CHKERRQ(ierr);   // t- parameter for Vertex B
    
    /* Check that (t) parameter gives the same (x,y,z) coordinates stored as part of the DMPlex */
    /* If not, correct (t) parameter. This sometimes happens in the EGADS library.              */
    double range[4];    // [umin, umax, vmin, vmax]
    int    peri;
    ierr = EG_getRange(edge, range, &peri); CHKERRQ(ierr);
    
    /* Verify right (t) parameter for Point A*/
    ierr = EG_evaluate(edge, paramsA, eval);
    
    delta = 0.;
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzA[ii]) * (eval[ii]-xyzA[ii]))/xyzA[ii];
    }
    delta = sqrt(delta);
    
    if (delta > errTolr){
      /* Output Some Diagnostic Data */
      double dx, dy, dz, dt;
      double dtdx, dtdy, dtdz;
             
      for (int jj=0; jj<6; ++jj){
        ierr = PetscPrintf(PETSC_COMM_SELF, " EDGE Point A doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "  jj = %d \n", jj); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "  delta = %.6e \n", delta); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "               (t) = (%lf) \n", paramsA[0]); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    range (tmin, tmax) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    peri = %d \n", peri);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    faceId = %d \n", faceID);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    edgeId = %d \n", edgeID);CHKERRQ(ierr);
      
        dx = eval[0] - xyzA[0];
        dy = eval[1] - xyzA[1];
        dz = eval[2] - xyzA[2];
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dx = %lf \n", dx); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dy = %lf \n", dy); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dz = %lf \n", dz); CHKERRQ(ierr);
        
        if (eval[3] != 0.){
          dtdx = 1. / eval[3];
        } else {
          dtdx = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dtdx = %lf \n", dtdx); CHKERRQ(ierr);
        
        if (eval[4] != 0.){
          dtdy = 1. / eval[4];
        } else {
          dtdy = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dtdy = %lf \n", dtdy); CHKERRQ(ierr);
        
        if (eval[5] != 0.){
          dtdz = 1. / eval[5];
        } else {
          dtdz = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dtdz = %lf \n", dtdz); CHKERRQ(ierr);
        
        dt = dtdx*dx + dtdy*dy + dtdz*dz;
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dt = %.6e \n", dt); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[3] = %.6e \n", eval[3]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[4] = %.6e \n", eval[4]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[5] = %.6e \n", eval[5]); CHKERRQ(ierr);
        
        if (paramsA[0] - dt < range[0]){
          //paramsA[0] = (paramsA[0] + range[0]) / 2.;
          //paramsA[0] = paramsA[0]; // - (dt/abs(dt))*1.e-4;
          paramsA[0] = range[0];
        } else if (paramsA[0] - dt > range[1]) {
          //paramsA[0] = (paramsA[0] + range[1]) / 2.;
          //paramsA[0] = paramsA[0]; // - (dt/abs(dt))*1.e-4;
          paramsA[0] = range[0];
        } else {
          paramsA[0] = paramsA[0] - dt;
        }
        
        ierr = EG_evaluate(edge, paramsA, &eval); CHKERRQ(ierr);
        
        delta = 0.;
        for (int ii=0; ii<3; ++ii){
          delta = delta + ((eval[ii]-xyzA[ii]) * (eval[ii]-xyzA[ii]));
        }
        delta = sqrt(delta);
      }
    } 
    
    /* Verify right (t) parameter for Point B */    
    ierr = EG_evaluate(edge, paramsB, eval);
    
    delta = 0.;
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzB[ii]) * (eval[ii]-xyzB[ii]))/xyzB[ii];
    }
    delta = sqrt(delta);
    
    if (delta > errTolr){
      double dx, dy, dz, dt;
      double dtdx, dtdy, dtdz;
       
      for (int jj=0; jj<6; ++jj){
       ierr = PetscPrintf(PETSC_COMM_SELF, " EDGE Point B doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "  jj = %d \n", jj); CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "    delta = %.6e \n", delta); CHKERRQ(ierr);
       //ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]);CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "               (t) = (%lf) \n", paramsB[0]); CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "    range (tmin, tmax) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
       //ierr = PetscPrintf(PETSC_COMM_SELF, "    range (vmin, vmax) = (%lf, %lf) \n", range[2], range[3]);CHKERRQ(ierr);
       //ierr = PetscPrintf(PETSC_COMM_SELF, "    peri = %d \n", peri);CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "    faceId = %d \n", faceID);CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "    edgeId = %d \n", edgeID);CHKERRQ(ierr);
      
       dx = eval[0] - xyzB[0];
       dy = eval[1] - xyzB[1];
       dz = eval[2] - xyzB[2];
        
       ierr = PetscPrintf(PETSC_COMM_SELF, "    dx = %lf \n", dx); CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "    dy = %lf \n", dy); CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_SELF, "    dz = %lf \n", dz); CHKERRQ(ierr);
        
        if (eval[3] != 0.){
          dtdx = 1. / eval[3];
        } else {
          dtdx = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dtdx = %lf \n", dtdx); CHKERRQ(ierr);
        
        if (eval[4] != 0.){
          dtdy = 1. / eval[4];
        } else {
          dtdy = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dtdy = %lf \n", dtdy); CHKERRQ(ierr);
        
        if (eval[5] != 0.){
          dtdz = 1. / eval[5];
        } else {
          dtdz = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dtdz = %lf \n", dtdz); CHKERRQ(ierr);
        
        dt = dtdx*dx + dtdy*dy + dtdz*dz;
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dt = %.6e \n", dt); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[3] = %.6e \n", eval[3]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[4] = %.6e \n", eval[4]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[5] = %.6e \n", eval[5]); CHKERRQ(ierr);
        
        if (paramsB[0] - dt < range[0]){
          //paramsB[0] = (paramsB[0] + range[0]) / 2.;
          //paramsB[0] = paramsB[0]; // - (dt/abs(dt))*1.e-4;
          paramsB[0] = range[0];
        } else if (paramsB[0] - dt > range[1]) {
          //paramsB[0] = (paramsB[0] + range[1]) / 2.;
          //paramsB[0] = paramsB[0]; // - (dt/abs(dt))*1.e-4;
          paramsB[0] = range[0];
        } else {
          paramsB[0] = paramsB[0] - dt;
        }
        
        ierr = EG_evaluate(edge, paramsB, &eval); CHKERRQ(ierr);
        
        delta = 0.;
        for (int ii=0; ii<3; ++ii){
          delta = delta + ((eval[ii]-xyzB[ii]) * (eval[ii]-xyzB[ii]));
        }
        delta = sqrt(delta);
      }
    } // End of Check
    
    /* Calculate (t) parameter for refinement point (average of 2 endpoints) */
    paramsNew[0] = (paramsA[0] + paramsB[0]) / 2.;          // t- parameter for New Vertex at midpoint    
    ierr = EG_evaluate(edge, paramsNew, eval); CHKERRQ(ierr);
    
    /* Set Correct BodyID, FaceID & edgeID for refinement point */
    ierr = DMLabelSetValue(bodyLabel, p, bodyID);CHKERRQ(ierr);
    //ierr = DMLabelSetValue(faceLabel, p, faceID);CHKERRQ(ierr);
    ierr = DMLabelSetValue(edgeLabel, p, edgeID);CHKERRQ(ierr);
  }

  /* Output Results */
  ierr = PetscPrintf(PETSC_COMM_SELF, " After Refinement \n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "  p = %d \n", p); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "   (faceID, edgeID, vertexID) = (%d, %d, %d) \n", faceID, edgeID, vertexID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, " ------------------------ \n"); CHKERRQ(ierr);

  /* Set Coordinates for Refinement Point */
  for (d = 0; d < cdim; ++d) coords[d] = eval[d];
#else
  for (d = 0; d < cdim; ++d) coords[d] = initCoords[d];
#endif
  PetscFunctionReturn(0);
}
