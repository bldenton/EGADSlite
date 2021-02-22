#!/bin/bash

# Compile Code
echo "---------------------------"
echo "     Recompiling Code"
echo "---------------------------"
make clean
make testNozzleNew

# Run code with all options
# Sphere and Nozzle Multisolid Development Files
#./testNozzleNew -filename ../examples/sphere_and_nozzle_tryB.egadslite -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_sphere_and_nozzle.h5 -dm_view2 hdf5:mesh_vol_sphere_and_nozzle.h5 -dm_view3 hdf5:mesh_vol_pSaN.h5 -dm_view4 hdf5:mesh_vol_r1_SaN.h5 -dm_view5 hdf5:mesh_vol_r2_SaN.h5 -dm_view6 hdf5:mesh_vol_r3_SaN.h5 -dm_view7 hdf5:mesh_vol_r4_SaN.h5 -dm_view8 hdf5:mesh_vol_r5_SaN.h5 -dm_refine 1

#  Mulitsolid Mesh Creation
#./testNozzleNew -filename ../examples/multisolids_v2.egadslite -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_multisolids.h5 -dm_view2 hdf5:mesh_vol_multisolids.h5 -dm_view3 hdf5:mesh_vol_pMultisolids.h5 -dm_view4 hdf5:mesh_vol_r1_Multisolids.h5 -dm_view5 hdf5:mesh_vol_r2_Multisolids.h5 -dm_view6 hdf5:mesh_vol_r3_Multisolids.h5 -dm_view7 hdf5:mesh_vol_r4_Multisolids.h5 -dm_view8 hdf5:mesh_vol_r5_Multisolids.h5 -dm_refine 1

#  Sphere Mesh Creation
#./testNozzleNew -filename ../examples/sphere_example.egadslite -dm_refine 0 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_Sphere.h5 -dm_view2 hdf5:mesh_vol_Sphere.h5 -dm_view3 hdf5:mesh_vol_pSphere.h5 -dm_view4 hdf5:mesh_vol_r1_Sphere.h5 -dm_view5 hdf5:mesh_vol_r2_Sphere.h5 -dm_view6 hdf5:mesh_vol_r3_Sphere.h5 -dm_view7 hdf5:mesh_vol_r4_Sphere.h5 -dm_view8 hdf5:mesh_vol_r5_Sphere.h5

#  Cylinder Mesh Creation
#./testNozzleNew -filename ../examples/cylinder_example.egadslite -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_Cylinder.h5 -dm_view2 hdf5:mesh_vol_Cylinder.h5 -dm_view3 hdf5:mesh_vol_pCylinder.h5 -dm_view4 hdf5:mesh_vol_r1_Cylinder.h5 -dm_view5 hdf5:mesh_vol_r2_Cylinder.h5 -dm_view6 hdf5:mesh_vol_r3_Cylinder.h5 -dm_view7 hdf5:mesh_vol_r4_Cylinder.h5 -dm_view8 hdf5:mesh_vol_r5_Cylinder.h5

#  Cone Mesh Creation using Universal Mesh Approach
#./testNozzleNew -filename ../examples/cone_example.egadslite -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_uni_cone.h5 -dm_view2 hdf5:mesh_uni_vol_cone.h5 -dm_view3 hdf5:mesh_vol_pCone.h5 -dm_view4 hdf5:mesh_vol_r1_Cone.h5 -dm_view5 hdf5:mesh_vol_r2_Cone.h5 -dm_view6 hdf5:mesh_vol_r3_Cone.h5 -dm_view7 hdf5:mesh_vol_r4_Cone.h5 -dm_view8 hdf5:mesh_vol_r5_Cone.h5

#  Cone Mesh Creation using Tessellation
#./testNozzleNew -filename ../examples/cone_example.egadslite -dm_plex_egads_with_tess 1 -dm_refine 0 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_tess_cone_r0.h5 -dm_view2 hdf5:mesh_tess_vol_cone_r0.h5 -dm_view3 hdf5:mesh_vol_pCone.h5 -dm_view4 hdf5:mesh_vol_r1_Cone.h5 -dm_view5 hdf5:mesh_vol_r2_Cone.h5 -dm_view6 hdf5:mesh_vol_r3_Cone.h5 -dm_view7 hdf5:mesh_vol_r4_Cone.h5 -dm_view8 hdf5:mesh_vol_r5_Cone.h5
#./testNozzleNew -filename ../examples/cone_example.egadslite -dm_plex_egads_with_tess 1 -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_tess_cone_r1.h5 -dm_view2 hdf5:mesh_tess_vol_cone_r1.h5 -dm_view3 hdf5:mesh_vol_pCone.h5 -dm_view4 hdf5:mesh_vol_r1_Cone.h5 -dm_view5 hdf5:mesh_vol_r2_Cone.h5 -dm_view6 hdf5:mesh_vol_r3_Cone.h5 -dm_view7 hdf5:mesh_vol_r4_Cone.h5 -dm_view8 hdf5:mesh_vol_r5_Cone.h5

# Complex Development Files
#./testNozzleNew -filename ../examples/complex.egadslite -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_complex.h5 -dm_view2 hdf5:mesh_vol_complex.h5 -dm_view3 hdf5:mesh_vol_complex.h5 -dm_view4 hdf5:mesh_vol_complex_r1.h5 -dm_view5 hdf5:mesh_vol_complex_r2.h5 -dm_view6 hdf5:mesh_vol_complex_r3.h5 -dm_view7 hdf5:mesh_vol_complex_r4.h5 -dm_view8 hdf5:mesh_vol_complex_r5.h5 -dm_refine 1

# Nozzle Mesh Creation using Universal Mesh Approach
#./testNozzleNew -filename ../examples/Nozzle_example.egadslite -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_uni_nozzle.h5 -dm_view2 hdf5:mesh_uni_vol_nozzle.h5 -dm_view3 hdf5:mesh_vol_pRR.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r10.h5 -dm_refine 10
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r25.h5 -dm_refine 25

# Nozzle Mesh Creation using Tessellation
./testNozzleNew -filename ../examples/Nozzle_example.egadslite -dm_plex_egads_with_tess 1 -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_tess_nozzle.h5 -dm_view2 hdf5:mesh_tess_vol_nozzle.h5 -dm_view3 hdf5:mesh_tess_vol_nozzle_postRefine.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5


#${PETSC_DIR}/lib/petsc/bin/petsc_gen_xdmf.py mesh_nozzle.h5

