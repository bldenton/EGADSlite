#!/bin/bash

# Compile Code
echo "---------------------------"
echo "     Recompiling Code"
echo "---------------------------"
make clean
make testPlexEx45

# Run code with all options
# Sphere and Nozzle Multisolid Development Files
#./testNozzle -filename ../examples/sphere_and_nozzle_tryB.egadslite -dm_view hdf5:mesh_sphere_and_nozzle.h5 -dm_view2 hdf5:mesh_vol_sphere_and_nozzle.h5 -dm_view3 hdf5:mesh_vol_pSaN.h5 -dm_view4 hdf5:mesh_vol_r1_SaN.h5 -dm_view5 hdf5:mesh_vol_r2_SaN.h5 -dm_view6 hdf5:mesh_vol_r3_SaN.h5 -dm_view7 hdf5:mesh_vol_r4_SaN.h5 -dm_view8 hdf5:mesh_vol_r5_SaN.h5 -dm_refine 1

#  Mulitsolid Mesh Creation
#./testNozzle -filename ../examples/multisolids_v2.egadslite -dm_view hdf5:mesh_multisolids.h5 -dm_view2 hdf5:mesh_vol_multisolids.h5 -dm_view3 hdf5:mesh_vol_pMultisolids.h5 -dm_view4 hdf5:mesh_vol_r1_Multisolids.h5 -dm_view5 hdf5:mesh_vol_r2_Multisolids.h5 -dm_view6 hdf5:mesh_vol_r3_Multisolids.h5 -dm_view7 hdf5:mesh_vol_r4_Multisolids.h5 -dm_view8 hdf5:mesh_vol_r5_Multisolids.h5 -dm_refine 1

#  Sphere Mesh Creation
#./testNozzle -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_Sphere.h5 -dm_view2 hdf5:mesh_vol_Sphere.h5 -dm_view3 hdf5:mesh_vol_pSphere.h5 -dm_view4 hdf5:mesh_vol_r1_Sphere.h5 -dm_view5 hdf5:mesh_vol_r2_Sphere.h5 -dm_view6 hdf5:mesh_vol_r3_Sphere.h5 -dm_view7 hdf5:mesh_vol_r4_Sphere.h5 -dm_view8 hdf5:mesh_vol_r5_Sphere.h5 -dm_refine 1

#  Cylinder Mesh Creation
#./testNozzle -filename ../examples/cylinder_example.egadslite -dm_view hdf5:mesh_Cylinder.h5 -dm_view2 hdf5:mesh_vol_Cylinder.h5 -dm_view3 hdf5:mesh_vol_pCylinder.h5 -dm_view4 hdf5:mesh_vol_r1_Cylinder.h5 -dm_view5 hdf5:mesh_vol_r2_Cylinder.h5 -dm_view6 hdf5:mesh_vol_r3_Cylinder.h5 -dm_view7 hdf5:mesh_vol_r4_Cylinder.h5 -dm_view8 hdf5:mesh_vol_r5_Cylinder.h5 -dm_refine 1

#  Cone Mesh Creation
#./testNozzle -filename ../examples/cone_example.egadslite -dm_view hdf5:mesh_Cone.h5 -dm_view2 hdf5:mesh_vol_Cone.h5 -dm_view3 hdf5:mesh_vol_pCone.h5 -dm_view4 hdf5:mesh_vol_r1_Cone.h5 -dm_view5 hdf5:mesh_vol_r2_Cone.h5 -dm_view6 hdf5:mesh_vol_r3_Cone.h5 -dm_view7 hdf5:mesh_vol_r4_Cone.h5 -dm_view8 hdf5:mesh_vol_r5_Cone.h5 -dm_refine 1

# Complex Development Files
#./testNozzle -filename ../examples/complex.egadslite -dm_view hdf5:mesh_complex.h5 -dm_view2 hdf5:mesh_vol_complex.h5 -dm_view3 hdf5:mesh_vol_complex.h5 -dm_view4 hdf5:mesh_vol_complex_r1.h5 -dm_view5 hdf5:mesh_vol_complex_r2.h5 -dm_view6 hdf5:mesh_vol_complex_r3.h5 -dm_view7 hdf5:mesh_vol_complex_r4.h5 -dm_view8 hdf5:mesh_vol_complex_r5.h5 -dm_refine 1

# Nozzle Development Files
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle.h5 -dm_view2 hdf5:mesh_vol_SM.h5
./testPlexEx45 -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r1.h5 -dm_view2 hdf5:mesh_vol_r1_RR.h5 -dm_view3 hdf5:mesh_vol_pRR.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5 -dm_view9 hdf5:mesh_vol_TempSol.h5 -dim 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor -sol_vec_view hdf5:mesh_vol_TempSol.h5::append 
./testPlexEx45 -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r1.h5 -dm_view2 hdf5:mesh_vol_r1_RR.h5 -dm_view3 hdf5:mesh_vol_pRR.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5 -dm_view9 hdf5:mesh_vol_r1_TempSol.h5 -dm_refine 1 -dim 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor -sol_vec_view hdf5:mesh_vol_r1_TempSol.h5::append 
./testPlexEx45 -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r1.h5 -dm_view2 hdf5:mesh_vol_r1_RR.h5 -dm_view3 hdf5:mesh_vol_pRR.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5 -dm_view9 hdf5:mesh_vol_r2_TempSol.h5 -dm_refine 2 -dim 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor -sol_vec_view hdf5:mesh_vol_r2_TempSol.h5::append 
./testPlexEx45 -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r1.h5 -dm_view2 hdf5:mesh_vol_r1_RR.h5 -dm_view3 hdf5:mesh_vol_pRR.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5 -dm_view9 hdf5:mesh_vol_r3_TempSol.h5 -dm_refine 3 -dim 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor -sol_vec_view hdf5:mesh_vol_r3_TempSol.h5::append 
./testPlexEx45 -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r1.h5 -dm_view2 hdf5:mesh_vol_r1_RR.h5 -dm_view3 hdf5:mesh_vol_pRR.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5 -dm_view9 hdf5:mesh_vol_r4_TempSol.h5 -dm_refine 4 -dim 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor -sol_vec_view hdf5:mesh_vol_r4_TempSol.h5::append 
./testPlexEx45 -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r1.h5 -dm_view2 hdf5:mesh_vol_r1_RR.h5 -dm_view3 hdf5:mesh_vol_pRR.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5 -dm_view9 hdf5:mesh_vol_r5_TempSol.h5 -dm_refine 5 -dim 3 -temp_petscspace_degree 1 -ts_type beuler -ts_max_steps 10 -ts_dt 0.1 -pc_type lu -ksp_monitor_short -ksp_converged_reason -snes_monitor_short -snes_converged_reason -ts_monitor -sol_vec_view hdf5:mesh_vol_r5_TempSol.h5::append 
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r2.h5 -dm_view2 hdf5:mesh_vol_r2_SM.h5 -dm_view3 hdf5:mesh_vol_r2_pSM.h5 -dm_refine 2
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r3.h5 -dm_view2 hdf5:mesh_vol_r3_SM.h5 -dm_view3 hdf5:mesh_vol_r3_pSM.h5 -dm_refine 3
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r4.h5 -dm_view2 hdf5:mesh_vol_r4_SM.h5 -dm_view3 hdf5:mesh_vol_r4_pSM.h5 -dm_refine 4
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r5.h5 -dm_view2 hdf5:mesh_vol_r5_SM.h5 -dm_view3 hdf5:mesh_vol_r5_pSM.h5 -dm_refine 5
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r10.h5 -dm_view2 hdf5:mesh_vol_r10.h5 -dm_view3 hdf5:mesh_vol_r1_pSM.h5 -dm_refine 10
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r10.h5 -dm_refine 10
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r25.h5 -dm_refine 25

#${PETSC_DIR}/lib/petsc/bin/petsc_gen_xdmf.py mesh_nozzle.h5

