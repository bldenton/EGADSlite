#!/bin/bash

# Compile Code
echo "---------------------------"
echo "     Recompiling Code"
echo "---------------------------"
make clean
make testNozzle

# Run code with all options
# Sphere and Nozzle Multisolid Development Files
#./testNozzle -filename ../examples/sphere_and_nozzle_tryB.egadslite -dm_view hdf5:mesh_sphere_and_nozzle.h5 -dm_view2 hdf5:mesh_vol_sphere_and_nozzle.h5 -dm_view3 hdf5:mesh_vol_pSaN.h5 -dm_view4 hdf5:mesh_vol_r1_SaN.h5 -dm_view5 hdf5:mesh_vol_r2_SaN.h5 -dm_view6 hdf5:mesh_vol_r3_SaN.h5 -dm_view7 hdf5:mesh_vol_r4_SaN.h5 -dm_view8 hdf5:mesh_vol_r5_SaN.h5 -dm_refine 1

#  Mulitsolid Mesh Creation
./testNozzle -filename ../examples/multisolids_v2.egadslite -dm_view hdf5:mesh_multisolids.h5 -dm_view2 hdf5:mesh_vol_multisolids.h5 -dm_view3 hdf5:mesh_vol_pMultisolids.h5 -dm_view4 hdf5:mesh_vol_r1_Multisolids.h5 -dm_view5 hdf5:mesh_vol_r2_Multisolids.h5 -dm_view6 hdf5:mesh_vol_r3_Multisolids.h5 -dm_view7 hdf5:mesh_vol_r4_Multisolids.h5 -dm_view8 hdf5:mesh_vol_r5_Multisolids.h5 -dm_refine 1

# Nozzle Development Files
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle.h5 -dm_view2 hdf5:mesh_vol_SM.h5
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r1.h5 -dm_view2 hdf5:mesh_vol_r1_RR.h5 -dm_view3 hdf5:mesh_vol_pRR.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5 -dm_refine 1
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r2.h5 -dm_view2 hdf5:mesh_vol_r2_SM.h5 -dm_view3 hdf5:mesh_vol_r2_pSM.h5 -dm_refine 2
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r3.h5 -dm_view2 hdf5:mesh_vol_r3_SM.h5 -dm_view3 hdf5:mesh_vol_r3_pSM.h5 -dm_refine 3
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r4.h5 -dm_view2 hdf5:mesh_vol_r4_SM.h5 -dm_view3 hdf5:mesh_vol_r4_pSM.h5 -dm_refine 4
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r5.h5 -dm_view2 hdf5:mesh_vol_r5_SM.h5 -dm_view3 hdf5:mesh_vol_r5_pSM.h5 -dm_refine 5
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r10.h5 -dm_view2 hdf5:mesh_vol_r10.h5 -dm_view3 hdf5:mesh_vol_r1_pSM.h5 -dm_refine 10
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r10.h5 -dm_refine 10
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r25.h5 -dm_refine 25

