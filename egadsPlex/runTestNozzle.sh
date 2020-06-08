#!/bin/bash

# Compile Code
echo "---------------------------"
echo "     Recompiling Code"
echo "---------------------------"
make clean
make testNozzle

# Run code with all options
./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle.h5 -dm_view2 hdf5:mesh_vol_SM.h5
./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r1.h5 -dm_view2 hdf5:mesh_vol_r1_SM.h5 -dm_view3 hdf5:mesh_vol_r1_pSM.h5 -dm_refine 1
./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r2.h5 -dm_view2 hdf5:mesh_vol_r2_SM.h5 -dm_view3 hdf5:mesh_vol_r2_pSM.h5 -dm_refine 2
./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r3.h5 -dm_view2 hdf5:mesh_vol_r3_SM.h5 -dm_view3 hdf5:mesh_vol_r3_pSM.h5 -dm_refine 3
./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r4.h5 -dm_view2 hdf5:mesh_vol_r4_SM.h5 -dm_view3 hdf5:mesh_vol_r4_pSM.h5 -dm_refine 4
./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r5.h5 -dm_view2 hdf5:mesh_vol_r5_SM.h5 -dm_view3 hdf5:mesh_vol_r5_pSM.h5 -dm_refine 5
#./testNozzle -filename ../examples/Nozzle_fromStep.egadslite -dm_view hdf5:mesh_nozzle_r10.h5 -dm_view2 hdf5:mesh_vol_r10.h5 -dm_view3 hdf5:mesh_vol_r1_pSM.h5 -dm_refine 10
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r10.h5 -dm_refine 10
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r25.h5 -dm_refine 25

