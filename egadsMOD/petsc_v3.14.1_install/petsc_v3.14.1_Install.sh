#!/bin/bash

# Unload petsc module if loaded
echo "Unloading Petsc if Loaded"
module unload petsc

# Go to build directory
cd $UBCESLAB_SWENV_PREFIX/builddir

# Remove current unpacked folder
echo "Removing Previous Install"
rm -rf petsc-3.14.1

# Unpack tarball
echo "Unpacking Original Petsc Tarball"
tar -xzf $UBCESLAB_SWENV_PREFIX/sourcesdir/petsc/petsc-lite-3.14.1.tar.gz

# Petsc v3.14.1 doesn't need the patch. Only Petsc v 3.12.0 does
## Copy patch for hdf5
#cp /kif1/data/shared/software/patch/hdf5.patch /kif1/data/shared/software/builddir/petsc-3.12.0/hdf5.patch
#
## Apply the patch
#cd /kif1/data/shared/software/builddir/petsc-3.12.0
#patch -p1 < hdf5.patch

echo "Inserting Modified Code before Install"

# Copy Modified plexrefine  (Already in v3.14.1 release)
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/plexrefine.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/impls/plex/plexrefine.c

# Copy Modified petscvariables
#cp /kif1/data/shared/software/egadsMOD/petscvariables /kif1/data/shared/software/builddir/petsc-3.12.0/lib/petsc/conf/petscvariables

# Copy egads.py	(Alredy in v3.14.1 release)
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/egads.py $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/config/BuildSystem/config/packages

# Copy makefile with EGADS library	(Already in v3.14.1 release)
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/makefile $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/impls/plex/makefile

# Copy plexegads.c to correct location	(Modified version in v3.14.1 release -- MAY WANT TO CHANGE BACK TO RUN MY CODES -- NEED TO TRY THE RELEASED VERSION FIRST)
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/plexegads.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/impls/plex/plexegads.c

# Copy code related to branch Knepley/feature-tetgen-labels
ORIG_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD
cp $ORIG_DIR/petsc_v3.14.1_install/petscdm.h $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/include/petscdm.h
cp $ORIG_DIR/petsc_v3.14.1_install/petscdmtypes.h $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/include/petscdmtypes.h
cp $ORIG_DIR/petsc_v3.14.1_install/dmimpl.h $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/include/petsc/private/dmimpl.h
cp $ORIG_DIR/petsc_v3.14.1_install/dm.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/interface/dm.c
cp $ORIG_DIR/petsc_v3.14.1_install/tetgenerate.cxx $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/impls/plex/generators/tetgen/tetgenerate.cxx
cp $ORIG_DIR/petsc_v3.14.1_install/ctetgenerate.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/impls/plex/generators/ctetgen/ctetgenerate.c
cp $ORIG_DIR/petsc_v3.14.1_install/ex11.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/impls/plex/tests/ex11.c
cp $ORIG_DIR/petsc_v3.14.1_install/ex11_univ.out $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/impls/plex/tests/output/ex11_univ.out
cp $ORIG_DIR/petsc_v3.14.1_install/ex11_univ_egads_ball.out $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/impls/plex/tests/output/ex11_univ_egads_ball.out
cp $ORIG_DIR/petsc_v3.14.1_install/ex11_univ_egads_sphere.out $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1/src/dm/impls/plex/tests/output/ex11_univ_egads_sphere.out

# Install the modified petsc code
echo "Start Petsc Build"
cd ~/ubaceslab/build_scripts
./build_petsc_egads_version.sh

