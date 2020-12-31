#!/bin/bash
ORIG_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD

# Unload petsc module if loaded
echo "[Unloading Petsc if Loaded]"
module unload petsc

# Go to build directory
cd $UBCESLAB_SWENV_PREFIX/builddir

# Remove current unpacked folder
echo "[Removing Previous Install]"
rm -rf petsc-3.14.1-clone

# Copy current cloned Petsc Repository
echo "[Copying Repository from /sourcesdir/ to /builddir/]"
cp -r $UBCESLAB_SWENV_PREFIX/sourcesdir/petsc/petsc $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1-clone
cp $ORIG_DIR/petsc_v3.14.1_install/plexegads_20201114.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1-clone/src/dm/impls/plex/plexegads.c
cp $ORIG_DIR/petsc_v3.14.1_install/egads_20201121.py $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1-clone/config/BuildSystem/config/packages/egads.py
#mkdir petsc-3.14.1-clone

# Install the modified petsc code
echo "[Start Petsc Build]"
cd ~/ubaceslab/build_scripts
./build_petsc_clone_egads_version.sh

# If you want to just update egads.py by itself
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/petsc_v3.14.1_install/egads_20201115.py $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1-clone/config/BuildSystem/config/packages/egads.py
