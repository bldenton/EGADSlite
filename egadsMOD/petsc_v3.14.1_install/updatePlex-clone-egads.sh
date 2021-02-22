#!/bin/bash
#PETSC_ARCH=gcc-8.3.0-mpich-3.2.1-openblas-0.2.20-opt

echo "[Updating plexegads.c in Petsc v3.14.1]"
# Copy plexegads.c to correct location
ORIG_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADS-dev/egadsMOD
cp $ORIG_DIR/petsc_v3.14.1_install/plexegads_20210219.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1-clone/src/dm/impls/plex/plexegads.c
#cp $ORIG_DIR/petsc_v3.14.1_install/plexegads_20201114.c $UBCESLAB_SWENV_PREFIX/sourcesdir/petsc/petsc/src/dm/impls/plex/plexegads.c
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/plexegads.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/plexegads.c
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/plexegads_old_20200607.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/plexegads.c

# Copy ctetgenerate.c to correct location
#echo "Updating ctetgenerate.c in Petsc"
#echo ""
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/ctetgenerate.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/generators/ctetgen/ctetgenerate.c

#echo "Updating tetgenerate.cxx in Petsc"
#echo ""
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/tetgenerate_EGADs.cxx $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/generators/tetgen/tetgenerate.cxx
#cp $ORIG_DIR/petsc_v3.14.1_install/tetgenerate.cxx $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1-clone/src/dm/impls/plex/generators/tetgen/tetgenerate.cxx

# Copy updated egads.py file
echo "[Updating egads.py]"
cp $ORIG_DIR/petsc_v3.14.1_install/egads_20201121.py $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1-clone/config/BuildSystem/config/packages/egads.py

# Update Petsc
echo "[Updating Petsc]"
cd $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1-clone

make -f ./gmakefile
make -j8 libs
make -j8 install

cd $goEGADSlite

#--prefix=/mnt/c/Users/Brandon/software/libs/petsc/3.14.1-clone/gcc/8.3.0/mpich/3.2.1/openblas/0.2.20/opt
