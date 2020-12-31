#!/bin/bash
echo "[Updating plexegads.c in Petsc v3.14.1]"
echo ""
# Copy plexegads.c to correct location
ORIG_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD
cp $ORIG_DIR/petsc_v3.14.1_install/plexegads.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.14.1-clone/src/dm/impls/plex/plexegads.c
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/plexegads.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/plexegads.c
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/plexegads_old_20200607.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/plexegads.c

# Copy ctetgenerate.c to correct location
#echo "Updating ctetgenerate.c in Petsc"
#echo ""
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/ctetgenerate.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/generators/ctetgen/ctetgenerate.c

#echo "Updating tetgenerate.cxx in Petsc"
#echo ""
#cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/tetgenerate_EGADs.cxx $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/generators/tetgen/tetgenerate.cxx

# Update Petsc
echo "[Updating Petsc]"
cd /mnt/c/Users/Brandon/software/builddir/petsc-3.14.1-clone
make -f ./gmakefile
make install

cd $goEGADSlite
