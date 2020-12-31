#!/bin/bash
UPDATE_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/ESP_v1.18/EngSketchPad/src/EGADS
DESTIN_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/EGADSlite

### Update EGADSlite Source Code ###
echo ""
echo "[Updating EGADSlite version for Petsc]"

# Update include folder
echo "    [Updating EGADSlite /include/ folder]"
cp $UPDATE_DIR/include/egads.h $DESTIN_DIR/include/egads_lite.h
cp $UPDATE_DIR/include/egadsErrors.h $DESTIN_DIR/include/egadsErrors_lite.h
cp $UPDATE_DIR/src/egadsInternals.h $DESTIN_DIR/include/egadsInternals_lite.h
cp $UPDATE_DIR/src/egadsTris.h $DESTIN_DIR/include/egadsTris_lite.h
cp $UPDATE_DIR/include/egadsTypes.h $DESTIN_DIR/include/egadsTypes_lite.h
cp $UPDATE_DIR/include/emp.h $DESTIN_DIR/include/emp_lite.h
cp $UPDATE_DIR/lite/liteClasses.h $DESTIN_DIR/include/liteClasses.h
cp $UPDATE_DIR/lite/cudaUtil.h $DESTIN_DIR/include/cudaUtil_lite.h
cp $UPDATE_DIR/lite/liteDevice.h $DESTIN_DIR/include/liteDevice.h
cp $UPDATE_DIR/lite/liteString.h $DESTIN_DIR/include/liteString.h
cp $UPDATE_DIR/util/regQuads.h $DESTIN_DIR/include/regQuads_lite.h
cp $UPDATE_DIR/include/DARWIN64 $DESTIN_DIR/include/DARWIN64
cp $UPDATE_DIR/include/LINUX64 $DESTIN_DIR/include/LINUX64
cp $UPDATE_DIR/include/WIN64.2015 $DESTIN_DIR/include/WIN64.2015
cp $UPDATE_DIR/include/WIN64.2017 $DESTIN_DIR/include/WIN64.2017
cp $UPDATE_DIR/include/WIN64.2019 $DESTIN_DIR/include/WIN64.2019

# Update src folder
echo "    [Updating EGADSlite /src/ folder]"
cp $UPDATE_DIR/lite/Makefile $DESTIN_DIR/src/Makefile
cp $UPDATE_DIR/lite/NMakefile $DESTIN_DIR/src/NMakefile
cp $UPDATE_DIR/src/egadsQuads.c $DESTIN_DIR/src/egadsQuads_lite.c
cp $UPDATE_DIR/src/egadsRobust.c $DESTIN_DIR/src/egadsRobust_lite.c
cp $UPDATE_DIR/src/egadsTess.c $DESTIN_DIR/src/egadsTess_lite.c
cp $UPDATE_DIR/src/egadsTessInp.c $DESTIN_DIR/src/egadsTessInp_lite.c
cp $UPDATE_DIR/src/egadsTris.c $DESTIN_DIR/src/egadsTris_lite.c
cp $UPDATE_DIR/lite/egadslite.def $DESTIN_DIR/src/egadslite.def
cp $UPDATE_DIR/lite/egadslite.rc $DESTIN_DIR/src/egadslite.rc
cp $UPDATE_DIR/util/emp.c $DESTIN_DIR/src/emp_lite.c
cp $UPDATE_DIR/util/evaluate.c $DESTIN_DIR/src/evaluate_lite.c
cp $UPDATE_DIR/lite/liteAttrs.c $DESTIN_DIR/src/liteAttrs.c
cp $UPDATE_DIR/lite/liteBase.c $DESTIN_DIR/src/liteBase.c
cp $UPDATE_DIR/lite/liteGeom.c $DESTIN_DIR/src/liteGeom.c
cp $UPDATE_DIR/lite/liteImport.c $DESTIN_DIR/src/liteImport.c
cp $UPDATE_DIR/lite/liteMemory.c $DESTIN_DIR/src/liteMemory.c
cp $UPDATE_DIR/lite/liteTess.mak $DESTIN_DIR/src/liteTess.mak
cp $UPDATE_DIR/lite/liteTess.make $DESTIN_DIR/src/liteTess.make
cp $UPDATE_DIR/lite/liteTest.c $DESTIN_DIR/src/liteTest.c
cp $UPDATE_DIR/lite/liteTest.mak $DESTIN_DIR/src/liteTest.mak
cp $UPDATE_DIR/lite/liteTest.make $DESTIN_DIR/src/liteTest.make
cp $UPDATE_DIR/lite/liteTopo.c $DESTIN_DIR/src/liteTopo.c
cp $UPDATE_DIR/util/rational.c $DESTIN_DIR/src/rational_lite.c
cp $UPDATE_DIR/lite/NR.make $DESTIN_DIR/src/NR.make
cp $UPDATE_DIR/lite/liteDevice.c $DESTIN_DIR/src/liteDevice.c
cp $UPDATE_DIR/lite/liteString.c $DESTIN_DIR/src/liteString.c
cp $UPDATE_DIR/lite/liteTessNR.make $DESTIN_DIR/src/liteTessNR.make

# Update liteBase.c with modified version containing EG_getInfo() -- For EGADS v1.18 only
echo "    [Updating EGADSlite /src/liteBase.c file with EG_getInfo() function]"
cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADS-dev/special-litebase.c/liteBase.c $DESTIN_DIR/src/liteBase.c

# Previously in util folder -- now moved to src folder for Petsc EGADSlite v1.18
echo "    [Updating EGADSlite /src/ folder with files originally stored in /util/ folder]"
cp $UPDATE_DIR/util/README.txt $DESTIN_DIR/src/README.txt
cp $UPDATE_DIR/util/SurrealD1_btest.cpp $DESTIN_DIR/src/SurrealD1_btest_lite.cpp
cp $UPDATE_DIR/util/SurrealD4_btest.cpp $DESTIN_DIR/src/SurrealD4_btest_lite.cpp
cp $UPDATE_DIR/util/SurrealS1_btest.cpp $DESTIN_DIR/src/SurrealS1_btest_lite.cpp
cp $UPDATE_DIR/util/SurrealS4_btest.cpp $DESTIN_DIR/src/SurrealS4_btest_lite.cpp
cp $UPDATE_DIR/util/ThreadTest.c $DESTIN_DIR/src/ThreadTest_lite.c
cp $UPDATE_DIR/util/emp.c $DESTIN_DIR/src/emp_lite.c
cp $UPDATE_DIR/util/evaluate.c $DESTIN_DIR/src/evaluate_lite.c
cp $UPDATE_DIR/util/limitTessBody.c $DESTIN_DIR/src/limitTessBody_lite.c
cp $UPDATE_DIR/util/limits.mak $DESTIN_DIR/src/limits_lite.mak
cp $UPDATE_DIR/util/limits.make $DESTIN_DIR/src/limits_lite.make
cp $UPDATE_DIR/util/rational.c $DESTIN_DIR/src/rational_lite.c
cp $UPDATE_DIR/util/retessFaces.c $DESTIN_DIR/src/retessFaces_lite.c
cp $UPDATE_DIR/util/regQuads.c $DESTIN_DIR/src/regQuads_lite.c

# Update util folder - All files have been moved to either include or src or Petsc integration
# They continue to be updated here as backup
#cp $UPDATE_DIR/util/README.txt $DESTIN_DIR/util/README.txt
#cp $UPDATE_DIR/util/SurrealD1_btest.cpp $DESTIN_DIR/util/SurrealD1_btest.cpp
#cp $UPDATE_DIR/util/SurrealD4_btest.cpp $DESTIN_DIR/util/SurrealD4_btest.cpp
#cp $UPDATE_DIR/util/SurrealS1_btest.cpp $DESTIN_DIR/util/SurrealS1_btest.cpp
#cp $UPDATE_DIR/util/SurrealS4_btest.cpp $DESTIN_DIR/util/SurrealS4_btest.cpp
#cp $UPDATE_DIR/util/ThreadTest.c $DESTIN_DIR/util/ThreadTest.c
#cp $UPDATE_DIR/util/emp.c $DESTIN_DIR/util/emp.c
#cp $UPDATE_DIR/util/evaluate.c $DESTIN_DIR/util/evaluate.c
#cp $UPDATE_DIR/util/limitTessBody.c $DESTIN_DIR/util/limitTessBody.c
#cp $UPDATE_DIR/util/limits.mak $DESTIN_DIR/util/limits.mak
#cp $UPDATE_DIR/util/limits.make $DESTIN_DIR/util/limits.make
#cp $UPDATE_DIR/util/rational.c $DESTIN_DIR/util/rational.c
#cp $UPDATE_DIR/util/retessFaces.c $DESTIN_DIR/util/retessFaces.c
#cp $UPDATE_DIR/util/regQuads.h $DESTIN_DIR/util/regQuads.h
#cp $UPDATE_DIR/util/regQuads.c $DESTIN_DIR/util/regQuads.c

# Update Variables, Parameters and Function Names for a independent EGADSlite library
echo "    [Updating EGADSlite Source files with unique Variables, Parameters, Function & File Names to work alongside EGADS]"
for f in $DESTIN_DIR/*/*; do sed -i 's/EGADS_H/EGADS_LITE_H/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/CUDA_UTIL_H/CUDA_UTIL_LITE_H/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/EGADSERRORS_H/EGADSERRORS_LITE_H/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/EGADSINTERNALS_H/EGADSINTERNALS_LITE_H/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/EGADSTYPES_H/EGADSTYPES_LITE_H/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/EMP_/EMP_LITE_/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/_EGADS_STRING_H/_EGADS_STRING_LITE_H/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/EG_/EGlite_/g' $f; done;

for f in $DESTIN_DIR/*/*; do sed -i 's/cudaUtil.h/cudaUtil_lite.h/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/egads.h/egads_lite.h/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/egadsErrors.h/egadsErrors_lite.h/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/egadsInternals.h/egadsInternals_lite.h/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/egadsTris.h/egadsTris_lite.h/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/egadsTypes.h/egadsTypes_lite.h/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/emp.h/emp_lite.h/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/regQuads.h/regQuads_lite.h/g' $f; done;

for f in $DESTIN_DIR/*/*; do sed -i 's/egadsQuads.c/egadsQuads_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/egadsRobust.c/egadsRobust_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/egadsTess.c/egadsTess_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/egadsTessInp.c/egadsTessInp_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/egadsTris.c/egadsTris_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/emp.c/emp_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/evaluate.c/evaluate_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/limits.mak/limits_lite.mak/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/limits.make/limits_lite.make/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/limitTessBody.c/limitTessBody_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/rational.c/rational_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/regQuads.c/regQuads_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/retessFaces.c/retessFaces_lite.c/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/SurrealD1_btest.cpp/SurrealD1_btest_lite.cpp/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/SurrealD4_btest.cpp/SurrealD4_btest_lite.cpp/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/SurrealS1_btest.cpp/SurrealS1_btest_lite.cpp/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/SurrealS4_btest.cpp/SurrealS4_btest_lite.cpp/g' $f; done;
for f in $DESTIN_DIR/*/*; do sed -i 's/ThreadTest.c/ThreadTest_lite.c/g' $f; done;

echo "[COMPLETED EGADSlite Codebase Update - Need to Still Modify MakeFile]"
echo ""

### Update EGADS Source Code ###
echo "[Updating EGADSlite version for Petsc]"

