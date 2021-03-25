#EGADS_LIBDIR = /PETSc3/cad/EGADSlite/src
#EGADS_LIBDIR2= ../src
#EGADS_INCLUDEDIR = $(ESP_ROOT)/include
#CFLAGS = -I../include

#//Use this when using indpendent EGADS libraries
#EGADS_LIBDIR = $(ESP_ROOT)/lib

# Note: I replaced -legadslite with -legads

testNozzleIges: egadsNozzle_new.o
	${CLINKER} -o $@ $< -legads ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR}

testNozzleStep: egadsNozzle_new.o
	${CLINKER} -o $@ $< -legads ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR}

testNozzleEgads: egadsNozzle_new.o
	${CLINKER} -o $@ $< -legads ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR}


clean::
	${RM} testNozzleIges testNozzleStep testNozzleEgads egadsNozzle_new.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
