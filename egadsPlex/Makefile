#EGADS_LIBDIR = /PETSc3/cad/EGADSlite/src
#EGADS_LIBDIR2= ../src
EGADS_LIBDIR = $(ESP_ROOT)/lib
#EGADS_INCLUDEDIR = $(ESP_ROOT)/include
#CFLAGS = -I../include

# Note: I replaced -legadslite with -legads

testPlex: egadsPlex.o
	${CLINKER} -o $@ $< -legadslite ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -L$(EGADS_LIBDIR2)

testNozzle: egadsNozzle.o
	${CLINKER} -o $@ $< -legads ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -L$(EGADS_LIBDIR2)

testNozzleNew: egadsNozzle_new.o
	${CLINKER} -o $@ $< -legads ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR}

test_egads_examples: egads_example.o
	${CLINKER} -o $@ $< -legadslite ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -L$(EGADS_LIBDIR2)
	
testTess: egadsTess.o
	${CLINKER} -o $@ $< -legadslite ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -L$(EGADS_LIBDIR2)

clean::
	${RM} testPlex egadsPlex.o egadsNozzle.o testNozzle egadsTess.o testTess testNozzleNew egadsNozzle_new.o test_egads_examples egads_example.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
