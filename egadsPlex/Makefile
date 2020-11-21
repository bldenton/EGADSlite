#EGADS_LIBDIR = /PETSc3/cad/EGADSlite/src
#EGADS_LIBDIR2= ../src
#EGADS_LIBDIR2=$ESP_ROOT/lib
#CFLAGS = -I../include

# Note: I replaced -legadslite with -legads

testPlex: egadsPlex.o
	${CLINKER} -o $@ $< -legadslite ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -L$(EGADS_LIBDIR2)

testNozzle: egadsNozzle.o
	${CLINKER} -o $@ $< -legadslite ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -L$(EGADS_LIBDIR2)

testNozzleNew: egadsNozzle_new.o
	${CLINKER} -o $@ $< -legadslite ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -L$(EGADS_LIBDIR2)
	
testTess: egadsTess.o
	${CLINKER} -o $@ $< -legadslite ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -L$(EGADS_LIBDIR2)

clean::
	${RM} testPlex egadsPlex.o egadsNozzle.o testNozzle egadsTess.o testTess testNozzleNew egadsNozzle_new.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
