#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif

LDIR  = $(ESP_ROOT)/lib

.PHONY: all egadslite egads clean cleanall

all: egads egadslite
#	(cd EGADS/src; make)
#	(cd EGADSlite/src; make)
#	@echo " "
#	@echo " *** Build Completed! ***"

egadslite:
	(cd EGADSlite/src; make)
	@echo " "
	@echo " *** EGADSlite Build Completed! ***"
	@echo " "
	
egads:
	(cd EGADS/src; make)
	@echo " "
	@echo " *** EGADS Build Completed! ***"
	@echo " "

clean:
	(cd EGADSlite/src; make clean)
	(cd EGADS/src; make clean)

cleanall:
	(cd EGADSlite/src; make cleanall)
	(cd EGADS/src; make cleanall)
