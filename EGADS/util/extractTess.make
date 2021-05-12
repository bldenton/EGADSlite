#
include ../include/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = $(ESP_BLOC)/obj
TDIR = $(ESP_BLOC)/test
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
endif

default:	$(TDIR)/extractTess

$(TDIR)/extractTess:	$(ODIR)/tessExtract.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/extractTess $(ODIR)/tessExtract.o \
		-L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/tessExtract.o:	extractTess.c ../src/egadsTris.h ../include/egads.h \
			../include/egadsTypes.h ../include/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -DSTANDALONE -DREPORT \
		extractTess.c -o $(ODIR)/tessExtract.o

clean:
	-rm $(ODIR)/tessExtract.o

cleanall:	clean
	-rm $(TDIR)/extractTess
