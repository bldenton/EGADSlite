#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = $(ESP_BLOC)/obj
TDIR = $(ESP_BLOC)/test
SFIL = ../test/vTess.c
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
SFIL = ../examples/vTess.c
endif

$(TDIR)/liteVTess:	$(ODIR)/liteVTess.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/liteVTess $(ODIR)/liteVTess.o \
		-L$(LDIR) -lwsserver -legadslite -lpthread -lz $(RPATH) -lm

$(ODIR)/liteVTess.o:	$(SFIL) $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I$(IDIR) $(SFIL) \
		-o $(ODIR)/liteVTess.o

clean:
	-rm $(ODIR)/liteTess.o

cleanall:	clean
	-rm $(TDIR)/liteTess
