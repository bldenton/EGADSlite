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

$(TDIR)/liteVTessNR:	$(ODIR)/liteVTess.o $(LDIR)/libwsserver.a
	$(CXX) -o $(TDIR)/liteVTessNR $(ODIR)/liteVTess.o \
		-L$(LDIR) -lwsserver -legadsliteNR -lpthread -lz $(RPATH) -lm

$(ODIR)/liteVTess.o:	$(SFIL) $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I$(IDIR) $(SFIL) \
		-o $(ODIR)/liteVTess.o

clean:
	-rm $(ODIR)/liteVTess.o

cleanall:	clean
	-rm $(TDIR)/liteVTessNR
