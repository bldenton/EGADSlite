#
# allow for clang static analyzer
CCSAV = $(CC)
CXSAV = $(CXX)
include ../include/$(ESP_ARCH)
ifneq (,$(findstring analyzer,$(CCSAV)))
CC     = $(CCSAV)
CXX    = $(CXSAV)
endif

LDIR  = $(ESP_ROOT)/lib
IDIR  = $(ESP_ROOT)/include
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

VPATH = $(ODIR)

OBJS  = liteBase.o liteMemory.o liteAttrs.o liteImport.o liteString.o


ifeq ($(ESP_ARCH),DARWIN64)
default:	$(LDIR)/libegadsliteNR.dylib
else
default:	$(LDIR)/libegadsliteNR.so
endif

$(LDIR)/libegadsliteNR.so:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o \
				ratLite.o liteRegQuads.o \
				evaluateNR.o liteGeomNR.o liteTopoNR.o
	touch $(LDIR)/libegadsliteNR.so
	rm $(LDIR)/libegadsliteNR.so
	(cd $(ODIR); $(CC) -shared -Wl,-no-undefined \
		-o $(LDIR)/libegadsliteNR.so $(OBJS) liteTess.o liteTris.o \
		liteQuads.o liteTessInp.o egadsRobust.o emp.o ratLite.o \
		liteRegQuads.o evaluateNR.o liteGeomNR.o \
		liteTopoNR.o -lpthread -lm )

$(LDIR)/libegadsliteNR.dylib:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o ratLite.o \
				liteRegQuads.o evaluateNR.o \
				liteGeomNR.o liteTopoNR.o
	touch $(LDIR)/libegadsliteNR.dylib
	rm $(LDIR)/libegadsliteNR.dylib
	(cd $(ODIR); $(CC) -dynamiclib -o $(LDIR)/libegadsliteNR.dylib \
                $(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o \
		egadsRobust.o emp.o ratLite.o liteRegQuads.o \
		evaluateNR.o liteGeomNR.o liteTopoNR.o \
		-undefined error -install_name '@rpath/libegadsliteNR.dylib' \
                -current_version $(EGREV) )

$(OBJS): %.o:	%.c ../include/egadsErrors.h ../src/egadsInternals.h \
		../include/egadsTypes.h liteClasses.h
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -I../include -I. -I../src $< \
		-o $(ODIR)/$@

evaluateNR.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h ../util/evaluateNR.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -I../include -I. -I../src \
                ../util/evaluateNR.c -o $(ODIR)/evaluateNR.o

liteGeomNR.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h liteGeom.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -DCUDA -I../include -I. \
		-I../src liteGeom.c -o $(ODIR)/liteGeomNR.o

liteTopoNR.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h liteTopo.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -DCUDA -I../include -I. \
		-I../src liteTopo.c -o $(ODIR)/liteTopoNR.o

liteTris.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h ../src/egadsTris.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -I../include -I. -I../src \
                ../src/egadsTris.c -o $(ODIR)/liteTris.o

liteQuads.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h ../src/egadsQuads.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -I../include -I. -I../src \
                ../src/egadsQuads.c -o $(ODIR)/liteQuads.o

liteTessInp.o:	../include/egadsErrors.h ../src/egadsInternals.h \
                ../include/egadsTypes.h liteClasses.h ../src/egadsTessInp.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -I../include -I. -I../src \
                -I../util ../src/egadsTessInp.c -o $(ODIR)/liteTessInp.o

emp.o:		../util/emp.c
	$(CC) -c $(COPTS) $(DEFINE) -D$(ESP_ARCH) -I../include ../util/emp.c \
		-o $(ODIR)/emp.o

ratLite.o:	../util/rational.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE ../util/rational.c \
		-o $(ODIR)/ratLite.o

liteRegQuads.o:	../util/regQuads.c ../util/regQuads.h
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I../util \
		../util/regQuads.c -o $(ODIR)/liteRegQuads.o

egadsRobust.o:	../src/egadsRobust.c
	$(CC) -c $(COPTS) $(DEFINE) ../src/egadsRobust.c \
		-o $(ODIR)/egadsRobust.o

clean:
	-(cd $(ODIR); rm $(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o evaLite.o ratLite.o liteRegQuads.o)

cleanall:	clean
	touch $(LDIR)/libegadslite.dylib $(LDIR)/libegadslite.so
	-rm $(LDIR)/libegadslite*
