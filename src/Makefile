#
# allow for clang static analyzer
CCSAV = $(CC)
CXSAV = $(CXX)
#include ../include/$(ESP_ARCH)
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

OBJS  = liteBase.o liteMemory.o liteGeom.o liteTopo.o liteAttrs.o liteImport.o


ifeq ($(ESP_ARCH),DARWIN64)
default:	$(LDIR)/libegadslite.dylib $(LDIR)/libegadslitestatic.a
else
default:	$(LDIR)/libegadslite.so    $(LDIR)/libegadslitestatic.a
endif

$(LDIR)/libegadslite.so:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o evaLite.o \
				rational.o
	touch $(LDIR)/libegadslite.so
	rm $(LDIR)/libegadslite.so
	(cd $(ODIR); $(CC) -shared -Wl,-no-undefined -o $(LDIR)/libegadslite.so \
		$(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o \
		egadsRobust.o emp.o evaLite.o rational.o -lpthread -lm )

$(LDIR)/libegadslite.dylib:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o evaLite.o \
				rational.o
	touch $(LDIR)/libegadslite.dylib
	rm $(LDIR)/libegadslite.dylib
	(cd $(ODIR); $(CC) -dynamiclib -o $(LDIR)/libegadslite.dylib \
                $(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o \
		egadsRobust.o emp.o evaLite.o rational.o \
		-undefined error \
                -install_name '@rpath/libegadslite.dylib' \
                -current_version $(EGREV) )

$(LDIR)/libegadslitestatic.a:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o evaLite.o \
				rational.o
	touch $(LDIR)/libegadslitestatic.a
	rm $(LDIR)/libegadslitestatic.a
	(cd $(ODIR); ar $(LOPTS) $(LDIR)/libegadslitestatic.a $(OBJS) \
                liteTess.o liteTris.o liteQuads.o liteTessInp.o egadsRobust.o \
		emp.o evaLite.o rational.o; $(RANLB) )

$(OBJS): %.o:	%.c ../include/egadsErrors.h ../include/egadsInternals.h \
		../include/egadsTypes.h ../include/liteClasses.h
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -I../include -I. -I../src $< \
		-o $(ODIR)/$@

liteTess.o:	../include/egadsErrors.h ../include/egadsInternals.h \
                ../include/egadsTypes.h ../include/liteClasses.h ../src/egadsTess.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -I../include -I. -I../src \
                ../src/egadsTess.c -o $(ODIR)/liteTess.o

liteTris.o:	../include/egadsErrors.h ../include/egadsInternals.h \
                ../include/egadsTypes.h ../include/liteClasses.h ../src/egadsTris.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -I../include -I. -I../src -I$(PETSC_DIR)/include \
                ../src/egadsTris.c -o $(ODIR)/liteTris.o

liteQuads.o:	../include/egadsErrors.h ../include/egadsInternals.h \
                ../include/egadsTypes.h ../include/liteClasses.h ../src/egadsQuads.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -I../include -I. -I../src \
                ../src/egadsQuads.c -o $(ODIR)/liteQuads.o

liteTessInp.o:	../include/egadsErrors.h ../include/egadsInternals.h \
                ../include/egadsTypes.h ../include/liteClasses.h ../src/egadsTessInp.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -I../include -I. -I../src \
                ../src/egadsTessInp.c -o $(ODIR)/liteTessInp.o

emp.o:		../util/emp.c
	$(CC) -c $(COPTS) $(DEFINE) -D$(ESP_ARCH) -I../include ../util/emp.c \
		-o $(ODIR)/emp.o

evaLite.o:	../util/evaluate.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
		../util/evaluate.c -o $(ODIR)/evaLite.o

rational.o:	../util/rational.c
	$(CC) -c $(COPTS) $(DEFINE) ../util/rational.c -o $(ODIR)/rational.o

egadsRobust.o:	../src/egadsRobust.c
	$(CC) -c $(COPTS) $(DEFINE) ../src/egadsRobust.c -o $(ODIR)/egadsRobust.o

clean:
	-(cd $(ODIR); rm $(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o)

cleanall:	clean
	touch $(LDIR)/libegadslite.dylib $(LDIR)/libegadslite.so
	-rm $(LDIR)/libegadslite*
