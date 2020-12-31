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

OBJS  = liteBase.o liteMemory.o liteGeom.o liteTopo.o liteAttrs.o liteImport.o \
	liteString.o

ifeq ($(ESP_ARCH),DARWIN64)
default:	$(LDIR)/libegadslite.dylib $(LDIR)/libegadslite.a
else
default:	$(LDIR)/libegadslite.so    $(LDIR)/libegadslite.a \
			$(IDIR)/cudaUtil_lite.h  $(IDIR)/egads_lite.h  $(IDIR)/egadsErrors_lite.h \
			$(IDIR)/egadsInternals_lite.h  $(IDIR)/egadsTris_lite.h  $(IDIR)/egadsTypes_lite.h \
			$(IDIR)/emp_lite.h  $(IDIR)/liteClasses.h  $(IDIR)/liteDevice.h \
			$(IDIR)/liteString.h  $(IDIR)/regQuads_lite.h  $(IDIR)/$(ESP_ARCH)
endif

$(LDIR)/libegadslite.so:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobustLite.o empLite.o evaLite.o \
				ratLite.o liteRegQuads.o
	touch $(LDIR)/libegadslite.so
	rm $(LDIR)/libegadslite.so
	(cd $(ODIR); $(CC) -shared -Wl,-no-undefined \
		-o $(LDIR)/libegadslite.so $(OBJS) liteTess.o liteTris.o \
		liteQuads.o liteTessInp.o egadsRobustLite.o empLite.o evaLite.o \
		ratLite.o liteRegQuads.o -lpthread -lm )

$(LDIR)/libegadslite.dylib:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobustLite.o empLite.o evaLite.o \
				ratLite.o liteRegQuads.o
	touch $(LDIR)/libegadslite.dylib
	rm $(LDIR)/libegadslite.dylib
	(cd $(ODIR); $(CC) -dynamiclib -o $(LDIR)/libegadslite.dylib \
                $(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o \
		egadsRobustLite.o empLite.o evaLite.o ratLite.o liteRegQuads.o \
		-undefined error -install_name '@rpath/libegadslite.dylib' \
                -current_version $(EGREV) )

$(LDIR)/libegadslite.a:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobustLite.o empLite.o evaLite.o \
				ratLite.o liteRegQuads.o
	touch $(LDIR)/libegadslite.a
	rm $(LDIR)/libegadslite.a
	(cd $(ODIR); ar $(LOPTS) $(LDIR)/libegadslite.a $(OBJS) \
                liteTess.o liteTris.o liteQuads.o liteTessInp.o egadsRobustLite.o \
		empLite.o evaLite.o ratLite.o liteRegQuads.o ; $(RANLB) )

$(OBJS): %.o:	%.c ../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
		../include/egadsTypes_lite.h ../include/liteClasses.h ../include/liteDevice.h
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -I../include -I. -I../src $< \
		-o $(ODIR)/$@

liteTess.o:	../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
                ../include/egadsTypes_lite.h ../include/liteClasses.h ../src/egadsTess_lite.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -I../include \
                ../src/egadsTess_lite.c -o $(ODIR)/liteTess.o

liteTris.o:	../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
                ../include/egadsTypes_lite.h ../include/liteClasses.h ../src/egadsTris_lite.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -I../include \
                ../src/egadsTris_lite.c -o $(ODIR)/liteTris.o

liteQuads.o:	../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
                ../include/egadsTypes_lite.h ../include/liteClasses.h ../src/egadsQuads_lite.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -I../include \
                ../src/egadsQuads_lite.c -o $(ODIR)/liteQuads.o

liteTessInp.o:	../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
                ../include/egadsTypes_lite.h ../include/liteClasses.h ../src/egadsTessInp_lite.c
	$(CC) -c $(COPTS) $(DEFANL) $(DEFINE) -DLITE -I../include \
                ../src/egadsTessInp_lite.c -o $(ODIR)/liteTessInp.o

empLite.o:		../src/emp_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -D$(ESP_ARCH) -I../include ../src/emp_lite.c \
		-o $(ODIR)/empLite.o

evaLite.o:	../src/evaluate_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include \
		../src/evaluate_lite.c -o $(ODIR)/evaLite.o

ratLite.o:	../src/rational_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE ../src/rational_lite.c \
		-o $(ODIR)/ratLite.o

liteRegQuads.o:	../src/regQuads_lite.c ../include/regQuads_lite.h
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include \
		../src/regQuads_lite.c -o $(ODIR)/liteRegQuads.o

egadsRobustLite.o:	../src/egadsRobust_lite.c
	$(CC) -c $(COPTS) $(DEFINE) ../src/egadsRobust_lite.c \
		-o $(ODIR)/egadsRobustLite.o

$(IDIR)/cudaUtil_lite.h:	../include/cudaUtil_lite.h
	touch $(IDIR)/cudaUtil_lite.h
	rm $(IDIR)/cudaUtil_lite.h
	cp -p ../include/cudaUtil_lite.h $(IDIR)/cudaUtil_lite.h

$(IDIR)/egads_lite.h:	../include/egads_lite.h
	touch $(IDIR)/egads_lite.h
	rm $(IDIR)/egads_lite.h
	cp -p ../include/egads_lite.h $(IDIR)/egads_lite.h

$(IDIR)/egadsErrors_lite.h:	../include/egadsErrors_lite.h
	touch $(IDIR)/egadsErrors_lite.h
	rm $(IDIR)/egadsErrors_lite.h
	cp -p ../include/egadsErrors_lite.h $(IDIR)/egadsErrors_lite.h

$(IDIR)/egadsInternals_lite.h:	../include/egadsInternals_lite.h
	touch $(IDIR)/egadsInternals_lite.h
	rm $(IDIR)/egadsInternals_lite.h
	cp -p ../include/egadsInternals_lite.h $(IDIR)/egadsInternals_lite.h

$(IDIR)/egadsTris_lite.h:	../include/egadsTris_lite.h
	touch $(IDIR)/egadsTris_lite.h
	rm $(IDIR)/egadsTris_lite.h
	cp -p ../include/egadsTris_lite.h $(IDIR)/egadsTris_lite.h

$(IDIR)/egadsTypes_lite.h:	../include/egadsTypes_lite.h
	touch $(IDIR)/egadsTypes_lite.h
	rm $(IDIR)/egadsTypes_lite.h
	cp -p ../include/egadsTypes_lite.h $(IDIR)/egadsTypes_lite.h

$(IDIR)/emp_lite.h:	../include/emp_lite.h
	touch $(IDIR)/emp_lite.h
	rm $(IDIR)/emp_lite.h
	cp -p ../include/emp_lite.h $(IDIR)/emp_lite.h

$(IDIR)/liteClasses.h:	../include/liteClasses.h
	touch $(IDIR)/liteClasses.h
	rm $(IDIR)/liteClasses.h
	cp -p ../include/liteClasses.h $(IDIR)/liteClasses.h

$(IDIR)/liteDevice.h:	../include/liteDevice.h
	touch $(IDIR)/liteDevice.h
	rm $(IDIR)/liteDevice.h
	cp -p ../include/liteDevice.h $(IDIR)/liteDevice.h

$(IDIR)/liteString.h:	../include/liteString.h
	touch $(IDIR)/liteString.h
	rm $(IDIR)/liteString.h
	cp -p ../include/liteString.h $(IDIR)/liteString.h

$(IDIR)/regQuads_lite.h:	../include/regQuads_lite.h
	touch $(IDIR)/regQuads_lite.h
	rm $(IDIR)/regQuads_lite.h
	cp -p ../include/regQuads_lite.h $(IDIR)/regQuads_lite.h

$(IDIR)/$(ESP_ARCH):	../include/$(ESP_ARCH)
	touch $(IDIR)/$(ESP_ARCH)
	rm $(IDIR)/$(ESP_ARCH)
	cp -p ../include/$(ESP_ARCH) $(IDIR)/$(ESP_ARCH)

clean:
	-(cd $(ODIR); rm $(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o evaLite.o ratLite.o liteRegQuads.o )

cleanall:	clean
	touch $(LDIR)/libegadslite.dylib $(LDIR)/libegadslite.so
	-rm $(LDIR)/libegadslite*