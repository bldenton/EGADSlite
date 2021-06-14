#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
TDIR  = $(ESP_BLOC)/test
else
ODIR  = .
TDIR  = .
endif

# allow for clang static analyzer
ifneq (,$(findstring analyzer,$(CC)))
DEFINE += -DDEBUG
endif
ifneq (,$(findstring intercept,$(CC)))
DEFINE += -DDEBUG
endif

VPATH = $(ODIR)

$(TDIR)/uvmap.exe:	$(ODIR)/uvmap.o $(LDIR)/libuvmap.a
	$(CC) -o $(TDIR)/uvmap.exe $(ODIR)/uvmap.o \
		-L$(LDIR) -luvmap $(RPATH) -lm

$(ODIR)/uvmap.o:	uvmap/main/uvmap.c 
	$(CC) -c $(COPTS) $(DEFINE) uvmap/main/uvmap.c -Iuvmap \
		-o $(ODIR)/uvmap.o

INC  = $(wildcard uvmap/*.h)
SRC  = $(wildcard uvmap/*.c)
OBJS = $(subst uvmap/,,$(SRC:.c=.o))

EG_%.o : uvmap/EG_%.c $(INC)
	$(CC) -c $(COPTS) $(DEFINE) -I. $< -o $(ODIR)/$@

uvmap_%.o : uvmap/uvmap_%.c $(INC)
	$(CC) -c $(COPTS) $(DEFINE) -I. $< -o $(ODIR)/$@

$(ODIR)/egadsUVmap.o:	egadsUVmap.c 
	$(CC) -c $(COPTS) $(DEFINE) egadsUVmap.c -Iuvmap \
		-o $(ODIR)/egadsUVmap.o

$(LDIR)/libuvmap.a:	 $(OBJS) $(ODIR)/egadsUVmap.o
	touch $(LDIR)/libuvmap.a
	rm $(LDIR)/libuvmap.a
	(cd $(ODIR); ar $(LOPTS) $(LDIR)/libuvmap.a $(ODIR)/egadsUVmap.o $(OBJS); $(RANLB) )

clean:
	-(cd $(ODIR); rm -f uvmap.o egadsUVmap.o $(OBJS))

cleanall:	clean
	-rm -f $(TDIR)/uvmap.exe $(LDIR)/libuvmap.a
