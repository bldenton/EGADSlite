#
!include ..\include\$(ESP_ARCH).$(MSVC)
IDIR = $(ESP_ROOT)\include
LDIR = $(ESP_ROOT)\lib
!ifdef ESP_BLOC
ODIR = $(ESP_BLOC)\obj
TDIR = $(ESP_BLOC)\test
!else
ODIR = .
TDIR = $(ESP_ROOT)\bin
!endif


$(TDIR)\extractTess.exe:	$(ODIR)\tessExtract.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\extractTess.exe $(ODIR)\tessExtract.obj \
		$(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\extractTess.exe.manifest \
		/outputresource:$(TDIR)\extractTess.exe;1

$(ODIR)\tessExtract.obj:	extractTess.c $(IDIR)\egadsTris.h \
		$(IDIR)\egads.h $(IDIR)\egadsTypes.h $(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -DSTANDALONE -DREPORT \
		extractTess.c /Fo$(ODIR)\tessExtract.obj

clean:
	-del $(ODIR)\tessExtract.obj

cleanall:	clean
	-del $(TDIR)\extractTess.exe
