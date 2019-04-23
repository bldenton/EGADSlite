#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
SDIR  = $(MAKEDIR)
UDIR  = $(SDIR)\..\util
LDIR  = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR  = $(ESP_BLOC)\obj
TDIR  = $(ESP_BLOC)\test
SFIL  = ..\test\vTess.c
!ELSE
ODIR  = .
TDIR  = $(ESP_ROOT)\bin
SFIL  = ..\examples\vTess.c
!ENDIF

$(TDIR)\liteVTess.exe:	$(ODIR)\liteVTess.obj $(LDIR)\egadslite.lib \
			$(LDIR)\wsserver.lib $(LDIR)\z.lib
	cl /Fe$(TDIR)\liteVTess.exe $(ODIR)\liteVTess.obj /link \
		/LIBPATH:$(LDIR) wsserver.lib z.lib egadslite.lib ws2_32.lib
	$(MCOMP) /manifest $(TDIR)\liteVTess.exe.manifest \
		/outputresource:$(TDIR)\liteVTess.exe;1

$(ODIR)\liteVTess.obj:	$(SFIL) $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
			$(IDIR)\egadsErrors.h $(IDIR)\wsss.h $(IDIR)\wsserver.h
	cl /c $(COPTS) $(DEFINE) -DLITE -I$(IDIR) -I$(IDIR)\winhelpers \
		$(SFIL) /Fo$(ODIR)\liteVTess.obj

clean:
	-del $(ODIR)\liteVTess.obj

cleanall:	clean
	-del $(TDIR)\liteVTess.exe $(TDIR)\liteVTess.exe.manifest
