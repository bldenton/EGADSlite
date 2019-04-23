#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
LDIR  = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR  = $(ESP_BLOC)\obj
TDIR  = $(ESP_BLOC)\test
!ELSE
ODIR  = .
TDIR  = $(ESP_ROOT)\bin
!ENDIF


$(TDIR)\liteTest.exe:	$(ODIR)\liteTest.obj $(LDIR)\egadslite.lib
	cl /Fe$(TDIR)\liteTest.exe $(ODIR)\liteTest.obj $(LIBPTH) egadslite.lib
	$(MCOMP) /manifest $(TDIR)\liteTest.exe.manifest \
		/outputresource:$(TDIR)\liteTest.exe;1

$(ODIR)\liteTest.obj:	liteTest.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) liteTest.c /Fo$(ODIR)\liteTest.obj

clean:
	-del $(ODIR)\liteTest.obj

cleanall:	clean
	-rm $(TDIR)\liteTest.exe $(TDIR)\liteTest.exe.manifest
