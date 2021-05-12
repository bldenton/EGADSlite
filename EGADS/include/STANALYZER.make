
ifeq ("$(ESP_ARCH)","DARWIN64")
# Python version of scan-build
#
# Install with: pip install scan-build
#
SCANBUILD=intercept-build --override-compiler --cdb=$(SCANDIR)compile_commands.json $(MAKE) CC=intercept-cc CXX=intercept-c++ && analyze-build --cdb=$(SCANDIR)compile_commands.json -o $(SCANDIR) $(SCANEXCLUDE)

.PHONY: scan-build
scan-build:
	@bash -c '[ -d $(SCANDIR) ] && rm -rf $(SCANDIR)' || true
	$(SCANBUILD) --status-bugs

.PHONY: scan-view
scan-view:
	@bash -c '[ -d $(SCANDIR) ] && rm -rf $(SCANDIR)' || true
	$(SCANBUILD)
	-@open $(SCANDIR)/*/index.html || true
endif

ifeq ("$(ESP_ARCH)","LINUX64")
# Perl version of scan-build
#
# Install on ubuntu: sudo apt-get install clang-tools
#
.PHONY: scan-build
scan-build:
	@bash -c '[ -d $(SCANDIR) ] && rm -rf $(SCANDIR)' || true
	scan-build -o $(SCANDIR) $(SCANEXCLUDE) --status-bugs $(MAKE)

.PHONY: scan-view
scan-view:
	@bash -c '[ -d $(SCANDIR) ] && rm -rf $(SCANDIR)' || true
	scan-build -o $(SCANDIR) $(SCANEXCLUDE) $(MAKE)
	-@xdg-open $(SCANDIR)/*/index.html || true
endif
