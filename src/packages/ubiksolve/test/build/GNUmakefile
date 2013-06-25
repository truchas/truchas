LIBNAME = UbikSolve
PROG    = UbikTest

# Directories.

dot := $(shell pwd)
top = ../..

include $(top)/config.mak
makeutilsdir := $(dot)/../../$(makeutilsdir)
libbuilddir = $(top)/build

# info about environment

include $(makeutilsdir)/env.mak

# Fortran setup (needed only to build the $(platform) variable)

include $(makeutilsdir)/f90.mak

default:
	@echo "To check compile settings:"
	@echo "  make check"
	@echo "edit" $(subst /test/build,,$(dot))/config.mak "to override"
	@echo
	@echo "To compile the test program and library:"
	@echo "  make UbikTest"
	@echo
	@echo "To remove trash:"
	@echo "  make cleanup"
	@echo
	@echo "To remove all binaries (libraries, object files, and executables):"
	@echo "  make clean"
	@echo "(files created by preprocessing are left)"
	@echo
	@echo "To remove everything that wasn't in the original distribution:"
	@echo "  make distclean"

# Test program target.

$(PROG): ensure_builddir
	$(MAKE) --directory=$(libbuilddir) $(LIBNAME)
	$(MAKE) --directory=$(platform) -f ../GNUmakefile.sub $@

# Echo what commands and flags will be used for the build.
check: ensure_builddir
	@echo
	@echo ">>>>>>> in" $(dot) "<<<<<<<<"
	@echo
	@echo "path to Makeutils: " $(makeutilsdir)
	@echo "platform:          " $(platform)
	$(MAKE) --directory=$(platform) -f ../GNUmakefile.sub $@

ensure_builddir:
        # Create the output directory if it doesn't exist.
	@$(ECHO)
	if [ ! -d $(platform) ];  \
        then  \
         $(MKDIR) $(platform);  \
        fi

clean cleanup mostlyclean:
	@$(MAKE) --directory=$(platform) -f ../GNUmakefile.sub $@

cleaner:
	$(RM) -rf *~ $(platform)

cleanest distclean realclean:
	$(RM) -rf *~ *-*

### Local Variables: ###
### tab-width: 8 ###
### End: ###
