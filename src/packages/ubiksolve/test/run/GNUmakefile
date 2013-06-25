#
# Author: John A. Turner
#         e-mail: john.turner@pobox.com
#         http://john.turner.org/
#
# This makefile requires GNU make, available for free from
# the Free Software Foundation, Inc. (http://www.gnu.org/).
# It could probably be modified for use with other versions of make.

LIBNAME = UbikSolve
PROG    = UbikTest

# Directories.

dot := $(shell pwd)
include ../../config.mak
makeutilsdir := $(dot)/../../$(makeutilsdir)

# Info about environment.

include $(makeutilsdir)/env.mak

# Fortran setup (needed only to build the $(platform) variable)

include $(makeutilsdir)/f90.mak

default:
	@$(ECHO) "To check compile settings:"
	@$(ECHO) "  make check"
	@$(ECHO) "edit" $(subst /test,,$(dot))/config.mak "to override"
	@$(ECHO)
	@$(ECHO) "To run a subset of the test problems and display a one-line"
	@$(ECHO) "summary for each problem:"
	@$(ECHO) "  make test              for the optimized version"
	@$(ECHO)
	@$(ECHO) "To run all the test problems and display a one-line summary for"
	@$(ECHO) "each problem:"
	@$(ECHO) "  make testall           for the optimized version"
	@$(ECHO)
	@$(ECHO) "To view the results of the last run of the test suite:"
	@$(ECHO) "  make summary"
	@$(ECHO)
	@$(ECHO) "To run problems of increasing size to test performance:"
	@$(ECHO) "  make cputest"
	@$(ECHO)
	@$(ECHO) "To view the results of the last performance test:"
	@$(ECHO) "  make cpusummary"
	@$(ECHO)
	@$(ECHO) "To remove trash:"
	@$(ECHO) "  make cleanup"
	@$(ECHO)
	@$(ECHO) "To remove all binaries (libraries, object files, and executables):"
	@$(ECHO) "  make clean"
	@$(ECHO) "(files created by preprocessing are left)"
	@$(ECHO)
	@$(ECHO) "To remove everything that wasn't in the original distribution:"
	@$(ECHO) "  make distclean"

# Test suite targets.

test testall cputest summary cpusummary: ensure_rundir
	cd $(platform); $(MAKE) -f ../GNUmakefile.sub $@

check: ensure_rundir
	@$(ECHO)
	@$(ECHO) ">>>>>>> in" $(dot) "<<<<<<<<"
	@$(ECHO)
	@$(ECHO) "path to Makeutils: " $(makeutilsdir)
	@$(ECHO) "platform:          " $(platform)
	$(MAKE) --directory=$(platform) -f ../GNUmakefile.sub $@

ensure_rundir:
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
