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
include ../config.mak
makeutilsdir := $(dot)/../$(makeutilsdir)

# Info about environment.

include $(makeutilsdir)/env.mak

# Fortran setup (needed only to build the $(platform) variable)

include $(makeutilsdir)/f90.mak

trash = *~ .*~ core

# Force these targets to be made, even if the target seems up-to-date.

.PHONY: $(PROG)$(exe)

default:
	@$(ECHO) "To check compile settings:"
	@$(ECHO) "  make check"
	@$(ECHO) "edit" $(subst /test,,$(dot))/config.mak "to override"
	@$(ECHO)
	@$(ECHO) "To compile the test program and library:"
	@$(ECHO) "  make $(PROG)"
	@$(ECHO)
	@$(ECHO) "To run a subset of the test problems and display a one-line"
	@$(ECHO) "summary for each problem:"
	@$(ECHO) "  make test"
	@$(ECHO)
	@$(ECHO) "To run all the test problems and display a one-line summary for"
	@$(ECHO) "each problem:"
	@$(ECHO) "  make testall"
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

# Test program targets.

$(PROG):
	$(MAKE) --directory=build $@

# Test suite targets.

test cputest summary cpusummary:
	$(MAKE) --directory=run $@

# tests

check:
	@$(ECHO)
	@$(ECHO) ">>>>>>> in" $(dot) "<<<<<<<"
	@$(ECHO)
	@$(ECHO) "path to Makeutils: " $(makeutilsdir)
	@$(ECHO) "operating system:  " $(os)$(os_ver)
	@$(ECHO) "architecture:      " $(arch)$(arch_ver)
	$(MAKE) --directory=build $@
	$(MAKE) --directory=run $@

# Cleanup.

cleanup:
	$(MAKE) --directory=build $@
	$(MAKE) --directory=run $@
	$(RM) $(RMFLAGS) $(trash)

clean mostlyclean:
	$(MAKE) --directory=build $@
	$(MAKE) --directory=run $@
	$(RM) $(RMFLAGS) $(trash)

cleaner cleanest distclean realclean:
	$(MAKE) --directory=build $@
	$(MAKE) --directory=run $@
	$(RM) $(RMFLAGS) $(trash)

### Local Variables: ###
### tab-width: 8 ###
### End: ###
