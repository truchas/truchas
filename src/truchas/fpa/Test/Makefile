VPATH=.:..:../../utilities

# Linux with the Lahey compiler
ifeq (Linux, $(shell uname))
  ifeq (lahey, $(compiler))
    F90= lf95
    F90FLAGS= -g $(INC)
  else
    F90= f95
    F90FLAGS= -g -gline $(INC)
  endif
endif

INC= -I.
OBJ= f90_assert.o fpa_kinds.o fixed_point_accelerator.o fpa_unit_test.o
PROGS= unit_test

unit_test: $(OBJ)
	$(F90) $(F90FLAGS) -o $@ $(OBJ)

CLEANFILES = *.o *.mod *.MOD
clean:; rm -rf $(CLEANFILES)
cleanest: clean; rm -f $(PROGS)

.SUFFIXES: .f90 .F90
.f90.o:
	$(F90) $(F90FLAGS) -c $<
.F90.o:
	$(F90) $(F90FLAGS) -c $<
