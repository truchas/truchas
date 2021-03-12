VPATH=.:..:../../exodus

INC= -I$(NETCDF_MOD_DIR)
LIB= -L$(NETCDF_LIB_DIR) -lnetcdf

# Linux with the Lahey (default) or NAG compiler
ifeq (Linux, $(shell uname))
  ifeq (nag, $(compiler))
    F90= nagfor
    F90FLAGS= -g -gline -C=all $(INC) -DNAG
  else ifeq (nagpre52, $(compiler))
    F90= f95
    F90FLAGS= -g -gline -C=all $(INC) -DNAG
  else
    F90= lf95
    F90FLAGS= -g $(INC)
  endif
endif
                                                                                
PROGS= maptest

EXOBJ= string_utilities.o exodus_errors.o exodus_mesh_type.o \
       exodus_mesh_reader.o exodus_mesh_writer.o exodus_truchas_hack.o exodus.o

OBJ= $(EXOBJ) hpsort.o overlap_module.o gm_mesh_type.o grid_mapping_utils.o \
	grid_mapping_exodus.o grid_mapping_module.o driver.o

maptest: $(OBJ)
	$(F90) $(F90FLAGS) -o $@ $(OBJ) $(LIB)

CLEANFILES= *.o *.mod int_vols_file *.gmv
clean:; rm -f $(CLEANFILES)
cleanest: clean; rm -f $(PROGS)

.SUFFIXES: .F90
.F90.o:
	$(F90) $(F90FLAGS) -c $<
