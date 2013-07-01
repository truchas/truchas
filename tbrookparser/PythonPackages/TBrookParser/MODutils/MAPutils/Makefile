TRUCHAS_SOURCE= ../../../../../src
VPATH= .:$(TRUCHAS_SOURCE)/exodus:$(TRUCHAS_SOURCE)/grid_mapping:$(TRUCHAS_SOURCE)/utilities

INC= -I$(NETCDF_MOD_DIR)
LIB= -L$(NETCDF_LIB_DIR) -lnetcdf

# Linux with the Lahey, NAG, Pathscale compiler
ifeq (Linux, $(shell uname))
    F90      = lf95
    F90FLAGS = -O $(INC)
  ifeq ($(compiler), nag)
    F90      = f95
    F90FLAGS = -O3 -PIC $(INC) -DNAG
  endif  
  ifeq ($(compiler), pathscale)
    F90      = pathf90
    F90FLAGS = -O -fPIC $(INC) 
  endif
  ifeq ($(compiler), g95)
    F90      = g95
    F90FLAGS =  -O2 -fno-second-underscore -ffree-line-length-huge -fimplicit-none -static -fPIC $(INC) -DG95
  endif  
endif
                                                                                
PROGS   = wrapper

EXOBJ   = f90_assert.o string_utilities.o \
	exodus_errors.o exodus_mesh_type.o exodus_mesh_utilities.o \
	exodus_mesh_reader.o exodus_mesh_writer.o exodus_truchas_hack.o exodus.o

OBJ     = hpsort.o overlap_module.o gm_mesh_type.o grid_mapping_utils.o \
	grid_mapping_module.o map_wrapper.o

MAINOBJ = main.o

ifeq ($(mode), main)
  OBJ = $(EXOBJ) hpsort.o overlap_module.o gm_mesh_type.o grid_mapping_utils.o \
	grid_mapping_module.o map_wrapper.o $(MAINOBJ)
endif

wrapper : $(OBJ)
	$(F90) $(F90FLAGS) -o $@ $(OBJ) $(LIB)

CLEANFILES= *.o *.mod int_vols_file
clean:; rm -f $(CLEANFILES)
cleanest: clean; rm -f $(PROGS)

.SUFFIXES: .F90
.F90.o:
	$(F90) $(F90FLAGS) -c $<
