VPATH=.:..:../../utilities

os=$(shell uname)

FC=nagfor
FFLAGS=-g90 -gline -C=all -DNAG $(INC)

NETCDF_MOD_DIR=../../../packages/netcdf/linux-x86_64-nag/include
NETCDF_LIB_DIR=../../../packages/netcdf/linux-x86_64-nag/lib

INC= -I$(NETCDF_MOD_DIR)
LIB= -L$(NETCDF_LIB_DIR) -lnetcdff -lnetcdf

PROGS= test1

OBJ1= 	string_utilities.o exodus_errors.o exodus_mesh_type.o \
	exodus_mesh_reader.o exodus_mesh_writer.o exodus_truchas_hack.o \
	exodus_mesh_utilities.o exodus.o test1.o

CLEANFILES= *.o *.mod *.MOD *.copy *.g90

test1: $(OBJ1)
	$(F90) -o $@ $(OBJ1) $(LIB)

clean:; rm -f $(CLEANFILES)
cleanest: clean; rm -f $(PROGS)

.SUFFIXES: .F90
.F90.o:
	$(FC) $(FFLAGS) -c $<
