# Set this to the path of the Truchas root directory
TRUCHAS_ROOT = ../../../..

include ${TRUCHAS_ROOT}/configuration

LIB = -L${BUILD_DIR} -ltruchas ${PGSLIB_LIBRARIES} ${CHACO_LIBRARIES} $(NETCDF_LIBRARIES) ${MPI_LD_FLAGS} -lmpi -ldl

INC = $(NETCDF_INCLUDES) $(PGSLIB_INCLUDES) -I$(BUILD_DIR)

FFLAGS += -u -gline -C -DNAG $(INC)

example: example.o
	$(FL) $(FFLAGS) -o $@ example.o  $(LIB)

CLEANFILES= *.mod *.MOD *.o

clean:; rm -f $(CLEANFILES)
cleanest: clean; rm -f $(PROGS) libre.a

.SUFFIXES: .F90
%.o: %.F90
	$(FC) $(FFLAGS) -c $<
