TRUCHAS_ROOT = ../../..

-include $(TRUCHAS_ROOT)/configuration

TRUCHAS_INCLUDES = -I$(BUILD_DIR)
TRUCHAS_LIBRARIES = -L$(BUILD_DIR) -ltruchas

LIB = $(TRUCHAS_LIBRARIES) $(PGSLIB_LIBRARIES) $(HYPRE_LIBRARIES) $(MPI_LIBRARIES)

FFLAGS += $(TRUCHAS_INCLUDES) $(PGSLIB_INCLUDES)

default:
	$(FC) $(FFLAGS) -o a.out test_scalar_func_copy.F90 $(LIB)
	$(FC) $(FFLAGS) -o a.out test_pressure_poisson.F90 $(LIB)
