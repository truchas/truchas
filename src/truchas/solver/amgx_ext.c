#include "amgx_c.h"
#include "mpi.h"

// device, double matrix, double vectors
#define MODE AMGX_mode_dDDI

int AMGX_resources_create_ext(AMGX_resources_handle *resources, AMGX_config_handle config,
                              int ndevices, int *gpu_ids)
{
   return AMGX_resources_create(resources, config, MPI_COMM_WORLD, ndevices, gpu_ids);
}

int AMGX_matrix_create_ext(AMGX_matrix_handle *A, AMGX_resources_handle resources)
{
  return AMGX_matrix_create(A, resources, MODE);
}

int AMGX_vector_create_ext(AMGX_vector_handle *x, AMGX_resources_handle resources)
{
  return AMGX_vector_create(x, resources, MODE);
}

int AMGX_solver_create_ext(AMGX_solver_handle *solver, AMGX_resources_handle resources,
                           AMGX_config_handle config)
{
  return AMGX_solver_create(solver, resources, MODE, config);
}
