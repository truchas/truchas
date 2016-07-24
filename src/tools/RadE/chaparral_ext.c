#include "vf_api.h"

/* An MPI communicator is not easily available on the Fortran
   side so we resort to using this C wrapper for VF_Setup. */

void VF_SetupDefault()
{
  VF_Setup(0,1,MPI_COMM_WORLD);
}
