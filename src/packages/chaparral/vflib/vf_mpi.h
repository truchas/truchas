#ifndef _VF_MPI_H_
#define _VF_MPI_H_

#ifdef VF_NO_MPI
typedef int MPI_Comm;
#define MPI_COMM_WORLD 0
#else
#include "mpi.h"
#endif

#endif
