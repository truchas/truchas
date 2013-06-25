/* C routines to support physical shift routines */

/* $Id: shift-c-log.c,v 1.1.1.1 2000/10/11 22:44:23 ferrell Exp $ */

#include "mpi.h"
#include "pgslib-include-c.h"

#define PGSLIB_DATA_TYPE     C_LOG_TYPE
#define SIZE_MULTIPLIER      BYTES_PER_LOG
#define PGSLIB_ROUTINE_NAME  physical_c_shift_up_ ## log_c
#define PGSLIB_MPI_DATA_TYPE MPI_LOG_TYPE

#define _DIRECTION_SIGN_ +

#include "shift-c.h"

/**********************************************************************/

#define PGSLIB_DATA_TYPE     C_LOG_TYPE
#define SIZE_MULTIPLIER      BYTES_PER_LOG
#define PGSLIB_ROUTINE_NAME  physical_c_shift_down_ ## log_c
#define PGSLIB_MPI_DATA_TYPE MPI_LOG_TYPE

#define _DIRECTION_SIGN_ -

#include "shift-c.h"


