/* C routines to support PREFIX and SUFFIX operations */
/* These routines perform a global scan,
   the routines assume only a single item per processor.  No 
   on-pe scans are done at this level. */

/* $Id: scan-c-double.c,v 1.1.1.1 2000/10/11 22:44:23 ferrell Exp $ */

#include "mpi.h"
#include "pgslib-include-c.h"

#define PGSLIB_DATA_TYPE     double
#define PGSLIB_ROUTINE_NAME off_node_sum_prefix_ ## double_c
#define PGSLIB_MPI_DATA_TYPE MPI_DOUBLE_INT

#define _PREFIX_    1
#undef  _SUFFIX_
#define  _SRC_SIGN_ -
#define _DEST_SIGN_ +

#include "scan-c.h"

/**********************************************************************/

#define PGSLIB_DATA_TYPE     double
#define PGSLIB_ROUTINE_NAME off_node_sum_suffix_ ## double_c
#define PGSLIB_MPI_DATA_TYPE MPI_DOUBLE_INT

#define _SUFFIX_    1
#undef  _PREFIX_
#define  _SRC_SIGN_ +
#define _DEST_SIGN_ -

#include "scan-c.h"


