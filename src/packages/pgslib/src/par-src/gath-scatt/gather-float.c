/* This routine moves data from nodeCommBufferData to offPEnodesData.
   gsENdata must be setup before this call.  That means that indexing
   and init_comm have been done.
   
The input is:
   nnodesThisPE
   nodeCommBufferData
   gsENdata (which knows how to do the communication)
   offPEnodesData (already allocated, but not loaded)

The output is node data at the elements
   offPEnodesData

   
   
*/

/* $Id: gather-float.c,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "pgslib-include-c.h"

#define PGSLIB_DATA_TYPE float
#define PGSLIB_ROUTINE_TYPE_POSTFIX float
#define PGSLIB_ROUTINE_NAME(Base_Name) Base_Name ## float_c
#define PGSLIB_TYPE_NAME(Base_Name) Base_Name ## float
#define PGSLIB_MPI_DATA_TYPE MPI_FLOAT

#include "pgslib-types.h"

#include "gather.h"
