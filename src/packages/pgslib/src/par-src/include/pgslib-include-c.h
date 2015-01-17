/* General include file for all c routines in PGSLib */

/* $Id: pgslib-include-c.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

#ifndef PGSLIB_INCLUDE_H__
#define PGSLIB_INCLUDE_H__

#define USE_MPI

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "pgslib-types.h"

/* Global defines used for all PGSLib C routines */
/* These need to be consistent with the definitions in pgslib_types_module.F */
#define TRUE  1
#define FALSE 0
#define PGSLIB_TRUE  1
#define PGSLIB_FALSE 0

/* Need to know how F90 stores characters. This will be different on each system. */
#define BYTES_PER_CHAR 1

/* Need to know how F90 stores logicals.  These will be different on each system. */
#define BYTES_PER_LOG  4
#define C_LOG_TYPE     char
#define MPI_LOG_TYPE   MPI_BYTE

extern parallel_state pgslib_state; /* Declared in utility-c.c */

#ifdef USE_SGI_SHMEM_LIB
extern int * sm_space;
#endif

#endif
