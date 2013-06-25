/* *************************************************************************** *
*                                                                              *
*                                                                              *
*                             Copyright  (C) 20xx,                             *
*                      Los Alamos National Security, LLC                       *
*                                                                              *
*                             LA-CC-xxxxxx                                     *
*                                                                              *
* **************************************************************************** */

/*
* danu_types.h
*
*  DANU simple data types
*
*/

#ifndef DANU_TYPES_H
#define DANU_TYPES_H

#include <hdf5.h>

/* Data types I find useful */
typedef int           dbool_t;       /* Boolean 0=FALSE, 1=TRUE */
typedef uint64_t      didx_t;        /* Array index types */
typedef int           dcyc_t;        /* Cycle number */


typedef int           ddim_t;          /* Array, dataset dimensions */

typedef uint64_t      dsize_t;        /* Always want 64-bit addressing */

/* Useful MACROS */
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif



#endif
