/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * danu_memory.h
 *
 *  DANU Memory Manager
 *
 */

#ifndef DANU_MEMORY_H
#define DANU_MEMORY_H


#include <stdlib.h>

/* MALLOC and FREE wrappers */
void * danu_malloc(size_t num);
void * danu_calloc(size_t num);
void * danu_realloc(void * ptr, size_t num);
void   danu_free(void * ptr);

/* MALLOC macros
 *
 *  Usage:
 *
 *  DANU_MALLOC(real,8)
 *  is equivalent to
 *
 *  (real *) malloc(sizeof(real)*8);
 */ 

#define DANU_MALLOC(type,num)           (type *) danu_malloc(num*sizeof(type)) 
#define DANU_CALLOC(type,num)           (type *) danu_calloc(num*sizeof(type)) 
#define DANU_REALLOC(type,ptr,num)      (type *) danu_realloc(ptr,num*sizeof(type)) 
#define DANU_FREE(old)                  danu_free((old))


#endif
