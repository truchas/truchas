/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * memory.c
 *
 *  DANU Memory Manager
 *
 *
 *  Purpose:
 *           This source file defines wrappers to malloc, calloc
 *           and realloc. The macros defined in danu_memory.h
 *           should be used instead of calling these routines
 *           directly. The macro simplify malloc statements by
 *           accepting the pointer type and size relative to the type
 *           and returning the correct pointer.
 *           Example:
 *
 *           ptr = DANU_MALLOC(int,128)
 *
 *           is equivalent to
 *
 *           ptr = (int *) malloc(sizeof(int)*128);
 *
 */

#if HAVE_CONFIG_H
# include <danu_config.h>
#endif

#include <stdlib.h>
#include <string.h>

#include <danu_error.h>
#include <danu_memory.h>

/*
 * Routine: void* ptr = danu_malloc(size_t num)
 * Purpose: Allocate specified bytes of memory 
 * Description: Pass in the number of bytes and return a void pointer
 *              to the new memory location. Routine will print an error
 *              message if malloc returns a NULL pointer. Calling functions
 *              or code should use the macros defined in danu_memory.h instead
 *              of calling this routine directly.
 *
 * Parameters:
 *           num		IN		Number of bytes
 *                              
 * Returns: Returns a void type pointer.
 * Errors: The returning pointer will be set to NULL if the malloc fails.
 *
 */ 
void * danu_malloc(size_t num)
{
    void * ptr = NULL;

    if ( num <= 0 ) {
      DANU_ERROR_MESS("Invalid memory size requested");
    }
    else {
      ptr = malloc(num);
      if ( ptr == NULL ) {
	danu_debug_printf("Attempted to allocated num=%lu bytes\n",
			  (unsigned long) num);
	DANU_ERROR_MESS("Malloc failed ... memory exhausted?");
      }
    }

    return ptr;
}
/*
 * Routine: void* ptr = danu_calloc(size_t num)
 * Purpose: Allocate specified bytes of memory and initialize the memory 
 *          with \0 
 * Description: Pass in the number of bytes and return a void pointer
 *              to the new memory location. The memory is initialized with
 *              \0 values before returning. Routine will print an error
 *              message if malloc returns a NULL pointer. Calling functions
 *              or code should use the macros defined in danu_memory.h instead
 *              of calling this routine directly.
 *
 * Parameters:
 *           num		IN		Number of bytes
 *                              
 * Returns: Returns a void type pointer.
 * Errors: The returning pointer will be set to NULL if the malloc fails.
 *
 */ 
void * danu_calloc(size_t num)
{
    void * ptr = danu_malloc(num);

    memset(ptr,'\0', num);

    return ptr;
}
/*
 * Routine: void* new = danu_realloc(void * ptr, size_t num)
 * Purpose: Attempts to resize the data allocation pointed to by pointer ptr 
 * Description: Pass in the new number of bytes and a pointer, if the pointer is NULL
 *              then a standard malloc is called otherwise realloc is called. 
 *              If the space can be expended to the new size num bytes
 *              then the memory is expanded and realloc returns ptr. If there is 
 *              not sufficient space to increase the memory then the behavior will
 *              depend on the system implementation of realloc. Please consult the system
 *              documentation (man realloc). 
 *
 * Parameters:
 *           ptr		IN/OUT	void pointer pointing to a memory address
 *           num		IN		Resize memory space ptr points to to num bytes 
 *                              
 * Returns: Returns a void type pointer.
 * Errors: The returning pointer will be set to NULL if the malloc fails or if the
 *         realloc cannot increase the memory location size. Warning messages will
 *         be triggered if ptr is changed to NULL indicating that the current location
 *         could not be resized.
 *
 */ 
void * danu_realloc(void * ptr, size_t num)
{
    void * new;
    
    if ( ptr == NULL ) {

        new = danu_malloc(num);

    }
    else {

        new = realloc(ptr,num);

        if ( new == NULL ) {
            danu_error("Realloc failed ... memory exhausted?");
        }

		if ( new != ptr && new != NULL ) {
			danu_warn("Realloc required new pointer");
		}

		if ( ptr == NULL ) {
			danu_warn("Pointer passed into danu_realloc no longer valid");
		}

    }

    return new;
}
/*
 * Routine: danu_free(void * ptr)
 * Purpose: Release the memory
 * Description: Release the memory ptr is assigned to if ptr is
 *              not NULL.
 *              
 * Parameters:
 *           ptr		IN/OUT		Pointer
 *                              
 * Returns: Nothing
 * Errors: This routine never fails.
 *
 */ 
void danu_free(void * ptr)
{
    if ( ptr != NULL ) {
        free(ptr);
    }
    else {
        danu_warn("Attempting to free NULL pointer");
    }

}
