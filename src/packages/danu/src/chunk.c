/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * chunk_params.c
 *
 *  DANU Chunking Parameters
 *
 *
 *  Purpose:
 *
 *          The HDF5 library allows users to group raw data into chunks
 *          which can improve I/O performance. Chunking data is also
 *          required if a user would like to create an extend-able dataset.
 *          This module creates a data structure that holds all of the 
 *          parameters that control chunk cache behavior and chunk dimensions
 *          and sizes.
 *
 */

#include <hdf5.h>

#include <danu_error.h>
#include <danu_types.h>
#include <danu_memory.h>

#include <danu_chunk.h>


/*
 * Routine: dchunk_params_t * danu_chunk_params_alloc(int ndim)
 * Purpose: Allocates a chunk parameter data struct 
 * Description: Pass in the dimension of the dataset that will use the
 *              chunk parameters. Routine allocates a chunk parameter
 *              struct and a pointer that will point to the dimensions
 *              of the chunk data. 
 *
 * Parameters:
 *           ndim                IN              Number of dimensions of the dataset to chunk
 *                              
 * Returns: Returns a dchunk_param_t pointer
 * Errors: The returning pointer will be set to NULL if an error occurs.
 *
 */ 
dchunk_params_t * danu_chunk_params_alloc(int ndim)
{
	dchunk_params_t * ptr;

	if ( ndim <= 0 ) {
		DANU_ERROR_MESS("Invalid dimension size (<=0)");
		ptr = NULL;
	}
	else {
		ptr = DANU_MALLOC(dchunk_params_t,1);

		ptr->cache_nslots = H5D_CHUNK_CACHE_NSLOTS_DEFAULT; 
		ptr->cache_nbytes = H5D_CHUNK_CACHE_NBYTES_DEFAULT; 
		ptr->cache_w0     = H5D_CHUNK_CACHE_W0_DEFAULT; 

		ptr->size = DANU_MALLOC(hsize_t,ndim);

		if ( ptr->size == NULL ) {
			DANU_ERROR_MESS("Invalid pointer returned from malloc");
			DANU_FREE(ptr);
		}
	}

	return ptr;
}
/*
 * Routine: void danu_chunk_params_free(dchunk_params_t * ptr)
 * Purpose: Frees a chunk parameter data struct 
 * Description: Pass in a pointer to a chunk parameter data struct
 *              and release all the memory associated with that data
 *              struct. If ptr is NULL, then routine does nothing.
 *
 * Parameters:
 *           ptr                IN              Pointer to chunk parameters data struct
 *                              
 * Returns: Returns nothing
 * Errors: This routine should never fail since nothing is freed if the ptr == NULL
 *
 */ 
void danu_chunk_params_free(dchunk_params_t * ptr)
{

	if ( ptr != NULL ) {
		if ( ptr->size != NULL ) {
			DANU_FREE(ptr->size);
		}

		DANU_FREE(ptr);
	}

}







