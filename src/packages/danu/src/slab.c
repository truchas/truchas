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
 * slab.c
 *
 *  DANU Hyperslabs
 *
 *
 *  Purpose:
 *
 *
 */

#if HAVE_CONFIG_H
# include <danu_config.h>
#endif

#include <hdf5.h>

#include <danu_error.h>
#include <danu_types.h>
#include <danu_memory.h>

#include <danu_slab.h>

/*
 * Routine: dslab_t danu_slab_alloc(hid_t dataspace)
 * Purpose: Create a HDF5 Hyperslab associated with HDF5 dataspace
 * Description: Given a defined HDF5 dataspace, allocate a hyperslab
 *              struct that contains pointers for offset, count, stride
 *              and block size. The dimensions of these array will be 1xdim
 *              where dim is the dimension of the dataspace.
 *              
 *
 * Parameters:
 *           dataspace                IN   HDF5 dataspace           
 *                              
 * Returns: Returns a pointer to a slab structure
 * Errors: Will return immediately if dataspace is not a valid HDF5 id. If 
 *         an error occurs the returning pointer is NULL.
 */
 dslab_t * danu_slab_alloc(hid_t id)
 {
     dslab_t *ptr = NULL;
     H5I_type_t type;
     hid_t dataspace;
     int dim, close_space;

     /* Check the id */
     if ( H5_ISA_INVALID_ID(id) ) {
         DANU_ERROR_MESS("Invalid HDF5 identifier");
         return ptr;
     }

     /* Id maybe a dataspace or a dataset */
     type = H5Iget_type(id);

     if ( type == H5I_DATASET ) {
         close_space = TRUE;
         dataspace = H5Dget_space(id);
     }
     else {
         close_space = FALSE;
         dataspace = id;
     }

     /* Get the dimensions */ 
     dim = H5Sget_simple_extent_ndims(dataspace);

     if ( dim < 0 ) {
         DANU_ERROR_MESS("Failed to read the dimensions from HDf5 dataspace");
         return ptr;
     }

     ptr = DANU_MALLOC(dslab_t,1);

     ptr->dim = dim;
     ptr->offset = DANU_MALLOC(hsize_t,dim);
     ptr->block = DANU_MALLOC(hsize_t,dim);
     ptr->stride = DANU_MALLOC(hsize_t,dim);
     ptr->count = DANU_MALLOC(hsize_t,dim);

     /* Close the dataspace if opened here */
     if ( close_space ) {
         H5Sclose(dataspace);
     }

     return ptr;
 }
/*
 * Routine: danu_slab_free(dslab_t *ptr)
 * Purpose: Free dslab_t struct memory 
 * Description: Free all the memory associated with hyper slab struct
 *
 * Parameters:
 *           ptr                IN   pointer to hyperslab struct          
 *                              
 * Returns: Nothing
 * Errors:  No errors are raised. The DANU_FREE macro will not free NULL pointers.
 */

 void danu_slab_free(dslab_t * ptr)
 {
     DANU_FREE(ptr->offset);
     DANU_FREE(ptr->block);
     DANU_FREE(ptr->stride);
     DANU_FREE(ptr->count);

     DANU_FREE(ptr);
 }
/*
 * Routine: herr_t danu_slab_select(hid_t dataspace,dslab_t *slab)
 * Purpose:  Select the hyperslab in dataspace
 * Description: Selects the hyperslab defined in the slab data struct slab. The call
 *              replaces any hyperslab currently defined in dataspace and will not
 *              allow overlapping blocks.
 *
 * Parameters:
 *           dataspace          IN   HDF5 dataspace
 *           slab               IN   pointer to hyperslab struct          
 *                              
 * Returns: The returning value of the HDF5 select hyperslab call
 * Errors: Checks dataspace is a valid HDF5 id and the dimensions of dataspace and slab
 *         will return immediately if these checks raise an error. 
 */
 herr_t danu_slab_select(hid_t id,dslab_t *slab)
 {
     int    dim,close_space;
     herr_t status;
     hid_t  dataspace;
     H5I_type_t type;

     /* Check input */
     if ( H5_ISA_INVALID_ID(id) ) {
         DANU_ERROR_MESS("Invalid dataspace identifier");
         return DANU_FAILURE;
     }

     if ( DANU_BAD_PTR(slab) ) {
         DANU_ERROR_MESS("Invalid slab pointer");
         return DANU_FAILURE;
     }

     /* Determine if id is a dataspace or dataset */
     type = H5Iget_type(id);
     if ( type == H5I_DATASET ) {
         dataspace = H5Dget_space(id);
         close_space = TRUE;
     }
     else {
         dataspace = id;
         close_space = FALSE;
     }

     dim = H5Sget_simple_extent_ndims(dataspace);
     if ( dim < 0 ) {
         DANU_ERROR_MESS("Failed to read dimensions from dataspace");
         return DANU_FAILURE;
     }

     if ( dim != slab->dim ) {
         DANU_ERROR_MESS("Mismatch between hyperslab dimension and dataspace dimension");
         return DANU_FAILURE;
     }

     status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,slab->offset,slab->stride,slab->count,slab->block);

     if ( close_space ) {
         H5Sclose(dataspace);
     }

     return status;


}
/*
 * Routine: herr_t danu_slab_contiguous(dslab_t slab, hsize_t *offset, hsize_t *size)
 * Purpose:  Define the a contiguous hyperslab
 * Description: Selects defines a contiguous hyperslab. Both offset and size  must be
 *              linear arrays 1 x slab->dim. This routine sets the stride and the 
 *              block arrays to 1. HDF5 allows users to pass NULL pointers for those 
 *              arrays and assumes that a NULL pointer is equivalent to an array of 1's,
 *              however the only place a slab struct is free'd is in danu_slab_free. So
 *              this routine populates the stride and block arrays with 1's to avoid 
 *              freeing memory here.
 *
 * Parameters:
 *           slab               IN   pointer to hyperslab struct          
 *           offset             IN   pointer to hyperslab offset
 *           size               IN   size of hyperslab
 *                              
 * Returns: The returning value is nonnegative if no errors are raised, negative value
 *          returned if an error occurs.
 * Errors: Checks dataspace is a valid HDF5 id and the dimensions of dataspace and slab
 *         will return immediately if these checks raise an error. 
 */


herr_t danu_slab_contiguous(dslab_t *ptr, hsize_t *offset, hsize_t *size)
{
    int dim;
    didx_t i;
    
    /* Check input */
    if ( DANU_BAD_PTR(ptr) ) {
        DANU_ERROR_MESS("Invalid pointer to hyperslab");
        return DANU_FAILURE;
    }

    if ( DANU_BAD_PTR(offset) ) {
        DANU_ERROR_MESS("Invalid pointer to offset");
        return DANU_FAILURE;
    }

    if ( DANU_BAD_PTR(size) ) {
        DANU_ERROR_MESS("Invalid pointer to size");
        return DANU_FAILURE;
    }

    /* HDF5 interpets NULL's for stride and count as arrays of 1's 
       Since I do not want to free the memory anywhere else but in danu_slab_free
       I am forced to put 1's in these arrays.Sigh.
    */
    dim = ptr->dim;
    for(i=0;i<dim;i++) {
        ptr->block[i]  = 1;
        ptr->stride[i] = 1;
        ptr->offset[i] = offset[i];
        ptr->count[i]  = size[i];
    }

    return DANU_SUCCESS;
}



    


        
    
    


