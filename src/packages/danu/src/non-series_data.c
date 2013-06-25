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
 * non_series_data.c
 *
 *  DANU Non-series Datasets
 *
 *
 *  Purpose:
 *           This source file defines functions that create non_series datasets. 
 *
 *
 */

#if HAVE_CONFIG_H
#include <danu_config.h>
#endif

#include <stdarg.h>
#include <string.h>

#include <hdf5.h>

#include <danu_error.h>
#include <danu_memory.h>
#include <danu_utils.h>
#include <danu_group.h>
#include <danu_dataset.h>

#include <danu_sim.h>
#include <danu_non-series.h>

/*
 * Routine: 
 * Purpose:
 * Description:             
 *
 *
 * Parameters:
 *           mid                IN              HDF5 identifier for existing mesh group 
 *           nnodes             IN              Number of nodes
 *           x                  IN              X (dim 0) nodal coordinate data
 *           y                  IN              Y (dim 1) nodal coordinate data [OPTIONAL]
 *           z                  IN              Z (dim 2) nodal coordinate data [OPTIONAL]
 *                              
 * Returns:
 *          
 * Errors:s
 *        
 */

/*
 * Routine: data_rank(hid_t sid, const char *name, int *rank) 
 * Purpose: Return the rank of a non series dataset.
 * Description: Find the rank of dataset name found under simulation group
 *              sid. Routine checks input and returns a newgative value if
 *              an error is detected. Return code should be checked before
 *              rank value is used in the calling routine.
 *
 *
 * Parameters:
 *           sid        IN      HDF5 identifier for simulation group 
 *           name       IN      Dataset name
 *           rank       OUT     Dataset rank
 *                              
 * Returns: Returns a negative value if an error is detected.
 *          
 * Errors: Errors will be raised if sid is an invalid id or dataset does not
 *         exist.
 */
herr_t data_rank(hid_t sid, const char * name, int *rank)
{
  herr_t status = DANU_FAILURE;
  hid_t gid;
  int exists;

  if ( DANU_RETURN_FAIL(data_exists(sid,name,&exists)) ) {
    DANU_ERROR_MESS("Can not stat dataset");
    return status;
  }

  if ( exists ) {
    gid = data_open_group(sid);
    if ( H5_ISA_VALID_ID(gid) ) {
      *rank = danu_dataset_rank(gid,name);
      if ( *rank > 0 ) {
        status = DANU_SUCCESS;
      }
      else {
        DANU_ERROR_MESS("Failed to find dataset rank");
      }
      danu_group_close(gid);
    }
    else {
      DANU_ERROR_MESS("Failed to open non-series group");
    }
  }
  else {
    DANU_ERROR_MESS("Dataset does not exist");
  }

  return status;

}
/*
 * Routine: data_dimensions(hid_t sid,
                            const char *name, int rank, int *dimensions)
 * Purpose: Return the dimensions of a non-series dataset.
 * Description: Find the dimensions of dataset name under the non-series
 *              group. Routine checks the input and will return a negative 
 *              number if an error is raised. Pointer to dimensions is
 *              a 1D array. This array is initialized to zero. Calling
 *              routine should check the return code before using the
 *              values in dimension.
 *
 * Parameters:
 *           sid          IN      HDF5 identifier for simulation group 
 *           name         IN      Dataset name
 *           rank         IN      Size of dimension array
 *          *dimensions   OUT     Dimensions of the dataset
 *                              
 * Returns: Returns a negative value if an error is detected.
 *          
 * Errors: Errors will be raised if sid is an invalid id or dataset does not
 *         exist.
 */
herr_t data_dimensions(hid_t sid,
		       const char * name,
		       int rank,
		       int *dimensions)
{
  herr_t status = DANU_FAILURE;
  hid_t gid;
  int i,exists;
  hsize_t *dims;

  /* Check input */
  if ( DANU_RETURN_FAIL(data_exists(sid,name,&exists)) ) {
    DANU_ERROR_MESS("Can not stat dataset");
    return status;
  }

  if ( rank <=0 )  {
    DANU_ERROR_MESS("Invalid dimension array size");
    return status;
  }


  /* Initialize the array */
  dims=DANU_CALLOC(hsize_t,rank);

  if ( exists ) {
    gid = data_open_group(sid);
    if ( H5_ISA_VALID_ID(gid) ) {
      status = danu_dataset_dimensions(gid,name,rank,dims);
    }
    else {
      DANU_ERROR_MESS("Failed to open non-series group");
    }
  }
  else {
    DANU_ERROR_MESS("Dataset does not exist");
  }

  /* Release the memory */
  for(i=0;i<rank;i++)
    dimensions[i] = (int) dims[i];
  DANU_FREE(dims);

  return status;

}
/*
 * Routine: data_type(hid_t sid, const char *name, int *type) 
 * Purpose: Return the typecode of non-series dataset name.
 * Description: Define the typecode of dataset name found under simulation group
 *              sid. Typecodes are defined in daun_dataset_types.h.
 *              Routine checks input and returns a negative value if
 *              an error is detected. Return code should be checked before
 *              type value is used in calling routines.
 *
 *
 * Parameters:
 *           sid        IN      HDF5 identifier for simulation group 
 *           name       IN      Dataset name
 *           type       OUT     Dataset typecode
 *                              
 * Returns: Returns a negative value if an error is detected.
 *          
 * Errors: Errors will be raised if sid is an invalid id or dataset does not
 *         exist.
 */
herr_t data_type(hid_t sid, const char * name, int *type)
{
  herr_t status = DANU_FAILURE;
  hid_t gid;
  int exists;

  if ( DANU_RETURN_FAIL(data_exists(sid,name,&exists)) ) {
    DANU_ERROR_MESS("Can not stat dataset");
    return status;
  }

  if ( exists ) {
    gid = data_open_group(sid);
    if ( H5_ISA_VALID_ID(gid) ) {
      status = danu_dataset_type(gid,name,type);
      danu_group_close(gid);
    }
    else {
      DANU_ERROR_MESS("Failed to open non-series group");
    }
  }
  else {
    DANU_ERROR_MESS("Dataset does not exist");
  }

  return status;

}
/*
 * Routine: data_create_group(hid_t sid) 
 * Purpose: Create the non_series sub group under simulation sid
 * Description: Create the non_series subgroup under simuation sid
 *              The non_series group will be closed before returning.
 *              Errors are raised if sid is not a valid HDf5 id or
 *              fail to create the subgroup.            
 *
 *
 * Parameters:
 *           sid                IN              HDF5 identifier for simulation group 
 *                              
 * Returns: Returns a negative value if an error is detected
 *          
 * Errors: Errors will be raised if sid is an invalid id.
 *         Failure to create the subgroup is also an error. 
 *        
 */
herr_t data_create_group(hid_t sid)
{
    herr_t status = DANU_FAILURE;

    char non_series_group_name[] = NON_SERIES_GROUP_NAME;
    hid_t gid = danu_group_create(sid,non_series_group_name);

    if ( H5_ISA_VALID_ID(gid) ) {
        status = 0;
        danu_group_close(gid);
    }
    else {
        DANU_ERROR_MESS("Failed to create the non_series group");
    }

    return status;
}

/*
 * Routine: data_open_group(hid_t sid) 
 * Purpose: Open the non_series sub group under simulation sid
 * Description: Open the non_series subgroup under simuation sid
 *
 *
 * Parameters:
 *           sid                IN              HDF5 identifier for simulation group 
 *                              
 * Returns: Returns a negative value if an error is detected
 *          
 * Errors: Errors will be raised if sid is an invalid id.
 *        
 */
hid_t data_open_group(hid_t sid)
{
    char non_series_group_name[] = NON_SERIES_GROUP_NAME;
    hid_t gid = danu_group_open(sid,non_series_group_name);

    if ( H5_ISA_INVALID_ID(gid) ) {
        printf("gid=%lX sid=%lx\n", (long) gid, (long) sid);
        DANU_ERROR_MESS("Failed to open the non_series group");
    }

    return gid;
}
/*
 * Routine: data_open_dataset(hid_t sid, const char * dataname) 
 * Purpose: Open dataset dataname under the the non_series sub group under simulation sid
 * Description: Open the non_series dataset dataname under simuation sid and
 *              return the HDF5 identifier for that dataset
 *
 *
 * Parameters:
 *           sid                IN              HDF5 identifier for simulation group 
 *           dataname           IN              Dataset name
 *                              
 * Returns: Returns a negative value if an error is detected
 *          
 * Errors: Errors will be raised if sid is an invalid id or the routine fails to open
 *         the dataset.
 *        
 */
hid_t data_open_dataset(hid_t sid,const char *dataname)
{
    hid_t gid;
    hid_t id = H5I_BADID;
    int   exists;

    if ( H5_RETURN_OK(data_exists(sid,dataname,&exists)) ) {
        if ( exists ) {
            gid = data_open_group(sid);
            id  = danu_dataset_open(gid,dataname);
            danu_group_close(gid);
        }
        else {
            DANU_ERROR_MESS("Non-series dataset does not exist");
        }
    }
    else {
        DANU_ERROR_MESS("Failed to stat the non_series dataset");
    }

    return id;
}






/*
 * Routine: data_exists(hid_t sid, const char * dataname, int *exists) 
 * Purpose: Checks the existence of dataset dataname under the Non-series subgroup
 * Description: Checks the existence of dataname under a Non-series subgroup in a
 *              Simulation group identified by sid. 'exists' is intially set to 
 *              FALSE. If the dataset is found exists is set to TRUE. 
 *
 * Parameters:
 *           sid           IN              HDF5 identifier for the Simulation group
 *          *dataname      IN              Dataset name to search
 *          *int           OUT             Flag set to TRUE or FALSE 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *          
 * Errors: Error is raised if the input is not valid, routine fails to open the 
 *         non_series data group or the existence check returns an error. 
 *         The danu_group_open routine checks the sid identifier. Data name
 *         input is checked in the danu_dataset_exists call. 
 *        
 */
 herr_t data_exists(hid_t sid, const char *dataname, int *exists)
 {
     herr_t status = DANU_FAILURE;

     hid_t gid = data_open_group(sid);

     *exists = FALSE;

     if ( H5_ISA_INVALID_ID(gid) ) {
         DANU_ERROR_MESS("Failed to open non-series data group");
     }
     else {
         *exists = danu_dataset_exists(gid,dataname);
         status = 0;
         danu_group_close(gid);
     }

     return status;

}
/*
 * Routine: data_count(hid_t sid,int *ndata) 
 * Purpose: Determine the number of non_series datasets found in simulation group sid
 * Description: Counts the number of datasets found in the non_series subgroup under
 *              the simulation group identified with sid. Routine assumes that only
 *              datasets are found under the non_series group. This save memory and
 *              time. ndata is initally set to zero.
 *
 * Parameters:
 *           sid      IN              HDF5 identifier for the simulation group
 *          *ndata    OUT             Number of datasets found the in non_series subgroup
 *                              
 * Returns: A negative value if an error occurs.
 *          
 * Errors: Input is checked in routine to open the non_series group. Error occurs if input
 *         is not valid or error occurs while querying the non_series data group. 
 *        
 */
herr_t data_count(hid_t sid, int *ndata)
{
    herr_t status = DANU_FAILURE;
    hsize_t tmp  = 0;

    hid_t gid = data_open_group(sid);

    *ndata = 0;
    if ( H5_ISA_INVALID_ID(gid) ){
        DANU_ERROR_MESS("Failed to open non_series data group");
        return status;
    }

    status = danu_group_get_nlinks(gid,&tmp);
    *ndata = (int)tmp;

    danu_group_close(gid);

    return status;
}
/*
 * Routine: data_list(hid_t sid, int num, char **datanames, int * num_found) 
 * Purpose: Find the dataset names in the non_series data group
 * Description: Returns the names of non_series datasets available under the
 *              non_series subgroup found in the simulation group sid. Warning
 *              message appears if the number found does not match num. Callers
 *              should use the data_count routine to determine the correct number
 *              of datasets available. The datanames strings are initialized to
 *              \0 chars. This will be overwritten when a dataset name is inserted. 
 *              This routine will allocate memory for each pointer in datanames.
 *              The calling routine will be responsible for freeing this memory.
 *
 * Parameters:
 *           sid      IN              HDF5 identifier for the simulation group
 *           num      IN              Number of string pointers in datanames
 *         **datanames  OUT           Array of string pointers pointing to the dataset names
 *          *num_found                Pointer to int containing number of found datasets
 *                              
 * Returns: A negative value if an error occurs.
 *          
 * Errors: Input is checked in routine to open the non_series group. Error occurs if input
 *         is not valid or error occurs while querying the non_series data group.  
 *        
 */
herr_t data_list(hid_t sid, int num, char **datanames, int *num_found)
{
    herr_t status = DANU_FAILURE;

    hid_t gid = data_open_group(sid);

    if ( H5_ISA_VALID_ID(gid) ) {
        status = danu_group_get_datasets(gid,num,datanames,num_found);
        danu_group_close(gid);
    }
    else {
        DANU_ERROR_MESS("Failed to open non_series data subgroup");
    }


    return status;
}
/*
 * Routine: data_write_int(hid_t sid, const char *name, int dim, int *dimensions, int *data) 
 * Purpose: Write non_series integer data to the non_series subgroup if simulation sid
 * Description: Creates a dataset of dimensions simensions under non_series data subgroup
 *              with label name.
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for the simulation group
 *          *name        IN      Dataset name 
 *           dim         IN      Dimension of the dataset
 *          *dimensions  IN      1D array of length dim that stores the size of each dimension
 *          *data        IN      Pointer pointing to the data to be written 
 *                              
 * Returns: A negative value if an error occurs
 *          
 * Errors: Error occurs if the input is not valid, if dataset already exists, or an error is raised
 *         when writing the data.
 *        
 */
herr_t data_write_int(hid_t sid, const char *name, int dim, const int *dimensions, int *data)
{
    herr_t status = DANU_FAILURE;

    hsize_t * sizes;

    hid_t  gid;
    int    exists;

    if ( H5_RETURN_OK(data_exists(sid,name,&exists)) ) {
      
        if ( !exists ) {
            gid = data_open_group(sid);
            sizes = DANU_MALLOC(hsize_t, dim);
            convert_int_to_hsize(dim,dimensions,sizes);
            status = danu_data_write_int(gid,name,dim,sizes,data);
            DANU_FREE(sizes);
        }
        else {
            DANU_ERROR_MESS("Non series dataset already exists will not overwrite");
        }
    }
    else {
        DANU_ERROR_MESS("Failed to stat the non-series dataset");
    }

    return status;

}
/*
 * Routine: data_read_int(hid_t sid, const char *name, int dim, int *dimensions, int *data) 
 * Purpose: Read non_series integer data found in the non_series subgroup of simulation sid
 * Description: Read a dataset of dimensions dimensions under non_series data subgroup
 *              with label name.
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for the simulation group
 *          *name        IN      Dataset name 
 *           dim         IN      Dimension of the dataset
 *          *dimensions  IN      1D array of length dim that stores the size of each dimension
 *          *data        IN      Pointer pointing to the buffer to be read in to 
 *                              
 * Returns: A negative value if an error occurs
 *          
 * Errors: Error occurs if the input is not valid, if dataset does not exists, or an error is raised
 *         when writing the data.
 *        
 */
herr_t data_read_int(hid_t sid, const char *name, int dim, const int *dimensions, int *data)
{
    herr_t status = DANU_FAILURE;

    hid_t    gid;

    hsize_t * sizes;
    int      exists;

    if ( H5_RETURN_OK(data_exists(sid,name,&exists) ) ) {

        if ( exists ) {
            gid = data_open_group(sid);
            sizes = DANU_MALLOC(hsize_t,dim);
            convert_int_to_hsize(dim,dimensions,sizes);
            status =  danu_data_read_int(gid,name,dim,sizes,data);
            danu_group_close(gid);
            DANU_FREE(sizes);
        }
        else {
            DANU_ERROR_MESS("Non-series dataset does not exist");
        }

    }

    return status;

}
/*
 * Routine: data_write_double(hid_t sid, const char *name, int dim, int *dimensions, double *data) 
 * Purpose: Write non_series double data to the non_series subgroup if simulation sid
 * Description: Creates a dataset of dimensions simensions under non_series data subgroup
 *              with label name.
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for the simulation group
 *          *name        IN      Dataset name 
 *           dim         IN      Dimension of the dataset
 *          *dimensions  IN      1D array of length dim that stores the size of each dimension
 *          *data        IN      Pointer pointing to the data to be written 
 *                              
 * Returns: A negative value if an error occurs
 *          
 * Errors: Error occurs if the input is not valid, if dataset already exists, or an error is raised
 *         when writing the data.
 *        
 */
herr_t data_write_double(hid_t sid, const char *name, int dim, const int *dimensions, double *data)
{
    herr_t status = DANU_FAILURE;

    hsize_t * sizes;

    hid_t  gid;
    int    exists;

    if ( H5_RETURN_OK(data_exists(sid,name,&exists)) ) {
      
        if ( !exists ) {
            gid = data_open_group(sid);
            sizes = DANU_MALLOC(hsize_t,dim);
            convert_int_to_hsize(dim,dimensions,sizes);
            status = danu_data_write_double(gid,name,dim,sizes,data);
            DANU_FREE(sizes);
        }
        else {
            DANU_ERROR_MESS("Non series dataset already exists will not overwrite");
        }
    }
    else {
        DANU_ERROR_MESS("Failed to stat the non-series dataset");
    }

    return status;

}

/*
 * Routine: data_read_double(hid_t sid, const char *name, int dim, int *dimensions, double *data) 
 * Purpose: Read non_series integer data found in the non_series subgroup of simulation sid
 * Description: Read a dataset of dimensions dimensions under non_series data subgroup
 *              with label name.
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for the simulation group
 *          *name        IN      Dataset name 
 *           dim         IN      Dimension of the dataset
 *          *dimensions  IN      1D array of length dim that stores the size of each dimension
 *          *data        IN      Pointer pointing to the buffer to be read in to 
 *                              
 * Returns: A negative value if an error occurs
 *          
 * Errors: Error occurs if the input is not valid, if dataset does not exists, or an error is raised
 *         when writing the data.
 *        
 */
herr_t data_read_double(hid_t sid, const char *name, int dim, const int *dimensions, double *data)
{
    herr_t status = DANU_FAILURE;

    hid_t    gid;

    hsize_t * sizes;
    int      exists;

    if ( H5_RETURN_OK(data_exists(sid,name,&exists) ) ) {

        if ( exists ) {
            gid = data_open_group(sid);
            sizes = DANU_MALLOC(hsize_t,dim);
            convert_int_to_hsize(dim,dimensions,sizes);
            status =  danu_data_read_double(gid,name,dim,sizes,data);
            danu_group_close(gid);
            DANU_FREE(sizes);
        }
        else {
            DANU_ERROR_MESS("Non-series dataset does not exist");
        }

    }

    return status;

}
/*
 * Routine: data_write_float(hid_t sid, const char *name, int dim, int *dimensions, float *data) 
 * Purpose: Write non_series float data to the non_series subgroup if simulation sid
 * Description: Creates a dataset of dimensions simensions under non_series data subgroup
 *              with label name.
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for the simulation group
 *          *name        IN      Dataset name 
 *           dim         IN      Dimension of the dataset
 *          *dimensions  IN      1D array of length dim that stores the size of each dimension
 *          *data        IN      Pointer pointing to the data to be written 
 *                              
 * Returns: A negative value if an error occurs
 *          
 * Errors: Error occurs if the input is not valid, if dataset already exists, or an error is raised
 *         when writing the data.
 *        
 */
herr_t data_write_float(hid_t sid, const char *name, int dim, const int *dimensions, float *data)
{
    herr_t status = DANU_FAILURE;

    hsize_t * sizes;

    hid_t  gid;
    int    exists;

    if ( H5_RETURN_OK(data_exists(sid,name,&exists)) ) {
      
        if ( !exists ) {
            gid = data_open_group(sid);
            sizes = DANU_MALLOC(hsize_t,dim);
            convert_int_to_hsize(dim,dimensions,sizes);
            status = danu_data_write_float(gid,name,dim,sizes,data);
            DANU_FREE(sizes);
        }
        else {
            DANU_ERROR_MESS("Non series dataset already exists will not overwrite");
        }
    }
    else {
        DANU_ERROR_MESS("Failed to stat the non-series dataset");
    }

    return status;

}

/*
 * Routine: data_read_float(hid_t sid, const char *name, int dim, int *dimensions, float *data) 
 * Purpose: Read non_series integer data found in the non_series subgroup of simulation sid
 * Description: Read a dataset of dimensions dimensions under non_series data subgroup
 *              with label name.
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for the simulation group
 *          *name        IN      Dataset name 
 *           dim         IN      Dimension of the dataset
 *          *dimensions  IN      1D array of length dim that stores the size of each dimension
 *          *data        IN      Pointer pointing to the buffer to be read in to 
 *                              
 * Returns: A negative value if an error occurs
 *          
 * Errors: Error occurs if the input is not valid, if dataset does not exists, or an error is raised
 *         when writing the data.
 *        
 */
herr_t data_read_float(hid_t sid, const char *name, int dim, const int *dimensions, float *data)
{
    herr_t status = DANU_FAILURE;

    hid_t    gid;

    hsize_t * sizes;
    int      exists;

    if ( H5_RETURN_OK(data_exists(sid,name,&exists) ) ) {

        if ( exists ) {
            gid = data_open_group(sid);
            sizes = DANU_MALLOC(hsize_t,dim);
            convert_int_to_hsize(dim,dimensions,sizes);
            status =  danu_data_read_float(gid,name,dim,sizes,data);
            danu_group_close(gid);
            DANU_FREE(sizes);
        }
        else {
            DANU_ERROR_MESS("Non-series dataset does not exist");
        }

    }

    return status;

}
