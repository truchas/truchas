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
 * series_data.c
 *
 *  TOUT time series datasets
 *
 *
 *  Purpose:
 *           This source file defines functions that create time series datasets. 
 *
 *
 */

#if HAVE_CONFIG_H
#include <danu_config.h>
#endif

#include <stdarg.h>
#include <string.h>
#include <math.h>

#include <hdf5.h>

#include <danu_error.h>
#include <danu_memory.h>
#include <danu_utils.h>
#include <danu_linked_list.h>
#include <danu_link.h>
#include <danu_group.h>
#include <danu_dataset.h>
#include <danu_dataset_types.h>

#include <danu_sim.h>
#include <danu_series-group.h>
#include <danu_series-data.h>
/*
 * Routine: simulation_open_data(hid_t nsid, const char *datanamei)  
 * Purpose: Open and return the identifier for the dataset dataname under the
 *          Series/<seriesname> group.
 * Description: Open the dataset named dataname under the series group nsid
 *              and return the HDF5 identifier. Invalid identifier is returned
 *              if an error is detected.
 *
 * Parameters:
 *           nsid       IN              HDF5 identifier for existing sequence group
 *           dataname   IN              Name of the dataset to open
 *
 * Returns: The return code is an HDF5 identifier. An invalid identifier is
 *          returned if an error occurs.
 *          
 * Errors: An error is raised if the input values are not valid or if the
 *         dataset does not exist.
 *        
 */
hid_t simulation_open_data(hid_t nsid, const char *dataname)
{
  hid_t hid = H5I_INVALID_HID;
  
  hid = danu_dataset_open(nsid,dataname);

  return hid;

}
/*
 * Routine: simulation_data_exists(hid_t nsid, const char *dataname, int * exists 
 * Purpose: Query the existence of a dataset under series subgroup nsid
 * Description: Compares the name dataname against the existing datasets under
 *              Series/<seriesname> identified with nsid. Exist is initially set
 *              to FALSE and only set to TRUE if the dataset is found. The return
 *              of this routine indicates a successful search or failure. Return
 *              will be negative if an error is encountered with the input arguments,
 *              or traversing the datasets in the group. 
 *
 * Parameters:
 *           nsid       IN              HDF5 identifier for existing sequence group
 *           dataname   IN              Name of the dataset to search 
 *           exists     OUT             Flag set to TRUE or FALSE 
 *
 * Returns: The return code is negative if an error occurs.
 *          
 * Errors: An error is raised if the input values are not valid or the routine can not
 *         traverse the series group datasets. Not finding the dataset by itself is not
 *         considered an error. Calling routine should check the return code for error
 *         checking not the exist flag.
 *        
 */
herr_t simulation_data_exists(hid_t nsid, const char * dataname, int * exist)
{
    herr_t ret_value = DANU_FAILURE;

    hbool_t flag;

    ret_value = danu_link_exists(nsid,dataname,&flag);

    if ( H5_RETURN_OK(ret_value) ) {
        if (flag) {
            *exist = TRUE;
        }
        else {
            *exist = FALSE;
        }
    }
    else {
        DANU_ERROR_MESS("Failed to stat the dataset");
        *exist = FALSE;
    }
    
    return ret_value;

}
/*
 * Routine: simulation_data_rank(hid_t nsid, const char * dataname)  
 * Purpose: Find the rank (ndims) of a data set dataname under series group
 *          nsid
 * Description: Return the rank of dataname dataset found under the series
 *              group nsid. Will return a negative value if an error is
 *              detected.
 *
 * Parameters:
 *           nsid       IN         HDF5 identifier for Sereies #n group
 *           dataname   IN         Data set name 
 *                              
 * Returns: Returns a negative value if an error occurs, otherwise a positive
 *          number.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected.
 *        
 */
int simulation_data_rank(hid_t nsid, const char * dataname)
{
  int rank=-1;
  int exists;

  if ( H5_ISA_INVALID_ID(nsid) ) {
    DANU_ERROR_MESS("Invalid series group identifier");
    return rank;
  }

  if ( DANU_RETURN_OK(simulation_data_exists(nsid,dataname,&exists) ) ) {
    if (exists) {
      rank=danu_dataset_rank(nsid,dataname);
    }
    else {
      DANU_ERROR_MESS("Dataset does not exist");
    }
  }
  else {
    DANU_ERROR_MESS("Failed to stat dataset");
  }

  return rank;

}
/*
 * Routine: simulation_data_dimensions(hid_t nsid, const char * dataname,
 *                                     int size, hsize_t *dimensions)
 *
 * Purpose: Find the current size of dataname dataset under a series group.
 *              group nsid. Will return a negative value if an error is
 *              detected.
 *
 * Description: Return the current size in each dimension in array dimensions
 *              Routine checks input and will return a negative value if
 *              an error is raised. Return value should be checked before 
 *              accessing data in dimensions.
 *
 * Parameters:
 *           nsid        IN         HDF5 identifier for Sereies #n group
 *           dataname    IN         Data set name
 *           size        IN         Length of array dimension
 *           *dimensions OUT        Array containing the size of the dataset
 *                                  in each dimension.
 *                              
 * Returns: Returns a negative value if an error occurs, otherwise a zero.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected.
 *        
 */
herr_t simulation_data_dimensions(hid_t nsid, const char * dataname,
                               int size, hsize_t *dimensions)
{
  herr_t stat=DANU_FAILURE;
  int exists;

  if ( H5_ISA_INVALID_ID(nsid) ) {
    DANU_ERROR_MESS("Invalid series group identifier");
    return stat;
  }

  if ( size <= 0 ) {
    DANU_ERROR_MESS("Invalid array size value.");
    return stat;
  }

  if ( DANU_RETURN_OK(simulation_data_exists(nsid,dataname,&exists) ) ) {
    if (exists) {
      stat=danu_dataset_dimensions(nsid,dataname,size,dimensions);
    }
    else {
      DANU_ERROR_MESS("Dataset does not exist");
    }
  }
  else {
    DANU_ERROR_MESS("Failed to stat dataset");
  }

  return stat;

}
/*
 * Routine: simulation_data_type(hid_t nsid, const char * dataname,
 *                                   int *typecode)
 *
 * Purpose: Return the type of dataset dataname under the series group nsid.
 *
 * Description: Return either DANU_DATASET_INT,DOUBLE,FLOAT,STRING,UNKNOWN
 *              for dataset dataname. Input is checked and error is raised
 *              Error is also raised if the dataset type is not INT, DOUBLE,
 *              FLOAT or STRING.
 *
 * Parameters:
 *           nsid        IN         HDF5 identifier for Sereies #n group
 *           dataname    IN         Data set name
 *           *typecode   OUT        Data type code defined in
 *                                  danu_data_types.h 
 *                              
 * Returns: Returns a negative value if an error occurs, otherwise a zero.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected.
 *        
 */
herr_t simulation_data_type(hid_t nsid, const char * dataname,
                               int *typecode)
{
  herr_t stat=DANU_FAILURE;
  int exists;

  *typecode=DANU_DATASET_UNKNOWN;

  if ( H5_ISA_INVALID_ID(nsid) ) {
    DANU_ERROR_MESS("Invalid series group identifier");
    return stat;
  }

  if ( DANU_RETURN_OK(simulation_data_exists(nsid,dataname,&exists) ) ) {
    if (exists) {
      stat=danu_dataset_type(nsid,dataname,typecode);
    }
    else {
      DANU_ERROR_MESS("Dataset does not exist");
    }
  }
  else {
    DANU_ERROR_MESS("Failed to stat dataset");
  }

  return stat;

}

/*
 * Routine: simulation_data_count(hid_t nsid, int *ndata)  
 * Purpose: Find the number of series datasets under the series group nsid
 * Description: Returns the number of links found under the series group
 *              identified. Routine assumes all links found under the 
 *              series subgroup are series time-step datasets. Calling 
 *              routines should check the return status for errors, NOT nseries.
 *
 * Parameters:
 *           nsid       IN              HDF5 identifier for Sereies #n group
 *           *ndata     OUT             Number of series datasets found 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected.
 *        
 */

herr_t simulation_data_count(hid_t nsid, int *ndata)
{
    herr_t ret_value = DANU_FAILURE;

    hsize_t flag;

    ret_value = danu_group_get_nlinks(nsid,&flag);

    if ( H5_RETURN_OK(ret_value) ) {
        *ndata = (int) flag;
    }
    else {
        DANU_ERROR_MESS("Failed to determine link count");
    }

    return ret_value;
}
/*
 * Routine: simulation_data_list(hid_t nsid, int num, char **datanames)  
 * Purpose: Return the list of dataset names found under Series #n identified with nsid
 * Description: Returns the list of series sub group names found under Series #n subgroup
 *              Calling routine must provide the number
 *              (num) of names datanames can hold and the size of memory each pointer
 *              points to. The number of datasets can be found using the simulation_data_count.
 *              This routine allocates memory for each pointer in datanames. The 
 *              calling routine is responsible for freeing this memory.
 *
 * Parameters:
 *           nsid        IN      HDF5 identifier for Series #n group
 *           num         IN      Number of char pointers in datanames
 *          **datanames  OUT     Array holding pointers to the Sequence subgroups names 
 *          *num_found   OUT     Number of sequence subgroups found
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected.
 *        
*/
herr_t simulation_data_list(hid_t nsid, int num, char **datanames, int *num_found)
{
    herr_t ret_value = DANU_FAILURE;

    ret_value = danu_group_get_datasets(nsid,num,datanames,num_found);

    if ( ret_value < 0 ) {
        DANU_ERROR_MESS("Failed to read dataset names");
    }

    return ret_value;
}
/*
 * Routine: simulation_data_write_int(hid_t nsid, const char *dataname, 
                                        int num, const int *size, int *data)  
 * Purpose: Create and write a dataset of dimensions size under Series #n group with id nsid.
 * Description: Will create and write a dataset to Series #n group identified with nsid. Routine
 *              checks the existence of dataname before writing.
 *
 * Parameters:
 *           nsid        IN      HDF5 identifier for Series #n group
 *           dataname    IN      Dataset name
 *           num         IN      Dim of dataset
 *          *size        IN      Array holding the size of each dimension
 *          *data        OUT     Pointer to the data 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: An error is raised if the input is not valid, the dataset exists or the routine
 *         fails to write the data.
 *        
*/
herr_t simulation_data_write_int(hid_t nsid,
                                 const char * dataname,
                                 int num,
                                 const int *size,
                                 const int * data)
{
    herr_t ret_value = DANU_FAILURE;

    int exists;
    hsize_t *hsize;

    /* Check the existence of the dataset */
    ret_value = simulation_data_exists(nsid,dataname,&exists);
    if ( H5_RETURN_FAIL(ret_value) ) {
        DANU_ERROR_MESS("Failed to stat the existence of dataset");
        return ret_value;
    }

    if ( exists ) {
        DANU_ERROR_MESS("Will not overwrite a dataset that exists");
        return ret_value;
    }
   
    /* Convert the int to hsize_t */
    if ( NULL != ( hsize = DANU_MALLOC(hsize_t,num)) ) {
        convert_int_to_hsize(num,size,hsize);
        /* Write out data ... the remaining args are checked in this routine */
        ret_value = danu_data_write_int(nsid,dataname,num,hsize,data);

        if ( H5_RETURN_FAIL(ret_value) ) {
            DANU_ERROR_MESS("Failed to write series dataset");
        }

        DANU_FREE(hsize);
    }
    else {
        DANU_ERROR_MESS("Failed to allocate size array .. memory exhausted");
    }

    return ret_value;
}
/*
 * Routine: simulation_data_read_int(hid_t nsid, const char *dataname, 
                                        int num, const int *size, int *data)  
 * Purpose: Create and read a dataset of dimensions size under Series #n group with id nsid.
 * Description: Will create and read a dataset to Series #n group identified with nsid. Routine
 *              checks the existence of dataname before writing.
 *
 * Parameters:
 *           nsid        IN      HDF5 identifier for Series #n group
 *           dataname    IN      Dataset name
 *           num         IN      Dim of dataset
 *          *size        IN      Array holding the size of each dimension
 *          *data        OUT     Pointer to the data 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: An error is raised if the input is not valid, the dataset exists or the routine
 *         fails to read the data.
 *        
*/
herr_t simulation_data_read_int(hid_t nsid,
                                 const char * dataname,
                                 int num,
                                 const int *size,
                                 int * data)
{
    herr_t ret_value = DANU_FAILURE;

    int exists;
    hsize_t *hsize;

    /* Check the existence of the dataset */
    ret_value = simulation_data_exists(nsid,dataname,&exists);
    if ( H5_RETURN_FAIL(ret_value) ) {
        DANU_ERROR_MESS("Failed to stat the existence of dataset");
        return ret_value;
    }

    if ( ! exists ) {
        DANU_ERROR_MESS("Dataset does not exist in this series group");
        return ret_value;
    }
   
    /* Convert the int to hsize_t */
    if ( NULL != ( hsize = DANU_MALLOC(hsize_t,num) ) ) {
        convert_int_to_hsize(num,size,hsize);
        /* Read data ... the remaining args are checked in this routine */
        ret_value = danu_data_read_int(nsid,dataname,num,hsize,data);

        if ( H5_RETURN_FAIL(ret_value) ) {
            DANU_ERROR_MESS("Failed to read series dataset");
        }

        DANU_FREE(hsize);
    }

    return ret_value;
}
/*
 * Routine: simulation_data_write_double(hid_t nsid, const char *dataname, 
                                        int num, const int *size, double *data)  
 * Purpose: Create and write a dataset of dimensions size under Series #n group with id nsid.
 * Description: Will create and write a dataset to Series #n group identified with nsid. Routine
 *              checks the existence of dataname before writing.
 *
 * Parameters:
 *           nsid        IN      HDF5 identifier for Series #n group
 *           dataname    IN      Dataset name
 *           num         IN      Dim of dataset
 *          *size        IN      Array holding the size of each dimension
 *          *data        OUT     Pointer to the data 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: An error is raised if the input is not valid, the dataset exists or the routine
 *         fails to write the data.
 *        
*/
herr_t simulation_data_write_double(hid_t nsid,
                                    const char * dataname,
                                    int num,
                                    const int *size,
                                    const double * data)
{
    herr_t ret_value = DANU_FAILURE;

    int exists;
    hsize_t *hsize;

    /* Check the existence of the dataset */
    ret_value = simulation_data_exists(nsid,dataname,&exists);
    if ( H5_RETURN_FAIL(ret_value) ) {
        DANU_ERROR_MESS("Failed to stat the existence of dataset");
        return ret_value;
    }

    if ( exists ) {
        DANU_ERROR_MESS("Will not overwrite a dataset that exists");
        return ret_value;
    }
   
    /* Convert the int to hsize_t */
    if ( NULL != ( hsize = DANU_MALLOC(hsize_t, num) ) ) {
        convert_int_to_hsize(num,size,hsize);
        /* Write out data ... the remaining args are checked in this routine */
        ret_value = danu_data_write_double(nsid,dataname,num,hsize,data);

        if ( H5_RETURN_FAIL(ret_value) ) {
            DANU_ERROR_MESS("Failed to write series dataset");
        }

        DANU_FREE(hsize);
    }
    else {
        DANU_ERROR_MESS("Failed to allocate size array .. memory exhausted");
    }

    return ret_value;
}
/*
 * Routine: simulation_data_read_double(hid_t nsid, const char *dataname, 
                                        int num, const int *size, double *data)  
 * Purpose: Create and read a dataset of dimensions size under Series #n group with id nsid.
 * Description: Will create and read a dataset to Series #n group identified with nsid. Routine
 *              checks the existence of dataname before writing.
 *
 * Parameters:
 *           nsid        IN      HDF5 identifier for Series #n group
 *           dataname    IN      Dataset name
 *           num         IN      Dim of dataset
 *          *size        IN      Array holding the size of each dimension
 *          *data        OUT     Pointer to the data 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: An error is raised if the input is not valid, the dataset exists or the routine
 *         fails to read the data.
 *        
*/
herr_t simulation_data_read_double(hid_t nsid,
                                   const char * dataname,
                                   int num,
                                   const int *size,
                                   double * data)
{
    herr_t ret_value = DANU_FAILURE;

    int exists;
    hsize_t *hsize;

    /* Check the existence of the dataset */
    ret_value = simulation_data_exists(nsid,dataname,&exists);
    if ( H5_RETURN_FAIL(ret_value) ) {
        DANU_ERROR_MESS("Failed to stat the existence of dataset");
        return ret_value;
    }

    if ( ! exists ) {
        DANU_ERROR_MESS("Dataset does not exist in this series group");
        return ret_value;
    }
   
    /* Convert the int to hsize_t */
    if ( NULL != ( hsize = DANU_MALLOC(hsize_t,num) ) ) {
        convert_int_to_hsize(num,size,hsize);
        /* Read data ... the remaining args are checked in this routine */
        ret_value = danu_data_read_double(nsid,dataname,num,hsize,data);

        if ( H5_RETURN_FAIL(ret_value) ) {
            DANU_ERROR_MESS("Failed to read series dataset");
        }

        DANU_FREE(hsize);
    }

    return ret_value;
}

/*
 * Routine: simulation_data_write_float(hid_t nsid, const char *dataname, 
                                        int num, const int *size, float *data)  
 * Purpose: Create and write a dataset of dimensions size under Series #n group with id nsid.
 * Description: Will create and write a dataset to Series #n group identified with nsid. Routine
 *              checks the existence of dataname before writing.
 *
 * Parameters:
 *           nsid        IN      HDF5 identifier for Series #n group
 *           dataname    IN      Dataset name
 *           num         IN      Dim of dataset
 *          *size        IN      Array holding the size of each dimension
 *          *data        OUT     Pointer to the data 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: An error is raised if the input is not valid, the dataset exists or the routine
 *         fails to write the data.
 *        
*/
herr_t simulation_data_write_float(hid_t nsid,
                                    const char * dataname,
                                    int num,
                                    const int *size,
                                    const float * data)
{
    herr_t ret_value = DANU_FAILURE;

    int exists;
    hsize_t *hsize;

    /* Check the existence of the dataset */
    ret_value = simulation_data_exists(nsid,dataname,&exists);
    if ( H5_RETURN_FAIL(ret_value) ) {
        DANU_ERROR_MESS("Failed to stat the existence of dataset");
        return ret_value;
    }

    if ( exists ) {
        DANU_ERROR_MESS("Will not overwrite a dataset that exists");
        return ret_value;
    }
   
    /* Convert the int to hsize_t */
    if ( NULL != ( hsize = DANU_MALLOC(hsize_t, num) ) ) {
        convert_int_to_hsize(num,size,hsize);
        /* Write out data ... the remaining args are checked in this routine */
        ret_value = danu_data_write_float(nsid,dataname,num,hsize,data);

        if ( H5_RETURN_FAIL(ret_value) ) {
            DANU_ERROR_MESS("Failed to write series dataset");
        }

        DANU_FREE(hsize);
    }
    else {
        DANU_ERROR_MESS("Failed to allocate size array .. memory exhausted");
    }

    return ret_value;
}
/*
 * Routine: simulation_data_read_float(hid_t nsid, const char *dataname, 
                                        int num, const int *size, float *data)  
 * Purpose: Create and read a dataset of dimensions size under Series #n group with id nsid.
 * Description: Will create and read a dataset to Series #n group identified with nsid. Routine
 *              checks the existence of dataname before writing.
 *
 * Parameters:
 *           nsid        IN      HDF5 identifier for Series #n group
 *           dataname    IN      Dataset name
 *           num         IN      Dim of dataset
 *          *size        IN      Array holding the size of each dimension
 *          *data        OUT     Pointer to the data 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: An error is raised if the input is not valid, the dataset exists or the routine
 *         fails to read the data.
 *        
*/
herr_t simulation_data_read_float(hid_t nsid,
                                   const char * dataname,
                                   int num,
                                   const int *size,
                                   float * data)
{
    herr_t ret_value = DANU_FAILURE;

    int exists;
    hsize_t *hsize;

    /* Check the existence of the dataset */
    ret_value = simulation_data_exists(nsid,dataname,&exists);
    if ( H5_RETURN_FAIL(ret_value) ) {
        DANU_ERROR_MESS("Failed to stat the existence of dataset");
        return ret_value;
    }

    if ( ! exists ) {
        DANU_ERROR_MESS("Dataset does not exist in this series group");
        return ret_value;
    }
   
    /* Convert the int to hsize_t */
    if ( NULL != ( hsize = DANU_MALLOC(hsize_t,num) ) ) {
        convert_int_to_hsize(num,size,hsize);
        /* Read data ... the remaining args are checked in this routine */
        ret_value = danu_data_read_float(nsid,dataname,num,hsize,data);

        if ( H5_RETURN_FAIL(ret_value) ) {
            DANU_ERROR_MESS("Failed to read series dataset");
        }

        DANU_FREE(hsize);
    }

    return ret_value;
}

