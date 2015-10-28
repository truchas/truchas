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
 * dataset.c
 *
 *  DANU Datasets
 *
 *
 *  Purpose:
 *
 *          The HDF5 library requires a location and a dataspace definition
 *          for each dataset. The Danu library only sees groups and datasets.
 *          The extra object required to create a dataset are stored in a 
 *          single struct ddataset_t. Each dataset struct will contain 
 *          HDF5 identifiers for the location and dataspace associated with 
 *          the dataset. It will also contain default attributes that will be
 *          defined for any data set.
 *
 *          To avoid internal copies, the raw data pointers are copies of
 *          the pointers to the raw data. Thus the 
 *
 *
 */

#if HAVE_CONFIG_H
# include <danu_config.h>
#endif

#include <hdf5.h>

#include <danu_error.h>
#include <danu_types.h>
#include <danu_utils.h>
#include <danu_memory.h>
#include <danu_link.h>

#include <danu_dataset_types.h>
#include <danu_dataset.h>

/* Global define's */

#define DATASET_SMALL  32768             /* Dataset smaller than 32k are stored in a COMPACT style, with meta-data */
#define DATASET_BIG    104857600         /* Datasets larger than this value are chunked (Not supported) */
#define DATASET_ORDER  H5T_ORDER_LE      /* All data stored on-disk is order little endian */


/* Private functions */
hsize_t compute_dataset_size(int dim, const hsize_t *dims, hid_t type);




hsize_t compute_dataset_size(int dim, const hsize_t *dims, hid_t type)
{
    hsize_t i;
    hsize_t num;
    size_t tbytes = H5Tget_size(type);
    
    num = 1;
    for(i=0; i<dim; i++) {
        num*=dims[i];
    }

    num*= (hsize_t)tbytes;

    return num;
}

/*
 * Routine: int danu_dataset_rank(hid_t loc,const char *name)
 * Purpose:  Return the dimensions of a dataset
 * Description:
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *
 * Returns: Returns a nonzero postive value or negative number if 
 *          an error is detected.
 * Errors: Checks input and returns a -1 if an error is detected.
 *
*/
int danu_dataset_rank(hid_t loc, const char *name)
{

  hid_t data, space_id;
  int   ndim;

  /* Error checking by opening the dataset and dataspace
     Return with invalid value if either fails */
  ndim = -1;
  data = danu_dataset_open(loc,name);
  if ( H5_ISA_INVALID_ID(data) ) {
    DANU_ERROR_MESS("Failed to open dataset");
    return ndim;
  }

  space_id = H5Dget_space(data);
  if ( H5_ISA_INVALID_ID(data) ) {
    DANU_ERROR_MESS("Failed to retrieve dataspace id");
    return ndim;
  }

  /* This is for simple datasets only! */
  if ( H5Sis_simple(space_id) ) {
    ndim = H5Sget_simple_extent_ndims(space_id);
  }
  else {
    DANU_ERROR_MESS("Complex HDF5 datasets are not supported at this time");
  }

  /* Close all the HIDs */
  H5Sclose(space_id);
  H5Dclose(data);

  return ndim;

}
/*
 * Routine: int danu_dataset_rank(hid_t id)
 * Purpose:  Return the dimensions of a dataset
 * Description:
 *
 * Parameters:
 *           id           IN              Dataset HDF5 identifier
 *
 * Returns: Returns a nonzero postive value or negative number if 
 *          an error is detected.
 * Errors: Checks input and returns a -1 if an error is detected.
 *
*/
int danu_dataset_rank2(hid_t id)
{

  hid_t space_id;
  int   ndim;

  /* Error checking by opening the dataset and dataspace
     Return with invalid value if either fails */
  ndim = -1;
  if ( H5_ISA_INVALID_ID(id) ) {
    DANU_ERROR_MESS("Invalid HDF5 identifer");
    return ndim;
  }

  space_id = H5Dget_space(id);
  if ( H5_ISA_INVALID_ID(space_id) ) {
    DANU_ERROR_MESS("Failed to retrieve dataspace id");
    return ndim;
  }

  /* This is for simple datasets only! */
  if ( H5Sis_simple(space_id) ) {
    ndim = H5Sget_simple_extent_ndims(space_id);
  }
  else {
    DANU_ERROR_MESS("Complex HDF5 datasets are not supported at this time");
  }

  /* Close all the HIDs */
  H5Sclose(space_id);

  return ndim;

}
/*
 * Routine: herr_t danu_dataset_dimensions(hid_t loc, const char *name, 
                                           int size, hsize_t *dimensions)
 * Purpose:  Return the dimensions of a dataset
 * Description: Every HDF5 dataset has two 1D dimension arrays 
 *              associated with it. One is the current size of each dimension
 *              and the other is the maximum size in each dimension. This
 *              routine returns the current dimension size array. Calling routine must 
 *              also pass in the size of the returning array. An error
 *              is raised if the rank of the dataset is larger than the 
 *              length of the dimension array passed in.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           size          IN              Length of the 1D array dimensions
 *          *dimensions    OUT             Array holding the dimensions of the 
 *                                         dataset.
 *
 * Returns: Returns a postive value if no error is detected.
 * Errors: Checks input and returns immediately if an error is detected.
 *         Only simple datasets are supported at this time.
 *
*/
herr_t danu_dataset_dimensions(hid_t loc, const char *name,
			       int size, hsize_t *dimensions)
{

  hid_t data, space_id;
  herr_t status = DANU_FAILURE;

  int ndim;

  /* Error checking by opening the dataset and dataspace
     Return with invalid value if either fails */
  data = danu_dataset_open(loc,name);
  if ( H5_ISA_INVALID_ID(data) ) {
    DANU_ERROR_MESS("Failed to open dataset");
    return status;
  }

  space_id = H5Dget_space(data);
  if ( H5_ISA_INVALID_ID(data) ) {
    DANU_ERROR_MESS("Failed to retrieve dataspace id");
    return status;
  }

  /* Check the pointer */
  if ( DANU_BAD_PTR(dimensions) ) {
    DANU_ERROR_MESS("Invalid pointer to dimension array");
    return status;
  }

  /* Only simple datasets are supported */
  if ( ! H5Sis_simple(space_id) ) {
    DANU_ERROR_MESS("Complex datasets are not supported at this time");
    return status;
  }

  /* Check the number of dimensions */
  ndim = H5Sget_simple_extent_ndims(space_id);

  if ( ndim < 0 ) {
    DANU_ERROR_MESS("Failed to retrieve the number of dimensions");
    return status;
  }

  if ( ndim > size ) {
    DANU_ERROR_MESS("The dimensions array is too small");
    return status;
  }

  /* Now populate the dimensions array */
  ndim = H5Sget_simple_extent_dims(space_id,dimensions,NULL);

  /* Set the return status on the ndim value, >0 is OK */
  if ( ndim > 0 ) {
    status = DANU_SUCCESS;
  }
  else {
    DANU_ERROR_MESS("Failed to define the dimension array");
  }

  /* Close the space id */
  H5Sclose(space_id);

  return status;

}
/*
 * Routine: herr_t danu_dataset_dimensions2(hid_t id, 
                                           int size, hsize_t *dimensions)
 * Purpose:  Return the dimensions of a dataset
 * Description: Every HDF5 dataset has two 1D dimension arrays 
 *              associated with it. One is the current size of each dimension
 *              and the other is the maximum size in each dimension. This
 *              routine returns the current dimension size array. Calling routine must 
 *              also pass in the size of the returning array. An error
 *              is raised if the rank of the dataset is larger than the 
 *              length of the dimension array passed in.
 *
 * Parameters:
 *           id            IN              Dataset HDF5 identifier
 *           size          IN              Length of the 1D array dimensions
 *          *dimensions    OUT             Array holding the dimensions of the 
 *                                         dataset.
 *
 * Returns: Returns a postive value if no error is detected.
 * Errors: Checks input and returns immediately if an error is detected.
 *         Only simple datasets are supported at this time.
 *
*/
herr_t danu_dataset_dimensions2(hid_t id,
			        int size, hsize_t *dimensions)
{

  hid_t space_id;
  herr_t status = DANU_FAILURE;

  int ndim;

  /* Error checking by opening the dataset and dataspace
     Return with invalid value if either fails */
  if ( H5_ISA_INVALID_ID(id) ) {
    DANU_ERROR_MESS("Invalid HDF5 identifier");
    return status;
  }

  space_id = H5Dget_space(id);
  if ( H5_ISA_INVALID_ID(space_id) ) {
    DANU_ERROR_MESS("Failed to retrieve dataspace id");
    return status;
  }

  /* Check the pointer */
  if ( DANU_BAD_PTR(dimensions) ) {
    DANU_ERROR_MESS("Invalid pointer to dimension array");
    return status;
  }

  /* Only simple datasets are supported */
  if ( ! H5Sis_simple(space_id) ) {
    DANU_ERROR_MESS("Complex datasets are not supported at this time");
    return status;
  }

  /* Check the number of dimensions */
  ndim = H5Sget_simple_extent_ndims(space_id);

  if ( ndim < 0 ) {
    DANU_ERROR_MESS("Failed to retrieve the number of dimensions");
    return status;
  }

  if ( ndim > size ) {
    DANU_ERROR_MESS("The dimensions array is too small");
    return status;
  }

  /* Now populate the dimensions array */
  ndim = H5Sget_simple_extent_dims(space_id,dimensions,NULL);

  /* Set the return status on the ndim value, >0 is OK */
  if ( ndim > 0 ) {
    status = DANU_SUCCESS;
  }
  else {
    DANU_ERROR_MESS("Failed to define the dimension array");
  }

  /* Close the space id */
  H5Sclose(space_id);

  return status;

}
/*
 * Routine: herr_t danu_dataset_max_dimensions(hid_t loc, const char *name, 
                                           int size, hsize_t *max_dimensions)
 * Purpose:  Return the dimensions of a dataset
 * Description: Every HDF5 dataset has two 1D dimension arrays 
 *              associated with it. One is the current size of each dimension
 *              and the other is the maximum size in each dimension. This
 *              routine returns the maximum dimension size array. Entries in this
 *              array that have the value of H5S_UNLIMITED correspond to
 *              dimensions with unlimited size. 
 *              Calling routines must also pass in the size of the 
 *              returning array. An error is raised if the rank of the 
 *              dataset is larger than the length of the dimension array passed in.
 *
 * Parameters:
 *           loc               IN              HDF5 location identifier
 *           name              IN              Dataset name
 *           size              IN              Length of the 1D array dimensions
 *          *max_dimensions    OUT             Array holding the max dimensions of the 
 *                                              dataset.
 *
 * Returns: A positive value if no error is detected. 
 * Errors: Checks input and returns immediately if an error in the input
 *         is detected. Only supports simple datasets at this time.
 *
*/
herr_t danu_dataset_max_dimensions(hid_t loc, const char *name,
			       int size, hsize_t *max_dimensions)
{

  hid_t data, space_id;
  herr_t status = DANU_FAILURE;

  int ndim;

  /* Error checking by opening the dataset and dataspace
     Return with invalid value if either fails */
  data = danu_dataset_open(loc,name);
  if ( H5_ISA_INVALID_ID(data) ) {
    DANU_ERROR_MESS("Failed to open dataset");
    return status;
  }

  space_id = H5Dget_space(data);
  if ( H5_ISA_INVALID_ID(data) ) {
    DANU_ERROR_MESS("Failed to retrieve dataspace id");
    return status;
  }

  /* Check the pointer */
  if ( DANU_BAD_PTR(max_dimensions) ) {
    DANU_ERROR_MESS("Invalid pointer to max dimension array");
    return status;
  }

  /* Only simple datasets are supported */
  if ( ! H5Sis_simple(space_id) ) {
    DANU_ERROR_MESS("Complex datasets are not supported at this time");
    return status;
  }

  /* Check the number of dimensions */
  ndim = H5Sget_simple_extent_ndims(space_id);

  if ( ndim < 0 ) {
    DANU_ERROR_MESS("Failed to retrieve the number of dimensions");
    return status;
  }

  if ( ndim > size ) {
    DANU_ERROR_MESS("The max dimensions array is too small");
    return status;
  }

  /* Now populate the dimensions array */
  ndim = H5Sget_simple_extent_dims(space_id,NULL,max_dimensions);

  /* Set the return status on the ndim value, >0 is OK */
  if ( ndim > 0 ) {
    status = DANU_SUCCESS;
  }
  else {
    DANU_ERROR_MESS("Failed to define the max dimension array");
  }

  /* Close the space id */
  H5Sclose(space_id);

  return status;

}
/*
 * Routine: herr_t danu_dataset_max_dimensions2(hid_t id, 
                                                int size, hsize_t *max_dimensions)
 * Purpose:  Return the dimensions of a dataset
 * Description: Every HDF5 dataset has two 1D dimension arrays 
 *              associated with it. One is the current size of each dimension
 *              and the other is the maximum size in each dimension. This
 *              routine returns the maximum dimension size array. Entries in this
 *              array that have the value of H5S_UNLIMITED correspond to
 *              dimensions with unlimited size. 
 *              Calling routines must also pass in the size of the 
 *              returning array. An error is raised if the rank of the 
 *              dataset is larger than the length of the dimension array passed in.
 *
 * Parameters:
 *           id                IN              Dataset HDF5 location identifier
 *           size              IN              Length of the 1D array dimensions
 *          *max_dimensions    OUT             Array holding the max dimensions of the 
 *                                              dataset.
 *
 * Returns: A positive value if no error is detected. 
 * Errors: Checks input and returns immediately if an error in the input
 *         is detected. Only supports simple datasets at this time.
 *
*/
herr_t danu_dataset_max_dimensions2(hid_t id,
			            int size, hsize_t *max_dimensions)
{

  hid_t space_id;
  herr_t status = DANU_FAILURE;

  int ndim;

  if ( H5_ISA_INVALID_ID(id) ) {
    DANU_ERROR_MESS("Invalid HDF5 identifier");
    return status;
  }

  space_id = H5Dget_space(id);
  if ( H5_ISA_INVALID_ID(space_id) ) {
    DANU_ERROR_MESS("Failed to retrieve dataspace id");
    return status;
  }

  /* Check the pointer */
  if ( DANU_BAD_PTR(max_dimensions) ) {
    DANU_ERROR_MESS("Invalid pointer to max dimension array");
    return status;
  }

  /* Only simple datasets are supported */
  if ( ! H5Sis_simple(space_id) ) {
    DANU_ERROR_MESS("Complex datasets are not supported at this time");
    return status;
  }

  /* Check the number of dimensions */
  ndim = H5Sget_simple_extent_ndims(space_id);

  if ( ndim < 0 ) {
    DANU_ERROR_MESS("Failed to retrieve the maximum number of dimensions");
    return status;
  }

  if ( ndim > size ) {
    DANU_ERROR_MESS("The max dimensions array is too small");
    return status;
  }

  /* Now populate the dimensions array */
  ndim = H5Sget_simple_extent_dims(space_id,NULL,max_dimensions);

  /* Set the return status on the ndim value, >0 is OK */
  if ( ndim > 0 ) {
    status = DANU_SUCCESS;
  }
  else {
    DANU_ERROR_MESS("Failed to define the max dimension array");
  }

  /* Close the space id */
  H5Sclose(space_id);

  return status;

}

/*
 * Routine: herr_t danu_dataset_type(hid_t loc, const char * dataname, int* type)
 * Purpose:  Return a dataset type located under object loc with name dataname
 * Description: Given a HDF5 object loc and and dataset name, return the 
 *              dataset type. Type is set to one of the
 *              dataset types defined in danu_dataset_types.h.
 *              Will raise an error if id is invalid or a type
 *              can not be termined. Returns a negative value
 *              if an error is detected. Routine only recognizes
 *              int, float, double and string dataset types.
 *
 * Parameters:
 *           loc           IN             HDF5 object identifier
 *           dataname      IN             Dataset name
 *           *type        OUT             Dataset type code. See
 *                                        danu_dataset_types.h for possible
 *                                        values.
 *
 * Returns: Returns a negative number if an error is detected.
 *        
 * Errors: Error is raised if loc is not a valid id or the internal
 *         HDF5 calls return an error, or the type is not supported. 
 *
*/
herr_t danu_dataset_type(hid_t loc, const char *dataname, int *typecode)
{
  hid_t id;
  herr_t status=DANU_FAILURE;
  hid_t data_type;
  H5T_class_t class;
  size_t bytes;

  *typecode=DANU_DATASET_UNKNOWN;

  /* 
     HDF5 does not have a simple datatype query function. 
     Here we open the datatype and then check it against the 
     class type. We assume INTEGER class is always int and 
     Split float/double types based on size.
  */
  id=danu_dataset_open(loc,dataname);
  if ( H5_ISA_DATASET_ID(id) ) {
    data_type=H5Dget_type(id);
    if ( H5_ISA_VALID_ID(data_type) ) {
      class=H5Tget_class(data_type);
      if ( class == H5T_STRING ) {
	*typecode=DANU_DATASET_STRING;
      }
      else if ( class == H5T_INTEGER ) {
	*typecode=DANU_DATASET_INT;
      }
      else if ( class == H5T_FLOAT ) {
        bytes=H5Tget_size(data_type);
	if ( bytes == sizeof(float) ) {
	  *typecode=DANU_DATASET_FLOAT;
	}
	else if ( bytes == sizeof(double) ) {
	  *typecode=DANU_DATASET_DOUBLE;
	}
	else {
	  DANU_ERROR_MESS("Unsupported float class datatype");
	}
      }
      else {
	DANU_ERROR_MESS("Failed to determine class datatype.");
      }
      H5Tclose(data_type);
    }
    else {
      DANU_ERROR_MESS("Failed to open dataset type.");
    }
  }
  else {
    DANU_ERROR_MESS("Invalid dataset type.");
  }

  /* Set the status */
  if ( *typecode == DANU_DATASET_UNKNOWN ) {
    DANU_ERROR_MESS("Failed to determine dataset type.");
  }
  else {
    status=DANU_SUCCESS;
  }

  return status;

}
/*
 * Routine: herr_t danu_dataset_type2(hid_t id, int* type)
 * Purpose:  Return a dataset type located under object loc with name dataname
 * Description: Given a HDF5 object loc and and dataset name, return the 
 *              dataset type. Type is set to one of the
 *              dataset types defined in danu_dataset_types.h.
 *              Will raise an error if id is invalid or a type
 *              can not be termined. Returns a negative value
 *              if an error is detected. Routine only recognizes
 *              int, float, double and string dataset types.
 *
 * Parameters:
 *           loc           IN             Dataset HDF5  identifier
 *           *type        OUT             Dataset type code. See
 *                                        danu_dataset_types.h for possible
 *                                        values.
 *
 * Returns: Returns a negative number if an error is detected.
 *        
 * Errors: Error is raised if loc is not a valid id or the internal
 *         HDF5 calls return an error, or the type is not supported. 
 *
*/
herr_t danu_dataset_type2(hid_t id, int *typecode)
{
  herr_t status=DANU_FAILURE;
  hid_t data_type;
  H5T_class_t class;
  size_t bytes;

  *typecode=DANU_DATASET_UNKNOWN;

  /* 
     HDF5 does not have a simple datatype query function. 
     Here we open the datatype and then check it against the 
     class type. We assume INTEGER class is always int and 
     Split float/double types based on size.
  */
  if ( H5_ISA_DATASET_ID(id) ) {
    data_type=H5Dget_type(id);
    if ( H5_ISA_VALID_ID(data_type) ) {
      class=H5Tget_class(data_type);
      if ( class == H5T_STRING ) {
	*typecode=DANU_DATASET_STRING;
      }
      else if ( class == H5T_INTEGER ) {
	*typecode=DANU_DATASET_INT;
      }
      else if ( class == H5T_FLOAT ) {
        bytes=H5Tget_size(data_type);
	if ( bytes == sizeof(float) ) {
	  *typecode=DANU_DATASET_FLOAT;
	}
	else if ( bytes == sizeof(double) ) {
	  *typecode=DANU_DATASET_DOUBLE;
	}
	else {
	  DANU_ERROR_MESS("Unsupported float class datatype");
	}
      }
      else {
	DANU_ERROR_MESS("Failed to determine class datatype.");
      }
      H5Tclose(data_type);
    }
    else {
      DANU_ERROR_MESS("Failed to open dataset type.");
    }
  }
  else {
    DANU_ERROR_MESS("Invalid dataset type.");
  }

  /* Set the status */
  if ( *typecode == DANU_DATASET_UNKNOWN ) {
    DANU_ERROR_MESS("Failed to determine dataset type.");
  }
  else {
    status=DANU_SUCCESS;
  }

  return status;

}

/*
 * Routine: hbool_t danu_dataset_exists(hid_t loc,const char *name)
 * Purpose:  Check the existence of dataset
 * Description:
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *
 * Returns: Returns a TRUE if dataset exists, FALSE otherwise.
 * Errors: Checks input and returns a FALSE if input is not valid.
 *
 */ 
hbool_t danu_dataset_exists(hid_t loc_id, const char * name)
{
    hbool_t  flag;
    herr_t ret;

    ret = danu_link_exists(loc_id,name,&flag);

    if ( ! H5_RETURN_OK(ret) ) {
        DANU_ERROR_MESS("Failed to status dataset");
    }
   
    return flag;
}


/*
 * Routine: hid_t danu_dataset_create(loc,name,type,size,dim,extend)
 * Purpose: Create a HDF5 dataset
 * Description: 
 *
 * Parameters:
 *           loc                IN              HDF5 id for data set location (file,group)
 *           name               IN              Name of the dataset
 *           type               IN              HDF5 data memory type
 *           size               IN              Pointer to array of dataset sizes
 *           dim                IN              Dimension of the dataset
 *           extend             IN              Boolean flag indicating if dataset is extendable
 *                              
 * Returns: Returns a HDF5 identifier
 * Errors:
 *
 */ 
hid_t danu_dataset_create(hid_t loc, const char * name, 
                                                  hid_t type,
                                                  int dim,
                                                  const hsize_t * size,
                                                  dbool_t extend)

{
        hid_t dataset;
        hid_t dataspace, datatype;
        hid_t lpl, cpl, apl;
        hsize_t * maxsize;
        hsize_t   nbytes;
        didx_t i;

        /* Check input */
        if ( H5_ISA_INVALID_ID(loc) ) {
                DANU_ERROR_MESS("Invalid HDF5 identifier");
                return H5I_INVALID_HID;
        }

        if ( DANU_BAD_PTR(name) ) {
                DANU_ERROR_MESS("Invalid name pointer");
                return H5I_INVALID_HID;
        }

        if ( DANU_EMPTY_STRING(name) ) {
                DANU_ERROR_MESS("Invalid name argument");
                return H5I_INVALID_HID;
        }

        if ( DANU_BAD_PTR(size) ) {
                DANU_ERROR_MESS("Invalid size array pointer");
                return H5I_INVALID_HID;
        }

        if ( dim <= 0 ) {
                DANU_ERROR_MESS("invalid dimension value");
                return H5I_INVALID_HID;
        }
        
        /* Create the dataspace */
        maxsize = DANU_MALLOC(hsize_t, dim);

        if ( extend ) {
                for(i=0;i<dim;i++)
                   maxsize[i] = H5S_UNLIMITED;
        }
        else {
                for(i=0;i<dim;i++)
                   maxsize[i] = size[i];
        }
        
        dataspace = H5Screate_simple(dim,size,maxsize);

        /* Create the parameter lists */
        lpl = H5Pcreate(H5P_LINK_CREATE);
        cpl = H5Pcreate(H5P_DATASET_CREATE);
        apl = H5Pcreate(H5P_DATASET_ACCESS);

        /* Set the property lists */ 
        if ( extend ) {
                H5Pset_layout(cpl,H5D_CHUNKED);
                H5Pset_chunk(cpl,dim,size);
        }
        else {
                nbytes = compute_dataset_size(dim,size,type);
                if ( nbytes <= DATASET_SMALL ) {
                   H5Pset_layout(cpl,H5D_COMPACT);
                }
                else {
                   H5Pset_layout(cpl,H5D_CONTIGUOUS);
                }
        }

        /* Create the file datatype (order) */
        datatype = H5Tcopy(type);
        H5Tset_order(datatype,DATASET_ORDER);

        /* Create the HDF5 dataset */
        dataset = H5Dcreate(loc,name,datatype,dataspace,lpl,cpl,apl);

        /* Close all the objects */
        H5Pclose(lpl);
        H5Pclose(cpl);
        H5Pclose(apl);
        H5Tclose(datatype);
        H5Sclose(dataspace);

        /* Free pointers */
        DANU_FREE(maxsize);


    return dataset;
}       
/*
 * Routine: hid_t danu_dataset_open(loc,name)
 * Purpose: Open an existing HDF5 dataset
 * Description: 
 *
 * Parameters:
 *           loc                IN              HDF5 id for data set location (file,group)
 *           name               IN              Name of the dataset
 *                              
 * Returns: Returns a HDF5 identifier
 * Errors: Checks the input and will raise an error if input is not valid. Will return a
 *         negative value if fails to open the dataset
 *
 */ 
hid_t danu_dataset_open(hid_t loc,const char *name)
{
        hid_t dataset = H5I_INVALID_HID;

        /* Check Input */
        if ( H5_ISA_INVALID_ID(loc) ) {
           DANU_ERROR_MESS("Invalid HDF5 location identifier");
           return dataset;
        }

        if ( DANU_BAD_STRING(name) ) {
           DANU_ERROR_MESS("Invalid name string argument");
           return dataset;
        }


        if ( danu_dataset_exists(loc,name) ) {
            dataset = H5Dopen(loc,name,H5P_DEFAULT);
        }
        else {
            DANU_ERROR_MESS("Dataset does not exist");
        }


        return dataset;
}

/*
 * Routine: hid_t danu_dataset_close(id)
 * Purpose: Close and existing HDf5 dataset
 * Description: Closes an HDF5 dataset object and the other objects
 *              associated with that dataset. At this time, HDF5
 *              does not provide a call to retrieve the link
 *              property list, thus this property list can not
 *              be closed here.
 *
 * Parameters:
 *           id                IN              HDF5 dataset id 
 *                              
 * Returns: Returns a negative value if not successful
 * Errors: Checks for a valid id and returns immediately 
 *
 */ 
herr_t danu_dataset_close(hid_t id)
{
        herr_t ret;

        if ( ! H5_ISA_DATASET_ID(id) ) {
            DANU_ERROR_MESS("Invalid HDF5 id argument");
            return DANU_FAILURE;
        }

        ret = H5Dclose(id);

        return ret;
}
/*
 * Routine: herr_t danu_dataset_read(hid_t id,hid_t mem_type, hid_t mem_space,void * buf)
 * Purpose: Reads an existing HDF5 dataset into a buffer
 * Description: Read a dataset into a buffer. The dataset's on-file 
 *              dataspace and file data type will be used. The
 *              default transfer properties are used.
 *
 *
 * Parameters:
 *           id                IN              HDF5 dataset id 
 *           mem_type          IN              HDF5 native datatype
 *           buf               IN              Data buffer
 *                              
 * Returns: Returns a negative value if not successful
 * Errors: Checks for a valid id and mem_type. Memory type must be a NATIVE
 *         to correctly preform the byte swapping. If either id or mem_type
 *         are not correct routine returns immediately. 
 *
 */ 
herr_t danu_dataset_read(hid_t id, 
                         dslab_t *slab, 
                         hid_t buf_type, 
                         int dim, 
                         const hsize_t * buf_size,
                         void *buf)
{
        herr_t status;
        hid_t  filespace,memspace;
        hid_t  memtype; 

        /* Check input */
        if ( H5_ISA_INVALID_ID(id) ) {
            DANU_ERROR_MESS("Invalid HDF5 identifier argument");
            return DANU_FAILURE;
        }
        
        if ( dim <= 0 ) {
            DANU_ERROR_MESS("Invalid dimension argument");
            return DANU_FAILURE;
        }

        if ( DANU_BAD_PTR(buf_size) ) {
            DANU_ERROR_MESS("Invalid buf_size pointer");
            return DANU_FAILURE;
        }

        if ( DANU_BAD_PTR(buf) ) {
            DANU_ERROR_MESS("Invalid buf pointer");
            return DANU_FAILURE;
        }

        /* Define the buffer memory space */
        memspace = H5Screate_simple(dim,buf_size,NULL);
        memtype  = H5Tget_native_type(buf_type,H5T_DIR_ASCEND);

        /* Define the file space 
        *  If slab is a valid pointer, then the data is
        *  a hyperslab
        */
        filespace = H5Dget_space(id);

        status = DANU_SUCCESS;
        if ( slab != NULL ) {
            status = danu_slab_select(filespace,slab);
        }

        if ( H5_RETURN_OK(status) || status == DANU_SUCCESS ) {
            status = H5Dread(id,memtype,memspace,filespace,H5P_DEFAULT,buf);
        }

        /* Close the HDF5 objects */
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Tclose(memtype);

        return status;
}
/*
 * Routine: herr_t danu_dataset_write(hid_t id,hid_t mem_type, void * buf)
 * Purpose: Write data in buffer to an existing HDF5 dataset
 * Description: Read a dataset into a buffer. The dataset's on-file 
 *              dataspace and file data type will be used. The
 *              default transfer properties are used.
 *
 *
 * Parameters:
 *           id                IN              HDF5 dataset id 
 *           mem_type          IN              The buffer's HDF5 datatype
 *           buf               IN              Data buffer
 *                              
 * Returns: Returns a negative value if not successful
 * Errors: Checks for a valid id and buffer pointer. Will return immediately
 *         if either are not valid.
 *
 */ 
herr_t danu_dataset_write(hid_t id,
                          dslab_t *slab,
                          hid_t buf_type,
                          int dim,
                          const hsize_t *buf_size,
                          const void *buf)
{
        herr_t status;
        hid_t filespace,memspace;
        hid_t memtype;
        
        /* Check input */
        if ( H5_ISA_INVALID_ID(id) ) {
            DANU_ERROR_MESS("Invalid HDF5 identifier");
            return DANU_FAILURE;
        }

        if ( dim <=0 ) {
            DANU_ERROR_MESS("Invalid dimension");
            return DANU_FAILURE;
        }

        if ( DANU_BAD_PTR(buf_size) ) {
            DANU_ERROR_MESS("Invalid size pointer");
            return DANU_FAILURE;
        }

        if ( DANU_BAD_PTR(buf) ) {
            DANU_ERROR_MESS("Invalid buffer pointer");
            return DANU_FAILURE;
        }


        /* Define the memory space */
        memspace = H5Screate_simple(dim,buf_size,buf_size);
        memtype  = H5Tget_native_type(buf_type, H5T_DIR_ASCEND);

        
        /* 
        *  Define the file space from the id passed in
        *  If slab is a valid pointer, then the data is a 
        *  hyperslab
        */
        filespace = H5Dget_space(id);

        status = DANU_SUCCESS;
        if ( slab != NULL ) {
            status = danu_slab_select(filespace,slab);
        }

        if ( H5_RETURN_OK(status) || status == DANU_SUCCESS ) {
            status = H5Dwrite(id,memtype,memspace,filespace,H5P_DEFAULT,buf);
        }

        /* Free the HDF5 objects */
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Tclose(memtype);


        return status;
}
/*
 * Routine: hid_t danu_dataset_create_byte(hid_t loc,const char *name, int dim, hsize_t *size )
 * Purpose: Create a dataset of opaque byte values
 * Description: Create a dataset of opaque bytes located in HDF5 object loc. Routine is 
 *              essentially a wrapper to danu_dataset_create.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input checked in danu_dataset_create. Error raised if input is not valid.
 *
 */ 
 hid_t danu_dataset_create_byte(hid_t loc, const char *name, int dim, const hsize_t *size)
 {
     return danu_dataset_create(loc,name,H5T_NATIVE_OPAQUE,dim,size,FALSE);
 }
/*
 * Routine: hid_t danu_dataset_create_double(hid_t loc,const char *name, int dim, hsize_t *size )
 * Purpose: Create a dataset of double values 
 * Description: Create a dataset of doubles located in HDF5 object loc. Routine is 
 *              essentially a wrapper to danu_dataset_create.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input checked in danu_dataset_create. Error raised if input is not valid.
 *
 */ 
 hid_t danu_dataset_create_double(hid_t loc, const char *name, int dim, const hsize_t *size)
 {
     return danu_dataset_create(loc,name,H5T_NATIVE_DOUBLE,dim,size,FALSE);
 }
/*
 * Routine: hid_t danu_dataset_create_float(hid_t loc,const char *name, int dim, hsize_t *size )
 * Purpose: Create a dataset of float values 
 * Description: Create a dataset of floats located in HDF5 object loc. Routine is 
 *              essentially a wrapper to danu_dataset_create.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input checked in danu_dataset_create. Error raised if input is not valid.
 *
 */ 
 hid_t danu_dataset_create_float(hid_t loc, const char *name, int dim, const hsize_t *size)
 {
     return danu_dataset_create(loc,name,H5T_NATIVE_FLOAT,dim,size,FALSE);
 }

/*
 * Routine: hid_t danu_dataset_create_int(hid_t loc,const char *name, int dim, hsize_t *size )
 * Purpose: Create a dataset of integer values 
 * Description: Create a dataset of integers located in HDF5 object loc. Routine is 
 *              essentially a wrapper to danu_dataset_create.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input checked in danu_dataset_create. Error raised if input is not valid.
 *
 */ 
 hid_t danu_dataset_create_int(hid_t loc, const char *name, int dim, const hsize_t *size)
 {
     return danu_dataset_create(loc,name,H5T_NATIVE_INT,dim,size,FALSE);
 }
/*
 * Routine: hid_t danu_dataset_create_string(hid_t loc,const char *name, int dim, hsize_t *size )
 * Purpose: Create a dataset of  strings
 * Description: Create a dataset of strings located in HDF5 object loc. Routine is 
 *              essentially a wrapper to danu_dataset_create with a defined atomic datatype.
 *              Routine computes the number of bytes required from the size array and includes
 *              the NULL character string.
 *              
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input checked in danu_dataset_create. Error raised if input is not valid.
 *
 */ 
 hid_t danu_dataset_create_string(hid_t loc, const char *name, int dim, const hsize_t *size)
 {
     int i;
     hsize_t bytes;
     hid_t string_type, id;

     /* Compute total number of bytes with NULL's*/
     bytes = 0;
     for(i=0;i<dim;i++) {
         bytes+=(size[i]+1);
     }

     /* Create the data type */
     string_type = H5Tcopy(H5T_C_S1);
     H5Tset_size(string_type,bytes);
     H5Tset_strpad(string_type,H5T_STR_NULLTERM);

     id = danu_dataset_create(loc,name,string_type,dim,size,FALSE);

     H5Tclose(string_type);

     return id;
 }
/*
 * Routine: herr_t danu_data_write_byte(hid_t loc,const char *name, int dim, hsize_t *size, int8_t *buf)
 * Purpose: Write a double dataset containing data in buf
 * Description: Create a dataset of bytes located in HDF5 object loc and write the data found in buf to that
 *              dataset. Routine does NOT return the dataset identifier, only the status of the write. The 
 *              is closed before the routine returns. Calling routine must re-open dataset to retrieve the 
 *              identifier.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Some of the input is checked in danu_dataset_create. Buffer pointer checked in this routine
 *         Return immediately if input is not valid.
 *
 */ 
 herr_t danu_data_write_byte(hid_t loc, 
                               const char *name,
                               int dim,
                               const hsize_t *size,
                               const int8_t *buf)
 {
     hid_t id;
     herr_t status;

     /* Creaate the dataset */
     id = danu_dataset_create_byte(loc,name,dim,size);

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_write(id,NULL,H5T_NATIVE_OPAQUE,dim,size,buf);
         danu_dataset_close(id);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
}
/*
 * Routine: herr_t danu_data_write_double(hid_t loc,const char *name, int dim, hsize_t *size,double *buf )
 * Purpose: Write a double dataset containing data in buf
 * Description: Create a dataset of doubles located in HDF5 object loc and write the data found in buf to that
 *              dataset. Routine does NOT return the dataset identifier, only the status of the write. The 
 *              is closed before the routine returns. Calling routine must re-open dataset to retrieve the 
 *              identifier.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Some of the input is checked in danu_dataset_create. Buffer pointer checked in this routine
 *         Return immediately if input is not valid.
 *
 */ 
 herr_t danu_data_write_double(hid_t loc, 
                               const char *name,
                               int dim,
                               const hsize_t *size,
                               const double *buf)
 {
     hid_t id;
     herr_t status;

     /* Creaate the dataset */
     id = danu_dataset_create_double(loc,name,dim,size);

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_write(id,NULL,H5T_NATIVE_DOUBLE,dim,size,buf);
         danu_dataset_close(id);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
}
/*
 * Routine: herr_t danu_data_write_double2(hid_t id, int dim, hsize_t *size,double *buf )
 * Purpose: Write a double dataset containing data in buf
 * Description: Given a dataset id  and write the data found in buf to that
 *              dataset.  
 *
 * Parameters:
 *           id            IN              Dataset HDF5 identifier
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Some of the input is checked in danu_dataset_create. Buffer pointer checked in this routine
 *         Return immediately if input is not valid.
 *
 */ 
 herr_t danu_data_write_double2(hid_t id, 
                               int dim,
                               const hsize_t *size,
                               const double *buf)
 {
     herr_t status;

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_write(id,NULL,H5T_NATIVE_DOUBLE,dim,size,buf);
     }
     else {
         DANU_ERROR_MESS("Invalid HDF5 identifer");
         status = DANU_FAILURE;
     }


     return status;
}

/*
 * Routine: herr_t danu_data_write_float(hid_t loc,const char *name, int dim, hsize_t *size,float *buf )
 * Purpose: Write a float dataset containing data in buf
 * Description: Create a dataset of floats located in HDF5 object loc and write the data found in buf to that
 *              dataset. Routine does NOT return the dataset identifier, only the status of the write. The 
 *              is closed before the routine returns. Calling routine must re-open dataset to retrieve the 
 *              identifier.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Some of the input is checked in danu_dataset_create. Buffer pointer checked in this routine
 *         Return immediately if input is not valid.
 *
 */ 
 herr_t danu_data_write_float(hid_t loc, 
                               const char *name,
                               int dim,
                               const hsize_t *size,
                               const float *buf)
 {
     hid_t id;
     herr_t status;

     /* Creaate the dataset */
     id = danu_dataset_create_float(loc,name,dim,size);

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_write(id,NULL,H5T_NATIVE_FLOAT,dim,size,buf);
         danu_dataset_close(id);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
 }
/*
 * Routine: herr_t danu_data_write_float2(hid_t id, int dim, hsize_t *size,double *buf )
 * Purpose: Write a float dataset containing data in buf
 * Description: Given a dataset id  and write the data found in buf to that
 *              dataset.  
 *
 * Parameters:
 *           id            IN              Dataset HDF5 identifier
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Some of the input is checked in danu_dataset_create. Buffer pointer checked in this routine
 *         Return immediately if input is not valid.
 *
 */ 
 herr_t danu_data_write_float2(hid_t id, 
                               int dim,
                               const hsize_t *size,
                               const float *buf)
 {
     herr_t status;

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_write(id,NULL,H5T_NATIVE_FLOAT,dim,size,buf);
     }
     else {
         DANU_ERROR_MESS("Invalid HDF5 identifer");
         status = DANU_FAILURE;
     }


     return status;
}


/*
 * Routine: herr_t danu_data_write_int(hid_t loc,const char *name, int dim, hsize_t *size,int *buf )
 * Purpose: Write a integer dataset containing data in buf
 * Description: Create a dataset of integers located in HDF5 object loc and write the data found in buf to that
 *              dataset. Routine does NOT return the dataset identifier, only the status of the write. The dataset 
 *              is closed before the routine returns. Calling routine must re-open dataset to retrieve the 
 *              identifier.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_create and danu_dataset_write.
 *         
 *
 */ 
 herr_t danu_data_write_int(hid_t loc,
                            const char *name,
                            int dim,
                            const hsize_t *size,
                            const int *buf)
 {
     hid_t id;
     herr_t status;

     /* Creaate the dataset */
     id = danu_dataset_create_int(loc,name,dim,size);

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_write(id,NULL,H5T_NATIVE_INT,dim,size,buf);
         danu_dataset_close(id);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
 }
/*
 * Routine: herr_t danu_data_write_int2(hid_t id, int dim, hsize_t *size,double *buf )
 * Purpose: Write a integer dataset containing data in buf
 * Description: Given a dataset id  and write the data found in buf to that
 *              dataset.  
 *
 * Parameters:
 *           id            IN              Dataset HDF5 identifier
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Some of the input is checked in danu_dataset_create. Buffer pointer checked in this routine
 *         Return immediately if input is not valid.
 *
 */ 
 herr_t danu_data_write_int2(hid_t id, 
                               int dim,
                               const hsize_t *size,
                               const int *buf)
 {
     herr_t status;

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_write(id,NULL,H5T_NATIVE_INT,dim,size,buf);
     }
     else {
         DANU_ERROR_MESS("Invalid HDF5 identifer");
         status = DANU_FAILURE;
     }


     return status;
}
/*
 * Routine: herr_t danu_data_write_string(hid_t loc,const char *name, int dim, hsize_t *size,char *buf )
 * Purpose: Write a character dataset containing strings stored in buf
 * Description: Create a dataset of characters located in HDF5 object loc and write the data found in buf to that
 *              dataset. Create a data type to hold all the strings. Size of this datatype includes NULL characters
 *              and pads the strings with NULL terminators. 
 *              Routine does NOT return the dataset identifier, only the status of the write. The dataset 
 *              is closed before the routine returns. Calling routine must re-open dataset to retrieve the 
 *              identifier.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_create and danu_dataset_write.
 *         
 *
*/ 
herr_t danu_data_write_strings(hid_t loc, const char *name, int num, const char **buf)
{
     hid_t dataset, dataspace, strspace;
     hid_t type;
     hid_t cparams;
     hsize_t chunk_dims[1] = {2}; /* One char plus a NULL */
     hsize_t dims[1];

     herr_t status;

     /* Check the input */
     if ( H5_ISA_INVALID_ID(loc) ) {
         DANU_ERROR_MESS("Invalid HDF5 location identifier");
         return DANU_FAILURE;
     }

     if ( DANU_BAD_PTR(name) ) {
         DANU_ERROR_MESS("invalid pointer for dataset name");
         return DANU_FAILURE;
     }

     if ( DANU_EMPTY_STRING(name) ) {
         DANU_ERROR_MESS("Empty string name");
         return DANU_FAILURE;
     }

     if ( num <= 0 ) {
         DANU_ERROR_MESS("Invalid number of strings");
         return DANU_FAILURE;
     }

     /* Create the datatype */
     type = H5Tcopy(H5T_C_S1);
     H5Tset_size(type,H5T_VARIABLE);
     H5Tset_strpad(type,H5T_STR_NULLTERM);

     /* Create the dataspace */
     dims[0] = (hsize_t) num;
     dataspace = H5Screate_simple(1,dims,NULL);

     /* Create parameter list MUST have chunk turned on for variable length */
     cparams = H5Pcreate(H5P_DATASET_CREATE);
     H5Pset_chunk(cparams,1,chunk_dims);

     /* Create the dataset */
     dataset = H5Dcreate(loc,name,type,dataspace,H5P_DEFAULT,cparams,H5P_DEFAULT);

     if ( H5_ISA_VALID_ID(dataset) ) {
         
         /* Memory space for each string */
         strspace =  H5Screate_simple(1,dims,NULL);
         status = H5Dwrite(dataset,type,strspace,dataspace,H5P_DEFAULT,buf);

         H5Sclose(strspace);
     }
     else {
         status = DANU_FAILURE;
     }

     /* Close up objects */
     H5Tclose(type);
     H5Pclose(cparams);
     H5Sclose(dataspace);
     H5Dclose(dataset);

     return status;
}
/*
 * Routine: herr_t danu_data_write_string2(hid_t id, int dim, hsize_t *size,char *buf )
 * Purpose: Write a character dataset containing strings stored in buf
 * Description: Write character strings found in buf the data found in buf to 
 *              dataset with id. Size of this datatype includes NULL characters
 *              and pads the strings with NULL terminators. 
 *
 * Parameters:
 *           id            IN              Dataset HDF5 location identifier
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_create and danu_dataset_write.
 *         
 *
*/ 
herr_t danu_data_write_strings2(hid_t id, int num, const char **buf)
{
     hid_t dataspace, strspace;
     hid_t type;
     hid_t cparams;
     hsize_t chunk_dims[1] = {2}; /* One char plus a NULL */
     hsize_t dims[1];

     herr_t status;

     /* Check the input */
     if ( H5_ISA_INVALID_ID(id) ) {
         DANU_ERROR_MESS("Invalid HDF5 location identifier");
         return DANU_FAILURE;
     }

     if ( num <= 0 ) {
         DANU_ERROR_MESS("Invalid number of strings");
         return DANU_FAILURE;
     }

     /* Create the datatype */
     type = H5Tcopy(H5T_C_S1);
     H5Tset_size(type,H5T_VARIABLE);
     H5Tset_strpad(type,H5T_STR_NULLTERM);

     /* Create the dataspace */
     dims[0] = (hsize_t) num;
     dataspace = H5Dget_space(id);

     /* Create parameter list MUST have chunk turned on for variable length */
     cparams = H5Pcreate(H5P_DATASET_CREATE);
     H5Pset_chunk(cparams,1,chunk_dims);

         
     /* Memory space for each string */
     strspace =  H5Screate_simple(1,dims,NULL);
     status = H5Dwrite(id,type,strspace,dataspace,H5P_DEFAULT,buf);

     H5Sclose(strspace);

     /* Close up objects */
     H5Tclose(type);
     H5Pclose(cparams);
     H5Sclose(dataspace);

     return status;
}
/*
 * Routine: herr_t danu_data_read_byte(hid_t loc,const char *name, int dim, hsize_t *size, int8_t *buf )
 * Purpose: Read a opaque byte dataset and store data in buf
 * Description: Open a dataset of doubles located in HDF5 object loc and read the data into buf.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_open and danu_dataset_read. Errors are raised in these routines.
 *
 */ 
 herr_t danu_data_read_byte(hid_t loc, const char *name, int dim, const hsize_t *size, int8_t *buf)
 {
     hid_t id;
     herr_t status;

     /* Open the dataset */
     id = danu_dataset_open(loc,name);

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_read(id,NULL,H5T_NATIVE_OPAQUE,dim,size,buf);
         danu_dataset_close(id);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
 }
/*
 * Routine: herr_t danu_data_read_double(hid_t loc,const char *name, int dim, hsize_t *size,double *buf )
 * Purpose: Read a double dataset and store data in buf
 * Description: Open a dataset of doubles located in HDF5 object loc and read the data into buf.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_open and danu_dataset_read. Errors are raised in these routines.
 *
 */ 
 herr_t danu_data_read_double(hid_t loc, const char *name, int dim, const hsize_t *size, double *buf)
 {
     hid_t id;
     herr_t status;

     /* Open the dataset */
     id = danu_dataset_open(loc,name);

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_read(id,NULL,H5T_NATIVE_DOUBLE,dim,size,buf);
         danu_dataset_close(id);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
 }
/*
 * Routine: herr_t danu_data_read_double2(hid_t id, int dim, hsize_t *size,double *buf )
 * Purpose: Read a double dataset and store data in buf
 * Description: Given dataset id of doubles, read the data into buf.
 *
 * Parameters:
 *           id            IN              Dataset HDF5 identifier
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_open and danu_dataset_read. Errors are raised in these routines.
 *
 */ 
 herr_t danu_data_read_double2(hid_t id, int dim, const hsize_t *size, double *buf)
 {
     herr_t status;

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_read(id,NULL,H5T_NATIVE_DOUBLE,dim,size,buf);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
}
/*
 * Routine: herr_t danu_data_read_float(hid_t loc,const char *name, int dim, hsize_t *size,float *buf )
 * Purpose: Read a float dataset and store data in buf
 * Description: Open a dataset of floats located in HDF5 object loc and read the data into buf.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_open and danu_dataset_read. Errors are raised in these routines.
 *
 */ 
 herr_t danu_data_read_float(hid_t loc, const char *name, int dim, const hsize_t *size, float *buf)
 {
     hid_t id;
     herr_t status;

     /* Open the dataset */
     id = danu_dataset_open(loc,name);

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_read(id,NULL,H5T_NATIVE_FLOAT,dim,size,buf);
         danu_dataset_close(id);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
 }
/*
 * Routine: herr_t danu_data_read_float2(hid_t id, int dim, hsize_t *size,double *buf )
 * Purpose: Read a float dataset and store data in buf
 * Description: Given dataset id of floats, read the data into buf.
 *
 * Parameters:
 *           id            IN              Dataset HDF5 identifier
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_open and danu_dataset_read. Errors are raised in these routines.
 *
 */ 
 herr_t danu_data_read_float2(hid_t id, int dim, const hsize_t *size, float *buf)
 {
     herr_t status;

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_read(id,NULL,H5T_NATIVE_FLOAT,dim,size,buf);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
}
/*
 * Routine: herr_t danu_data_read_int(hid_t loc,const char *name, int dim, hsize_t *size,int *buf )
 * Purpose: Read an integer dataset and store data in buf
 * Description: Read a dataset of integers located in HDF5 object loc and store the data in buf.
 *              Routine does NOT return the dataset identifier, only the status of the read. The dataset 
 *              is closed before the routine returns. Calling routine must re-open dataset to retrieve the 
 *              identifier.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_create and danu_dataset_read.
 *         
 *
 */ 
 hid_t danu_data_read_int(hid_t loc, const char *name, int dim, const hsize_t *size, int *buf)
 {
     hid_t id;
     herr_t status;

     /* Creaate the dataset */
     id = danu_dataset_open(loc,name);

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_read(id,NULL,H5T_NATIVE_INT,dim,size,buf);
         danu_dataset_close(id);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
 }
/*
 * Routine: herr_t danu_data_read_int2(hid_t id, int dim, hsize_t *size, int *buf )
 * Purpose: Read an integer dataset and store data in buf
 * Description: Given dataset id of floats, read the data into buf.
 *
 * Parameters:
 *           id            IN              Dataset HDF5 identifier
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_open and danu_dataset_read. Errors are raised in these routines.
 *
 */ 
 herr_t danu_data_read_int2(hid_t id, int dim, const hsize_t *size, int *buf)
 {
     herr_t status;

     if ( H5_ISA_VALID_ID(id) ) {
         status = danu_dataset_read(id,NULL,H5T_NATIVE_INT,dim,size,buf);
     }
     else {
         status = DANU_FAILURE;
     }


     return status;
 }
/*
 * Routine: herr_t danu_data_read_string(hid_t loc,const char *name, int dim, hsize_t *size,char *buf )
 * Purpose: Read a character dataset and store the strings in buf
 * Description: Read a dataset of characters located in HDF5 object loc and store the data in buf.
 *              Routine does NOT return the dataset identifier, only the status of the write. The dataset 
 *              is closed before the routine returns. Calling routine must re-open dataset to retrieve the 
 *              identifier.
 *
 * Parameters:
 *           loc           IN              HDF5 location identifier
 *           name          IN              Dataset name
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_create and danu_dataset_write.
 *         
 *
*/ 
hid_t danu_data_read_strings(hid_t loc, const char *name, int num, char **buf)
{
     hid_t dataset, type;
     herr_t status;

     /* Check input */
     if ( H5_ISA_INVALID_ID(loc) ) {
         DANU_ERROR_MESS("Invalid HDF5 location identifier");
         return DANU_FAILURE;
     }

     if ( DANU_BAD_PTR(name) ) {
         DANU_ERROR_MESS("Bad pointer for name");
         return DANU_FAILURE;
     }

     if ( DANU_EMPTY_STRING(name) ) {
         DANU_ERROR_MESS("Empty name string");
         return DANU_FAILURE;
     }

     if ( num <= 0 ) {
         DANU_ERROR_MESS("Invalid number of strings");
         return DANU_FAILURE;
     }

     /* Create datatype */
     type = H5Tcopy(H5T_C_S1);
     H5Tset_size(type,H5T_VARIABLE);
     H5Tset_strpad(type,H5T_STR_NULLTERM);

     /* Open dataset */
     dataset = H5Dopen(loc,name,H5P_DEFAULT);

     if ( H5_ISA_VALID_ID(dataset) ) {
         status = H5Dread(dataset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,buf);
     }
     else {
         status = DANU_FAILURE;
     }

     /* Close objects */
     H5Tclose(type);
     H5Dclose(dataset);

     return status;
}
/*
 * Routine: herr_t danu_data_read_string2(hid_t id, int dim, hsize_t *size,char *buf )
 * Purpose: Read a character dataset and store the strings in buf
 * Description: Read a dataset of characters identified by id into buf.
 *
 * Parameters:
 *           id            IN              Dataset HDF5 location identifier
 *           dim           IN              Dataset dimension
 *           size          IN              1D array of length dim that defines the dataset size
 *           buf           IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked in danu_dataset_create and danu_dataset_write.
 *         
 *
*/ 
hid_t danu_data_read_strings2(hid_t id, int num, char **buf)
{
     hid_t type;
     herr_t status;

     /* Check input */
     if ( H5_ISA_INVALID_ID(id) ) {
         DANU_ERROR_MESS("Invalid HDF5 identifier");
         return DANU_FAILURE;
     }

     if ( num <= 0 ) {
         DANU_ERROR_MESS("Invalid number of strings");
         return DANU_FAILURE;
     }

     /* Create datatype */
     type = H5Tcopy(H5T_C_S1);
     H5Tset_size(type,H5T_VARIABLE);
     H5Tset_strpad(type,H5T_STR_NULLTERM);

     status = H5Dread(id,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,buf);

     /* Close objects */
     H5Tclose(type);

     return status;
}
#if 0
/*
 * Routine: herr_t danu_dataset_append(hid_t id, hid_t buf_type, const hsize_t * size, const void * buffer )
 * Purpose: Append data contained in buffer to an existing dataset
 * Description: Append data contained in buffer to an existing dataset. Dataset must be configured
 *              at it's creation to allow extensions. See HDF5 documentation for details.
 *              Existing dataset will be extended to accommodate the existing data and buffer.
 *              The routine appends the existing dataset to the first dimension.
 *
 *              Example: Current dataset is 3x3, buffer is 5x1. Dataset after this routine
 *              will be 8x3 with data from buffer written to position (3,0) through (7,0).
 *              x x x       x x x
 *              x x x   ->  x x x
 *              x x x       x x x
 *                          y 0 0
 *                          y 0 0
 *                          y 0 0
 *                          y 0 0
 *                          y 0 0
 *
 *             HDF5 zero-fills datasets. 
 *             Since there are several steps that must complete successfully for the data to be
 *             written, the routine's logic structure is a goto to improve the readability of the 
 *             code. The FAIL_EXIT label will return a negative value if an error is raised.
 *
 * Parameters:
 *           id           IN              HDF5 dataset identifier
 *           buf_type     IN              Buffer type
 *           size         IN              1D array that holds the dimensions of buffer
 *           buffer       IN              Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked and errors are raised if the input is not valid. Errors can also occur
 *         if dataspaces are not created correctly or if the write fails. HDF5 has restrictions
 *         on extending external datasets. In this case ONLY the first dimension is allowed to extend.
 *         This routine does not have the capability to check this before attempting to extend the
 *         dataset. See HDF5 error messages when this error occurs.
 *
*/
herr_t danu_dataset_append(hid_t id, hid_t buf_type, const hsize_t * size, const void * buffer)
{
    herr_t status = DANU_FAILURE;

    int i,rank;
    hsize_t  *old_size, *max_size, *new_size;
    hsize_t *offset=NULL;
    hid_t   mspace, dspace;

    /* Initialize the pointers to NULL */
    old_size=NULL;
    new_size=NULL;
    max_size=NULL;

    /* Check the input */
    if ( ! H5_ISA_DATASET_ID(id) ) {
        DANU_ERROR_MESS("Invalid dataset id");
        goto FAIL_EXIT;
    }

    if ( ! H5_ISA_DATATYPE_ID(buf_type) ) {
        DANU_ERROR_MESS("Invalid buffer data type id");
        goto FAIL_EXIT;
    }

    if ( DANU_BAD_PTR(size) ) {
        DANU_ERROR_MESS("Invalid size pointer");
        goto FAIL_EXIT;
    }

    if ( DANU_BAD_PTR(buffer) ) {
        DANU_ERROR_MESS("Invalid data buffer pointer");
        goto FAIL_EXIT;
    }

    /* Get the current space id */
    dspace = H5Dget_space(id);
    if ( H5_ISA_INVALID_ID(dspace) ) {
        DANU_ERROR_MESS("Failed to define the dataspace id");
        goto FAIL_EXIT;
    }

    /* Define the current dataspace size, close the space once we have the info we need */
    rank = H5Sget_simple_extent_ndims(dspace);
    old_size = DANU_MALLOC(hsize_t,rank);
    max_size = DANU_MALLOC(hsize_t,rank);
    new_size = DANU_MALLOC(hsize_t,rank);
    if ( H5Sget_simple_extent_dims(dspace,old_size,max_size) < 0 ) {
        DANU_ERROR_MESS("Failed to find current dataspace dimensions");
        goto FAIL_EXIT;
    }
    H5Sclose(dspace);

    /* Define the new_size */
    new_size[0] = old_size[0] + size[0];
    for(i=1;i<rank;i++) {
        if ( old_size[i] != H5S_UNLIMITED ) {
            new_size[i] = MAX(size[i],old_size[i]);
        }
        else {
            new_size[i] = H5S_UNLIMITED;
        }
    }

    /* Now check the new size */
    for(i=0;i<rank;i++) {
        if ( new_size[i] != H5S_UNLIMITED && new_size[i] > max_size[i] ) {
            DANU_ERROR_MESS("Can not extend the dataset past the configured bounds");
            goto FAIL_EXIT;
        }
    }

    /* Extend the dataset to the new size and open a new space with this new size */
    if ( H5Dset_extent(id,new_size) < 0 ) {
        DANU_ERROR_MESS("Failed to extend the dataset");
        goto FAIL_EXIT;
    }
    dspace = H5Dget_space(id);

    /* Create a hyperslab */
    offset = DANU_MALLOC(hsize_t,rank);
    offset[0] = old_size[0];
    for(i=1;i<rank;i++) {
        offset[i] = 0;
    }
    if ( H5Sselect_hyperslab(dspace,H5S_SELECT_SET,offset,NULL,size,NULL) < 0 ) {
        DANU_ERROR_MESS("Failed to select a hyperslab");
        goto FAIL_EXIT;
    }

    /* Create a memory space for the buffer */
    mspace = H5Screate_simple(rank,size,NULL);

    /* Now write the buffer */
    status = H5Dwrite(id,buf_type,mspace,dspace,H5P_DEFAULT,buffer);

    /* Clean up the HDF5 objects created */
    H5Sclose(mspace);
    H5Sclose(dspace);

    /* Clean up the memory alloc'd here */
    DANU_FREE(offset);
    DANU_FREE(old_size);
    DANU_FREE(new_size);
    DANU_FREE(max_size);

    return status;

FAIL_EXIT:

    /* Clean up the arrays */
    if(old_size) DANU_FREE(old_size);
    if(new_size) DANU_FREE(new_size);
    if(max_size) DANU_FREE(max_size);
    if(offset)   DANU_FREE(offset);

    return DANU_FAILURE;
}
/*
 * Routine: herr_t danu_data_append_int(hid_t loc, const char * name, int rank, const hsize_t * size, const int * buffer )
 * Purpose: Append integer data to an existing dataset
 * Description: This is a wrapper routine for danu_dataset_append for integer buffers.
 *              Routine checks the existence of the dataset and if it does not exist will simply create it with
 *              extend flag set to TRUE. Input values are checked in the subroutines called to create or
 *              append the data. See danu_dataset_write and danu_dataset_append for more information. 
 *
 * Parameters:
 *           loc           IN             HDF5 identifier for the location of the dataset
 *           name          IN             Name of the dataset
 *           rank          IN             Rank of the dataset
 *           size          IN             1D array of length rank containing dataset dimensions
 *           buffer        IN             Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked and errors are raised if the input is not valid. Most input checked in routines called
 *         here. will call danu_dataset_create,write with extend flag set to TRUE if dataset does not exist
 *
*/
herr_t danu_data_append_int(hid_t loc, const char *name, int rank, const hsize_t *size, const int *buffer)
{
   herr_t status = DANU_FAILURE;

   hid_t dataset;
   hbool_t append = danu_dataset_exists(loc,name);

   /* Check the location id  and name by defining the dataset id */
   if ( append ) {
       dataset = danu_dataset_open(loc,name);
   }
   else {
       dataset = danu_dataset_create(loc,name,H5T_NATIVE_INT,rank,size,TRUE);
   }

   if ( H5_ISA_INVALID_ID(dataset) ) {
       DANU_ERROR_MESS("Could not open/create dataset");
       return status;
   }

   
   /* Either write or append the data */
   if ( append ) {
       status = danu_dataset_append(dataset,H5T_NATIVE_INT,size,buffer);
   }
   else {
       status = danu_dataset_write(dataset,NULL,H5T_NATIVE_INT,rank,size,buffer);
   }

   danu_dataset_close(dataset);

   return status;
}
/*
 * Routine: herr_t danu_data_append_double(hid_t loc, const char * name, double rank, const hsize_t * size, const double * buffer )
 * Purpose: Append double data to an existing dataset
 * Description: This is a wrapper routine for danu_dataset_append for doubleeger buffers.
 *              Routine checks the existence of the dataset and if it does not exist will simply create it with
 *              extend flag set to TRUE. Input values are checked in the subroutines called to create or
 *              append the data. See danu_dataset_write and danu_dataset_append for more information. 
 *
 * Parameters:
 *           loc           IN             HDF5 identifier for the location of the dataset
 *           name          IN             Name of the dataset
 *           rank          IN             Rank of the dataset
 *           size          IN             1D array of length rank containing dataset dimensions
 *           buffer        IN             Buffer containing data
 *
 * Returns: Returns a negative value if not successful
 * Errors: Input is checked and errors are raised if the input is not valid. Most input checked in routines called
 *         here. will call danu_dataset_create,write with extend flag set to TRUE if dataset does not exist
 *
*/
herr_t danu_data_append_double(hid_t loc, const char *name, double rank, const hsize_t *size, const double *buffer)
{
   herr_t status = DANU_FAILURE;

   hid_t dataset;
   hbool_t append = danu_dataset_exists(loc,name);

   /* Check the location id  and name by defining the dataset id */
   if ( append ) {
       dataset = danu_dataset_open(loc,name);
   }
   else {
       dataset = danu_dataset_create(loc,name,H5T_NATIVE_DOUBLE,rank,size,TRUE);
   }

   if ( H5_ISA_INVALID_ID(dataset) ) {
       DANU_ERROR_MESS("Could not open/create dataset");
       return status;
   }

   
   /* Either write or append the data */
   if ( append ) {
       status = danu_dataset_append(dataset,H5T_NATIVE_INT,size,buffer);
   }
   else {
       status = danu_dataset_write(dataset,NULL,H5T_NATIVE_INT,rank,size,buffer);
   }

   danu_dataset_close(dataset);

   return status;
}
#endif
