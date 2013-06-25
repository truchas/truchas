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
 * attribute.c
 *
 *  DANU Attributes
 *
 *
 *  Purpose:
 *
 *          The HDF5 library allows attributes to be assigned to any HDF5 object
 *          This source file contains wrapper routines to write scalar attributes 
 *          to any valid HDF5 object.
 *
 *
 */

#if HAVE_CONFIG_H
# include <danu_config.h>
#endif

#include <string.h>

#include <hdf5.h>

#include <danu_utils.h>
#include <danu_memory.h>
#include <danu_h5_object.h>
#include <danu_linked_list.h>

#include <danu_attribute.h>

/* Global defines */
#define ATTRIBUTE_ORDER H5T_ORDER_LE  /* Attributes are stored in little endian format */

/* private function proto-types */
typedef struct attr_data_t {
    char * name;
    hid_t  type_id;
    hid_t  space_id;
} attr_data_t;

typedef struct attr_names_t {
    char  **names;
    int     idx;
    int     num;
} attr_names_t;


attr_data_t * create_attr_data(size_t num);
danu_err_t    delete_attr_data(void * data);

/* Functions use over the list iterator */
danu_err_t copy_attr_names(danu_node_t *node, void *op_data);

/* Functions used with H5Aiterate */
herr_t attr_simple_search(hid_t loc, const char * attr_name, const H5A_info_t * info, void * op_data);
herr_t attr_simple_count(hid_t loc, const char * attr_name, const H5A_info_t * info, void * op_data);


/* Private fnctions */
attr_data_t * create_attr_data(size_t num)
{
    attr_data_t *ptr = DANU_MALLOC(attr_data_t,1);

    if ( ptr != NULL ) {
        ptr->name = DANU_MALLOC(char,num);
    }
    
    return ptr;
}

danu_err_t delete_attr_data(void * data)
{
    attr_data_t *a_data = (attr_data_t *)data;

    DANU_FREE(a_data->name);
    DANU_FREE(a_data);

    return DANU_SUCCESS;
}
 
danu_err_t copy_attr_names(danu_node_t *node, void * op_data)
{
    attr_names_t * ndata = (attr_names_t *) op_data;
    attr_data_t  * adata = (attr_data_t *) node->data;

    size_t name_len = strlen(adata->name) + 1;

    ndata->names[ndata->idx] = DANU_MALLOC(char,name_len);
    strcpy(ndata->names[ndata->idx],adata->name);
    ndata->idx++;

    if ( ndata->idx < ndata->num ) {
	return DANU_CONTINUE;
    }
    else {
	return DANU_STOP;
    }
}

herr_t attr_simple_search(hid_t loc, const char * attr_name, const H5A_info_t *info, void *op_data)
{
    herr_t status = 0;

    danu_list_t * list = (danu_list_t *) op_data;
    attr_data_t * attr_data;
    size_t len;

    len = strlen(attr_name);
    attr_data = create_attr_data(len+1);
    sprintf(attr_data->name,attr_name);

    danu_list_append(list,attr_data,&delete_attr_data);

    return status;
}

herr_t attr_simple_count(hid_t loc, const char * attr_name, const H5A_info_t *info, void *op_data)
{
    herr_t status = 0;

    int *count = (int *) op_data;

    (*count)++;

    return status;
}


herr_t danu_attr_find_all(hid_t loc, danu_list_t *list)
{
    herr_t status = DANU_FAILURE;

    hsize_t n;

    /* Check the input */
    if ( DANU_BAD_PTR(list) ) {
        DANU_ERROR_MESS("Invalid Danu Linked List pointer");
        return status;
    }

    if ( H5_ISA_INVALID_ID(loc) ) {
        DANU_ERROR_MESS("Invalid HDF5 object identifier");
        return status;
    }

    n=0;
    status = H5Aiterate(loc,H5_INDEX_CRT_ORDER,H5_ITER_NATIVE,&n,&attr_simple_search,list);

    return status;
}

/*
 * Routine: herr_t danu_attr_exists(hid_t loc, const char *attr_name, int *exists)
 * Purpose: Test the existence of attribute associated with HDF5 object loc
 * Description: Query the existence of attribute with name attr_name
 *              associated with HDF5 object loc. Routine checks loc and
 *              attr_name as valid input. Will set return value to an error
 *              code if either are invalid. Flag is set to  TRUE if
 *              attribute does exist or FALSE if it does not. Default is
 *              FALSE.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           attr_name          IN              Attribute name
 *           exists             OUT             Flag set to TRUE or FALSE.
 *                                               Default value is FALSE.
 *                              
 * Returns: Returns a negative value if error detected, otherwise returns 0
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. 
 *
*/
herr_t danu_attr_exists(hid_t loc, const char * attr_name, int *exists)
{
  herr_t status = DANU_FAILURE;
  htri_t tri_val;

  *exists = FALSE;
 
  if ( H5_ISA_INVALID_ID(loc) ) {
    DANU_ERROR_MESS("Invalid HDF5 identifier");
    return status;
  }

  if ( DANU_BAD_STRING(attr_name) ) {
    DANU_ERROR_MESS("Invalid string argument");
    return status;
  }

  if ( DANU_BAD_PTR(exists) ) {
    DANU_ERROR_MESS("Invalid pointer argument");
    return status;
  }

  tri_val = H5Aexists(loc,attr_name);
  if ( tri_val > 0 ) {
    *exists = TRUE;
    status = DANU_SUCCESS;
  } else if ( tri_val == 0 ) {
    *exists = FALSE;
    status = DANU_SUCCESS;
  } else {
    DANU_ERROR_MESS("H5Aexists returned a negative value");
  }

  return status;

}

/*
 * Routine: herr_t danu_attr_count(hid_t loc,int *num_found)
 * Purpose: Find the number of attributes attached to HDF5 object loc
 * Description: Iterate over all attributes attached to object
 *              loc and set num_found to the number of attributes found.
 *              Return code indicates if an error was raised. num_found is
 *              initialized to -1 before the iterator is called.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           num_found          OUT             Number of attributes found
 *                              
 * Returns: Returns a negative value if error detected, otherwise returns 0
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. 
 *
*/
herr_t danu_attr_count(hid_t loc, int * num_found)
{
    herr_t status = DANU_FAILURE;
    hsize_t n;
    int count;

    /* Check the input */
    if ( DANU_BAD_PTR(num_found) ) {
        DANU_ERROR_MESS("Invalid pointer");
        return status;
    }
    *num_found = -1;

    if ( H5_ISA_INVALID_ID(loc) ) {
        DANU_ERROR_MESS("Invalid HDF5 identifier");
        return status;
    }

    n = 0;
    count = 0;
    status = H5Aiterate(loc,H5_INDEX_CRT_ORDER,H5_ITER_NATIVE,&n,&attr_simple_count,&count);

    if ( status != 0 ) {
        DANU_ERROR_MESS("Attribute iterator stopped before processing all the attributes");
    }
    else {
        *num_found = count;
    }

    return status;
}

/*
 * Routine: herr_t danu_attr_names(hid_t loc,int num, char **names, int *num_found)
 * Purpose: Find the names of attributes attached to HDF5 object loc
 * Description: Iterate over all attributes attached to object
 *              loc and copy the names found into array names. Will safely copy names to names
 *              and print warning messages if the name is too long for the pointer or if
 *              more attributes were found that could not be copied into names. Calling routines
 *              used call danu_attr_count to determine the number of attributes.
 *              Return code indicates if an error was raised.
 *
 * Parameters:
 *           loc          IN              HDF5 id for location (file,group,dataset)
 *           num          IN              Number of pointers in names array
 *           names        OUT             Array of char pointers to each attribute name found
 *           num_found    OUT             Number of attributes found
 *                              
 * Returns: Returns a negative value if error detected, otherwise returns 0
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Results, such as the names array can not hold all of the attribute
 *         names found or a pointer does not have sufficient memory for a name, are
 *         treated as error conditions. 
 *
*/
herr_t danu_attr_names(hid_t loc, int num, char **names, int *num_found)
{
    herr_t status = DANU_FAILURE;
    danu_list_t *list;
    attr_names_t ndata;
    int i;

    if ( H5_ISA_INVALID_ID(loc) ) {
      DANU_ERROR_MESS("Invalid HDF5 identifier");
      return status;
    }

    if ( num <= 0 ) {
      DANU_ERROR_MESS("Invalid num value (<=0)");
      return status;
    }

    if ( DANU_BAD_PTR(names) ) {
      DANU_ERROR_MESS("Invalid names pointer");
      return status;
    }

    if ( DANU_BAD_PTR(num_found) ) {
      DANU_ERROR_MESS("Invalid num_found pointer");
      return status;
    }
    *num_found = -1;

    /* Find all the attributes  */
    if ( NULL != ( list = danu_list_create() ) ) {
        if ( H5_RETURN_OK(danu_attr_find_all(loc,list)) ) {
          
           *num_found = danu_list_node_count(list);


           /* Now copy the names attached to each node to names */
           ndata.names = names;
           ndata.num   = num;
           ndata.idx = 0;
           if ( DANU_FAILURE != danu_list_iterate(list, &copy_attr_names, &ndata) ) {

	     if ( *num_found > num ) { /* Warn the calling routine */
	       DANU_WARN_MESS("Names array sent to danu_attr_names too small") ; 
	     }
	     else {
	       for(i=*num_found; i<num; i++)
		 names[i] = DANU_CALLOC(char,1);
	       status = DANU_SUCCESS;
	     }

	   }
	   else {
	     DANU_ERROR_MESS("Failed to copy attribute names");
	   }

           danu_list_delete(list);
	}
	else {
	  DANU_ERROR_MESS("Failed to search attribute names");
	}
    }
    else {
        DANU_ERROR_MESS("Failed to allocate linked list");
    }

    return status;

}

/*
 * Routine: char danu_attr_get_type(hid_t loc,char *name)
 * Purpose: Return the data type of attribute name
 * Description: Determine the data type of an existsing attribute named
 *              name attached to HDF5 object loc and return a single
 *              char character to indicate types (int, double, float, string). 
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *                              
 * Returns: Returns a single char. Possible values 'i' (int), 'f' (float), 'd' (double)
 *          'c' (string) or \0 (error)
 * Errors: The input is checked and if an error is detected in the input the routine
 *         will return the NULL character \0. If the attribute does not exist
 *         the routine will return \0.
 *
*/
char danu_attr_get_type(hid_t loc, const char *name)
{
    hid_t attr_id, dtype_id;
    H5T_class_t type;
    size_t type_size;
    char r = '\0';
    int exists;
    herr_t err;

    /* Check input */
    if ( H5_ISA_INVALID_ID(loc) ) {
        DANU_ERROR_MESS("Invalid HDF5 identifier");
        return r;
    }

    if ( DANU_BAD_STRING(name) ) {
        DANU_ERROR_MESS("Invalid name input");
        return r;
    }

    err = danu_attr_exists(loc,name,&exists);
    if (H5_RETURN_OK(err) && exists ) {

        attr_id = H5Aopen(loc, name, H5P_DEFAULT);
        dtype_id = H5Aget_type(attr_id);
        type = H5Tget_class(dtype_id);
        type_size = H5Tget_size(dtype_id);

        if ( type == H5T_INTEGER ) {
	    if ( type_size == sizeof(int) ) {
		r = 'i';
	    }
	    else if ( type_size == sizeof(short) ) {
		r = 's';
	    }
	    else if ( type_size == sizeof(long) ) {
		r = 'l';
	    }
	    else {
		DANU_ERROR_MESS("Unknown integer size");
	    }
        }
        else if ( type == H5T_FLOAT ) {
            if ( type_size == sizeof(float) ) {
                r = 'f';
            }
            else if ( type_size == sizeof(double) ) {
                r = 'd';
            }
            else {
                DANU_ERROR_MESS("Unknown FLOAT class size");
            }
        }
        else if ( type == H5T_STRING ) {
            r = 'c';
        }
        else if ( type == H5T_NO_CLASS ) {
            DANU_ERROR_MESS("Unknown attribute data type");
        }
        else {
            DANU_ERROR_MESS("Undefined error condition");
        }

        H5Tclose(dtype_id);
        H5Aclose(attr_id);

    }
    else {

        if ( H5_RETURN_FAIL(err) ) {
            DANU_ERROR_MESS("Failed to check existence of attribute");
        }
        else {
            DANU_ERROR_MESS("Attribute does not exists");
        }

    }

    return r;
}
/*
 * Routine: size_t danu_attr_get_size(hid_t loc,char *name)
 * Purpose: Return the size of attribute name in bytes
 * Description: Determine the the size of an existsing attribute named
 *              name attached to HDF5 object loc. Size is the number of bytes.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *                              
 * Returns: Returns the size of attribute in number of bytes. 
 *          'c' (string) or \0 (error)
 * Errors: The input is checked and if an error is detected in the input the routine
 *         will return a zero. Will also return zero if any of the internal HDF5
 *         calls return an err condition, otherwise the return value is strictly
 *         positive.
 *
*/
size_t danu_attr_get_size(hid_t loc, const char *name)
{
    size_t size = 0;
    hid_t aid, dtype_id;

    if ( H5_ISA_INVALID_ID(loc) ) {
	DANU_ERROR_MESS("Invalid HDF5 identifier input");
	return size;
    }

    if ( DANU_BAD_PTR(name) ) {
	DANU_ERROR_MESS("Invlaid pointer");
	return size;
    }

    aid = H5Aopen(loc,name,H5P_DEFAULT);
    if ( H5_ISA_VALID_ID(aid) ) {
        dtype_id = H5Aget_type(aid);
	if ( H5_ISA_VALID_ID(dtype_id) ) {
	    size = H5Tget_size(dtype_id);
	    H5Tclose(dtype_id);
	}
	else {
	    DANU_ERROR_MESS("Failed to open data type");
	}
	H5Aclose(aid);
    }
    else {
	DANU_ERROR_MESS("Failed to open attribute");
    }

    return size;
}

    

/*
 * Routine: herr_t danu_attr_write(hid_t loc,char *name, void * value, hid_t type)
 * Purpose: Write a scalar attribute
 * Description: Write a scalar attribute for HDF5 object
 *              assigned to loc identifier. Routine checks that loc is
 *              a valid HDF5 object. Create all the neccessary objects
 *              to write the atrribute to object and close all of the 
 *              objects created in the routine. The return value of the
 *              H5Awrite is passed back to the caller.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           value              IN              Pointer to value
 *           type               IN              Datatype of value
 *                              
 * Returns: Returns the returning value of H5Awrite
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Returns a neagtive value if an error occurs.
 *
*/
herr_t danu_attr_write(hid_t loc, const char * name, void * value, hid_t type)
{
    herr_t status;
    hid_t  attr, dataspace, datatype;
    int attr_exists;
    hsize_t  dims[1] = {1};
    size_t   rank = 1;

    if ( H5_ISA_INVALID_ID(loc) ) {
        DANU_ERROR_MESS("Invalid HDF5 identifier argument");
        return DANU_FAILURE;
    }


    if ( DANU_BAD_STRING(name) ) {
        DANU_ERROR_MESS("Invalid name argument");
        return DANU_FAILURE;
    }

    if ( DANU_BAD_PTR(value) ) {
        DANU_ERROR_MESS("Invalid value pointer");
        return DANU_FAILURE;
    }

    if ( ! H5_ISA_DATATYPE_ID(type) ) {
        DANU_ERROR_MESS("Invalid datatype identifier");
        return DANU_FAILURE;
    }

    status = danu_attr_exists(loc, name, &attr_exists);
    if ( H5_RETURN_FAIL(status) ) {
      DANU_ERROR_MESS("Failed to query existence of attribute");
      return status;
    }


    attr=H5I_INVALID_HID;
    datatype=H5I_INVALID_HID;
    dataspace=H5I_INVALID_HID;
    if ( attr_exists == TRUE ) {
      attr = H5Aopen(loc,name,H5P_DEFAULT);

      if ( H5_ISA_INVALID_ID(attr) ) {
        DANU_ERROR_MESS("Failed to open existing attribute");
      }
    }
    else {
      /* Define the datatype */
      datatype = H5Tcopy(type);
      if ( H5_ISA_INVALID_ID(datatype) ) {
        DANU_ERROR_MESS("Failed to create a copy datatype");
        return DANU_FAILURE;
      }

      /* Set the order if NOT char data type */
      if ( H5Tget_class(datatype) != H5T_STRING ) {
          H5Tset_order(datatype,ATTRIBUTE_ORDER);
      }

      /* Create attribute dataspace */
      dataspace = H5Screate_simple(rank,dims,NULL);

      /* Create the HDF5 attribute object */
      attr = H5Acreate(loc,name,datatype,dataspace,H5P_DEFAULT,H5P_DEFAULT);
    }

    if ( H5_ISA_VALID_ID(attr) ) {
        status = H5Awrite(attr,type,value);
        H5Aclose(attr);
    }
    else {
        DANU_ERROR_MESS("Failed to create HDF5 attribute object");
        status = DANU_FAILURE;
    }

    /* Be tidy and close all the objects created in this routine */
    if ( H5_ISA_VALID_ID(datatype) ) {
      H5Tclose(datatype);
    }

    if ( H5_ISA_VALID_ID(dataspace) ) {
      H5Sclose(dataspace);
    }

    return status;
}
/*
 * Routine: herr_t danu_attr_read(hid_t loc,char *name, void * buffer, hid_t type)
 * Purpose: Read a scalar attribute
 * Description: Read a scalar attribute for HDF5 object
 *              assigned to loc identifier. Routine checks that loc is
 *              a valid HDF5 object. Create all the neccessary objects
 *              to read the atrribute to object and close all of the 
 *              objects created in the routine. The routine assumes that loc
 *              completely describes the HDF5 object and hence the relative
 *              object name defaults to ".". The return value of the
 *              H5Aread is passed back to the caller.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           buffer             OUT             Pointer to buffer attribute is read into
 *           type               IN              Datatype of attribute
 *                              
 * Returns: Returns the returning value of H5Aread
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Returns a neagtive value if an error occurs.
 *
*/
herr_t danu_attr_read(hid_t loc, const char * name, void * buffer, hid_t mem_type)
{
    herr_t status;
    htri_t test;
    hid_t attr;

    /* Check input */
    if ( H5_ISA_INVALID_ID(loc) ) {
        DANU_ERROR_MESS("Invalid HDF5 identifier argument");
        return DANU_FAILURE;
    }

    if ( DANU_BAD_STRING(name) ) {
        DANU_ERROR_MESS("Invalid name argument");
        return DANU_FAILURE;
    }

    if ( DANU_BAD_PTR(buffer) ) {
        DANU_ERROR_MESS("Invalid buffer pointer");
        return DANU_FAILURE;
    }

    if ( ! H5_ISA_DATATYPE_ID(mem_type) ) {
        DANU_ERROR_MESS("Invalid memory datatype");
        return DANU_FAILURE;
    }

    /* Check that the attribute exists */
    test = H5Aexists(loc,name);
    if ( test > 0 ) {
        attr = H5Aopen(loc,name,H5P_DEFAULT);
        DANU_CHECK_H5_RETURN(attr);

        /* Read data inot buffer if the open was successful */
        if ( H5_ISA_VALID_ID(attr) ) {
            status = H5Aread(attr,mem_type,buffer);
            H5Aclose(attr);
        }
        else {
            DANU_ERROR_MESS("Failed to to open attribute");
            status = DANU_FAILURE;
        }
    
    }
    else {
        status = DANU_FAILURE;
        DANU_ERROR_MESS("Attribute does not exist");
    }

    return status;
}

/*
 * Routine: herr_t danu_attr_write_int(hid_t loc,char *name, int value)
 * Purpose: Write a scalar float attribute
 * Description: Write a scalar float attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the 
 *              general write function danu_attr_write. Since that routine
 *              checks the input, no checking is performed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           value              IN              Integer attribute value
 *                              
 * Returns: Returns the returning value danu_attr_write
 * Errors: Errors are possible, but are checked in danu_atrr_write. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_write_int(hid_t loc, const char * name, int value)
{
    return danu_attr_write(loc,name,&value,H5T_NATIVE_INT);
}
/*
 * Routine: herr_t danu_attr_read_int(hid_t loc,char *name, int * buffer)
 * Purpose: Read a scalar float attribute
 * Description: Read a scalar float attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the
 *              general read function danu_attr_read. Since that routine
 *              checks the input, no checking is preformed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           buffer             OUT             Pointer to float buffer attribute is read into
 *                              
 * Returns: Returns the returning value of danu_attr_read
 * Errors: Errors are possible, but are checked in danu_attr_read. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_read_int(hid_t loc, const char * name, int * buffer)
{
    return danu_attr_read(loc,name,buffer,H5T_NATIVE_INT);
}
/*
 * Routine: herr_t danu_attr_write_uint(hid_t loc,char *name, unsigned int value)
 * Purpose: Write a scalar unsigned float attribute
 * Description: Write a scalar unsigned float attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the 
 *              general write function danu_attr_write. Since that routine
 *              checks the input, no checking is performed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           value              IN              Unsigned float attribute value
 *                              
 * Returns: Returns the returning value danu_attr_write
 * Errors: Errors are possible, but are checked in danu_atrr_write. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_write_uint(hid_t loc, const char * name, unsigned int value)
{
    return danu_attr_write(loc,name,&value,H5T_NATIVE_UINT);
}
/*
 * Routine: herr_t danu_attr_read_uint(hid_t loc,char *name, unsigned int * buffer)
 * Purpose: Read a scalar unsigned float attribute
 * Description: Read a scalar unsigned float attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the
 *              general read function danu_attr_read. Since that routine
 *              checks the input, no checking is preformed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           buffer             OUT             Pointer to unsigned float buffer attribute is read into
 *                              
 * Returns: Returns the returning value of danu_attr_read
 * Errors: Errors are possible, but are checked in danu_attr_read. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_read_uint(hid_t loc, const char * name, unsigned int * buffer)
{
    return danu_attr_read(loc,name,buffer,H5T_NATIVE_UINT);
}
/*
 * Routine: herr_t danu_attr_write_float(hid_t loc,char *name, float value)
 * Purpose: Write a scalar float attribute
 * Description: Write a scalar float attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the 
 *              general write function danu_attr_write. Since that routine
 *              checks the input, no checking is performed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           value              IN              Float attribute value
 *                              
 * Returns: Returns the returning value danu_attr_write
 * Errors: Errors are possible, but are checked in danu_atrr_write. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_write_float(hid_t loc, const char * name, float value)
{
    return danu_attr_write(loc,name,&value,H5T_NATIVE_FLOAT);
}
/*
 * Routine: herr_t danu_attr_read_float(hid_t loc,char *name, float * buffer)
 * Purpose: Read a scalar float attribute
 * Description: Read a scalar float attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the
 *              general read function danu_attr_read. Since that routine
 *              checks the input, no checking is preformed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           buffer             OUT             Pointer to float buffer attribute is read into
 *                              
 * Returns: Returns the returning value of danu_attr_read
 * Errors: Errors are possible, but are checked in danu_attr_read. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_read_float(hid_t loc, const char * name, float * buffer)
{
    return danu_attr_read(loc,name,buffer,H5T_NATIVE_FLOAT);
}
/*
 * Routine: herr_t danu_attr_write_double(hid_t loc,char *name, double value)
 * Purpose: Write a scalar double attribute
 * Description: Write a scalar double attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the 
 *              general write function danu_attr_write. Since that routine
 *              checks the input, no checking is performed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           value              IN              Float attribute value
 *                              
 * Returns: Returns the returning value danu_attr_write
 * Errors: Errors are possible, but are checked in danu_atrr_write. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_write_double(hid_t loc, const char * name, double value)
{
    return danu_attr_write(loc,name,&value,H5T_NATIVE_DOUBLE);
}
/*
 * Routine: herr_t danu_attr_read_double(hid_t loc,char *name, double * buffer)
 * Purpose: Read a scalar double attribute
 * Description: Read a scalar double attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the
 *              general read function danu_attr_read. Since that routine
 *              checks the input, no checking is preformed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           buffer             OUT             Pointer to double buffer attribute is read into
 *                              
 * Returns: Returns the returning value of danu_attr_read
 * Errors: Errors are possible, but are checked in danu_attr_read. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_read_double(hid_t loc, const char * name, double * buffer)
{
    return danu_attr_read(loc,name,buffer,H5T_NATIVE_DOUBLE);
}
/*
 * Routine: herr_t danu_attr_write_string(hid_t loc,char *name, char * string)
 * Purpose: Write a scalar char attribute
 * Description: Write a scalar char attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the 
 *              general write function danu_attr_write. Since that routine
 *              checks the input, no checking is performed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           string             IN              Char string attribute value
 *                              
 * Returns: Returns the returning value danu_attr_write
 * Errors: Errors are possible, but are checked in danu_atrr_write. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_write_string(hid_t loc, const char * name, const char * string)
{
    hid_t string_type;
    herr_t status;

    string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type,strlen(string)+1);
    H5Tset_strpad(string_type,H5T_STR_NULLTERM);
    H5Tset_order(string_type,H5T_ORDER_NONE);
    status = danu_attr_write(loc,name,(void*)string,string_type);
    H5Tclose(string_type);
    return status;
}
/*
 * Routine: herr_t danu_attr_read_string(hid_t loc,char *name, char * buffer)
 * Purpose: Read a scalar char attribute
 * Description: Read a scalar char attribute for HDF5 object
 *              assigned to loc identifier. This is a wrapper for the
 *              general read function danu_attr_read. Since that routine
 *              checks the input, no checking is preformed here.
 *
 * Parameters:
 *           loc                IN              HDF5 id for location (file,group,dataset)
 *           name               IN              Name of attribute
 *           buffer             OUT             Pointer to char buffer attribute is read into
 *                              
 * Returns: Returns the returning value of danu_attr_read
 * Errors: Errors are possible, but are checked in danu_attr_read. If an error has occured
 *         a negative is returned.
 *
*/
herr_t danu_attr_read_string(hid_t loc, const char * name, char * buffer, size_t buf_len)
{
    hid_t string_type;
    herr_t status;

    string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type,buf_len);
    H5Tset_strpad(string_type,H5T_STR_NULLTERM);
    status = danu_attr_read(loc,name,(void*)buffer,string_type);
    H5Tclose(string_type);

    return status;
}







