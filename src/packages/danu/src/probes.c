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
 * simulation.c
 *
 *  DANU Simulation Datasets, Groups and Files
 *
 *
 *  Purpose:
 *           This source file defines functions that create simulation datasets, 
 *           groups and files.
 *
 *
 */

#if HAVE_CONFIG_H
#include <danu_config.h>
#endif

#include <hdf5.h>

#include <danu_error.h>
#include <danu_memory.h>
#include <danu_utils.h>
#include <danu_group.h>
#include <danu_dataset.h>
#include <danu_dataset_types.h>
#include <danu_attribute.h>

#include <danu_probes.h>

/* PRIVATE define's */


/* PRIVATE FUNCTIONS */

 
/*
 * Routine: probe_exists(hid_t sid, const char *probename, int *exists)
 * Purpose: Check the existence of probe dataset under the simulation group sid
 * Description: Check the existence of probe dataset probename. The flag exists
 *              indicates if the dataset exists. The return code indicates the status
 *              of the call. Routine will return a negative number if an error occurs.
 *              
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *           exists    OUT             Pointer to flag TRUE or FALSE 
 *                              
 * Returns: A negative return value indicates an error otherwise 0 is returned.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         is not a valid string. Error will also be raised if the routine can not open
 *         the probe subgroup. 
 *         
 *        
 */
herr_t probe_exists(hid_t sid, const char *probename, int *exists)
{
    herr_t status = DANU_FAILURE;

    hid_t pid = probe_open_group(sid);

    if ( DANU_BAD_PTR(exists) ) {
        DANU_ERROR_MESS("Invalid flag pointer");
        return status;
    }

    if ( H5_ISA_VALID_ID(pid) ) {
        status = 0;
        if ( danu_dataset_exists(pid,probename) ) {
            *exists = TRUE;
        }
        else {
            *exists = FALSE;
        }

    }

    return status;
}
/*
 * Routine: probe_data_type(hid_t sid, const char *probename, int *typecode)
 * Purpose: Return the datatype of the probe dataset probename
 * Description: Find the data type of the probe dataset probename. Valid data
 *              types are DANU_DATASET_INT,FLOAT,DOUBLE,STRING. If an error
 *              occurs the typecode is set to DANU_DATASET_UNKNOWN and the
 *              return value is negative.
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *           probename IN              Probe dataset named probename
 *           typecode  OUT             Dataset type code
 *                                     DANU_DATASET_{INT,FLOAT,DOUBLE,STRING} 
 *                              
 * Returns: A negative return value indicates an error otherwise 0 is returned.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         is not a valid string. Error will also be raised if the routine can not open
 *         the probe subgroup. 
 *         
 *        
 */
herr_t probe_data_type(hid_t sid, const char * probename, int *typecode)
{

  herr_t stat = DANU_FAILURE;
  hid_t pid = probe_open_group(sid);

  *typecode=DANU_DATASET_UNKNOWN;
  if ( H5_ISA_VALID_ID(pid) ) {
    stat=danu_dataset_type(pid,probename,typecode);
    danu_group_close(pid);
  }
  else {
    DANU_ERROR_MESS("Failed to open probe dataset group");
  }

  return stat;

}
/*
 * Routine: probe_data_type2(hid_t pid, int *typecode)
 * Purpose: Return the datatype of the probe dataset 
 * Description: Find the data type of the probe dataset identified with id.
 *              Valid datatypes are DANU_DATASET_INT,FLOAT,DOUBLE,STRING. If an error
 *              occurs the typecode is set to DANU_DATASET_UNKNOWN and the
 *              return value is negative.
 *
 * Parameters:
 *           pid       IN              Probe dataset HDF5 identifier
 *           typecode  OUT             Dataset type code
 *                                     DANU_DATASET_{INT,FLOAT,DOUBLE,STRING} 
 *                              
 * Returns: A negative return value indicates an error otherwise 0 is returned.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         is not a valid string. Error will also be raised if the routine can not open
 *         the probe subgroup. 
 *         
 *        
 */
herr_t probe_data_type2(hid_t pid, int *typecode)
{

  herr_t stat = DANU_FAILURE;

  *typecode=DANU_DATASET_UNKNOWN;
  if ( H5_ISA_VALID_ID(pid) ) {
    stat=danu_dataset_type2(pid,typecode);
  }
  else {
    DANU_ERROR_MESS("Failed to open probe dataset group");
  }

  return stat;

}

/*
 * Routine: probe_create_group(hid_t sid)
 * Purpose: Create the probe subgroup under simulation sid.
 * Description: Create the probe subgroup under simulation sid. Routine will close the
 *              HDF5 resource before returning.
 *              
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *                              
 * Returns: A negative return value indicates an error otherwise 0 is returned.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not created.
 *        
 */
herr_t probe_create_group(hid_t sid)
{
    herr_t status = DANU_FAILURE;

    char  probe_grp_name[] = PROBES_GROUP_NAME;
    hid_t pid;

    if ( H5_ISA_INVALID_ID(sid) ) {
        DANU_ERROR_MESS("Invalid simulation HDF5 id");
        return status;
    }

    pid = danu_group_create(sid,probe_grp_name);

    if ( H5_ISA_INVALID_ID(pid) ) {
        DANU_ERROR_MESS("Failed to create the probe subgroup");
    }
    else {
        status = 0;
        danu_group_close(pid);
    }

    return pid;
}
/*
 * Routine: probe_open_group(hid_t sid)
 * Purpose: Open the probe subgroup under simulation sid and return the HDF5 identifier.
 * Description: Open the probe subgroup under simulation sid and return the HDF5
 *              identifier. Input is checked and will return immediatly if
 *              input is not valid.
 *              
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened.
 *        
 */
hid_t probe_open_group(hid_t sid)
{
    char probe_grp_name[] = PROBES_GROUP_NAME;
    hid_t pid = danu_group_open(sid,probe_grp_name);

    return pid;
}

/*
 * Routine: probe_create_data(hid_t sid, const char *probename,
 *                           h5id_t type, int len, int num, void * data)
 *
 * Purpose: Create the probe dataset under the simulation group sid
 *
 * Description: Create the probe dataset under simulation sid. Routine will
 *              return the probe dataset HDF5 identifier 
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *           probename IN              Probe dataset name
 *           type      IN              Dataset type
 *           len       IN              Cannonical length of the 1D buffers
 *           num       IN              Number of 1D array buffers
 *           data      IN              Pointer to data
 *                              
 * Returns: A negative return value indicates an error otherwise valid id is returned.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not created.
 *        
 */
herr_t probe_create_data(hid_t sid, const char *probename,
                         hid_t type, int len, int num, const void *data,
			 hid_t *probe)
{
    herr_t status = DANU_FAILURE;

    char  probe_grp_name[] = PROBES_GROUP_NAME;
    hid_t group,space,filespace;
    hid_t cpl;
    int rank = 2;
    hsize_t hsize[rank], max_size[rank], chunk_size[rank];

    *probe=H5I_INVALID_HID;
    group = danu_group_open(sid,probe_grp_name);

    if ( H5_ISA_INVALID_ID(group) ) {
        DANU_ERROR_MESS("Failed to open the probe subgroup");
    }
    else {
        /* memory space */
	hsize[PROBE_LEN_IDX]=(hsize_t) len;
	hsize[PROBE_NUM_IDX]= (hsize_t) num;
	/* Limits on the dataset */
	max_size[PROBE_LEN_IDX] = (hsize_t) len;
	max_size[PROBE_NUM_IDX] = H5S_UNLIMITED;
	space=H5Screate_simple(rank,hsize,max_size);
	/* Chunk spacing .... how the data is laid out in the file */
	chunk_size[PROBE_LEN_IDX]=(hsize_t) len;
	chunk_size[PROBE_NUM_IDX]= (hsize_t) 1;
	cpl = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk(cpl,rank,chunk_size);
        /* Create the probe dataset */
	*probe=H5Dcreate(group,probename,type,space,
			H5P_DEFAULT,cpl,H5P_DEFAULT);
	/* Write the data */
	filespace=H5Dget_space(*probe);
        status = H5Dwrite(*probe,type,space,filespace,H5P_DEFAULT,data);
	H5Sclose(filespace);
	H5Sclose(space);
	H5Pclose(cpl);
	danu_group_close(group);
    }

    return status; 
}
/*
 * Routine: probe_create_data_int(hid_t sid, const char *probename, int len)
 *
 * Purpose: Create the probe integer dataset under the simulation group sid
 *
 * Description: Create the probe integer dataset under simulation sid. Routine will
 *              return the probe dataset HDF5 identifier. Routine is a wrapper
 *              probe_create_data.
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *           probename IN              Probe dataset name
 *           len       IN              Cannonical length of the 1D buffers
 *           num       IN              Number of 1D arrays
 *           data      IN              Pointer to integer data to write
 *           pid       OUT             Probe data set ID
 *                              
 * Returns: A negative return value indicates an error otherwise valid id is returned.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not created.
 *        
 */

herr_t probe_create_data_int(hid_t sid, const char *probename,
			     int len, int num, const int *data,
			     hid_t *pid)
{
  return probe_create_data(sid,probename,H5T_NATIVE_INT,len,num,data,pid);
}
/*
 * Routine: probe_create_data_floatt(hid_t sid, const char *probename, int len)
 *
 * Purpose: Create the probe float dataset under the simulation group sid
 *
 * Description: Create the probe integer dataset under simulation sid. Routine will
 *              return the probe dataset HDF5 identifier. Routine is a wrapper
 *              probe_create_data.
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *           probename IN              Probe dataset name
 *           len       IN              Cannonical length of the 1D buffers
 *           num       IN              Number of 1D arrays
 *           data      IN              Pointer to float data to write
 *           pid       OUT             Probe data set ID
 *                              
 * Returns: A negative return value indicates an error otherwise valid id is returned.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not created.
 *        
 */

herr_t probe_create_data_float(hid_t sid, const char *probename,
			     int len, int num, const float *data,
			     hid_t *pid)
{
  return probe_create_data(sid,probename,H5T_NATIVE_FLOAT,len,num,data,pid);
}
/*
 * Routine: probe_create_data_doublet(hid_t sid, const char *probename, int len)
 *
 * Purpose: Create the probe double dataset under the simulation group sid
 *
 * Description: Create the probe integer dataset under simulation sid. Routine will
 *              return the probe dataset HDF5 identifier. Routine is a wrapper
 *              probe_create_data.
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *           probename IN              Probe dataset name
 *           len       IN              Cannonical length of the 1D buffers
 *           num       IN              Number of 1D arrays
 *           data      IN              Pointer to double data to write
 *           pid       OUT             Probe data set ID
 *                              
 * Returns: A negative return value indicates an error otherwise valid id is returned.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not created.
 *        
 */

herr_t probe_create_data_double(hid_t sid, const char *probename,
			        int len, int num, const double *data,
			        hid_t *pid)
{
  return probe_create_data(sid,probename,H5T_NATIVE_DOUBLE,len,num,data,pid);
}


/*
 * Routine: probe_open_data(hid_t sid)
 * Purpose: Open the probe dataset under simulation sid and return the HDF5 identifier.
 * Description: Open the probe dataset under simulation sid and return the HDF5
 *              identifier. Input is checked and will return immediately if
 *              input is not valid.
 *              
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *           probename IN              Probe dataset name
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened.
 *        
 */
hid_t probe_open_data(hid_t sid, const char * probename)
{
    char probe_grp_name[] = PROBES_GROUP_NAME;
    hid_t pid = danu_group_open(sid,probe_grp_name);
    hid_t data = H5I_INVALID_HID;

    if ( H5_ISA_VALID_ID(pid) ) {
        data = danu_dataset_open(pid,probename);
        danu_group_close(pid);
    }

    return data;
}

/*
 * Routine: probe_count(hid_t sid,int * num)
 * Purpose: Find the number of probe datasets under simulation sid
 * Description: Open the probe subgroup under simulation sid and find
 *              the number of datasets in the probe subgroup. num is initialized
 *              to -1 and the return code should be checked for errors.
 *              Routine assumes that ALL links under the probe
 *              subgroup are datasets.
 *              
 *
 * Parameters:
 *           sid       IN              HDF5 identifier for  Simulation group
 *           num       OUT             Number of datasets under probes group
 *                              
 * Returns: A negative return value indicates an error, other wise return code is 0.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or an error occurs while traversing the probe
 *         subgroup.
 *        
 */
herr_t probe_count(hid_t sid, int * num)
{
    herr_t status = DANU_FAILURE;

    hid_t pid = probe_open_group(sid);
    hsize_t nlinks;

    /* Check the input */
    if ( DANU_BAD_PTR(num) ) {
        DANU_ERROR_MESS("Invalid pointer");
        return status;
    }

    if ( H5_ISA_INVALID_ID(pid) ) {
        DANU_ERROR_MESS("Failed to open probes group");
        return status;
    }

    status = danu_group_get_nlinks(pid,&nlinks);
    *num = -1;
    if ( H5_RETURN_OK(status) ) {
        *num = (int) nlinks;
    }

    danu_group_close(pid);

    return status;
}


/*
 * Routine: probe_list(hid_t sid,int num, char **probenames)
 * Purpose: Return the probe dataset names under simulation sid. 
 * Description: Open the probe subgroup under simulation sid and return the dataset
 *              probe names. The input variables num and size, define the number
 *              and size of each memory pointed to in probenames. The number of
 *              probe datasets can be found by calling probe_count. This routine
 *              allocates memory for each pointer in probenames. The calling
 *              routine is responsible for freeing this memory. 
 *              A negative return code indicates an error.
 *              
 *              
 *
 * Parameters:
 *           sid        IN              HDF5 identifier for  Simulation group
 *           num        IN              Number of pointers in probenames
 *           probenames OUT             Array holding char pointers to probe 
 *                                       dataset names                           
 *           num_found  OUT             Number of probe datasets found
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or if an error is encountered while traversing
 *         the probe datasets. Finding 
 *        
 */
herr_t probe_list(hid_t sid, int num, char **probenames, int *num_found)
{
    herr_t status = DANU_FAILURE;

    hid_t pid = probe_open_group(sid);

    /* Check the input */
    if ( H5_ISA_INVALID_ID(pid) ) {
        DANU_ERROR_MESS("Failed to open the probe subgroup");
        return status;
    }

    status = danu_group_get_datasets(pid,num,probenames,num_found);

    if ( H5_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to retrieve probe dataset names");
    }
    else {
        if ( *num_found > num ) {
            DANU_WARN_MESS("Missing probe datasets, char pointer array was not large enough");
            status = -1;
        }
    }

    danu_group_close(pid);

    return status;
}

/*
 * Routine: probe_attribute_list(hid_t sid,
 *                               const char *probename, 
 *                               int num, char **attrnames, int *num_found)
 * Purpose: Return the attribute names associated with probename under simulation sid.
 * Description: Open the probename dataset under the simulation sid subgroup Probes.
 *              The attribute names are copied to the attrnames. This copy will
 *              not overwrite the pointer bounds.  
 *              
 *
 * Parameters:
 *           sid        IN              HDF5 identifier for  Simulation group
 *           probename  IN              Probe dataset name
 *           num        IN              Number of pointers in probenames
 *           probenames OUT             Array holding char pointers to probe 
 *                                       dataset names                           
 *           num_found  OUT             Number of attribute names found
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or if an error is encountered while traversing
 *         the probe attributes. 
 *        
 */
herr_t probe_attribute_list(hid_t sid, 
	                    const char * probename,
			    int num,
			    char **attrnames,
			    int *num_found)
{
    herr_t status = DANU_FAILURE;

    hid_t pid = probe_open_group(sid);
    hid_t data;
    int i;

    if ( H5_ISA_INVALID_ID(pid) ) {
        DANU_ERROR_MESS("Failed to open probe group");
        return status;
    }

    data = danu_dataset_open(pid,probename);

    if ( H5_ISA_VALID_ID(data) ) {
       status = danu_attr_names(data,num,attrnames,num_found);
       danu_dataset_close(data);

       if ( *num_found < num ) {
	   for(i=*num_found;i<num;i++) {
	       attrnames[i] = DANU_MALLOC(char,1);
	       *attrnames[i] = '\0';
	   }
       }
    }

    danu_group_close(pid);

    return status;
}
/*
 * Routine: probe_attribute_list(hid_t pid,
 *                               int num, char **attrnames, int *num_found)
 * Purpose: Return the attribute names associated with pid.
 * Description: Read the attribute names associated with pid.
 *              The attribute names are copied to the attrnames. This copy will
 *              not overwrite the pointer bounds.  
 *              
 *
 * Parameters:
 *           pid        IN              Probe dataset HDF5 identifier 
 *           num        IN              Number of pointers in probenames
 *           probenames OUT             Array holding char pointers to probe 
 *                                       dataset names                           
 *           num_found  OUT             Number of attribute names found
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or if an error is encountered while traversing
 *         the probe attributes. 
 *        
 */
herr_t probe_attribute_list2(hid_t pid, 
			    int num,
			    char **attrnames,
			    int *num_found)
{
    herr_t status = DANU_FAILURE;
    int i;

    if ( H5_ISA_INVALID_ID(pid) ) {
        DANU_ERROR_MESS("Invalid HDf5 identifer");
        return status;
    }

    status = danu_attr_names(pid,num,attrnames,num_found);

    if ( *num_found < num ) {
      for(i=*num_found;i<num;i++) {
        attrnames[i] = DANU_MALLOC(char,1);
	*attrnames[i] = '\0';
      }
    }

    return status;
}

/*
 * Routine: probe_data_num(hid_t sid, const char *pname, int *num)
 *
 * Purpose: Return the number of probe data arrays.
 * Description: Probe datasets are num x N datasets that allow appending
 *              1D N-length datasets. This routine returns the 
 *              number of 1D arrays.
 *
 * Parameters:
 *           sid        IN              HDF5 identifier for  Simulation group
 *           probename  IN              Probe dataset name
 *           num        OUT             Number of 1D probe data arrays
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or if an error is encountered while reading
 *         the dataset dimension.
 *        
 */
herr_t probe_data_num(hid_t sid, const char *pname, int *num)
{
  herr_t status = DANU_FAILURE;
  hid_t group;
  hsize_t size[2];
  int rank = 2;

  /* Check input */
  if ( ! H5_ISA_GROUP_ID(sid) ) {
    DANU_ERROR_MESS("Invalid simulation group id");
    return status;
  }

  if ( DANU_BAD_STRING(pname) ) {
    DANU_ERROR_MESS("Invalid probe name");
    return status;
  }

  if ( DANU_BAD_PTR(num) ) {
    DANU_ERROR_MESS("Invalid length pointer");
    return status;
  }

  size[PROBE_LEN_IDX]=-1;
  size[PROBE_NUM_IDX]=-1;
  group=probe_open_group(sid);
  if ( H5_ISA_VALID_ID(group) ) {
    status = danu_dataset_dimensions(group,pname,rank,size);
    danu_group_close(group);
  }
  *num=size[PROBE_NUM_IDX];

  return status;
}
/*
 * Routine: probe_data_num2(hid_t pid, int *num)
 *
 * Purpose: Return the number of probe data arrays.
 * Description: Probe datasets are num x N datasets that allow appending
 *              1D N-length datasets. This routine returns the 
 *              number of 1D arrays.
 *
 * Parameters:
 *           pid        IN              Probe dataset HDF5 identifier
 *           num        OUT             Number of 1D probe data arrays
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or if an error is encountered while reading
 *         the dataset dimension.
 *        
 */
herr_t probe_data_num2(hid_t pid, int *num)
{
  herr_t status = DANU_FAILURE;
  hid_t group;
  hsize_t size[2];
  int rank = 2;

  /* Check input */
  if ( H5_ISA_INVALID_ID(pid) ) {
    DANU_ERROR_MESS("Invalid HDF5 identifier");
    return status;
  }

  if ( DANU_BAD_PTR(num) ) {
    DANU_ERROR_MESS("Invalid length pointer");
    return status;
  }

  size[PROBE_LEN_IDX]=-1;
  size[PROBE_NUM_IDX]=-1;
  status = danu_dataset_dimensions2(pid,rank,size);
  *num=size[PROBE_NUM_IDX];

  return status;
}

/*
 * Routine: probe_data_length(hid_t sid, const char *pname, int *len)
 *
 * Purpose: Return the canonical probe length.
 * Description: Probe datasets are num x N datasets that allow appending
 *              1D N-length datasets. This routine returns the canonical 
 *              length N.
 *
 * Parameters:
 *           sid        IN              HDF5 identifier for  Simulation group
 *           probename  IN              Probe dataset name
 *           len        OUT             Allowed 1D array size
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or if an error is encountered while reading
 *         the dataset dimension.
 *        
 */
herr_t probe_data_length(hid_t sid, const char *pname, int *len)
{
  herr_t status = DANU_FAILURE;
  hid_t group;
  hsize_t size[2];
  int rank = 2;

  /* Check input */
  if ( ! H5_ISA_GROUP_ID(sid) ) {
    DANU_ERROR_MESS("Invalid simulation group id");
    return status;
  }

  if ( DANU_BAD_PTR(len) ) {
    DANU_ERROR_MESS("Invalid length pointer");
    return status;
  }

  /* Get the probe dimensions */
  group = probe_open_group(sid);
  size[PROBE_LEN_IDX]=-1;
  size[PROBE_NUM_IDX]=-1;
  if ( H5_ISA_VALID_ID(group) ) {
    status = danu_dataset_dimensions(group,pname,rank,size);
    danu_group_close(group);
  }
  *len=size[PROBE_LEN_IDX];

  return status;
}

/*
 * Routine: probe_data_length2(hid_t pid, int *len)
 *
 * Purpose: Return the canonical probe length.
 * Description: Probe datasets are num x N datasets that allow appending
 *              1D N-length datasets. This routine returns the canonical 
 *              length N.
 *
 * Parameters:
 *           pid        IN              Probe dataset HDF5 identifier
 *           len        OUT             Allowed 1D array size
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or if an error is encountered while reading
 *         the dataset dimension.
 *        
 */
herr_t probe_data_length2(hid_t pid, int *len)
{
  herr_t status = DANU_FAILURE;
  hsize_t size[2];
  int rank = 2;

  /* Check input */
  if ( H5_ISA_INVALID_ID(pid) ) {
    DANU_ERROR_MESS("Invalid HDF5 identifier");
    return status;
  }

  if ( DANU_BAD_PTR(len) ) {
    DANU_ERROR_MESS("Invalid length pointer");
    return status;
  }

  /* Get the probe dimensions */
  size[PROBE_LEN_IDX]=-1;
  size[PROBE_NUM_IDX]=-1;
  status = danu_dataset_dimensions2(pid,rank,size);
  *len=size[PROBE_LEN_IDX];

  return status;
}

/*
 * Routine: probe_data_read(hid_t pid,
 *                              hid_t type, 
 *                              void * data)
 * Purpose: Read data contained in probe dataset associated with pid.
 * Description: This is the main probe data read routine.
 *              The specific data type read calls are wrappers for
 *              this routine. The input is checked if accessed in this routine, otherwise
 *              the routines called in this routine will check the input. If an error
 *              is encountered the return value will be negative. 
 *              
 *
 * Parameters:
 *           pid        IN              Porbe dataset HDF5 identifier
 *           type       IN              Data type (HDF5 type)
 *           data       OUT             Array holding probe data
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or if an error is encountered while traversing
 *         the probe attributes. 
 *        
 */
herr_t probe_data_read(hid_t pid, hid_t type, void *data)
{
    herr_t status = DANU_FAILURE;

    int len, num;
    int rank = 2;
    hsize_t hsize[rank];

    /* Check the input values */
    if ( H5_ISA_INVALID_ID(pid) ) {
      DANU_ERROR_MESS("Invalid HDF5 identifier");
      return status;
    }

    if ( H5Tget_class(type) == H5T_NO_CLASS ) {
      DANU_ERROR_MESS("Invalid HDF5 dataset type");
      return status;
    }

    if ( DANU_BAD_PTR(data) ) {
      DANU_ERROR_MESS("Invalid data pointer");
      return status;
    }

    if ( DANU_RETURN_FAIL(probe_data_length2(pid,&len)) ||
	 DANU_RETURN_FAIL(probe_data_num2(pid,&num)) ) {
      DANU_ERROR_MESS("Failed to stat the size of the probe dataset");
      return status;
    }


    hsize[PROBE_LEN_IDX] = (hsize_t) len;
    hsize[PROBE_NUM_IDX] = (hsize_t) num;
    status = danu_dataset_read(pid,NULL,type,rank,hsize,data);

    return status;
}
/*
 * Routine: probe_data_read_int(hid_t pid,
 *                              int * data)
 * Purpose: Read integer data contained in dataset pid.
 * Description: This is a wrapper function for probe_data_read for integer datasets.
 *
 * Parameters:
 *           sid        IN              Probe daataset HDF5 identifier
 *           data       OUT             Integer data array
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         dataset does not exist. 
 *        
 */
herr_t probe_data_read_int(hid_t pid, int *data)
{
    return probe_data_read(pid,H5T_NATIVE_INT,data);
}
/*
 * Routine: probe_data_read_double(hid_t pid,
 *                              double * data)
 * Purpose: Read double data contained in dataset pid
 * Description: This is a wrapper function for probe_data_read for double datasets.
 *
 * Parameters:
 *           pid        IN              Probe dataset HDF5 identifier
 *           data       OUT             Double data array
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         dataset does not exist. 
 *        
 */
herr_t probe_data_read_double(hid_t pid, double *data)
{
    return probe_data_read(pid,H5T_NATIVE_DOUBLE,data);
}

/*
 * Routine: probe_data_read_float(hid_t pid,
 *                                float * data)
 * Purpose: Read float data contained in dataset pid
 * Description: This is a wrapper function for probe_data_read for float datasets.
 *
 * Parameters:
 *           pid        IN              Probe dataset HDF5 identifier
 *           data       OUT             Float data array
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         dataset does not exist. 
 *        
 */
herr_t probe_data_read_float(hid_t pid, float *data)
{
    return probe_data_read(pid,H5T_NATIVE_FLOAT,data);
}
/*
 * Routine: probe_data_write(hid_t pid,
 *                           hid_t type,
 *                           int num
 *                           void * data)
 * Purpose: Write data contained in data to probename dataset pid 
 * Description: WRite data to probe dataset associated with ppid. Routine assumes
 *              the data is shaped num x N where N is the canonical 1D array
 *              length. This is the main probe write routine. All other
 *              routines are wrappers to this function. 
 *
 * Parameters:
 *           pid        IN              Probe dataset HDF5 identifier
 *           type       IN              Data type (HDF5 type)
 *           num        IN              Number of 1D arrays
 *           data       IN              Array holding data to write
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         subgroup is not opened, or if an error is encountered while traversing
 *         the probe attributes. 
 *        
 */
herr_t probe_data_write(hid_t pid, hid_t type, int num, const void *data)
{
    herr_t status = DANU_FAILURE;
    int i, exists, cur_num, len;
    hid_t filespace, memspace, slab;
    int rank = 2;
    hsize_t hsize[rank], new_size[rank], offset[rank];

    /* Check the input values */
    if ( H5_ISA_INVALID_ID(pid) ) {
      DANU_ERROR_MESS("Invalid HDF5 identifier");
      return status;
    }

    if ( H5Tget_class(type) == H5T_NO_CLASS ) {
      DANU_ERROR_MESS("Invalid HDF5 dataset type");
      return status;
    }

    if ( num <= 0 ) {
      DANU_ERROR_MESS("Invalid probe number value");
      return status;

    }

    if ( DANU_BAD_PTR(data) ) {
      DANU_ERROR_MESS("Invalid data pointer");
      return status;
    }

    if ( DANU_RETURN_FAIL(probe_data_length2(pid,&len)) ) {
      DANU_ERROR_MESS("Failed to determine canonical length");
      return status;
    }

    if ( DANU_RETURN_FAIL(probe_data_num2(pid,&cur_num)) ) {
      DANU_ERROR_MESS("Failed to determine probe data number");
      return status;
    }

  
    /* Size array and HDF5 space associated with the data */
    hsize[PROBE_LEN_IDX] = (hsize_t) len;
    hsize[PROBE_NUM_IDX] = (hsize_t) num;
    memspace = H5Screate_simple(rank,hsize,NULL);

    offset[PROBE_LEN_IDX] = (hsize_t) 0;
    offset[PROBE_NUM_IDX] = (hsize_t) cur_num; /* Zero-based offset */
    new_size[PROBE_LEN_IDX] = (hsize_t) len;
    new_size[PROBE_NUM_IDX] = (hsize_t) (num + cur_num);
    if ( DANU_RETURN_FAIL(H5Dextend(pid,new_size) ) ) {
        DANU_ERROR_MESS("Failed to extend dataset");
    }


    filespace = H5Dget_space(pid);
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset,
                                 NULL, hsize, NULL);
    /* Write the data */
    status = H5Dwrite(pid,type,memspace,filespace,H5P_DEFAULT,data);
    
    /* Close all the dataspaces */
    H5Sclose(filespace);
    H5Sclose(memspace);

    return status;
}
/*
 * Routine: probe_data_write_int(hid_t pid,
 *                              int dim, const int * size, int * data)
 * Purpose: Write integer data to dataset probename pid
 * Description: This is a wrapper function for probe_data_write for integer datasets.
 *
 * Parameters:
 *           pid        IN              Probe dataset HDF5 identifier
 *           dim        IN              Dimension of the dataset
 *           size       IN              1D array containing the size of each dimension
 *           data       OUT             Integer data array
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         dataset does not exist. 
 *        
 */
herr_t probe_data_write_int(hid_t pid, int num, const int *data)
{
    return probe_data_write(pid,H5T_NATIVE_INT,num,data);
}
/*
 * Routine: probe_data_write_double(hid_t pid,
 *                                  int num, double * data)
 * Purpose: Write double data to probe dataset associated with pid.
 * Description: This is a wrapper function for probe_data_write for double datasets.
 *
 * Parameters:
 *           pid        IN              HDF5 identifier for  Simulation group
 *           num        IN              Number of 1D n-length arrays 
 *           data       OUT             Double data array
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         dataset does not exist. 
 *        
 */
herr_t probe_data_write_double(hid_t pid, int num, const double *data)
{
    return probe_data_write(pid,H5T_NATIVE_DOUBLE,num,data);
}
/*
 * Routine: probe_data_write_float(hid_t pid,
 *                                 int num, const float * data)
 * Purpose: Write float data to dataset probename associatesd with pid
 * Description: This is a wrapper function for probe_data_write for float datasets.
 *
 * Parameters:
 *           pid        IN              Probe dataset HDF5 identifier
 *           num        IN              Number of 1D n-length arrays
 *           data       OUT             Float data array
 *                              
 * Returns: A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the probe
 *         dataset does not exist. 
 *        
 */
herr_t probe_data_write_float(hid_t pid, int num, const float *data)
{
    return probe_data_write(pid,H5T_NATIVE_FLOAT,num,data);
}
