/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * series_group.c
 *
 *  TOUT time series groups
 *
 *
 *  Purpose:
 *           This source file defines functions that create time series groups. 
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
#include <danu_group.h>
#include <danu_attribute.h>
#include <danu_sim.h>
#include <danu_series-group.h>


 /* Private functions */

/*
 * Routine: sequence_create_root_group(hid_t)
 * Purpose: Create the base group Series Data under simulation sid
 * Description: Create the base group Series Data under simulation sid. The group 
 *              is closed before the routine returns. Error will occur if sid
 *              is not a valid HDF5 identifier or if the Series Data group is
 *              not created.
 *
 * Parameters:
 *           sid                IN              HDF5 identifier for the Simulation group
 *                              
 * Returns: A negative value is returned if an error occurs.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the routine fails
 *         to create the group. 
 *        
 */
herr_t sequence_create_root_group(hid_t sid)
{
    herr_t status = DANU_FAILURE;

    const char seq_root_name[] = SERIES_GROUP_NAME;
    hid_t gid = danu_group_create(sid,seq_root_name);

    if ( H5_ISA_VALID_ID(gid) ) {
        status = 0;
        danu_group_close(gid);
    }
    else {
        DANU_ERROR_MESS("Failed to create the series data group");
    }

    return status;
}

/*
 * Routine: sequence_open_root_group(hid_t)
 * Purpose: Open the base group Series Data under simulation sid
 * Description: Open the base group Series Data under simulation sid.
 *              Error will occur if sid is not a valid HDF5 identifier or
 *              if the Series Data group does not exist.
 *
 * Parameters:
 *           sid                IN              HDF5 identifier for the Simulation group
 *                              
 * Returns: A negative value is returned if an error occurs.
 *          
 * Errors: An error is raised if sid is not a valid HDF5 identifier or the routine fails
 *         to open the group.
 *        
 */
hid_t sequence_open_root_group(hid_t sid)
{
    const char seq_root_name[] = SERIES_GROUP_NAME;
    hid_t gid = danu_group_open(sid,seq_root_name);

    if ( H5_ISA_INVALID_ID(gid) ) {
        DANU_ERROR_MESS("Failed to open the series data group");
    }

    return gid;
}



/*
 * Routine: sequence_exists(hid_t sid, const char *seriesname, int *exists)  
 * Purpose: Checks the existence of seriesname subgroup under the Simulation group sid 
 * Description: Checks the seriesname against the list of series subgroups under
 *              a Simualtion Series group. The Simulation handle is sid. 
 *              The 'exists' is initall yet to FALSE. If the series subgroup is found
 *              it is set to TRUE. While building the series name, the routine checks
 *              the input to see if the series already includes the Simulation group name.
 *              It will pre-append this name if it does not.
 *               
 *
 * Parameters:
 *           sid                IN              HDF5 identifier for Simulation group
 *           seriesname         IN              Series subgroup name to search
 *           *exists            OUT             Flag set to TRUE or FALSE
 *                              
 * Returns: Returns a negative value if an error occurs.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routien will not
 *         accept a NULL pointer for seriesname. Other errors may be raised while 
 *         retrieving the series subgroup names.
 *        
 */
herr_t sequence_exists(hid_t sid, const char *seriesname, int *exists)
{
    herr_t status = DANU_FAILURE;

    hid_t gid = sequence_open_root_group(sid);

    if ( DANU_BAD_PTR(exists) ) {
        DANU_ERROR_MESS("Invalid pointer");
        return status;
    }
    *exists = FALSE;

    if ( H5_ISA_VALID_ID(gid) ) {
        if ( danu_group_exists(gid,seriesname) ) {
            *exists = TRUE;
        }
        status = 0;
    }

    return status;
}

/*
 * Routine: sequence_count(hid_t sid, int *nseries)  
 * Purpose: Find the number of series subgroups under the simulation group sid
 * Description: Returns the number of links found under the simulation group
 *              identified. Routine assumes all links found under the 
 *              series subgroup are series time-step (sequence) groups. Calling 
 *              routines should check the return status for errors, NOT nseries.
 *
 * Parameters:
 *           sid                IN              HDF5 identifier for Simulation group
 *           *nseries           OUT             Number of series subgroups found 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected.
 *        
 */
herr_t sequence_count(hid_t sid, int *nseries)
{
    herr_t status = DANU_FAILURE;

    hid_t gid = sequence_open_root_group(sid); 
    hsize_t nlinks;

    if (DANU_BAD_PTR(nseries)){
        DANU_ERROR_MESS("Invalid pointer");
        return status;
    }
    *nseries = 0;

    if ( H5_RETURN_OK(danu_group_get_nlinks(gid,&nlinks) ) ) {
        *nseries = (int) nlinks;
        status = 0;
    }
    else {
        DANU_ERROR_MESS("Failed to count links in Series Data group");
    }

    return status;
}
/*
 * Routine: sequence_name_max_bytes(hid_t sid, size_t *max_bytes)  
 * Purpose: Determine the max number of bytes that will hold a sequence group name
 * Description: Routine queries the Simulation group for the number of sequence
 *              subgroups. It then uses this number to build the last sequence
 *              group name. Since the format of these groups is
 *              "SEQUENCE_BASE_NAME nnn"
 *              the group with the largest number will also have the longest name.
 *              Note the number of bytes returned INCLUDES the '\0' character.
 *
 * Parameters:
 *           sid                IN              HDF5 identifier for the Simulation group
 *           *nseries           OUT             Number of series subgroups found 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected.
 *        
 */

herr_t sequence_name_max_bytes(hid_t sid, size_t*max_bytes)
{
    herr_t status = DANU_FAILURE;

    int ngrps;

    char  base_grp_name[] = SEQUENCE_BASE_NAME;

    if ( DANU_BAD_PTR(max_bytes) ) {
        DANU_ERROR_MESS("Invalid pointer");
        return status;
    }
    *max_bytes = 0;

    if (H5_RETURN_OK(sequence_count(sid,&ngrps) ) ) {
        *max_bytes = strlen(base_grp_name);
        (*max_bytes)++;
        if ( ngrps ) {
            *max_bytes+=2;
            *max_bytes+= (int) log10( (double) ngrps);
        }
        status = 0;
    }


    return status;
}
/*
 * Routine: sequence_list(hid_t sid, int num, const size_t*size, char **datanames)  
 * Purpose: Return the list of series sub-group names in datanames
 * Description: Returns the list of series sub group names found under the Simulation
 *              subgroup SERIES_GROUP_NAME. Calling routine must provide the number
 *              (num) of names datanames can hold. This routine will allocate memory
 *              for each pointer in datanames. The calling routine is responsible
 *              for freeing this memory.
 *              
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for Simulation group
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
herr_t sequence_list(hid_t sid, int num, char **datanames, int *num_found)
{
    herr_t status = DANU_FAILURE;
    hid_t gid = sequence_open_root_group(sid);

    if ( H5_ISA_VALID_ID(gid) ) {
        status = danu_group_get_subgroups(gid,num,datanames,num_found);
        danu_group_close(gid);
        status = 0;
    }

    return status;
}
/*
 * Routine: sequence_get_handle(hid_t sid, const char *seriesname, hid_t * nsid)  
 * Purpose: Return the HDF5 handle for sequence subgroup seriesname.
 * Description: Returns the handle of a sequence series subgroup named seriesname
 *              Function will append the series data group name if it is not
 *              included in seriesname. Calling routine should check the return
 *              code for errors. In teh event an error is detected id
 *              will be set to an invalid HDF5 handle.
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for Simulation group
 *          *seriesname  IN      Sequence series name  
 *          *nsid        OUT     HDF5 identifier for the Sequence series group 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected. The nsid handle will be set to
 *         an invalid HDF5 identifier if an error occurs.
 *        
*/
herr_t sequence_get_handle(hid_t sid, const char *seriesname, hid_t * nsid)  
{
    herr_t status = DANU_FAILURE;
    hid_t gid = sequence_open_root_group(sid);

    if ( DANU_BAD_PTR(nsid) ) {
        DANU_ERROR_MESS("Invalid pointer");
        return status;
    }
    *nsid = H5I_BADID;

    if ( H5_ISA_VALID_ID(gid) ) {
        *nsid = danu_group_open(gid,seriesname);
        status = ( H5_ISA_VALID_ID(*nsid) ) ? 0 : DANU_FAILURE;
        danu_group_close(gid);
    }
    else {
        DANU_ERROR_MESS("Failed to open the root Series Data group");
    }

    return status;
}
/*
 * Routine: sequence_get_handle_byid(hid_t sid, int id, hid_t * nsid)  
 * Purpose: Return the HDF5 handle for sequence subgroup by series id.
 * Description: Returns the handle of a sequence series subgroup by series id.
 *              Function builds the series name from the id and then calls
 *              sequence_get_handle with this name.
 *              Calling routine should check the return
 *              code for errors. In the event an error is detected id
 *              will be set to an invalid HDF5 handle.
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for Simulation group
 *           id          IN      Sequence series id 
 *          *nsid        OUT     HDF5 identifier for the Sequence series group 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected. The nsid handle will be set to
 *         an invalid HDF5 identifier if an error occurs.
 *        
*/
herr_t sequence_get_handle_byid(hid_t sid, int id, hid_t *nsid)
{

  herr_t stat = DANU_FAILURE;
  char *name;

  *nsid=H5I_INVALID_HID;
  if ( id < 0 ) {
    DANU_ERROR_MESS("Invalid series group identifier");
    return stat;
  }
  name=sequence_get_name(id);
  stat=sequence_get_handle(sid,name,nsid);
  DANU_FREE(name);

  return stat;

}
/*
 * Routine: sequence_getNextID(hid_t sid, int cycle, double time, hid_t * nsid)  
 * Purpose: Return the identifier for the next sequence series group 
 * Description: Returns the handle of the next sequence series group. The sequence
 *              number and time are stored as attributes in this group.
 *
 * Parameters:
 *           sid         IN      HDF5 identifier for Simulation group
 *           cycle       IN      Cycle number of the timestep
 *           time        IN      Time value for this sequence  
 *          *nsid        OUT     HDF5 identifier for the sequence series group 
 *                              
 * Returns: Returns a negative value if an error occurs, zero (0) otherwise.
 *         
 * Errors: Input is checked for valid HDF5 handle and pointers. Routine returns a
 *         negative value if an error is detected. The nsid handle will be set to
 *         an invalid HDF5 identifier if an error occurs.
 *        
*/
herr_t sequence_getNextID(hid_t sid, int cycle, double time, hid_t *nsid)
{
    herr_t status = DANU_FAILURE;

    hid_t gid = sequence_open_root_group(sid);

    char time_attr_name[]   = SEQ_TIME_ATTR_NAME;
    char cycle_attr_name[]  = SEQ_CYCNUM_ATTR_NAME;
    char seqnum_attr_name[] = SEQ_SEQNUM_ATTR_NAME; 

    int nseries;
    char * new_series;

    if ( DANU_BAD_PTR(nsid) ) {
        DANU_ERROR_MESS("Invalid pointer");
        return status;
    }
    *nsid = -1;

    if ( H5_ISA_VALID_ID(gid) ) {
        sequence_count(sid,&nseries);
        nseries++;
	new_series = sequence_get_name(nseries);
        if ( new_series != NULL ) {
            *nsid = danu_group_create(gid,new_series);
            if ( H5_ISA_VALID_ID(*nsid) ) {
                danu_attr_write_double(*nsid,time_attr_name,time);
                danu_attr_write_int(*nsid,seqnum_attr_name,nseries);
                danu_attr_write_int(*nsid,cycle_attr_name,cycle);
                status = 0;
            }
            else {
                DANU_ERROR_MESS("Failed to create the next series group");
            }

            DANU_FREE(new_series);
        }
        else {
            DANU_ERROR_MESS("Failed to allocate new_series name...memory exhausted");
        }
        danu_group_close(gid);
    }
    else {
        DANU_ERROR_MESS("Failed to open the root Series Data group");
    }

    return status;

}
/*
 * Routine: sequence_get_name(num)
 * Purpose: Return pointer to string containing the series name for series num
 * Description: Given a series number, return the series group associated with
 *              that number. Calling routine is responsible for freeing the 
 *              pointer.
 *
 * Parameters:
 *           num         IN      Integer series number
 *          *name        OUT     Pointer to series name
 *                              
 * Returns: Returns a NULL pointer if error occurs.
 *         
 * Errors: Input is checked for a non-negative integer type. Routine does not
 *         check the returning pointer.
 *        
*/
char * sequence_get_name(int num)
{
  char *ptr;
  size_t bytes;
  char base_series_name[] = SEQUENCE_BASE_NAME;

  ptr = NULL;
  if ( num > 0 ) {
    bytes = strlen(base_series_name);
    bytes++;
    bytes+=  ( (size_t) ( log10((double)num) ) + 1);
    bytes++;
    ptr = DANU_MALLOC(char, bytes);
    if ( ptr != NULL)
      sprintf(ptr,"%s %d",base_series_name,num);
  } else {
    DANU_ERROR_MESS("Invalid sequence number");
  }
  

  return ptr;
}
/*
 * Routine: sequence_get_time(nsid, double *time)
 * Purpose: Return the time attribute associated with sequence group nsid
 * Description: Given a identifier for an existing sequence group, return the 
 *              time attribute associated with this groups. Returns a
 *              non-negative value if an error is detected. 
 *
 * Parameters:
 *           nsid         IN     HDF5 identifier for the sequence group 
 *          *time         OUT    Time attribute value 
 *                              
 * Returns: Non-negative number if an error occurs, zero otherwise.
 *         
 * Errors: Identifier is checked against valid id values, check the existence
 *         of the attribute name.
 *        
*/
herr_t sequence_get_time(hid_t nsid, double *time)
{
  herr_t stat = DANU_FAILURE;

  if ( H5_ISA_VALID_ID(nsid) ) {
    stat=danu_attr_read_double(nsid,SEQ_TIME_ATTR_NAME,time);
  }
  else {
    DANU_ERROR_MESS("Invalid sequence group identifier.");
  }

  return stat;

}
/*
 * Routine: sequence_get_cycle(nsid, int * cycle)
 * Purpose: Return the cycle attribute associated with sequence group nsid
 * Description: Given a identifier for an existing sequence group, return the 
 *              cycle attribute associated with this groups. Returns a
 *              non-negative value if an error is detected. 
 *
 * Parameters:
 *           nsid         IN     HDF5 identifier for the sequence group 
 *          *cycle        OUT    Cycle attribute value 
 *                              
 * Returns: Non-negative number if an error occurs, zero otherwise.
 *         
 * Errors: Identifier is checked against valid id values, check the existence
 *         of the attribute name.
 *        
*/
herr_t sequence_get_cycle(hid_t nsid, int *cycle)
{
  herr_t stat = DANU_FAILURE;

  if ( H5_ISA_VALID_ID(nsid) ) {
    stat=danu_attr_read_int(nsid,SEQ_CYCNUM_ATTR_NAME,cycle);
  }
  else {
    DANU_ERROR_MESS("Invalid sequence group identifier.");
  }

  return stat;

}
/*
 * Routine: sequence_get_id(nsid, int * id)
 * Purpose: Return the sequence id attribute associated with sequence group nsid
 * Description: Given a identifier for an existing sequence group, return the 
 *              sequence id attribute associated with this groups. Returns a
 *              non-negative value if an error is detected. 
 *
 * Parameters:
 *           nsid      IN     HDF5 identifier for the sequence group 
 *          *id        OUT    Sequence id attribute value 
 *                              
 * Returns: Non-negative number if an error occurs, zero otherwise.
 *         
 * Errors: Identifier is checked against valid id values, check the existence
 *         of the attribute name.
 *        
*/
herr_t sequence_get_id(hid_t nsid, int *id)
{
  herr_t stat = DANU_FAILURE;

  if ( H5_ISA_VALID_ID(nsid) ) {
    stat=danu_attr_read_int(nsid,SEQ_SEQNUM_ATTR_NAME,id);
  }
  else {
    DANU_ERROR_MESS("Invalid sequence group identifier.");
  }

  return stat;

}





