
/*
 * group.c
 *
 *  DANU Groups
 *
 *
 *  Purpose:
 *
 *          The HDF5 library
 *
 *
 */

#if HAVE_CONFIG_H
# include <danu_config.h>
#endif

#include <string.h>

#include <hdf5.h>

#include <danu_error.h>
#include <danu_h5_error.h>
#include <danu_h5_object.h>
#include <danu_types.h>
#include <danu_memory.h>
#include <danu_utils.h>
#include <danu_link.h>
#include <danu_file.h>

#include <danu_group.h>

/* Global define's */

/*
 * Routine: hid_t danu_group_create(hid_t loc,char *name)
 * Purpose: Create a HDF5 group in HDF5 object identified by loc
 * Description: Creates a an HDF5 group in location identified by loc. Will
 *              create intermediate groups if those groups do not exist and
 *              the index of links to order created. The default create and
 *              access property lists are used.
 *
 * Parameters:
 *           loc                IN              HDF5 id for group location (file,group)
 *           name               IN              Name of the group
 *                              
 * Returns: Returns a HDF5 identifier for the new group.
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately. Returns a neagtive value if an error occurs.
 *
 */ 
hid_t danu_group_create(hid_t loc, const char *name)
{
	hid_t group;

	hid_t lpl,cpl;

	/* Check input */
	if ( H5_ISA_INVALID_ID(loc) ) {
           DANU_ERROR_MESS("Invalid HDF5 identifier argument");
           return H5I_INVALID_HID;
	}

	if ( DANU_BAD_PTR(name) ) {
	   DANU_ERROR_MESS("Invalid group name argument");
           return H5I_INVALID_HID;
	}

	if ( DANU_EMPTY_STRING(name) ) {
	   DANU_ERROR_MESS("Invalid group name argument");
           return H5I_INVALID_HID;
	}

	/* Define the create property list ... create intermediate groups */

	cpl = H5Pcreate(H5P_GROUP_CREATE);
        lpl = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(lpl,TRUE);

	group = H5Gcreate(loc,name,lpl,cpl,H5P_DEFAULT);

	/* Release the link and create property list */
	H5Pclose(cpl);
	H5Pclose(lpl);

	return group;

}
/*
 * Routine: hid_t danu_group_open(hid_t loc,char *name)
 * Purpose: Opens an existing HDF5 group in HDF5 object identified by loc
 * Description: Opens an HDF5 group in location identified by loc. Uses the
 *              default access proerties. 
 *
 * Parameters:
 *           loc                IN              HDF5 id for group location (file,group)
 *           name               IN              Name of the group
 *                              
 * Returns: Returns a HDF5 identifier for the group.
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Returns a negative value if an error occurs.
 *
 */ 
hid_t danu_group_open(hid_t loc, const char * name )
{
	hid_t group = H5I_INVALID_HID;

	/* Check input */
        if ( H5_ISA_INVALID_ID(loc) ) {
           DANU_ERROR_MESS("Invalid HDF5 identifier argument");
           return group;
	}

	if (DANU_BAD_STRING(name) ) {
           DANU_ERROR_MESS("Invalid group name argument");
           return group;
	}


	/* Now open the group */

        if ( danu_group_exists(loc,name) ) {
	    group = H5Gopen(loc,name,H5P_DEFAULT);
        }
        else {
            DANU_ERROR_MESS("Group does not exist");
        }

	return group;

}
/*
 * Routine: hid_t danu_group_close(hid_t id)
 * Purpose: Closes an HDF5 group identified by id 
 * Description: Close an HDF5 group object identified with id. The input
 *              is checked before the close is called. If the id is not valid the 
 *              routine returns immediately. 
 *
 * Parameters:
 *           loc                id              HDF5 id for group
 *                              
 * Returns: Returns the returning value of H5Gclose
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Returns a negative value if an error occurs.
 *
 */ 
herr_t danu_group_close(hid_t id )
{
	/* Check input */
	if ( H5_ISA_INVALID_ID(id) ) {
           DANU_ERROR_MESS("Invalid HDF5 identifier argument");
	   return H5I_INVALID_HID;
	}

	return H5Gclose(id);

}
/*
 * Routine: herr_t danu_group_get_nlinks(hid_t gid, hsize_t *nlinks)
 * Purpose: Return the number of links found under gid
 * Description: Calls H5Gget_info and returns the nlinks value in the group 
 *              info struct. Will raise an error if input is not valid.
 *              Return is negative if an error occurs.
 *
 * Parameters:
 *           gid                IN              HDF5 id for group
 *          *nlinks             OUT             Number of links
 *                              
 * Returns: Returns the returning value of H5Gget_info
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Returns a negative value if an error occurs.
 *
 */
herr_t danu_group_get_nlinks(hid_t gid, hsize_t *nlinks)
{
    herr_t status = DANU_FAILURE;

    H5G_info_t info;
    
    /* Check input */
    if ( ! H5_ISA_GROUP_ID(gid) ) {
        DANU_ERROR_MESS("Invalid HDf5 group identifier");
        return status;
    }

    if ( DANU_BAD_PTR(nlinks) ) {
        DANU_ERROR_MESS("Invalid pointer to nlinks");
        return status;
    }

    status = H5Gget_info(gid,&info);

    if ( status >= 0 ) {
        *nlinks = info.nlinks;
    }

    return status;
}
/*
 * Routine: herr_t danu_group_get_nlinks_by_name(hid_t loc_id,const char *name, hsize_t *nlinks)
 * Purpose: Return the number of links found under group named 'name' found in loc_id
 * Description: Calls H5Gget_info_by_name and returns the nlinks value in the group 
 *              info struct. Will raise an error if input is not valid.
 *              Return is negative if an error occurs.
 *
 * Parameters:
 *           loc_id             IN              HDF5 id for location
 *           name               IN              Name of group
 *          *nlinks             OUT             Number of links
 *                              
 * Returns: Returns the returning value of H5Gget_info
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Returns a negative value if an error occurs.
 *
 */
herr_t danu_group_get_nlinks_by_name(hid_t loc_id, const char *name, hsize_t *nlinks)
{
    herr_t status = DANU_FAILURE;

    H5G_info_t info;
    
    /* Check input */
    if ( H5_ISA_INVALID_ID(loc_id) ) {
        DANU_ERROR_MESS("Invalid HDF5 location identifier");
        return status;
    }
    
    if ( DANU_BAD_STRING(name) ) {
        DANU_ERROR_MESS("Invalid pointer to group name");
        return status;
    }


    if ( DANU_BAD_PTR(nlinks) ) {
        DANU_ERROR_MESS("Invalid pointer to nlinks");
        return status;
    }

    status = H5Gget_info_by_name(loc_id,name,&info,H5P_DEFAULT);

    if ( status >= 0 ) {
        *nlinks = info.nlinks;
    }

    return status;
}

/*
 * Routine: herr_t danu_group_get_subgroups(hid_t gid, int num,  char **subgroups, int*num_found)
 * Purpose: Search for subgroups under gid and return the name of those groups and the number found.
 * Description: Searchs the group (gid) for objects that are groups. This search is NOT recursive. Returns
 *              the subgroups array with the names of the subgroups found and num_found is the number of
 *              groups found. Routine allocates memory for each subgroup name. The calling
 *              routine is responsible for freeing this memory.
 *
 * Parameters:
 *           loc_id             IN              HDF5 id for location
 *           num                IN              Number of pointers in subgroups array
 *           **subgroups        INOUT           Array of num pointers
 *          *num_found          OUT             Number of groups found
 *                              
 * Returns: Returns the returning value of H5Gget_info
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Returns a negative value if an error occurs.
 *
 */
herr_t danu_group_get_subgroups(hid_t gid, int num, char **subgroups, int *num_found)
{
    herr_t status = DANU_FAILURE;

    h5o_node_list_t *list;
    int               num_cpy;

    /* Check input */
    if ( ! H5_ISA_GROUP_ID(gid) ) {
        DANU_ERROR_MESS("Invalid group identifier");
        return status;
    }

    if ( num <= 0 ) {
        DANU_ERROR_MESS("Invalid num (size of subgroups array) value");
        return status;
    }

    if ( DANU_BAD_PTR(subgroups) ) {
        DANU_ERROR_MESS("Invalid subgroups pointer");
        return status;
    }

    if ( DANU_BAD_PTR(num_found) ) {
        DANU_ERROR_MESS("Invalid subgroups pointer");
        return status;
    }
    *num_found = -1;


    /* Initialize the list */
    if ( NULL != ( list = h5_object_list_create(H5O_TYPE_GROUP) ) ) {
       
        status = h5_object_search(gid,H5O_TYPE_GROUP,list);

        if ( status >= 0 ) {

	  *num_found = (int) h5_object_list_nnodes(list);
	  if ( *num_found > 0 ) {
	    num_cpy = h5_object_list_copy_names(list,num,subgroups);
            if ( num_cpy <= 0 ) {
                DANU_ERROR_MESS("Copy of linked list failed");
                status = DANU_FAILURE;
            }
            else {
                status = DANU_SUCCESS;
            }
	  }
	  else if ( *num_found == 0 ) {
	    memset(subgroups,0,sizeof(char*)*num);
            status = DANU_SUCCESS;
	  }
	  else {
	    DANU_ERROR_MESS("Unknown error state....should not be possible");
	    status = DANU_FAILURE;
	  }

	}  

        h5_object_list_delete(list);
    }

    return status;
}

/*
 * Routine: herr_t danu_group_get_datasets(hid_t gid, int num,  char **datasets, int*num_found)
 * Purpose: Search for datasets under gid and return the name of those datasets and the number found.
 * Description: Searchs the group (gid) for objects that are datasets. This search is NOT recursive. Returns
 *              the datasets array with the names of the subgroups found and num_found is the number of
 *              datasets found. The routine allocates memory for each pointer in datasets. The calling
 *              routine is responsible for freeing this memory, 
 *              Use danu_group_get_nlinks to create arrays of the correct size.
 *
 * Parameters:
 *           gid          IN    HDF5 group identifier.          
 *           num          IN    num length of the 1D array size and the number of string pointers dataset can hold
 *          *datasets     OUT   Array of string pointers
 *          *num_found    OUT   Number of datasets found
 *                              
 * Returns: Returns a negative value if the input is not valid. If the input is valid returns the 
 *          returning value of H5Literate.
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Returns a negative value if an error occurs.
 *
 */
herr_t danu_group_get_datasets(hid_t gid, int num, char **datasets, int *num_found)
{
    herr_t status = DANU_FAILURE;

    h5o_node_list_t *list;
    int              num_cpy;
    char             null_char = '\0';

    /* Check input */
    if ( ! H5_ISA_GROUP_ID(gid) ) {
        DANU_ERROR_MESS("Invalid group identifier");
        return status;
    }

    if ( num <= 0 ) {
        DANU_ERROR_MESS("Invalid num (size of datasets array) value");
        return status;
    }

    if ( DANU_BAD_PTR(datasets) ) {
        DANU_ERROR_MESS("Invalid datasets pointer");
        return status;
    }

    if ( DANU_BAD_PTR(num_found) ) {
        DANU_ERROR_MESS("Invalid num_found pointer");
        return status;
    }
    *num_found = -1;


    /* Initialize the list */
    if ( NULL != ( list = h5_object_list_create(H5O_TYPE_DATASET) ) ) {

        status = h5_object_search(gid,H5O_TYPE_DATASET,list);

        if ( status >= 0 ) {

	  *num_found = h5_object_list_nnodes(list);
          if ( *num_found > 0 ) {
	    num_cpy = h5_object_list_copy_names(list,num,datasets);
            if ( num_cpy <= 0 ) {
                DANU_ERROR_MESS("Copy of dataset names from linked list failed");
            }
            else {
                status = DANU_SUCCESS;
            }
                
	  }
	  else if ( *num_found == 0 ) {
	    memcpy(datasets,&null_char,sizeof(char*)*num);
            status = DANU_SUCCESS;
	  }
	  else {
	    DANU_ERROR_MESS("Unknown error state....should not be possible");
	    status = DANU_FAILURE;
	  }

	}
           
        h5_object_list_delete(list);
    }
                    

    return status;
}


hbool_t danu_group_exists(hid_t loc_id, const char * name)
{
    hbool_t flag;
    herr_t ret;

    ret = danu_link_exists(loc_id,name,&flag);

    if ( ! H5_RETURN_OK(ret) ) {
        DANU_ERROR_MESS("Failed to status group");
    }
   
    return flag;
}
/*
 * Routine: herr_t danu_group_find_target(hid_t gid, const char * target_name, H5O_type_t type, int *found)
 * Purpose: Search for H5 object named target_name under group identified by gid
 * Description: Searchs the group (gid) for objects of type 'type' with name target_name. The search name is
 *              relative to the group name used to retrieve gid. Flag found will be set to either TRUE or FALSE.
 *
 * Parameters:
 *           gid          IN    HDF5 group identifier.          
 *           target_name  IN    Target name to search
 *           type         IN    Target H5 object type, only group and dataset supported
 *          *found        OUT   Flag set to TRUE or FALSE
 *                              
 * Returns: Returns a negative value if an error occurs. Default setting for found is FALSE. Calling code 
 *          should check the return status for errors NOT the found flag.
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Return value will be negative if an error occured.
 *
 */
herr_t danu_group_find_target(hid_t loc, const char * target, H5O_type_t type, int *found)
{
    herr_t status = DANU_FAILURE;
 
    h5o_node_list_t *list;

    /* Initialize the pointer */
    if ( DANU_BAD_PTR(found) ) {
        DANU_ERROR_MESS("Invalid pointer input");
        return status;
    }

    *found = FALSE;

    /* Define the list */
    if ( type >= H5O_TYPE_NTYPES || type < 0) {
        DANU_ERROR_MESS("Invalid HDF5 object type");
        return status;
    }

    if ( NULL != ( list = h5_object_list_create(type) ) ) {

        status = h5_object_search_by_name(loc,target,type,list);

        if ( status >= 0 ) {
            
            if ( list->root != NULL ) {
                *found = TRUE;
            }

        }

        h5_object_list_delete(list);
    }
    else {
        DANU_ERROR_MESS("Failed to allocate memory for object list");
        status = -1;
    }

    return status;

}
/*
 * Routine: herr_t danu_group_find_target_subgroup(hid_t gid, const char * target_name, int *found)
 * Purpose: Search for H5 subgroup named target_name under group identified by gid
 * Description: Wrapper for the danu_group_find_target subroutine for subgroups.
 *
 * Parameters:
 *           gid          IN    HDF5 group identifier.          
 *           target_name  IN    Target subgroup name to search
 *          *found        OUT   Flag set to TRUE or FALSE
 *                              
 * Returns: Returns a negative value if an error occurs. Default setting for found is FALSE. Calling code 
 *          should check the return status for errors NOT the found flag.
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Return value will be negative if an error occured.
 *
 */
herr_t danu_group_find_target_subgroup(hid_t loc, const char * target, int *found)
{
    return danu_group_find_target(loc,target,H5O_TYPE_GROUP,found);
}
/*
 * Routine: herr_t danu_group_find_target_dataset(hid_t gid, const char * target_name, int *found)
 * Purpose: Search for H5 dataset named target_name under group identified by gid
 * Description: Wrapper for the danu_group_find_target subroutine for datasets.
 *
 * Parameters:
 *           gid          IN    HDF5 group identifier.          
 *           target_name  IN    Target dataset name to search
 *          *found        OUT   Flag set to TRUE or FALSE
 *                              
 * Returns: Returns a negative value if an error occurs. Default setting for found is FALSE. Calling code 
 *          should check the return status for errors NOT the found flag.
 * Errors: The input is checked and if an error is detected in the input the routine
 *         returns immediately. Return value will be negative if an error occured.
 *
 */
herr_t danu_group_find_target_dataset(hid_t loc, const char * target, int *found)
{
    return danu_group_find_target(loc,target,H5O_TYPE_DATASET,found);
}






            



    

/*
 * Routine: hid_t danu_group_create_external(hid_t link_loc, const char *group_link, const char *file, const char *obj_name)
 * Purpose: Create an group that is an external link.
 * Description: Create an external link located under HDF5 object link_loc. External file must exist before this
 *              routine is called. If link is created successfully, the
 *              new link is opened and the HDF5 identifier for that new link is returned.
 *
 * Parameters:
 *           loc               IN              HDF5 id for link location
 *           group_link        IN              Group link name under link location
 *           file              IN              External file link points to
 *           obj_name          IN              External object name
 *                              
 * Returns: Returns the HDF5 identifier for the link if the link creation is successful. 
 * Errors: The input is check in danu_link_create_external.
 *
 */ 
hid_t danu_group_create_external(hid_t link_loc, const char *group_link, const char *file, const char *obj_name)
{
    hid_t group = H5I_INVALID_HID;  /* By default return an invalid group id */

    hid_t fid, gid;

    herr_t status;

    /* Check file name input */
    if ( DANU_BAD_PTR(file) ) {
        DANU_ERROR_MESS("Invliad pointer to file name");
        return H5I_INVALID_HID;
    }

    if (DANU_EMPTY_STRING(file) ) {
        DANU_ERROR_MESS("Empty string file name");
        return H5I_INVALID_HID;
    }

    /* Check the object name input */ 
    if ( DANU_BAD_PTR(obj_name) ) {
        DANU_ERROR_MESS("Invalid pointer to object name");
        return H5I_INVALID_HID;
    }

    if (DANU_EMPTY_STRING(obj_name) ) {
        DANU_ERROR_MESS("Empty string object name");
        return H5I_INVALID_HID;
    }

    /* Create the file if it does not exist */
    if ( FILE_EXISTS(file) ) {
        fid = danu_file_open(file,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_APPEND);
    }
    else {
        fid = danu_file_create(file);
    }

    /* Now create the object group */
    if ( H5_ISA_VALID_ID(fid) ) {
        gid = danu_group_create(fid,obj_name);
        if ( H5_ISA_INVALID_ID(gid) ) {
            DANU_ERROR_MESS("Failed to create object group");
        }
    }
    else {
        DANU_ERROR_MESS("Failed to create/open file");
        gid = H5I_INVALID_HID;
    }

    /* Create the linke if file and group are OK */
    if ( H5_ISA_VALID_ID(fid) && H5_ISA_VALID_ID(gid) ) {
        status = danu_link_create_external(link_loc,group_link,file,obj_name);
    }
    else {
        status = DANU_FAILURE;
    }

    /* Close the external file and object */
    //if ( H5_ISA_VALID_ID(fid) ) danu_file_close(fid);
    //if ( H5_ISA_VALID_ID(gid) ) danu_group_close(gid);


    if ( H5_RETURN_OK(status) ) {
        group = H5Gopen(link_loc,group_link,H5P_DEFAULT);
    }

    return group;
}





