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
 * link.c
 *
 *  DANU Links
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

#include <hdf5.h>

#include <danu_error.h>
#include <danu_h5_error.h>
#include <danu_types.h>
#include <danu_memory.h>
#include <danu_utils.h>
#include <danu_group.h>

#include <danu_link.h>

/* Global define's */

/*
 * Routine: hid_t danu_link_create_soft(hid_t link_loc, const char * link_name, const char * target)
 * Purpose: Create a soft link name link_name that points to target.
 * Description: Creates a HDF5 soft link named link_name located in link_loc pointing
 *              to target. The target object named target may not exist. Soft links
 *              are allowed to dangle. Default property lists used in the create call.
 *
 * Parameters:
 *           link_loc   IN      HDF5 id for link location
 *           link_name  IN      name of the link        
 *           target     IN      Target name that link refers to
 *                              
 * Returns: Returns a HDF5 return value of the H5Lcreate_soft call. Will return a
 *          negative value if an error is raised.
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately. Returns a neagtive value if an error occurs.
 *
 */ 
hid_t danu_link_create_soft(hid_t link_loc, const char *link_name, const char *target )
{


    /* Check input */
    if ( H5_ISA_INVALID_ID(link_loc) ) {
       DANU_ERROR_MESS("Invalid HDF5 identifier argument");
       return DANU_FAILURE;
    }

    if ( DANU_BAD_PTR(link_name)) {
        DANU_ERROR_MESS("Invalid pointer to link_name");
        return DANU_FAILURE;
    }
    
    if ( DANU_BAD_PTR(target)) {
        DANU_ERROR_MESS("Invalid pointer to target");
        return DANU_FAILURE;
    }

    if ( DANU_EMPTY_STRING(link_name)) {
        DANU_ERROR_MESS("link_name is an empty string");
        return DANU_FAILURE;
    }
    
    if ( DANU_EMPTY_STRING(target)) {
        DANU_ERROR_MESS("target is an empty string");
        return DANU_FAILURE;
    }

    return H5Lcreate_soft(target, link_loc,link_name, H5P_DEFAULT,H5P_DEFAULT);
}

/*
 * Routine: hid_t danu_link_create_hard(hid_t link_loc, const char * link_name, hid_t obj_loc, const char * obj_name)
 * Purpose: Create a hard link name link_name that points to target defined by obj_loc and obj_name.
 * Description: Creates a HDF5 hard link named link_name located in link_loc pointing
 *              to target defined by obj_loc and obj_name. Names are interpeted relative to the 
 *              locations. The target MUST be located in the same file as link_loc. Default
 *              property list parameters are used in the create call. An obj_name pointer that is
 *              is interpeted as '.', that is the object location is the target link.
 *
 * Parameters:
 *           link_loc   IN      HDF5 id for link location
 *           link_name  IN      name of the link       
 *           obj_id     IN      Target HDF5 location identifier
 *           obj_name   IN      Target name relative to the object location
 *                              
 * Returns: Returns a HDF5 return value of the H5Lcreate_hard call. Will return a
 *          negative value if an error is raised.
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately. Returns a neagtive value if an error occurs.
 *
 */ 
hid_t danu_link_create_hard(hid_t link_loc, const char *link_name, hid_t obj_loc, const char *obj_name )
{
    /* Check input */
    if ( H5_ISA_INVALID_ID(link_loc) ) {
       DANU_ERROR_MESS("Invalid HDF5 identifier argument for link_loc");
       return DANU_FAILURE;
    }

    if ( DANU_BAD_PTR(link_name)) {
        DANU_ERROR_MESS("Invalid pointer to link_name");
        return DANU_FAILURE;
    }
    
    if ( H5_ISA_INVALID_ID(obj_loc)) {
        DANU_ERROR_MESS("Invalid HDF5 identifier argument for obj_loc");
        return DANU_FAILURE;
    }

    if ( DANU_EMPTY_STRING(obj_name) ) {
        return H5Lcreate_hard(obj_loc,".",link_loc,link_name, H5P_DEFAULT,H5P_DEFAULT);
    }
    else {
        return H5Lcreate_hard(obj_loc,obj_name,link_loc,link_name, H5P_DEFAULT,H5P_DEFAULT);
    }
}
/*
 * Routine: hid_t danu_link_create_hard(hid_t link_loc, const char * link_name, const name *file, const char * obj_name)
 * Purpose: Create an external link name link_name that points to target named obj_name in file.
 * Description: Creates an HDF5 external link named link_name located in link_loc pointing
 *              to target name obj_name in external file. Names are interpeted relative to the 
 *              locations. The file and obj_name in file  MUST exist. Default
 *              property list parameters are used in the create call. An obj_name pointer that is
 *
 * Parameters:
 *           link_loc   IN      HDF5 id for link location
 *           link_name  IN      name of the link       
 *           file       IN      External file name
 *           obj_name   IN      Target named obj_name in file
 *                              
 * Returns: Returns a HDF5 return value of the H5Lcreate_external call. Will return a
 *          negative value if an error is raised.
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately. Returns a neagtive value if an error occurs.
 *
 */ 
hid_t danu_link_create_external(hid_t link_loc, const char *link_name, const char *file, const char *obj_name )
{
    /* Check input */
    if ( H5_ISA_INVALID_ID(link_loc) ) {
       DANU_ERROR_MESS("Invalid HDF5 identifier argument for link_loc");
       return DANU_FAILURE;
    }

    if ( DANU_BAD_PTR(link_name)) {
        DANU_ERROR_MESS("Invalid pointer to link_name");
        return DANU_FAILURE;
    }

    if ( DANU_EMPTY_STRING(link_name) ) {
        DANU_ERROR_MESS("Empty string for link_name");
        return DANU_FAILURE;
    }
    
    if ( DANU_BAD_PTR(file) ) {
        DANU_ERROR_MESS("Invalid pointer to file");
        return DANU_FAILURE;
    }

    if ( DANU_EMPTY_STRING(file) ) {
        DANU_ERROR_MESS("Empty string for file");
        return DANU_FAILURE;
    }

    if ( DANU_BAD_PTR(obj_name) ) {
        DANU_ERROR_MESS("Invalid pointer to obj_name");
        return DANU_FAILURE;
    }

    if ( DANU_EMPTY_STRING(obj_name) ) {
        DANU_ERROR_MESS("Empty string for obj_name");
        return DANU_FAILURE;
    }

    if ( ! FILE_EXISTS(file) ) {
        DANU_ERROR_MESS("Attempting to create link to external file that does not exist");
        return DANU_FAILURE;
    }

    return H5Lcreate_external(file,obj_name,link_loc,link_name,H5P_DEFAULT,H5P_DEFAULT);

}

herr_t danu_link_exists(hid_t loc_id, const char * name, hbool_t * flag)
{
    herr_t status = DANU_FAILURE;
    htri_t ret;
    hid_t  id;

    /* Default is FALSE */
    if ( DANU_BAD_PTR(flag) ) {
        DANU_ERROR_MESS("Invalid flag pointer");
        return status;
    }
    *flag = FALSE;

    /* Check input */
    if ( ! H5_ISA_GROUP_ID(loc_id) && ! H5_ISA_FILE_ID(loc_id) ) {
        danu_print_hid_info(loc_id);
        DANU_ERROR_MESS("Invalid HDF5 location identifier");
        return status;
    }

    if ( DANU_BAD_STRING(name) ) {
        DANU_ERROR_MESS("Invalid name string pointer");
        return status;
    }

    /* Want to query on file id's 
       to prevent a recursive call use the HDF5 group open call 
       here. Assume loc id is a group id, otherwise open the  
       root group in a file.
     */
    id = loc_id; 
    if ( H5_ISA_FILE_ID(loc_id) ) {
        id = H5Gopen(loc_id,"/",H5P_DEFAULT);
        if ( H5_ISA_INVALID_ID(id) ) {
            DANU_ERROR_MESS("Failed to open file root group");
            return status;
        }
    }

    ret = H5Lexists(id,name,H5P_DEFAULT);

    /* Set the flag and return status */
    if ( ret >= 0 ) {
        *flag = ret;
        status = 0;
    }
    else {
        *flag = FALSE;
        status = DANU_FAILURE;
    }

    /* Close the root group if needed */
    if ( id != loc_id ) {
        danu_group_close(id);
    }

    return status;
}

    


    
    
    
        



