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
 * output.c
 *
 *  TOUT Root output file
 *
 *
 *  Purpose:
 *           This source file creates, opens and closes the root output file. 
 *
 */
#include <danu_error.h>
#include <danu_memory.h>
#include <danu_file.h>
#include <danu_group.h>
#include <danu_dataset.h>

#include <danu_mesh.h>
#include <danu_sim.h>

#include <danu_output.h>

/* Private Defines */

/* Private prototypes */
herr_t output_file_initial_setup(hid_t fid);

/*
* Routine: herr_t output_file_initial_setup(hid_t fid)
* Purpose: Create the simulations and mesh root groups
* Description: Create the root simulation and mesh groups in output file
*              fid. Both groups are closed before returning.
*
* Parameters:
*           
*           fid     IN     HDF5 identifier for output file
*                              
* Returns: A negative value if an error occurs, otherwise returns a zero. 
* Errors: Possible error conditions are failure to create the groups or
*          invalid file id.
*     
*/
herr_t output_file_initial_setup(hid_t fid)
{
    herr_t stat1, stat2;

    stat1 = simulations_create_root_group(fid);
    stat2 = mesh_create_root_group(fid);

    if ( H5_RETURN_FAIL(stat1) ) {
        DANU_ERROR_MESS("Failed to create the simulation root group");
    }

    if ( H5_RETURN_FAIL(stat2) ) {
        DANU_ERROR_MESS("Failed to create the mesh root group");
    }

    return DANU_RETURN_SELECT(stat1,stat2);
}
/*
* Routine: hbool_t output_file_is_valid(hid_t fid)
* Purpose: Check root output file structure.
* Description: Check the existence of the simulation and mesh root groups
*              If either are not present the file is not valid.
*
* Parameters:
*           
*           fid     IN     HDF5 identifier for output file
*                              
* Returns: TRUE if file passes all tests, otherwise returns FALSE. 
* Errors: Possible error conditions are failure to stat the groups, in this
*         case the return value is always FALSE.
*     
*/
hbool_t output_file_is_valid(hid_t fid)
{
    hbool_t flag = FALSE;
    hbool_t f1,f2;
    hid_t sim = simulations_open_root_group(fid);
    hid_t mesh = mesh_open_root_group(fid);

    f1 = FALSE;
    if ( H5_ISA_INVALID_ID(sim) ) {
        DANU_ERROR_MESS("Failed to open simulations root group");
    }
    else {
        danu_group_close(sim);
        f1 = TRUE;
    }
 
    f2 = FALSE;
    if ( H5_ISA_INVALID_ID(mesh) ) {
        DANU_ERROR_MESS("Failed to open mesh root group");
    }
    else {
        danu_group_close(mesh);
        f2 = TRUE;
    }

    if ( f1 == TRUE && f2 == TRUE ) {
        flag = TRUE;
    }

    return flag;
}


/*
* Routine: hid_t output_file_open(const char * filename, unsigned access, unsigned action)
* Purpose: Create,open and possibly clobber the root output file
* Description: This is the driver routine for all *create and *open* calls. The action
*              flag is either open (DANU_FILE_ACT_OPEN) or create (DANU_FILE_ACT_CREATE).
*              The access flag is either read only (DANU_FILE_ACC_RDONLY), 
*              read/write (DANU_FILE_ACC_RDWR) or append (DANU_FILE_ACC_APPEND). 
*              If the file action is successful then the HDF5 identifier for the file
*              is returned. An invalid id will be returned if an error occurs. If the 
*              calling routine is a create action and the file exists, then the file will
*              be deleted with a warning message.
*
* Parameters:
*           filename      IN         File name
*           access        IN         File access type: rdonly,rdwr or append.    
*           action        IN         File action type: open, create EXITING FILE WILL BE CLOBBERED
*                              
* Returns: An HDF5 identifier for the file, this will be an invalif value if an error occurs. 
* Errors: Possible error conditions are bad string pointers, bad flag values for access, action
*         and failed to perform the requestion file access and action.
*/ 
hid_t output_file_open(const char *filename, unsigned access, unsigned action)
{
    hid_t fid = H5I_INVALID_HID;
    herr_t setup;


    if ( ACCESS_IS_INVALID(access) ) {
        DANU_ERROR_MESS("Invalid access flag");
        return fid;
    }

    if ( ACTION_IS_INVALID(action) ) {
        DANU_ERROR_MESS("Invalid action flag");
        return fid;
    }

    if ( DANU_BAD_STRING(filename) ) {
        DANU_ERROR_MESS("Invalid file name pointer");
        return fid;
    }


    if ( action == DANU_FILE_ACT_CREATE ) {

	/* Create the file and setup the groups */
	fid = danu_file_create(filename);
	setup = output_file_initial_setup(fid); /* fid is checked here! */
        if ( H5_RETURN_FAIL(setup) ) {
	    DANU_ERROR_MESS("Failed to create and initialize output file");
	    if ( H5_ISA_VALID_ID(fid) ) {
		danu_file_close(fid);
		fid = H5I_INVALID_HID;
	    }
	}
    }
    else {

	fid = danu_file_open(filename,access,action);
	if ( H5_ISA_INVALID_ID(fid) || ! output_file_is_valid(fid) ) {
	    if ( H5_ISA_INVALID_ID(fid) ) {
		DANU_ERROR_MESS("Failed to open output file");
	    }
	    else {
		danu_file_close(fid);
		DANU_ERROR_MESS("Invalid output file");
		fid = H5I_INVALID_HID;
	    }
	}
    }

    return fid;
}
/*
* Routine: herr_t output_file_create(const char * filename, hid_t * fid)
* Purpose: Create the root output file
* Description: Create and clobber, if the file exists, a root output file.
*              The file id is intitially set to an invalid id. Calling routine
*              should check return flag before using fid.
*
* Parameters:
*           filename      IN         File name
*           fid           OUT        HDF5 file id
*                              
* Returns: A negative value if an error occurs, otherwise returns a zero. 
* Errors: See output_file_open for a full list of possible error conditions.
*/
herr_t output_file_create(const char * filename, hid_t *fid)
{
    herr_t status = DANU_FAILURE;

    if ( DANU_BAD_PTR(fid) ) {
        DANU_ERROR_MESS("Invalid file id pointer");
        return status;
    }

    *fid = output_file_open(filename,DANU_FILE_ACC_RDWR,DANU_FILE_ACT_CREATE);

    if ( H5_ISA_VALID_ID(*fid) ) {
        status = 0;
    }

    return status;
}
/*
* Routine: herr_t output_file_open_rdonly(const char * filename, hid_t * fid)
* Purpose: Open root output file with read only access
* Description: Open a root output file with read only access. The simulation and 
*              mesh groups are checked before the routine returns. If either
*              group is not found the file is closed and the return value
*              is set to a neagtive value.
*
* Parameters:
*           filename      IN         File name
*           fid           OUT        HDF5 file id
*                              
* Returns: A negative value if an error occurs, otherwise returns a zero. 
* Errors: See output_file_open for a full list of possible file access error conditions.
*         If either the simulation or mesh root group does not exist, an error is raised.
*/
herr_t output_file_open_rdonly(const char * filename, hid_t *fid)
{
    herr_t status = DANU_FAILURE;

    if ( DANU_BAD_PTR(fid) ) {
        DANU_ERROR_MESS("Invalid file id pointer");
        return status;
    }

    *fid = output_file_open(filename,DANU_FILE_ACC_RDONLY,DANU_FILE_ACT_OPEN);

    if ( H5_ISA_VALID_ID(*fid) ) {

        if ( output_file_is_valid(*fid) ) {
            status = 0;
        }
        else {
            DANU_ERROR_MESS("Bad file format");
            danu_file_close(*fid);
        }
    }

    return status;
}
/*
* Routine: herr_t output_file_open_rdwr(const char * filename, hid_t * fid)
* Purpose: Open root output file with read and write access
* Description: Open a root output file with read and write access. The simulation and 
*              mesh groups are checked before the routine returns. If either
*              group is not found the file is closed and the return value
*              is set to a neagtive value.
*
* Parameters:
*           filename      IN         File name
*           fid           OUT        HDF5 file id
*                              
* Returns: A negative value if an error occurs, otherwise returns a zero. 
* Errors: See output_file_open for a full list of possible file access error conditions.
*         If either the simulation or mesh root group does not exist, an error is raised.
*/
herr_t output_file_open_rdwr(const char * filename, hid_t *fid)
{
    herr_t status = DANU_FAILURE;

    if ( DANU_BAD_PTR(fid) ) {
        DANU_ERROR_MESS("Invalid file id pointer");
        return status;
    }

    *fid = output_file_open(filename,DANU_FILE_ACC_RDWR,DANU_FILE_ACT_OPEN);

    if ( H5_ISA_VALID_ID(*fid) ) {

        if ( output_file_is_valid(*fid) ) {
            status = 0;
        }
        else {
            DANU_ERROR_MESS("Bad file format");
            danu_file_close(*fid);
        }
    }

    return status;
}
/*
* Routine: herr_t output_file_open_append(const char * filename, hid_t * fid)
* Purpose: Open root output file to apend data
* Description: Open a root output file with read and write access. The simulation and 
*              mesh groups are checked before the routine returns. If either
*              group is not found the file is closed and the return value
*              is set to a neagtive value.
*
* Parameters:
*           filename      IN         File name
*           fid           OUT        HDF5 file id
*                              
* Returns: A negative value if an error occurs, otherwise returns a zero. 
* Errors: See output_file_open for a full list of possible file access error conditions.
*         If either the simulation or mesh root group does not exist, an error is raised.
*/
herr_t output_file_open_append(const char * filename, hid_t *fid)
{
    herr_t status = DANU_FAILURE;

    if ( DANU_BAD_PTR(fid) ) {
        DANU_ERROR_MESS("Invalid file id pointer");
        return status;
    }

    *fid = output_file_open(filename,DANU_FILE_ACC_APPEND,DANU_FILE_ACT_OPEN);

    if ( H5_ISA_VALID_ID(*fid) ) {

        if ( output_file_is_valid(*fid) ) {
            status = 0;
        }
        else {
            DANU_ERROR_MESS("Bad file format");
            danu_file_close(*fid);
        }
    }

    return status;
}
/*
* Routine: herr_t output_file_close(hid_t fid)
* Purpose: Close root output file
* Description:  Close root output file. Once the routine returns any access using
*               fid will not be valid. File is flushed before the close is called.
*
* Parameters:
*           fid           IN        HDF5 file id
*                              
* Returns: A negative value if an error occurs, otherwise returns a zero. 
* Errors: Possible error conditions are invalid file identifier, fail to flush file
*         or fail to close the file.
*        
*/
herr_t output_file_close(hid_t *fid)
{
    herr_t status = DANU_FAILURE;

    if ( H5_ISA_INVALID_ID(*fid) ) {
        DANU_ERROR_MESS("Invalid file identifier");
        return status;
    }

    if ( H5_RETURN_FAIL(danu_file_flush_global(*fid)) ) {
        DANU_ERROR_MESS("Failed to flush output file");
        return status;
    }

    return danu_file_close(*fid);
}
