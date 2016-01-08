/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

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

#include <stdarg.h>
#include <string.h>

#include <hdf5.h>

#include <danu.h>

#include <danu_non-series.h>
#include <danu_series-group.h>
#include <danu_series-data.h>
#include <danu_probes.h>
#include <danu_mesh.h>

#include <danu_sim.h>

/* PRIVATE define's */
#define SIM_ROOT_NAME "Simulations"


/*
 * Routine: simulations_create_root_group(hid_t fid)
 * Purpose: Create the Simulations root group under the root file 
 * Description: Create the Simulations group under the root file identified by
 *              fid. Will return error if fid is an invalid or if the group is not created
 *              successfully. Group will be closed before returning.
 *
 * Parameters:
 *           fid       IN              HDF5 identifier for root file 
 *                              
 * Returns: A flag (herr_t) is returned indicating the status. A negative return value 
 *          indicates an error.
 *          
 * Errors: An error is raised if fid is not a valid file identifier or if the group
 *         is not created.
 *        
 */
 herr_t simulations_create_root_group(hid_t fid)
 {
     herr_t ret = DANU_FAILURE;

     char sim_root_name[] = SIM_ROOT_NAME;
     hid_t gid = danu_group_create(fid,sim_root_name);

     if ( H5_ISA_GROUP_ID(gid) ) {
         danu_group_close(gid);
         ret = DANU_SUCCESS;
     }
     else {
         DANU_ERROR_MESS("Failed to create the root simulation group");
     }

     return ret;
 }
/*
 * Routine: simulations_open_root_group(hid_t fid)
 * Purpose: Open the Simulations root group under the root file 
 * Description: Open the Simulations group under the root file identified by
 *              fid. Will return error if fid is an invalid or if the group is not opened
 *              successfully.
 *
 * Parameters:
 *           fid       IN              HDF5 identifier for root file 
 *                              
 * Returns: The HDF5 identifier for the group is returned. A negative return value 
 *          indicates an error.
 *          
 * Errors: An error is raised if fid is not a valid file identifier or if the group
 *         is not opened.
 *        
 */
 hid_t simulations_open_root_group(hid_t fid)
 {
     char sim_root_name[] = SIM_ROOT_NAME;
     hid_t gid = danu_group_open(fid,sim_root_name);

     if ( ! H5_ISA_GROUP_ID(gid) ) {
         DANU_ERROR_MESS("Failed to open the root simulation group");
         gid = H5I_BADID;
     }

     return gid;
 }

 /*
 * Routine: simulation_exist(hid_t fid, const char * sim_name, int *exist)
 * Purpose: Test the existence of a simulation
 * Description: Tests the existence of the simulation group under HDF5 object fid.
 *              Since HDf5 does not provide existence functions, H5Lexist is used.
 *              This function returns a TRUE, FALSE or negative value if an error
 *              occurs. The exist flag will be set to either TRUE or FALSE. If the
 *              H5Lexists call returns an error exist is be set to FALSE. This 
 *              routine will not guarantee that an open will succeed. This call
 *              ONLY checks the existence of the link in fid it does not
 *              resolve the link. In a future HDF5 release, a function that
 *              resolves links will be implemented and this function will be updated. 
 *              
 *
 * Parameters:
 *           fid       IN              HDF5 identifier for root file 
 *           sim_name  IN              Name of simulation to search
 *           *exist    OUT             Flag set to TRUE (simulation exists) or 
 *                                      FALSE (simulation does not exist)   
 *                              
 * Returns: A flag (herr_t) is returned indicating the status. A negative return value 
 *          indicates an error otherwise 0 is returned.
 *          
 * Errors: An error is raised if fid is not a valid file identifier or sim_name
 *         is not a a valid char pointer.
 *      
 *        
 */
 herr_t simulation_exists(hid_t fid, const char * sim_name, int *exists)
 {
     herr_t ret_value = DANU_FAILURE;
     hbool_t flag;
     hid_t rid = simulations_open_root_group(fid);



     if ( H5_ISA_VALID_ID(rid) ) {
         flag = danu_group_exists(rid,sim_name);
         ret_value = DANU_SUCCESS;
         if (flag) {
             *exists = TRUE;
         }
         else {
             *exists = FALSE;
         }

     }

     return ret_value;
 }
/*
 * Routine: simulation_count(hid_t fid, int *cnt) 
 * Purpose: Return the number of simulation groups under the Simulations (root) group
 * Description: Opens the Simulations (root) group and returns the number of subgroups
 *              under that group. Function assumes each subgroup is a simulation group.
 *              Error is raised if the input is not valid. An error is also raised
 *              if the sizeof int and hsize_t are mismatched. In this case, the function
 *              will cast the number of links (hsize_t) into an int and print a wanring
 *              message. The return value is set to indicate an error.
 *              Return value is negative if an error is raised and cnt is set to zero.
 *
 * Parameters:
 *           fid      IN              HDF5 identifier for  output file
 *           *cnt     OUT             Number of subgroups under Simulations (root) group
 *                              
 * Returns: Returns a negative value if an error is encountered.
 *          
 * Errors: Input is checked when the number of links to Simulations (root) group is queried.
 *         Error also raised when there is a size mismatch between int and hsize_t.
 *        
 */
 herr_t simulation_count(hid_t fid, int *cnt)
 {
     herr_t ret_value = DANU_FAILURE;

     hid_t rid = simulations_open_root_group(fid);

     hsize_t nlinks;


     if ( H5_ISA_VALID_ID(rid) ) {
         if( H5_RETURN_OK(danu_group_get_nlinks(rid, &nlinks) ) ) {
             *cnt = (int) nlinks;
             ret_value = DANU_SUCCESS;
         }
         else {
             DANU_ERROR_MESS("Failed to count the number of simulations");
         }
         danu_group_close(rid);
     }
     else {
         DANU_ERROR_MESS("Failed to open the root simulation group");
     }

     return ret_value;
}
/*
 * Routine: simulation_list(hid_t fid, 
                            int num,
                            char **sim_names,
                            int *num_found)
 * Purpose: Return the simulation group names in sim_names array
 * Description: Given an output file descriptor (fid), return the simulation names
 *               stored in the sim_names array. Calling routine should call 
 *               simulation_count to determine the size of sim_names. This routine
 *               will allocate memory for each pointer in sim_name. The calling
 *               routine is responsible for freeing this memory. An error
 *               is raised if the input is not valid. 
 *
 * Parameters:
 *           fid        IN              HDF5 identifier for the output file 
 *           num        IN              Number of sim_name strings
 *           sim_names  OUT             Array of string pointers pointing to the 
 *                                       simulation names.
 *           *num_found OUT             Pointer to integer that returns the
 *                                       number of found simulations 
 *                              
 * Returns: Returns a negative value if an error occurs. 
 *          
 * Errors: Invalid input values will raise an error. An error will be raised if the
 *         Simulations group does not exist or the search for simulation groups encounters
 *         an error. 
 *        
 */
herr_t  simulation_list(hid_t fid, int num, char **sim_names, int *num_found)
{
    herr_t status = DANU_FAILURE;

    hid_t rid = simulations_open_root_group(fid);

    if ( H5_ISA_VALID_ID(rid) ) {
        status = danu_group_get_subgroups(rid,num,sim_names,num_found);
        if ( H5_RETURN_OK(status) ) {
            if ( *num_found > num ) {
                DANU_WARN_MESS("Number of simulation groups found exceeds the size of the array");
            }
        }
        else {
            DANU_ERROR_MESS("Failed to search subgroups in Simulation root group");
        }
        danu_group_close(rid);
    }
    else {
        DANU_ERROR_MESS("Failed to open the Simulations root group");
    }

    return status;
}
/*
 * Routine: simulation_add(hid_t fid, const char * sim_name, hid_t * sid)
 * Purpose: Create a simulation group structure under the root file fid
 * Description: Create a simulation group under the root file identified by
 *              fid. Will return error if fid is an invalid or if the group is not created
 *              successfully. Functions checks the existence of the simulation group and
 *              will return an error if the group already exists.
 *              All subgroups such as Non Series, Series and Probes will be created
 *              and closed in this routine. Failure to create any of these subgroups
 *              will result in an error.
 *
 *
 * Parameters:
 *           fid       IN              HDF5 identifier for root file 
 *           sim_name  IN              Simulation name
 *           *sid      OUT             HDF5 group identifier (pointer)
 *                              
 * Returns: A flag (herr_t) is returned indicating the status. A negative return value 
 *          indicates an error.
 *          
 * Errors: An error is raised if fid is not a valid file identifier or if the group
 *         is not created or if any of the sub groups are not created.
 *        
 */
 herr_t simulation_add(hid_t fid, const char * sim_name, hid_t *sid)
 {
     herr_t ret_value = DANU_FAILURE;
     herr_t num_errs  = 0;
     int exists = FALSE;
     hid_t rid = simulations_open_root_group(fid);

     /* Check input */
     if ( H5_RETURN_FAIL(simulation_exists(fid,sim_name,&exists) ) ) {
         DANU_ERROR_MESS("Failed to stat existence of simulation");
         return ret_value;
     }

     if ( H5_ISA_INVALID_ID(rid) ) {
         DANU_ERROR_MESS("Failed to open the root Simulation group");
         return ret_value;
     }

     if ( exists ) {
         DANU_ERROR_MESS("Simulation already exists, will not over write exiting group");
         return ret_value;
     }

     if ( DANU_BAD_PTR(sid) ) {
         DANU_ERROR_MESS("Invalid pointer");
         return ret_value;
     }


     /* Now create a simulation group structure */
     
     /* Create the top group  */
     *sid = H5I_BADID;
     *sid = danu_group_create(rid,sim_name);

     if ( H5_ISA_VALID_ID(*sid) ) {
         num_errs = 0;
         if ( H5_RETURN_FAIL(data_create_group(*sid) ) ) {
             DANU_ERROR_MESS("Failed to create non-series data group");
             num_errs++;
         }

         if ( H5_RETURN_FAIL(sequence_create_root_group(*sid) ) ) {
             DANU_ERROR_MESS("Failed to create series data group");
             num_errs++;
         }

         if ( H5_RETURN_FAIL(probe_create_group(*sid) ) ) {
             DANU_ERROR_MESS("Failed to create probe data group");
             num_errs++;
         }

         danu_group_close(rid);

     }
     else {
         DANU_ERROR_MESS("Failed to create the simulation group");
     }

     ret_value = (num_errs == 0 ) ? num_errs : DANU_FAILURE;

     return ret_value;


}   
/*
 * Routine: simulation_open(hid_t fid, const char * sim_name, hid_t * sid)
 * Purpose: Open an existing simulation group
 * Description: Open an existing Simulations group under the root file identified by
 *              fid. Will return error if fid is an invalid or if the Simulation group
 *              does not exist.
 *
 * Parameters:
 *           fid       IN              HDF5 identifier for root file 
 *           sim_name  IN              Simulation name
 *           *sid      OUT             HDF5 group identifier (pointer)
 *                              
 * Returns: A flag (herr_t) is returned indicating the status. A negative return value 
 *          indicates an error.
 *          
 * Errors: An error is raised if fid is not a valid file identifier or if the group
 *         is not created.
 *        
 */
herr_t simulation_open(hid_t fid, const char *sim_name, hid_t *sid)
{

    herr_t status = DANU_FAILURE;
    
    hid_t rid;
    int exists;

    if (DANU_BAD_PTR(sid) ) {
        DANU_ERROR_MESS("Invalid pointer");
        return status;
    }
    
    if ( H5_RETURN_OK(simulation_exists(fid,sim_name,&exists) ) ) {
        if (exists) {
            rid = simulations_open_root_group(fid);
            *sid = H5I_BADID;
            *sid = danu_group_open(rid,sim_name);
            if ( H5_ISA_INVALID_ID(*sid) ) {
                DANU_ERROR_MESS("Failed to open the simulation group");
            }
            danu_group_close(rid);
            status = 0;
        }
        else {
            DANU_ERROR_MESS("Can not simulation group; does not exist");
        }
    }
    else {
        DANU_ERROR_MESS("Failed to stat the existence of simulation group");
    }


    return status;
}
/*
 * Routine: herr_t simulation_link_mesh(hid_t fid, hid_t sid, const char * meshname)
 * Purpose: Link meshname to simulation sid. 
 * Description: Create a soft link to mesh group meshname in fid, under the simulation group
 *              labeled SIM_MESH_LINK_NAME. Since the link is soft, only a warning
 *              will occur if the meshname does not exist. 
 *
 * Parameters:
 *           fid       IN     HDF5 identifier for file
 *           sid       IN     Simulation identifier
 *           meshname  IN     Mesh name found in file fid 
 *                              
 * Returns: A flag (herr_t) is returned indicating the status. A negative return value 
 *          indicates an error.
 *          
 * Errors: An error is raised if fid or sid are not a valid, or if the link fails.
 *         Since the link is a soft link, only a warning will be issued if meshname
 *         does not exist.
 *        
 */
herr_t simulation_link_mesh(hid_t fid, hid_t sid, const char * meshname)
{
    herr_t status = DANU_FAILURE;
    char  mesh_link_name[] = SIM_MESH_LINK_NAME;
    char  mesh_root_name[] = MESH_ROOT_GROUP_NAME;
    char *mesh_fullname;
    size_t nchar;
    int exists;

    if ( H5_ISA_INVALID_ID(sid) ) {
        DANU_ERROR_MESS("Invalid simulation id");
        return status;
    }

    if ( H5_RETURN_OK(mesh_exists(fid,meshname,&exists) ) ) {

        if ( ! exists ) {
            DANU_WARN_MESS("Linkinng to a mesh that does not exist");
        }

        /* Build the full mesh name */
        nchar = 1;
        nchar = strlen(mesh_root_name);
        nchar++;
        nchar+= strlen(meshname);
        nchar++;
        mesh_fullname = DANU_MALLOC(char,nchar);

        sprintf(mesh_fullname,"/%s/%s",mesh_root_name,meshname);

        /* Create a soft link if the mesh does not exist, hard otherwise */
	if ( exists ) {
	  status = danu_link_create_hard(sid,mesh_link_name,fid,mesh_fullname);
        }
	else {
          status = danu_link_create_soft(sid,mesh_link_name,mesh_fullname);
	}

        DANU_FREE(mesh_fullname);

    }

    return status;
}

/*
 * Routine: hid_t simulation_open_mesh_link(hid_t sid)
 * Purpose: Open soft-link mesh group. 
 * Description: Open a soft linked mesh group under the simulation id sid
 *              Since the link is soft, only a warning
 *              will occur if the mesh does not exist.  
 *
 * Parameters:
 *           sid       IN     Simulation identifier
 *           mid       OUT    HDF5 identifier for the mesh
 *                              
 * Returns: An identifier (hid_t) is returned. A negative return value 
 *          indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid, or if the link fails.
 *         Since the link is a soft link, only a warning will be issued if meshname
 *         does not exist.
 *        
 */

hid_t simulation_open_mesh_link(hid_t sid)
{
  hid_t mid = H5I_BADID;
  char  mesh_link_name[] = SIM_MESH_LINK_NAME;

  if ( danu_group_exists(sid,mesh_link_name) ) {
    mid = danu_group_open(sid,mesh_link_name);
  }
  else {
    DANU_WARN_MESS("Link mesh group does not exist");
  }

  return mid;

}

/*
 * Routine: hid_t simulation_mesh_link_exists(hid_t sid)
 * Purpose: Check the existence of the soft link mesh group
 * Description: Check the existence of a soft linked mesh group
 *              under the simulation id sid.
 *
 * Parameters:
 *           sid       IN     Simulation identifier
 *           *exists   OUT    Flag TRUE (!=0) or FALSE(=0) if
 *                            mesh link group exists
 *                              
 * Returns: A flag (hid_t) is returned to indicate status.
 *           A negative return value indicates an error.
 *          
 * Errors: An error is raised if sid is not a valid.
 *        
 */
herr_t simulation_mesh_link_exists(hid_t sid, int *exists)
{

  herr_t status = DANU_FAILURE;
  char  mesh_link_name[] = SIM_MESH_LINK_NAME;
  hid_t gid;
  hbool_t flag;

  if ( H5_ISA_INVALID_ID(sid) ) {
     DANU_ERROR_MESS("Invalid simulation id");
     return status;
  }

  status = danu_link_exists(sid,mesh_link_name,&flag);
  if ( flag == TRUE ) {
    *exists=1;
  }
  else {
    *exists=0;
  }

  return status;

}








