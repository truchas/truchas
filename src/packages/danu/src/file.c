/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * file.c
 *
 *  DANU Files
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
#include <danu_utils.h>
#include <danu_memory.h>
#include <danu_types.h>

#include <danu_file.h>

/* Global defines for cache (chunking) */
#define FILE_CACHE_SIZE 1048576
#define FILE_CACHE_NUM  50
#define FILE_CACHE_W0   0.75

/* Global defines for on-disk data types and data order */
#define FILE_DATA_ORDER H5T_ORDER_LE

/*
 * Routine: hid_t danu_file_create(const char *name)
 * Purpose: Create an HDF5 file
 * Description: Creates an HDF5 file named name. If file exists, it will be 
 *              clobbered. There will be a warning message if this occurs.
 *              Checks the input and will return immediately if error is
 *              detected.
 *
 * Parameters:
 *           name               IN              Name of the file
 *
 * Returns: Returns a HDF5 identifier for the new file.
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately. Returns a neagtive value if an error occurs.
 *
 */
hid_t danu_file_create(const char * name)
{
     hid_t fid;
     hid_t apl;

     /* Check input */
     if ( DANU_BAD_PTR(name) ) {
         DANU_ERROR_MESS("Invalid string pointer argument");
         return H5I_INVALID_HID;
     }

     if ( DANU_EMPTY_STRING(name)) {
         DANU_ERROR_MESS("Invalid empty string argument");
         return H5I_INVALID_HID;
     }

     /* Check if the file exists */
     if (  FILE_EXISTS(name) ) {
         danu_warn_printf("File %s exists ... will clobber", name);
     }

     /* Create the access property list 
      *  enable caching 
     */
     apl = H5Pcreate(H5P_FILE_ACCESS);
     H5Pset_cache(apl,0,FILE_CACHE_NUM,FILE_CACHE_SIZE,FILE_CACHE_W0);

     /* Create the file */
     fid = H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,apl);

     /* Close the access property list */
     H5Pclose(apl);

     return fid;
}

/*
 * Routine: hid_t danu_file_open(const char *name, unsigned action, unsigned access)
 * Purpose: Open an HDF5 file
 * Description: Open an HDF5 file named name. If file exists, routine checks
 *              if it is an HDf5 file. An error is raised if it is not an HDF5
 *              file. The following table describes the behavior for each combination
 *              of access or action, and when errors are raised. The CREATE access
 *              option ALWAYS clobbers any existing file. 
 *              
 *              X = NOT POSSIBLE
 *
 *                 ACCESS
 *              RDONLY          RDWR        CLOBBER         APPEND     IF_FILE_EXISTS
 *    ACTION    
 *    OPEN      OK              OK          OK              OK         File clobbered with CLOBBER
 *
 *    CREATE    X               OK          OK              X          File is always clobbered
 *
 * Parameters:
 *           name               IN              Name of the file
 *           action             IN              Action flag (DANU_FILE_ACT_CREATE|OPEN)
 *           access             IN              Access flag (HDF5 Access flag)
 *
 * Returns: Returns a HDF5 identifier for the file.
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately. Also checks if an existing file is an HDF5 file. 
 *         will Returns a neagtive value if an error occurs.
 *
 */
 hid_t danu_file_open(const char * name, unsigned access, unsigned action)
 {
     hid_t fid, apl;
     htri_t  file_status;


     /* Check input */
     if ( DANU_BAD_PTR(name) ) {
         DANU_ERROR_MESS("Invalid string pointer argument");
         return H5I_INVALID_HID;
     }

     if ( DANU_EMPTY_STRING(name) ) {
         DANU_ERROR_MESS("Invalid string name argument");
         return H5I_INVALID_HID;
     }

     /* Check if the file exists and is an HDF5 file. Save the value to
     *  avoid another query
     *  Three possible values
     *   1 = File exists and is HDF5 file
     *   0 = File does exists but is not HDF5 file
     *   <0 = File does not exist
     */
     file_status = H5Fis_hdf5(name);

     if ( action == DANU_FILE_ACT_CREATE ) {
         fid = danu_file_create(name);
     }
     else if ( action == DANU_FILE_ACT_OPEN ) {

         if ( file_status > 0 ) {

             /* Create the file access parameter list ... enable caching */
             apl = H5Pcreate(H5P_FILE_ACCESS);
             H5Pset_cache(apl,0, FILE_CACHE_NUM,FILE_CACHE_SIZE,FILE_CACHE_W0);

             fid = H5Fopen(name,access,apl);

             /* Release the access list once file is open */
             H5Pclose(apl);

         }
         else {
             fid = H5I_INVALID_HID;
             danu_error_printf("H5Fis_hdf5 returned %d\n", file_status);
             if ( file_status == 0 ) {
                 danu_error_printf("File %s is not an HDF5 file\n",name);
             }
             else {
                 danu_error_printf("File %s does not exist. Can not create file with open action.\n",name);
             }
         }

     }
     else {
         fid = H5I_INVALID_HID;
         danu_error_printf("Unknown file action (%d)",action);
     }

     return fid;
 }

/*
 * Routine: herr_t danu_file_flush(hid_t fid)
 * Purpose: Locally (only the file fid) flush an opened HDF5 file.
 * Description: Flush all the internal HDF5 buffers associated with this
 *              file. Will request that the OS flush the buffers and write to disk.
 *              Check the input for a valid id. Will return immediately if fid is
 *              not valid. 
 *              
 *
 * Parameters:
 *           fid               IN              HDF5 file identifier
 *
 * Returns: Returns the return code from H5Fflush
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately.
 *
 */
herr_t danu_file_flush(hid_t fid)
{
     /* Check input */
     if (! H5_ISA_VALID_ID(fid) ) {
         DANU_ERROR_MESS("Invalif HDf5 file identifier argument");
         return DANU_FAILURE;
     }

     return H5Fflush(fid,H5F_SCOPE_LOCAL);
}
/*
 * Routine: herr_t danu_file_flush_global(hid_t fid)
 * Purpose: Global (file fid and all buffers for the entire virtual file) flush an opened HDF5 file.
 * Description: Flush all the internal HDF5 buffers associated with this
 *              file and the any buffers associated with the virtual file.
 *              Will request that the OS flush the buffers and write to disk.
 *              Check the input for a valid id. Will return immediately if fid is
 *              not valid. 
 *              
 *
 * Parameters:
 *           fid               IN              HDF5 file identifier
 *
 * Returns: Returns the return code from H5Fflush
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately.
 *
 */
herr_t danu_file_flush_global(hid_t fid)
{
     /* Check input */
     if (! H5_ISA_VALID_ID(fid) ) {
         DANU_ERROR_MESS("Invalif HDf5 file identifier argument");
         return DANU_FAILURE;
     }

     return H5Fflush(fid,H5F_SCOPE_GLOBAL);
}
/*
 * Routine: herr_t danu_file_close(hid_t fid)
 * Purpose: Close an HDF5 file.
 * Description: Close file fid and all objects associated with the file.
 *              Once this routine is called any dataspace, link or group 
 *              associated with this file will also be closed
 *
 * Parameters:
 *           fid               IN              HDF5 file identifier
 *
 * Returns: Returns the return code from H5close
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately.
 *
 */
herr_t danu_file_close(hid_t fid)
{
    ssize_t  num_obj;
    hid_t *  objects;
    didx_t   i;


    /* Check input */
    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("Invalid HDF5 file identifier argument");
        return DANU_FAILURE;
    }

    if (! H5_ISA_FILE_ID(fid) ) {
        DANU_ERROR_MESS("HDF5 identifer is not a file id");
        return DANU_FAILURE;
    }

    /* Flush the all the buffers. This will speed up the close calls. */
    danu_file_flush_global(fid);

    /* Determine the number of objects */
    num_obj = H5Fget_obj_count(fid,H5F_OBJ_ALL);
    objects = DANU_MALLOC(hid_t,num_obj);

    /* Close all the attributes */
    num_obj = H5Fget_obj_count(fid,H5F_OBJ_ATTR);
    if ( num_obj > 0 ) {
        H5Fget_obj_ids(fid,H5F_OBJ_ATTR,num_obj,objects);
        for(i=0;i<num_obj;i++)
            H5Aclose(objects[i]);
    }

    /* Close all the datatypes */
    num_obj = H5Fget_obj_count(fid,H5F_OBJ_DATATYPE);
    if ( num_obj > 0 ) {
        H5Fget_obj_ids(fid,H5F_OBJ_DATATYPE,num_obj,objects);
        for(i=0;i<num_obj;i++)
            H5Tclose(objects[i]);
    }

    /* Close all the datasets */
    num_obj = H5Fget_obj_count(fid,H5F_OBJ_DATASET);
    if ( num_obj > 0 ) {
        H5Fget_obj_ids(fid,H5F_OBJ_DATASET,num_obj,objects);
        for(i=0;i<num_obj;i++)
            H5Dclose(objects[i]);
    }

    /* Close all the groups */
    num_obj = H5Fget_obj_count(fid,H5F_OBJ_GROUP);
    if ( num_obj > 0 ) {
        H5Fget_obj_ids(fid,H5F_OBJ_GROUP,num_obj,objects);
        for(i=0;i<num_obj;i++)
            H5Gclose(objects[i]);
    }

    /* Close all the files, this call will return the current file id */
    num_obj = H5Fget_obj_count(fid,H5F_OBJ_FILE);
    if ( num_obj > 1 ) {
        H5Fget_obj_ids(fid,H5F_OBJ_FILE,num_obj,objects);
        for(i=0;i<num_obj;i++) {
            if ( objects[i] != fid ) H5Fclose(objects[i]);
        }
    }

    /* Free memory */
    DANU_FREE(objects);

    return 0;
    return H5Fclose(fid);
}
/*
 * Routine: hid_t danu_file_datatype(hid_t mem_type)
 * Purpose: Define the on-disk data type for memory type mem_type
 * Description: All on-disk data is little endian. This routine provides
 *              a uniform data type order for all data on a file. The calling
 *              function is responsible for closing this data type to control
 *              resource usage.
 *
 * Parameters:
 *           mem_type               IN              HDF5 memory type
 *
 * Returns: Returns a HDF5 data type with order little endian
 * Errors: The input is checked and if an error is deteced in the input the routine
 *         returns immediately.
 *
 */
 hid_t danu_file_dataype(hid_t memtype)
 {
     hid_t dt;

     if ( H5_ISA_VALID_ID(memtype) ) {
         dt = H5Tcopy(memtype);
         H5Tset_order(dt,FILE_DATA_ORDER);
     }
     else {
         DANU_ERROR_MESS("Invalid HDF5 memory identifier");
         dt = H5I_INVALID_HID;
     }

     return dt;
 }

hid_t danu_file_open_rdonly(const char *name)
{
    return danu_file_open(name,DANU_FILE_ACC_RDONLY,DANU_FILE_ACT_OPEN);
}

hid_t danu_file_open_append(const char *name)
{
    return danu_file_open(name,DANU_FILE_ACC_APPEND,DANU_FILE_ACT_OPEN);
}

hid_t danu_file_open_rdwr(const char *name)
{
    return danu_file_open(name,DANU_FILE_ACC_RDWR,DANU_FILE_ACT_OPEN);
}
