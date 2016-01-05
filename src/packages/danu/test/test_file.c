/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
/*
*
*  Sample File
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_utils.h>
#include <danu_error.h>
#include <danu_memory.h>
#include <danu_file.h>

#include <stdio.h>


int main(int argc, char ** argv)
{

    derr_t status;
    didx_t i;
    hid_t fid;

    const char test_file[] = "deleteme.h5";

    /* Get rid of old copies of the test file */
    if ( FILE_EXISTS(test_file) ) {
        printf("Test file exists will delete\n");
        danu_file_delete(test_file);
    }

    /* Simple create */
    printf("Creating test file:%s\n", test_file);

    fid = danu_file_create(test_file);

    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("Failed to create a test file");
        goto FAIL_EXIT;
    }

    /* Simple close */
    if ( DANU_RETURN_FAIL(danu_file_close(fid)) ) {
        DANU_ERROR_MESS("Failed to close test file");
        goto FAIL_EXIT;
    }

    /* Open existing file for readonly */
    fid = danu_file_open(test_file,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_RDONLY);
    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("Failed to open a test file in read only mode");
        goto FAIL_EXIT;
    }
    danu_file_close(fid); 

    /* Open existing file for read write */
    fid = danu_file_open(test_file,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_RDWR);
    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("Failed to open a test file in read/write mode");
        goto FAIL_EXIT;
    }
    danu_file_close(fid); 

    
    /* Open existing file for append */
    fid = danu_file_open(test_file,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_APPEND);
    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("Failed to open a test file in append mode");
        goto FAIL_EXIT;
    }
    danu_file_close(fid); 



    return 0;

    FAIL_EXIT:
             return 1;
}







