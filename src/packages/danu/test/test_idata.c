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
#include <danu_h5_error.h>
#include <danu_memory.h>
#include <danu_file.h>
#include <danu_attribute.h>
#include <danu_group.h>

#include <danu_dataset.h>


#define TEST_FILE "deleteme.h5"

#define GROUP_NAME "DataGoesHere"

#define DATA_INT_NAME "My INT Data"
#define DATA_INT_DIM  1
#define DATA_INT_SIZE 10000

#define DATA_ATTR_NAME "DUMMY ATTRIBUTE"
#define DATA_ATTR_VALUE 0.98765432110985763829

#define STR_ATTR_LEN 128

int main(int argc, char ** argv)
{
    hid_t  fid, gid;
    hid_t  data;

    didx_t  i;

    hsize_t   *size;
    dsize_t    num;
    int      *int_data, *ref;
    double    value,vref;

    /* Get rid of old copies of the test file */
    if ( FILE_EXISTS(TEST_FILE) ) {
        danu_file_delete(TEST_FILE);
    }

    /* Generate the data */
    size = DANU_MALLOC(hsize_t, DATA_INT_DIM);
    num = 1;
    for(i=0;i<DATA_INT_DIM;i++) { 
        size[i] = DATA_INT_SIZE;
        num*=DATA_INT_SIZE;
    }

    int_data = DANU_MALLOC(int,num);
    ref      = DANU_MALLOC(int,num);
    danu_rand_data_int(0,DATA_INT_SIZE,num,int_data);
    memcpy(ref,int_data,num*sizeof(int));

    /* Simple file and group create */
    fid = danu_file_create(TEST_FILE);
    gid = danu_group_create(fid,GROUP_NAME);

    data = danu_dataset_create(gid,DATA_INT_NAME,H5T_STD_I32LE,DATA_INT_DIM,size,FALSE);

    if ( H5_ISA_INVALID_ID(data) ) {
        DANU_ERROR_MESS("Failed to create the dataspace");
        goto FAIL_EXIT;
    }

    /* Write the data */
    if ( DANU_RETURN_FAIL(danu_dataset_write(data,NULL,H5T_NATIVE_INT,DATA_INT_DIM,size,int_data) ) ) {
        DANU_ERROR_MESS("Failed to write data to file");
        goto FAIL_EXIT;
    }

    /* Write data attribute */
    value = DATA_ATTR_VALUE; 
    danu_attr_write_double(data,DATA_ATTR_NAME,value);

    /* Close the file */
    danu_file_close(fid);

    /* Zero out the int_data array */
    for(i=0;i<DATA_INT_DIM*DATA_INT_SIZE;i++)
        int_data[i] = 0;

    /*Open the file to read */
    fid = danu_file_open(TEST_FILE,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_RDONLY);
    gid = danu_group_open(fid,GROUP_NAME);

    
    data = danu_dataset_open(gid,DATA_INT_NAME);
    if ( H5_ISA_INVALID_ID(data) ) {
        DANU_ERROR_MESS("Failed to open data set for reading");
        goto FAIL_EXIT;
    }

    if ( DANU_RETURN_FAIL(danu_dataset_read(data,NULL,H5T_NATIVE_INT,DATA_INT_DIM,size,int_data)) ) {
        DANU_ERROR_MESS("Failed to read data set");
        goto FAIL_EXIT;
    }

    /* Check the data */
    for(i=0;i<num;i++) {
        if ( ref[i] != int_data[i] ) {
            danu_error_printf("MISMATCH FOUND index=%d IN=0x%lX OUT=0x%lX",i,ref[i],int_data[i]);
            goto FAIL_EXIT;
        }
    }

    /* Read the attribute and check it */
    vref = value;
    value = 0.0;
    danu_attr_read_double(data,DATA_ATTR_NAME,&value);
    if ( vref != value ) {
        DANU_ERROR_MESS("Failed to read attribute");
        goto FAIL_EXIT;
    }

    /* Close the objects */
    danu_dataset_close(data);
    danu_group_close(gid);
    danu_file_close(fid);

    /* Free Memory */
    DANU_FREE(size);
    DANU_FREE(ref);
    DANU_FREE(int_data);

    goto SUCCESS_EXIT;

    FAIL_EXIT:
             return DANU_FAIL;
    SUCCESS_EXIT:
             printf("Test PASS\n");
             return DANU_SUCCESS;
}







