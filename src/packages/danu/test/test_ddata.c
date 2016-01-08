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

#define DATA_DOUBLE_NAME "My DOUBLE Data"
#define DATA_DOUBLE_DIM  3
#define DATA_DOUBLE_SIZE 100

#define DATA_ATTR_NAME "DUMMY ATTRIBUTE"
#define DATA_ATTR_VALUE 0xFABC1234

#define STR_ATTR_LEN 128

int main(int argc, char ** argv)
{
    hid_t  fid, gid;
    hid_t  data;

    hsize_t   *size;
    dsize_t    i,num;
    double    *double_data, *ref;
    int        value,vref;

    /* Get rid of old copies of the test file */
    if ( FILE_EXISTS(TEST_FILE) ) {
        danu_file_delete(TEST_FILE);
    }

    /* Generate the data */
    size = DANU_MALLOC(hsize_t, DATA_DOUBLE_DIM);
    num = 1;
    for(i=0;i<DATA_DOUBLE_DIM;i++) { 
        size[i] = DATA_DOUBLE_SIZE;
        num*=DATA_DOUBLE_SIZE;
    }

    double_data = DANU_MALLOC(double,num);
    ref      = DANU_MALLOC(double,num);
    danu_rand_data_double(-5.0,5.0,num,double_data);
    memcpy(ref,double_data,num*sizeof(double));

    /* Simple file and group create */
    fid = danu_file_create(TEST_FILE);
    gid = danu_group_create(fid,GROUP_NAME,FALSE);

    data = danu_dataset_create(gid,DATA_DOUBLE_NAME,H5T_IEEE_F64LE,DATA_DOUBLE_DIM,size,FALSE,FALSE);

    if ( H5_ISA_INVALID_ID(data) ) {
        DANU_ERROR_MESS("Failed to create the dataspace");
        goto FAIL_EXIT;
    }

    /* Write the data */
    if ( DANU_RETURN_FAIL(danu_dataset_write(data,NULL,H5T_NATIVE_DOUBLE,DATA_DOUBLE_DIM,size,double_data) ) ) {
        DANU_ERROR_MESS("Failed to write data to file");
        goto FAIL_EXIT;
    }

    /* Write data attribute */
    value = DATA_ATTR_VALUE; 
    danu_attr_write_int(data,DATA_ATTR_NAME,value);

    /* Close the file */
    danu_file_close(fid);

    /* Zero out the double_data array */
    for(i=0;i<num;i++)
        double_data[i] = 0;

    /*Open the file to read */
    fid = danu_file_open(TEST_FILE,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_RDONLY);
    gid = danu_group_open(fid,GROUP_NAME);

    
    data = danu_dataset_open(gid,DATA_DOUBLE_NAME);
    if ( H5_ISA_INVALID_ID(data) ) {
        DANU_ERROR_MESS("Failed to open data set for reading");
        goto FAIL_EXIT;
    }

    if ( DANU_RETURN_FAIL(danu_dataset_read(data,NULL,H5T_NATIVE_DOUBLE,DATA_DOUBLE_DIM,size,double_data)) ) {
        DANU_ERROR_MESS("Failed to read data set");
        goto FAIL_EXIT;
    }

    /* Check the data */
    for(i=0;i<num;i++) {
        if ( ref[i] != double_data[i] ) {
            danu_error_printf("MISMATCH FOUND index=%d IN=0x%lX OUT=0x%lX",i,ref[i],double_data[i]);
            goto FAIL_EXIT;
        }
    }

    /* Read the attribute and check it */
    vref = value;
    value = 0;
    danu_attr_read_int(data,DATA_ATTR_NAME,&value);
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
    DANU_FREE(double_data);

    return DANU_SUCCESS;

    FAIL_EXIT:
             printf("Test %s FAILED\n", argv[0]); 
             return DANU_FAIL;
}







