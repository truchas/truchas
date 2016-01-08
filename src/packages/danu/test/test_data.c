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
#define DATA_INT_SIZE 10

#define DATA_DOUBLE_NAME "My DOUBLE Data"
#define DATA_DOUBLE_DIM  2
#define DATA_DOUBLE_SIZE 10

#define DATA_CHAR_NAME "My CHAR Data"
#define DATA_CHAR_NUM  3
#define DATA_CHAR_LEN 128

int main(int argc, char ** argv)
{
    hid_t  fid;
    herr_t status;

    didx_t  i;

    hsize_t   *int_size, *d_size;
    dsize_t    num;
    int       *int_data, *int_ref;
    double    *d_data, *d_ref;

    char *strings[] = {"This is a string", "And another", "Ditto", "Hello World"};
    int   str_cnt = 4;

    /* Get rid of old copies of the test file */
    if ( FILE_EXISTS(TEST_FILE) ) {
        danu_file_delete(TEST_FILE);
    }

    /* Generate the data */
    int_size = DANU_MALLOC(hsize_t, DATA_INT_DIM);
    num = 1;
    for(i=0;i<DATA_INT_DIM;i++) { 
        int_size[i] = DATA_INT_SIZE;
        num*=DATA_INT_SIZE;
    }

    int_data = DANU_MALLOC(int,num);
    int_ref  = DANU_MALLOC(int,num);
    danu_rand_data_int(0,DATA_INT_SIZE,num,int_data);
    memcpy(int_ref,int_data,num*sizeof(int));

    d_size = DANU_MALLOC(hsize_t, DATA_DOUBLE_DIM);
    num = 1;
    for(i=0;i<DATA_DOUBLE_DIM;i++) { 
        d_size[i] = DATA_DOUBLE_SIZE;
        num*=DATA_DOUBLE_SIZE;
    }

    d_data = DANU_MALLOC(double,num);
    d_ref  = DANU_MALLOC(double,num);
    danu_rand_data_double(-100.0,100,num,d_data);
    memcpy(d_ref,d_data,num*sizeof(double));

    /* Simple file create */
    fid = danu_file_create(TEST_FILE);

    /* Write the data to disk */
    status = danu_data_write_int(fid,DATA_INT_NAME,DATA_INT_DIM,int_size,int_data);

    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to write integer data");
        goto FAIL_EXIT;
    }

    status = danu_data_write_double(fid,DATA_DOUBLE_NAME,DATA_DOUBLE_DIM,d_size,d_data);

    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to write double data");
        goto FAIL_EXIT;
    }

    status = danu_data_write_strings(fid,DATA_CHAR_NAME,str_cnt,strings);

    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to write char data");
        goto FAIL_EXIT;
    }

    /* Close the file */
    danu_file_close(fid);

    /* Zero out the int_data array */
    for(i=0;i<DATA_INT_DIM*DATA_INT_SIZE;i++)
        int_data[i] = 0;

    for(i=0;i<DATA_DOUBLE_DIM*DATA_DOUBLE_SIZE;i++)
        d_data[i] = 0.0;

    for(i=0;i<str_cnt;i++) 
        strings[i] = '\0';


    /*Open the file to read */
    fid = danu_file_open(TEST_FILE,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_RDONLY);
    
    status = danu_data_read_int(fid,DATA_INT_NAME,DATA_INT_DIM,int_size,int_data);
    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to open int data set for reading");
        goto FAIL_EXIT;
    }

    status = danu_data_read_double(fid,DATA_DOUBLE_NAME,DATA_DOUBLE_DIM,d_size,d_data);
    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to open double data set for reading");
        goto FAIL_EXIT;
    }

    status = danu_data_read_strings(fid,DATA_CHAR_NAME,num,strings);
    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to read string data");
        goto FAIL_EXIT;
    }

    /* Check the data */
    for(i=0;i<DATA_INT_DIM*DATA_INT_SIZE;i++) {
        if ( int_ref[i] != int_data[i] ) {
            danu_error_printf("MISMATCH FOUND index=%d IN=0x%lX OUT=0x%lX",i,int_ref[i],int_data[i]);
            goto FAIL_EXIT;
        }
    }

    for(i=0;i<DATA_DOUBLE_DIM*DATA_DOUBLE_SIZE;i++) {
        if ( d_ref[i] != d_data[i] ) {
            danu_error_printf("MISMATCH FOUND index=%d IN=0x%lX OUT=0x%lX",i,d_ref[i],d_data[i]);
            goto FAIL_EXIT;
        }
    }

    printf("DATA:\n");
    for(i=0;i<str_cnt;i++) {
        printf("<%s>\n",strings[i]);
    }

    /* Close the objects */
    danu_file_close(fid);

    /* Free Memory */
    DANU_FREE(int_size);
    DANU_FREE(int_ref);
    DANU_FREE(int_data);

    DANU_FREE(d_size);
    DANU_FREE(d_ref);
    DANU_FREE(d_data);

    printf("Test PASS\n");
    return DANU_SUCCESS;

    FAIL_EXIT:
             return DANU_FAILURE;
}







