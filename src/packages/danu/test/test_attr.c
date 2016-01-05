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
#include <danu_group.h>
#include <danu_attribute.h>

#define TEST_FILE "deleteme.h5"

#define INT_ATTR_NAME     "Integer Attribute"
#define INT_ATTR_VALUE     0x1234ABDC
#define DOUBLE_ATTR_NAME  "DOUBLE Attribute"
#define DOUBLE_ATTR_VALUE  1.23456789e-10
#define STR_ATTR_NAME     "STRING Attribute"
#define STR_ATTR_VALUE     "This is a string attribute"
#define STR_ATTR_LEN       128

int main(int argc, char ** argv)
{
    hid_t fid,gid;
    double dtest;
    int    itest;
    char   stest[STR_ATTR_LEN]; 
    int    i, num_found;

    char  **names;
    size_t  *size;
    char    group_name[] = "Group"; 

    /* Get rid of old copies of the test file */
    if ( FILE_EXISTS(TEST_FILE) ) {
        danu_file_delete(TEST_FILE);
    }

    /* Simple create */
    fid = danu_file_create(TEST_FILE);
    gid = danu_group_create(fid,group_name);

    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("Failed to create a test file");
        goto FAIL_EXIT;
    }

    /* Now Write the attributes */
    if ( DANU_RETURN_FAIL(danu_attr_write_int(fid, INT_ATTR_NAME, INT_ATTR_VALUE) ) ){
        DANU_ERROR_MESS("Failed to write int attribute");
        goto FAIL_EXIT;
    }

    if ( DANU_RETURN_FAIL(danu_attr_write_double(fid, DOUBLE_ATTR_NAME, DOUBLE_ATTR_VALUE) ) ){
        DANU_ERROR_MESS("Failed to write double attribute");
        goto FAIL_EXIT;
    }

    
    if ( DANU_RETURN_FAIL(danu_attr_write_string(fid, STR_ATTR_NAME, STR_ATTR_VALUE) ) ){
        DANU_ERROR_MESS("Failed to write char attribute");
        goto FAIL_EXIT;
    }

    danu_file_close(fid); 

    /* Re-open file and read the attributes back */
    fid = danu_file_open(TEST_FILE,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_RDONLY);

    /* Int attribute */
    if ( DANU_RETURN_FAIL(danu_attr_read_int(fid, INT_ATTR_NAME, &itest) ) ){
        DANU_ERROR_MESS("Failed to read int attribute");
        goto FAIL_EXIT;
    }

    /* Check the value */
    printf("INTEGER IN=%d OUT=%d\n", INT_ATTR_VALUE, itest);
    if ( itest != INT_ATTR_VALUE ) {
        danu_error_printf("Failed to read the correct value!\n");
        goto FAIL_EXIT;
    }
    
    /* Double attribute */
    if ( DANU_RETURN_FAIL(danu_attr_read_double(fid, DOUBLE_ATTR_NAME, &dtest) ) ){
        DANU_ERROR_MESS("Failed to read double attribute");
        goto FAIL_EXIT;
    }

    /* Check the value */
    printf("Double IN=%1.9e OUT=%1.9e\n", DOUBLE_ATTR_VALUE, dtest);
    if ( dtest != DOUBLE_ATTR_VALUE ) {
        danu_error_printf("Failed to read the correct value!\n");
        goto FAIL_EXIT;
    }


    /* String attribute */
    if ( DANU_RETURN_FAIL(danu_attr_read_string(fid, STR_ATTR_NAME, stest,STR_ATTR_LEN) ) ){
        DANU_ERROR_MESS("Failed to read string attribute");
        goto FAIL_EXIT;
    }

    /* Check the value */
    printf("string IN='%s' OUT='%s'\n", STR_ATTR_VALUE, stest);
    if ( strncmp(STR_ATTR_VALUE,stest,STR_ATTR_LEN) != 0 ) {
        danu_error_printf("Failed to read the correct value!\n");
        goto FAIL_EXIT;
    }

    /* Find the count n the file*/
    if ( DANU_RETURN_FAIL(danu_attr_count(fid,&num_found)) ) {
        DANU_ERROR_MESS("Failed to count the number of attributes");
        goto FAIL_EXIT;
    }
    printf("Found %d attributes in file\n",num_found);

    names = DANU_MALLOC(char *, num_found);
    size  = DANU_MALLOC(size_t, num_found);
    for(i=0;i<num_found;i++) {
        size[i] = 128;
        names[i] = DANU_MALLOC(char,128);
    }

    /* Grab the names */
    if ( DANU_RETURN_FAIL(danu_attr_names(fid,num_found,size,names)) ) {
        DANU_ERROR_MESS("Failed to read the attribute names");
        goto FAIL_EXIT;
    }
    else {
        printf("Found the attributes\n");
        for(i=0; i<num_found; i++) {
            printf("I=%d\tName=%s\n",i,names[i]);
        }
    }

    /* Find the count in the group*/
    gid = danu_group_open(fid,group_name);
    if ( DANU_RETURN_FAIL(danu_attr_count(gid,&num_found)) ) {
        DANU_ERROR_MESS("Failed to count the number of attributes");
        goto FAIL_EXIT;
    }
    printf("Found %d attributes in group\n",num_found);


    /* Free memory */
    DANU_FREE(size);
    for(i=0;i<3;i++) {
        DANU_FREE(names[i]);
    }
    DANU_FREE(names);


 
    goto SUCCESS_EXIT;



    FAIL_EXIT:
             DANU_ERROR_MESS("Test FAIL");
             return DANU_FAILURE;
    SUCCESS_EXIT:
             //danu_file_delete(TEST_FILE);
             return DANU_SUCCESS;
}







