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


#define TEST_FILE "deleteme.h5"

#define GROUP_1_NAME "GroupA"
#define GROUP_2_NAME "Group B"
#define GROUP_3_NAME "GROUP/C"

#define GRP_ATTR1_NAME "GRP ATTRIBUTE1"
#define GRP_ATTR1_VALUE 0xf1234

#define GRP_ATTR2_NAME "GRP ATTRIBUTE2"
#define GRP_ATTR2_VALUE "This is a group with no data"

#define STR_ATTR_LEN 128

int main(int argc, char ** argv)
{
    hid_t fid, gid;
    int    itest,dum;
    char   stest[STR_ATTR_LEN]; 

    /* Get rid of old copies of the test file */
    if ( FILE_EXISTS(TEST_FILE) ) {
        danu_file_delete(TEST_FILE);
    }

    /* Simple create */
    fid = danu_file_create(TEST_FILE);

    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("Failed to create a test file");
        goto FAIL_EXIT;
    }

    /* Create each test group */
    if DANU_RETURN_FAIL(danu_group_create(fid,GROUP_1_NAME,FALSE) ) {
        DANU_ERROR_MESS("Failed to create group 1");
        goto FAIL_EXIT;
    }

    /* Create each test group */
    if DANU_RETURN_FAIL(danu_group_create(fid,GROUP_2_NAME,FALSE) ) {
        DANU_ERROR_MESS("Failed to create group 2");
        goto FAIL_EXIT;
    }

    /* Create each test group */
    gid = danu_group_create(fid,GROUP_3_NAME,FALSE);
    if ( H5_ISA_INVALID_ID(gid) ) {
        DANU_ERROR_MESS("Failed to create group 3");
        goto FAIL_EXIT;
    }


    /* Write attributes to the group 3 */
    if ( DANU_RETURN_FAIL(danu_attr_write_int(gid,GRP_ATTR1_NAME,GRP_ATTR1_VALUE)) ){
        danu_error_printf("Failed to write an attribute to group3");
        goto FAIL_EXIT;
    }

    /* Write attributes to the group 3 */
    if ( DANU_RETURN_FAIL(danu_attr_write_string(gid,GRP_ATTR2_NAME,GRP_ATTR2_VALUE)) ) {
        danu_error_printf("Failed to write an attribute to group3");
        goto FAIL_EXIT;
    }

    danu_file_close(fid); 

    /* Re-open file and open groups */
    fid = danu_file_open(TEST_FILE,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_RDONLY);

    /* Open group 1 */
    if ( DANU_RETURN_FAIL(danu_group_open(fid, GROUP_1_NAME) ) ){
        DANU_ERROR_MESS("Failed to open group 1");
        goto FAIL_EXIT;
    }

    /* Open group 2 */
    if ( DANU_RETURN_FAIL(danu_group_open(fid, GROUP_2_NAME) ) ){
        DANU_ERROR_MESS("Failed to open group 2");
        goto FAIL_EXIT;
    }

    /* Open group 3 */
    gid = danu_group_open(fid,GROUP_3_NAME);
    if ( H5_ISA_INVALID_ID(gid) ){
        DANU_ERROR_MESS("Failed to open group 3");
        goto FAIL_EXIT;
    }
    
    /* String attribute */
    if ( DANU_RETURN_FAIL(danu_attr_read_string(gid, GRP_ATTR2_NAME, stest,STR_ATTR_LEN) ) ){
        DANU_ERROR_MESS("Failed to read string attribute from group 3");
        goto FAIL_EXIT;
    }

    /* Check the value */
    printf("string IN='%s' OUT='%s'\n", GRP_ATTR2_VALUE, stest);
    if ( strncmp(GRP_ATTR2_VALUE,stest,STR_ATTR_LEN) != 0 ) {
        danu_error_printf("Failed to read the correct value!\n");
        goto FAIL_EXIT;
    }

    /* Integer attribute */
    itest = 0x0;
    if ( DANU_RETURN_FAIL(danu_attr_read_int(gid, GRP_ATTR1_NAME, &itest) ) ){
        DANU_ERROR_MESS("Failed to read int attribute from group 3");
        goto FAIL_EXIT;
    }

    /* Check the value */
    printf("Integer IN=%lx and OUT=%lx\n", GRP_ATTR1_VALUE, itest);
    if ( GRP_ATTR1_VALUE != itest ) {
        danu_error_printf("Failed to read the correct value!\n");
        goto FAIL_EXIT;
    }
 
    goto SUCCESS_EXIT;



    FAIL_EXIT:
             printf("Test FAILED\n");
             fflush(stdout);
             return DANU_FAIL;
    SUCCESS_EXIT:
             //danu_file_delete(TEST_FILE);
             return DANU_SUCCESS;
}







