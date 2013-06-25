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
#define MESH_FILE "mesh001.h5"

#define LINK_NAME  "Mesh"
#define OBJ_NAME   "Mesh001"

#define DATA_DOUBLE_NAME "My DOUBLE Data"
#define DATA_DOUBLE_DIM  2
#define DATA_DOUBLE_SIZE 10

int main(int argc, char ** argv)
{
    hid_t  fid, gid;
    hid_t  id;
    herr_t status;

    didx_t  i;

    dsize_t    num;
    hsize_t   *d_size;
    double    *d_data, *d_ref;
    char      link_name[128];

    ssize_t    onum,obj;
    hid_t      objects;
    H5O_info_t 

    /* Get rid of old copies of the test file */
    if ( FILE_EXISTS(TEST_FILE) ) {
        danu_file_delete(TEST_FILE);
    }

    if ( FILE_EXISTS(MESH_FILE) ) {
        danu_file_delete(MESH_FILE);
    }



    /* Generate the data */
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

    /* Create an external group in mesh file */
    gid = danu_group_create_external(fid,LINK_NAME,MESH_FILE,OBJ_NAME);
    if ( H5_ISA_INVALID_ID(gid) ) {
        printf("Return = %d LINK= %s\n", gid, LINK_NAME);
        DANU_ERROR_MESS("Failed to create the link");
        goto FAIL_EXIT;
    }
#if 0    
    id = danu_dataset_create_double(gid,DATA_DOUBLE_NAME,DATA_DOUBLE_DIM,d_size);
    printf("Create a DOUBLE returned %d\n",id);
    status = danu_dataset_write(id,NULL,H5T_NATIVE_DOUBLE,DATA_DOUBLE_DIM,d_size,d_data);
    printf("write DOUBLE returned %d\n",status);
    status = danu_dataset_close(id);
    printf("close returned %d\n",status);
#endif

    /* Create a dataset ... this data physically exists in the mesh file */
    status = danu_data_write_double(gid,DATA_DOUBLE_NAME,DATA_DOUBLE_DIM,d_size,d_data);
    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Falied to write double dataset");
        goto FAIL_EXIT;
    }

    /* Close the group and file */
    danu_group_close(gid);
    printf("After closing group number of objects open (global) %d\n", (int) H5Fget_obj_count(fid,H5F_OBJ_ALL));
    printf("After closing group number of objects open (local) %d\n", (int) H5Fget_obj_count(fid,H5F_OBJ_ALL|H5F_OBJ_LOCAL));
    danu_file_close(fid);

    printf("Close file %s\n", TEST_FILE);

    onum = H5Fget_obj_count(H5F_OBJ_ALL,H5F_OBJ_ALL);
    if ( onum > 0 ) {
        printf("After file close ... number of objects open %d\n", (int) onum);
        objects = DANU_MALLOC(hid_t,onum);
        for(obj=0;obj<onum;obj++) {
            printf("Object ID = %d NAME=%s\n", objects[obj], 

   // printf("CALLING HDF5 CLOSE! Return %d\n", H5close());

    /* Zero out the data to read */
    for(i=0;i<DATA_DOUBLE_DIM*DATA_DOUBLE_SIZE;i++)
        d_data[i] = 0.0;


    /*Open the file to read */
    fid = danu_file_open(TEST_FILE,DANU_FILE_ACT_OPEN,H5F_ACC_RDONLY);
    //fid = H5Fopen(TEST_FILE,H5F_ACC_RDONLY,H5P_DEFAULT);
    sprintf(link_name,"/%s",LINK_NAME);
    printf("Opening %s in %s\n", link_name,TEST_FILE);
    gid = danu_group_open(fid,LINK_NAME);
    if ( H5_ISA_INVALID_ID(gid) ) {
        DANU_ERROR_MESS("Failed to open link for reading");
        goto FAIL_EXIT;
    }

    
    id = danu_dataset_open(gid,DATA_DOUBLE_NAME);
    printf("Opening dataset returned %d\n", id);
    status = danu_dataset_read(id,NULL,H5T_NATIVE_DOUBLE,DATA_DOUBLE_DIM,d_size,d_data);
    printf("reading  dataset returned %d\n", status);
    status = danu_dataset_close(id);
    printf("closing  dataset returned %d\n", status);

#if 0
    status = danu_data_read_double(gid,DATA_DOUBLE_NAME,DATA_DOUBLE_DIM,d_size,d_data);
    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to open double data set for reading");
        goto FAIL_EXIT;
    }

    for(i=0;i<DATA_DOUBLE_DIM*DATA_DOUBLE_SIZE;i++) {
        if ( d_ref[i] != d_data[i] ) {
            danu_error_printf("MISMATCH FOUND index=%d IN=0x%lX OUT=0x%lX",i,d_ref[i],d_data[i]);
            goto FAIL_EXIT;
        }
    }
#endif

    /* Close the objects */
    danu_group_close(gid);
    danu_file_close(fid);

    /* Free Memory */
    DANU_FREE(d_size);
    DANU_FREE(d_ref);
    DANU_FREE(d_data);

    printf("Test PASS\n");
    return DANU_SUCCESS;

    FAIL_EXIT:
             return DANU_FAIL;
}







