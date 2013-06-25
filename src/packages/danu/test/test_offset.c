/*
*
*  Sample File
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

#define GROUP_NAME "Mesh/Mesh1/Coordinates"

#define DATA_NAME "Nodal Coordinates"
#define DIM       3
#define COORD_IDX 1
#define NCELLS    10000

int main(int argc, char ** argv)
{
    hid_t  fid, gid;
    hid_t  dsid;
    herr_t status;

    didx_t  i;
    hsize_t size[2];
    hsize_t offset[2], coord_size[2];
    dslab_t *slab;
    double   *coord,*ref;

    /* Get rid of old copies of the test file */
    if ( FILE_EXISTS(TEST_FILE) ) {
        danu_file_delete(TEST_FILE);
    }


    /* Create the memory sapce and generate the data */
    coord    = DANU_MALLOC(double,NCELLS);
    ref      = DANU_MALLOC(double,NCELLS);
    danu_rand_data_double(0.0,10.0,NCELLS,coord);
    memcpy(ref,coord,sizeof(double)*NCELLS);
    //for(i=0;i<NCELLS;i++)
    //    printf("i=%d COORD=%1.5f\n", i, coord[i]);

    /* Simple file and group create */
    fid = danu_file_create(TEST_FILE);
    gid = danu_group_create(fid,GROUP_NAME,FALSE);

    /* Create the file data space */
    size[0] = 3;
    size[1] = NCELLS;
    dsid = danu_dataset_create(gid,DATA_NAME,H5T_IEEE_F64LE,2,size,FALSE,FALSE);


    /* Create the Hyperslab */
    slab = danu_slab_alloc(dsid);
    offset[0] = COORD_IDX;
    offset[1] = 0;
    coord_size[0] = 1;
    coord_size[1] = NCELLS;
    danu_slab_contiguous(slab,2,offset,coord_size);

    /* Now Write the coordinate data */
    danu_dataset_write(dsid,slab,H5T_NATIVE_DOUBLE,2,coord_size,coord);


    /* Close the file */
    status = H5Dclose(dsid);
    status = danu_group_close(gid);
    status = danu_file_close(fid);

    /* Zero out the coordinate array */
    for(i=0;i<NCELLS;i++)
        coord[i] = 0.0;

    /*Open the file to read */
    fid = danu_file_open(TEST_FILE,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_RDONLY);
    gid = danu_group_open(fid,GROUP_NAME);

    
    dsid = danu_dataset_open(gid,DATA_NAME);
    if ( H5_ISA_INVALID_ID(dsid) ) {
        DANU_ERROR_MESS("Failed to open data set for reading");
        goto FAIL_EXIT;
    }

    if ( DANU_RETURN_FAIL(danu_dataset_read(dsid,slab,H5T_NATIVE_DOUBLE,2,coord_size,coord)) ) {
        DANU_ERROR_MESS("Failed to read data set");
        goto FAIL_EXIT;
    }

    /* Check the data */
    for(i=0;i<NCELLS;i++) {
        if ( ref[i] != coord[i] ) {
            danu_error_printf("MISMATCH FOUND index=%d IN=0x%lX OUT=0x%lX",i,ref[i],coord[i]);
            goto FAIL_EXIT;
        }
    }


    /* Close the objects */
    danu_dataset_close(dsid);
    danu_group_close(gid);
    danu_file_close(fid);


    /* Free Memory */
    DANU_FREE(ref);
    DANU_FREE(coord);


    return DANU_SUCCESS;


    FAIL_EXIT:
             return DANU_FAIL;
}







