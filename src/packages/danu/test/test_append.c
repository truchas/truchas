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

void print_array(int num, const hsize_t * array)
{
    int i;

    for(i=0;i<num;i++) {
        printf("%d %d\n", i, array[i]);
    }

}

int main(int argc, char ** argv)
{
    hid_t  fid,id;

    unsigned  i;

    int        dim = 2;
    int        ncell = 8;
    int        num;
    hsize_t   *isize;
    int       *idata;

    int       data3[3] = { -1, -1, -1 };
    hsize_t   size[2] = { 3, 1 }; 


    const char  test_file[] = "deleteme.h5";
    const char  data_name[] = "Extend this data";

    /* Get rid of old copies of the test file */
    if ( FILE_EXISTS(test_file) ) {
        danu_file_delete(test_file);
    }

    /* Generate the data */
    isize = DANU_MALLOC(hsize_t, dim);
    num = 1;
    for(i=0;i<dim;i++) { 
        isize[i] = ncell;
        num*=ncell;
    }

    idata = DANU_MALLOC(int,num);
    danu_rand_data_int(0,num,num,idata);

    /* Simple file create */
    fid = danu_file_create(test_file);

    /* Create the dataset with the extend flag set to TRUE */
    id = danu_dataset_create(fid,data_name,H5T_NATIVE_INT,dim,isize,TRUE);

    if ( H5_ISA_INVALID_ID(id) ) {
        DANU_ERROR_MESS("Failed to create dataset");
        goto FAIL_EXIT;
    }

    /* First write ..... */
    if ( danu_dataset_write(id,NULL,H5T_NATIVE_INT,dim,isize,idata) < 0 ) {
        DANU_ERROR_MESS("Failed the first data write");
        goto FAIL_EXIT;
    }

    /* Now the second ..... */
    if ( danu_dataset_append(id,H5T_NATIVE_INT,isize,idata) < 0 ) {
        DANU_ERROR_MESS("Failed to write the second dataset");
        goto FAIL_EXIT;
    }

    /* And the third .... */
    if ( danu_dataset_append(id,H5T_NATIVE_INT,size,data3) < 0 ) {
        DANU_ERROR_MESS("Failed to write the third dataset");
        goto FAIL_EXIT;
    }

    /* And try out the append_int call */
    if ( danu_data_append_int(fid,data_name,dim,size,data3) < 0 ) {
        DANU_ERROR_MESS("append_int call failed");
        goto FAIL_EXIT;
    }
    

    /* Close the objects */
    danu_dataset_close(id);
    danu_file_close(fid);

    /* Free Memory */
    DANU_FREE(isize);
    DANU_FREE(idata);

    printf("Test PASS\n");
    return 0;

    FAIL_EXIT:
             return DANU_FAIL;
}







