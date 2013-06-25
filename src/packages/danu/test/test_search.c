#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include <hdf5.h>

#include <danu_utils.h>
#include <danu_memory.h>
#include <danu_file.h>
#include <danu_group.h>
#include <danu_dataset.h>

#include <danu_h5_object.h>


#define FILE "deleteme.h5"

#define GRPA "GRPA"
#define GRPB "GRPB"
#define GRPC "GRPC"
#define GRPD "GRPD"

#define DATASET_NAME  "Data"
#define DATASET_DIM   2
#define DATASET_SIZE  2
#define DATASET_VALUE 0xABCD0123


herr_t create_test_file(void);


herr_t create_test_file(void)
{
    herr_t status;

    hid_t fid, grp_a, grp_b, grp_c, grp_d;

    hsize_t *size;
    dsize_t  num;
    int      i,dim;
    int     *data;

    int  grp_num;
    char grp_name[16];
    int  ds_num;
    char ds_name[16];


    fid = danu_file_create(FILE);

    grp_a = danu_group_create(fid,GRPA);
    grp_b = danu_group_create(fid,GRPB);
    grp_c = danu_group_create(grp_b,GRPC);
    grp_d = danu_group_create(fid,GRPD);

    /* Generate a dataset */
    dim = DATASET_DIM;
    size = DANU_MALLOC(hsize_t,dim);
    num = 1;
    for(i=0;i<dim;i++) {
        num*=DATASET_SIZE;
        size[i] = DATASET_SIZE;
    }

    data = DANU_MALLOC(int,num);

    for(i=0;i<num;i++)
        data[i] = DATASET_VALUE;

    /* Write the data */ 
    status = danu_data_write_int(grp_a,DATASET_NAME,DATASET_DIM,size,data);
    status = danu_data_write_int(grp_b,DATASET_NAME,DATASET_DIM,size,data);
    status = danu_data_write_int(grp_c,DATASET_NAME,DATASET_DIM,size,data);

    /* Create a bunch of groups under GRPB */
    grp_num = 10;
    for(i=0;i<grp_num;i++) {
        sprintf(grp_name,"Group%04d",i);
        danu_group_create(grp_b,grp_name);
    }

    /* Create a bunch of datasets under GRPA */
    ds_num = 13;
    for(i=0;i<ds_num;i++) {
        sprintf(ds_name,"DataSet%04d",i);
        status = danu_data_write_int(grp_a,ds_name,DATASET_DIM,size,data);
    }
    

    /* Free memory */
    DANU_FREE(size);
    DANU_FREE(data);

    /* Close up file */
    danu_file_close(fid);


    return status;
}
    

/* Main */
int main(int argc, char ** argv)
{
    herr_t status;

    hid_t fid, gid;
    
    int     found;

    /* Delete the test file */
    if ( FILE_EXISTS(FILE) ) {
        danu_file_delete(FILE);
    }

    create_test_file();

    fid = danu_file_open(FILE,DANU_FILE_ACC_RDONLY,DANU_FILE_ACT_OPEN);
    gid = danu_group_open(fid,GRPB);

    status = danu_group_find_target_subgroup(gid,GRPC,&found);

    if ( status >= 0 && found == TRUE ) {
        printf("Found group %s in group %s\n",GRPC,GRPB);
    }
    else {
        DANU_ERROR_MESS("Failed to find group");
        goto TEST_FAIL;
    }

    status = danu_group_find_target_subgroup(gid,GRPA,&found);

    if ( status >= 0 || found == FALSE ) {
        printf("Did not find group %s in group %s\n",GRPA,GRPB);
    }
    else {
        DANU_ERROR_MESS("Failed to find group that did not exist");
        goto TEST_FAIL;
    }

   status = danu_group_find_target_dataset(gid,"Data",&found);

    if ( status >= 0 && found == TRUE ) {
        printf("Found dataset %s in group %s\n","Data",GRPB);
    }
    else {
        DANU_ERROR_MESS("Failed to find dataset");
        goto TEST_FAIL;
    }

    status = danu_group_find_target_dataset(gid,"Dummy That Does Not Exist",&found);

    if ( status >= 0 || found == FALSE ) {
        printf("Did not find dataset %s in group %s\n","Dummy That Does Not Exist",GRPB);
    }
    else {
        DANU_ERROR_MESS("Failed to find dataset that did not exist");
        goto TEST_FAIL;
    }

    
    //danu_file_delete(FILE);

    return 0;

TEST_FAIL:
    
    return 1;
}
