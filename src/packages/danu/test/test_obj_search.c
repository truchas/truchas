/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
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

#define DATASET_NAME  "Data"
#define DATASET_DIM   2
#define DATASET_SIZE  2
#define DATASET_VALUE 0xABCD0123


herr_t create_test_file(void);


herr_t create_test_file(void)
{
    herr_t status;

    hid_t fid, grp_a, grp_b, grp_c;

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
    ds_num = 50;
    for(i=0;i<50;i++) {
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

    h5o_node_list_t *list;
    h5o_node_t      *ptr;


    /* Delete the test file */
    if ( FILE_EXISTS(FILE) ) {
        danu_file_delete(FILE);
    }

    create_test_file();

    fid = danu_file_open(FILE,DANU_FILE_ACC_RDONLY,DANU_FILE_ACT_OPEN);
    gid = danu_group_open(fid,GRPA);

    list = h5_object_list_create(H5O_TYPE_DATASET);

    status = h5_object_search(gid,H5O_TYPE_DATASET,list);


    /* print out */
    ptr = list->root;
    printf("Search for ");
    if ( list->type == H5O_TYPE_GROUP ) {
        printf("groups\n");
    }
    else if ( list->type == H5O_TYPE_DATASET ) {
        printf("datasets\n");
    }
    else {
        printf("Invalid serach type\n");
    }
    printf("Number found %d\n", h5_object_list_nnodes(list));
    while(ptr != NULL) {
        printf("Object id: %d\n", ptr->id);
        printf("Object Name: %s\n", ptr->name);
        printf("Address: 0x%lX\n",ptr->addr);
        ptr = ptr->next;
    }

    /* Free memory */
    h5_object_list_delete(list);

    //danu_file_delete(FILE);

    return 0;
}
