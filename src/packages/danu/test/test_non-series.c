#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include <hdf5.h>

#include <danu.h>

#include<danu_sim.h>

#include <danu_non-series.h>


#define FILE     "deleteme.h5"
#define SIM_NAME "Test Sim"

#define DATA_INT_NAME "NUMPE"
#define DATA_DBL_NAME "Dummy Data"

#define DUMMY_INT 8
#define DIM 2
#define NUMCELLS  16
herr_t create_test_file(const char *name);


herr_t create_test_file(const char *name)
{
    herr_t status = DANU_FAILURE;

    hid_t fid;

    fid = danu_file_create(name);
    if ( H5_ISA_VALID_ID(fid) ) {
        status = simulations_create_root_group(fid);
        danu_file_close(fid);
    }

    return status;
}

int main(int argc, char ** argv)
{
    int num_errs;

    char root_file[] = FILE;
    char test_sim[]  = SIM_NAME;

    char data_name[128];
    char **names;

    int numpe;
    int i,dim, *size;
    size_t *nsizes;
    int ncnt;
    dsize_t num;
    double *data,*ref_data;
    hid_t fid, sid;


    if ( FILE_EXISTS(root_file) ) {
        danu_file_delete(root_file);
    }

    if ( H5_RETURN_FAIL(create_test_file(root_file) ) ) {
        DANU_ERROR_MESS("Failed to create the test file");
        num_errs++;
        goto EXIT_NOW;
    }

    fid = danu_file_open_append(root_file);
    if ( H5_RETURN_FAIL(simulation_add(fid,test_sim,&sid) ) ) {
        DANU_ERROR_MESS("Failed to add simulation");
        goto EXIT_NOW;
    }

    num_errs = 0;

    /* Arrays */
    size = DANU_MALLOC(int,DIM);

    /* Write the data */
    size[0] = 1;
    numpe = DUMMY_INT;
    sprintf(data_name,DATA_INT_NAME);
    if ( H5_RETURN_FAIL(data_write_int(sid,data_name,1,size,&numpe)) ) {
        DANU_ERROR_MESS("Failed to write int data");
        num_errs++;
        goto EXIT_NOW;
    }

    dim = DIM;
    num = 1;
    for(i=0;i<dim;i++) {
       size[i] = NUMCELLS;
       num*=NUMCELLS;
    }
    data = DANU_MALLOC(double,num);
    sprintf(data_name,DATA_DBL_NAME);
    danu_rand_data_double(-5.0,5.0,num,data);
    if ( H5_RETURN_FAIL(data_write_double(sid,data_name,dim,size,data)) ) {
        DANU_ERROR_MESS("Failed to write double data");
        num_errs++;
        goto EXIT_NOW;
    }
  
    /* Check data by reading */ 
    ref_data = DANU_MALLOC(double,num);
    if ( H5_RETURN_FAIL(data_read_double(sid,data_name,dim,size,ref_data) ) ) {
        num_errs++;
        goto EXIT_NOW;
    }
    for(i=0;i<num;i++) {
        if ( ref_data[i] != data[i] ) {
            num_errs++;
        }
    }

    /* Find the number of datasets */
    data_count(sid,&ncnt);
    printf("Found %d datasets\n",ncnt);
    names = DANU_MALLOC(char *,ncnt);
    nsizes = convert_int_to_size(ncnt,size);
    for(i=0;i<ncnt;i++) {
        size[i] = 128;
        names[i] = DANU_MALLOC(char,128);
    }
    data_list(sid,ncnt,nsizes,names);
    printf("Found the following datasets\n");
    for(i=0;i<ncnt;i++) {
        printf("\t<%s>\n",names[i]);
    }



    
    /* Free memory */
    for(i=0;i<ncnt;i++) {
        DANU_FREE(names[i]);
    }
    DANU_FREE(names);
    DANU_FREE(ref_data);
    DANU_FREE(data);    
    DANU_FREE(size);

    /* Free HDF5 resources */
    danu_group_close(sid);
    danu_file_close(fid);


EXIT_NOW:
    printf("Found %d errors\n",num_errs);
    return num_errs;

}
