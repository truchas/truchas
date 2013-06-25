#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include <hdf5.h>

#include <danu.h>

#include<danu_sim.h>

#include <danu_probes.h>


#define FILE     "deleteme.h5"
#define SIM_NAME "Test Sim"

#define DATA_INT_NAME "NUMPE"
#define DATA_DBL_NAME "Dummy Data"

#define DUMMY_INT 8
#define DIM 3
#define NUMCELLS 128 
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
    char probeName[128] = "Viscosity";

    char data_name[128];
    char **names;

    int numpe;
    int i,dim, *size;
    size_t *nsizes;
    int ncnt, exists;
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

    /* Check existence of probe data */
    if ( H5_RETURN_FAIL(probe_exists(sid,probeName,&exists) ) ) {
        DANU_ERROR_MESS("Failed to stat probe dataset");
        num_errs++;
        goto EXIT_NOW;
    }
    
    if ( exists ) {
        printf("Probe name = %s does exist\n", probeName);
    }
    else {
        printf("Probe name = %s does not exist\n", probeName);
    }
        

    /* Arrays */
    size = DANU_MALLOC(int,DIM);

    /* Write the data */
    size[0] = 1;
    numpe = DUMMY_INT;
    sprintf(data_name,DATA_INT_NAME);
    if ( H5_RETURN_FAIL(probe_data_write_int(sid,data_name,1,size,&numpe)) ) {
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
    danu_rand_data_double(-10.0,10.0,num,data);
    if ( H5_RETURN_FAIL(probe_data_write_double(sid,data_name,dim,size,data)) ) {
        DANU_ERROR_MESS("Failed to write double data");
        num_errs++;
        goto EXIT_NOW;
    }
  
    /* Check data by reading */ 
    ref_data = DANU_MALLOC(double,num);
    if ( H5_RETURN_FAIL(probe_data_read_double(sid,data_name,dim,size,ref_data) ) ) {
        DANU_ERROR_MESS("Failed to read the double data");
        num_errs++;
        goto EXIT_NOW;
    }
    for(i=0;i<num;i++) {
        if ( ref_data[i] != data[i] ) {
            num_errs++;
        }
    }

    /* Append probe data */
    i = 0;
    size[0] = 1;
    while ( i< 10 ) {
        if ( H5_RETURN_FAIL(probe_data_write_int(sid,probeName,1,size,&numpe) ) ) {
            DANU_ERROR_MESS("Failed to write int data");
            num_errs++;
        }
        i++;
    }
        

    /* Find the number of datasets */
    probe_count(sid,&ncnt);
    printf("Found %d probe datasets\n",ncnt);
    names = DANU_MALLOC(char *,ncnt);
    nsizes = DANU_MALLOC(size_t,ncnt);
    for(i=0;i<ncnt;i++) {
        nsizes[i] = 128;
        names[i] = DANU_MALLOC(char,128);
    }
    probe_list(sid,ncnt,nsizes,names);
    printf("Found the following probe datasets\n");
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
    DANU_FREE(nsizes);

    /* Free HDF5 resources */
    danu_group_close(sid);
    danu_file_close(fid);


EXIT_NOW:
    printf("Found %d errors\n",num_errs);
    return num_errs;

}
