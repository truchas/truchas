/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include <hdf5.h>

#include <danu.h>

#include<danu_sim.h>

#include <danu_series-group.h>
#include <danu_series-data.h>


#define FILE     "deleteme.h5"
#define SIM_NAME "Test Sim"

#define DATA_INT_NAME "NUMPE"
#define DATA_DBL_NAME "Dummy Data"

#define DUMMY_INT 0xabcd0123 
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
    char seriesname[128] = "Series 1";

    char data_name[128];
    char **names;

    int numpe;
    int i,dim, *size;
    size_t *nsizes;
    int ncnt,nfound,save;
    dsize_t num;
    double *data,*ref_data;
    double t0 = 0.0;
    double time;
    hid_t fid, sid, nsid;


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

    /* Add a series group to the simulation */
    if ( H5_RETURN_FAIL(sequence_getNextID(sid,t0,&nsid) ) ) {
        DANU_ERROR_MESS("Failed to get the next id handle");
        goto EXIT_NOW;
    }

    num_errs = 0;

    /* Arrays */
    size = DANU_MALLOC(int,DIM);

    /* Write the data */
    size[0] = 1;
    numpe = DUMMY_INT;
    sprintf(data_name,DATA_INT_NAME);
    if ( H5_RETURN_FAIL(simulation_data_write_int(nsid,data_name,1,size,&numpe)) ) {
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
    if ( H5_RETURN_FAIL(simulation_data_write_double(nsid,data_name,dim,size,data)) ) {
        DANU_ERROR_MESS("Failed to write double data");
        num_errs++;
        goto EXIT_NOW;
    }
  
    /* Check data by reading */ 
    ref_data = DANU_MALLOC(double,num);
    if ( H5_RETURN_FAIL(simulation_data_read_double(nsid,data_name,dim,size,ref_data) ) ) {
        num_errs++;
        goto EXIT_NOW;
    }
    for(i=0;i<num;i++) {
        if ( ref_data[i] != data[i] ) {
            num_errs++;
        }
    }

    /* Now loop through and create a buch of dummy series datasets */
    danu_group_close(nsid);
    i=0;
    time = t0;
    while ( i<10 ) {
        time+=0.05;
        if ( H5_RETURN_FAIL(sequence_getNextID(sid,time,&nsid) ) ) {
            DANU_ERROR_MESS("Failed to create a series group");
            num_errs++;
        }
        if (H5_RETURN_FAIL(simulation_data_write_double(nsid,data_name,dim,size,data) ) ) {
            DANU_ERROR_MESS("Failed to write double data set");
            num_errs++;
        }
        danu_group_close(nsid);
        i++;
    }

    /* Find the number of datasets */
    sequence_count(sid,&ncnt);
    printf("Found %d series groups\n",ncnt);
    names = DANU_MALLOC(char *,ncnt);
    nsizes = DANU_MALLOC(size_t, ncnt); 
    save = ncnt;
    for(i=0;i<ncnt;i++) {
        nsizes[i] = 128;
        names[i] = DANU_MALLOC(char,128);
    }
    sequence_list(sid,ncnt,nsizes,names,&nfound);
    printf("Found the following series groups\n");
    for(i=0;i<ncnt;i++) {
        printf("\t<%s>\n",names[i]);
    }

    /* Open the first series group */
    sequence_get_handle(sid,seriesname,&nsid);
    simulation_data_count(nsid,&ncnt);
    printf("Found %d datasets in %s\n", ncnt, seriesname);
    simulation_data_list(nsid,ncnt,nsizes,names,&nfound);
    for(i=0; i<nfound; i++ ) {
        printf("\t<%s>\n", names[i]);
    }



    
    /* Free memory */
    for(i=0;i<save;i++) {
        DANU_FREE(names[i]);
    }
    DANU_FREE(names);
    DANU_FREE(ref_data);
    DANU_FREE(data);
    DANU_FREE(nsizes);
    DANU_FREE(size);

    /* Free HDF5 resources */
    danu_group_close(sid);
    danu_file_close(fid);


EXIT_NOW:
    printf("Found %d errors\n",num_errs);
    return num_errs;

}
