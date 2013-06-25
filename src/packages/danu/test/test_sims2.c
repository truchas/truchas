#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include <hdf5.h>

#include <danu.h>

#include <danu_sim.h>
#include <danu_non-series.h>
#include <danu_series-group.h>
#include <danu_series-data.h>


#define FILE "deleteme.h5"
#define NUM_SIM   3
#define SIM1_NAME "Simulation 1"
#define SIM2_NAME "Simulation 2"
#define SIM3_NAME "EM Simulation"

#define INT_DATA_NAME1 "Int Data Set 1"
#define INT_DATA_NAME2 "Int Data Set 2"

#define DATA_DIM 2
#define DATA_SIZE 8


herr_t create_test_file(void)
{
    herr_t status, num_fails;
    hid_t fid,sid;

    fid = danu_file_create(FILE);

    if ( H5_ISA_VALID_ID(fid) ) {
        num_fails = simulations_create(fid,&sid,FALSE);
        danu_group_close(sid);
    }

    status = simulation_add(fid,SIM1_NAME,&sid,FALSE);
    danu_group_close(sid);

    status = simulation_add(fid,SIM2_NAME,&sid,FALSE);
    danu_group_close(sid);

    danu_file_close(fid);

    return status;
}

int * generate_size_array(int dim)
{
    int *ptr = DANU_MALLOC(int,dim);
    int i;
    
    for(i=0; i< dim; i++)
        ptr[i] = DATA_SIZE;

    return ptr;
}

int * generate_int_data(int dim, int *size, int *num)
{
    int * data;
    int  i;
    int tmp_num;

    tmp_num=1;
    for(i=0;i<dim;i++) {
        (tmp_num)*=size[i];
    }

    data = DANU_MALLOC(int,tmp_num);

    danu_rand_data_int(0,10000,(dsize_t)tmp_num,data);

    *num=tmp_num;


    return data;
}
double * generate_double_data(int dim, const int *size, int *num)
{
    double * data;
    int  i;
    int tmp_num;

    tmp_num=1;
    for(i=0;i<dim;i++) {
        (tmp_num)*=size[i];
    }

    data = DANU_MALLOC(double,tmp_num);

    danu_rand_data_double(-10.0,10.0,(dsize_t)tmp_num,data);

    *num=tmp_num;


    return data;
}

/* Main */
int main(int argc, char ** argv)
{
    herr_t status, num_fails;

    hid_t fid,sid,nsid;

    int     i, data_num,*data_size, *idata, *rd_buf;

    /* Delete the test file */
    if ( FILE_EXISTS(FILE) ) {
        danu_file_delete(FILE);
    }

    create_test_file();

    fid = danu_file_open(FILE,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_APPEND);

    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("Failed to open test file");
        goto FAIL_EXIT;
    }

    if ( simulation_open(fid,SIM1_NAME,&sid) < 0 ) {
        DANU_ERROR_MESS("Failed to open simulation");
        goto FAIL_EXIT;
    }
    
    if ( sequence_getNextID(sid,1.0,&nsid) < 0 ) {
        DANU_ERROR_MESS("Failed to get the next ID");
        goto FAIL_EXIT;
    }

    /* Generate garbage data */
    data_size = generate_size_array(DATA_DIM);
    idata = generate_int_data(DATA_DIM,data_size,&data_num);

    if ( simulation_data_write_int(nsid,INT_DATA_NAME1,DATA_DIM,data_size,idata) < 0 ) {
        DANU_ERROR_MESS("Failed to write data set 1");
        goto FAIL_EXIT;
    }
    
    if ( simulation_data_write_int(nsid,INT_DATA_NAME2,DATA_DIM,data_size,idata) < 0 ) {
        DANU_ERROR_MESS("Failed to write data set 2");
        goto FAIL_EXIT;
    }
    

    /* Close the file */
    danu_group_close(nsid);
    danu_group_close(sid);
    danu_file_close(fid);

    /* Now reopen and read the data */
    fid = danu_file_open(FILE,DANU_FILE_ACT_OPEN,DANU_FILE_ACC_RDONLY);

    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("Failed to open file for reading");
        goto FAIL_EXIT;
    }

    if ( simulation_open(fid,SIM1_NAME,&sid) < 0 ) {
        DANU_ERROR_MESS("Failed to open simulation group");
        goto FAIL_EXIT;
    }

    if ( sequence_get_handle(sid,"Series 1",&nsid) < 0 ) {
        DANU_ERROR_MESS("Failed to open the series group");
        goto FAIL_EXIT;
    }

    /* Malloc */
    rd_buf = DANU_MALLOC(int,data_num);
    if ( simulation_data_read_int(nsid,INT_DATA_NAME1,DATA_DIM,data_size,rd_buf) < 0 ) {
        DANU_ERROR_MESS("Failed to read data set");
        goto FAIL_EXIT;
    }

    /* Now check the data */
    for(i=0;i<data_num;i++) {
        if ( rd_buf[i] != idata[i] ) {
            DANU_ERROR_MESS("CHECK FAILED");
            printf("INDEX=%d READ=%d OUT=%d\n", i,rd_buf[i], idata[i]);
            goto FAIL_EXIT;
        }
    }

    /* WE PASSSED! */
    printf("YES wrote and read in data set no errors!\n");

    /* Memory cleanup */
    DANU_FREE(rd_buf);
    DANU_FREE(idata);
    DANU_FREE(data_size);


    /* Close the file....again */
    danu_group_close(nsid);
    danu_group_close(sid);
    danu_file_close(fid);

    


PASS_EXIT:
    return DANU_SUCCESS;
FAIL_EXIT:
    return DANU_FAIL;
}
