/*
*
*  Sample 2D Mesh 
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_utils.h>
#include <danu_error.h>
#include <danu_memory.h>
#include <danu_file.h>
#include <danu_mesh.h>

#define TEST_FILE "test-mesh.h5"

#define NCELLS  10
#define DIM     2
#define X_MIN   0.0
#define X_MAX   1.0
#define Y_MIN   0.0
#define Y_MAX   1.0

void print_buffer_double(int rank, size_t * size, double * buffer); 
void print_buffer_int(int rank, size_t * size, int * buffer); 


void print_buffer_double(int rank, size_t * size, double * buffer)
{
    didx_t i,j;
    double * p=buffer;
    for(i=0;i<rank;i++) {
        for(j=0;j<size[i];j++) {
            printf("BUF[%d,%d] = %1.5e\n",i,j,*p);
            p++;
        }
    }

}

void print_buffer_int(int rank, size_t * size, int * buffer)
{
    didx_t i,j;
    int * p=buffer;
    for(i=0;i<rank;i++) {
        for(j=0;j<size[i];j++) {
            printf("BUF[%d,%d] = %d\n",i,j,*p);
            p++;
        }
    }

}


int main(int argc, char ** argv)
{
    derr_t status;
    didx_t i;
    size_t size[DIM];
    size_t tmp;
    double *x, *y;

    hid_t fid;

    /* Initialize the x and y coordinate buffers */
    x = DANU_MALLOC(double,NCELLS+1);
    y = DANU_MALLOC(double,NCELLS+1);

    for(i=0;i<DIM;i++)
        size[i] = NCELLS+1;

    danu_rand_data_double(X_MIN,X_MAX,1,size,x);
    danu_rand_data_double(Y_MIN,Y_MAX,1,size,y);

    printf("X Coordinate Vector (x address 0x%0lX)\n",x);
    print_buffer_double(1,size,x);
    printf("Y Coordinate Vector (y address 0x%0lX)\n",y);
    print_buffer_double(1,size,y);

    /* Create a file */
    fid = danu_file_create(TEST_FILE);

    if ( H5_ISA_INVALID_ID(fid) ) {
        DANU_ERROR_MESS("ACK CAN NOT OPEN A TEST FILE");
        exit(DANU_FAIL);
    }

    return status;
}







