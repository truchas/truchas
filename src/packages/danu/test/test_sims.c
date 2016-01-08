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

#include <danu_sim.h>


#define FILE "deleteme.h5"
#define SIM1_NAME "Simulation 1"
#define SIM2_NAME "EM Simulation"
#define MAX_NAME_LEN 128


herr_t create_test_file(void);


herr_t create_test_file(void)
{
    herr_t status;


    return status;
}
    

/* Main */
int main(int argc, char ** argv)
{
    herr_t num_fails;

    hid_t fid,sid,sid1,sid2;
    int   cnt, i;
    char * names[2];
    size_t * size;
    
    /* Delete the test file */
    if ( FILE_EXISTS(FILE) ) {
        danu_file_delete(FILE);
    }

    fid = danu_file_create(FILE);
    DANU_DEBUG_MESS("Create test file");

    if ( H5_ISA_VALID_ID(fid) ) {

        num_fails = simulations_create(fid,&sid,FALSE);
        DANU_DEBUG_MESS("Create Simulations groups");
        num_fails += simulation_add(fid,SIM1_NAME,&sid1,FALSE);
        DANU_DEBUG_MESS("Add a simulation");
        num_fails += simulation_add(fid,SIM2_NAME,&sid2,FALSE);
        DANU_DEBUG_MESS("Add another simulation");

        danu_group_close(sid2);
        danu_group_close(sid1);
        danu_group_close(sid);

        num_fails += simulation_count(fid,&cnt);
        printf("There are %d simulations in %s\n",cnt,FILE);
        size = DANU_MALLOC(size_t,2);
        for(i=0;i<2;i++) {
            size[i] = MAX_NAME_LEN;
            names[i] = DANU_MALLOC(char,MAX_NAME_LEN);
        }
        num_fails += simulation_list(fid,cnt,size,names);
        printf("The simulations are:\n");
        for(i=0;i<2;i++) 
            printf("\t%s\n",names[i]);

        DANU_FREE(size);
        for(i=0;i<2;i++) {
            DANU_FREE(names[i]);
        }
    }


    return num_fails;
}
