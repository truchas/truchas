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
#include <danu_mesh.h>
#include <danu_output.h>

#include <danu_probes.h>


#define FILE       "deleteme.h5"
#define MESH_NAME1 "TestMesh"
#define MESH_NAME2 "EM Mesh"

#define DIM 3 
#define NNODES 128

int main(int argc, char ** argv)
{
    int num_errs = 0;

    hid_t fid, mid;

    char root_file[] = FILE;
    char meshname[] = MESH_NAME1;
    tmesh_t mesh_type = UNSTRUCTURED_MESH;
    telem_t elem_type = HEX_ELEM;
    int elem_order    = HEX_ELEM_ORDER;

    double *coordinates[DIM], *rd_coordinates[DIM];

    int status;
    int i;



    if ( DANU_RETURN_FAIL(output_file_create(root_file,&fid)) ) {
        DANU_ERROR_MESS("Failed to create the output file");
        num_errs++;
        goto EXIT_NOW;
    }


    /* Create a mesh */
    if ( DANU_RETURN_FAIL(mesh_create(fid,meshname,DIM,mesh_type,elem_type)) ) {
        DANU_ERROR_MESS("Failed to create the mesh group");
        num_errs++;
        goto EXIT_NOW;
    }

    mid = mesh_open(fid,meshname);

    /* Generate coordinate data and write to a file */
    coordinates[0] = NULL;
    coordinates[1] = NULL;
    coordinates[2] = NULL;
    rd_coordinates[0] = NULL;
    rd_coordinates[1] = NULL;
    rd_coordinates[2] = NULL;
    for(i=0; i< DIM; i++ ) {
        coordinates[i] = DANU_MALLOC(double, NNODES);
        rd_coordinates[i] = DANU_MALLOC(double, NNODES);
        danu_rand_data_double(-1.0,1.0,NNODES,coordinates[i]);
    }

    status = mesh_write_coordinates(mid,NNODES,coordinates[0],coordinates[1],coordinates[2]);
    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to write coordinates");
        num_errs++;
        goto EXIT_NOW;
    }
    output_file_close(fid);

    /* Reopen file and  read coordinates */
    if ( DANU_RETURN_FAIL(output_file_open_rdonly(root_file,&fid) ) ) {
        DANU_ERROR_MESS("Failed to open the file for reading");
        num_errs++;
        goto EXIT_NOW;
    }

    mid = mesh_open(fid,meshname);

    status = mesh_read_coordinates(mid,rd_coordinates[0], rd_coordinates[1], rd_coordinates[2]);
    if ( DANU_RETURN_FAIL(status) ) {
        DANU_ERROR_MESS("Failed to read coordinates");
        num_errs++;
        goto EXIT_NOW;
    }

    for(i=0; i<DIM; i++ ) {
        if ( memcmp(coordinates[i],rd_coordinates[i],sizeof(double)*NNODES) != 0 ) {
            DANU_ERROR_MESS("Mismatch between read/write data");
            num_errs++;
        }
    }

    /* Close the file */
    output_file_close(fid);


    /* Free memory */
    for(i=0; i< DIM; i++ ) {
        DANU_FREE(coordinates[i]);
        DANU_FREE(rd_coordinates[i]);
    }


    goto EXIT_NOW;


EXIT_NOW:
    printf("Found %d errors\n",num_errs);
    return num_errs;

}
