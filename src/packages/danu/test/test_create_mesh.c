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

    double *coordinates[DIM];

    int exists;
    int nmesh;
    char **meshnames;
    int i,*size;



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

    /* Query the existence of the mesh */
    if ( DANU_RETURN_FAIL(mesh_exists(fid,meshname,&exists) ) ) {
        DANU_ERROR_MESS("Failed to query mesh");
        num_errs++;
        goto EXIT_NOW;
    }

    if ( exists ) {
        printf("Mesh existence query successful\n");
    }

    /* Query a mesh that does not exist */
    if ( DANU_RETURN_FAIL(mesh_exists(fid,"DUHMesh",&exists) ) ) {
        DANU_ERROR_MESS("Failed to query fake mesh");
        num_errs++;
        goto EXIT_NOW;
    }

    if ( ! exists ) {
        printf("Query of fake mesh passed\n");
    }

    /* Create another mesh */
    sprintf(meshname,MESH_NAME2);
    mesh_type = STRUCTURED_MESH;
    elem_type = HEX_ELEM;
    if ( DANU_RETURN_FAIL(mesh_create(fid,meshname,DIM,mesh_type,elem_type)) ) {
        DANU_ERROR_MESS("Failed to create the second mesh group");
        num_errs++;
        goto EXIT_NOW;
    }

    /* Count the number of meshes */
    if ( DANU_RETURN_FAIL(mesh_count(fid,&nmesh) ) ) {
        DANU_ERROR_MESS("Failed to count meshes");
        num_errs++;
        goto EXIT_NOW;
    }


    printf("Found %d meshes in %s\n", nmesh, root_file);

    /* Now get the names */
    meshnames = DANU_MALLOC(char *, nmesh);
    size      = DANU_MALLOC(int,nmesh);
    for(i=0;i<nmesh;i++) {
        size[i] = 64;
        meshnames[i] = DANU_MALLOC(char, size[i]);
    }

    if ( DANU_RETURN_FAIL(mesh_list(fid,nmesh,size,meshnames) ) ) {
        DANU_ERROR_MESS("Failed to get mesh list");
        num_errs++;
        goto EXIT_NOW;
    }

    DANU_FREE(size);
    for(i=0;i<nmesh;i++) {
        printf("\t-->%s\n",meshnames[i]);
        DANU_FREE(meshnames[i]);
    }
    DANU_FREE(meshnames);
       

    danu_file_close(fid);

    goto EXIT_NOW;


EXIT_NOW:
    printf("Found %d errors\n",num_errs);
    return num_errs;

}
