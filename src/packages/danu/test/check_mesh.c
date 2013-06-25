/* ****************************************************************************
*                                                                             *
*  Danu Unit Test                                                             *
*  Unit: File                                                                  *
*                                                                             *
* Requires Check Unit Test software package                                   *
* http://check.sourceforge.net/                                               *
*                                                                             *
* *************************************************************************** */ 
#include <stdlib.h>
#include <check.h>

#include <hdf5.h>

#include "unit_test_utils.h"
#include "random_generators.h"

#include <danu_error.h>
#include <danu_memory.h>

#include <danu_output.h>

#include <danu_mesh.h>

#define FILENAME   "deleteme-mesh.h5"
#define MESHNAME   "Test Mesh"


/* Test Utilities */

void create_random_meshes(const char *filename, int num_meshes);
void generate_random_coordinates(int nnodes, double *x, double *y, double *z);
void generate_random_connectivity(int size, int *connectivity);


void create_random_meshes(const char *filename,int num_meshes)
{
    hid_t fid,mid;
    herr_t status = output_file_create(filename,&fid);
    char mesh_name[64];
    telem_t elem;
    int seed = 0x123456;
    int m;

    if ( status != DANU_SUCCESS ) {
	printf("FAILED TO CREATE DUMMY TEST FILE\n");
	fail_exit_now();
    }

    for(m=1; m <= num_meshes; m++) {
        elem = generate_random_bound_int(LINE_ELEM,HEX_ELEM,&seed);
	sprintf(mesh_name,"Mesh %02d",m);
	status = mesh_create(fid,mesh_name,UNSTRUCTURED_MESH,elem,&mid);
	if ( status != DANU_SUCCESS ) {
	    printf("FAILED TO CREATE DUMMY MESH\n");
	    fail_exit_now();
	}
	danu_group_close(mid);
    }

    output_file_close(&fid);
}

void generate_random_coordinates(int nnodes, double *x, double *y, double *z)
{
  int iseed = 0xFF01;
  int n;

  for(n = 0; n < nnodes; n++)
    x[n] = generate_random_double(&iseed);
  
  if ( y != NULL ) {
    for(n = 0; n < nnodes; n++)
      y[n] = generate_random_double(&iseed);
  }

  if ( z != NULL ) {
    for(n = 0; n < nnodes; n++)
      z[n] = generate_random_double(&iseed);
  }

}

void generate_random_connectivity(int size, int *connectivity)
{
    int iseed = 0x1456;
    int n;

    for(n=0; n < size; n++)
	connectivity[n] = generate_random_int(&iseed);


}
/* BEGIN TESTS */

START_TEST (test_unstruct_mesh_exists)
{
    const char test_file[] = FILENAME;
    hid_t  fid,mid;
    herr_t status = output_file_create(test_file,&fid);
    char mesh_name[64];
    int exists;
 
    fail_unless ( status == DANU_SUCCESS,
                  "Failed to create output file");

    if ( status == DANU_SUCCESS ) {
	sprintf(mesh_name,"Mesh DNE");
	status = mesh_exists(fid,mesh_name,&exists);
	fail_unless(status == DANU_SUCCESS,
		    "Mesh Exists call failed");
	fail_unless(exists == FALSE,
		    "Return incorrect flag value (FALSE)");

	sprintf(mesh_name,"Mesh Test");
	status = mesh_create(fid,mesh_name,UNSTRUCTURED_MESH,HEX_ELEM,&mid);
	fail_unless(status == DANU_SUCCESS,
		    "Failed to create a test mesh");
	output_file_close(&fid);

	status = output_file_open_rdonly(test_file,&fid);
	fail_unless(status == DANU_SUCCESS,
		    "Failed to open test file");

	status = mesh_exists(fid,mesh_name,&exists);
	fail_unless(status == DANU_SUCCESS,
		    "Mesh Exists call failed");
	fail_unless(exists == TRUE,
		    "Failed to return correct exists flag (TRUE)");

    }

    output_file_close(&fid);
    danu_file_delete(test_file);

}
END_TEST


START_TEST (test_unstruct_mesh_open)
{ 
    const char test_file[] = FILENAME;
    hid_t  fid, mid;
    herr_t status;
    int seed = 0xABCD0123;
    int num_meshes = generate_random_bound_int(1,10,&seed);
    int i,m;
    char mesh_name[64];


    create_random_meshes(test_file,num_meshes);

    status = output_file_open_rdonly(test_file,&fid);
    fail_unless ( status == DANU_SUCCESS,
                  "Failed to open  output file");

    for(m=1;m<=num_meshes;m++) {
	sprintf(mesh_name,"Mesh %02d",m);
	mid = mesh_open(fid,mesh_name);
	fail_unless(H5_ISA_VALID_ID(mid),
		    "Failed to open existing mesh");
    }

    output_file_close(&fid);
    danu_file_delete(test_file);


}
END_TEST

START_TEST (test_unstruct_mesh_create)
{
    const char test_file[] = FILENAME;
    hid_t  fid,mid;
    herr_t status = output_file_create(test_file,&fid);
    char mesh_name[64];
    telem_t elem;



    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to create output file. Invalid H5 ID");

    if ( H5_ISA_VALID_ID(fid) ) {

	elem = -2;
	sprintf(mesh_name,"Failed Mesh");
	status = mesh_create(fid,mesh_name,UNSTRUCTURED_MESH,elem,&mid);
	fail_unless( status != DANU_SUCCESS,
		     "Failed to flag bad elem type");

	/* Create mesh for each element type */
	for(elem=LINE_ELEM; elem <= HEX_ELEM; elem++) {
	    sprintf(mesh_name,"%s elemtype=%d",MESHNAME, elem);
	    status = mesh_create(fid,mesh_name,UNSTRUCTURED_MESH,elem,&mid);
            fail_unless( status == DANU_SUCCESS,
	                 "Failed to create 2D unstructured mesh");
	}

	output_file_close(&fid);
	danu_file_delete(test_file);

    }

}
END_TEST

START_TEST(test_unstruct_mesh_add)
{ 
    const char test_file[] = FILENAME;
    hid_t  fid,mid;
    herr_t status;
    int d, m, e;
    char mesh_name[32] = "Test Unstructured";

    status = output_file_create(test_file,&fid);
    fail_unless(status == DANU_SUCCESS,
	        "Faile to create output file");

    /* These should fail */

    d = 1;
    for( e=3; e<= HEX_ELEM_ORDER ; e++) {
      status = mesh_add_unstructured(fid,mesh_name,e,d,&mid);
      fail_unless(status == DANU_FAILURE,
	          "Failed to return bad status with garbage input");
    }

    /* These should pass */
    d = 1;
    sprintf(mesh_name, "Test %dD",d);
    status = mesh_add_unstructured(fid,mesh_name,LINE_ELEM_ORDER,d,&mid);
    fail_unless( status == DANU_SUCCESS,
	         "Failed to add 1D mesh");

    d = 2;
    sprintf(mesh_name, "Quad Test %dD",d);
    status = mesh_add_unstructured(fid,mesh_name,QUAD_ELEM_ORDER,d,&mid);
    fail_unless( status == DANU_SUCCESS,
	         "Failed to add 2D quad mesh");

    sprintf(mesh_name, "TRI Test %dD",d);
    status = mesh_add_unstructured(fid,mesh_name,TRI_ELEM_ORDER,d,&mid);
    fail_unless( status == DANU_SUCCESS,
	         "Failed to add 2D tri mesh");

    d = 3;
    sprintf(mesh_name, "TET Test %dD",d);
    status = mesh_add_unstructured(fid,mesh_name,TET_ELEM_ORDER,d,&mid);
    fail_unless( status == DANU_SUCCESS,
	         "Failed to add 3D tet mesh");

    sprintf(mesh_name, "HEX Test %dD",d);
    status = mesh_add_unstructured(fid,mesh_name,HEX_ELEM_ORDER,d,&mid);
    fail_unless( status == DANU_SUCCESS,
	         "Failed to add 3d HEX mesh");



    output_file_close(&fid);
    danu_file_delete(test_file);

}
END_TEST

START_TEST(test_unstruct_mesh_count)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;
    int seed = 0xABCD0123;
    int num_meshes = generate_random_bound_int(1,10,&seed);
    int count;

    create_random_meshes(test_file,num_meshes);

    status = output_file_open_rdonly(test_file,&fid);
    fail_unless ( status == DANU_SUCCESS,
                  "Failed to open  output file");

    if ( status == DANU_SUCCESS ) {
	status = mesh_count(fid,&count);
	fail_unless(status == DANU_SUCCESS,
		"Failed to count the number of meshes");
	fail_unless(count == num_meshes,
		    "Number of meshes created does not match number count");
    }

    output_file_close(&fid);
    danu_file_delete(test_file);
}
END_TEST

START_TEST(test_unstruct_mesh_list)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;
    int seed = 0xABCD0123;
    int num_meshes = generate_random_bound_int(1,10,&seed);
    int i,m,num_found;
    char **mesh_names;
    char test_name[64];
    size_t len;

    create_random_meshes(test_file,num_meshes);

    status = output_file_open_rdonly(test_file,&fid);
    fail_unless ( status == DANU_SUCCESS,
                  "Failed to open  output file");

    if ( status == DANU_SUCCESS ) {

	mesh_names = DANU_MALLOC(char *, num_meshes);
	status = mesh_list(fid,num_meshes,mesh_names,&num_found);
	fail_unless(status == DANU_SUCCESS,
		    "Failed to read the list of meshes");

	for(m=1; m<=num_meshes; m++) {
	    i = m -1;
	    sprintf(test_name,"Mesh %02d",m);
	    len = strlen(test_name);
	    printf("test=%s read=%s\n",test_name,mesh_names[i]);
	    fail_unless( strncmp(test_name,mesh_names[i],len) == 0,
		         "Incorrect mesh name returned");
	    DANU_FREE(mesh_names[i]);
	}
	DANU_FREE(mesh_names);

    }

    output_file_close(&fid);
    danu_file_delete(test_file);
}
END_TEST

START_TEST(test_unstruct_mesh_coordinate_3d)
{
  const char test_file[] = FILENAME;
  const char mesh_name[] = "Test 3D Mesh";
  hid_t  fid, mid;
  herr_t status;
  int iseed = 0xABCD0123;
  int n,nnodes,i;
  double *x,*y,*z;
  double *rx,*ry,*rz;
  double **coord;

  status = output_file_create(test_file,&fid);
  fail_unless ( status == DANU_SUCCESS,
                "Failed to create output file");

  if ( status == DANU_SUCCESS ) {
    status = mesh_add_unstructured(fid,mesh_name,TET_ELEM_ORDER,3,&mid);
    fail_unless(status == DANU_SUCCESS,
		"Failed to create a test mesh");
    
    /* Generate the coordinates */
    nnodes = generate_random_bound_int(512,4096,&iseed);
    x = DANU_MALLOC(double,nnodes);
    rx = DANU_MALLOC(double,nnodes);
    y = DANU_MALLOC(double,nnodes);
    ry = DANU_MALLOC(double,nnodes);
    z = DANU_MALLOC(double,nnodes);
    rz = DANU_MALLOC(double,nnodes);
    generate_random_coordinates(nnodes,x,y,z);

    /* Calls to fail */
    status = mesh_write_coordinates(mid,nnodes,NULL,y,z);
    fail_unless(status != DANU_SUCCESS,
		"Failed to flag bad arguments (x)");
    
    status = mesh_write_coordinates(mid,nnodes,x,NULL,z);
    fail_unless(status != DANU_SUCCESS,
		"Failed to flag bad arguments (y)");

    status = mesh_write_coordinates(mid,nnodes,x,y,NULL);
    fail_unless(status != DANU_SUCCESS,
		"Failed to flag bad arguments (z)");


    /* Calls that should pass */
    status = mesh_write_coordinates(mid,nnodes,x,y,z);
    fail_unless(status == DANU_SUCCESS,
		"Failed to write 3D mesh coordinates");

    output_file_close(&fid);

    /* Open file for read only */
    status = output_file_open_rdonly(test_file,&fid);
    fail_unless( status == DANU_SUCCESS,
		 "Failed to open test file read only");

    mid = mesh_open(fid,mesh_name);
    fail_unless(H5_ISA_VALID_ID(mid),
		"Failed to open existing mesh group");

    /* Calls designed to fail */
    status = mesh_read_coordinates(mid,rx,NULL,NULL);
    fail_unless(status == DANU_FAILURE,
		"Failed to flag bad arguments (yz)");

    status = mesh_read_coordinates(mid,rx,ry,NULL);
    fail_unless(status == DANU_FAILURE,
		"Failed to flag bad arguments (z)");

    status = mesh_read_coordinates_byindex(mid,-1,rx);
    fail_unless(status == DANU_FAILURE,
		"Failed to flag bad arguments invalid index");

    status = mesh_read_coordinates_byindex(mid,3,rx);
    fail_unless(status == DANU_FAILURE,
		"Failed to flag bad arguments invalid index");
    
    status = mesh_read_coordinates_byindex(mid,0,NULL);
    fail_unless(status == DANU_FAILURE,
		"Failed to flag bad arguments invalid pointer");


    /* Calls designed to pass */
    status = mesh_read_coordinates(mid,rx,ry,rz);
    fail_unless(status == DANU_SUCCESS,
		"Failed to read 3D mesh coordinates");


    for(n=0; n < nnodes; n++) {
      fail_unless(x[n] == rx[n],
		  "Failed to read x-coordinate correctly");
      fail_unless(y[n] == ry[n],
		  "Failed to read y-coordinate correctly");
      fail_unless(z[n] == rz[n],
		  "Failed to read z-coordinate correctly");
    }

    /* Now zero out the data and read by coordinate */
    bzero(rx,sizeof(double)*nnodes);
    bzero(ry,sizeof(double)*nnodes);
    bzero(rz,sizeof(double)*nnodes);
    coord = DANU_MALLOC(double *,3);
    coord[0] = rx;
    coord[1] = ry;
    coord[2] = rz;
    for( i = 0; i < 3; i ++ ) {
      status = mesh_read_coordinates_byindex(mid,i,coord[i]);
      fail_unless(status == DANU_SUCCESS,
		  "Failed to read coordinates by index");
    }

    for(n=0; n < nnodes; n++) {
      fail_unless(x[n] == rx[n],
		  "Failed to read x-coordinate correctly");
      fail_unless(y[n] == ry[n],
		  "Failed to read y-coordinate correctly");
      fail_unless(z[n] == rz[n],
		  "Failed to read z-coordinate correctly");
    }


    DANU_FREE(coord);
    DANU_FREE(x);
    DANU_FREE(rx);
    DANU_FREE(y);
    DANU_FREE(ry);
    DANU_FREE(z);
    DANU_FREE(rz);

    output_file_close(&fid);

  }

  danu_file_delete(test_file);

}
END_TEST
START_TEST(test_unstruct_mesh_connectivity)
{
  const char test_file[] = FILENAME;
  const char mesh_name[] = "Test 3D Mesh";
  hid_t  fid, mid,cid;
  herr_t status;
  int iseed = 0xABCD0123;
  int nelem;
  int n, num,i,r_num;
  int *data, *rdata;

  status = output_file_create(test_file,&fid);
  fail_unless ( status == DANU_SUCCESS,
                "Failed to create output file");

  if ( status == DANU_SUCCESS ) {
    /* Generate the connectivity data */
    //num = generate_random_bound_int(2048,10000,&iseed);
    num=10;
    data = DANU_MALLOC(int,num*HEX_ELEM_ORDER);
    generate_random_connectivity(num*HEX_ELEM_ORDER,data);

    /* Calls that should fail */
#if 0
    status = mesh_write_connectivity(mid,num,data);
    fail_unless(status != DANU_SUCCESS,
	        "Failed flag bad mesh HDF5 id");
#endif    

    status = mesh_add_unstructured(fid,mesh_name,HEX_ELEM_ORDER,3,&mid);
    fail_unless(status == DANU_SUCCESS,
		"Failed to create a test mesh");

    status = mesh_write_connectivity(mid,-1,data);
    fail_unless(status != DANU_SUCCESS,
	        "Failed flag bad data size");

    status = mesh_write_connectivity(mid,num,NULL);
    fail_unless(status != DANU_SUCCESS,
	        "Failed flag bad data size");
    
    /* Write the connectivity data should pass */
    status=mesh_write_connectivity(mid,num,data);
    fail_unless(status == DANU_SUCCESS,
	        "Failed to write connectivity data");



    output_file_close(&fid);

    /* Now open and read data */
    status = output_file_open_rdonly(test_file,&fid);
    fail_unless(status == DANU_SUCCESS,
	        "Failed to open test file for read only");

    /* Calls that should fail */
    status = mesh_read_connectivity(mid,rdata);
    fail_unless(status != DANU_SUCCESS,
	        "Failed to flag bad mesh HDF5 id");

    mid=mesh_open(fid,mesh_name);
    status = mesh_read_connectivity(mid,NULL);
    fail_unless(status != DANU_SUCCESS,
	        "Failed to flag data pointer");

    rdata = DANU_CALLOC(int,num*HEX_ELEM_ORDER);
    status = mesh_read_connectivity(mid, rdata);
    fail_unless(status == DANU_SUCCESS,
	        "Failed to read connectivity data");

    for(n=0;n < num*HEX_ELEM_ORDER; n++) {
	fail_unless(rdata[n] == data[n],
		    "Failed to read connectivity data correctly");
    }

    /* Cleanup */
    output_file_close(&fid);

    DANU_FREE(data);
    DANU_FREE(rdata);

  }

}
END_TEST

START_TEST(test_unstruct_mesh_coordinate_2d)
{
  const char test_file[] = FILENAME;
  char mesh_name[32];
  hid_t  fid, mid, threed_mid;
  herr_t status;
  int iseed = 0xABCD0123;
  int n,nnodes,i;
  double *x,*y;

  status = output_file_create(test_file,&fid);
  fail_unless ( status == DANU_SUCCESS,
                "Failed to create output file");

  if ( status == DANU_SUCCESS ) {
    sprintf(mesh_name, "Test 3D Mesh");
    status = mesh_add_unstructured(fid,mesh_name,TET_ELEM_ORDER,3,&threed_mid);
    fail_unless(status == DANU_SUCCESS,
		"Failed to create a test 3D mesh");
    
    sprintf(mesh_name, "Test 2D Mesh");
    status = mesh_add_unstructured(fid,mesh_name,QUAD_ELEM_ORDER,2,&mid);
    fail_unless(status == DANU_SUCCESS,
		"Failed to create a test 2D mesh");
    
    /* Generate the coordinates */
    nnodes = generate_random_bound_int(512,4096,&iseed);
    x = DANU_MALLOC(double,nnodes);
    y = DANU_MALLOC(double,nnodes);
    generate_random_coordinates(nnodes,x,y,NULL);

    /* Calls to fail */
    status = mesh_write_coordinates_2d(threed_mid,nnodes,x,y);
    fail_unless(status != DANU_SUCCESS,
		"Failed to flag bad call to 3D Mesh");
    
    status = mesh_write_coordinates_2d(mid,nnodes,NULL,y);
    fail_unless(status != DANU_SUCCESS,
		"Failed to flag bad arguments (x)");

    status = mesh_write_coordinates_2d(mid,nnodes,x,NULL);
    fail_unless(status != DANU_SUCCESS,
		"Failed to flag bad arguments (y)");

    /* Calls that should pass */
    status = mesh_write_coordinates_2d(mid,nnodes,x,y);
    fail_unless(status == DANU_SUCCESS,
		"Failed to write 2D mesh coordinates");

    output_file_close(&fid);
    DANU_FREE(x);
    DANU_FREE(y);

  }

  danu_file_delete(test_file);

}
END_TEST
START_TEST(test_unstruct_mesh_coordinate_1d)
{
  const char test_file[] = FILENAME;
  char mesh_name[32];
  hid_t  fid, mid, threed_mid;
  herr_t status;
  int iseed = 0xABCD0123;
  int n,nnodes;
  double *x;

  status = output_file_create(test_file,&fid);
  fail_unless ( status == DANU_SUCCESS,
                "Failed to create output file");

  if ( status == DANU_SUCCESS ) {
    sprintf(mesh_name, "Test 3D Mesh");
    status = mesh_add_unstructured(fid,mesh_name,TET_ELEM_ORDER,3,&threed_mid);
    fail_unless(status == DANU_SUCCESS,
		"Failed to create a test 3D mesh");
    
    sprintf(mesh_name, "Test 2D Mesh");
    status = mesh_add_unstructured(fid,mesh_name,LINE_ELEM_ORDER,1,&mid);
    fail_unless(status == DANU_SUCCESS,
		"Failed to create a test 2D mesh");
    
    /* Generate the coordinates */
    nnodes = generate_random_bound_int(512,4096,&iseed);
    x = DANU_MALLOC(double,nnodes);
    generate_random_coordinates(nnodes,x,NULL,NULL);

    /* Calls to fail */
    status = mesh_write_coordinates_1d(threed_mid,nnodes,x);
    fail_unless(status != DANU_SUCCESS,
		"Failed to flag bad call to 3D Mesh");
    
    status = mesh_write_coordinates_1d(mid,nnodes,NULL);
    fail_unless(status != DANU_SUCCESS,
		"Failed to flag bad arguments (x)");

    /* Calls that should pass */
    status = mesh_write_coordinates_1d(mid,nnodes,x);
    fail_unless(status == DANU_SUCCESS,
		"Failed to write 1D mesh coordinates");

    output_file_close(&fid);
    DANU_FREE(x);

  }

  danu_file_delete(test_file);

}
END_TEST


Suite *
mesh_suite (void)
{
    Suite *s = suite_create("Danu Mesh");


    TCase *mesh_unstruct = tcase_create("Unstructured Mesh Create/Open");
    tcase_add_test(mesh_unstruct,test_unstruct_mesh_create);
    tcase_add_test(mesh_unstruct,test_unstruct_mesh_open);
    tcase_add_test(mesh_unstruct,test_unstruct_mesh_add);
    suite_add_tcase(s,mesh_unstruct);

    
    TCase *mesh_unstruct_query = tcase_create("Unstructured Mesh Query");
    tcase_add_test(mesh_unstruct_query,test_unstruct_mesh_exists);
    tcase_add_test(mesh_unstruct_query,test_unstruct_mesh_count);
    tcase_add_test(mesh_unstruct_query,test_unstruct_mesh_list);
    suite_add_tcase(s,mesh_unstruct_query);

    TCase *mesh_unstruct_coord = tcase_create("Unstructured Mesh Coordinates");
    tcase_add_test(mesh_unstruct_coord,test_unstruct_mesh_coordinate_3d);
    tcase_add_test(mesh_unstruct_coord,test_unstruct_mesh_coordinate_2d);
    tcase_add_test(mesh_unstruct_coord,test_unstruct_mesh_coordinate_1d);
    suite_add_tcase(s,mesh_unstruct_coord);

    TCase *mesh_unstruct_con = tcase_create("Unstructured Mesh Conectivity");
    tcase_add_test(mesh_unstruct_con,test_unstruct_mesh_connectivity);
    suite_add_tcase(s,mesh_unstruct_con);

    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = mesh_suite ();
  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
