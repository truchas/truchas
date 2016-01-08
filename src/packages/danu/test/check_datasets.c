/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
/* ****************************************************************************
*                                                                             *
*  Danu Unit Test                                                             *
*  Unit: Datasets                                                             *
*                                                                             *
* Requires Check Unit Test software package                                   *
* http://check.sourceforge.net/                                               *
*                                                                             *
* *************************************************************************** */ 
#include <stdlib.h>
#include <string.h>

#include <check.h>

#include <hdf5.h>

#include "unit_test_utils.h"
#include "random_generators.h"

#include <danu_memory.h>
#include <danu_error.h>
#include <danu_dataset.h>

#define INT_DATA_NAME    "Integer Data"
#define DOUBLE_DATA_NAME "Double Data"
#define FLOAT_DATA_NAME  "Float Data"

void create_dummy_data_int(hid_t loc,  int ndim, hsize_t *size);

void create_dummy_data_int(hid_t loc, int ndim, hsize_t *size)
{
  int *data;
  int i,num;
  herr_t status;
  const char data_name[] = INT_DATA_NAME;
  static int iseed = 0;

  num = 1;
  for(i=0;i<ndim;i++) {
    num*=(int)size[i];
  }

  data = DANU_MALLOC(int,num);
  generate_random_int_data(data,num,&iseed);
  status = danu_data_write_int(loc,data_name,ndim,size,data);

  if ( H5_RETURN_FAIL(status) ) {
    DANU_ERROR_MESS("Failed to create dummy data set");
    exit_now(1);
  }

  DANU_FREE(data);

}



START_TEST (test_dataset_dim)
{
    hid_t  fid, did;
    herr_t status;
    int iseed = 0;
    hsize_t *dims;
    int idx,rank,ndim,max_size,num;
    char data_name[] = INT_DATA_NAME;

    create_test_h5_file();
    fid = open_test_file();

    ndim = generate_random_bound_int(1,3,&iseed);
    max_size = generate_random_bound_int(1,128,&iseed);

    dims = DANU_MALLOC(hsize_t,ndim);
    num = 1;
    for(idx=0; idx < ndim; idx++) {
      dims[idx] = (hsize_t) generate_random_bound_int(1,max_size,&iseed);
      num*=(int)dims[idx];
    }

    create_dummy_data_int(fid,ndim,dims);

    /* Find the rank */
    rank = danu_dataset_rank(fid,data_name);
    fail_unless(rank > 0,
		"Failed to return positive rank value");
    fail_unless(rank == ndim,
		"Failed to return correct rank value");

     
}
END_TEST

START_TEST (test_dataset_dimension)
{
    hid_t  fid, did;
    herr_t status;
    int iseed = 0;
    hsize_t *dims,*size;
    int idx,rank,ndim,max_size,num;
    char data_name[] = INT_DATA_NAME;

    create_test_h5_file();
    fid = open_test_file();

    ndim = generate_random_bound_int(1,3,&iseed);
    max_size = generate_random_bound_int(1,128,&iseed);

    dims = DANU_MALLOC(hsize_t,ndim);
    size = DANU_MALLOC(hsize_t,ndim);
    num = 1;
    for(idx=0; idx < ndim; idx++) {
      dims[idx] = (hsize_t) generate_random_bound_int(1,max_size,&iseed);
      num*=(int)dims[idx];
    }

    create_dummy_data_int(fid,ndim,dims);

    /* Find the dimensions */
    status = danu_dataset_dimensions(fid,data_name,ndim,size);
    fail_unless(H5_RETURN_OK(status),
		"Failed to return pass status for dimension check");

    for(idx=0;idx < ndim; idx++) {
      fail_unless(dims[idx] == size[idx],
		  "Failed to return correct dimension size");
    }

     
}
END_TEST

START_TEST (test_dataset_max_dimension)
{
    hid_t  fid, did;
    herr_t status;
    int iseed = 0;
    hsize_t *dims,*size;
    int idx,rank,ndim,max_size,num;
    char data_name[] = INT_DATA_NAME;

    create_test_h5_file();
    fid = open_test_file();

    ndim = generate_random_bound_int(1,3,&iseed);
    max_size = generate_random_bound_int(1,128,&iseed);

    dims = DANU_MALLOC(hsize_t,ndim);
    size = DANU_MALLOC(hsize_t,ndim);
    num = 1;
    for(idx=0; idx < ndim; idx++) {
      dims[idx] = (hsize_t) generate_random_bound_int(1,max_size,&iseed);
      num*=(int)dims[idx];
    }

    create_dummy_data_int(fid,ndim,dims);

    /* Find the dimensions */
    status = danu_dataset_max_dimensions(fid,data_name,ndim,size);
    fail_unless(H5_RETURN_OK(status),
		"Failed to return pass status for max dimension check");

    for(idx=0;idx < ndim; idx++) {
      fail_unless(dims[idx] == size[idx],
		  "Failed to return correct dimension size");
    }

     
}
END_TEST

/* Define the test suite */
Suite *
dataset_suite (void)
{
    Suite *s = suite_create("Danu Datasets");

    TCase *dataset_utils = tcase_create("Dataset Utils");
    tcase_add_test(dataset_utils,test_dataset_dim);
    tcase_add_test(dataset_utils,test_dataset_dimension);
    tcase_add_test(dataset_utils,test_dataset_max_dimension);
    suite_add_tcase(s,dataset_utils);

    return s;
}


      
int main(void)
{
  int number_failed;
  Suite *s = dataset_suite ();
  SRunner *sr = srunner_create (s);
  
  srunner_set_log(sr, "dataset-results.log");
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}






