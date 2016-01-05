/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
/* ****************************************************************************
*                                                                             *
*  Danu Unit Test                                                             *
*  Unit:  Non-series datasets                                                 *
*                                                                             *
* Requires Check Unit Test software package                                   *
* http://check.sourceforge.net/                                               *
*                                                                             *
* *************************************************************************** */ 
#include <stdlib.h>

#include <hdf5.h>
#include <check.h>

#include "unit_test_utils.h"

#include <danu_types.h>
#include <danu_h5_error.h>
#include <danu_memory.h>
#include <danu_error.h>
#include <danu_utils.h>
#include <danu_file.h>

#include <danu_output.h>
#include <danu_sim.h>

#include <danu_non-series.h>

#define FILENAME "deleteme_non-series.h5"
#define SIMNAME  "Test Simulation"

void generate_output_file(hid_t * fid, hid_t *sid)
{
  char test_file[] = FILENAME;
  char sim_name[]  = SIMNAME;

  herr_t err = output_file_create(test_file,fid);

  if ( H5_RETURN_FAIL(err) ) {
    danu_error_printf("Failed to create test output file");
    exit_now(err);
  }

  err = simulation_add(*fid,sim_name,sid);

  if ( H5_RETURN_FAIL(err) ) {
    danu_error_printf("Failed to create test sim group");
    exit_now(err);
  }

}

void open_file(hid_t *fid, hid_t *sid)
{
  char test_file[] = FILENAME;
  char sim_name[]  = SIMNAME;

  herr_t err = output_file_open_rdonly(test_file,fid);

  if ( H5_RETURN_FAIL(err) ) {
    danu_error_printf("Failed to create test output file");
    exit_now(err);
  }

  err = simulation_open(*fid,sim_name,sid);

  if ( H5_RETURN_FAIL(err) ) {
    danu_error_printf("Failed to create test sim group");
    exit_now(err);
  }

}

void file_delete()
{
  danu_file_delete(FILENAME);
}


/* BEGIN TESTS */

START_TEST (test_nsdata_write_int)
{
  hid_t fid, sid;
  char data_name[] = "Test Non-Series INT";
  int iseed = 2090;
  herr_t err;
  int dim = generate_random_bound_int(1,3,&iseed);
  int *dimensions;
  int *data;
  int data_size,idx;

  generate_output_file(&fid,&sid);

  dimensions = DANU_MALLOC(int,dim);
  data_size = 1;
  for(idx=0; idx < dim; idx++) {
    dimensions[idx] = generate_random_bound_int(10,1024,&iseed);
    data_size*=dimensions[idx];
  }
  data = DANU_MALLOC(int,data_size);

  generate_random_int_data(data,data_size,&iseed);


  err = data_write_int(sid,data_name,dim, dimensions, data);
  fail_unless(H5_RETURN_OK(err), 
	      "Failed to write non-series int data");


  DANU_FREE(dimensions);
  DANU_FREE(data);

  file_delete();

}
END_TEST
START_TEST (test_nsdata_write_double)
{
  hid_t fid, sid;
  char data_name[] = "Test Non-Series DOUBLE";
  int iseed = 2090;
  herr_t err;
  int dim = generate_random_bound_int(1,3,&iseed);
  int *dimensions;
  double *data;
  int data_size,idx;

  generate_output_file(&fid,&sid);

  dimensions = DANU_MALLOC(int,dim);
  data_size = 1;
  for(idx=0; idx < dim; idx++) {
    dimensions[idx] = generate_random_bound_int(10,1024,&iseed);
    data_size*=dimensions[idx];
  }
  data = DANU_MALLOC(double,data_size);

  generate_random_double_data(data,data_size,&iseed);


  err = data_write_double(sid,data_name,dim, dimensions, data);
  fail_unless(H5_RETURN_OK(err), 
	      "Failed to write non-series int data");


  DANU_FREE(dimensions);
  DANU_FREE(data);

  file_delete();

}
END_TEST
START_TEST (test_nsdata_write_float)
{
  hid_t fid, sid;
  char data_name[] = "Test Non-Series FLOAT";
  int iseed = 10;
  herr_t err;
  int dim = generate_random_bound_int(1,3,&iseed);
  int *dimensions;
  float *data;
  int data_size,idx;

  generate_output_file(&fid,&sid);

  dimensions = DANU_MALLOC(int,dim);
  data_size = 1;
  for(idx=0; idx < dim; idx++) {
    dimensions[idx] = generate_random_bound_int(10,1024,&iseed);
    data_size*=dimensions[idx];
  }
  data = DANU_MALLOC(float,data_size);

  generate_random_float_data(data,data_size,&iseed);


  err = data_write_float(sid,data_name,dim, dimensions, data);
  fail_unless(H5_RETURN_OK(err), 
	      "Failed to write non-series int data");


  DANU_FREE(dimensions);
  DANU_FREE(data);

  file_delete();

}
END_TEST

START_TEST (test_nsdata_read_int)
{
  hid_t fid, sid;
  char data_name[] = "Test Non-Series INT";
  int iseed = 2090;
  herr_t err;
  int dim = generate_random_bound_int(1,3,&iseed);
  int *dimensions;
  int *data, *read_data;
  int data_size,idx;

  generate_output_file(&fid,&sid);

  dimensions = DANU_MALLOC(int,dim);
  data_size = 1;
  for(idx=0; idx < dim; idx++) {
    dimensions[idx] = generate_random_bound_int(10,1024,&iseed);
    data_size*=dimensions[idx];
  }
  data = DANU_MALLOC(int,data_size);
  read_data = DANU_MALLOC(int,data_size);

  generate_random_int_data(data,data_size,&iseed);


  data_write_int(sid,data_name,dim, dimensions,data);

  output_file_close(&fid);
  
  open_file(&fid,&sid);

  data_read_int(sid,data_name,dim,dimensions,read_data);

  fail_unless(memcmp(read_data,data,sizeof(int)*data_size) == 0,
              "Failed to read correct INT data");


  DANU_FREE(dimensions);
  DANU_FREE(data);
  DANU_FREE(read_data);

  file_delete();

}
END_TEST

START_TEST (test_nsdata_read_double)
{
  hid_t fid, sid;
  char data_name[] = "Test Non-Series DOUBLE";
  int iseed = 2090;
  herr_t err;
  int dim = generate_random_bound_int(1,3,&iseed);
  int *dimensions;
  double *data, *read_data;
  int data_size,idx;

  generate_output_file(&fid,&sid);

  dimensions = DANU_MALLOC(int,dim);
  data_size = 1;
  for(idx=0; idx < dim; idx++) {
    dimensions[idx] = generate_random_bound_int(10,1024,&iseed);
    data_size*=dimensions[idx];
  }
  data = DANU_MALLOC(double,data_size);
  read_data = DANU_MALLOC(double,data_size);

  generate_random_double_data(data,data_size,&iseed);


  data_write_double(sid,data_name,dim, dimensions, data);

  err = data_read_double(sid,data_name,dim,dimensions,read_data);
  fail_unless(memcmp(read_data,data,sizeof(double)*data_size) == 0,
              "Failed to read correct DOUBLE data");

  DANU_FREE(dimensions);
  DANU_FREE(data);
  DANU_FREE(read_data);

  file_delete();

}
END_TEST
START_TEST (test_nsdata_read_float)
{
  hid_t fid, sid;
  char data_name[] = "Test Non-Series FLOAT";
  int iseed = 10;
  herr_t err;
  int dim = generate_random_bound_int(1,3,&iseed);
  int *dimensions;
  float *data, *read_data;
  int data_size,idx;

  generate_output_file(&fid,&sid);

  dimensions = DANU_MALLOC(int,dim);
  data_size = 1;
  for(idx=0; idx < dim; idx++) {
    dimensions[idx] = generate_random_bound_int(10,1024,&iseed);
    data_size*=dimensions[idx];
  }
  data = DANU_MALLOC(float,data_size);
  read_data = DANU_MALLOC(float,data_size);

  generate_random_float_data(data,data_size,&iseed);

  data_write_float(sid,data_name,dim, dimensions, data);

  err = data_read_float(sid,data_name,dim,dimensions,read_data);
  fail_unless(memcmp(read_data,data,sizeof(float)*data_size) == 0,
              "Failed to read correct FLOAT data");
 
  DANU_FREE(dimensions);
  DANU_FREE(data);
  DANU_FREE(read_data);

  file_delete();

}
END_TEST

START_TEST(test_nsdata_exists)
{
  hid_t fid, sid;
  herr_t err;
  char data_name[] = "Data DNE";
  int exists;
  int iseed = 0x44508;
  int data;
  int size;

  generate_output_file(&fid,&sid);

  err = data_exists(sid,data_name,&exists);
  fail_unless(H5_RETURN_OK(err),
              "Failed to return success status");
  fail_unless(exists == 0,
              "Failed to return correct exists flag (F)");

  data = generate_random_int(&iseed);
  size = 1;
  data_write_int(sid,data_name,1,&size,&data);

  err = data_exists(sid,data_name,&exists);
  fail_unless(H5_RETURN_OK(err),
              "Failed to return success status");
  fail_unless(exists != 0,
              "Failed to return correct exists flag (T)");

  danu_group_close(sid);
  output_file_close(&fid);

  file_delete();

}
END_TEST

START_TEST(test_nsdata_count)
{
  hid_t fid, sid;
  herr_t err;
  int iseed = 0xaa44;
  char data_name[32];
  int cnt;
  int gen_cnt, idx;
  int data;
  int size;

  generate_output_file(&fid,&sid);
  fail_unless(H5_ISA_VALID_ID(sid) && H5_ISA_VALID_ID(fid),
	      "failed to generate test file");

  printf("line=%d sid=0x%lx\n", __LINE__,sid);
  err = data_count(sid,&cnt);
  fail_unless(H5_RETURN_OK(err),
              "Failed to return success status");
  fail_unless(cnt == 0,
              "Failed to return correct count (0)");

  printf("line=%d sid=0x%lx\n", __LINE__,sid);
  data = generate_random_int(&iseed);
  size = 1;
  //gen_cnt = generate_random_bound_int(10,1024,&iseed);
  gen_cnt=1;
  for(idx=0; idx < gen_cnt; idx++) {
    sprintf(data_name, "Data Mu %04d",idx);
    data_write_int(sid,data_name,1,&size,&data);
  }

  printf("line=%d sid=0x%lx\n", __LINE__,sid);
  err = data_count(sid,&cnt);
  fail_unless(H5_RETURN_OK(err),
              "Failed to return success status");
  fail_unless(cnt == gen_cnt,
              "Failed to return correct count (<0)");

  printf("line=%d sid=0x%lx\n", __LINE__,sid);
  danu_group_close(sid);
  output_file_close(&fid);

  //file_delete();

}
END_TEST

START_TEST(test_nsdata_list)
{
  hid_t fid, sid;
  herr_t err;
  int iseed = 0xaa44;
  char **data_names;
  char **read_names;
  int gen_cnt, read_size, idx, num_found;
  int data;
  int size;

  generate_output_file(&fid,&sid);

  data = generate_random_int(&iseed);
  size = 1;

  gen_cnt = generate_random_bound_int(10,1024,&iseed);
  data_names = DANU_MALLOC(char *, gen_cnt);
  for(idx=0; idx < gen_cnt; idx++) {
    data_names[idx] = DANU_MALLOC(char, 32);  
    sprintf(data_names[idx], "Data Mu %04d",idx);
    data_write_int(sid,data_names[idx],1,&size,&data);
  }
  
  read_names = DANU_MALLOC(char *, gen_cnt);

  err = data_list(sid,gen_cnt,read_names,&num_found);
  fail_unless(H5_RETURN_OK(err),
              "Failed to return success status");

  fail_unless(num_found == gen_cnt,
              "Failed to find the correct number in data_list");

  fail_unless(char_array_cmp(read_names,data_names,gen_cnt) == 0,
              "Failed to read data names");

  
  DANU_FREE(read_names);
  read_size = gen_cnt / 2;
  read_names = DANU_MALLOC(char *, read_size);

  err = data_list(sid,read_size,read_names, &num_found);
  fail_unless(H5_RETURN_FAIL(err),
              "Failed to return fail status when array too small");

  fail_unless(num_found == gen_cnt,
              "Failed to find correct number in data_list");
               
  fail_unless(char_array_cmp(read_names,data_names,read_size) == 0,
              "Failed to read data names");

  

  danu_group_close(sid);
  output_file_close(&fid);

  file_delete();

}
END_TEST

Suite *
nsdata_suite (void)
{
    Suite *s = suite_create("Danu Non-Series Data");

#if 0
    TCase *nsdata_write = tcase_create("Non-Series Data Write");
    tcase_add_test(nsdata_write,test_nsdata_write_int);
    tcase_add_test(nsdata_write,test_nsdata_write_float);
    tcase_add_test(nsdata_write,test_nsdata_write_double);
    suite_add_tcase(s,nsdata_write);
#endif    
    TCase *nsdata_utils = tcase_create("Non-Series Data Utils");
//    tcase_add_test(nsdata_utils,test_nsdata_exists);
    tcase_add_test(nsdata_utils,test_nsdata_count);
//    tcase_add_test(nsdata_utils,test_nsdata_list);
    suite_add_tcase(s,nsdata_utils);
#if 0
    TCase *nsdata_read = tcase_create("Non-Series Data Read");
    tcase_add_test(nsdata_read,test_nsdata_read_int);
    tcase_add_test(nsdata_read,test_nsdata_read_float);
    tcase_add_test(nsdata_read,test_nsdata_read_double);
    suite_add_tcase(s,nsdata_read);
#endif
    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = nsdata_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_log(sr, "results.log");

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}



