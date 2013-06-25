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
#include "random_generators.h"

#include <danu_types.h>
#include <danu_h5_error.h>
#include <danu_memory.h>
#include <danu_error.h>
#include <danu_utils.h>
#include <danu_file.h>
#include <danu_dataset_types.h>

#include <danu_output.h>
#include <danu_sim.h>

#include <danu_series-group.h>
#include <danu_series-data.h>

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

START_TEST (test_data_exists)
{
  hid_t fid, sid, nid;
  herr_t err;
  int iseed = 0xeeee4444;
  double time = generate_random_double(&iseed);
  int cycle = generate_random_bound_int(1,1000,&iseed);
  char data_name[] = "Pressure";
  int flag;

  /* Setup for the test */
  generate_output_file(&fid,&sid);
  sequence_getNextID(sid,cycle,time,&nid);

  err = simulation_data_exists(nid,data_name,&flag);
  fail_unless(H5_RETURN_OK(err),
	      "Failed to return correctly");
  fail_unless(flag == 0,
	      "failed to return correct flag(0)");



  file_delete();

}
END_TEST

START_TEST (test_data_count)
{
  hid_t fid, sid, nid;
  herr_t err;
  int iseed = 0xeeee4444;
  double time = generate_random_double(&iseed);
  int cycle = generate_random_bound_int(0,1000,&iseed);
  int count;

  /* Setup for the test */
  generate_output_file(&fid,&sid);
  sequence_getNextID(sid,cycle,time,&nid);

  /* Zero count test */
  err = simulation_data_count(nid,&count);
  fail_unless(H5_RETURN_OK(err),
	      "Failed to return correctly");
  fail_unless(count == 0,
	      "Failed to return correct flag(0)");


  file_delete();

}
END_TEST

START_TEST (test_data_list)
{
  hid_t fid, sid, nid;
  herr_t err;
  int iseed = 0xeeee4444;
  double time = generate_random_double(&iseed);
  int cycle=generate_random_int(&iseed);
  int count;
 
  /* Setup for the test */
  generate_output_file(&fid,&sid);
  sequence_getNextID(sid,cycle,time,&nid);

  file_delete();

}
END_TEST

START_TEST(test_data_description)
{
  hid_t fid, sid, nid;
  herr_t stat;
  int iseed = 0xFFFFF;
  double time = generate_random_double(&iseed);
  int cycle = generate_random_int(&iseed);
  int ndims = generate_random_bound_int(1,3,&iseed);
  int *dims;
  int *idata;
  double *data8;
  float *data4;
  int idx, num;
  hsize_t *read_dims;
  int rank, read_code;

  /* Setup for the test */
  generate_output_file(&fid,&sid);
  sequence_getNextID(sid,cycle,time,&nid);

  dims=DANU_MALLOC(int,ndims);
  num=1;
  for(idx=0; idx<ndims; idx++) {
    dims[idx]=generate_random_bound_int(10,20,&iseed);
    num*=dims[idx];
  }

  idata=DANU_MALLOC(int,num);
  data8=DANU_MALLOC(double,num);
  data4=DANU_MALLOC(float,num);
  for(idx=0;idx<num;idx++) {
    idata[idx]=generate_random_int(&iseed);
    data8[idx]=generate_random_double(&iseed);
    data4[idx]=generate_random_float(&iseed);
  }

  simulation_data_write_int(nid, "Integer Data", ndims, dims, idata);
  simulation_data_write_double(nid, "Double Data", ndims, dims, data8);
  simulation_data_write_float(nid, "Float Data", ndims, dims, data4);

  /* Check the integer dataset size and rank */
  rank=simulation_data_rank(nid,"Integer Data");
  fail_unless(rank == ndims,
              "Returned incorrect rank value");
   
  read_dims=DANU_MALLOC(hsize_t,ndims);
  stat=simulation_data_dimensions(nid, "Integer Data", ndims, read_dims);

  fail_unless( DANU_RETURN_OK(stat),
               "Failed to find data set size");

  for(idx=0;idx<ndims;idx++) {
    fail_unless(read_dims[idx] == dims[idx],
                "Failed to rad correct data size");
  }

  /* Now check the type codes */
  stat = simulation_data_type(nid,"Integer Data",&read_code);
  fail_unless( DANU_RETURN_OK(stat),
               "Failed to read data type");
  fail_unless(read_code == DANU_DATASET_INT,
              "Failed to read correct data type");


  stat = simulation_data_type(nid,"Double Data",&read_code);
  fail_unless( DANU_RETURN_OK(stat),
               "Failed to read data type");
  fail_unless(read_code == DANU_DATASET_DOUBLE,
              "Failed to read correct data type");

  stat = simulation_data_type(nid,"Float Data",&read_code);
  fail_unless( DANU_RETURN_OK(stat),
               "Failed to read data type");
  fail_unless(read_code == DANU_DATASET_FLOAT,
              "Failed to read correct data type");
  

  /* Release memory */
  DANU_FREE(idata);
  DANU_FREE(dims);
  DANU_FREE(read_dims);
  DANU_FREE(data8);
  DANU_FREE(data4);

  file_delete(); 

}
END_TEST
Suite *
data_suite (void)
{
    Suite *s = suite_create("Danu Series Data");

    TCase *data_utils = tcase_create("Series Data Utils");
    tcase_add_test(data_utils,test_data_exists);
    tcase_add_test(data_utils,test_data_count);
    tcase_add_test(data_utils,test_data_list);
    tcase_add_test(data_utils,test_data_description);
    suite_add_tcase(s,data_utils);

    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = data_suite();
  SRunner *sr = srunner_create (s);

  srunner_set_log(sr, "series-data-results.log");

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}



