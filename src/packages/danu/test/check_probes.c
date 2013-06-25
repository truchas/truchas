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

#include <danu_probes.h>

#define FILENAME "deleteme_probes.h5"
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

START_TEST (test_probe_create)
{
  hid_t fid, sid, pid;
  char probe_name[] = "Dummy Probe";
  herr_t err;
  int len;
  int num;
  double *data;
  hid_t type = H5T_NATIVE_DOUBLE;

  // This should fail
  len=4;
  num=10;
  err = probe_create_data(sid,probe_name,type,len,num,NULL,&pid);
  fail_unless(DANU_RETURN_FAIL(err),
	      "Failed to return fail status");

  // This should pass
  generate_output_file(&fid,&sid);
  data=DANU_MALLOC(double,len*num);
  err = probe_create_data(sid,probe_name,type,len,num,data,&pid);
  fail_unless(H5_ISA_VALID_ID(pid),
	      "Failed to return pass status");
  DANU_FREE(data);

  file_delete();
}
END_TEST

START_TEST ( test_probe_open )
{
  hid_t fid, sid, pid;
  char probe_name[] = "Dummy Probe";
  herr_t err;
  int len = 2;
  int num = 3;
  double *data;
  hid_t type = H5T_NATIVE_DOUBLE;

  // This should fail
  pid = probe_open_data(sid,probe_name);
  fail_unless(H5_ISA_INVALID_ID(pid),
	      "Failed to return correct status (F) 1");

  generate_output_file(&fid,&sid);

  // This should also fail since probe DNE
  pid = probe_open_data(sid,probe_name);
  fail_unless(H5_ISA_INVALID_ID(pid),
	      "Failed to return correct status (F) 2");

  // Generate the probe
  data = DANU_CALLOC(double,len*num);
  err = probe_create_data(sid,probe_name,type,len,num,data,&pid);
  fail_unless(H5_ISA_VALID_ID(pid),
	      "Failed to create dummy probe");
  danu_dataset_close(pid);

  pid = probe_open_data(sid,probe_name);
  fail_unless(H5_ISA_VALID_ID(pid),
	      "Failed to open existing probe");


  file_delete();

}
END_TEST

START_TEST(test_probe_count)
{
  hid_t fid, sid, pid;
  char probe_base_name[] = "Dummy Probe";
  char probe_name[64];
  herr_t err;
  int i,cnt,num;
  int iseed = 0x1234FFFF;
  int len = 1;
  hid_t type = H5T_NATIVE_INT;
  int * data;

  // This should fail
  err = probe_count(sid,&cnt);
  fail_unless(H5_RETURN_FAIL(err),
	      "Failed to return correct status (F)");

  generate_output_file(&fid,&sid);

  err = probe_count(sid,&cnt);
  fail_unless(H5_RETURN_OK(err),
	      "Fail to return correct status (OK)");
  fail_unless(cnt == 0,
	      "Fail to return correct number 0");


  // Generate dummy probe datasets
  num = generate_random_bound_int(32,128,&iseed);
  for(i=0;i<num; i ++ ) {
    data = DANU_CALLOC(int,len*num);
    sprintf(probe_name,"%s %d", probe_base_name,i);
    err = probe_create_data(sid,probe_name,type,len,num,data,&pid);
    DANU_FREE(data);
    danu_dataset_close(pid);
  }

  // This should pass
  err = probe_count(sid,&cnt);
  fail_unless(H5_RETURN_OK(err),
	      "Fail to return correct status (OK)");
  fail_unless(cnt == num,
	      "Fail to return correct number num");

  file_delete();

}
END_TEST

START_TEST(test_probe_data_int)
{
  hid_t fid, sid, pid;
  char probe_name[64];
  herr_t err;
  int i,cnt;
  int iseed = 0x1234FFFF;
  int len = generate_random_bound_int(2,10,&iseed);
  int num = generate_random_bound_int(10,50,&iseed);
  int array[num][len];
  int read_array[num][len];
  int *data, *rdata;

  // These calls should all fail
  data=DANU_CALLOC(int,num*len);
  err = probe_data_write_int(pid,num,data);
  fail_unless(H5_RETURN_FAIL(err),
	      "Failed to return correct status (F)");
  
  err = probe_data_write_int(pid,num,data);
  fail_unless(H5_RETURN_FAIL(err),
	      "Failed to return correct status (F)");


  err = probe_data_write_int(pid,-2,data);
  fail_unless(H5_RETURN_FAIL(err),
	      "Failed to return correct status (F)");

  err = probe_data_write_int(pid,num,NULL);
  fail_unless(H5_RETURN_FAIL(err),
	      "Failed to return correct status (F)");
  DANU_FREE(data);

  /* Create the file and data */
  generate_output_file(&fid,&sid);
  sprintf(probe_name,"%s", "Integer Probe Data");
  for(cnt=0; cnt < num; cnt++) {
    for(i=0; i<len; i++ ) {
      array[cnt][i]=cnt+1;
    }
  }

  data = &array[0][0];
  err = probe_create_data_int(sid,probe_name,len,num,data,&pid);
  fail_unless(DANU_RETURN_OK(err),
              "Failed to create and write probe dataset");
  danu_dataset_close(pid);
  danu_file_close(fid);

  /* Reopen and read data */
  open_file(&fid,&sid);
  rdata=&read_array[0][0];
  pid = probe_open_data(sid,probe_name);
  err = probe_data_read_int(pid,rdata);
  fail_unless(DANU_RETURN_OK(err),
              "Failed to read probe data");

  /* Now check */
  data = &array[0][0];
  rdata = &read_array[0][0];
  fail_unless(memcmp(data,rdata,(len*num)*sizeof(int)) == 0,
              "Read in data does not match written data");

  file_delete();

}
END_TEST

START_TEST(test_probe_data_double)
{
  hid_t fid, sid, pid;
  char probe_name[64];
  herr_t err;
  int i,cnt;
  int iseed = 0x1234FFFF;
  int len = generate_random_bound_int(2,10,&iseed);
  int num = generate_random_bound_int(10,50,&iseed);
  double array[num][len];
  double read_array[num][len];
  double *data, *rdata;

  /* Create the file and data */
  generate_output_file(&fid,&sid);
  sprintf(probe_name,"%s", "Double Probe Data");
  for(cnt=0; cnt < num; cnt++) {
    for(i=0; i<len; i++ ) {
      array[cnt][i]= (cnt+1)*10.1;
    }
  }

  data = &array[0][0];
  err = probe_create_data_double(sid,probe_name,len,num,data,&pid);
  fail_unless(DANU_RETURN_OK(err),
              "Failed to write double probe data");

  danu_file_close(fid);

  /* Reopen and read data */
  open_file(&fid,&sid);
  rdata=&read_array[0][0];
  pid = probe_open_data(sid,probe_name);
  err = probe_data_read_double(pid,rdata);
  fail_unless(DANU_RETURN_OK(err),
              "Failed to read probe data");

  /* Now check */
  data = &array[0][0];
  rdata = &read_array[0][0];
  fail_unless(memcmp(data,rdata,(len*num)*sizeof(double)) == 0,
              "Read in data does not match written data");


  file_delete();
}
END_TEST

START_TEST(test_probe_data_float)
{
  hid_t fid, sid, pid;
  char probe_name[64];
  herr_t err;
  int i,cnt;
  int iseed = 0x1234FFFF;
  int len = generate_random_bound_int(2,10,&iseed);
  int num = generate_random_bound_int(10,50,&iseed);
  float array[num][len];
  float read_array[num][len];
  float *data, *rdata;

  /* Create the file and data */
  generate_output_file(&fid,&sid);
  sprintf(probe_name,"%s", "Float Probe Data");
  for(cnt=0; cnt < num; cnt++) {
    for(i=0; i<len; i++ ) {
      array[cnt][i]= (cnt+1)*10.1;
    }
  }

  data = &array[0][0];
  err = probe_create_data_float(sid,probe_name,len,num,data,&pid);
  fail_unless(DANU_RETURN_OK(err),
              "Failed to write float probe data");

  danu_file_close(fid);

  /* Reopen and read data */
  open_file(&fid,&sid);
  rdata=&read_array[0][0];
  pid = probe_open_data(sid,probe_name);
  err = probe_data_read_float(pid,rdata);
  fail_unless(DANU_RETURN_OK(err),
              "Failed to read probe data");

  /* Now check */
  data = &array[0][0];
  rdata = &read_array[0][0];
  fail_unless(memcmp(data,rdata,(len*num)*sizeof(float)) == 0,
              "Read in data does not match written data");


  file_delete();
}
END_TEST

START_TEST(test_probe_data_append)
{
  hid_t fid, sid, pid;
  char probe_name[64];
  herr_t err;
  int i,cnt;
  int iseed = 0x1234FFFF;
  int len = generate_random_bound_int(2,10,&iseed);
  int num = generate_random_bound_int(10,50,&iseed);
  int new_num = generate_random_bound_int(1,5,&iseed);
  double array[num][len], new_array[new_num][len];
  double read_array[num+new_num][len];
  double *data, *rdata;

  /* Create the file and data */
  generate_output_file(&fid,&sid);
  sprintf(probe_name,"%s", "Double Probe Append Data");
  for(cnt=0; cnt < num; cnt++) {
    for(i=0; i<len; i++ ) {
      array[cnt][i]= (cnt+1)*10.1;
    }
  }

  data = &array[0][0];
  err = probe_create_data_double(sid,probe_name,len,num,data,&pid);
  fail_unless(DANU_RETURN_OK(err),
              "Failed to write double probe data");

  /* Append data */
  data = &new_array[0][0];
  for(i=0;i<new_num*len;i++) {
      *data=-9999.9;
      data++;
  }
  data=&new_array[0][0];
  err= probe_data_write_double(pid,new_num,data);
  fail_unless(DANU_RETURN_OK(err),
              "Failed to append double data");

  danu_file_close(fid);

  /* Reopen and read data */
  open_file(&fid,&sid);
  pid=probe_open_data(sid,probe_name);
  rdata=&read_array[0][0];
  err = probe_data_read_double(pid,rdata);
  fail_unless(DANU_RETURN_OK(err),
              "Failed to read probe data");

  /* Now check */
  data = &array[0][0];
  rdata = &read_array[0][0];
  fail_unless(memcmp(data,rdata,(len*num)*sizeof(double)) == 0,
              "Read in data does not match written data");

  data = &new_array[0][0];
  rdata = &read_array[num+new_num-1][0];
  fail_unless(memcmp(data,rdata,(len*new_num)*sizeof(double)) == 0,
              "Read in data does not match append data");


  file_delete();
}
END_TEST

Suite *
probe_suite (void)
{
    Suite *s = suite_create("Danu Probe Data");

    TCase *probe_utils = tcase_create("Probe Utils");
    tcase_add_test(probe_utils,test_probe_create);
    tcase_add_test(probe_utils,test_probe_open);
    tcase_add_test(probe_utils,test_probe_count);
    //tcase_add_test(probe_utils,test_probe_list);
    suite_add_tcase(s,probe_utils);


    TCase *probe_data = tcase_create("Probe Data");
    tcase_add_test(probe_data,test_probe_data_int);
    tcase_add_test(probe_data,test_probe_data_double);
    tcase_add_test(probe_data,test_probe_data_float);
    tcase_add_test(probe_data,test_probe_data_append);

    suite_add_tcase(s,probe_data);

    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = probe_suite ();
  SRunner *sr = srunner_create (s);

  srunner_set_log(sr, "probe-test-results.log");

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}



