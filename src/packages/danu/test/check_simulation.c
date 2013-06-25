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

#include <danu_sim.h>

#define FILENAME      "deleteme-simulation.h5"
#define MESHNAME      "Test Mesh"
#define SIMULATION_ROOT_NAME "Simulation Number"


/* BEGIN TESTS */
START_TEST (test_sim_add)
{
  const char test_file[] = FILENAME;
  const char sim_name[]  = SIMULATION_ROOT_NAME;
  hid_t  fid,sid;
  herr_t status;
  int exists;

  status = simulation_add(fid,sim_name,&sid);
  fail_unless(H5_RETURN_FAIL(status),
              "Failed to return fail status with bad HID");

  output_file_create(test_file, &fid);

  status = simulation_add(fid,sim_name,&sid);
  fail_unless(H5_RETURN_OK(status),
              "failed to add simulation group");

  output_file_close(&fid);
  danu_file_delete(test_file);

}
END_TEST

START_TEST (test_sim_count)
{
  const char test_file[] = FILENAME;
  char sim_name[32];
  int idx;
  int iseed = 0x0FFF123;
  int num = generate_random_bound_int(1,64,&iseed);
  hid_t  fid, sid;
  int count;
  herr_t status;

  status = simulation_count(fid,&count);
  fail_unless(H5_RETURN_FAIL(status),
              "Failed to return fail status with bad HID");

  output_file_create(test_file,&fid);
  
  status = simulation_count(fid,&count);
  fail_unless(H5_RETURN_OK(status),
              "Failed to return valid status");

  fail_unless(count == 0,
              "Failed to return correct count (0)");

  for(idx = 1; idx <=num; idx++) {
   sprintf(sim_name,"%s %04d", SIMULATION_ROOT_NAME,idx);
   simulation_add(fid,sim_name,&sid);
   danu_group_close(sid);
  }
  
  status = simulation_count(fid,&count);
  fail_unless(H5_RETURN_OK(status),
              "Failed to return valid status");
  fail_unless(count == num,
              "Failed to return correct count (>0)");

  output_file_close(&fid);
  danu_file_delete(test_file);

}
END_TEST

START_TEST (test_sim_exists)
{
    const char test_file[] = FILENAME;
    const char sim_name[]  = SIMULATION_ROOT_NAME;
    hid_t  fid, sid;
    herr_t status;
    int exists;

    status = simulation_exists(fid,sim_name,&exists);
    fail_unless( H5_RETURN_FAIL(status),
                 "Failed to return fail status with bad HID");
    
    output_file_create(test_file, &fid);

    status = simulation_exists(fid,sim_name,&exists);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to return valid status");

    fail_unless(! exists,
                "Failed to return correct status (FALSE)");

    simulation_add(fid,sim_name,&sid);
    danu_group_close(sid);

    status = simulation_exists(fid,sim_name,&exists);
    fail_unless(H5_RETURN_OK(status),
                "Failed to return valid status");

    fail_unless(exists,
                "Failed to return correct status (TRUE)");

    

    output_file_close(&fid);
    danu_file_delete(test_file);

}
END_TEST

START_TEST (test_sim_open)
{
    const char test_file[] = FILENAME;
    const char sim_name[]  = SIMULATION_ROOT_NAME;
    hid_t  fid, sid;
    herr_t status;

    status = simulation_open(fid,sim_name,&sid);
    fail_unless( H5_RETURN_FAIL(status),
                 "Failed to return fail status with bad HID");
    
    output_file_create(test_file, &fid);

    status = simulation_open(fid,sim_name,&sid);
    fail_unless( H5_RETURN_FAIL(status),
                 "Failed to return valid status");
    fail_unless(H5_ISA_INVALID_ID(sid),
                "Failed to return bad HID when sim DNE");

    simulation_add(fid,sim_name,&sid);
    danu_group_close(sid);

    status = simulation_open(fid,sim_name,&sid);
    fail_unless(H5_RETURN_OK(status),
                "Failed to return valid status");

    fail_unless(H5_ISA_VALID_ID(sid),
                "Failed to return valid HID");

    output_file_close(&fid);
    danu_file_delete(test_file);

}
END_TEST

START_TEST (test_sim_list)
{
  const char test_file[] = FILENAME;
  char sim_name[32];
  int idx;
  int iseed = 0x0FFF123;
  int num = generate_random_bound_int(2,64,&iseed);
  hid_t  fid, sid;
  int num_found;
  int num_ok;
  size_t len;
  char **read_names;
  herr_t status;

  read_names = DANU_MALLOC(char *, num);

  status = simulation_list(fid,num, read_names, &num_found);
  fail_unless(H5_RETURN_FAIL(status),
              "Failed to return fail status with bad HID");

  output_file_create(test_file,&fid);
  
  status = simulation_list(fid,num, read_names, &num_found);
  fail_unless(H5_RETURN_OK(status),
              "Failed to return valid status");
  fail_unless(num_found == 0,
              "failed to return correct num_found (0)");

  /* Create the simulations */
  for(idx = 1; idx <=num; idx++) {
   sprintf(sim_name,"%s %04d", SIMULATION_ROOT_NAME,idx);
   simulation_add(fid,sim_name,&sid);
   danu_group_close(sid);
  }

  status = simulation_list(fid,num,read_names,&num_found);
  fail_unless(H5_RETURN_OK(status),
              "Failed to return valid status");
  fail_unless(num_found == num,
              "Failed to return correct num_found (>0)");
  num_ok = 0;
  for(idx=1; idx <= num; idx++) {
    sprintf(sim_name,"%s %04d", SIMULATION_ROOT_NAME,idx);
    num_ok += strcmp(sim_name, read_names[idx-1]);
    if ( num_ok != 0 ) {
      printf("%s %s num_ok=%d\n", sim_name, read_names[idx-1]);
    }
    DANU_FREE(read_names[idx-1]);
  }

  fail_unless(num_ok == 0,
              "Failed to return correct simulation names");

  status = simulation_list(fid,num-1,read_names,&num_found);
  fail_unless(H5_RETURN_OK(status),
              "Failed to return fail status");
  num_ok = 0;
  for(idx=1;idx<=num-1; idx++) {
    sprintf(sim_name,"%s %04d", SIMULATION_ROOT_NAME,idx);
    num_ok += strcmp(sim_name, read_names[idx-1]);
    if ( num_ok != 0 ) {
      printf("%s %s num_ok=%d\n", sim_name, read_names[idx-1],num_ok);
    }
    DANU_FREE(read_names[idx-1]);
  }

  fail_unless(num_ok == 0,
              "Failed to return correct simulation names (num<num_found)");
  
  read_names = DANU_REALLOC(char*,read_names, 2*num);
  status = simulation_list(fid,2*num,read_names,&num_found);
  fail_unless(H5_RETURN_OK(status),
              "Failed to return fail status");
  num_ok = 0;
  for(idx=1;idx<=num_found; idx++) {
    sprintf(sim_name,"%s %04d", SIMULATION_ROOT_NAME,idx);
    num_ok += strcmp(sim_name, read_names[idx-1]);
    if ( num_ok != 0 ) {
      printf("<%s> <%s> num_ok=%d\n", sim_name, read_names[idx-1],num_ok);
    }
    DANU_FREE(read_names[idx-1]);
  }

  for(idx=num_found+1;idx<=2*num;idx++) {
    if ( ! DANU_BAD_STRING(read_names[idx-1]) )  {
      printf("Found a non-NULL pointer %d 0x%lx <%s>\n",idx-1,read_names[idx-1],read_names[idx-1]);
      num_ok++;
    }
    DANU_FREE(read_names[idx-1]);
  }

  fail_unless(num_ok == 0,
              "Failed to return correct simulation names (num>num_found)");

   
  output_file_close(&fid);
  danu_file_delete(test_file);

  DANU_FREE(read_names);

}
END_TEST

START_TEST (test_sim_link_mesh)
{
  const char test_file[] = FILENAME;
  char sim_name[] = SIMULATION_ROOT_NAME;
  char mesh_name[] = "Test Mesh";
  hid_t  fid, sid, mesh;
  herr_t status;

  status = simulation_link_mesh(fid,sid,mesh_name);
  fail_unless(H5_RETURN_FAIL(status),
              "Failed to return fail status with invalid file HID");

  output_file_create(test_file,&fid);
  
  status = simulation_link_mesh(fid,sid,mesh_name);
  fail_unless(H5_RETURN_FAIL(status),
              "Failed to return fail status with invalid sim HID");
  
  simulation_add(fid,sim_name,&sid);

  status = simulation_link_mesh(fid,sid,mesh_name);
  fail_unless(H5_RETURN_OK(status),
              "Failed to return valid status with mesh link");

  danu_group_close(sid);
  output_file_close(&fid);
  danu_file_delete(test_file);

}
END_TEST

Suite *
sim_suite (void)
{
    Suite *s = suite_create("Danu Simulation");

    TCase *sim_create = tcase_create("Simulation Create");
    tcase_add_test(sim_create,test_sim_add);
    suite_add_tcase(s,sim_create);

    TCase *sim_open = tcase_create("Simulation Open");
    tcase_add_test(sim_open,test_sim_open);
    suite_add_tcase(s,sim_open);


    TCase *sim_utils = tcase_create("Simulation Utils");
    tcase_add_test(sim_utils,test_sim_exists);
    tcase_add_test(sim_utils,test_sim_count);
    tcase_add_test(sim_utils,test_sim_list);
    tcase_add_test(sim_utils,test_sim_link_mesh);
    suite_add_tcase(s,sim_utils);


    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = sim_suite ();
  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
