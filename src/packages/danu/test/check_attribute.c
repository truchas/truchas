/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
/* ****************************************************************************
*                                                                             *
*  Danu Unit Test                                                             *
*  Unit: Attribute                                                            *
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

#include <danu_attribute.h>

#define INT_ATTR_NAME    "Integer Attribute"
#define DOUBLE_ATTR_NAME "Double Attribute"
#define FLOAT_ATTR_NAME  "Float Attribute"
#define UINT_ATTR_NAME   "Unisigned Integer Attribute"
#define STRING_ATTR_NAME "String Attribute"

#define INT_ATTR_VALUE    0xF1234568
#define DOUBLE_ATTR_VALUE 1.98765432e-33
#define FLOAT_ATTR_VALUE  0.1234567891234
#define UINT_ATTR_VALUE   0x01012020
#define STRING_ATTR_VALUE "The quick brown fox jumped over the lazy dog"

START_TEST (test_attr_error_badid)
{
    hid_t  dummy_id;
    herr_t status;
    int data = INT_ATTR_VALUE;
    int result;
    char attr_name[] = INT_ATTR_NAME;


    status = danu_attr_write(dummy_id, attr_name, &data, H5T_NATIVE_INT);
    fail_unless(H5_RETURN_FAIL(status), 
                "danu_attr_write failed to return an error with invalid id");

    status = danu_attr_read(dummy_id, attr_name, &result, H5T_NATIVE_INT);
    fail_unless(H5_RETURN_FAIL(status), 
                "danu_attr_read failed to return an error with invalid id");
}
END_TEST

START_TEST(test_attr_error_badtype)
{
    hid_t fid;
    hid_t type = TEST_HID;
    herr_t status;
    int data = INT_ATTR_VALUE;
    int result;
    char attr_name[] = INT_ATTR_NAME;

    /* Create a test file so the loc id is valid */
    create_test_h5_file();
    fid = open_test_file();

    status = danu_attr_write(fid, attr_name, &data, type);
    fail_unless(H5_RETURN_FAIL(status), 
                "danu_attr_write failed to return an error with invalid type");

    status = danu_attr_read(fid, attr_name, &result, type);
    fail_unless(H5_RETURN_FAIL(status), 
                "danu_attr_read failed to return an error with invalid type");

    close_test_file(fid);
    delete_test_file();
}
END_TEST

START_TEST (test_attr_error_baddata)
{
    hid_t fid;
    herr_t status;
    char attr_name[] = INT_ATTR_NAME;

    /* Create a test file */
    create_test_h5_file();
    fid = open_test_file();

    /* Bad data pointers */
    status = danu_attr_write(fid,attr_name,NULL,H5T_NATIVE_INT);
    fail_unless( H5_RETURN_FAIL(status),
                 "danu_attr_write Failed to return an error with NULL data pointer" );

    status = danu_attr_read(fid,attr_name,NULL,H5T_NATIVE_INT);
    fail_unless( H5_RETURN_FAIL(status),
                 "danu_attr_read Failed to return an error with NULL data pointer" );

    close_test_file(fid);
    delete_test_file();
}
END_TEST

START_TEST (test_attr_error_badname)
{
    hid_t fid;
    herr_t status;
    char attr_name[] = "This attribute does not exist";
    char dummy = '\0';
    int data = INT_ATTR_VALUE;
    int result;

    /* Create test file */
    create_test_h5_file();
    fid = open_test_file();

    status = danu_attr_write(fid,NULL,&data,H5T_NATIVE_INT);
    fail_unless( H5_RETURN_FAIL(status),
                        "danu_attr_write Failed to return fail status with null name pointer");

    status = danu_attr_write(fid,&dummy,&data,H5T_NATIVE_INT);
    fail_unless( H5_RETURN_FAIL(status),
                 "danu_attr_write Failed to return fail status with empty name pointer");

    status = danu_attr_read(fid,NULL,&result,H5T_NATIVE_INT);
    fail_unless( H5_RETURN_FAIL(status),
                 "danu_attr_read Failed to return fail status with null name pointer");

    status = danu_attr_read(fid,&dummy,&result,H5T_NATIVE_INT);
    fail_unless( H5_RETURN_FAIL(status),
                 "danu_attr_read Failed to return fail status with empty name pointer");

    status = danu_attr_read(fid,attr_name,&result,H5T_NATIVE_INT);
    fail_unless( H5_RETURN_FAIL(status),
                 "danu_attr_read Failed to return fail status for attribute that is DNE");


    close_test_file(fid);
    delete_test_file();
}
END_TEST

START_TEST (test_attr_error_count)
{
  hid_t fid;
  herr_t status;
  int num;

  status = danu_attr_count(fid, &num);
  fail_unless(H5_RETURN_FAIL(status),
              "danu_attr_count did not fail with bad id");

  create_test_h5_file();
  fid = open_test_file();

  status = danu_attr_count(fid,NULL);
  fail_unless(H5_RETURN_FAIL(status),
              "danu_attr_count did not fail with NULL pointer");

  close_test_file(fid);
  delete_test_file();
}
END_TEST

START_TEST (test_attr_error_names)
{
  hid_t fid;
  herr_t status;
  int num = 1;
  char **names;
  int num_found;

  names = (char **)malloc(sizeof(char*)*16);

  status = danu_attr_names(fid,num,names,&num_found);
  fail_unless( H5_RETURN_FAIL(status),
               "danu_attr_names failed to return an error with bad id");

  fid = create_open_test_h5_file();

  status = danu_attr_names(fid,-1,names,&num_found);
  fail_unless( H5_RETURN_FAIL(status),
               "danu_attr_names failed to return an error with bad num value");

  status = danu_attr_names(fid,num,NULL,&num_found);
  fail_unless( H5_RETURN_FAIL(status),
               "danu_attr_names failed to return an error with bad names");

  status = danu_attr_names(fid,num,names,NULL);
  fail_unless( H5_RETURN_FAIL(status),
               "danu_attr_names failed to return an error with bad num_found");


  close_delete_test_file(fid);
}
END_TEST

START_TEST (test_attr_exists)
{

  hid_t fid;
  char attr_name[] = INT_ATTR_NAME;
  int data = INT_ATTR_VALUE;

  int exists;
  herr_t status;

  /* Pass a bad HID */
  status = danu_attr_exists(fid,attr_name,&exists);
  fail_unless( H5_RETURN_FAIL(status),
	       "Failed to raise error on bad id");

  /* Create test file */
  fid = create_open_test_h5_file();

  /* Another fail */
  status = danu_attr_exists(fid,NULL,&exists);
  fail_unless( H5_RETURN_FAIL(status),
	       "Failed to raise error on bad pointer");

  /* And another fail .. */
  status = danu_attr_exists(fid,attr_name,&exists);
  fail_unless( H5_RETURN_OK(status),
	       "Unexpected return condition");
  fail_unless(exists == 0,
	      "Failed to return corrected flag expected FALSE");

  status = danu_attr_write_int(fid,attr_name,data);
  fail_unless(H5_RETURN_OK(status),
	      "Failed to write dummy attribute in exist test");

  close_test_file(fid);

  /* Now THIS should pass */
  fid = open_test_file();
  status = danu_attr_exists(fid,attr_name,&exists);
  fail_unless( H5_RETURN_OK(status),
	       "Unexpected return condition");
  fail_unless(exists > 0,
	      "Failed to return corrected flag expected TRUE");

	       
  /* Remove the test file */
  close_delete_test_file(fid);

}
END_TEST

START_TEST (test_attr_int)
{
    hid_t fid;
    char attr_name[] = INT_ATTR_NAME;
    int data = INT_ATTR_VALUE;
    int result;
    herr_t status;
    char type;

    /* Create the test file */
    create_test_h5_file();
    fid = open_test_file();

    /* Write the integer value */
    status = danu_attr_write(fid,attr_name,&data,H5T_NATIVE_INT);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to write an int attribute");
    
    /* Now read the value after closing the file....this forces HDF5 to flush memory */
    close_test_file(fid);
    fid = open_test_file();
    status = danu_attr_read(fid,attr_name,&result,H5T_NATIVE_INT);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to read int attribute");

    fail_unless( data == result,
                 "Data written and does not match read result");

    type = danu_attr_get_type(fid,attr_name);
    fail_unless(type == 'i',
	        "Failed to return correct datatype (int)");

    close_test_file(fid);
    delete_test_file();

}
END_TEST

START_TEST (test_attr_int_wrapper)
{
    hid_t fid;
    char attr_name[] = INT_ATTR_NAME;
    int data = INT_ATTR_VALUE;
    int result;
    herr_t status;

    /* Create the file */
    create_test_h5_file();
    fid = open_test_file();

    /* Test the wrappers */
    status = danu_attr_write_int(fid,attr_name,data);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to write int attribute through wrapper");

    close_test_file(fid);
    fid = open_test_file();
    status = danu_attr_read_int(fid,attr_name,&result);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to read int attribute through wrapper");

    fail_unless( data == result,
                 "Data written through wrapper did not match wrapper read result");

    close_test_file(fid);
    delete_test_file();

}
END_TEST

START_TEST (test_attr_double)
{
    hid_t fid;
    char attr_name[] = DOUBLE_ATTR_NAME;
    double data = DOUBLE_ATTR_VALUE;
    double result;
    char type;
    herr_t status;

    /* Create the test file */
    create_test_h5_file();
    fid = open_test_file();

    /* Write the double value */
    status = danu_attr_write(fid,attr_name,&data,H5T_NATIVE_DOUBLE);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to write an double attribute");
    
    /* Now read the value after closing the file....this forces HDF5 to flush memory */
    close_test_file(fid);
    fid = open_test_file();
    status = danu_attr_read(fid,attr_name,&result,H5T_NATIVE_DOUBLE);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to read double attribute");

    fail_unless( data == result,
                 "Data written and does not match read result");

    type = danu_attr_get_type(fid,attr_name);
    fail_unless(type == 'd',
	        "Failed to return correct datatype (double)");

    close_test_file(fid);
    delete_test_file();

}
END_TEST

START_TEST (test_attr_overwrite)
{
  hid_t fid = create_open_test_h5_file();
  int iseed = 0x34AAEEFF;
  int data = generate_random_int(&iseed);;
  int data_read;
  char attr_name[] = INT_ATTR_NAME;
  herr_t status;

  status = danu_attr_write_int(fid,attr_name,data);
  fail_unless ( H5_RETURN_OK(status),
		"Unexpected failed return");

  /* Now close and reopen  the file */
  close_test_file(fid);
  fid = open_test_file();

  /* Overwrite the data */
  data = generate_random_int(&iseed);
  status = danu_attr_write_int(fid,attr_name,data);
  fail_unless ( H5_RETURN_OK(status),
		"Unexpected fail return when overwriting attribute");

  /* reopen again to read */
  close_test_file(fid);
  fid = open_test_file();

  status = danu_attr_read_int(fid,attr_name, &data_read);
  fail_unless ( H5_RETURN_OK(status),
		"Unexpected fail return when reading attribute");
  fail_unless( data == data_read,
	       "Failed to read the new data value");

  close_delete_test_file(fid);
}
END_TEST

START_TEST (test_attr_double_wrapper)
{
    hid_t fid;
    char attr_name[] = DOUBLE_ATTR_NAME;
    double data = DOUBLE_ATTR_VALUE;
    double result;
    herr_t status;

    /* Create the file */
    create_test_h5_file();
    fid = open_test_file();

    /* Test the wrappers */
    status = danu_attr_write_double(fid,attr_name,data);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to write double attribute through wrapper");

    close_test_file(fid);
    fid = open_test_file();
    status = danu_attr_read_double(fid,attr_name,&result);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to read double attribute through wrapper");

    fail_unless( data == result,
                 "Data written through wrapper did not match wrapper read result");

    close_test_file(fid);
    delete_test_file();

}
END_TEST
START_TEST (test_attr_float)
{
    hid_t fid;
    char attr_name[] = FLOAT_ATTR_NAME;
    float data = FLOAT_ATTR_VALUE;
    float result;
    herr_t status;

    /* Create the test file */
    create_test_h5_file();
    fid = open_test_file();

    /* Write the floateger value */
    status = danu_attr_write(fid,attr_name,&data,H5T_NATIVE_FLOAT);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to write an float attribute");
    
    /* Now read the value after closing the file....this forces HDF5 to flush memory */
    close_test_file(fid);
    fid = open_test_file();
    status = danu_attr_read(fid,attr_name,&result,H5T_NATIVE_FLOAT);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to read float attribute");

    fail_unless( data == result,
                 "Data written and does not match read result");

    close_test_file(fid);
    delete_test_file();

}
END_TEST

START_TEST (test_attr_float_wrapper)
{
    hid_t fid;
    char attr_name[] = FLOAT_ATTR_NAME;
    float data = FLOAT_ATTR_VALUE;
    float result;
    herr_t status;

    /* Create the file */
    create_test_h5_file();
    fid = open_test_file();

    /* Test the wrappers */
    status = danu_attr_write_float(fid,attr_name,data);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to write float attribute through wrapper");

    close_test_file(fid);
    fid = open_test_file();
    status = danu_attr_read_float(fid,attr_name,&result);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to read float attribute through wrapper");

    fail_unless( data == result,
                 "Data written through wrapper did not match wrapper read result");

    close_test_file(fid);
    delete_test_file();

}
END_TEST

START_TEST (test_attr_uint)
{
    hid_t fid;
    char attr_name[] = UINT_ATTR_NAME;
    unsigned int data = UINT_ATTR_VALUE;
    unsigned int result;
    herr_t status;

    /* Create the test file */
    create_test_h5_file();
    fid = open_test_file();

    /* Write the uinteger value */
    status = danu_attr_write(fid,attr_name,&data,H5T_NATIVE_UINT);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to write an uint attribute");
    
    /* Now read the value after closing the file....this forces HDF5 to flush memory */
    close_test_file(fid);
    fid = open_test_file();
    status = danu_attr_read(fid,attr_name,&result,H5T_NATIVE_UINT);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to read uint attribute");

    fail_unless( data == result,
                 "Data written and does not match read result");

    close_test_file(fid);
    delete_test_file();

}
END_TEST

START_TEST (test_attr_uint_wrapper)
{
    hid_t fid;
    char attr_name[] = UINT_ATTR_NAME;
    unsigned int data = UINT_ATTR_VALUE;
    unsigned int result;
    herr_t status;

    /* Create the file */
    create_test_h5_file();
    fid = open_test_file();

    /* Test the wrappers */
    status = danu_attr_write_uint(fid,attr_name,data);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to write uint attribute through wrapper");

    close_test_file(fid);
    fid = open_test_file();
    status = danu_attr_read_uint(fid,attr_name,&result);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to read uint attribute through wrapper");

    fail_unless( data == result,
                 "Data written through wrapper did not match wrapper read result");

    close_test_file(fid);
    delete_test_file();

}
END_TEST

START_TEST (test_attr_string)
{
    hid_t fid;
    char attr_name[] = STRING_ATTR_NAME;
    char *data;
    int data_size;
    char result[128];
    herr_t status;
    int iseed = 123456;

    /* Generate random string data */
    data_size =  generate_random_bound_int(64,128,&iseed);
    data = (char *)malloc(sizeof(char)*data_size);
    generate_random_string(data,data_size,&iseed);

    /* Create the file */
    create_test_h5_file();
    fid = open_test_file();

    /* Test the wrappers */
    status = danu_attr_write_string(fid,attr_name,data);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to write string attribute");
    
    close_test_file(fid);

    /* RE-open the file */
    fid = open_test_file();
    status = danu_attr_read_string(fid,attr_name,result,128);
    fail_unless( H5_RETURN_OK(status),
                 "Failed to read string attribute");

    printf("data=%s result=%s\n",data,result);

    fail_unless( strncmp(data,result,data_size) == 0,
                 "Data written through wrapper did not match read result");

    close_delete_test_file(fid);

}
END_TEST

START_TEST (test_attr_count_zero)
{
  hid_t fid;
  int num;
  herr_t status;

  create_test_h5_file();
  fid = open_test_file();

  status = danu_attr_count(fid,&num);
  fail_unless ( H5_RETURN_OK(status),
                "Call to danu_attr_count failed");

  fail_unless( num == 0,
               "Number of found attributes is not zero");

  close_test_file(fid);
  delete_test_file();
}
END_TEST

START_TEST (test_attr_count)
{
  hid_t fid;
  herr_t status;
  int count;
  int inum = INT_ATTR_VALUE;
  float fnum = FLOAT_ATTR_VALUE;
  double dnum = DOUBLE_ATTR_VALUE;

  fid = create_open_test_h5_file();

  H5LTset_attribute_int(fid,"/", INT_ATTR_NAME,&inum,1);
  H5LTset_attribute_float(fid,"/", FLOAT_ATTR_NAME,&fnum,1);
  H5LTset_attribute_double(fid,"/", DOUBLE_ATTR_NAME,&dnum,1);

  close_test_file(fid);

  fid = open_test_file();

  status = danu_attr_count(fid,&count);
  fail_unless(H5_RETURN_OK(status),
              "Failed to count attributes");

  fail_unless(count == 3,
              "danu_attr_count returned the incorrect number of attributes");

  close_delete_test_file(fid);
}
END_TEST

START_TEST (test_attr_names)
{
  hid_t fid;
  herr_t status;

  int inum = INT_ATTR_VALUE;
  float fnum = FLOAT_ATTR_VALUE;
  double dnum = DOUBLE_ATTR_VALUE;

  int i,num_found;
  char   **names;
  size_t len;

  fid = create_open_test_h5_file();
  H5LTset_attribute_int(fid,"/", INT_ATTR_NAME,&inum,1);
  H5LTset_attribute_float(fid,"/", FLOAT_ATTR_NAME,&fnum,1);
  H5LTset_attribute_double(fid,"/", DOUBLE_ATTR_NAME,&dnum,1);
  close_test_file(fid);

  names = (char**) malloc(3*sizeof(char *));

  fid = open_test_file();


  status = danu_attr_names(fid,3,names,&num_found);
  fail_unless(H5_RETURN_OK(status),
              "Failed to read attribute names");

  fail_unless(num_found == 3,
	      "Number of attribute names found incorrect");

  len = strlen(INT_ATTR_NAME);
  fail_unless ( strncmp(INT_ATTR_NAME,names[0],len) == 0,
                "Incorrect attribute name read");       
  
  len = strlen(FLOAT_ATTR_NAME);
  fail_unless ( strncmp(FLOAT_ATTR_NAME,names[1],len) == 0,
                "Incorrect attribute name read");       

  len = strlen(DOUBLE_ATTR_NAME);
  fail_unless ( strncmp(DOUBLE_ATTR_NAME,names[2],len) == 0,
                "Incorrect attribute name read"); 

  status = danu_attr_names(fid,2,names,&num_found);
  fail_unless(H5_RETURN_FAIL(status),
	      "Did not fail with array too small");

  fail_unless(num_found == 3,
	      "Did not return correct number of attributes found");

  for(i=0;i<3;i++) {
    free(names[i]);
  }
  free(names);

  names = (char **)malloc(4*sizeof(char*));
  status = danu_attr_names(fid,4,names,&num_found);
  fail_unless(H5_RETURN_OK(status),
	      "Failed to return correct status");

  fail_unless(strlen(names[3]) == 0,
	      "Last attribute entry was not NULL");

  for(i=0;i<4;i++) {
      free(names[i]);
  }
  free(names);



  close_delete_test_file(fid);

}
END_TEST


/* Define the test suite */
Suite *
attribute_suite (void)
{
    Suite *s = suite_create("Danu Attribute");

    TCase *attribute_error = tcase_create("Attribute Errors");
    tcase_add_test(attribute_error,test_attr_error_badid);
    tcase_add_test(attribute_error,test_attr_error_badtype);
    tcase_add_test(attribute_error,test_attr_error_baddata);
    tcase_add_test(attribute_error,test_attr_error_badname);
    tcase_add_test(attribute_error,test_attr_error_count);
    tcase_add_test(attribute_error,test_attr_error_names);
    suite_add_tcase(s,attribute_error);

    TCase *attribute_counter = tcase_create("Attribute Counter");
    tcase_add_test(attribute_counter,test_attr_count_zero);
    tcase_add_test(attribute_counter,test_attr_count);
    tcase_add_test(attribute_counter,test_attr_names);
    suite_add_tcase(s,attribute_counter);

    TCase *attribute = tcase_create("Attribute");
    tcase_add_test(attribute,test_attr_exists);
    tcase_add_test(attribute,test_attr_int);
    tcase_add_test(attribute,test_attr_int_wrapper);
    tcase_add_test(attribute,test_attr_double);
    tcase_add_test(attribute,test_attr_double_wrapper);
    tcase_add_test(attribute,test_attr_float);
    tcase_add_test(attribute,test_attr_float_wrapper);
    tcase_add_test(attribute,test_attr_uint);
    tcase_add_test(attribute,test_attr_uint_wrapper);
    tcase_add_test(attribute,test_attr_string);
    tcase_add_test(attribute,test_attr_overwrite);
    suite_add_tcase(s,attribute);

    return s;
}


      
int main(void)
{
  int number_failed;
  Suite *s = attribute_suite ();
  SRunner *sr = srunner_create (s);
  
  srunner_set_log(sr, "results.log");
  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}






