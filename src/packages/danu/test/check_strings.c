/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
/* ****************************************************************************
*                                                                             *
*  Danu Unit Test                                                             *
*  Unit: Fortran Strings                                                      *
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

#include <danu_memory.h>
#include <danu_types.h>
#include <danu_h5_error.h>
#include <danu_error.h>

#include <danu_fort_strings.h>

#define MAX_STRING_LEN 128
#define MAX_STRING_NUM 8

#define RANDOM_SEED 0x0FF456789


/* BEGIN TESTS */


START_TEST(test_trim_len)
{
  char *fort_string;
  int  len_bound = MAX_STRING_LEN/2;
  int  len, trim_len;
  int seed = RANDOM_SEED;

  len = generate_random_bound_int(1,len_bound,&seed);
  fort_string = DANU_MALLOC(char,MAX_STRING_LEN);
  memset(fort_string,' ',MAX_STRING_LEN);
  generate_random_fortran_string(fort_string,len,&seed);

  trim_len = compute_fortran_trim_len(fort_string,MAX_STRING_LEN);
 
  fail_unless( len == trim_len,
	       "Failed to compute the correct trim length");

   DANU_FREE(fort_string);
}
END_TEST

START_TEST( test_convert_f2c)
{
    danu_err_t status;
    char *fort_string, *c_string;
    int len_bound; 
    int len, trim_len, c_len;
    int seed = RANDOM_SEED;

    len_bound = generate_random_bound_int(1,MAX_STRING_LEN/2,&seed);
    len = generate_random_bound_int(len_bound,MAX_STRING_LEN,&seed);

    fort_string = DANU_MALLOC(char, len);
    memset(fort_string,' ',len);
    generate_random_fortran_string(fort_string,len,&seed);

    trim_len = compute_fortran_trim_len(fort_string,len);

    /* Tests that should fail */
    c_len = trim_len - 1;
    c_string = DANU_MALLOC(char, c_len);
    status = convert_string_f2c(fort_string,len,c_string, c_len);
    fail_unless(status == DANU_FAILURE,
	        "Failed to flag too short c_string");
    DANU_FREE(c_string);

    /* This should pass */
    c_len = len + 5;
    c_string = DANU_MALLOC(char,c_len);
    status = convert_string_f2c(fort_string,len,c_string, c_len);
    fail_unless(status == DANU_SUCCESS,
	        "Failed to copy fortran string to c string");

    /* And this verifies that each is the same up to the trim length */
    fail_unless(strncmp(fort_string,c_string,trim_len) == 0,
	        "Failed to copy fortran string correctly");

    DANU_FREE(c_string);
    DANU_FREE(fort_string);

}
END_TEST

START_TEST( test_convert_c2f)
{
    danu_err_t status;
    char *fort_string, *c_string;
    int len_bound; 
    int len, trim_len, c_len;
    int seed = RANDOM_SEED;

    len_bound = generate_random_bound_int(1,MAX_STRING_LEN/2,&seed);
    c_len = generate_random_bound_int(len_bound,MAX_STRING_LEN,&seed);

    c_string = DANU_MALLOC(char, c_len);
    generate_random_string(c_string,c_len,&seed);

    /* Tests that should fail */
    len = c_len - 5;
    fort_string = DANU_MALLOC(char, len);
    status = convert_string_c2f(c_string, fort_string, len);
    fail_unless(status == DANU_FAILURE,
	        "Failed to flag too short fortran string");
    DANU_FREE(fort_string);

    /* This should pass */
    len = c_len + 5;
    fort_string = DANU_MALLOC(char,len);
    status = convert_string_c2f(c_string, fort_string, len);
    fail_unless(status == DANU_SUCCESS,
	        "Failed to copy C string to fortran string");

    /* And this verifies that each is the same up to the trim length */
    fail_unless(strncmp(fort_string,c_string,c_len-1) == 0,
	        "Failed to copy C string correctly");

    trim_len = compute_fortran_trim_len(fort_string,len);
    printf("trim_len=%d, strlen=%d\n", trim_len, strlen(c_string));
    print_fortran_string_in_c(fort_string,len);
    printf("c_string='%s'\n",c_string);
    fail_unless(trim_len == strlen(c_string),
		"Failed to add appropriate white space to fortran string");

    DANU_FREE(c_string);
    DANU_FREE(fort_string);

}
END_TEST

Suite *
string_suite (void)
{
    Suite *s = suite_create("Danu Fortran Strings");

    TCase *string_utils = tcase_create("Fortran String Utils");
    tcase_add_test(string_utils,test_trim_len);
    tcase_add_test(string_utils,test_convert_f2c);
    tcase_add_test(string_utils,test_convert_c2f);
    suite_add_tcase(s,string_utils);

    return s;
}

int main(void)
{
  int number_failed = 0;
  
  Suite *s = string_suite ();
  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
