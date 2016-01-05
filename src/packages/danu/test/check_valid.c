/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
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

#include <danu_types.h>
#include <danu_h5_error.h>
#include <danu_error.h>
#include <danu_utils.h>
#include <danu_file.h>
#include <danu_h5_error.h>

#define FILENAME "deleteme-valid.h5"


/* BEGIN TESTS */

START_TEST (test_file_valid)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;
    htri_t res;

    /* Create a file */
    fid = danu_file_create(test_file);
    res = H5Iis_valid(fid);
    fail_unless(res == TRUE,
	        "Simple create failed to return valid HID");
    fail_unless(H5_ISA_VALID_ID(fid),
	        "Result failed to match VALID macro");
    fail_unless(H5_ISA_FILE_ID(fid),
	        "File id macro failed");
    danu_file_close(fid);

    /* Check ID after a close */
    res = H5Iis_valid(fid);
    fail_unless(res != TRUE,
	        "Failed to return !TRUE on valid query");
    fail_unless(H5_ISA_INVALID_ID(fid),
	        "INVALID macro failed");
    fail_unless(!H5_ISA_FILE_ID(fid),
	        "File id macro failed");

    /* Open Read Only */
    fid = danu_file_open_rdonly(test_file);
    res = H5Iis_valid(fid);
    fail_unless(res == TRUE,
	        "Simple read only open failed to return a valid HID");
    fail_unless(H5_ISA_VALID_ID(fid),
	        "Result failed to match macro");
    fail_unless(H5_ISA_FILE_ID(fid),
	        "File id macro failed");

    danu_file_close(fid);

    danu_file_delete(test_file);
}
END_TEST

Suite *
file_suite (void)
{
    Suite *s = suite_create("Danu Valid Id");

    TCase *file_utils = tcase_create("File ID");
    tcase_add_test(file_utils,test_file_valid);
    suite_add_tcase(s,file_utils);

    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = file_suite ();
  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}



