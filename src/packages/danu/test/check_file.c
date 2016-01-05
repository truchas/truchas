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
#include <danu_memory.h>
#include <danu_h5_error.h>
#include <danu_error.h>
#include <danu_utils.h>
#include <danu_file.h>

#define FILENAME "deleteme-file.h5"


/* Test Utilities */
void create_dummy_file(void);
void create_dummy_h5_file(void);
void delete_dummy_file(void);

void create_dummy_file()
{
    const char name[] = FILENAME;
    char buffer[] = {'A', 'B', 'C'};
    FILE *fh;
    
    
    fh = fopen(name,"wb");
    fwrite(buffer,1,sizeof(buffer),fh);
    fclose(fh);

}

void create_dummy_h5_file()
{
    const char name[] = FILENAME;

    hid_t id = danu_file_create(name);
    herr_t status = danu_file_close(id);

}

void delete_dummy_file()
{
    const char name[] = FILENAME;
    int stat = remove(name);
}


/* BEGIN TESTS */

START_TEST (test_file_create)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    fid = danu_file_create(test_file);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to create file. Invalid H5 ID");

    status = danu_file_close(fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");

    delete_dummy_file();
}
END_TEST

START_TEST (test_file_rdonly)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    create_dummy_h5_file();

    fid = danu_file_open_rdonly(test_file);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to open file rdonly mode. Invalid H5 ID");

    status = danu_file_close(fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");

    delete_dummy_file();
}
END_TEST

START_TEST (test_file_append)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    create_dummy_h5_file();

    fid = danu_file_open_append(test_file);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to open file in append mode. Invalid H5 ID");

    status = danu_file_close(fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");

    delete_dummy_file();
}
END_TEST

START_TEST (test_file_rdwr)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    create_dummy_h5_file();

    fid = danu_file_open_rdwr(test_file);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to open file in rdwr mode. Invalid H5 ID");

    status = danu_file_close(fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");

    delete_dummy_file();
}
END_TEST

START_TEST (test_file_flush)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    create_dummy_h5_file();

    fid = danu_file_open_rdwr(test_file);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to open file in rdwr mode. Invalid H5 ID");

    status = danu_file_flush(fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to flush file");
 
    status = danu_file_close(fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");

    delete_dummy_file();

}
END_TEST

START_TEST (test_file_flush_global)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    create_dummy_h5_file();

    fid = danu_file_open_rdwr(test_file);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to open file in rdwr mode. Invalid H5 ID");

    status = danu_file_flush_global(fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed a global flush");
 
    status = danu_file_close(fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");

    delete_dummy_file();

}
END_TEST

START_TEST (test_file_create_errors)
{
    char *test_file;
    hid_t  fid;
    herr_t status;

    fid = danu_file_create(NULL);

    fail_unless (  H5_ISA_INVALID_ID(fid),
                  "NULL pointer failed to return bad id");

    test_file = DANU_MALLOC(char,8);
    fid = danu_file_create(test_file);

    fail_unless ( H5_ISA_INVALID_ID(fid),
                  "Uninitialized pointer failed to return bad id");

    DANU_FREE(test_file);

}
END_TEST

START_TEST (test_file_open_errors)
{
    char *test_file;
    hid_t  fid;
    herr_t status;
    FILE *fh;

    fid = danu_file_open_rdonly(NULL);
    fail_unless ( H5_ISA_INVALID_ID(fid),
                  "NULL pointer failed to return bad id");

    test_file = DANU_MALLOC(char,8);
    fid = danu_file_open_rdonly(test_file);
    fail_unless ( H5_ISA_INVALID_ID(fid),
                  "Uninitialized pointer failed to return bad id");
    DANU_FREE(test_file);

    fid = danu_file_open_rdonly(FILENAME);
    fail_unless( H5_ISA_INVALID_ID(fid) ,
                 "Open of file that does not exist failed to return a bad id");

    create_dummy_file();
    fid = danu_file_open_rdonly(FILENAME);
    fail_unless(  H5_ISA_INVALID_ID(fid) ,
                 "Open of file that is not an HDF5 file failed to return a bad id");

    delete_dummy_file();

}
END_TEST

START_TEST (test_file_exists)
{
    const char test_file[] = FILENAME;

    create_dummy_h5_file();
    fail_unless( FILE_EXISTS(test_file),
               "Test file does not exist");
    delete_dummy_file();


}
END_TEST

START_TEST (test_file_delete)
{
    const char test_file[] = FILENAME;
    danu_err_t err;

    create_dummy_h5_file();
    err = danu_file_delete(test_file);
    fail_unless ( DANU_RETURN_SUCCESS(err),
	              "Failed to delete existing file");

    if ( ! DANU_RETURN_SUCCESS(err) ) {
	delete_dummy_file();
    }

}
END_TEST

Suite *
file_suite (void)
{
    Suite *s = suite_create("Danu File");

    TCase *file_utils = tcase_create("File Utils");
    tcase_add_test(file_utils,test_file_exists);
    tcase_add_test(file_utils,test_file_delete);
    suite_add_tcase(s,file_utils);

    TCase *file_errors = tcase_create("File Errors");
    tcase_add_test(file_errors, test_file_create_errors);
    tcase_add_test(file_errors, test_file_open_errors);
    suite_add_tcase(s,file_errors);

    TCase *file_create = tcase_create("File Create");
    tcase_add_test(file_create, test_file_create);
    suite_add_tcase(s,file_create);

    TCase *file_open = tcase_create("File Open");
    tcase_add_test(file_open, test_file_rdonly);
    tcase_add_test(file_open, test_file_append);
    tcase_add_test(file_open, test_file_rdwr);
    suite_add_tcase(s,file_open);

    TCase *file_flush = tcase_create("File Flush");
    tcase_add_test(file_flush, test_file_flush);
    tcase_add_test(file_flush, test_file_flush_global);
    suite_add_tcase(s,file_flush);

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



