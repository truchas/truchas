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

#include <danu.h>

#include <danu_output.h>

#define FILENAME "output-deleteme.h5"

/* Test Utils */
void create_dummy_file();
void create_dummy_h5_file();
void create_dummy_output_file();
void delete_dummy_output_file();

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
    hid_t id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    H5Fclose(id);
}


void create_dummy_output_file()
{
    const char name[] = FILENAME;
    hid_t fid; 
    output_file_create(name,&fid);
    output_file_close(&fid);

}

void delete_dummy_output_file()
{
    const char name[] = FILENAME;
    int stat = remove(name);
}

/* Begin Tests */
START_TEST (test_ofile_valid)
{
    const char test_file[] = FILENAME;
    hbool_t flag;
    hid_t id;

    create_dummy_h5_file();
    id  = H5Fopen(test_file,H5F_ACC_RDONLY,H5P_DEFAULT);
    flag = output_file_is_valid(id);
    fail_unless( flag == FALSE,
                 "Failed to flag invalid output file");
    H5Fclose(id);
    delete_dummy_output_file();

    create_dummy_output_file();
    id  = H5Fopen(test_file,H5F_ACC_RDONLY,H5P_DEFAULT);
    flag = output_file_is_valid(id);
    fail_unless ( flag != FALSE,
                  "Failed to flag valid output file");
    H5Fclose(id);
    delete_dummy_output_file();
}
END_TEST

START_TEST (test_ofile_flags)
{
    const char test_file[] = FILENAME;
    hid_t id;

    /* Test all possible open combinations */
    create_dummy_output_file();

    id = output_file_open(test_file,DANU_FILE_ACC_RDONLY,DANU_FILE_ACT_OPEN);
    fail_unless (H5_ISA_VALID_ID(id),
                 "Failed to open file with RDONLY, OPEN flags");
    H5Fclose(id);

    id = output_file_open(test_file,DANU_FILE_ACC_RDWR,DANU_FILE_ACT_OPEN);
    fail_unless (H5_ISA_VALID_ID(id),
                 "Failed to open file with RDWR, OPEN flags");
    H5Fclose(id);

    id = output_file_open(test_file,DANU_FILE_ACC_APPEND,DANU_FILE_ACT_OPEN);
    fail_unless (H5_ISA_VALID_ID(id),
                 "Failed to open file with APPEND, OPEN flags");
    H5Fclose(id);
    delete_dummy_output_file();

    /* Test all possible create combinations */
    id = output_file_open(test_file,DANU_FILE_ACC_RDONLY,DANU_FILE_ACT_CREATE);
    fail_unless (H5_ISA_VALID_ID(id),
                 "Failed to create file with RDONLY, CREATE flags");
    H5Fclose(id);

    id = output_file_open(test_file,DANU_FILE_ACC_RDWR,DANU_FILE_ACT_CREATE);
    fail_unless (H5_ISA_VALID_ID(id),
                 "Failed to create file with RDONLY, CREATE flags");
    H5Fclose(id);

    id = output_file_open(test_file,DANU_FILE_ACC_APPEND,DANU_FILE_ACT_CREATE);
    fail_unless (H5_ISA_VALID_ID(id),
                 "Failed to create file with APPEND, CREATE flags");
    H5Fclose(id);
    delete_dummy_output_file();

}
END_TEST
		
START_TEST (test_ofile_flag_errors)
{
    const char test_file[] = FILENAME;
    hid_t id;
    unsigned int bad_flag = 0x111;

    /* Test all possible access combinations */
    id = output_file_open(test_file,DANU_FILE_ACC_RDONLY,bad_flag);
    fail_unless (H5_ISA_INVALID_ID(id),
                 "Failed to catch bad ACTION FLAG and RDONLY");

    id = output_file_open(test_file,DANU_FILE_ACC_RDWR,bad_flag);
    fail_unless (H5_ISA_INVALID_ID(id),
                 "Failed to catch bad ACTION FLAG and RDWR");

    id = output_file_open(test_file,DANU_FILE_ACC_APPEND,bad_flag);
    fail_unless (H5_ISA_INVALID_ID(id),
                 "Failed to catch bad ACTION FLAG and APPEND");

    /* Test all possible create/open combinations */
    id = output_file_open(test_file,bad_flag,DANU_FILE_ACT_CREATE);
    fail_unless (H5_ISA_INVALID_ID(id),
                 "Failed to catch bad ACCESS FLAG and CREATE");

    id = output_file_open(test_file,bad_flag,DANU_FILE_ACT_OPEN);
    fail_unless (H5_ISA_INVALID_ID(id),
                 "Failed to catch bad ACCESS FLAG and OPEN");

}
END_TEST


START_TEST (test_ofile_open_errors)
{
    const char test_file[] = FILENAME;
    hid_t id;

    id = output_file_open(test_file,DANU_FILE_ACC_RDONLY,DANU_FILE_ACT_OPEN);
    fail_unless ( H5_ISA_INVALID_ID(id),
                  "Failed to flag an open of file that does not exist");

    create_dummy_file();
    id = output_file_open(test_file, DANU_FILE_ACC_RDONLY, DANU_FILE_ACT_OPEN);
    fail_unless ( H5_ISA_INVALID_ID(id),
                  "Failed to flag an open of invalid output file");
    delete_dummy_output_file();

    create_dummy_h5_file();
    id = output_file_open(test_file, DANU_FILE_ACC_RDONLY, DANU_FILE_ACT_OPEN);
    fail_unless ( H5_ISA_INVALID_ID(id),
                   "Failed to return invalid id for an HDF5 file as invalid");
    delete_dummy_output_file();

}
END_TEST

START_TEST (test_ofile_create)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    status = output_file_create(test_file,&fid);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to create output file. Invalid H5 ID");

    status = output_file_close(&fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");
}
END_TEST

START_TEST (test_ofile_rdonly)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    create_dummy_output_file();

    status = output_file_open_rdonly(test_file,&fid);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to open output file rdonly mode. Invalid H5 ID");

    status = output_file_close(&fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");

    delete_dummy_output_file();
}
END_TEST

START_TEST (test_ofile_append)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    create_dummy_output_file();

    status = output_file_open_append(test_file,&fid);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to open output file in append mode. Invalid H5 ID");

    status = output_file_close(&fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");

    delete_dummy_output_file();
}
END_TEST

START_TEST (test_ofile_rdwr)
{
    const char test_file[] = FILENAME;
    hid_t  fid;
    herr_t status;

    create_dummy_output_file();

    status = output_file_open_rdwr(test_file,&fid);

    fail_unless ( H5_ISA_VALID_ID(fid),
                  "Failed to open output file in rdwr mode. Invalid H5 ID");

    status = output_file_close(&fid);

    fail_unless( DANU_RETURN_OK(status),
                 "Failed to close file");

    delete_dummy_output_file();
}
END_TEST


Suite *
ofile_suite (void)
{
    Suite *s = suite_create("Output File");

    TCase *ofile_errors = tcase_create("Output File Errors");
    tcase_add_test(ofile_errors, test_ofile_valid);
    tcase_add_test(ofile_errors, test_ofile_open_errors);
    tcase_add_test(ofile_errors, test_ofile_flag_errors);
    suite_add_tcase(s,ofile_errors);

    TCase *ofile_flags = tcase_create("Ouptut File Flags");
    tcase_add_test(ofile_flags, test_ofile_flags);
    suite_add_tcase(s,ofile_flags);

    TCase *ofile_create = tcase_create("File Create");
    tcase_add_test(ofile_create, test_ofile_create);
    suite_add_tcase(s,ofile_create);

    TCase *ofile_open = tcase_create("File Open");
    tcase_add_test(ofile_open, test_ofile_rdonly);
    tcase_add_test(ofile_open, test_ofile_append);
    tcase_add_test(ofile_open, test_ofile_rdwr);
    suite_add_tcase(s,ofile_open);

    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = ofile_suite ();
  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}



