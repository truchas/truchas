/* *************************************************************************** *
*                                                                              *
*                                                                              *
*                             Copyright  (C) 20xx,                             *
*                      Los Alamos National Security, LLC                       *
*                                                                              *
*                             LA-CC-xxxxxx                                     *
*                                                                              *
* **************************************************************************** */

/*
* unit_test_utils.h
*
*  DANU unit test utility functions
*
*/

#ifndef UTEST_UTILS_H
#define UTEST_UTILS_H

#include <hdf5.h>

#define TEST_PASS 0
#define TEST_FAIL 1 

#define TEST_FILENAME "testfile.h5"
#define TEST_GROUP    "Test Group"

#define TEST_HID 0x1


/* Test Utilities */
void   create_test_h5_file(void);
hid_t  create_open_test_h5_file(void);
hid_t  open_test_file(void);
void   close_test_file(hid_t id);
void   delete_test_file(void);
void   close_delete_test_file(hid_t id);

hid_t  test_group_create(hid_t loc, const char *name);
hid_t  test_group_open(hid_t loc, const char *name);
void   test_group_close(hid_t gid);

void   print_fortran_string_in_c(char *fort_string, int len);

int char_array_cmp(char **a, char **b, int size); 
int int_array_cmp(const int *a, const int *b, int num);
int double_array_cmp(const double *a, const double *b, int num);
int float_array_cmp(const float *a, const float *b, int num);

void   exit_now(int);
void   fail_exit_now(void);
void   pass_exit_now(void);

#endif

