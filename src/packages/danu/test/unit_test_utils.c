/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "hdf5.h"

#include <danu_memory.h>

#include "unit_test_utils.h"



void create_test_h5_file()
{
    hid_t fid = H5Fcreate(TEST_FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 

    H5Fclose(fid);
}

hid_t open_test_file()
{
    hid_t fid = H5Fopen(TEST_FILENAME, H5F_ACC_RDWR, H5P_DEFAULT); 

    return fid;
}

hid_t create_open_test_h5_file()
{
  hid_t fid = H5Fcreate(TEST_FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
  return fid;
}

void close_test_file(hid_t id)
{
    H5Fclose(id);
}

void delete_test_file()
{
    int stat = remove(TEST_FILENAME);
}
 
void close_delete_test_file(hid_t fid)
{
  int stat;
  close_test_file(fid);
  stat = remove(TEST_FILENAME);
}

hid_t test_group_create(hid_t loc, const char *name)
{
	hid_t gid;
	
	gid = H5Gcreate(loc,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	
	return gid;
}

hid_t test_group_open(hid_t loc, const char *name)
{
	hid_t gid;
	
	gid = H5Gopen(loc,name,H5P_DEFAULT);
	
	return gid;
}

void test_group_close(hid_t gid)
{
	H5Gclose(gid);
}

void print_fortran_string_in_c(char * fort_string, int len)
{
  char *p = DANU_MALLOC(char,len+1);
  memcpy(p,fort_string,len);
  p[len] = '\0';
  printf("FORTRAN_STRING=\"%s\"\n",p);
  DANU_FREE(p);
}

int char_array_cmp(char ** a, char **b, int size)
{
  int i;
  int ret;

  ret = 0;
  for(i=0;i<size;i++)
    ret+= strcmp(a[i],b[i]);

  return ret;

}

int int_array_cmp(const int *a, const int *b, int num)
{
    size_t bytes = sizeof(int)*num;
    return memcmp(a,b,bytes);
}

int double_array_cmp(const double *a, const double *b, int num)
{
    size_t bytes = sizeof(double)*num;
    return memcmp(a,b,bytes);
}

int float_array_cmp(const float *a, const float *b, int num)
{
    size_t bytes = sizeof(float)*num;
    return memcmp(a,b,bytes);
}
void exit_now(int code)
{
    printf("Exit Now %d\n",code);
    exit(code);
}

void fail_exit_now(void)
{
    exit_now(TEST_FAIL);
}

void pass_exit_now(void)
{
    exit_now(TEST_PASS);
}





