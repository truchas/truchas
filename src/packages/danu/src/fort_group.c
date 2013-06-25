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
*
*  DANU Fortran/C Attribute Interface
*
*/
#include<string.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_memory.h>
#include <danu_error.h>
#include <danu_group.h>
#include <danu_fort_hid.h>
#include <danu_fort_error.h>
#include <danu_fort_strings.h>
#include <danu_fort_group.h>

#ifndef MIN
#define MIN(a,b)   ( (a) < (b) ? (a) : (b) )
#endif

/* Group exists */
void danu_group_exists_f(hid_t_ptr *ptr, const char *name, int *flen, int *exists, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  hbool_t flag;
  int c_len = *flen + 1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(name,*flen,c_name,c_len);

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    flag = danu_group_exists(loc,c_name);
    if ( flag ) {
      *exists = (int) flag;
      err = DANU_SUCCESS;
    }
    else {
      *exists = FALSE;
    }
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert fortran to c string");
  }

  error_translate(err,ierr);

}

/* Basic group control */
void danu_group_create_f(hid_t_ptr *ptr, const char *name, int *flen, hid_t_ptr *group, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  hid_t gid;
  herr_t err;
  int c_len = *flen + 1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(name,*flen,c_name,c_len);

  err = DANU_FAILURE;
  *group = NULL;
  if ( status == DANU_SUCCESS ) {
    gid = danu_group_create(loc,c_name);
    if ( H5_ISA_VALID_ID(gid) ) {
      err = DANU_SUCCESS;
      *group = create_hid_struct(gid);
    }
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert Fortran string to C string");
  }

  error_translate(err,ierr);

}

void danu_group_open_f(hid_t_ptr *ptr, const char *name, int *flen, hid_t_ptr *group, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  hid_t gid;
  herr_t err;
  int c_len = *flen + 1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(name,*flen,c_name,c_len);

  err = DANU_FAILURE;
  *group = NULL;
  if ( status == DANU_SUCCESS ) {
    *group = create_hid_struct(0);  
    gid = danu_group_open(loc,c_name);
    if ( H5_ISA_VALID_ID(gid) ) {
      err = DANU_SUCCESS;
      *group = create_hid_struct(gid);
    }
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert Fortran string to C string");
  }

  error_translate(err,ierr);

}

void danu_group_close_f(hid_t_ptr *group, int *ierr)
{
  
  hid_t gid = GET_HID_VALUE(group);
  herr_t err = danu_group_close(gid);

  destroy_hid_struct(*group);

  error_translate(err,ierr);

} 

