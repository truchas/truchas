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
#include <danu_attribute.h>
#include <danu_fort_hid.h>
#include <danu_fort_error.h>
#include <danu_fort_strings.h>
#include <danu_fort_attribute.h>

#ifndef MIN
#define MIN(a,b)   ( (a) < (b) ? (a) : (b) )
#endif

/* Attribute exists */
void danu_attr_exists_f(hid_t_ptr *ptr, const char *name, int *flen, int *exists, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  int c_len = *flen + 1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(name,*flen,c_name,c_len);

  if ( status == DANU_SUCCESS ) {
    err = danu_attr_exists(loc,c_name,exists);
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert fortran to c string");
  }

  error_translate(err,ierr);

}

/* Attribute count */
void danu_attr_count_f(hid_t_ptr *ptr, int *num_found, int *ierr)
{

  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err = danu_attr_count(loc,num_found);

  error_translate(err,ierr);

}

/* Attribute list */
void danu_attr_names_f(hid_t_ptr *ptr, char *names, int *flen, int *num, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  char **attr_names, *p;
  size_t str_len;
  int  i, num_found, cp_len;

  attr_names = DANU_MALLOC(char *, *num);
  err = danu_attr_names(loc, *num, attr_names, &num_found);

  p = names;
  for(i=0; i< *num; i++) {
      str_len = strlen(attr_names[i]);
      cp_len = MIN(str_len,*flen);
      memcpy(p,attr_names[i],cp_len);
      p+=*flen;
      DANU_FREE(attr_names[i]);
  }
  DANU_FREE(attr_names);

  error_translate(err,ierr);

}

/* Attribute write*/
void danu_attr_write_int_f(hid_t_ptr *ptr, const char *name, int *name_flen, int *value, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  int  c_len  = *name_flen+1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status  = convert_string_f2c(name,*name_flen,c_name,c_len);

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    err = danu_attr_write_int(loc, c_name, *value);
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert Fortran string to C string");
  }

  error_translate(err,ierr);

}

void danu_attr_write_real4_f(hid_t_ptr *ptr, const char *name, int *name_flen, float *value, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  int  c_len  = *name_flen+1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status  = convert_string_f2c(name,*name_flen,c_name,c_len);

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    err = danu_attr_write_float(loc, c_name, *value);
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert Fortran string to C string");
  }

  error_translate(err,ierr);
}

void danu_attr_write_real8_f(hid_t_ptr *ptr, const char *name, int *name_flen, double *value, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  int  c_len  = *name_flen+1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status  = convert_string_f2c(name,*name_flen,c_name,c_len);

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    err = danu_attr_write_double(loc, c_name, *value);
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert Fortran string to C string");
  }

  error_translate(err,ierr);
}

void danu_attr_write_char_f(hid_t_ptr *ptr,
			    const char *name,
			    int *name_flen,
			    const char *string,
			    int *str_flen,
			    int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  int  c_len  = *name_flen+1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status  = convert_string_f2c(name,*name_flen,c_name,c_len);
  int c_string_len = *str_flen + 1;
  char *c_string = DANU_MALLOC(char, c_string_len);


  err = DANU_FAILURE;
  if ( (status == DANU_SUCCESS) && (c_string != NULL) ) {
    status = convert_string_f2c(string,*str_flen,c_string, c_string_len);
    if ( status == DANU_SUCCESS ) {
      err = danu_attr_write_string(loc,c_name, c_string);
      DANU_FREE(c_string);
    }
    else { 
      DANU_ERROR_MESS("Failed to convert fortran character data to C string");
    }
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert Fortran string to C string");
  }

  error_translate(err,ierr);
}

void danu_attr_read_int_f(hid_t_ptr *ptr, const char *name, int *name_flen, int *buffer, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  int  c_len  = *name_flen+1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status  = convert_string_f2c(name,*name_flen,c_name,c_len);

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    err = danu_attr_read_int(loc, c_name, buffer);
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert Fortran string to C string");
  }

  error_translate(err,ierr);
}
void danu_attr_read_real4_f(hid_t_ptr *ptr, const char *name, int *name_flen, float *buffer, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  int  c_len  = *name_flen+1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status  = convert_string_f2c(name,*name_flen,c_name,c_len);

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    err = danu_attr_read_float(loc, c_name, buffer);
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert Fortran string to C string");
  }

  error_translate(err,ierr);
}

void danu_attr_read_real8_f(hid_t_ptr *ptr, const char *name, int *name_flen, double *buffer, int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  int  c_len  = *name_flen+1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status  = convert_string_f2c(name,*name_flen,c_name,c_len);

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    err = danu_attr_read_double(loc, c_name, buffer);
    DANU_FREE(c_name);
  }
  else {
    DANU_ERROR_MESS("Failed to convert Fortran string to C string");
  }

  error_translate(err,ierr);
}

void danu_attr_read_char_f(hid_t_ptr *ptr,
			    const char *name,
			    int *name_flen,
			    char *string,
			    int *str_flen,
			    int *ierr)
{
  hid_t loc = GET_HID_VALUE(ptr);
  herr_t err;
  int  c_len  = *name_flen+1;
  char *c_name = DANU_MALLOC(char, c_len);
  danu_err_t status  = convert_string_f2c(name,*name_flen,c_name,c_len);
  size_t c_string_len = *str_flen + 1;
  char *c_string = DANU_MALLOC(char, c_string_len);


  err = DANU_FAILURE;
  if ( (status == DANU_SUCCESS) && (c_string != NULL) ) {
    err = danu_attr_read_string(loc, c_name, c_string, c_string_len);
    status = convert_string_c2f(c_string, string, *str_flen);
    if ( status != DANU_SUCCESS ) {
      err = DANU_FAILURE;
    }
    DANU_FREE(c_name);
    DANU_FREE(c_string);
  }
  else {
    err = DANU_FAILURE;
    if ( status != DANU_SUCCESS ) {
     DANU_ERROR_MESS("Failed to convert attribute name to C string");
    }
    else {
     DANU_ERROR_MESS("Failed to create C buffer for read");
    }
  } 

  error_translate(err,ierr);
}

