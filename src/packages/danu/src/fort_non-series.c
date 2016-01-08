/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
*
*  DANU Fortran/C Non-series Data Interface
*
*/
#include<string.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_memory.h>
#include <danu_error.h>
#include <danu_non-series.h>
#include <danu_sim.h>
#include <danu_fort_error.h>
#include <danu_fort_transpose.h>
#include <danu_fort_strings.h>
#include <danu_fort_non-series.h>

/* Dataset open */
void data_open_dataset_f(const hid_t_ptr *fptr, const char *data_name, const int *len, hid_t_ptr *hid, int *ierr)
{
    hid_t sid = GET_HID_VALUE(fptr);
    hid_t data;
    herr_t err;
    int  c_len        = *len+1;
    char *c_data       = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(data_name,*len,c_data,c_len);
    
    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
        data = data_open_dataset(sid,c_data);
        DANU_FREE(c_data);
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }    

    if ( H5_ISA_VALID_ID(data) ) {
      *hid = create_hid_struct(data);
      err = DANU_SUCCESS;
    }
     
    error_translate(err,ierr);
    
}

/* Dataset query */

void data_exists_f(const hid_t_ptr *fptr, const char *data_name, const int *len, int *flag, int *ierr)
{
    hid_t sid = GET_HID_VALUE(fptr);
    herr_t err;
    int  c_len        = *len+1;
    char *c_data       = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(data_name,*len,c_data,c_len);
    
    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
        err = data_exists(sid,c_data,flag);
        DANU_FREE(c_data);
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }    
     
    error_translate(err,ierr);
    
}

void data_count_f(const hid_t_ptr *fptr, int *count, int *ierr)
{
    hid_t sid = GET_HID_VALUE(fptr);
    herr_t err;

    *count = 0;
    err = data_count(sid,count);

    error_translate(err,ierr);

}

void data_list_f(const hid_t_ptr *fptr, char *names, const int *len, const int *num, int *ierr)
{

    hid_t sid = GET_HID_VALUE(fptr);
    herr_t err, flag;
    danu_err_t copy_ok;
    int i, num_found;
    char **data_names;
    int min_free;

    data_names = DANU_MALLOC(char *, *num);
    
    err = data_list(sid,*num,data_names,&num_found);

    if ( H5_RETURN_OK(err) && num_found > 0) {
      copy_ok =
       copy_c_array_to_fortran_array((const char * const *)data_names,num_found,names,num,len);

      if ( DANU_RETURN_OK(copy_ok) && ( *num >= num_found ) ) {
	err = DANU_SUCCESS;
      }
      else {

	if ( ! DANU_RETURN_OK(copy_ok) ) {
	  DANU_ERROR_MESS("Failed to copy C strings to fortran array");
	}

	if ( *num < num_found ) {
	  DANU_ERROR_MESS("data_names array too small to hold all data"
			  " names");
	}

	err = DANU_FAILURE;
      }
   
      min_free = num_found < *num ? num_found : *num;
      for(i=0;i<min_free;i++)
        DANU_FREE(data_names[i]);

    }
    else {

      if ( H5_RETURN_FAIL(err) ) {
        DANU_ERROR_MESS("Failed to read data name list");
      }

      if ( num_found == 0 ) {
	DANU_WARN_MESS("No non-series data found");
      }

      err = DANU_FAILURE;
    }


    DANU_FREE(data_names);


    error_translate(err,ierr);
        
}

void data_type_f(const hid_t_ptr *sptr,
		 const char *data_name, const int *flen, int *type, int *ierr)
{
  hid_t sid = GET_HID_VALUE(sptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);

  if ( status == DANU_SUCCESS ) {
    err = data_type(sid,c_data,type);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void data_rank_f(const hid_t_ptr *sptr,
		 const char *data_name, const int *flen, int *rank, int *ierr)
{
  hid_t sid = GET_HID_VALUE(sptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);

  if ( status == DANU_SUCCESS ) {
    err = data_rank(sid,c_data,rank);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void data_dimensions_f(const hid_t_ptr *sptr,
		       const char *data_name, const int *flen,
		       const int *rank, int *dimensions, int *ierr)
{
  hid_t sid = GET_HID_VALUE(sptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dims;

  if ( status == DANU_SUCCESS ) {
    c_dims=DANU_MALLOC(int,*rank);
    err = data_dimensions(sid,c_data,*rank,c_dims);
    transpose_array_int(*rank,c_dims,dimensions);
    DANU_FREE(c_data);
    DANU_FREE(c_dims);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void data_write_byte0_f(const hid_t_ptr *fptr, const char *data_name, const int *flen,
		       int8_t *data, int *ierr)
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = data_write_byte(sid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void data_write_byte_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, int8_t *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim, i, j;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = data_write_byte(sid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void data_write_int0_f(const hid_t_ptr *fptr, const char *data_name, const int *flen,
		       int *data, int *ierr)
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = data_write_int(sid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void data_write_int_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, int *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim, i, j;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = data_write_int(sid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void data_write_float0_f(const hid_t_ptr *fptr, const char *data_name, const int *flen,
		         float *data, int *ierr)
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = data_write_float(sid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void data_write_float_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, float *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim, i, j;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = data_write_float(sid,c_data,*dim,c_dim,data);
    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void data_write_double0_f(const hid_t_ptr *fptr, const char *data_name, const int *flen,
		           double *data, int *ierr)
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = data_write_double(sid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void data_write_double_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, double *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = data_write_double(sid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void data_read_byte0_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      int8_t *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim = 1;
    size = 1;
    err = data_read_byte(sid,c_data,dim,&size,data);

    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}


void data_read_byte_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, int8_t *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = data_read_byte(sid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void data_read_int0_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      int *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim = 1;
    size = 1;
    err = data_read_int(sid,c_data,dim,&size,data);

    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}


void data_read_int_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, int *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = data_read_int(sid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void data_read_float0_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		         float *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim = 1;
    size = 1;
    err = data_read_float(sid,c_data,dim,&size,data);

    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}


void data_read_float_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, float *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = data_read_float(sid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void data_read_double0_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		         double *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim = 1;
    size = 1;
    err = data_read_double(sid,c_data,dim,&size,data);

    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

 
void data_read_double_f(const hid_t_ptr *fptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, double *data, int *ierr )
{
  hid_t sid = GET_HID_VALUE(fptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = data_read_double(sid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}
 



