/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
*
*  DANU Fortran/C Probe Data Interface
*
*/

#include <hdf5.h>

#include <danu_memory.h>
#include <danu_error.h>
#include <danu_probes.h>

#include <danu_fort_error.h>
#include <danu_fort_strings.h>
#include <danu_fort_probes.h>

/* Basic utilities */
void probe_data_exists_f(hid_t_ptr *sptr, 
			 const char *probe_name,
			 const int  *flen,
			 int *flag,
			 int *ierr)
{
    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;
    int  c_len        = *flen+1;
    char *c_data       = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(probe_name,*flen,c_data,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
      err = probe_exists(sid,c_data,flag);
      DANU_FREE(c_data);
    }

    error_translate(err,ierr);

}

void probe_data_dimensions_f(const hid_t_ptr *sptr, 
                             const char *probename, const int *flen,
                             int *len, int *num, int *ierr)
{
    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;
    int  c_len        = *flen+1;
    char *c_data      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(probename,*flen,c_data,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
      err =probe_data_length(sid,c_data,len);
      err &= probe_data_num(sid,c_data,num);
      DANU_FREE(c_data);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C");
    }

    error_translate(err,ierr);

}

void probe_data_dimensions2_f(const hid_t_ptr *pptr, 
                              int *len, int *num, int *ierr)
{
    hid_t pid = GET_HID_VALUE(pptr);
    herr_t err;

    err = DANU_FAILURE;
    err =probe_data_length2(pid,len);
    err &= probe_data_num2(pid,num);

    error_translate(err,ierr);

}


void probe_data_list_f(hid_t_ptr *sptr,
		       char *names,
		       const int *len,
		       const int *num,
		       int *ierr)
{
   
    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;
    danu_err_t copy_ok;
    int i, num_found, num_free;
    char **probe_names;

    probe_names = DANU_MALLOC(char *, *num);
    
    err = probe_list(sid,*num,probe_names,&num_found);

    if ( H5_RETURN_OK(err) && num_found > 0) {
      copy_ok =
       copy_c_array_to_fortran_array((const char * const *)probe_names,num_found,names,num,len);

      if ( DANU_RETURN_OK(copy_ok) && ( *num >= num_found ) ) {
	err = DANU_SUCCESS;
      }
      else {

	if ( ! DANU_RETURN_OK(copy_ok) ) {
	  DANU_ERROR_MESS("Failed to copy C strings to fortran array");
	}

	if ( *num < num_found ) {
	  DANU_ERROR_MESS("probes_names array too small to hold all probe"
			  " names");
	}

	err = DANU_FAILURE;
      }


    }
    else {

      if ( H5_RETURN_FAIL(err) ) {
        DANU_ERROR_MESS("Failed to read probe name list");
      }

      if ( num_found == 0 ) {
	DANU_WARN_MESS("No probes found");
      }

      err = DANU_FAILURE;
    }
  
    num_free = *num < num_found ? *num : num_found;
    for(i=0;i<num_free;i++)
      DANU_FREE(probe_names[i]);

    DANU_FREE(probe_names);


    error_translate(err,ierr);
        
}


void probe_data_count_f(hid_t_ptr *sptr, int * count, int *ierr)
{
    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;

    err = probe_count(sid,count);

    error_translate(err,ierr);

}

void probe_data_open_f(hid_t_ptr *sptr,
                       const char *probe_name,
                       const int *len,
                       hid_t_ptr *pptr,
                       int *ierr)
{
    hid_t sid = GET_HID_VALUE(sptr);
    hid_t pid;
    herr_t err;
    int  c_len        = *len+1;
    char *c_data       = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(probe_name,*len,c_data,c_len);
    
    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
      pid = probe_open_data(sid,c_data);
      if ( H5_ISA_VALID_ID(pid) ) {
	err = DANU_SUCCESS;
      }
      DANU_FREE(c_data);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }

    *pptr = create_hid_struct(pid);

    error_translate(err,ierr);

}

/* Create routines */
void probe_create_data_int_f(hid_t_ptr *sptr,
			     const char *probe_name,
			     const int  *flen,
			     const int  *len,
			     const int  *num,
			     const int  *data,
			     hid_t_ptr  *probe,
			     int        *ierr)

{
    hid_t sid = GET_HID_VALUE(sptr);
    hid_t pid;
    herr_t err;
    int  c_len        = *flen+1;
    char *c_data       = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(probe_name,*flen,c_data,c_len);

    if ( status == DANU_SUCCESS ) {
      err = probe_create_data_int(sid,c_data,*len,*num,data,&pid);
      *probe=create_hid_struct(pid);
      DANU_FREE(c_data);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }

    error_translate(err,ierr);

}

void probe_create_data_float_f(hid_t_ptr *sptr,
			       const char *probe_name,
			       const int  *flen,
			       const int  *len,
			       const int  *num,
			       const float  *data,
			       hid_t_ptr  *probe,
			       int        *ierr)

{
    hid_t sid = GET_HID_VALUE(sptr);
    hid_t pid;
    herr_t err;
    int  c_len        = *flen+1;
    char *c_data       = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(probe_name,*flen,c_data,c_len);

    if ( status == DANU_SUCCESS ) {
      err = probe_create_data_float(sid,c_data,*len,*num,data,&pid);
      *probe=create_hid_struct(pid);
      DANU_FREE(c_data);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }

    error_translate(err,ierr);

}

void probe_create_data_double_f(hid_t_ptr *sptr,
			       const char *probe_name,
			       const int  *flen,
			       const int  *len,
			       const int  *num,
			       const double  *data,
			       hid_t_ptr  *probe,
			       int        *ierr)

{
    hid_t sid = GET_HID_VALUE(sptr);
    hid_t pid;
    herr_t err;
    int  c_len        = *flen+1;
    char *c_data       = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(probe_name,*flen,c_data,c_len);

    if ( status == DANU_SUCCESS ) {
      err = probe_create_data_double(sid,c_data,*len,*num,data,&pid);
      *probe=create_hid_struct(pid);
      DANU_FREE(c_data);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }

    error_translate(err,ierr);

}






/* Write routines */

void probe_data_write_int_f(hid_t_ptr *pptr, 
			     const int *num,
			     const int *data,
			     int *ierr)
{
    hid_t pid = GET_HID_VALUE(pptr);
    herr_t err;

    err = probe_data_write_int(pid,*num,data);
    error_translate(err,ierr);

}

void probe_data_write_float_f(hid_t_ptr *pptr, 
			     const int *num,
			     const float *data,
			     int *ierr)
{
    hid_t pid = GET_HID_VALUE(pptr);
    herr_t err;
    
    err = probe_data_write_float(pid,*num,data);

    error_translate(err,ierr);

}

void probe_data_write_double_f(hid_t_ptr *pptr, 
			     const int *num,
			     const double *data,
			     int *ierr)
{
    hid_t pid = GET_HID_VALUE(pptr);
    herr_t err;

    
    err = probe_data_write_double(pid,*num,data);
    
    error_translate(err,ierr);

}

/* Read routines */
void probe_data_read_int_f(hid_t_ptr *pptr, 
			     int *data,
			     int *ierr)
{
    hid_t pid = GET_HID_VALUE(pptr);
    herr_t err;


    err = probe_data_read_int(pid,data);

    error_translate(err,ierr);

}

void probe_data_read_float_f(hid_t_ptr *pptr, 
			     float *data,
			     int *ierr)
{
    hid_t pid = GET_HID_VALUE(pptr);
    herr_t err;

    err = probe_data_read_float(pid,data);

    error_translate(err,ierr);

}

void probe_data_read_double_f(hid_t_ptr *pptr, 
			     double *data,
			     int *ierr)
{
    hid_t pid = GET_HID_VALUE(pptr);
    herr_t err;

    err = probe_data_read_double(pid,data);

    error_translate(err,ierr);

}

