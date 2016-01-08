/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
*
*  DANU Fortran/C File Interface
*
*/
#include<string.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_memory.h>
#include <danu_error.h>
#include <danu_file.h>
#include <danu_fort_hid.h>
#include <danu_fort_error.h>
#include <danu_fort_strings.h>
#include <danu_fort_file.h>

/* Basic file control */
void danu_file_create_f(const char *filename, const int *len, hid_t_ptr *ptr, int *ierr)
{
    hid_t id;
    herr_t err;
    int  c_len        = *len+1;
    char *c_file      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(filename,*len,c_file,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
        id = danu_file_create(c_file);

        if ( H5_ISA_VALID_ID(id) ) {
            err = DANU_SUCCESS;
        }
        
        DANU_FREE(c_file);
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }    
     
    *ptr  = create_hid_struct(id);

    error_translate(err,ierr);
    
}

void danu_file_close_f(hid_t_ptr *ptr, int *ierr)
{
    hid_t id = GET_HID_VALUE(ptr);
    herr_t err = danu_file_close(id);

    destroy_hid_struct(*ptr);
    error_translate(err,ierr);
}

void danu_file_open_rdonly_f(const char *filename, const int *len, hid_t_ptr *ptr, int *ierr)
{
    hid_t id;
    herr_t err;
    int  c_len        = *len+1;
    char *c_file      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(filename,*len,c_file,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
        id = danu_file_open_rdonly(c_file);

        if ( H5_ISA_VALID_ID(id) ) {
            err = DANU_SUCCESS;
        }
        
        DANU_FREE(c_file);
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }    
     
    *ptr  = create_hid_struct(id);

    error_translate(err,ierr);
    
}

void danu_file_open_rdwr_f(const char *filename, const int *len, hid_t_ptr *ptr, int *ierr)
{
    hid_t id;
    herr_t err;
    int  c_len        = *len+1;
    char *c_file      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(filename,*len,c_file,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
        id = danu_file_open_rdwr(c_file);

        if ( H5_ISA_VALID_ID(id) ) {
            err = DANU_SUCCESS;
        }
        
        DANU_FREE(c_file);
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }    
     
    *ptr  = create_hid_struct(id);

    error_translate(err,ierr);
    
}


