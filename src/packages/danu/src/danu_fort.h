/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_fort.h
*
*  DANU FORTRAN C interfaces
*
*/

#ifndef DANU_FORT_H
#define DANU_FORT_H

#include <danu.h>

/* Create a struct to hold the hid_t, hsize_t arrays and herr_t values */
typedef struct hid_t_struct_ {
    hid_t     id;
} hid_t_struct;


/* Subroutines to create, destroy and manipulate the structs above */
hid_t_struct *  create_hid_ptr(hid_t id);
void            destroy_hid_ptr(hid_t_struct *ptr);


/* Struct to handle Fortran Strings */
typedef struct f_str_t {
    char    *data;
    size_t   bytes;
    int *    trim_len;
    int      len;
    int      num;
} danu_fort_char_t;

danu_fort_char_t * create_fortran_char_data(char *fort_data, int len, const int num);
void               delete_fortran_char_data(danu_fort_char_t *ptr);


/* Use these defines to improve the readability of the code */

/* Basic file control */
#define danu_file_create_f        FORTRAN_FUNC_GLOBAL_(danu_file_create_f,DANU_FILE_CREATE_F)
#define danu_file_close_f         FORTRAN_FUNC_GLOBAL_(danu_file_close_test_f,DANU_FILE_CLOSE_TEST_F)


/* Function prototypes corresponding interface definition in module_danu_iface.f90 */
//void danu_file_create_f(const char *filename, const int *len, hid_t_ptr *ptr, int *ierr);
//void danu_file_close_f(hid_t_ptr *ptr, int *ierr);


#endif /* DANU_FORT_H */

