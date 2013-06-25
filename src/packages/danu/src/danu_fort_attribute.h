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
* danu_fort_file.h
*
*  DANU FORTRAN/C interfaces for files
*
*/

#ifndef DANU_FORT_ATTR_H
#define DANU_FORT_ATTR_H

#include <danu_fort_hid.h>

/* Use these defines to improve the readability of the code */
/* Basic file control */
#if 0
#include "Fortran2C.h"
#define danu_attr_exists_f            FORTRAN_FUNC_GLOBAL_(danu_attr_exists_f,DANU_ATTR_EXISTS_F)
#define danu_attr_count_f            FORTRAN_FUNC_GLOBAL_(danu_attr_count_f,DANU_ATTR_CREATE_F)
#define danu_attr_names_f            FORTRAN_FUNC_GLOBAL_(danu_attr_names_f,DANU_ATTR_NAMES_F)

#define danu_attr_write_int_f        FORTRAN_FUNC_GLOBAL_(danu_attr_write_int_f,DANU_ATTR_WRITE_INT_F)
#define danu_attr_write_real4_f      FORTRAN_FUNC_GLOBAL_(danu_attr_write_real4_f,DANU_ATTR_WRITE_REAL4_F)
#define danu_attr_write_real8_f      FORTRAN_FUNC_GLOBAL_(danu_attr_write_real8_f,DANU_ATTR_WRITE_REAL8_F)
#define danu_attr_write_char_f       FORTRAN_FUNC_GLOBAL_(danu_attr_write_char_f,DANU_ATTR_WRITE_CHAR_F)

#define danu_attr_read_int_f         FORTRAN_FUNC_GLOBAL_(danu_attr_read_int_f,DANU_ATTR_READ_INT_F)
#define danu_attr_read_real4_f       FORTRAN_FUNC_GLOBAL_(danu_attr_read_real4_f,DANU_ATTR_READ_REAL4_F)
#define danu_attr_read_real8_f       FORTRAN_FUNC_GLOBAL_(danu_attr_read_real8_f,DANU_ATTR_READ_REAL8_F)
#define danu_attr_read_char_f        FORTRAN_FUNC_GLOBAL_(danu_attr_read_char_f,DANU_ATTR_READ_CHAR_F)
#endif


/* Function prototypes corresponding interface definition in module_danu_iface.f90 */

/* Attribute exists, count and names */
void danu_attr_exists_f(hid_t_ptr *ptr, const char *name, int *flen, int *exists, int *ierr);
void danu_attr_count_f(hid_t_ptr *ptr, int *num_found, int *ierr); 
void danu_attr_names_f(hid_t_ptr *ptr, char *names, int *flen, int *num, int *ierr);


/* Attribute write */
void danu_attr_write_int_f(hid_t_ptr *ptr, const char *name, int *name_flen, int *value, int *ierr);
void danu_attr_write_real4_f(hid_t_ptr *ptr, const char *name, int *name_flen, float *value, int *ierr);
void danu_attr_write_real8_f(hid_t_ptr *ptr, const char *name, int *name_flen, double *value, int *ierr);
void danu_attr_write_char_f(hid_t_ptr *ptr, const char *name, int *name_flen, const char *string, int *str_flen, int *ierr);

/* Attribute read */
void danu_attr_read_int_f(hid_t_ptr *ptr, const char *name, int *name_flen, int *buffer, int *ierr);
void danu_attr_read_real4_f(hid_t_ptr *ptr, const char *name, int *name_flen, float *buffer, int *ierr);
void danu_attr_read_real8_f(hid_t_ptr *ptr, const char *name, int *name_flen, double *buffer, int *ierr);
void danu_attr_read_char_f(hid_t_ptr *ptr, const char *name, int *name_flen, char *buffer, int *buf_flen, int *ierr);



#endif /* DANU_FORT_ATTR_H */

