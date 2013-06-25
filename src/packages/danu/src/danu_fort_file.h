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

#ifndef DANU_FORT_FILE_H
#define DANU_FORT_FILE_H

#include <danu_fort_hid.h>

/* Use these defines to improve the readability of the code */
/* Basic file control */
#if 0
#include "Fortran2C.h"
#define danu_file_create_f        FORTRAN_FUNC_GLOBAL_(danu_file_create_f,DANU_FILE_CREATE_F)
#define danu_file_close_f         FORTRAN_FUNC_GLOBAL_(danu_file_close_f,DANU_FILE_CLOSE_F)
#define danu_file_open_rdonly_f   FORTRAN_FUNC_GLOBAL_(danu_file_open_rdonly_f,DANU_FILE_OPEN_RDONLY_F)
#define danu_file_open_rdwr_f     FORTRAN_FUNC_GLOBAL_(danu_file_open_rdwr_f,DANU_FILE_OPEN_RDWR_F)
#endif


/* Function prototypes corresponding interface definition in module_danu_iface.f90 */
void danu_file_create_f(const char *filename, const int *len, hid_t_ptr *ptr, int *ierr);
void danu_file_close_f(hid_t_ptr *ptr, int *ierr);
void danu_file_open_rdonly_f(const char *filename, const int *len, hid_t_ptr *ptr, int *ierr);
void danu_file_open_rdwr_f(const char *filename, const int *len, hid_t_ptr *ptr, int *ierr);

#endif /* DANU_FORT_FILE_H */

