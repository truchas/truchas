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
* danu_fort_non-series.h
*
*  Output Simulation Non-series data Fortran/C interface
*
*/

#ifndef DANU_FORT_NONSERIES_H
#define DANU_FORT_NONSERIES_H

#include <danu_fort_hid.h>
/* Use these defines to improve the readability of the code */
/* Basic simulation control */
#if 0
#include "Fortran2C.h"
#define data_exists_f            FORTRAN_FUNC_GLOBAL_(data_exists_f, DATA_EXISTS_F)
#define data_count_f             FORTRAN_FUNC_GLOBAL_(data_count_f, DATA_COUNT_F)
#define data_list_f              FORTRAN_FUNC_GLOBAL_(data_list_f, DATA_LIST_F)
#endif

/* Public Functions */
void data_exists_f(const hid_t_ptr *sid, const char *dataname, const int *flen, int *exists, int *ierr);

void data_count_f(const hid_t_ptr *sid, int *ndata, int *ierr);
void data_list_f(const hid_t_ptr *sid, char *datanames, const int *flen, const int *num, int *ierr);

void data_write_int_f(const hid_t_ptr *sid, const char * name, const int *flen,
                      const int *dim, const int *dimensions, int *data, int * ierr);
void data_write_float_f(const hid_t_ptr *sid, const char * name, const int *flen,
                        const int *dim, const int *dimensions, float *data, int *ierr);
void data_write_double_f(const hid_t_ptr *sid, const char * name, const int *flen,
                         const int *dim, const int *dimensions, double *data, int *ierr);
void data_write_string_f(const hid_t_ptr *sid, const char * name, const int *flen,
                         const int *dim, const int *dimensions, char *data, int *ierr);

void data_read_int_f(const hid_t_ptr *sid, const char * name, const int *flen,
                     const int *dim, const int * dimensions, int *data, int *ierr);
void data_read_float_f(const hid_t_ptr *sid, const char * name, const int *flen, const int *dim, const int * dimensions, float *data, int *ierr);
void data_read_double_f(const hid_t_ptr *sid, const char * name, const int *flen, const int *dim, const int * dimensions, double *data, int *ierr);
void data_read_string_f(const hid_t_ptr *sid, const char * name, const int *flen, const int *dim, const int * dimensions, char *data, int *ierr);



#endif
