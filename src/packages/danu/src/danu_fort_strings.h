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

#ifndef DANU_FORT_STRINGS_H
#define DANU_FORT_STRINGS_H

#include <danu.h>

int compute_fortran_trim_len(const char * f_str, const int f_len);

danu_err_t convert_string_f2c(const char *fort_str, int f_len, char *c_str, int c_len);    
danu_err_t convert_string_c2f(const char *c_str,  char *f_str, int f_len);    

danu_err_t copy_c_array_to_fortran_array(const char * const *c_array, int num,  
                                         char *fort_data, const int *fnum, const int *flen);

char ** copy_fortran_array_to_c_array(const char * fort_data, const int *fnum, const int *flen);

#endif /* DANU_FORT_STRINGS_H */

