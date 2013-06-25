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
* danu_utils.h
*
*  DANU General utility functions
*
*/

#ifndef DANU_UTILS_H
#define DANU_UTILS_H

#include <unistd.h>

#include<hdf5.h>

#include <danu_error.h>


int danu_file_access(const char * file, int access_mode);
danu_err_t danu_file_delete(const char *name);
danu_err_t danu_rand_data_int(int min, int max, size_t size, int * buffer); 
danu_err_t danu_rand_data_double(double min, double max, size_t  size, double * buffer); 

herr_t     convert_int_to_hsize(int size, const int *array, hsize_t *out);
danu_err_t convert_int_to_size(int size, const int *array, size_t *out);

void danu_print_hid_info(hid_t hid);

/*
 * Macros
 *
 * FILE_EXISTS(FILE)   Checks existence of FILE
 * FILE_READ_OK(FILE)  Checks read permission of FILE
 * FILE_WRITE_OK(FILE) Checks write access of FILE
 *
 * F_OK, R_OK and W_OK are macros defined in unistd.h
 *
 */
#define FILE_EXISTS(file)        ( danu_file_access((file),F_OK) )
#define FILE_READ_OK(file)       ( danu_file_access((file),R_OK) )
#define FILE_WRITE_OK(file)      ( danu_file_access((file),W_OK) )

/* Useful MACROS */
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define MAX(a,b)                 ( (a) <= (b) ? (b) : (a) )
#define MIN(a,b)                 ( (a) <= (b) ? (a) : (b) )

#endif

