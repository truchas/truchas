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
*
*
*/

#ifndef TOUT_OUTPUT_H
#define TOUT_OUTPUT_H

#include <hdf5.h>

/* Public defines */

/* Public Functions */
hbool_t output_file_is_valid(hid_t fid);

hid_t  output_file_open(const char *filename, unsigned access, unsigned action);
herr_t output_file_create(const char *filename, hid_t *fid);
herr_t output_file_open_rdonly(const char *filename, hid_t *fid);
herr_t output_file_open_append(const char *filename, hid_t *fid);
herr_t output_file_open_rdwr(const char * filename, hid_t *fid);
herr_t output_file_close(hid_t *fid);

#endif
