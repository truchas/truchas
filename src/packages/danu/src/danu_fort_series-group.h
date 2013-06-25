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
*  DANU Fortran/C Series Group Interface
*
*/

#ifndef DANU_FORT_SERIES_GROUP_H
#define DANU_FORT_SERIES_GROUP_H

#include <danu_fort_hid.h>
void sequence_exists_f(const hid_t_ptr *sptr, const char * sname, const int *len, int *flag, int *ierr);
void sequence_count_f(const hid_t_ptr *sptr, int * num, int *ierr);

void sequence_get_next_id_f(const hid_t_ptr *sptr, const int *cycle, const double *t, hid_t_ptr *nid, int *ierr); 
void sequence_list_f(const hid_t_ptr *sptr, char *names, const int *len, const int *num, int *ierr);

void sequence_get_id_f(const hid_t_ptr *sptr, const char *sname, const int *len, hid_t_ptr *nid, int *ierr);
void sequence_get_id_bynum_f(const hid_t_ptr *sptr, const int * num, hid_t_ptr *nptr, int *ierr);

#endif

