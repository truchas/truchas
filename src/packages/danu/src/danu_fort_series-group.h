/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

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

