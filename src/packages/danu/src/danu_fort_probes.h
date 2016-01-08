/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
*
*  DANU Fortran/C Probe Data Interface
*
*/

#ifndef DANU_FORT_PROBES_H
#define DANU_FORT_PROBES_H

#include <danu_fort_hid.h>

/* Probe Queries */
void probe_count_f(hid_t_ptr *sptr,
                 int *count,
                 int *ierr);

void probe_data_open_f(hid_t_ptr *sptr,
                     const char *probe_name,
                     const int *len,
                     hid_t_ptr *pptr,
                     int *ierr);

void probe_create_data_int_f(hid_t_ptr *sptr,
                             const char *probe_name,
                             const int  *flen,
                             const int  *len,
                             const int  *num,
                             const int  *data,
                             hid_t_ptr  *pptr,
                             int        *ierr);

void probe_data_write_int_f(hid_t_ptr *pptr,
                            const int  *num,
                            const int  *data,
                            int        *ierr);

#endif

