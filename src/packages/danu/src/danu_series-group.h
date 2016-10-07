/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_series-group.h
*
*  Telluride Output Simulation Module 
*
*/

#ifndef TOUT_SERIES_GROUP_H
#define TOUT_SERIES_GROUP_H

#include <hdf5.h>

/* Public define's */
#define SERIES_GROUP_NAME     "Series Data"

#define SEQUENCE_BASE_NAME      "Series"
#define SEQ_TIME_ATTR_NAME      "time"
#define SEQ_SEQNUM_ATTR_NAME    "sequence number"
#define SEQ_CYCNUM_ATTR_NAME    "cycle"

#define DANU_VALID_SEQNUM(a)       ( (a) > 0 ? 1 : 0 )
#define DANU_INVALID_SEQNUM        -1


/* Public Functions */

herr_t sequence_create_root_group(hid_t sid);
hid_t  sequence_open_root_group(hid_t sid);

herr_t sequence_exists(hid_t sid, const char *seriesname, int *exists);


herr_t sequence_name_max_bytes(hid_t sid, size_t*max_bytes);
herr_t sequence_count(hid_t sid, int *nseries);
herr_t sequence_list(hid_t sid, int num, char **datanames, int *num_found);

herr_t sequence_get_handle(hid_t sid, const char * seriesname, hid_t *nsid);
herr_t sequence_get_handle_byid(hid_t sid, int id, hid_t *nsid);
herr_t sequence_getNextID(hid_t sid, int cycle, double t, hid_t *nsid);

/* Return the series name given an id */
char * sequence_get_name(int num);

/* Sequence group attributes */
herr_t    sequence_get_time(hid_t nsid, double *time);
herr_t    sequence_get_cycle(hid_t nsid, int *cycle);      
herr_t    sequence_get_id(hid_t nsid, int *id);      

#endif
