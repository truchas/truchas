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
* danu_probe.h
*
*  Telluride Output Probe Module 
*
*/

#ifndef TOUT_PROBES_H
#define TOUT_PROBES_H

#include <hdf5.h>

/* Public define's */

#define PROBES_GROUP_NAME "Probes"
#define PROBE_NUM_IDX 0
#define PROBE_LEN_IDX 1

/* Public Functions */

herr_t probe_exists(hid_t sid, const char *probename, int *exists);

herr_t probe_create_group(hid_t sid);
hid_t  probe_open_group(hid_t sid);
hid_t  probe_open_data(hid_t sid, const char *probename);

herr_t  probe_create_data(hid_t sid, const char *probe_name, 
                         hid_t type, int len, int num, const void *data, hid_t *pid);
herr_t probe_create_data_int(hid_t sid, const char * probe_name,
                             int len, int num, const int * data, hid_t *pid);
herr_t probe_create_data_float(hid_t sid, const char * probe_name,
                               int len, int num, const float * data, hid_t *pid);
herr_t probe_create_data_double(hid_t sid, const char * probe_name,
                                int len, int num, const double * data, hid_t *pid);

herr_t probe_count(hid_t sid, int * num);
herr_t probe_list(hid_t sid, int num, char **probenames, int *num_found);

herr_t probe_data_length(hid_t sid, const char *pname, int *len);
herr_t probe_data_num(hid_t sid, const char *pname, int *num);
herr_t probe_data_type(hid_t sid, const char *pname, int *typecode);

herr_t probe_data_length2(hid_t pid, int *len);
herr_t probe_data_num2(hid_t pid, int *num);
herr_t probe_data_type2(hid_t pid, int *typecode);

herr_t probe_attribute_list(hid_t sid, const char *probename, int num, char **attribs, int *num_found);
herr_t probe_attribute_list2(hid_t pid, int num, char **attribs, int *num_found);

herr_t probe_data_write(hid_t pid, hid_t type, int num, const void * data);
herr_t probe_data_write_int(hid_t pid, int num, const int * data);
herr_t probe_data_write_double(hid_t pid, int num, const double * data);
herr_t probe_data_write_float(hid_t pid, int num, const float * data);

herr_t probe_data_read(hid_t pid, hid_t type, void * data);
herr_t probe_data_read_int(hid_t pid, int * data);
herr_t probe_data_read_double(hid_t pid, double * data);
herr_t probe_data_read_float(hid_t pid, float * data);


#endif
