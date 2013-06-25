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
 * group.c
 *
 *  DANU Groups 
 *
 *
 *  Purpose:
 *
 *          The HDF5 library
 *
 *
 */

#ifndef DANU_GROUP_H
#define DANU_GROUP_H

#include <hdf5.h>

hbool_t danu_group_exists(hid_t loc_id, const char * name);

hid_t danu_group_create(hid_t loc, const char * name);
hid_t danu_group_open(hid_t loc, const char *name);
herr_t danu_group_close(hid_t id);

herr_t  danu_group_get_info(hid_t gid, H5G_info_t *info);
herr_t  danu_group_get_info_by_name(hid_t loc, const char *group, H5G_info_t *info);

herr_t danu_group_get_nlinks(hid_t gid, hsize_t *nlinks);
herr_t danu_group_get_nlinks_by_name(hid_t loc, const char * group, hsize_t *nlinks);

herr_t danu_group_get_subgroups(hid_t gid, int num, char **subgroups, int *num_found);

herr_t danu_group_get_datasets(hid_t gid, int num, char **datasets, int *num_found);

herr_t danu_group_find_target(hid_t gid, const char * target, H5O_type_t type, int *found);
herr_t danu_group_find_target_subgroup(hid_t gid, const char * target, int *found);
herr_t danu_group_find_target_dataset(hid_t gid, const char * target, int *found);


/* THIS IS BROKEN DO NOT USE */
hid_t danu_group_create_external(hid_t link_loc, 
                                 const char * group_link,
                                 const char * file,
                                 const char * obj_name);

#endif
