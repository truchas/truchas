/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_h5_error.h
*
*  DANU Error Handler
*
*/

#ifndef DANU_H5_OBJECT_H
#define DANU_H5_OBJECT_H


#include <hdf5.h>

typedef struct h5o_node_t {
    unsigned     id;
    char *       name;
    haddr_t      addr;
    struct h5o_node_t * prev;
    struct h5o_node_t * next;
} h5o_node_t;

typedef struct h5o_node_list_t {
    struct h5o_node_t *  root;
    struct h5o_node_t *  tail;
    H5O_type_t    type;
} h5o_node_list_t;


h5o_node_list_t * h5_object_list_create(H5O_type_t type);
void              h5_object_list_delete( h5o_node_list_t *list);
h5o_node_t *      h5_object_list_append_node( h5o_node_list_t *list, const char *name, haddr_t addr);
unsigned          h5_object_list_nnodes(const h5o_node_list_t *list);
int               h5_object_list_copy_names(const h5o_node_list_t *list, int num, char **names);

h5o_node_t *      h5_object_node_create(const char * name, haddr_t addr);
void              h5_object_node_delete( h5o_node_t *node);

herr_t h5_object_search(hid_t loc_id, H5O_type_t type,  h5o_node_list_t *list);
herr_t h5_object_search_by_name(hid_t loc_id, const char * target, H5O_type_t type,  h5o_node_list_t *list);

#endif


