/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * h5_objects.c
 *
 *  Purpose:
 *
 *         Danu defines an HDF5 struct called h5_object_t. This struct is a
 *         a container for HDF5 identifiers and parameter list identifiers
 *         This code is the memory manger for these structs.
 */

#include <string.h>

#include <hdf5.h>

#include <danu_memory.h>
#include <danu_error.h>
#include <danu_types.h>
#include <danu_utils.h>

#include <danu_h5_object.h>


/* PRIVATE define's */

typedef struct {

    h5o_node_list_t *list;
    H5O_type_t       type;
    char *           target_name;
} op_search_data_t;



/* PRIVATE FUNCTION PROTOTYPES */

/* Operators */
herr_t obj_type_search(hid_t loc, const char *name, const H5L_info_t *link_info, void * data);
herr_t obj_type_search_by_name(hid_t loc, const char *name, const H5L_info_t *link_info, void * data);

h5o_node_t * h5_object_node_create(const char * name, haddr_t addr)
{
    h5o_node_t *ptr = NULL;

    size_t str_len;

    if ( DANU_BAD_STRING(name) ) {
        DANU_ERROR_MESS("Invalid pointer to name");
        return ptr;
    }

    if ( NULL == (ptr = DANU_MALLOC(h5o_node_t,1) ) ) {
        DANU_ERROR_MESS("Failed to allocate the node mempry");
        return ptr;
    }

    ptr->addr = addr;

    str_len = strlen(name);
    ptr->name = DANU_MALLOC(char,str_len+1);
    strcpy(ptr->name,name);

    ptr->id = 0;
    ptr->prev = NULL;
    ptr->next = NULL;

    return ptr;
}

void h5_object_node_delete(h5o_node_t *node)
{
    h5o_node_t *prev;
    h5o_node_t *next;
    h5o_node_t *ptr;

    if ( DANU_BAD_PTR(node) ) {
        DANU_ERROR_MESS("Invalid node pointer");
        return;
    }

    prev = node->prev;
    next = node->next;

    /* Update the id's */
    ptr = next;
    while ( ptr != NULL ) {
        ptr->id--;
        ptr = ptr->next;
    }
    
    /* Update the nodes bounding the current node */
    if ( prev != NULL ) {
        prev->next = next;
    }

    if ( next != NULL ) {
        next->prev = prev;
    }

    DANU_FREE(node->name);
    DANU_FREE(node);

}



h5o_node_list_t * h5_object_list_create(H5O_type_t type)
{
    h5o_node_list_t * ptr = NULL;

    /* Check input */
    if ( type < 0  || type >= H5O_TYPE_NTYPES ) {
        DANU_ERROR_MESS("Invalid HDF5 object type");
        return ptr;
    }

    if ( NULL == (ptr = DANU_MALLOC(h5o_node_list_t,1) ) ) {
        DANU_ERROR_MESS("Failed to allocate ptr to node list");
        return ptr;
    }

    ptr->type = type;
    ptr->root = NULL;
    ptr->tail = NULL;

    return ptr;
    
}

void h5_object_list_delete(h5o_node_list_t *list)
{
    h5o_node_t *ptr, *next;

    /* Starting at the root delete each node */
    ptr = list->root;
    
    while( ptr != NULL ) {
        next = ptr->next;
        h5_object_node_delete(ptr);
        ptr = next;
    }

    DANU_FREE(list);

}

h5o_node_t * h5_object_list_append_node(h5o_node_list_t *list, const char * name, haddr_t addr)
{
    h5o_node_t *ptr = NULL;

    h5o_node_t *tail, *root;

    /* Check input */
    if ( DANU_BAD_PTR(list) ) {
        DANU_ERROR_MESS("Invalid list pointer");
        return ptr;
    }

    if ( NULL == ( ptr = h5_object_node_create(name,addr) ) ) {
        DANU_ERROR_MESS("Node allocate failed");
        return ptr;
    }

    tail = list->tail;
    root = list->root;

    if ( root == NULL ) {
        list->root = ptr;
        
        ptr->id = 0;
        ptr->prev = NULL;
        ptr->next = tail;
    }
    else {

        if ( tail == NULL ) {

            ptr->id = 1;
            ptr->prev = root;
            ptr->next = NULL;

            root->next = ptr;

            list->tail = ptr;

        }
        else {

            ptr->id   = 1 + tail->id;
            ptr->prev = tail;
            ptr->next = NULL;

            tail->next = ptr;

            list->tail = ptr;

        }

    }

    return ptr;
}

unsigned h5_object_list_nnodes(const h5o_node_list_t *list)
{

    unsigned ret = 0;

    h5o_node_t *root,*tail;

    root = list->root;
    tail = list->tail;

    if ( root != NULL ) {

        if ( tail != NULL ) {

            ret = tail->id + 1;

        }
        else {

            ret = 1;

        }

    }

    return ret;

}

int h5_object_list_copy_names(const h5o_node_list_t *list, int num, char **names)
{
    int ncpy = -1;
    size_t str_len;
    int i;
    h5o_node_t *node, *next;

    /* Check input */
    if (DANU_BAD_PTR(list)) {
        DANU_ERROR_MESS("Invalid object list pointer");
        return ncpy;
    }

    if ( num <= 0 ) {
        DANU_ERROR_MESS("Invalid number of pointers");
        return ncpy;
    }

    if ( DANU_BAD_PTR(names) ) {
        DANU_ERROR_MESS("Invalid names pointer");
        return ncpy;
    }

    /* Loop through the list */
    node = list->root;
    i = 0;
    ncpy = 0;
    while( (i<num) && (node != NULL) ) {
      next = node->next;
      str_len = strlen(node->name);
      str_len++;
      names[i] = DANU_MALLOC(char,str_len);
      strcpy(names[i],node->name);
      ncpy++;
      i++;
      node = next;
    }

    if ( i<num ) {
      while ( i<num ) {
	names[i] = DANU_CALLOC(char,1);
	i++;
      }
    }

    return ncpy;
}


/* PRIVATE FUNCTIONS */

/* Search by types either groups or datasets */
herr_t obj_type_search(hid_t loc, const char *name, const H5L_info_t *link_info, void * data)
{
    op_search_data_t * op_data = (op_search_data_t *)data;

    h5o_node_list_t *list        = op_data->list;
    H5O_type_t       search_type = op_data->type; 

    H5O_info_t obj_info;
    herr_t status;


    status =  H5Oget_info_by_name(loc,name,&obj_info,H5P_DEFAULT);
    if ( status >= 0 ) {
        
        if ( obj_info.type == search_type ) {
            if ( NULL == h5_object_list_append_node(list,name,obj_info.addr) ) {
                DANU_ERROR_MESS("Failed to add node to search list");
                status = DANU_FAILURE;
            }
        }

    }

    return status;
}

/* Search by type AND name */
herr_t obj_type_search_by_name(hid_t loc, const char *name, const H5L_info_t *link_info, void * data)
{
    op_search_data_t * op_data = (op_search_data_t *)data;

    h5o_node_list_t *list        = op_data->list;
    H5O_type_t       search_type = op_data->type; 
    char *           target_name = op_data->target_name;

    H5O_info_t obj_info;
    herr_t status;

    status =  H5Oget_info_by_name(loc,name,&obj_info,H5P_DEFAULT);
    if ( status >= 0 ) {
        
        if ( ( obj_info.type == search_type ) && ( strcmp(name,target_name) == 0 ) ) {
            if ( NULL == h5_object_list_append_node(list,name,obj_info.addr) ) {
                DANU_ERROR_MESS("Failed to add node to search list");
                status = DANU_FAILURE;
            }

        }

    }

    return status;
}


/* PUBLIC FUNCTIONS */
herr_t h5_object_search(hid_t loc_id, H5O_type_t type, h5o_node_list_t *list) 
{
    herr_t status = DANU_FAILURE;
    
    op_search_data_t  data;
    hsize_t idx = 0;

    /* Check Input */
    if ( H5_ISA_INVALID_ID(loc_id) ) {
        DANU_ERROR_MESS("Invalid HDF5 identifier");
        return status;
    }

    if ( type >= H5O_TYPE_NTYPES || type < 0) {
        DANU_ERROR_MESS("Invalid HDF5 object type");
        return status;
    }

    if ( type != H5O_TYPE_GROUP && type != H5O_TYPE_DATASET ) {
        DANU_ERROR_MESS("Only groups and dataset are allowed in the search");
        return status;
    }

    if ( DANU_BAD_PTR(list) ) {
        DANU_ERROR_MESS("Invalid list pointer");
        return status;
    }
    

    /* Now iterate */
    data.list = list;
    data.type = type;
    status = H5Literate(loc_id, H5_INDEX_NAME, H5_ITER_NATIVE, &idx, obj_type_search,(void*)&data);

    return status;
}

herr_t h5_object_search_by_name(hid_t loc_id, const char *target_name, H5O_type_t type, h5o_node_list_t *list)
{
    herr_t status = DANU_FAILURE;
   
    size_t len;
    op_search_data_t  data;

    /* Check Input */
    if ( H5_ISA_INVALID_ID(loc_id) ) {
        DANU_ERROR_MESS("Invalid HDF5 identifier");
        return status;
    }

    if ( DANU_BAD_STRING(target_name) ) {
        DANU_ERROR_MESS("Invalid string pointer")
        return status;
    }

    if ( type >= H5O_TYPE_NTYPES || type < 0) {
        DANU_ERROR_MESS("Invalid HDF5 object type");
        return status;
    }

    if ( type != H5O_TYPE_GROUP && type != H5O_TYPE_DATASET ) {
        DANU_ERROR_MESS("Only groups and dataset are allowed in the search");
        return status;
    }

    if ( DANU_BAD_PTR(list) ) {
        DANU_ERROR_MESS("Invalid list pointer");
        return status;
    }
   
    /* Now iterate */
    data.list        = list;
    data.type        = type;

    len = strlen(target_name);
    data.target_name = DANU_MALLOC(char,len+1);
    strcpy(data.target_name,target_name);
    status = H5Literate(loc_id, H5_INDEX_NAME, H5_ITER_NATIVE, 
			NULL, obj_type_search_by_name,(void*)&data);
    DANU_FREE(data.target_name);

    return status;
}
 


