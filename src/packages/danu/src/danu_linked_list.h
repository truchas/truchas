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
 * danu_linked_list.h
 *
 *  DANU linked list manager
 *
 */

#ifndef DANU_LL_H
#define DANU_LL_H

#include <danu_error.h>

/* Basic data structures */

/* Node destructor function */
typedef danu_err_t (*danu_node_data_destruct_t)(void * data);

typedef struct danu_node_t {
    struct danu_node_t * next;
    struct danu_node_t * prev;

    danu_node_data_destruct_t destruct;
    void *data;
} danu_node_t;

typedef struct danu_list_t {
    struct danu_node_t * root;

} danu_list_t;

/* Function pointer that the iterator calls */
typedef danu_err_t (*danu_list_iter_t)(danu_node_t * node, void * op_data);

danu_node_t * danu_node_create(void * data, danu_node_data_destruct_t destruct);
void          danu_node_delete(danu_node_t *node);
void          danu_node_insert_next(danu_node_t * node, danu_node_t * new_node);
void          danu_node_insert_prev(danu_node_t * node, danu_node_t * new_node);


danu_list_t * danu_list_create(void);
void          danu_list_delete(danu_list_t * list);

int danu_list_iterate(danu_list_t *list, danu_list_iter_t func, void *op_data);

danu_node_t * danu_list_append(danu_list_t *list, void * data, danu_node_data_destruct_t func);
danu_node_t * danu_list_define_root(danu_list_t *list, void *data, danu_node_data_destruct_t func);

int danu_list_node_count(const danu_list_t *list);




#endif
