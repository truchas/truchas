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
 * linked_list.c
 *
 *  DANU linked list data structure
 *
 *
 *  Purpose:
 *
 */

#include <danu_memory.h>
#include <danu_error.h>

#include <danu_linked_list.h>

/* Private function prototypes */
danu_node_t * danu_node_create(void * data, danu_node_data_destruct_t func)
{
    danu_node_t * node = DANU_MALLOC(danu_node_t, 1);

    node->next = NULL;
    node->prev = NULL;
    node->data = data;
    node->destruct = func;

    return node;
}

void danu_node_delete(danu_node_t * node)
{
    danu_node_t *next, *prev;

    next  = node->next;
    prev  = node->prev;

    /* Hook up the next node to prev */
    if ( next != NULL ) {
        next->prev = prev;
    }

    if ( prev != NULL ) {
        prev->next = next;
    }

    /* Now delete the node */
    node->destruct(node->data);

    DANU_FREE(node);

}

void danu_node_insert_next(danu_node_t * node, danu_node_t * new_node)
{
    danu_node_t *next = node->next;

    node->next = new_node;
    new_node->prev = node;
    new_node->next = next;

    if ( next != NULL ) {
        next->prev = new_node;
    }

}

void danu_node_insert_prev(danu_node_t * node, danu_node_t * new_node)
{
    danu_node_t *prev = node->prev;

    node->prev = new_node;
    new_node->next = node;
    new_node->prev = prev;

    if ( prev != NULL ) {
        prev->next = new_node;
    }

}

danu_list_t * danu_list_create(void)
{
    danu_list_t *list = DANU_MALLOC(danu_list_t, 1);

    list->root = NULL;

    return list;
}


void danu_list_delete(danu_list_t *list)
{
    danu_node_t *root = list->root;
    danu_node_t *node, *tmp;

    if ( root != NULL ) {
        node = root->next;
        while( node != NULL ) {
            tmp = node->next;
            danu_node_delete(node);
            node = tmp;
        }

        danu_node_delete(root);
    }

    DANU_FREE(list);
}

int danu_list_iterate(danu_list_t *list, danu_list_iter_t func, void * op_data)
{
    int return_code = DANU_FAILURE;
    int loop_ctrl = DANU_CONTINUE;
    danu_node_t *node = list->root;
    danu_node_t *tmp;

    if ( node != NULL ) {
        loop_ctrl = DANU_CONTINUE;
        while ( node != NULL && loop_ctrl == DANU_CONTINUE ) {
            tmp = node->next;
            loop_ctrl = func(node, op_data);
            node = tmp;
        }

        return_code = loop_ctrl;
    }
    else {
        DANU_ERROR_MESS("Empty list can not iterate");
        return_code = DANU_FAILURE;
    }

    return return_code;
}

danu_node_t * danu_list_find_tail(danu_list_t *list)
{
    danu_node_t *node = list->root;
    danu_node_t *tail = NULL;

    while (node != NULL ) {
        if ( node->next == NULL ) {
            tail = node;
        }
        node=node->next;
    }

    return tail;
}

danu_node_t * danu_list_append(danu_list_t *list, void *data, danu_node_data_destruct_t func)
{
    danu_node_t * node;
    danu_node_t * tail;

    tail = danu_list_find_tail(list);
    if ( tail == NULL ) {
        node = danu_list_define_root(list,data,func);
    }
    else {
        if ( NULL != ( node = danu_node_create(data,func) ) ) {
            danu_node_insert_next(tail,node);
            node->next = NULL;
        }
        else {
            node = NULL;
        }
    }

    return node;
}

danu_node_t * danu_list_define_root(danu_list_t *list, void * data, danu_node_data_destruct_t func)
{
    danu_node_t * node = danu_node_create(data,func);
    danu_node_t * root = list->root;
   
   
    if ( root != NULL ) {
        danu_node_insert_prev(root,node);
    }

    list->root = node;

    return node;
}

int danu_list_node_count(const danu_list_t *list)
{
  danu_node_t *root, *ptr;
  int n = 0;

  root = list->root;
  if ( root != NULL ) {
    ptr = root;
    n = 1;
    while (ptr->next != NULL ) {
      n++;
      ptr = ptr->next;
    }
  }

  return n;
}




            
    





