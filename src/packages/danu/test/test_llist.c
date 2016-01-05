/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
#include <stdlib.h>
#include <stdio.h>

#include <danu_memory.h>
#include <danu_error.h>

#include <danu_linked_list.h>

extern danu_err_t find_tail(danu_node_t *node, void *op_data);

typedef struct person_t {
    char * name;
    int id;
} person_t;

/* Prototypes */
person_t * create_person(int id);
danu_err_t delete_person(void *data);
danu_err_t print_person(danu_node_t *node,void *op_data);

person_t * create_person(int id)
{
    char dummy_name[] = "John Q. Public";

    person_t * person = DANU_MALLOC(person_t, 1);

    person->name = DANU_MALLOC(char,128);
    sprintf(person->name,dummy_name);

    person->id = id;

    return person;
}

danu_err_t delete_person(void * data)
{
    person_t *person = (person_t *) data;
    DANU_FREE(person->name);
    DANU_FREE(person);

    return DANU_SUCCESS;
}

danu_err_t print_person(danu_node_t *node, void *op_data)
{
    person_t *person = (person_t *) node->data; 

    printf("ID=%d Name %s\n", person->id, person->name);

    return DANU_CONTINUE;
}



int main(int argc, char ** argv)
{
    danu_err_t flag;

    danu_list_t *list;
    person_t *person = create_person(0);
    danu_node_t *node;
    danu_node_t *a;
    int num_people = 10;
    int save_idx = 3;
    int n;

    if ( NULL == ( list = danu_list_create() ) ) {
        DANU_ERROR_MESS("Failed to allocate list");
        goto FAIL_EXIT;
    }

    if ( NULL == danu_list_define_root(list,person,&delete_person) ) {
        DANU_ERROR_MESS("Failed to define root node");
        goto FAIL_EXIT;
    }


    if ( argc > 1 ) {
        num_people = atoi(argv[1]);
    }
    n=1;
    while(n<=num_people) {
        person = create_person(n);
        node = danu_list_append(list,person,&delete_person);
        if ( n == save_idx ) {
            a = node;
        }
        n++;
    }

    person = create_person(n);
    node = danu_node_create(person,&delete_person);
    danu_node_insert_next(a,node);
    n++;
    person = create_person(n);
    node = danu_node_create(person,&delete_person);
    danu_node_insert_prev(a,node);

    danu_list_iterate(list,&print_person,NULL);


    danu_list_delete(list);

    return DANU_SUCCESS;

FAIL_EXIT:

     return DANU_FAILURE;
}



    


