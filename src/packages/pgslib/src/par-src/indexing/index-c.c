/* Support for the indexing routines */

/* $Id: index-c.c,v 1.2 2001/03/22 00:26:13 ferrell Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "pgslib-include-c.h"
#include "pgslib-types.h"

void pgslib_init_access_table_c(gsTable_ptrptr)
     GS_TABLE_STRUCT **gsTable_ptrptr;
{ GS_TABLE_STRUCT   *gsTable_ptr;

  gsTable_ptr = (GS_TABLE_STRUCT *)malloc(sizeof(GS_TABLE_STRUCT));
  pgslib_check_malloc_c(gsTable_ptr, "malloc failed in pgslib_init_access_table");

  *gsTable_ptrptr = gsTable_ptr;

  initialize_item_hash_table(&(gsTable_ptr -> Hash_Table), pgslib_state.nPE);

}

void pgslib_free_access_table_c(gsTable_ptrptr)
     GS_TABLE_STRUCT **gsTable_ptrptr;
{ GS_TABLE_STRUCT   *gsTable_ptr;

  gsTable_ptr = *gsTable_ptrptr;

  free_hash_table(&(gsTable_ptr -> Hash_Table), pgslib_state.nPE);

}

void pgslib_add_item_to_table_c(item_ptr, pe_ptr, gsTable_ptrptr, ierror_ptr)
     int *item_ptr, *pe_ptr, *ierror_ptr;
     GS_TABLE_STRUCT **gsTable_ptrptr;
{
  if ( ((*gsTable_ptrptr) -> Hash_Table) == NULL)
    { pgslib_fatal_error_c("Initialize access table before adding items to it");
      return;
    }
  add_item_to_hash_table( ((*gsTable_ptrptr) -> Hash_Table), pe_ptr, item_ptr);
  *ierror_ptr = 0;
  return;
}

void pgslib_count_items_in_table_c(count_ptr, gsTable_ptrptr)
     int *count_ptr;
     GS_TABLE_STRUCT ** gsTable_ptrptr;
{
  *count_ptr = hash_count_items_in_table( ((*gsTable_ptrptr) -> Hash_Table), pgslib_state.nPE);
  return;
}
  
void pgslib_items_from_table_c(items, pes, count_ptr, gsTable_ptrptr, ierror_ptr)
     int items[], pes[];
     int *count_ptr;
     GS_TABLE_STRUCT **gsTable_ptrptr;
     int *ierror_ptr;
{
  int nitems, local_index, pe, npe;
  int N_CmplOwnerPEs, peIndex, offset, ncmpls, c;
  HASH_TABLE_ENTRY *Hash_Table;
  GS_TABLE_STRUCT  *gsTable_ptr;
  char estr[512];


  local_index = 1;
  npe = pgslib_state.nPE;

  Hash_Table = (*gsTable_ptrptr)-> Hash_Table;

  /* Setup some of the fields in the trace */
  gsTable_ptr = *gsTable_ptrptr;
  pgslib_count_items_in_table_c( &nitems, &gsTable_ptr );
  if (*count_ptr != nitems)
    pgslib_fatal_error_c("wrong count in pgslib_items_from_table");

  /* This step does
     1. Pluck items out of hash table and put into items buffer.
     2. Put PE number of each item into correspnding slot of pes buffer.
     3. Put index of item extracted (index = offset from beginning of buffer)
        into redirect slot of hash table.  This is needed later.
*/

  peIndex = 0;
  offset  = 0;
  for(pe=0;pe< pgslib_state.nPE; pe++)
    if(Hash_Table[pe].nentries > 1) {
      
      ncmpls = Hash_Table[pe].nentries - 1;

      move_items_from_hash(Hash_Table[pe], ncmpls, offset, items );
      for (c=0; c<ncmpls; c++)
	pes[offset+c] = pe;

      offset += ncmpls;
      peIndex++;
    }

  *ierror_ptr = 0;
  return;
}	

void pgslib_item_index_from_table_c(index_ptr, item_ptr, pe_ptr, gsTable_ptrptr)
     int *index_ptr, *item_ptr, *pe_ptr;
     GS_TABLE_STRUCT **gsTable_ptrptr;
{
  int index;

  index = get_redirection(((*gsTable_ptrptr)-> Hash_Table), pe_ptr, item_ptr);
  *index_ptr = index+1;

  return;
}
