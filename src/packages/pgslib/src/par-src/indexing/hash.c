/* Hash table creation and maintenance functions */

/* $Id: hash.c,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include "pgslib-include-c.h" 
#include "pgslib-types.h"

static itemlink *get_new_slot(HASH_TABLE_ENTRY *);
static itemlink *get_item_from_hash_table(HASH_TABLE_ENTRY *, int);

void add_item_to_hash_table(Hash_Table, PEptr, itemptr)
     HASH_TABLE_ENTRY *Hash_Table;
     int              *PEptr, *itemptr;
{ 
    itemlink *thisslot, *nextslot, *newslot;
    HASH_TABLE_ENTRY *hash;
    int item;

/* Find the hash table entry to work with. */
    hash = & Hash_Table[*PEptr];
    item = *itemptr;
/* First find out if this item is already in the hash table entry.
   If it is, then we are done.  */

  thisslot=get_item_from_hash_table(hash, item);

/* If this item is not in the table, then add it */
/* thisslot is set from above. */
  if (item != thisslot -> items) {
    newslot = get_new_slot(hash);
    nextslot = thisslot -> links;
    thisslot -> links = newslot;
    newslot -> links = nextslot;
    newslot -> items = item;
  }

  return;
}

int get_redirection(Hash_Table, PEptr, itemptr)
     HASH_TABLE_ENTRY *Hash_Table;
     int              *PEptr, *itemptr;
{ 	
  HASH_TABLE_ENTRY *hash;
  int item, PEIndex;
  itemlink *thisslot;
  
  hash = &Hash_Table[*PEptr];
  item = *itemptr;
  thisslot = get_item_from_hash_table(hash, item);
  return(thisslot->redirect);
} 


/* Get_item_from_hash_table returns the itemlink in the hash
   table which satisfies:
   1. input_item > next_item_in_hash_table 
   2. if hash_table has more than a single item then
      returned_item >= input_item.
   Therefore, if input_item is in the hash table then it is 
   returned by get_item_from_hash_table.  
   */
itemlink *get_item_from_hash_table(hash, item)
     HASH_TABLE_ENTRY *hash;
     int              item;
{ itemlink *thisslot, *nextslot, *newslot;

  thisslot = hash->firstitemlinkchunk->itemlinkbuffer;
  nextslot = thisslot -> links;
  while (item <= nextslot -> items) {
    thisslot = nextslot;
    nextslot = nextslot -> links;
  }

  return(thisslot);
}


itemlink *get_new_slot(hash)
     HASH_TABLE_ENTRY *hash;
{ itemlinkchunk *thischunk, *newchunk;
  if (hash->nentries < hash->totalslots) {
    hash->nentries += 1;
    hash->newslot++;
  }
  else {
    /* Find that last chunk, so we can link to it */
    thischunk = hash->firstitemlinkchunk; 
    while (thischunk -> nextitemlinkchunk != NULL)
      thischunk = thischunk->nextitemlinkchunk;

    /* Malloc a new chunk */
    thischunk->nextitemlinkchunk = (itemlinkchunk *)malloc(sizeof(itemlinkchunk));
    pgslib_check_malloc_c(thischunk -> nextitemlinkchunk,
		      "malloc failed in get_new_slot");

    newchunk = thischunk -> nextitemlinkchunk;
    newchunk -> size = HASH_EXPANSION_RATE*thischunk->size;
    newchunk -> nextitemlinkchunk = NULL;
    newchunk->itemlinkbuffer = (itemlink *)malloc(newchunk->size * sizeof(itemlink));
    pgslib_check_malloc_c(newchunk -> itemlinkbuffer,
		      "malloc failed in get_new_slot");


    /* Now we have a new chunk of itemlinks. */
    hash->totalslots += newchunk->size;
    hash->nentries += 1;
    hash->newslot = newchunk->itemlinkbuffer;
  }

  return(hash->newslot);
}

void initialize_item_hash_table(Hash_Table_ptr,nproc)
     HASH_TABLE_ENTRY **Hash_Table_ptr;
     int nproc;
{ int p;
  itemlinkchunk *firstchunk;
/* Allocate space for the table.  */
  *Hash_Table_ptr = (HASH_TABLE_ENTRY *)malloc(nproc*sizeof(HASH_TABLE_ENTRY));

/* Initialize each entry */
  for (p=0;p<nproc;p++) {
    (*Hash_Table_ptr)[p].totalslots = HASH_INITIAL_SIZE;
    (*Hash_Table_ptr)[p].nentries = 1; /* first slot is used for list terminator */
 
/* Initialize first chunk for each table entry */
    (*Hash_Table_ptr)[p].firstitemlinkchunk = (itemlinkchunk *)malloc(sizeof(itemlinkchunk));
    pgslib_check_malloc_c((*Hash_Table_ptr)[p].firstitemlinkchunk,
		      "malloc failed in initialize_item_hash_table");
    firstchunk = (*Hash_Table_ptr)[p].firstitemlinkchunk;
    firstchunk->size = HASH_INITIAL_SIZE;
    firstchunk->nextitemlinkchunk = NULL;

/* Initialize first itemlinkbuffer for each chunk */
    firstchunk->itemlinkbuffer = 
      (itemlink *)malloc(HASH_INITIAL_SIZE*sizeof(itemlink));
    pgslib_check_malloc_c(firstchunk->itemlinkbuffer,
		      "malloc failed in initialize_item_hash_table");

    firstchunk->itemlinkbuffer->links = firstchunk->itemlinkbuffer;
    firstchunk->itemlinkbuffer->items = HASH_MINIMUM_ITEM;

    (*Hash_Table_ptr)[p].newslot = firstchunk->itemlinkbuffer->links;

  }
  return;
}

void free_hash_table(Hash_Table_ptr, nproc)
     HASH_TABLE_ENTRY **Hash_Table_ptr;
     int nproc;
{
  int p;
  HASH_TABLE_ENTRY *Hash_Table;	
  itemlinkchunk *thischunk, *nextchunk;

  Hash_Table = *Hash_Table_ptr;
/* First loop throug all entries, freeing each chunk and buffer.  Then free table */
  for(p=0;p<nproc;p++) {
    nextchunk = Hash_Table[p].firstitemlinkchunk;
    while(nextchunk != NULL) {
      thischunk = nextchunk;
      nextchunk = nextchunk->nextitemlinkchunk;
      free(thischunk->itemlinkbuffer);
      free(thischunk);
    }
  }

/* Now free the table itself. */
  free(Hash_Table);
  *Hash_Table_ptr = NULL;
}

void move_items_from_hash(hash, itembuflen, itemoffset, itembuffer)
    HASH_TABLE_ENTRY hash;
    int              itembuflen, itemoffset, *itembuffer;

/* This routine moves item global numbers from the hash table entry
   into the itembuffer array.  It also puts a rediretion pointer into
   the redirect slot of the hash table.*/

{ int *items;
  int itemn, itemsthischunk, bufferoffset;
  itemlinkchunk *currentchunk;
  itemlink      *currentitemlinkbuffer;
  itemlink      *thisslot;

/* Check that the itembuffer has the right number of entries. */
  if(hash.nentries-1 != itembuflen) {
      pgslib_error_c("Wrong number of entries for itembuflen in move_items_from_hash.");
      return;
  }

  /* Start filling from first position */
  itemn = 0; 
  thisslot = hash.firstitemlinkchunk->itemlinkbuffer -> links;
  /* Extract items with largest item first.  Walk pointers until end */
  while (thisslot -> items > HASH_MINIMUM_ITEM) {
    itembuffer[itemoffset+itemn] =  thisslot -> items;
    itemn++;
    thisslot = thisslot -> links;
  }
#ifdef USE_OLD_HASH_CODE
#error "You Shouldn't be using this code!"
/*Get pointer to the first buffer of items.*/
  currentchunk = hash.firstitemlinkchunk;
  itemsthischunk = currentchunk->size; 
  bufferoffset = 1; /* First one is a dummy entry */
  currentitemlinkbuffer = currentchunk->itemlinkbuffer;
  itemn = 0;
  while(itemn<itembuflen) {
      if(bufferoffset>=itemsthischunk) {  /* get a new chunk */
	  currentchunk = currentchunk->nextitemlinkchunk;
	  currentitemlinkbuffer = currentchunk->itemlinkbuffer;
	  itemsthischunk = currentchunk->size;
	  bufferoffset = 0;
      }
      itembuffer[itemoffset] = (currentitemlinkbuffer+bufferoffset)->items;
      (currentitemlinkbuffer+bufferoffset)->redirect = itemoffset;
      itemn++;
      itemoffset++;
      bufferoffset++;
  }
#endif
  return;
}

int hash_count_items_in_table(Hash_Table, nPE)
     HASH_TABLE_ENTRY *Hash_Table;
     int nPE;
{
  int pe, count;
  
  count = 0;
  for (pe=0;pe<nPE; pe++)
    count += Hash_Table[pe].nentries - 1;

  return(count);
}



