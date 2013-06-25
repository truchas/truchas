/* Header file for routines in the indexing directory */

/* $Id: index-c.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */
/* Prototypes for routines in index-c.c */

void pgslib_init_access_table_c(GS_TRACE_STRUCT **);
void pgslib_add_item_to_table_c(int *, int *, GS_TRACE_STRUCT **, int *);
void pgslib_count_items_in_table_c(int *, GS_TRACE_STRUCT **);
void pgslib_items_from_table_c(int[], int[], int *, GS_TRACE_STRUCT **, int *);
void pgslib_item_index_from_table_c(int *, int *, int *, GS_TRACE_STRUCT **);
