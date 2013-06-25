/* Header file for all files in gath-scatt */

/* $Id: gs-c.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

/* Prototypes for routines in gs-util-c.c */

void pgslib_gs_init_trace_c(GS_TRACE_STRUCT **);
void pgslib_setup_n_duplicate_c(int *, GS_TRACE_STRUCT **);
void pgslib_gs_release_trace_c(GS_TRACE_STRUCT **);
void gs_free_send_buffer(COMM_SEND *);
void gs_free_rcv_buffer(COMM_RCV *);
