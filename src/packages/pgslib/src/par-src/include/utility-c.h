/* Header file for utility routines */

/* $Id: utility-c.h,v 1.2 2002/09/12 20:52:35 lally Exp $ */

#ifndef UTILITY_H__
#define UTILITY_H__

#define DEFAULT_IO_ROOT_PE 0

/* Global variables */
extern int PGSLib_IO_ROOT_PE;		/* Root PE used during this particular run */

/* Prototypes for utility routines in utility-c.c */
void pgslib_initialize_c(int *, int *, int *, int *, char *);
void pgslib_finalize_c();
void pgslib_error_c(char *);
void pgslib_fatal_error_c(char *);
void pgslib_abort_c();
void pgslib_output_c(char *);
void pgslib_close_output_c();
void pgslib_check_malloc_c(void *, char *);

void pgslib_get_argc(int *);
void pgslib_get_argv(int *, int *, char *);
void pgslib_cleanup_argv();

#endif
