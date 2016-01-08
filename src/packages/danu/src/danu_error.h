/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_error.h
*
*  DANU Error Handler
*
*/

#ifndef DANU_ERROR_H
#define DANU_ERROR_H

#include <stdio.h>
#include <stdarg.h>

#include <danu_types.h>
#include <danu_h5_error.h>

void danu_fprint_mess(FILE *stream, const char *mode, const char *mess);
void danu_fprintf_mess(FILE *stream, const char *mode, const char *format, va_list ap);

void danu_error(const char *);
void danu_warn(const char *);
void danu_debug(const char *);

void danu_error_printf(const char * format, ... );
void danu_warn_printf(const char * format, ... );
void danu_debug_printf(const char * format, ... );

/* SUCCESS and FAIL return codes */
typedef enum _danu_err_t {
    DANU_FAILURE  = -2,  /* General fail condition */
    DANU_CONTINUE = -1,  /* Loop continue condition */
    DANU_SUCCESS  =  0,  /* Successfull return status */
    DANU_STOP     =  1,  /* Loop stop condition */
} danu_err_t;


#define DANU_RETURN_SUCCESS(a)       (  ( (a) == DANU_SUCCESS )  ? TRUE : FALSE )
#define DANU_RETURN_FAIL(a)          (  ( (a) == DANU_FAILURE    )  ? TRUE : FALSE )
#define DANU_RETURN_OK(a)          (  ( (a) >= 0    )  ? TRUE : FALSE )

#define DANU_RETURN_SELECT(a,b)  ( (a) != DANU_SUCCESS ? (a) : ( (b) != DANU_SUCCESS ? (b) : DANU_SUCCESS ))

#define DANU_BAD_PTR(a)          ( (a) == NULL ) ? TRUE  : FALSE
#define DANU_EMPTY_STRING(a)     ( *(a) == '\0' ) ? TRUE : FALSE
#define DANU_BAD_STRING(a)       ( ( (a) == NULL ) || ( *(a) == '\0' ) ) ? TRUE : FALSE




/* Useful macros to print error and debug messages */

/* DANU_ERROR_MESS: Print formatted error message the includes line number
 *                  and function name
 *
 */

#define DANU_ERROR_MESS(message) \
			( fprintf(stderr,"ERROR (%s %s,line=%d):%s\n", __func__, __FILE__,__LINE__, (message)) );

/* DANU_WARN_MESS: Print formatted warning message the includes line number
 *                  and function name
 *
 */

#define DANU_WARN_MESS(message) \
			( fprintf(stdout,"Warning (%s %s,line=%d): %s\n", __func__, __FILE__, __LINE__, (message)) );

/* DANU_DEBUG_MESS: Print formatted debug message the includes line number
 *                  and function name
 *
 */

#define DANU_DEBUG_MESS(message) \
			( fprintf(stdout,"Debug (%s %s,line=%d): %s\n",__func__, __FILE__, __LINE__, (message)) ); 

#endif
