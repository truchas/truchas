/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * error.c
 *
 *  DANU Error Handler
 *
 *
 *  Purpose:
 *           This source file defines several useful error, warning and debug
 *           printf routines that allow the package to report errors, warning 
 *           and debug messages in a standardized way.
*/

#if HAVE_CONFIG_H
#include <danu_config.h>
#endif

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

#include <danu_error.h>


/*
 * Routine: danu_fprint_mess(void)
 * Purpose: Print an ERROR, WARN or DEBUG message to stream pipe
 * Description: Print function for one of three possible modes
 *              ERROR, WARN or DEBUG. The stream pipe is
 *              flushed before exting the routine. This routine
 *              and danu_fprintf_mess are the main routines. All of the
 *              print calls in this file are wrappers to one or the other.
 *              
 * Parameters:
 *           stream		IN      Pointer to file stream
 *           mode		IN		const char pointer to string mode id (ERROR,WARN,DEBUG)
 *           message	IN		const char pointer to string with message           
 * Returns: Nothing
 * Errors: This routine should never fail
 *
 */ 
void danu_fprint_mess(FILE* stream, const char *mode, const char *message)
{
    fprintf(stream,"%s:%s\n",mode,message);

    fflush(stream);
}
/*
 * Routine: danu_fprintf_mess(void)
 * Purpose: Formatted ERROR, WARN or DEBUG message to a stream
 * Description: Print function for one of three possible modes
 *              ERROR, WARN or DEBUG. The pipe is
 *              flushed before exting the routine. This routine
 *              and danu_fprint_err_mess are the main routines. All of the
 *              print calls in this file are wrappers to one or the other.
 * Parameters:
 *           stream		IN		pointer to the file stream
 *           mode		IN		const char pointer to string mode id (ERROR,WARN,DEBUG)
 *           format		IN		const char pointer to string that formats the printf
 *           ...        IN      optional arguments. Similar to the optional arguments 
 *                              in printf.
 * Returns: Nothing
 * Errors: This routine should never fail
 *
 */ 
void danu_fprintf_mess(FILE *stream,const char * mode, const char * format, va_list ap )
{
    fprintf(stream, "%s:",mode);
    vfprintf(stream,format,ap);
    fprintf(stream,"\n");
	
    fflush(stream);

}
/*
 * Routine: danu_error(message)
 * Purpose: Print standardized error message to STDERR
 * Description: Print a formatted error message. Flush the 
 *              pipe to ensure message is not lost in a buffer
 *              if program exits.
 *              
 * Parameters:
 *           message	IN		const char pointer to message to print
 *
 * Returns: Nothing
 * Errors: This routine should never fail
 *
 */ 
void danu_error(const char * message)
{
    danu_fprint_mess(stderr,"ERROR", message);
}
/*
 * Routine: danu_warn(message)
 * Purpose: Print standardized warning message to STDOUT
 * Description: Print a formatted warning message to STDOUT. Flush the
 *              pipe to ensure message is not lost in a buffer 
 *              if program exits.
 *              
 * Parameters:
 *           message	IN		const char pointer to message to print
 *
 * Returns: Nothing
 * Errors: This routine should never fail
 *
 */ 
void danu_warn(const char *message)
{
    danu_fprint_mess(stdout,"Warning",message);
}
/*
 * Routine: danu_debug(message)
 * Purpose: Print standardized debug message to STDOUT
 * Description: Print a formatted debug message to STDOUT. Flush the
 *              pipe to ensure message is not lost in a buffer 
 *              if program exits.
 *              
 * Parameters:
 *           message	IN		const char pointer to message to print
 *
 * Returns: Nothing
 * Errors: This routine should never fail
 *
 */
void danu_debug(const char *message)
{
    danu_fprint_mess(stdout,"DEBUG", message);
}
/*
 * Routine: danu_error_printf(format,...)
 * Purpose: Print formatted error message to STDERR
 * Description: Print a formatted error message. Flush the 
 *              pipe to ensure message is not lost in a buffer
 *              if program exits.
 *              
 * Parameters:
 *           format	IN		const char pointer to a printf format string
 *           ...	IN		(optional) optional arguments
 *
 * Returns: Nothing
 * Errors: This routine should never fail
 *
 */ 
void danu_error_printf(const char * format, ... )
{
    va_list ap;

    va_start(ap,format);
    danu_fprintf_mess(stderr,"ERROR",format,ap);
    va_end(ap);
}
/*
 * Routine: danu_warn_printf(format,...)
 * Purpose: Print formatted warning message to STDERR
 * Description: Print a formatted warning message. Flush the 
 *              pipe to ensure message is not lost in a buffer
 *              if program exits.
 *              
 * Parameters:
 *           format	IN		const char pointer to a printf format string
 *           ...	IN		(optional) optional arguments
 *
 * Returns: Nothing
 * Errors: This routine should never fail
 *
 */    
void danu_warn_printf(const char * format, ... )
{
    va_list ap;

    va_start(ap,format);
    danu_fprintf_mess(stdout,"Warning",format,ap);
    va_end(ap);
}
/*
 * Routine: danu_debug_printf(format,...)
 * Purpose: Print formatted debug message to STDERR
 * Description: Print a formatted debug message. Flush the 
 *              pipe to ensure message is not lost in a buffer
 *              if program exits.
 *              
 * Parameters:
 *           format	IN		const char pointer to a printf format string
 *           ...	IN		(optional) optional arguments
 *
 * Returns: Nothing
 * Errors: This routine should never fail
 *
 */   
void danu_debug_printf(const char * format, ... )
{
    va_list ap;

    va_start(ap,format);
    danu_fprintf_mess(stdout,"DEBUG",format,ap);
    va_end(ap);
}
