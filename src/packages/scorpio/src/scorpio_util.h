/*******************************************************************************
* Copyright Notice
*  + 2010-2012 North Carolina State University
*  + 2010-2012 Pacific Northwest National Laboratory
* 
* This file is part of SCORPIO.
* 
* SCORPIO is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* any later version.
* 
* SCORPIO is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public License
* along with SCORPIO.  If not, see <http://www.gnu.org/licenses/>.
* 
*******************************************************************************/
#ifndef _SCORPIO_UTIL_H
#define _SCORPIO_UTIL_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdarg.h>
#include <string.h>

// Error is 1, Warning is 2 and so on.
typedef enum 
{
	SCORPIO_ERROR = 1,
	SCORPIO_WARNING,
	SCORPIO_INFO,
	SCORPIO_VERBOSE
} scorpio_loglevel_t ;

/** Length of string used to store file and function names. */
#define MAXSTR 512

typedef struct
{
	char function[MAXSTR];
	char file[MAXSTR];
	int line;
	int rank;
} scorpio_debuginfo_t ;

/** 
 * @brief Initialize debugging log
 * @param rank global rank that is used to print debug messages.
 */
void scorpio_debuglog_init(int rank);

/** 
 * @brief Helper function to log messages to stdout
 * 
 * @param scorpio_loglevel scorpio_loglevel applicable to this message
 * @param format format string similar to printf format
 * @param ... variable argument list that is passed on to printf
 */
void scorpio_printmsg(scorpio_loglevel_t scorpio_loglevel, const char *format, ...);

// Stores current debugging level to print messages. 
extern int scorpio_debug;

// Stores current function name, file and line number for priting debugging messages
extern scorpio_debuginfo_t scorpio_debuginfo;

// gcc.gnu.org docs:  __FUNCTION__ is another name for __func__. Older versions of GCC recognize only this name. However, it is not standardized. For maximum portability, we recommend you use __func__, but provide a fallback definition with the preprocessor: 	
#if __STDC_VERSION__ < 199901L
	#if __GNUC__ >= 2
		#define __func__ __FUNCTION__
	#else
		#define __func__ "<unknown>"
	#endif
#endif

/* Macro to help log messages to stdout. Only defined when SCORPIO_DEBUG is defined at compile time.
 * NOTE: Use (( )) to enclose the log message.
 * Example call to log message: PRINT_MSG((SCORPIO_INFO," Sample message: %d %f %s ", 10, 234.56, "sarat"));  
 * IMPLEMENTATION NOTE: Equivalent solution can be obtained using variadic macros (C99 standard) though it may not be supported on all platforms.
 * Whereas this solution works on all platforms.
 * __FILE__ , __LINE__ for current file and line number, __func__ for the current function (C99) which is equivalent to __FUNCTION__ in earlier versions. 
 */ 
#if defined(SCORPIO_DEBUG)
	#define PRINT_MSG(x) strcpy(scorpio_debuginfo.function, (char *) __func__); strcpy(scorpio_debuginfo.file, (char *)  __FILE__); scorpio_debuginfo.line = __LINE__; scorpio_printmsg x 
#else
	#define PRINT_MSG(x)
#endif

#ifdef __cplusplus
}
#endif

#endif 
