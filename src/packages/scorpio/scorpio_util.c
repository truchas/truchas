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
#include "scorpio_util.h"

// Set debugging level. 
#if defined(SCORPIO_DEBUG)
	// Stores current debugging level to print messages.
	int scorpio_debug = SCORPIO_DEBUG;
#else
	int scorpio_debug = 0;
#endif

scorpio_debuginfo_t scorpio_debuginfo; 

void scorpio_debuglog_init(int rank)
{
	// Set stderr to line buffering mode. Printing out once a '\n' is encountered. This is useful to separate log messages from different processes (parallel)
	setlinebuf(stderr);

	scorpio_debuginfo.rank = rank;
}

/** 
 * @brief Helper function to log messages to stdout
 * 
 * @param scorpio_loglevel scorpio_loglevel applicable to this message
 * @param format format string similar to printf format
 * @param ... variable argument list that is passed on to printf
 */
void scorpio_printmsg(scorpio_loglevel_t scorpio_loglevel, const char *format, ...)
{
	va_list varargs;

	// Sarat: if -DSCORPIO_DEBUG is given in CFLAGS, debug is initialized to 1. If -DSCORPIO_DEBUG=NUM is used, then debug=NUM
	/* printf("\n scorpio_debug level: %d, scorpio_loglevel: %d", scorpio_debug, scorpio_loglevel); */

	// If current message's level is greater (more verbose) than debug level (defined at compile-time) that message is ignored. 
	if (scorpio_loglevel > scorpio_debug)
	{ 
		/* fprintf(stderr, " Ignoring log message. \n"); */
		return;
	}

	fprintf(stderr, "\n %d:LOG:", scorpio_debuginfo.rank);
	
	switch(scorpio_loglevel)
	{
		case SCORPIO_VERBOSE: 
			fprintf(stderr, "SCORPIO_VERBOSE: ");
			break;
		case SCORPIO_INFO: 
			fprintf(stderr, "SCORPIO_INFO: ");
			break;
		case SCORPIO_WARNING:
			fprintf(stderr, "SCORPIO_WARNING: ");
			break;
		case SCORPIO_ERROR: 
			fprintf(stderr, "SCORPIO_ERROR: ");
			break;
	}	

	fprintf(stderr, "%s (%s:%d)--> ", scorpio_debuginfo.function, scorpio_debuginfo.file, scorpio_debuginfo.line );

	/* Start variable parameter list after format */
	va_start(varargs, format);
	vfprintf(stderr, format, varargs);
	va_end(varargs);
	fprintf(stderr, "\n");

	fflush(stderr);

	/* int nvarArgs = va_arg(varargs,int); */
	/* printf("\n %d\n", nvarArgs); */
}
