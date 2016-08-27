/* runinfo.C */

/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*------------------------------------------------------------------------------
 * get system information at runtime
 *
 * too hard (or impossible) to do in F90
 *
 * author: Bryan Lally, lally@lanl.gov
 *----------------------------------------------------------------------------*/

#include <sys/utsname.h>
#include <unistd.h>
#include <string.h>

static void f90strcpy(char *ostring, char *istring)
{
  while (*istring)
    *ostring++ = *istring++;
}

void getrunhostinfo(char *arch, char *host)
{
  char string[128];
  struct utsname data;

  /* fill the uname structure */
  uname(&data);

  /* build one arch string, as 'uname -a' would return */
  strcpy (string, data.sysname);
  strcat (string, " ");
  strcat (string, data.nodename);
  strcat (string, " ");
  strcat (string, data.release);
  strcat (string, " ");
  strcat (string, data.version);
  strcat (string, " ");
  strcat (string, data.machine);

  /* send to f90 */
  f90strcpy (arch, string);

  /* get the fully qualified domain name */
  gethostname (string, sizeof(string));

  /* send to f90 */
  f90strcpy (host, string);
}

/* runinfo.C end */
