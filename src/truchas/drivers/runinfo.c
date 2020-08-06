/* runinfo.C */

/*==============================================================================

  This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.

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
#include <stdio.h>

static void f90strcpy(char *ostring, char *istring)
{
  while (*istring)
    *ostring++ = *istring++;
}

void getrunhostinfo(int n, char *arch, char *host)
{
  char string[n];
  struct utsname data;

  /* fill the uname structure */
  uname(&data);

  /* build one arch string, as 'uname -a' would return */
  int sz = snprintf(string,
		    sizeof(string),
		    "%s %s %s %s %s",
		    data.sysname, data.nodename, data.release, data.version, data.machine);
  if (!sz) string[0] = 0;

  /* send to f90 */
  f90strcpy (arch, string);

  /* get the fully qualified domain name */
  int err = gethostname (string, sizeof(string));
  if (err) string[0] = 0;

  /* send to f90 */
  f90strcpy (host, string);
}

/* runinfo.C end */
