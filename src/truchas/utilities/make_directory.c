/* make_directory.C */

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


void make_directory_c(char*,int*);

void make_directory_hier_c(char *path, int *status)
{
    /* make a copy of our path */
    char *p = strdup(path);
    int n = strlen(path);
    char *c = p;

    *status = 0; /* Mark as success since we might not enter make_directory_c */
    errno = 0;
    /* Try to make the entire directory tree */
    while (c) {

	c = strchr(c+1, '/');

	if (!c) {
	    make_directory_c(p, status);
	    break;
	}

	*c = 0; /* end path string at current directory separator */

        /* only call make_directory_c with a non-existing directory */
	struct stat s;
	if (stat(p, &s)) {
	    if (errno != ENOENT) {
		perror(p);
		break;
	    }
	} else {
	    *c = '/';
	    continue; /* skip call to make_direcitory_c */
	}

	/* attempt to make the named directory */
	make_directory_c(p, status);

	*c = '/';

	/* Did an error occur in creating directory */
	if (*status || n > (c+1-p)) break;
    }

    /* Free memory and return */
    free(p);
    return;

}

void make_directory_c(char *path, int *status)
{
  struct stat buf;

  /* attempt to make the named directory */
  if (!(*status = mkdir(path, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))) return;

  /* status != 0, so we've got a problem, look at errno */

  /* if it's not EEXIST (already exists) return the non-zero status */
  if (errno != EEXIST) return;

  /* the path already exists, so we need to know if it's a directory */
  /* if stat fails, we can't read something, fail */
  if (stat(path,&buf)) return;

  /* if it's a directory, return success, otherwise just return */
  if (S_ISDIR(buf.st_mode)) *status = 0;
  return;
}



/* make_directory.C end */
