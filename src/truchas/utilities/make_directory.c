/* make_directory.C */

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

void make_directory_c(char*,int*);

void make_directory_hier_c(char *path, int *status)
{
  struct stat buf;
  char *tmppath;
  char *c;
  /* make a copy of our path */
  tmppath = (char *) malloc((10+strlen(path))*sizeof(char));

  /* Try to make the entire directory tree */

  c = path;
  while (c) {

    c = strchr(c, '/');

    (void) strcpy(tmppath, path);
    if ( c ) { 
      tmppath[(int)(c-path)+1] = '\0';
      c++; /* advance to next character */
    }

    /* attempt to make the named directory */
    (void) make_directory_c(tmppath, status);

    /* Did an error occur in creating directory */
    if ( *status) break; 
  }

  /* Free memory and return */
  free(tmppath);
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
