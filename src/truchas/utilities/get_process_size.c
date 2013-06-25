//
// Get various memory sizes for the process.
//

#ifdef LINUX

#include <unistd.h>
#include <stdio.h>

void get_process_size (int *vsize, int *rsize, int *dsize)
{
  // See 'man proc' for a description of the fields in the statm proc file.
  // Other possibilities are the stat and status proc files.
  char fname[128];
  FILE *fp;
  sprintf(fname,"/proc/%d/statm",getpid());
  if ((fp = fopen(fname,"r")) == NULL) {
    *vsize = 0;
    *rsize = 0;
    *dsize = 0;
  } else {
    fscanf(fp, "%d %d %*s %*s %*s %d", vsize, rsize, dsize);
    fclose(fp);
    int psize = getpagesize();
    *vsize = *vsize * psize / 1024;
    *rsize = *rsize * psize / 1024;
    *dsize = *dsize * psize / 1024;
  }
}

#else

void get_process_size (int *vsize, int *rsize, int *dsize)
{
  *vsize = 0;
  *rsize = 0;
  *dsize = 0;
}

#endif
