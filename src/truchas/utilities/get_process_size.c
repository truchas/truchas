//
// Get various memory sizes for the process.
//

#ifdef LINUX

#include <unistd.h>
#include <stdio.h>
#include <stdint.h>

void get_process_size (int64_t *vsize, int64_t *rsize, int64_t *dsize)
{
  // See 'man proc' for a description of the fields in the statm proc file.
  // Other possibilities are the stat and status proc files.
  char fname[128];
  FILE *fp;
  long pagesize;
  long long vpages, rpages, dpages;
  sprintf(fname,"/proc/%d/statm",getpid());
  if ((fp = fopen(fname,"r")) == NULL) {
    *vsize = 0;
    *rsize = 0;
    *dsize = 0;
  } else {
    if (fscanf(fp, "%lld %lld %*s %*s %*s %lld", &vpages, &rpages, &dpages) != 3) {
      vpages = 0;
      rpages = 0;
      dpages = 0;
    }
    fclose(fp);
    pagesize = sysconf(_SC_PAGESIZE);
    if (pagesize <= 0) pagesize = getpagesize();
    *vsize = vpages * pagesize / 1024;
    *rsize = rpages * pagesize / 1024;
    *dsize = dpages * pagesize / 1024;
  }
}

#else

#include <stdint.h>

void get_process_size (int64_t *vsize, int64_t *rsize, int64_t *dsize)
{
  *vsize = 0;
  *rsize = 0;
  *dsize = 0;
}

#endif
