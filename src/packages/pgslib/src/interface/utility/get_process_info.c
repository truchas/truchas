/* get_process_info.c */

/* provide procedures to get the pid and the virtual size of the current process */

/* $Id: get_process_info.c,v 1.1.1.1 2000/10/11 22:44:31 ferrell Exp $ */

#include <unistd.h>

/**********************************************************************/
/* Translator Macro                                                   */
/**********************************************************************/

/* FortranCInterface_names.h is created by CMake. It contians
   macros that manage the Fortran to C name mangling. The
   macro TR_ROUTINE_GLOBAL_ handles global routines with underscores
   in the name */
#include <FortranCInterface_names.h>


#define pgslib_get_process_id_c   TR_ROUTINE_GLOBAL_(pgslib_get_process_id_c,PGSLIB_GET_PROCESS_ID_C)
#define pgslib_get_vm_size_c      TR_ROUTINE_GLOBAL_(pgslib_get_vm_size_c,PGSLIB_GET_VM_SIZE_C)

/* get process id */
void pgslib_get_process_id_c (int *pid)
{
  *pid = getpid();
}

#ifdef LINUX
#ifndef GET_VM_SIZE_DEFINED
#define GET_VM_SIZE_DEFINED 1
void pgslib_get_vm_size_c(int *vmsize)
{
  /* Linux version not implemented yet */
  *vmsize = -1;
  return;
}
#endif
#endif /* LINUX */

#ifdef IRIX
#ifndef GET_VM_SIZE_DEFINED
#define GET_VM_SIZE_DEFINED 1

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/procfs.h>
#include <stdio.h>
#include <stropts.h>

/* get virtual size - works for IRIX, might work for other /proc filesystems */
void pgslib_get_vm_size_c (int *vmsize)
{
  char fname[128];
  int fd;
  struct prpsinfo info;

  sprintf(fname,"/proc/pinfo/%010d",getpid());
  if ((fd = open(fname,O_RDONLY)) == -1) {
    *vmsize = -1;		/* return -1 if we can't get it this way */
  } else {
    ioctl(fd, PIOCPSINFO, &info);
    close(fd);
    *vmsize = info.pr_size * getpagesize() / 1024;  /* return vsize in kb */
  }
}

#endif
#endif /* IRIX */

/* If we didn't pick up any of those, then we definitely don't have 
   a useful function */
#ifndef GET_VM_SIZE_DEFINED    /* Generic, undefined */
#define GET_VM_SIZE_DEFINED 1
void pgslib_get_vm_size_c(int *vmsize)
{
  /* not implemented yet */
  *vmsize = -1;
  return;
}
#endif  /* Generic, undefined */
