/*  dlwrap.c
 *
 *  Fortran-callable wrapper routines around the dynamic linking loader
 *  procedures dlopen, dlclose, dlsym, dlerror from libdl.  This solves
 *  several issues:
 *    o pass-by-value of integer flag to dlopen.
 *    o passing void pointer; we pass a pointer to an address-storing integer.
 *    o Fortran's hidden character length argument.
 *
 *  10 Apr 2006, Neil Carlson <nnc@lanl.gov>: Initial version
 *  17 Apr 2006, Neil Carlson <nnc@lanl.gov>: 
 *
 */

/* typical Fortran name mangling; your mileage may vary */
#include <FortranCInterface_names.h>

#define f_dlopen      TR_ROUTINE_GLOBAL_(f_dlopen,F_DLOPEN)
#define f_dlclose     TR_ROUTINE_GLOBAL_(f_dlclose,F_DLCLOSE)
#define f_dlsym       TR_ROUTINE_GLOBAL_(f_dlsym,F_DLSYM)
#define f_dlerror     TR_ROUTINE_GLOBAL_(f_dlerror,F_DLERROR)
#define c_string_to_f TR_ROUTINE_GLOBAL_(c_string_to_f,C_STRING_TO_F)

#include <dlfcn.h>
#include <string.h>

void *f_dlopen(const char *filename, int *flag, int len_filename)
{
  return dlopen(filename, *flag);
}

void *f_dlsym(void **handle, const char *symbol, int len_symbol)
{
  return dlsym(*handle, symbol);
}

int f_dlclose(void **handle)
{
  return dlclose(*handle);
}

char *f_dlerror(void)
{
  return dlerror();
}

/* Fortran can't directly dereference C pointers.  The following function
   takes a (pointer to a) string pointer returned by f_dlerror and copies
   the contents to a Fortran character variable. */
   
void c_string_to_f (char **s, char *f, int f_len)
{
  int i;
  for(i=0; i<f_len; i++) {
    if (i<strlen(*s))
      f[i] = (*s)[i];
    else
      f[i] = ' ';
  }
}
    
