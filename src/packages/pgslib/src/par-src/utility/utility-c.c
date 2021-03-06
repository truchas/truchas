/* These subroutines support the utilities in PGSLib.
   They are called by routines in the PGLSib_Utility_Module.
*/

/* $Id: utility-c.c,v 1.2 2002/09/12 20:52:35 lally Exp $ */

#include <stdlib.h>
#include <stdio.h>
#ifdef GNU_COMPILER
#include <string.h>
#else
#include <strings.h>
#endif
#include "pgslib-include-c.h"

#define USE_FILE_OUTPUT

#include "utility-c.h"

#ifdef USE_SGI_SHMEM_LIB
#include "sm.h"
int *sm_space;
long allocated;
#endif

#ifdef NAG_COMPILER
#  define _IARGC_   nag_iargc_
#  define _GETARG_  nag_getarg_
#elif GNU_COMPILER
#  define _IARGC_   _gfortran_iargc
#  define _GETARG_  _gfortran_getarg_i4
#else
#  define _IARGC_   iargc_
#  define _GETARG_  getarg_
#endif

/* Hack stuff to get some timing info */

#include "timing-c.h"
float barrier_time = 0.0;
float sr_time = 0.0;


/* General information struct about state of the system */
parallel_state pgslib_state;		/* Made available through pgslib-include-c.h */

int		PGSLib_IO_ROOT_PE;	/* Made available through utility-c.h	       */

static int pgslib_fperr_opened=0;
static FILE *fperr;
static char ferrname[2048] = "";
static int pgslib_fpout_opened=0;
static FILE *fpout;
static char foutname[2048] = "";

static int been_initialized = FALSE;
static int pgslib_doesnt_init_mpi;

/* supports mpi command line mess */
/* this must agree with what's on the F side - PGSLIB_CL_MAX_TOKEN_LENGTH */
#define MAX_TOKEN_LENGTH 1023

/* create argc/argv */
static int argc;
static char **argv;

/* We're also going to save pointers to the originals, as we allocate
   them and should deallocate them later.  MPI_Init is free to change
   them, so we can't rely on what comes back from that call */
static int argc_save;
static char **argv_save;
static char **argv_value_save;

/* end of mpi command line mess */


void pgslib_initialize_c(nPE, thisPE, IO_ROOT_PE, File_Per_PE, File_Prefix)
     int *nPE, *thisPE, *IO_ROOT_PE, *File_Per_PE;
     char *File_Prefix;
{ /* int ierror; */
  int smreturn;
  char ProgName[] = "PGSLib_MPI";
  char ErrorString[] = "ERROR in initialize.";

  if (! been_initialized)
    {	been_initialized = TRUE;

	MPI_Initialized(&pgslib_doesnt_init_mpi);
	if (! pgslib_doesnt_init_mpi) {
	  /* we now lift the restriction of calling the fortran mpi init */
	  /* Need to call the fortran version of init, since main is F90 code. */
	  /* pgslib_mpi_init(&ierror); */

	  char *p;
	  int i;

	  /* get the arg count, and save the original value */
	  argc_save = argc = _IARGC_() + 1;
  
	  /* allocate a vector of character pointers argc long, and save the original pointer */
	  argv_save = argv = (char **)malloc(argc*sizeof(char *));

	  /* allocate a vector of character pointers argc long to save the original pointers to the command line tokens in */
	  argv_value_save = (char **)malloc(argc*sizeof(char *));

	  for (i=0; i<argc; i++) {
	    argv_value_save[i] = argv[i] = (char *)malloc(MAX_TOKEN_LENGTH+1);
	    _GETARG_(&i, argv[i], MAX_TOKEN_LENGTH);

	    /* trim trailing blanks (arg came from fortran rtl */
	    p = argv[i] + MAX_TOKEN_LENGTH - 1;
	    while (p >= argv[i]) {
	      if (*p != ' ' && *p) {
		*(p+1) = '\0';
		break;
	      }
	      p--;
	    }
	  }

#define PASS_NULL_TO_MPI_INIT
#ifdef PASS_NULL_TO_MPI_INIT
          /* MPI-1.2 permits applications to pass NULL for &argc and &argv; */
          /* this is a workaround for an apparent bug in lam-7.1.1 mpi_init */
	  if (MPI_Init(NULL, NULL) != MPI_SUCCESS)
#else
	  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
#endif
	    fprintf(stderr, "ERROR: Could not initialize MPI");

	}     
	
	/*Rather than use MPI_COMM_WORLD everywhere, start using
	  my own communicator.  For the moment (7/27/00) most of PGSLib
	  assumes that pgslib_state.PGSLib_Comm is the same as MPI_COMM_WORLD.  Eventually that
	  assumption will vanish.*/

	MPI_Comm_dup(MPI_COMM_WORLD, &(pgslib_state.PGSLib_Comm));
	MPI_Comm_size( pgslib_state.PGSLib_Comm, nPE);
	MPI_Comm_rank( pgslib_state.PGSLib_Comm, thisPE);

#ifdef USE_SGI_SHMEM_LIB
        smreturn = sm_init( 0, 0, pgslib_state.PGSLib_Comm);
#endif	

	/* Establish file names */
	if (*File_Per_PE == PGSLIB_TRUE)
	  {
	    sprintf(foutname, "%s-out.%04d", File_Prefix, *thisPE);
	    sprintf(ferrname, "%s-err.%04d", File_Prefix, *thisPE);
	  }
	else
	  {
	    sprintf(foutname, "%s-out.ALL", File_Prefix);
	    sprintf(ferrname, "%s-err.ALL", File_Prefix);
	  }

	/* Use user request for IO_ROOT_PE, unless it is out of range */
	/* User request is 0 based */
	PGSLib_IO_ROOT_PE = *IO_ROOT_PE;
	if( *IO_ROOT_PE < 0) 
	  PGSLib_IO_ROOT_PE = DEFAULT_IO_ROOT_PE;
	if ( *IO_ROOT_PE > *nPE) 
	  PGSLib_IO_ROOT_PE = DEFAULT_IO_ROOT_PE;

	pgslib_state.initialized = TRUE;

	pgslib_state.nPE = *nPE;
	pgslib_state.thisPE = *thisPE;
	pgslib_state.io_pe = PGSLib_IO_ROOT_PE;
      }

  /* Need to initialize the tags */
  initializeTagBits();

  *nPE = pgslib_state.nPE;
  *thisPE = pgslib_state.thisPE;
  *IO_ROOT_PE = pgslib_state.io_pe;

}

void pgslib_finalize_c()
{ int ierror;
  pgslib_close_output_c();	
  been_initialized = FALSE;
  if (! pgslib_doesnt_init_mpi) {
#ifdef USE_SGI_SHMEM_LIB
/*	  sm_exit();*/
#endif	
    MPI_Finalize();	
  }
}


void pgslib_error_c(ErrorString)
     char *ErrorString;
{ int myrank;

  MPI_Comm_rank( pgslib_state.PGSLib_Comm, &myrank );

#ifdef USE_FILE_OUTPUT
  if( !pgslib_fperr_opened) {
    fperr = fopen(ferrname, "w");
    pgslib_fperr_opened = TRUE;
  }
  fprintf(fperr, " %3d: ERROR %s\n", myrank, ErrorString);
#else
  fprintf(stderr, " %3d: ERROR %s\n", myrank, ErrorString);
#endif
}


void pgslib_fatal_error_c(ErrorString)
     char *ErrorString;
{ int myrank;
  pgslib_error_c(ErrorString);
  pgslib_abort_c();
}

void pgslib_abort_c()
{
  pgslib_close_output_c();
  MPI_Abort( pgslib_state.PGSLib_Comm, 1);
}

void pgslib_output_c(Message)
     char *Message;
{ int myrank;

  MPI_Comm_rank( pgslib_state.PGSLib_Comm, &myrank );

#ifdef USE_FILE_OUTPUT
  if( !pgslib_fpout_opened) {
    fpout = fopen(foutname, "w");
    if( fpout == NULL)
      fprintf(stderr, "fopen failed in pgslib_output");
    pgslib_fpout_opened = TRUE;
  }
  fprintf(fpout, " %3d: %s\n", myrank, Message);
#else
  fprintf(stdout, " %3d: %s\n", myrank, Message);
#endif
}

void pgslib_flush_output_c()
{
#ifdef USE_FILE_OUTPUT
  if ( fpout != NULL)
  {
    fflush(fpout);
  }
  if ( fperr != NULL)
  {
    fflush(fperr);
  }
#else
  fflush(stdout);
  fflush(stderr);
#endif
}

void pgslib_close_output_c()
{
#ifdef USE_FILE_OUTPUT
  if ( fpout != NULL)
  {
    fclose(fpout);
    pgslib_fpout_opened = FALSE;
  }
  if ( fperr != NULL)
  {
    fclose(fperr);
    pgslib_fperr_opened = FALSE;
  }
#endif
}

void pgslib_check_malloc_c(Pointer, ErrorString)
     void *Pointer;
     char *ErrorString;
{
  if(Pointer == NULL) {
    pgslib_error_c(ErrorString);
  }
}

void pgslib_barrier_c()
{
  MPI_Barrier( pgslib_state.PGSLib_Comm );
}

void pgslib_barrier_time_c(float *bt)
{
  *bt = barrier_time;
  return;
}

void pgslib_sr_time_c(float *t)
{
  *t = sr_time;
  return;
}

/* return argc to F */
void pgslib_get_argc (int *n)
{
  *n = argc;
}

/* return argv strings to F one at a time, as requested */
void pgslib_get_argv(int *i, int *string_l, char *string)
{
  char *p;

  /* return the i'th arg */
  strncpy(string, argv[*i], *string_l);

  /* pad with spaces, as it's going back to the fortran rtl */
  for (p=string+strlen(string); p<string+(*string_l); p++)
    *p = ' ';
}

/* deallocate what we allocated */
void pgslib_cleanup_argv()
{
  int i;

  /* free what we allocated */
  for (i=0; i<argc_save; i++) {
    free(argv_value_save[i]);
  }
  free(argv_save);
  free(argv_value_save);
}
