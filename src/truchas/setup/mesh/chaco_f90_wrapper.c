/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/* this routine is an interface to chaco
   this facilitates calling Chaco (C code) from F90 */
/* Robert Ferrell, ferrell@cpca, Oct 24, 1996 */

/* If USE_CHACO is defined, then this code calls the "interface"routine 
   to Chaco and returns coloring in assignment_int.
   If USE_CHACO is not defined, then Chaco is not linked in and
   assignment_int is set to 0.
*/

#include <stdlib.h>
#include <stdio.h>

#include <FortranCInterface_names.h>

#define chaco_f90_wrapper  TR_ROUTINE_GLOBAL_(chaco_f90_wrapper,CHACO_F90_WRAPPER)
#define chaco_f90_wrapper2 TR_ROUTINE_GLOBAL_(chaco_f90_wrapper2,CHACO_F90_WRAPPER2)

void chaco_f90_wrapper(nvtxs_ptr, nPEs_ptr, start, adjacency, assignment_int, status)
int	*nvtxs_ptr;      /* number of vertices in full graph */
int	*start;          /* Array of start indices for edge segments */
int     *adjacency;      /* edge list data */
int     *assignment_int; /* set number for each vertex */
int     *nPEs_ptr;       /* Number of processors to partition for */
int     *status;         /* Return status */

{
  int	nvtxs, vtx;
  int mesh_dims[3];
  short *assignment;
  double eigtol;
  mesh_dims[0] = *nPEs_ptr;
  mesh_dims[1] = 0;
  mesh_dims[2] = 0;

  nvtxs = *nvtxs_ptr;
  if ((assignment = (short *)malloc(sizeof(short)*nvtxs)) == NULL)
    {
      *status = 1;
      return;
    }
  
  /* ifdef this out if we aren't using Chaco, to avoid link thime error */
#ifdef USE_CHACO

  eigtol = (double)0.0001;
    
  interface(nvtxs, start, adjacency,
	    /* vwgts  */        (int *)NULL,
	    /* ewgts  */        (float *)NULL,
	    /* coords*/         (float *)NULL, (float *)NULL, (float *)NULL,
	    /* outassignname */ (char *)NULL,
	    /* outfilename   */ (char *)NULL,
	                         assignment,
	    /* architecture */  1,
	    /* ndims_tot */     0,
	                        mesh_dims,
	    /* goal */          (double *)NULL,
	    /* global_method */ 1,
	    /* local_method */  1,
	    /* rqi_flag */      0,
	    /* vmax */          1000,
	    /* ndims */         1,
	    /* eigtol */        eigtol, 
	    /* seed  */         (long)1234567);

  for(vtx=0;vtx<nvtxs;vtx++)
    assignment_int[vtx] = (int) assignment[vtx];

#else
/* If we didn't call Chaco, then return a single color */
  for(vtx=0;vtx<nvtxs;vtx++)
    assignment_int[vtx] = 0;
#endif
  
  *status = 0;
  return;
}

void chaco_f90_wrapper2(nvtxs_ptr, nPEs_ptr, start, adjacency, ewgt, assignment_int, status)
int	*nvtxs_ptr;      /* number of vertices in full graph */
int	*start;          /* Array of start indices for edge segments */
int     *adjacency;      /* edge list data */
float   *ewgt;           /* edge weight data */
int     *assignment_int; /* set number for each vertex */
int     *nPEs_ptr;       /* Number of processors to partition for */
int     *status;         /* Return status */

{
  int	nvtxs, vtx;
  int mesh_dims[3];
  short *assignment;
  double eigtol;
  mesh_dims[0] = *nPEs_ptr;
  mesh_dims[1] = 0;
  mesh_dims[2] = 0;

  nvtxs = *nvtxs_ptr;
  if ((assignment = (short *)malloc(sizeof(short)*nvtxs)) == NULL)
    {
      *status = 1;
      return;
    }
  
  /* ifdef this out if we aren't using Chaco, to avoid link thime error */
#ifdef USE_CHACO

  eigtol = (double)0.0001;
    
  interface(nvtxs, start, adjacency,
	    /* vwgts  */        (int *)NULL,
	    /* ewgts  */        ewgt,
	    /* coords*/         (float *)NULL, (float *)NULL, (float *)NULL,
	    /* outassignname */ (char *)NULL,
	    /* outfilename   */ (char *)NULL,
	                         assignment,
	    /* architecture */  1,
	    /* ndims_tot */     0,
	                        mesh_dims,
	    /* goal */          (double *)NULL,
	    /* global_method */ 1,
	    /* local_method */  1,
	    /* rqi_flag */      0,
	    /* vmax */          1000,
	    /* ndims */         1,
	    /* eigtol */        eigtol, 
	    /* seed  */         (long)1234567);

  for(vtx=0;vtx<nvtxs;vtx++)
    assignment_int[vtx] = (int) assignment[vtx];

#else
/* If we didn't call Chaco, then return a single color */
  for(vtx=0;vtx<nvtxs;vtx++)
    assignment_int[vtx] = 0;
#endif
  
  *status = 0;
  return;
}
