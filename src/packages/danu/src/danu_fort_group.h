/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_fort_file.h
*
*  DANU FORTRAN/C interfaces for files
*
*/

#ifndef DANU_FORT_GROUP_H
#define DANU_FORT_GROUP_H

#include <danu_fort_hid.h>

/* Use these defines to improve the readability of the code */
/* Basic file control */
#if 0
#include "Fortran2C.h"
#define danu_group_exists_f            FORTRAN_FUNC_GLOBAL_(danu_group_exists_f,DANU_GROUP_EXISTS_F)

#define danu_group_create_f            FORTRAN_FUNC_GLOBAL_(danu_group_create_f,DANU_GROUP_CREATE_F)
#define danu_group_open_f            FORTRAN_FUNC_GLOBAL_(danu_group_open_f,DANU_GROUP_EXISTS_F)
#define danu_group_close_f            FORTRAN_FUNC_GLOBAL_(danu_group_close_f,DANU_GROUP_CLOSE_F)

$define danu_group_get_nlinks_f        FORTRAN_FUNC_GLOBAL_(danu_group_get_nlinks_f, DANU_GROUP_GET_NLINKS_F)
#define danu_group_get_subgroups_f     FORTRAN_FUNC_GLOBAL_(danu_group_get_subgroups_f,DANU_GROUP_GET_SUBGROUPS_F)
#define danu_group_get_datasets_f      FORTRAN_FUNC_GLOBAL_(danu_group_get_datasets_f,DANU_GROUP_GET_DATASETS_F)

#endif


/* Function prototypes corresponding interface definition in module_danu_iface.f90 */

/* Attribute exists, count and names */
void danu_group_exists_f(hid_t_ptr *ptr, const char *name, int *flen, int *exists, int *ierr);


void danu_group_create_f(hid_t_ptr *ptr, const char *name, int *flen, hid_t_ptr *gid, int *ierr); 
void danu_group_open_f(hid_t_ptr *ptr, const char *name, int *flen, hid_t_ptr *gid, int *ierr); 
void danu_group_close_f(hid_t_ptr *gid, int *ierr); 



#endif /* DANU_FORT_GROUP_H */

