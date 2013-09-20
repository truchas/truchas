/* *************************************************************************** *
*                                                                              *
*                                                                              *
*                             Copyright  (C) 20xx,                             *
*                      Los Alamos National Security, LLC                       *
*                                                                              *
*                             LA-CC-xxxxxx                                     *
*                                                                              *
* **************************************************************************** */

/*
* danu_fort_sim.h
*
*  DANU FORTRAN/C interfaces for simulations
*
*/

#ifndef DANU_FORT_SIM_H
#define DANU_FORT_SIM_H

#include <danu_fort_hid.h>

/* Use these defines to improve the readability of the code */
/* Basic simulation control */
#if 0
#include "Fortran2C.h"
#define simulation_exists_f            FORTRAN_FUNC_GLOBAL_(simulation_exists_f, SIMULATION_EXISTS_F)
#define simulation_count_f             FORTRAN_FUNC_GLOBAL_(simulation_count_f, SIMULATION_COUNT_F)
#define simulation_list_f              FORTRAN_FUNC_GLOBAL_(simulation_list_f, SIMULATION_LIST_F)
#define simulation_add_f               FORTRAN_FUNC_GLOBAL_(simulation_add_f, SIMULATION_ADD_F)
#define simulation_open_f              FORTRAN_FUNC_GLOBAL_(simulation_open_f, SIMULATION_OPEN_F)
#define simulation_link_mesh_f         FORTRAN_FUNC_GLOBAL_(simulation_link_mesh_f, SIMULATION_LINK_MESH_F)
#endif


/* Function prototypes corresponding interface definition in module_iface.f90 */
void simulation_exists_f(const hid_t_ptr *fid, const char *sim_name, const int *flen, int *flag, int *ierr);
void simulation_count_f(const hid_t_ptr *fid,  int *count, int *ierr);
void simulation_list_f(const hid_t_ptr *fid,  char *names, const int *flen, const int *num, int *ierr);

void simulation_add_f(const hid_t_ptr *fid, const char *name, const int *flen, hid_t_ptr *sid, int *ierr); 
void simulation_open_f(const hid_t_ptr *fid, const char *name, const int *flen, hid_t_ptr *sid, int *ierr); 

void simulation_link_mesh_f(const hid_t_ptr *fid, const hid_t_ptr *sid, const char * mesh_name, const int *flen, int *ierr); 
void simulation_open_mesh_link_f(const hid_t_ptr *sid, hid_t_ptr *mid, int *ierr); 
void simulation_mesh_link_exists_f(const hid_t_ptr *sid, int *flag, int *ierr); 

                              

#endif /* DANU_FORT_SIM_H */

