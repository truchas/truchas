/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_sim.h
*
*  Telluride Output Simulation Module 
*
*/

#ifndef TOUT_SIM_H
#define TOUT_SIM_H

#include <hdf5.h>

/* Public define's */
#define SIM_ROOT_NAME         "Simulations"
#define SIM_MESH_GROUP_NAME   "Sim Meshes"
#define SIM_MESH_LINK_NAME    "Mesh"


herr_t  simulations_create_root_group(hid_t fid);
hid_t   simulations_open_root_group(hid_t fid);

herr_t  simulation_add(hid_t fid, const char *sim_name, hid_t *sid);
herr_t  simulation_open(hid_t fid, const char *sim_name, hid_t *sid);
herr_t  simulation_count(hid_t fid, int *cnt);
herr_t  simulation_list(hid_t fid, int num, char **sim_names, int *num_found);
herr_t  simulation_exists(hid_t fid, const char * sim_name, int *exists);

herr_t  simulation_link_mesh(hid_t fid, hid_t sid, const char * meshname);
hid_t   simulation_open_mesh_link(hid_t sid);
herr_t  simulation_mesh_link_exists(hid_t sid, int *exists);


#endif
