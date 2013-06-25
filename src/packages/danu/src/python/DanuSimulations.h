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
 * Danu Simulation Objects
 *
 * Python Interface Simulation Objects
 *
 */

#ifndef DANU_PY_SIMS_H
#define DANU_PY_SIMS_H

#include <hdf5.h>

#include "DanuH5Object.h"
#include "DanuOutputObject.h"
#include "DanuSimulationObjects.h"

/* Object Macros */
#define GET_SIMOBJ_HID(a)           ( (a)->h5->hid )
#define GET_SIMOBJ_H5OBJECT(a)      ( (H5Obj *) (a)->h5 )

#define GET_SEQOBJ_HID(a)           ( (a)->h5->hid )
#define GET_SEQOBJ_H5OBJECT(a)      ( (H5Obj *) (a)->h5 )

#define GET_PRBOBJ_HID(a)           ( (a)->h5->hid )
#define GET_PRBOBJ_H5OBJECT(a)      ( (H5Obj *) (a)->h5 )

/* Create Simulation object */
Simulation * create_simulation_object(const Output *file, const char *sim_name, int new_obj);
void         destroy_simulation_object(Simulation *sim);
herr_t       simulation_data_write_iface(const Simulation *sim, 
                                         const char *data_name,
                                         int ndim,
                                         const int * dimensions,
                                         int typecode,
                                         const void *data);

herr_t       simulation_data_read_iface(const Simulation *sim, 
                                         const char *data_name,
                                         int ndim,
                                         const int * dimensions,
                                         int typecode,
                                         void *data);

/* Create Sequence object */
Sequence *  allocate_sequence_objecti(hid_t nsid, const char *seriesname, int id, int cycle, double time);
Sequence *  create_sequence_object(const Simulation *sim, const char *seriesname, int cycle, double time);
void        destroy_sequence_object(Sequence *seq);

herr_t       sequence_data_write_iface(const Sequence *seq, 
                                         const char *data_name,
                                         int ndim,
                                         const int * dimensions,
                                         int typecode,
                                         const void *data);

herr_t       sequence_data_read_iface(const Sequence *seq, 
                                         const char *data_name,
                                         int ndim,
                                         const int * dimensions,
                                         int typecode,
                                         void *data);

 

/* Probe Objects */

Probe * allocate_probe_object(hid_t pid, const char *name);
void    destroy_probe_object(Probe * probe);
Probe * create_probe_object(Simulation *sim,
                            const char *probename,
                            int len, int num, int typecode, void *data); 

herr_t       probe_data_write_iface(const Probe *probe, 
                                    int num,
                                    int typecode,
                                    const void *data);

herr_t       probe_data_read_iface(const Probe *probe, 
                                   int typecode,
                                   void *data);


#endif
