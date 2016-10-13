/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * Danu Simulations Objects
 *
 * Interface for Simulations Objects
 *
 */

/* System includes */
#include <string.h>

#include <hdf5.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/* Danu includes */
#include <danu.h>

#include "DanuH5Object.h" 
#include "DanuOutput.h"
#include "DanuPyException.h"

#include "DanuSimulations.h"


/* Simulation Object interfaces */
Simulation * create_simulation_object(const Output *file, const char *sim_name, int newsim)
{
  Simulation * sim = NULL;
  hid_t hid = GET_OUTPUT_OBJ_HID(file);
  hid_t sid;
  herr_t stat;

  if ( newsim ) {
    stat=simulation_add(hid,sim_name,&sid);
  }
  else {
    stat=simulation_open(hid,sim_name,&sid);
  }

  if ( DANU_RETURN_OK(stat) ) {
    sim=DANU_MALLOC(Simulation,1);
    sim->h5=allocate_h5_object(sid,sim_name);
  }
  else {
    throw_exception("Failed to create simulation object");
  }

  return sim;
}

void destroy_simulation_object(Simulation *sim)
{
  hid_t hid;
  H5Obj *h5;
  if (!sim) {
    hid=GET_SIMOBJ_HID(sim);
    h5=(H5Obj*)GET_SIMOBJ_H5OBJECT(sim);
    danu_group_close(hid);
    deallocate_h5_object(h5);
    DANU_FREE(sim);
  }
}

herr_t simulation_data_write_iface(const Simulation *sim, 
                                         const char *data_name,
                                         int ndim,
                                         const int * dimensions,
                                         int typecode,
                                         const void *data)
{
  hid_t sid;
  herr_t stat=DANU_FAILURE;

  if (!sim) {
    throw_exception("Invalid Simulation object pointer");
    return stat;
  }

  sid=GET_SIMOBJ_HID(sim);

  switch(typecode) {
    case PyArray_DOUBLE:
      stat=data_write_double(sid,data_name,ndim,dimensions,(double*)data);
      break;
    case PyArray_FLOAT:
      stat=data_write_float(sid,data_name,ndim,dimensions,(float*)data);
      break;
    case PyArray_INT:
      stat=data_write_int(sid,data_name,ndim,dimensions,(int*)data);
      break;
    default:
      throw_exception("Invalid simulation non-series data type");
  }
  
  return stat;

}

herr_t simulation_data_read_iface(const Simulation *sim, 
                                         const char *data_name,
                                         int ndim,
                                         const int * dimensions,
                                         int typecode,
                                         void *data)
{
  hid_t sid;
  herr_t stat=DANU_FAILURE;

  if (!sim) {
    throw_exception("Invalid Simulation object pointer");
    return stat;
  }

  sid=GET_SIMOBJ_HID(sim);

  switch(typecode) {
    case PyArray_DOUBLE:
      stat=data_read_double(sid,data_name,ndim,dimensions,(double*)data);
      break;
    case PyArray_FLOAT:
      stat=data_read_float(sid,data_name,ndim,dimensions,(float*)data);
      break;
    case PyArray_INT:
      stat=data_read_int(sid,data_name,ndim,dimensions,(int*)data);
      break;
    default:
      throw_exception("Invalid simulation non-series data type");
  }
  
  return stat;

}

/* Sequence Object API */
Sequence * allocate_sequence_object(hid_t nsid, const char *seriesname, int id, int cycle, double time)
{
  Sequence * seq=NULL;

  seq =  DANU_MALLOC(Sequence,1);
  if ( seq ) {
    seq->h5=allocate_h5_object(nsid,seriesname);
    seq->cycle=cycle;
    seq->id=id;
    seq->time=time;
  }
  else {
    throw_exception("Failed to allocate Sequence object");
  }

  return seq;

}
void destroy_sequence_object(Sequence *seq)
{
  hid_t nsid = GET_SEQOBJ_HID(seq);
  H5Obj *h5 = GET_SEQOBJ_H5OBJECT(seq);

  /* Close the group */
  danu_group_close(nsid);
  deallocate_h5_object(h5);

  DANU_FREE(seq);

}
Sequence * create_sequence_object(const Simulation *sim, const char *seriesname, int cycle, double time)
{
  hid_t sid=GET_SIMOBJ_HID(sim);
  hid_t nsid;
  herr_t stat;
  Sequence *seq=NULL;
  char * use_name;
  int id, use_cycle;
  double use_time;

  if ( H5_ISA_VALID_ID(sid) ) {

    if ( seriesname ) {
      stat = sequence_get_handle(sid,seriesname,&nsid);
      use_name=DANU_MALLOC(char,strlen(seriesname)+1);
      strcpy(use_name,seriesname);
      stat = stat || sequence_get_id(nsid,&id);
      stat = stat || sequence_get_cycle(nsid,&use_cycle);
      stat = stat || sequence_get_time(nsid,&use_time);
    }
    else {
      stat = sequence_getNextID(sid,cycle,time,&nsid);
      stat = stat || sequence_get_id(nsid,&id);
      use_name=sequence_get_name(id);
      use_cycle=cycle;
      use_time=time;
    }

    /* Define the name, cycle and time */
    if ( DANU_RETURN_OK(stat) ) {
      seq=allocate_sequence_object(nsid,use_name,id,use_cycle,use_time);
      DANU_FREE(use_name);
    }
    else {
      throw_exception("Failed to define series group");
    }
  }
  else {
    throw_exception("Invalid simulation id");
  }

  return seq;

}

herr_t sequence_data_write_iface(const Sequence *seq, 
                                 const char *data_name,
                                 int ndim,
                                 const int * dimensions,
                                 int typecode,
                                 const void *data)
{
  hid_t nsid;
  herr_t stat=DANU_FAILURE;

  if (!seq) {
    throw_exception("Invalid Sequence object pointer");
    return stat;
  }

  nsid=GET_SEQOBJ_HID(seq);

  switch(typecode) {
    case PyArray_DOUBLE:
      stat=simulation_data_write_double(nsid,data_name,ndim,dimensions,(double*)data);
      break;
    case PyArray_FLOAT:
      stat=simulation_data_write_float(nsid,data_name,ndim,dimensions,(float*)data);
      break;
    case PyArray_INT:
      stat=simulation_data_write_int(nsid,data_name,ndim,dimensions,(int*)data);
      break;
    default:
      throw_exception("Invalid simulation series data type");
  }
  
  return stat;

}

herr_t sequence_data_read_iface(const Sequence *seq, 
                                const char *data_name,
                                int ndim,
                                const int * dimensions,
                                int typecode,
                                void *data)
{
  hid_t nsid;
  herr_t stat=DANU_FAILURE;

  if (!seq) {
    throw_exception("Invalid Sequence object pointer");
    return stat;
  }

  nsid=GET_SEQOBJ_HID(seq);

  switch(typecode) {
    case PyArray_DOUBLE:
      stat=simulation_data_read_double(nsid,data_name,ndim,dimensions,(double*)data);
      break;
    case PyArray_FLOAT:
      stat=simulation_data_read_float(nsid,data_name,ndim,dimensions,(float*)data);
      break;
    case PyArray_INT:
      stat=simulation_data_read_int(nsid,data_name,ndim,dimensions,(int*)data);
      break;
    default:
      throw_exception("Invalid simulation series data type");
  }
  
  return stat;

}

/* Probe Objects */
Probe * allocate_probe_object(hid_t pid, const char *name)
{
  Probe * probe;

  probe = DANU_MALLOC(Probe,1);
  if ( probe ) {
    probe->h5=allocate_h5_object(pid,name);
  }

  return probe;

}

void destroy_probe_object(Probe * probe)
{
  hid_t pid = GET_PRBOBJ_HID(probe);
  H5Obj *h5 = GET_PRBOBJ_H5OBJECT(probe);

  danu_dataset_close(pid);
  deallocate_h5_object(h5);

  DANU_FREE(probe);

}

Probe * create_probe_object(Simulation *sim,
			    const char *pname,
			    int len,
			    int num, 
			    int typecode,
			    void *data)
{
  hid_t sid = GET_SIMOBJ_HID(sim);
  hid_t pid;
  herr_t stat;
  int exists = 0;
  Probe * probe = NULL;


  switch(typecode) {
    case PyArray_DOUBLE:
      stat=probe_create_data_double(sid,pname,len,num,(double*)data,&pid);
      break;
    case PyArray_FLOAT:
      stat=probe_create_data_float(sid,pname,len,num,(float*)data,&pid);
      break;
    case PyArray_INT:
      stat=probe_create_data_int(sid,pname,len,num,(int*)data,&pid);
      break;
    default:
      throw_exception("Invalid simulation probe data type");
  }

  if ( DANU_RETURN_OK(stat) ) {
    probe=allocate_probe_object(pid,pname);
  }
  else {
    throw_exception("Failed to create probe dataset");
  }

  return probe;
  
}

herr_t probe_data_write_iface(const Probe *probe, 
			      int num,
                              int typecode,
                              const void *data)
{
  hid_t pid;
  herr_t stat=DANU_FAILURE;

  if (!probe) {
    throw_exception("Invalid probe object pointer");
    return stat;
  }

  pid=GET_PRBOBJ_HID(probe);

  switch(typecode) {
    case PyArray_DOUBLE:
      stat=probe_data_write_double(pid,num,(double*)data);
      break;
    case PyArray_FLOAT:
      stat=probe_data_write_float(pid,num,(float*)data);
      break;
    case PyArray_INT:
      stat=probe_data_write_int(pid,num,(int*)data);
      break;
    default:
      throw_exception("Invalid simulation probe data type");
  }
  
  return stat;

}

herr_t probe_data_read_iface(const Probe *probe, 
                             int typecode,
                             void *data)
{
  hid_t pid;
  herr_t stat=DANU_FAILURE;

  if (!probe) {
    throw_exception("Invalid Simulation object pointer");
    return stat;
  }

  pid=GET_PRBOBJ_HID(probe);

  switch(typecode) {
    case PyArray_DOUBLE:
      stat=probe_data_read_double(pid,(double*)data);
      break;
    case PyArray_FLOAT:
      stat=probe_data_read_float(pid,(float*)data);
      break;
    case PyArray_INT:
      stat=probe_data_read_int(pid,(int*)data);
      break;
    default:
      throw_exception("Invalid simulation probe data type");
  }

  return stat;

}






    










        
        



