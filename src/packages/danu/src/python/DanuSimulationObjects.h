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

#ifndef DANU_PY_SIMOBJS_H
#define DANU_PY_SIMOBJS_H

#include "DanuH5Object.h"
#include "DanuOutputObject.h"

typedef struct {

  const H5Obj *h5;

} Simulation;  

typedef struct {

  const H5Obj *h5;
  int id;
  int cycle;
  double time;

} Sequence;

typedef struct {

  const H5Obj *h5;

} Probe;



#endif

