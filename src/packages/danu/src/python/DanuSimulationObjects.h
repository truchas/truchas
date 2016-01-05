/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

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

