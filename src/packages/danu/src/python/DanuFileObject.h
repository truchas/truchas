/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * Danu File Objects
 *
 * Python Interface File Objects
 *
 */

#ifndef DANU_PY_SWIG_FILE_H
#define DANU_PY_SWIG_FILE_H

#include "DanuH5Object.h"

typedef struct {

  H5Obj *h5;
  int  access;

} File;  


#endif

