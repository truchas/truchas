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

