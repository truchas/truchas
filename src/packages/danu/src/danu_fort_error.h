/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* danu_fort_error.h
*
*  DANU FORTRAN/C Error Translator
*
*/

#ifndef DANU_FORT_ERROR_H
#define DANU_FORT_ERROR_H

void error_translator(herr_t err, int *ierr);

#endif /* DANU_FORT_ERROR_H */

