/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
*
*  DANU Fortran/C Error Translator
*
*/
#include<string.h>

#include <hdf5.h>

#include <danu_error.h>
#include <danu_fort_error.h>

void error_translate(herr_t err, int * ierr)
{
    if ( H5_RETURN_OK(err) ) {
	*ierr = 0;
    }
    else {
	*ierr = 1;
    }
}

