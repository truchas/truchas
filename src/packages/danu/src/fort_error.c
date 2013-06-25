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

