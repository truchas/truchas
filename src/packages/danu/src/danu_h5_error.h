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
* danu_h5_error.h
*
*  DANU Error Handler
*
*/

#ifndef DANU_H5_ERROR_H
#define DANU_H5_ERROR_H

#include <stdio.h>

#include <hdf5.h>

/* HDF5 uses the defened data type hid_t for handles 
* This macro is useful to check for valid handles
*/
#define H5_ISA_VALID_ID(id)    ( ( H5Iis_valid((id)) == TRUE ) ? 1 : 0 )
#define H5_ISA_INVALID_ID(id)  (int)( H5Iis_valid((id)) != TRUE ? 1 : 0 )

#define H5_ISA_FILE_ID(id)        (int)( H5Iget_type((id)) == H5I_FILE ? 1 : 0 )
#define H5_ISA_GROUP_ID(id)       (int)( H5Iget_type((id)) == H5I_GROUP ? 1 : 0 )
#define H5_ISA_DATATYPE_ID(id)    (int)( H5Iget_type((id)) == H5I_DATATYPE ? 1 : 0 )
#define H5_ISA_DATASET_ID(id)     (int)( H5Iget_type((id)) == H5I_DATASET ? 1 : 0 )
#define H5_ISA_ATTR_ID(id)        (int)( H5Iget_type((id)) == H5I_ATTR ? 1 : 0 )

#define H5_RETURN_OK(status)   (int)( (status) >= 0 ? 1 : 0 )
#define H5_RETURN_FAIL(status) (int)( (status) < 0 ? 1 : 0 )

int  danu_check_h5_return(unsigned line, const char * file, const char * func, herr_t value);

void danu_h5_error_init(void);
void danu_h5_error_suppress(void);
void danu_h5_error_activate(void);
int  danu_h5_error_is_active(void);

/*
 * Macros
 *
 *    To check return value of an HDF5 call
 *
 *    DANU_CHECK_H5_RETURN(value)
 *
 */
#define DANU_CHECK_H5_RETURN(value) \
	( danu_check_h5_return(__LINE__,__FILE__,__func__,(herr_t) (value)) )



#endif


