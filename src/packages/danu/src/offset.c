/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
* 
*
*
*/

#include <hdf5.h>

#include <danu_error.h>
#include <danu_memory.h>
#include <danu_attribute.h>
#include <danu_offset.h>


/*
 * Routine: herr_t danu_set_offset(hid_t loc, int offset)
 * Purpose: Define an offset attribute and set it's value
 * Description: Creates, if needed, an integer attribute that
 *              that defines the index offset for dataset loc.
 *              Routine can be called multiple times to alter the
 *              offset value.
 *
 * Parameters:
 *           loc       IN              Dataset HDf5 identifier
 *           offset    IN              Offset value
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero 
 * Errors: Possible error conditions are invalid HDF5 dataset id or fails
 *         to create or update the offset value.
 */
herr_t danu_set_offset(hid_t loc, int offset)
{
  herr_t status = DANU_FAILURE;

  if ( ! H5_ISA_DATASET_ID(loc) ) {
    DANU_ERROR_MESS("Invalid dataset id when setting offset");
    return status;
  }

  status = danu_attr_write_int(loc,DANU_OFFSET_ATTR_NAME,offset);

  return status;
}

/*
 * Routine: herr_t danu_get_offset(hid_t loc, int *offset)
 * Purpose: Return the offset attribute value
 * Description: Queries the offset attribute value for dataset
 *              loc and returns the value if it exists.
 *
 * Parameters:
 *           loc       IN              Dataset HDf5 identifier
 *           offset    IN              Offset value
 *                              
 * Returns: A negative value if an error occurs, otherwise returns zero 
 * Errors: Possible error conditions are invalid HDF5 dataset id or fails
 *         to create or update the offset value.
 */
herr_t danu_get_offset(hid_t loc, int *offset)
{
  herr_t status = DANU_FAILURE;

  if ( ! H5_ISA_DATASET_ID(loc) ) {
    DANU_ERROR_MESS("Invalid dataset id when fetching offset");
    return status;
  }

  status = danu_attr_read_int(loc,DANU_OFFSET_ATTR_NAME,offset);

  return status;
}



