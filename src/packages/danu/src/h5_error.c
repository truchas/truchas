/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*
 * h5_error.c
 *
 *  DANU wrappers to HDF5 error handler calls
 *
 *  Purpose:
 *           The HDF5 error handling is not useful for most applications. The
 *           library builds a stack of error messages organized by major and 
 *           minor character strings and numbers. By default, this stack rolled
 *           out to STDOUT has the call stack is rolled back to the initial
 *           entry point in the library. For a calling program, this approach
 *           is not useful since the top of the stack will not include the
 *           program's calling routine. 
 *           DANU has it's own error handling and will check common
 *           errors such missing files, invalid permissions and return
 *           values from the HDF5 library calls. It will disable the
 *           HDF5 error handler by default. However calls to activate the 
 *           handler will be included along with calls to push messages onto
 *           the error stack and print out the current HDF5 error stack. A
 *           macro, DANU_H5EPUSH_MESS has been defined in danu_h5_error.h
 *           to simplify pushing messages to the error stack.
 *
 *
 */

#if HAVE_CONFIG_H
#include <danu_config.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_utils.h>
#include <danu_error.h>

#include<danu_h5_error.h>




/*
 * Static variables that indicate if the error stack printing
 * is active (H5_ErrorStack_Active) and pointers to store the
 * function and the function's data for the default error
 * stack H5E_DEFAULT. 
 */
static int         H5_ErrorStack_Active = 0;
static H5E_auto_t  H5_ErrorStack_FUNC;
static void*       H5_ErrorStack_Data = NULL;




/*
 * Routine: void danu_h5_error_suppress(void)
 * Purpose: Suppress HDF5 error messages
 * Description: This routine suppresses the error messages
 *              from the HDF5 library which can be prolific
 *              and useless. It assumes that only the default
 *              error handler is active and will only deactivate
 *              this handler. To reactivate the default handler
 *              call danu_h5_error_activate().
 * *
 * Parameters: None.
 * Returns: Nothing.
 * Errors: Will print an error message if the call returns a negative
 *         value.
 *
 */ 
void danu_h5_error_suppress(void)
{
    herr_t ret;

    ret = H5Eset_auto(H5E_DEFAULT,NULL,NULL);
    DANU_CHECK_H5_RETURN(ret);

    H5_ErrorStack_Active = 0;
}

/*
 * Routine: void danu_h5_error_activate(void)
 * Purpose: Allow the printing of the default HDF5 error stack
 * Description: This routine activates the error stack printing
 *              in HDF5. The error stack is the default error stack
 *              in HDF5. The routine restores the print function
 *              and data that is stored in the static variables
 *              H5_ErrorStack_FUNC and H5_ErrorStack_Data, if 
 *              the H5_ErrorStack_FUNC is not NULL. 
 *
 * Parameters: None.
 * Returns: Nothing.
 * Errors: Will print an error message if the call returns a negative
 *         value and will dump the HDF5 error stack if an error is
 *         detected. This is the ONLY function that will ignore the
 *         default behavior of suppressing the HDF5 error messages.
 *
 */ 
void danu_h5_error_activate(void)
{

        //if ( H5_ErrorStack_FUNC != NULL ) {
                H5Eset_auto(H5E_DEFAULT,H5_ErrorStack_FUNC,H5_ErrorStack_Data);
                H5_ErrorStack_Active = 1;
        //}

}

/*
 * Routine: int danu_h5_error_is_active(void)
 * Purpose: Returns the value of H5_ErrorStack_Active
 * Description: This routine returns the current value of H5_ErrorStack_Active
 *              A return value of 1 (TRUE) indicates that the default
 *              HDF5 Error Stack will be printed out if an error is detected
 *              in the HDF5 library. 
 *
 * Parameters: None.
 * Returns: Current value of H5_ErrorStack_Active
 * Errors: This routine should never fail.
 *
 */ 
int danu_h5_error_is_active(void)
{
        return H5_ErrorStack_Active;
}

/*
 * Routine: int danu_check_h5_return(unsigned line, char * file, char * func, value)
 * Purpose: Check the return status of an HDF5 call and return
 *          either DANU_SUCCESS or DANU_FAIL.
 * Description: Most HDF5 calls return a flag indicating if an error
 *              has been raised. A negative return value indicates
 *              an error has occurred. In this routine if value is
 *              negative an error message with the line number, file and function
 *              name is printed. User's should use the macro 
 *              DANU_CHECK_H5_RETURN(value) instead of calling this routine
 *              directly.
 * Parameters:
 *           line               IN              Line number
 *           file               IN              File name
 *           func               IN              Function name
 *           value              IN              HDF5 return value
 *
 * Returns: DANU_SUCCESS or DANU_FAIL
 *
 * Errors: This routine simply checks the return value
 *         and prints a message. Should not raise any errors
 *
 */
int danu_check_h5_return(unsigned line, const char * file, const char * func, herr_t value)
{
        int flag = DANU_SUCCESS;

        if (value < 0 ) {
           danu_error_printf("Error detected in HDF5 call in %s (File=%s,line=%d). HDF5 call returned %d\n",
                          func, file, line, value);
           flag = DANU_FAILURE;
        }

        return flag;
}
/*
 * Routine: void danu_h5_error_init(unsigned )
 * Purpose: Initialize the global static variables
 *          
 * Description: Initialize the variables H5_ErrorStack_FUNC and H5_ErrorStack_Data
 *              and then suppress the error stack printing.
 *
 * Parameters: None.
 *
 * Returns: Nothing
 *
 * Errors: Does not check the return value of H5Epush
 * 
 */
void danu_h5_error_init(void)
{
        herr_t ret;

        ret = H5Eget_auto(H5E_DEFAULT,&H5_ErrorStack_FUNC,H5_ErrorStack_Data);
        DANU_CHECK_H5_RETURN(ret);
                
        danu_h5_error_suppress();
}

/*
 * Routine: dbool_t h5_object_exists(hid_t loc, const char * object)
 * Purpose: Check the existence of object in loc
 *          
 * Description: Attempts to open object in location loc. Returns a
 *              boolean flag TRUE (exists) or FALSE (does not exist).
 *              Will not open object if loc is not a valid HDF5 id.
 *              Uses H5Oopen which opens a group, dataset or named
 *              data type.
 *           
 * Parameters: 
 *             loc              IN              HDF5 identifier location of object
 *             object           IN              String name of object
 *
 * Returns: A boolean flag 0 (FALSE) or 1 (TRUE)
 *
 * Errors: Since this routine only checks the existence of an object, no errors
 *         are raised
 * 
 */
dbool_t h5_object_exists(hid_t loc, const char * object)
{
    dbool_t flag;

    hid_t oid;

    if ( H5_ISA_INVALID_ID(loc) ) {
        flag = FALSE;
    }
    else {
        oid = H5Oopen(loc,object,H5P_DEFAULT);
        if ( H5_ISA_INVALID_ID(oid) ) {
            flag = FALSE;
        }
        else {
            flag = TRUE;
            H5Oclose(oid);
        }
    }

    return flag;
}

