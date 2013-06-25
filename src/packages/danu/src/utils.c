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
 * utils.c
 *
 *  DANU General Utilities 
 *
 *
 *  Purpose:
 *          This file contains routines that are useful utility functions
 *          that do not fit in other categories.
 *
 */

#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>

#include <hdf5.h>

#include <danu_error.h>
#include <danu_types.h>
#include <danu_memory.h>
//#include <danu_fort.h>
#include <danu_utils.h>

extern int errno;

/*
 * Routine: int = danu_file_access(const char * file, int * mode)
 * Purpose: Checks the access mode of file
 * Description: Return a 1 (true) or 0 (false) if the file access permissions
 *              match the mode passed in. The variable mode should be one of three
 *              macros (F_OK, R_OK, W_OK) or a bitwise OR of these macros. 
 *              F_OK => Check the existence of file
 *              R_OK => Check the read permission of file
 *              W_OK => Check the write permission of file
 *              See 'man access' for more documentation on file access calls.
 *              A file which fails the access test because the file does not 
 *              exist or does not have the correct permissions is handled silently.
 *              It is assumed that the calling program will handle these error
 *              conditions.
 *
 * Parameters:
 *           file		IN		File or pathname
 *           mode		IN		Permissions indicated by mode (F_OK, R_OK,W_OK)
 *                              
 * Returns: Returns a 0 (false) or 1 (true) result from the access test.
 * Errors: Will raise an error if mode is incorrect, an IO error occurs or
 *         file system is a read-only file system.
 *
 */ 
int danu_file_access(const char * file, int amode)
{
	int flag = FALSE;
	int ret;

	ret = access(file,amode);

	if ( ret == 0 ) {
		flag = TRUE;
	}
	else {
		switch (errno) {
			case EINVAL:
				danu_error("Invalid mode passed into danu_file_access");
				break;
			case EIO:
				danu_error_printf("An I/O error occurred while reading %s", file);
				break;
			case EROFS:
				danu_error("Write access is requested on a read-only file system");
				break;
		        case EACCES: /* Quiet failure */
                                break;
			case ENOENT: /*Quiet failure */
                                break;
			default:
				danu_error_printf("Unknown access error while accessing %s",file);
				break;
		}
                flag = FALSE;
      
	}

	return flag;
}

herr_t danu_rand_data_int(int min, int max, size_t size, int * buffer)
{
    time_t seconds;
    double scale;
    int shift;
    int r;
    didx_t i;
    int * p = buffer;

    /* Seed the generator */
    time(&seconds);
    srand((unsigned int)seconds); 

    if ( min>max ) {
        scale = (double)(min-max+1);
        shift = max;
    }
    else {
        scale = (double)(max-min+1);
        shift = min;
    }

    for(i=0;i<size;i++) {
        r =  (int) ( scale* rand()/ (RAND_MAX + 1.0) );
       *p = r + shift;
        p++;
    }

    return 0;
}



herr_t danu_rand_data_double( double min, double max, size_t  size, double * buffer)
{
    unsigned int seed = (unsigned int) time(NULL);
    didx_t i;

    double scale;
    double shift;
    double r;
    double * p = buffer;
	
    /* Compute the scale and shift */
    if ( min>max ) {
        scale = min - max;
        shift = max;
    }
    else {
        scale = max - min;
        shift = min;
    }

    srandom(seed);
    for(i=0;i<size;i++) {
        r = (double) random()/ ( (double) RAND_MAX + 1.0); /* Random number between 0.0 and 1.0 */
       *p = r*scale+shift;
        p++;
    }

    return 0;
}
/*
 * Routine: int = danu_file_delete(const char * name)
 * Purpose: Removes file (or directory) name
 * Description: Return DANU_SUCCESS or DAUN_FAIL if the is removed. Checks for the
 *              existence of the file and write permissions before deleting. Errors are possible
 *              in the unlink or rmdir call. Consult system documentation for an explanation
 *              for those error codes.
 *
 * Parameters:
 *           name       IN		File or pathname
 *                              
 * Returns: Returns a DANU_SUCCESS or DANU_FAIL
 * Errors: Will raise an error if file does not exist or caller does not have
 *         write permission. Other errors are raised in the unlink or rmdir
 *         calls. Consult system documentation for full explaination of thise errors.
 *
 */ 
danu_err_t danu_file_delete(const char * name )
{
    int    rm_stat = -1;

    if ( DANU_BAD_PTR(name) ) {
        DANU_ERROR_MESS("Invalid char pointer argument");
        return DANU_FAILURE;
    }

    if ( DANU_EMPTY_STRING(name) ) {
        DANU_ERROR_MESS("Empty string argument");
        return DANU_FAILURE;
    }

    if ( FILE_WRITE_OK(name) ) {
        rm_stat = remove(name);
    }

    if ( rm_stat != 0 ) {
        danu_error_printf("remove call failed errno set to %d",errno); 
        return DANU_FAILURE;
    }
    else {
        return DANU_SUCCESS;
    }
}



herr_t convert_int_to_hsize(int size, const int *array, hsize_t *out)
{
    herr_t    status = DANU_FAILURE;
    int       i;

    if ( out != NULL ) {
        status = DANU_SUCCESS;
        for(i=0;i<size;i++)
            out[i] = (hsize_t) array[i];
    }

    return status;

}

danu_err_t convert_int_to_size(int size, const int *array, size_t *out)
{
    danu_err_t status = DANU_FAILURE;
    int       i;


    if ( out != NULL ) {
        status = DANU_SUCCESS;
        for(i=0;i<size;i++)
            out[i] = (size_t) array[i];
    }

    return status;

}

void print_string(const char * string, int len)
{
    const char *ptr;
    int i;

    ptr = string;
    i = 0;
    printf("Printing a char array length %d\n", len);
    while(i < len) {
	if ( *ptr == 0 ) {
	    printf("Yep this is NULL===>");
	}
	if ( isspace(*ptr) ) {
	    printf("Yep this is blank===>");
	}
	printf("i=%d CHAR=<%c>\n",i, *ptr);
	ptr++;
	i++;
    }
}
    

void danu_print_hid_info(hid_t hid)
{
  char name[128];
  char type[64];
  H5I_type_t h5_type;
  
  name[0]='\0';
  H5Iget_name(hid,name,128);
  h5_type = H5Iget_type(hid);
  switch(h5_type) {
  case H5I_FILE:
    sprintf(type,"%s","FILE");
    break;
  case H5I_GROUP:
    sprintf(type,"%s","GROUP");
    break;
  case H5I_DATATYPE:
    sprintf(type,"%s","DATATYPE");
    break;
  case H5I_DATASPACE:
    sprintf(type,"%s","DATASPACE");
    break;
  case H5I_DATASET:
    sprintf(type,"%s","DATASET");
    break;
  case H5I_ATTR:
    sprintf(type,"%s","ATTR");
    break;
  case H5I_BADID:
    sprintf(type,"%s","BADID");
    break;
  default:
    sprintf(type,"%s","UNKNOWNN");
  }
  printf("HDF5 ID = 0x%08lx\n", hid);
  printf("\tname=%s\n",name);
  printf("\ttype=%s\n",type);

}




