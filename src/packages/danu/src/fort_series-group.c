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
*  DANU Fortran/C Series Group Interface
*
*/
#include<string.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_memory.h>
#include <danu_error.h>
#include <danu_series-group.h>
#include <danu_sim.h>
#include <danu_fort_error.h>
#include <danu_fort_strings.h>
#include <danu_fort_series-group.h>


/* Series group queries */
void sequence_exists_f(const hid_t_ptr *sptr, const char * sname, const int *len, int *flag, int *ierr)
{
    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;
    int  c_len        = *len+1;
    char *c_data      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(sname,*len,c_data,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
      err = sequence_exists(sid,c_data,flag);
      DANU_FREE(c_data);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C");
    }

    error_translate(err,ierr);

}

void sequence_count_f(const hid_t_ptr *sptr, int * num, int *ierr)
{
    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;

    *num = 0;
    err = sequence_count(sid,num);

    error_translate(err,ierr);

}

void sequence_get_next_id_f(const hid_t_ptr *sptr, const int *cycle, const double * t, hid_t_ptr *nptr, int * ierr)
{
    hid_t sid = GET_HID_VALUE(sptr);
    hid_t hid;
    herr_t err;

    err = sequence_getNextID(sid,*cycle,*t,&hid);
    *nptr = NULL;
    if ( err == DANU_SUCCESS ) {
      *nptr = create_hid_struct(hid);
    }
    else {
      DANU_ERROR_MESS("Failed to define next series group id");
    }

    error_translate(err,ierr);

}

void sequence_list_f(const hid_t_ptr *sptr, char *names, const int *len, const int *num, int *ierr) 
{

    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;
    danu_err_t copy_ok;
    int i, num_found;
    char **data_names;
    int min_free;

    data_names = DANU_MALLOC(char *, *num);
    
    err = sequence_list(sid,*num,data_names,&num_found);

    if ( H5_RETURN_OK(err) && num_found > 0 ) {
       copy_ok =
       copy_c_array_to_fortran_array((const char* const *)data_names,num_found,names,num,len);
       if ( DANU_RETURN_OK(copy_ok) && ( *num >= num_found ) ) {
	err = DANU_SUCCESS;
      }
      else {

	if ( ! DANU_RETURN_OK(copy_ok) ) {
	  DANU_ERROR_MESS("Failed to copy C strings to fortran array");
	}

	if ( *num < num_found ) {
	  DANU_ERROR_MESS("data_names array too small to hold all data"
			  " names");
	}

	err = DANU_FAILURE;
      }
   
      min_free = num_found < *num ? num_found : *num;
      for(i=0;i<min_free;i++)
        DANU_FREE(data_names[i]);

    }
    else {

      if ( H5_RETURN_FAIL(err) ) {
        DANU_ERROR_MESS("Failed to read data name list");
      }

      if ( num_found == 0 ) {
	DANU_WARN_MESS("No series groups found");
      }

      err = DANU_FAILURE;
    }


    DANU_FREE(data_names);


    error_translate(err,ierr);
        
}

void sequence_get_id_f(const hid_t_ptr *sptr, const char *sname, const int *len, hid_t_ptr *nid, int *ierr)
{
    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;
    hid_t hid;
    int  c_len        = *len+1;
    char *c_data      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(sname,*len,c_data,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
      err = sequence_get_handle(sid,c_data,&hid);
      *nid = create_hid_struct(hid);
      DANU_FREE(c_data);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }

    error_translate(err,ierr);

}

void sequence_get_id_bynum_f(const hid_t_ptr *sptr, const int * num, hid_t_ptr *nptr, int *ierr)
{
    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;
    hid_t hid;
    char * name;

    name = sequence_get_name(*num);

    if ( name != NULL ) {
      err = sequence_get_handle(sid,name,&hid);
      *nptr = create_hid_struct(hid);
      DANU_FREE(name);
    }
    else {
      DANU_ERROR_MESS("Failed to retrieve series group name");
    }

    error_translate(err,ierr);

}




