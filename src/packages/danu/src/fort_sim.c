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
*  DANU Fortran/C Simulation Interface
*
*/
#include<string.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_memory.h>
#include <danu_error.h>
#include <danu_sim.h>
#include <danu_sim.h>
#include <danu_fort_error.h>
#include <danu_fort_strings.h>
#include <danu_fort_sim.h>

/* Simulation query */

void simulation_exists_f(const hid_t_ptr *fptr, const char *sim_name, const int *len, int *flag, int *ierr)
{
    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err;
    int  c_len        = *len+1;
    char *c_sim       = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(sim_name,*len,c_sim,c_len);
    
    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
        err = simulation_exists(fid,c_sim,flag);
        DANU_FREE(c_sim);
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }    
     
    error_translate(err,ierr);
    
}

void simulation_count_f(const hid_t_ptr *fptr, int *count, int *ierr)
{
    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err;

    *count = 0;
    err = simulation_count(fid,count);

    error_translate(err,ierr);

}

void simulation_list_f(const hid_t_ptr *fptr, char *names, const int *len, const int *num, int *ierr)
{

    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err, flag;
    danu_err_t copy_ok;
    int i, num_found;
    char **simulation_names;
    int min_free;

    simulation_names = DANU_MALLOC(char *, *num);
    
    err = simulation_list(fid,*num,simulation_names,&num_found);

    if ( H5_RETURN_OK(err) && num_found > 0) {
      copy_ok =
       copy_c_array_to_fortran_array((const char * const *)simulation_names,num_found,names,num,len);

      if ( DANU_RETURN_OK(copy_ok) && ( *num >= num_found ) ) {
	err = DANU_SUCCESS;
      }
      else {

	if ( ! DANU_RETURN_OK(copy_ok) ) {
	  DANU_ERROR_MESS("Failed to copy C strings to fortran array");
	}

	if ( *num < num_found ) {
	  DANU_ERROR_MESS("simulation_names array too small to hold all simulation"
			  " names");
	}

	err = DANU_FAILURE;
      }
   
      min_free = num_found < *num ? num_found : *num;
      for(i=0;i<min_free;i++)
        DANU_FREE(simulation_names[i]);

    }
    else {

      if ( H5_RETURN_FAIL(err) ) {
        DANU_ERROR_MESS("Failed to read simulation name list");
      }

      if ( num_found == 0 ) {
	DANU_WARN_MESS("No simulations found");
      }

      err = DANU_FAILURE;
    }


    DANU_FREE(simulation_names);


    error_translate(err,ierr);
        
}
void simulation_add_f(const hid_t_ptr *fptr,
                      const char *simulation_name,
                      const int *flen,
                      hid_t_ptr *simulation,
                      int *ierr)
{
    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err;
    hid_t sid;
    int  c_len        = *flen+1;
    char *c_name      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(simulation_name,*flen,c_name,c_len);

    sid = H5I_INVALID_HID;
    if ( status == DANU_SUCCESS ) {
      err = simulation_add(fid,c_name,&sid);
      DANU_FREE(c_name);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C string");
      err = DANU_FAILURE;
    }
    *simulation = create_hid_struct(sid);

    error_translate(err,ierr);

}

void simulation_open_f(const hid_t_ptr *fptr,
                       const char *simulation_name,
                       const int *flen,
                       hid_t_ptr *simulation,
                       int *ierr)
{
    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err;
    hid_t sid;
    int  c_len        = *flen+1;
    char *c_name      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(simulation_name,*flen,c_name,c_len);

    sid = H5I_INVALID_HID;
    if ( status == DANU_SUCCESS ) {
        err = simulation_open(fid,c_name,&sid);
        DANU_FREE(c_name);
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
        err = DANU_FAILURE;
    }
    *simulation = create_hid_struct(sid);

    error_translate(err,ierr);

}

void simulation_link_mesh_f(const hid_t_ptr *fptr,
                            const hid_t_ptr *sptr,
                            const char *mesh_name,
                            const int  *flen,
                            int *ierr)
{
    hid_t fid = GET_HID_VALUE(fptr);
    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err;
    int  c_len        = *flen+1;
    char *c_name      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(mesh_name,*flen,c_name,c_len);

    if ( status == DANU_SUCCESS ) {
      err = simulation_link_mesh(fid,sid,c_name);
      DANU_FREE(c_name);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C string");
      err = DANU_FAILURE;
    }

    error_translate(err,ierr);

}

void simulation_open_mesh_link_f(const hid_t_ptr *sptr,
                                 hid_t_ptr *mptr,
                                 int *ierr)
{
   hid_t sid = GET_HID_VALUE(sptr);
   hid_t mid;
   herr_t err;

   mid = simulation_open_mesh_link(sid);
   if ( H5_ISA_VALID_ID(mid) ) {
     *mptr= create_hid_struct(mid);
     err = DANU_SUCCESS;
   }
   else {
     err = DANU_FAILURE;
   }

   error_translate(err,ierr);

}

void simulation_mesh_link_exists_f(const hid_t_ptr *sptr,
                                  int *flag,
                                  int *ierr)
{
   hid_t sid = GET_HID_VALUE(sptr);
   hid_t mid;
   herr_t err;

   err = simulation_mesh_link_exists(sid,flag);
   error_translate(err,ierr);

}






