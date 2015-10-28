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
#include <danu_fort_transpose.h>
#include <danu_fort_strings.h>
#include <danu_fort_series-data.h>


/* Series data queries */
void simulation_data_exists_f(const hid_t_ptr *nptr,
          	              const char * dataname, 
		              const int *len,
		              int *flag, int *ierr)
{
    hid_t nid = GET_HID_VALUE(nptr);
    herr_t err;
    int  c_len        = *len+1;
    char *c_data      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(dataname,*len,c_data,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
      err = simulation_data_exists(nid,c_data,flag);
      DANU_FREE(c_data);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C");
    }

    error_translate(err,ierr);

}

void simulation_data_count_f(const hid_t_ptr *nptr, int * count, int *ierr)
{
    hid_t nid = GET_HID_VALUE(nptr);
    herr_t err;

    *count = 0;
    err = simulation_data_count(nid,count);

    error_translate(err,ierr);

}

void simulation_data_rank_f(const hid_t_ptr *nptr, 
                            const char *dataname, const int *len,
                           int * rank, int *ierr)
{
    hid_t nid = GET_HID_VALUE(nptr);
    herr_t err;
    int  c_len        = *len+1;
    char *c_data      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(dataname,*len,c_data,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
      *rank = simulation_data_rank(nid,c_data);
      DANU_FREE(c_data);
      if ( *rank <= 0 ) {
        err=DANU_FAILURE;
      }
      else {
        err=DANU_SUCCESS;
      }
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C");
    }

    error_translate(err,ierr);

}

void simulation_data_dimensions_f(const hid_t_ptr *nptr, 
                                  const char *dataname, const int *len,
                                  const int * size, int *dimensions, int *ierr)
{
    hid_t nid = GET_HID_VALUE(nptr);
    herr_t err;
    int  c_len        = *len+1;
    char *c_data      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(dataname,*len,c_data,c_len);
    hsize_t * c_dims;
    int i,j,ndims;

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
      ndims = *size;
      c_dims=DANU_MALLOC(hsize_t,ndims);
      err=simulation_data_dimensions(nid,c_data,ndims,c_dims);
      transpose_array_hsize(ndims,c_dims,dimensions);
      DANU_FREE(c_data);
      DANU_FREE(c_dims);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C");
    }

    error_translate(err,ierr);

}

void simulation_data_type_f(const hid_t_ptr *nptr, 
                            const char *dataname, const int *len,
                           int * typecode, int *ierr)
{
    hid_t nid = GET_HID_VALUE(nptr);
    herr_t err;
    int  c_len        = *len+1;
    char *c_data      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(dataname,*len,c_data,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
      err = simulation_data_type(nid,c_data,typecode);
      DANU_FREE(c_data);
    }
    else {
      DANU_ERROR_MESS("Failed to convert Fortran string to C");
    }

    error_translate(err,ierr);

}


void simulation_data_list_f(const hid_t_ptr *sptr, char *names, const int *len, const int *num, int *ierr) 
{

    hid_t sid = GET_HID_VALUE(sptr);
    herr_t err, flag;
    danu_err_t copy_ok;
    int i, num_found;
    char **data_names;
    int min_free;

    data_names = DANU_MALLOC(char *, *num);
    
    err = simulation_data_list(sid,*num,data_names,&num_found);

    if ( H5_RETURN_OK(err) && num_found > 0 ) {
       copy_ok =
       copy_c_array_to_fortran_array((const char * const *)data_names,num_found,names,num,len);
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
	DANU_WARN_MESS("No series datasets found");
      }

      err = DANU_FAILURE;
    }


    DANU_FREE(data_names);


    error_translate(err,ierr);
        
}


/* Data write routines */

void simulation_data_write_byte0_f(const hid_t_ptr *nptr, const char *data_name, const int *flen,
				  const int8_t *data, int *ierr)
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = simulation_data_write_byte(nid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void simulation_data_write_byte_f(const hid_t_ptr *nptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, const int8_t *data, int *ierr )
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = simulation_data_write_byte(nid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void simulation_data_write_int0_f(const hid_t_ptr *nptr, const char *data_name, const int *flen,
				  const int *data, int *ierr)
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = simulation_data_write_int(nid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void simulation_data_write_int_f(const hid_t_ptr *nptr, const char * data_name, const int *flen,
		      const int *dim, const int * dimensions, const int *data, int *ierr )
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = simulation_data_write_int(nid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void simulation_data_write_float0_f(const hid_t_ptr *nptr, const char *data_name, const int *flen,
                                    const float *data, int *ierr)
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = simulation_data_write_float(nid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void simulation_data_write_float_f(const hid_t_ptr *nptr, const char * data_name, const int *flen,
                                    const int *dim, const int * dimensions, const float *data, int *ierr )
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = simulation_data_write_float(nid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void simulation_data_write_double0_f(const hid_t_ptr *nptr, const char *data_name, const int *flen,
                                    const double *data, int *ierr)
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = simulation_data_write_double(nid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void simulation_data_write_double_f(const hid_t_ptr *nptr, const char * data_name, const int *flen,
                                    const int *dim, const int * dimensions, const double *data, int *ierr )
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = simulation_data_write_double(nid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}



/* Data read routines */

void simulation_data_read_byte0_f(const hid_t_ptr *nptr, const char *data_name, const int *flen,
				  int8_t *data, int *ierr)
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = simulation_data_read_byte(nid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void simulation_data_read_byte_f(const hid_t_ptr *nptr, const char * data_name, const int *flen,
		                 const int *dim, const int * dimensions, int8_t *data, int *ierr )
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = simulation_data_read_byte(nid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void simulation_data_read_int0_f(const hid_t_ptr *nptr, const char *data_name, const int *flen,
				  int *data, int *ierr)
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = simulation_data_read_int(nid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void simulation_data_read_int_f(const hid_t_ptr *nptr, const char * data_name, const int *flen,
		                 const int *dim, const int * dimensions, int *data, int *ierr )
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = simulation_data_read_int(nid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void simulation_data_read_float0_f(const hid_t_ptr *nptr, const char *data_name, const int *flen,
                                   float *data, int *ierr)
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = simulation_data_read_float(nid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void simulation_data_read_float_f(const hid_t_ptr *nptr, const char * data_name, const int *flen,
		                 const int *dim, const int * dimensions, float *data, int *ierr )
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = simulation_data_read_float(nid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

void simulation_data_read_double0_f(const hid_t_ptr *nptr, const char *data_name, const int *flen,
                                   double *data, int *ierr)
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int dim, size;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    dim=1;
    size=1;
    err = simulation_data_read_double(nid,c_data,dim,&size,data);
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string C string");
  }

  error_translate(err,ierr);

}

void simulation_data_read_double_f(const hid_t_ptr *nptr, const char * data_name, const int *flen,
		                 const int *dim, const int * dimensions, double *data, int *ierr )
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);
  int *c_dim;

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    c_dim = DANU_MALLOC(int,*dim);
    transpose_array_int(*dim,dimensions,c_dim);
    err = simulation_data_read_double(nid,c_data,*dim,c_dim,data);

    DANU_FREE(c_data);
    DANU_FREE(c_dim);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran string to C string");
  }

  error_translate(err,ierr);

}

/* Simulation data open */

void simulation_open_data_f(const hid_t_ptr * nptr, const char * data_name, const int * flen, 
			    hid_t_ptr * hptr, int * ierr)
{
  hid_t nid = GET_HID_VALUE(nptr);
  herr_t err;
  hid_t hid;
  int  c_len        = *flen+1;
  char *c_data       = DANU_MALLOC(char, c_len);
  danu_err_t status = convert_string_f2c(data_name,*flen,c_data,c_len);

  err = DANU_FAILURE;
  if ( status == DANU_SUCCESS ) {
    hid = simulation_open_data(nid,c_data);
    if ( H5_ISA_VALID_ID(hid) ) {
      *hptr = create_hid_struct(hid);
      err = DANU_SUCCESS;
    }
    DANU_FREE(c_data);
  }
  else {
    DANU_ERROR_MESS("Failed to copy Fortran to C string");
  }

  error_translate(err,ierr);

}


