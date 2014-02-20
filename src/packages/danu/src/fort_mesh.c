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
*  DANU Fortran/C Mesh Interface
*
*/
#include<string.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_memory.h>
#include <danu_error.h>
#include <danu_mesh.h>
#include <danu_fort_error.h>
#include <danu_fort_strings.h>
#include <danu_fort_mesh.h>

/* Basic mesh control */

void mesh_add_unstructured_f(const hid_t_ptr *ptr, 
                             const char *name, 
                             int *flen,
                             int *elemorder,
                             int *dim,
                             hid_t_ptr *mesh,
                             int *ierr)
{
    hid_t fid = GET_HID_VALUE(ptr);
    herr_t err;
    hid_t mid;
    int  c_len        = *flen+1;
    char *c_name      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(name,*flen,c_name,c_len);

    mid = H5I_INVALID_HID;
    if ( status == DANU_SUCCESS ) {
        err = mesh_add_unstructured(fid,c_name,*elemorder,*dim,&mid);
        DANU_FREE(c_name);
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
        err = DANU_FAILURE;
    }
    *mesh = create_hid_struct(mid);
    error_translate(err,ierr);

}

void mesh_exists_f(const hid_t_ptr *fptr, const char *meshname, const int *len, int *flag, int *ierr)
{
    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err;
    int  c_len        = *len+1;
    char *c_mesh      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(meshname,*len,c_mesh,c_len);

    err = DANU_FAILURE;
    if ( status == DANU_SUCCESS ) {
        err = mesh_exists(fid,c_mesh,flag);
        DANU_FREE(c_mesh);
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
    }    
     
    error_translate(err,ierr);
    
}

void mesh_count_f(const hid_t_ptr *fptr, int *count, int *ierr)
{
    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err;

    *count = 0;
    err = mesh_count(fid,count);

    error_translate(err,ierr);

}

void mesh_list_f(const hid_t_ptr *fptr, char *names, const int *len, const int *num, int *ierr)
{

    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err, flag;
    danu_err_t copy_ok;
    int i, num_found;
    char **mesh_names;

    mesh_names = DANU_MALLOC(char *, *num);
    
    err = mesh_list(fid,*num,mesh_names,&num_found);

    if ( H5_RETURN_OK(err) && num_found > 0) {
      copy_ok =
       copy_c_array_to_fortran_array((const char * const *)mesh_names,num_found,names,num,len);

      if ( DANU_RETURN_OK(copy_ok) && ( *num >= num_found ) ) {
	err = DANU_SUCCESS;
      }
      else {

	if ( ! DANU_RETURN_OK(copy_ok) ) {
	  DANU_ERROR_MESS("Failed to copy C strings to fortran array");
	}

	if ( *num < num_found ) {
	  DANU_ERROR_MESS("mesh_names array too small to hold all mesh"
			  " names");
	}

	err = DANU_FAILURE;
      }


    }
    else {

      if ( H5_RETURN_FAIL(err) ) {
        DANU_ERROR_MESS("Failed to read mesh name list");
      }

      if ( num_found == 0 ) {
	DANU_WARN_MESS("No meshes found");
      }

      err = DANU_FAILURE;
    }

    for(i=0;i<*num;i++)
      DANU_FREE(mesh_names[i]);

    DANU_FREE(mesh_names);


    error_translate(err,ierr);
        
}

void mesh_open_f(const hid_t_ptr *fptr,
                 const char *mesh_name,
                 const int *flen,
                 hid_t_ptr *mesh,
                 int *ierr)
{
    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err;
    hid_t mid;
    int  c_len        = *flen+1;
    char *c_name      = DANU_MALLOC(char, c_len);
    danu_err_t status = convert_string_f2c(mesh_name,*flen,c_name,c_len);

    mid = H5I_INVALID_HID;
    if ( status == DANU_SUCCESS ) {
        mid = mesh_open(fid,c_name);
        DANU_FREE(c_name);
        if ( H5_ISA_VALID_ID(mid) ) {
          err = DANU_SUCCESS;
        }
        else {
          err = DANU_FAILURE;
        }
    }
    else {
        DANU_ERROR_MESS("Failed to convert Fortran string to C string");
        err = DANU_FAILURE;
    }
    *mesh = create_hid_struct(mid);

    error_translate(err,ierr);

}

void mesh_write_coordinates_f(const hid_t_ptr *mptr, 
                              int *nnodes,
                              double *x,
                              double *y,
                              double *z,
                              int *ierr)
{
    hid_t mid = GET_HID_VALUE(mptr);
    herr_t err;

    err = mesh_write_coordinates(mid,*nnodes,x,y,z);
  
    error_translate(err,ierr);

}

void mesh_write_coordinates_1d_f(const hid_t_ptr *mptr,
                                 int *nnodes,
                                 double *x,
                                 int *ierr)
{
    hid_t mid = GET_HID_VALUE(mptr);
    herr_t err;

    err = mesh_write_coordinates_1d(mid,*nnodes,x);
  
    error_translate(err,ierr);

}

void mesh_write_coordinates_2d_f(const hid_t_ptr *mptr,
                                 int *nnodes,
                                 double *x,
                                 double *y,
                                 int *ierr)
{
    hid_t mid = GET_HID_VALUE(mptr);
    herr_t err;

    err = mesh_write_coordinates_2d(mid,*nnodes,x,y);
  
    error_translate(err,ierr);

}


void mesh_read_coordinates_f(const hid_t_ptr *mptr, 
                             double *x,
                             double *y,
                             double *z,
                             int *ierr)
{
    hid_t mid = GET_HID_VALUE(mptr);
    herr_t err;

    err = mesh_read_coordinates(mid,x,y,z);
  
    error_translate(err,ierr);

}

void mesh_read_coordinates_byindex_f(const hid_t_ptr *mptr,
			             const int *f_idx,
				     double *buf,
				     int *ierr)
{
    hid_t mid = GET_HID_VALUE(mptr);
    herr_t err;
    int c_idx;

    c_idx = *f_idx - 1;
    err = mesh_read_coordinates_byindex(mid,c_idx,buf);
    error_translate(err,ierr);

}


void mesh_read_coordinates_1d_f(const hid_t_ptr *mptr,
                                 double *x,
                                 int *ierr)
{
    hid_t mid = GET_HID_VALUE(mptr);
    herr_t err;

    err = mesh_read_coordinates_1d(mid,x);
  
    error_translate(err,ierr);

}

void mesh_read_coordinates_2d_f(const hid_t_ptr *mptr,
                                 double *x,
                                 double *y,
                                 int *ierr)
{
    hid_t mid = GET_HID_VALUE(mptr);
    herr_t err;

    err = mesh_read_coordinates_2d(mid,x,y);
  
    error_translate(err,ierr);

}

/* Connectivity */

void mesh_write_connectivity_f(const hid_t_ptr *mptr,
                               const int *nelem,
                               const int *data,
                               int *ierr)
{

    hid_t mid = GET_HID_VALUE(mptr);
    herr_t err;

    err = mesh_write_connectivity(mid,*nelem,data);
    err &= mesh_connectivity_set_offset(mid,1);

    error_translate(err,ierr);

}

void mesh_read_connectivity_f(const hid_t_ptr *mptr,
                              int *data,
                              int *ierr)
{

    hid_t mid = GET_HID_VALUE(mptr);
    herr_t err;

    err = mesh_read_connectivity(mid,data);

    error_translate(err,ierr);

}

/* Mesh Attributes */
void mesh_get_type_f(const hid_t_ptr *mptr, tmesh_t *type, int *ierr)
{ 
  hid_t mid = GET_HID_VALUE(mptr);
  herr_t err = mesh_get_type(mid,type);

  error_translate(err,ierr);
}

void mesh_get_elementtype_f(const hid_t_ptr *mptr, telem_t *elem_type, int *ierr)
{ 
  hid_t mid = GET_HID_VALUE(mptr);
  herr_t err = mesh_get_elementtype(mid,elem_type);

  error_translate(err,ierr);
}

void mesh_get_dimension_f(const hid_t_ptr *mptr, int *dim, int *ierr)
{ 
  hid_t mid = GET_HID_VALUE(mptr);
  herr_t err = mesh_get_dimension(mid,dim);

  error_translate(err,ierr);
}

void mesh_get_nnodes_f(const hid_t_ptr *mptr, int *nnodes, int *ierr)
{ 
  hid_t mid = GET_HID_VALUE(mptr);
  herr_t err = mesh_get_nnodes(mid,nnodes);

  error_translate(err,ierr);
}

void mesh_get_nelem_f(const hid_t_ptr *mptr, int *nelem, int *ierr)
{ 
  hid_t mid = GET_HID_VALUE(mptr);
  herr_t err = mesh_get_nelem(mid,nelem);

  error_translate(err,ierr);
}

void mesh_get_elem_order_f(const hid_t_ptr *mptr, int *elem_order, int *ierr)
{ 
  hid_t mid = GET_HID_VALUE(mptr);
  herr_t err = mesh_get_elem_order(mid,elem_order);

  error_translate(err,ierr);
}

void mesh_create_hex_unstruct_f(const hid_t_ptr *fptr,
                                const char *name,
                                int * flen,
                                int * nnodes,
                                double *x,
                                double *y,
                                double *z,
                                int *nelem,
                                int *conn,
				hid_t_ptr *mptr,
				int *ierr)
{
    hid_t fid = GET_HID_VALUE(fptr);
    herr_t err;
    hid_t mid;
    int order = HEX_ELEM_ORDER;
    int dim = 3;
   
    mesh_add_unstructured_f(fptr,name,flen,&order,&dim,mptr,ierr);
    if ( *ierr == 0 ) {
      mesh_write_coordinates_f(mptr,nnodes,x,y,z,ierr);
      mesh_write_connectivity_f(mptr,nelem,conn,ierr);
    }
    else {
      DANU_ERROR_MESS("Failed to add unstructured hex mesh");
    }

}









