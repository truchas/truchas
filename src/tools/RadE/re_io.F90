#include "f90_assert.fpp"

module re_io

  use netcdf
  implicit none
  private
  
  public :: re_open_rw, re_close, re_init_vf, re_put_ia, re_put_vf_rows, re_put_ambient
  public :: re_open_ro, re_get_vf_dims, re_get_ia, re_get_vf_rows, re_get_ambient
  
contains
  
  subroutine re_open_rw (path, ncid)
    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid
    integer :: status
    status = nf90_open(path, NF90_WRITE, ncid)
    call handle_netcdf_error (status)
  end subroutine re_open_rw
  
  subroutine re_open_ro (path, ncid)
    character(len=*), intent(in) :: path
    integer, intent(out) :: ncid
    integer :: status
    status = nf90_open(path, NF90_NOWRITE, ncid)
    call handle_netcdf_error (status)
  end subroutine re_open_ro
  
  subroutine re_close (ncid)
    integer, intent(in) :: ncid
    integer :: status
    status = nf90_close (ncid)
    call handle_netcdf_error (status)
  end subroutine re_close
  
  subroutine re_init_vf (ncid, nnonz)
  
    integer, intent(in) :: ncid, nnonz
    
    integer :: status, dim1, dim2, varid, nface
    
    status = nf90_redef(ncid)
    call handle_netcdf_error (status)

    !! Get the number of faces and its dimension ID.
    status = nf90_inq_dimid(ncid, 'num_faces', dim1)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dim1, len=nface)
    call handle_netcdf_error (status)
    
    !! Define the VAL and JA arrays.
    status = nf90_def_dim(ncid, 'num_nonzero', nnonz, dim2)
    call handle_netcdf_error (status)
    status = nf90_def_var(ncid, 'val', NF90_FLOAT, (/dim2/), varid)
    call handle_netcdf_error (status)
    status = nf90_def_var(ncid, 'ja', NF90_INT, (/dim2/), varid)
    call handle_netcdf_error (status)
    
    !! Define the IA array.
    status = nf90_def_dim(ncid, 'num_ia', nface+1, dim2)
    call handle_netcdf_error (status)
    status = nf90_def_var(ncid, 'ia', NF90_INT, (/dim2/), varid)
    call handle_netcdf_error (status)
    
    !! Define the ambient view factor array.
    status = nf90_def_var(ncid, 'ambient', NF90_FLOAT, (/dim1/), varid)
    call handle_netcdf_error (status)
    
    status = nf90_enddef(ncid)
    call handle_netcdf_error (status)
    
  end subroutine re_init_vf
  
  subroutine re_get_vf_dims (ncid, nface, nnonz)
  
    integer, intent(in)  :: ncid
    integer, intent(out) :: nface, nnonz
    
    integer :: status, dimid
    
    !! Get the number of faces.
    status = nf90_inq_dimid(ncid, 'num_faces', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=nface)
    call handle_netcdf_error (status)
    
    !! Get the number of nonzeros in the VF matrix.
    status = nf90_inq_dimid(ncid, 'num_nonzero', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=nnonz)
    call handle_netcdf_error (status)

  end subroutine re_get_vf_dims
  
  subroutine re_put_ia (ncid, ia)
    integer, intent(in) :: ncid, ia(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'ia', varid)
    call handle_netcdf_error (status)
    status = nf90_put_var(ncid, varid, ia)
    call handle_netcdf_error (status)
  end subroutine re_put_ia
  
  subroutine re_get_ia (ncid, ia)
    integer, intent(in)  :: ncid
    integer, intent(out) :: ia(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'ia', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, ia)
    call handle_netcdf_error (status)
  end subroutine re_get_ia
  
  subroutine re_put_vf_rows (ncid, val, ja, start)
    integer, intent(in) :: ncid, ja(:), start
    real, intent(in) :: val(:)
    integer :: status, varid
    ASSERT( size(ja) == size(val) )
    if (size(ja) == 0) return
    !! Should check that we're not trying to write too much data
    status = nf90_inq_varid(ncid, 'val', varid)
    call handle_netcdf_error (status)
    status = nf90_put_var(ncid, varid, val, (/start/))
    call handle_netcdf_error (status)
    status = nf90_inq_varid(ncid, 'ja', varid)
    call handle_netcdf_error (status)
    status = nf90_put_var(ncid, varid, ja, (/start/))
    call handle_netcdf_error (status)
  end subroutine
  
  subroutine re_get_vf_rows (ncid, val, ja, start)
    integer, intent(in)  :: ncid, start
    integer, intent(out) :: ja(:)
    real,    intent(out) :: val(:)
    integer :: status, varid
    ASSERT( size(ja) == size(val) )
    if (size(ja) == 0) return
    !! Should check that we're not trying to read too much data
    status = nf90_inq_varid(ncid, 'val', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, val, (/start/))
    call handle_netcdf_error (status)
    status = nf90_inq_varid(ncid, 'ja', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, ja, (/start/))
    call handle_netcdf_error (status)
  end subroutine
  
  subroutine re_put_ambient (ncid, ambient)
    integer, intent(in) :: ncid
    real,    intent(in) :: ambient(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'ambient', varid)
    call handle_netcdf_error (status)
    status = nf90_put_var(ncid, varid, ambient)
    call handle_netcdf_error (status)
  end subroutine re_put_ambient
  
  subroutine re_get_ambient (ncid, ambient)
    integer, intent(in)  :: ncid
    real,    intent(out) :: ambient(:)
    integer :: status, varid
    !! Should check the length is correct
    status = nf90_inq_varid(ncid, 'ambient', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, ambient)
    call handle_netcdf_error (status)
  end subroutine re_get_ambient
  
  subroutine handle_netcdf_error (status)
    integer, intent(in) :: status
    if (status == NF90_NOERR) return
    print *, trim(nf90_strerror(status))
    stop
  end subroutine handle_netcdf_error
  
end module re_io
