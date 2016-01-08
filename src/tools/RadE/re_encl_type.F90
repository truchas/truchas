!!
!! RE_ENCL_TYPE
!!
!! This module provides a derived type for describing the enclosing surface of
!! a radiation enclosure and methods that operate on instances of this type.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 3 Apr 2008
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL DESTROY (THIS) deallocates the allocated storage components of THIS.
!!
!!  CALL BCAST_ENCL (THIS) broadcasts the contents of THIS on process rank 1 to
!!    all the other process ranks, replicating the enclosure on all processes.
!!
!!  CALL READ_ENCL (THIS, PATH[, HAS_VF]) reads the radiation enclosure
!!    dataset PATH and initializes the enclosure THIS with the data.  The
!!    optional logical arugment HAS_VF is assigned the value true if the
!!    dataset contains view factor data.  This is a collective procedure.
!!    Input occurs on process rank 1 and the enclosure data is then replicated
!!    on all other processes.
!!
!!  CALL WRITE_ENCL (THIS, PATH) writes the enclosure data contained in
!!    THIS to the radiation enclosure dataset PATH. Any previous file is
!!    overwritten.  This is a collective procedure, but only for the purposes
!!    of error handling.  Output occurs on process rank 1 using its enclosure
!!    data; THIS is ignored on all other processes.
!!

#include "f90_assert.fpp"

module re_encl_type

  use kinds, only: r8
  implicit none
  private

  public :: destroy, bcast_encl
  public :: read_encl, write_encl

  integer, parameter :: MAX_NAME_LEN = 32

  type, public :: encl
    character(len=MAX_NAME_LEN) :: name
    !! Bare surface representation
    integer :: nnode=0, nface=0
    integer,  pointer :: xface(:) => null()
    integer,  pointer :: fnode(:) => null()
    real(r8), pointer :: x(:,:) => null()
    !! Face group attributes
    integer,  pointer :: gnum(:) => null()
    integer,  pointer :: group_id_list(:) => null()
    !! Connection to the source mesh
    integer, pointer :: src_elem(:) => null()
    integer, pointer :: src_side(:) => null()
    !! Symmetries
    logical :: mirror(3)
    integer :: rot_axis, num_rot
  end type

  interface destroy
    module procedure destroy_encl
  end interface

contains

  subroutine destroy_encl (e)
    type(encl), intent(inout) :: e
    if (associated(e%xface)) deallocate(e%xface)
    if (associated(e%fnode)) deallocate(e%fnode)
    if (associated(e%x))     deallocate(e%x)
    if (associated(e%gnum))  deallocate(e%gnum)
    if (associated(e%group_id_list)) deallocate(e%group_id_list)
    if (associated(e%src_elem)) deallocate(e%src_elem)
    if (associated(e%src_side)) deallocate(e%src_side)
  end subroutine destroy_encl

  subroutine bcast_encl (this)

    use scl

    type(encl), intent(inout) :: this

    integer :: n, my_rank

    my_rank = scl_rank()

    call scl_bcast (this%name)
    call scl_bcast (this%nnode)
    call scl_bcast (this%nface)

    if (my_rank == 1) n = size(this%fnode)
    call scl_bcast (n)
    if (my_rank > 1) allocate(this%xface(this%nface+1), this%fnode(n), this%x(3,this%nnode))
    call scl_bcast (this%xface)
    call scl_bcast (this%fnode)
    call scl_bcast (this%x)

    if (my_rank == 1) n = size(this%group_id_list)
    call scl_bcast (n)
    if (my_rank > 1) allocate(this%group_id_list(n), this%gnum(this%nface))
    call scl_bcast (this%group_id_list)
    call scl_bcast (this%gnum)

    if (my_rank > 1) allocate(this%src_elem(this%nface), this%src_side(this%nface))
    call scl_bcast (this%src_elem)
    call scl_bcast (this%src_side)

    call scl_bcast (this%mirror)
    call scl_bcast (this%rot_axis)
    call scl_bcast (this%num_rot)

  end subroutine bcast_encl

  subroutine write_encl (this, path)

    use netcdf
    use re_io
    use scl

    type(encl),  intent(in) :: this
    character(len=*), intent(in) :: path

    integer :: status, ncid, cmode, dim1, dim2
    integer :: coord_id, xface_id, fnode_id, src_elem_id, src_side_id, gnum_id, group_id
    integer :: symmetry_id
    integer :: symmetry(5)

    if (scl_rank() > 1) return

    cmode = ior(NF90_CLOBBER, NF90_64BIT_OFFSET)
    status = nf90_create(path, cmode, ncid)
    call handle_netcdf_error (status)

    !! Define the coordinate variable.
    status = nf90_def_dim(ncid, 'num_dim', size(this%x,dim=1), dim1)
    call handle_netcdf_error (status)
    status = nf90_def_dim(ncid, 'num_nodes', size(this%x,dim=2), dim2)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'coord', NF90_DOUBLE, (/dim1,dim2/), coord_id)
    call handle_netcdf_error (status)

    !! Define the face connectivity variables.
    status = nf90_def_dim(ncid, 'num_xface', size(this%xface), dim1)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'xface', NF90_INT, (/dim1/), xface_id)
    call handle_netcdf_error (status)
    status = nf90_def_dim(ncid, 'num_fnode', size(this%fnode), dim1)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'fnode', NF90_INT, (/dim1/), fnode_id)
    call handle_netcdf_error (status)

    !! Define the gnum, src_elem and src_side variables.
    status = nf90_def_dim(ncid, 'num_faces', this%nface, dim1)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'src_elem', NF90_INT, (/dim1/), src_elem_id)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'src_side', NF90_INT, (/dim1/), src_side_id)
    call handle_netcdf_error (status)
    status = nf90_def_var (ncid, 'gnum', NF90_INT, (/dim1/), gnum_id)
    call handle_netcdf_error (status)

    !! Defined the group_id variable.
    status = nf90_def_dim(ncid, 'num_group', size(this%group_id_list), dim2)
    call handle_netcdf_error (status)
    status = nf90_def_var(ncid, 'group_ids', NF90_INT, (/dim2/), group_id)
    call handle_netcdf_error (status)

    !! Define the symmetry variable.
    status = nf90_def_dim(ncid, 'num_symmetry', size(symmetry), dim2)
    call handle_netcdf_error (status)
    status = nf90_def_var(ncid, 'symmetry', NF90_INT, (/dim2/), symmetry_id)
    call handle_netcdf_error (status)

    status = nf90_enddef (ncid)
    call handle_netcdf_error (status)

    status = nf90_put_var (ncid, coord_id, this%x)
    call handle_netcdf_error (status)

    status = nf90_put_var (ncid, xface_id, this%xface)
    call handle_netcdf_error (status)

    status = nf90_put_var (ncid, fnode_id, this%fnode)
    call handle_netcdf_error (status)

    status = nf90_put_var (ncid, gnum_id, this%gnum)
    call handle_netcdf_error (status)

    status = nf90_put_var (ncid, src_elem_id, this%src_elem)
    call handle_netcdf_error (status)

    status = nf90_put_var (ncid, src_side_id, this%src_side)
    call handle_netcdf_error (status)

    status = nf90_put_var (ncid, group_id, this%group_id_list)
    call handle_netcdf_error (status)

    symmetry = 0
    where (this%mirror) symmetry(1:3) = 1
    symmetry(4) = this%rot_axis
    if (this%rot_axis > 0) symmetry(5) = this%num_rot
    status = nf90_put_var (ncid, symmetry_id, symmetry)
    call handle_netcdf_error (status)

    call re_close (ncid)

  end subroutine write_encl

  subroutine read_encl (this, path, has_vf)

    use netcdf
    use re_io
    use scl

    type(encl),  intent(out) :: this
    character(len=*), intent(in) :: path
    logical, intent(out), optional :: has_vf

    integer :: status, ncid, dimid, varid
    integer :: n, num_dim
    integer, allocatable :: symmetry(:)

    if (scl_rank() == 1) then

    status = nf90_open(path, NF90_NOWRITE, ncid)
    call handle_netcdf_error (status)

    !! Get the number of nodes.
    status = nf90_inq_dimid(ncid, 'num_nodes', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=this%nnode)
    call handle_netcdf_error (status)

    !! Get the dimension of the coordinate space.
    status = nf90_inq_dimid(ncid, 'num_dim', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=num_dim)
    call handle_netcdf_error (status)

    !! Allocate and read the coordinate array.
    allocate(this%x(num_dim,this%nnode))
    status = nf90_inq_varid(ncid, 'coord', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, this%x)
    call handle_netcdf_error (status)


    status = nf90_inq_dimid(ncid, 'num_xface', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=n)
    call handle_netcdf_error (status)
    allocate(this%xface(n))
    status = nf90_inq_varid(ncid, 'xface', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, this%xface)
    call handle_netcdf_error (status)

    status = nf90_inq_dimid(ncid, 'num_fnode', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=n)
    call handle_netcdf_error (status)
    allocate(this%fnode(n))
    status = nf90_inq_varid(ncid, 'fnode', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, this%fnode)
    call handle_netcdf_error (status)

    status = nf90_inq_dimid(ncid, 'num_faces', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=n)
    call handle_netcdf_error (status)
    this%nface = n
    allocate(this%gnum(n), this%src_elem(n), this%src_side(n))
    status = nf90_inq_varid(ncid, 'gnum', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, this%gnum)
    call handle_netcdf_error (status)
    status = nf90_inq_varid(ncid, 'src_elem', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, this%src_elem)
    call handle_netcdf_error (status)
    status = nf90_inq_varid(ncid, 'src_side', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, this%src_side)
    call handle_netcdf_error (status)

    status = nf90_inq_dimid(ncid, 'num_group', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=n)
    call handle_netcdf_error (status)
    allocate(this%group_id_list(n))
    status = nf90_inq_varid(ncid, 'group_ids', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, this%group_id_list)
    call handle_netcdf_error (status)

    status = nf90_inq_dimid(ncid, 'num_symmetry', dimid)
    call handle_netcdf_error (status)
    status = nf90_inquire_dimension(ncid, dimid, len=n)
    call handle_netcdf_error (status)
    allocate(symmetry(n))
    status = nf90_inq_varid(ncid, 'symmetry', varid)
    call handle_netcdf_error (status)
    status = nf90_get_var(ncid, varid, symmetry)
    call handle_netcdf_error (status)
    this%mirror  = (symmetry(1:3) == 1)
    this%rot_axis = symmetry(4)
    this%num_rot  = symmetry(5)
    deallocate(symmetry)

    if (present(has_vf)) has_vf = (NF90_NOERR == nf90_inq_varid(ncid, 'val', varid))

    status = nf90_close (ncid)
    call handle_netcdf_error (status)

    end if

    call bcast_encl (this)
    if (present(has_vf)) call scl_bcast (has_vf)

  end subroutine read_encl

  subroutine handle_netcdf_error (status)
    use netcdf
    integer, intent(in) :: status
    if (status == NF90_NOERR) return
    print *, trim(nf90_strerror(status))
    stop
  end subroutine handle_netcdf_error

end module re_encl_type
