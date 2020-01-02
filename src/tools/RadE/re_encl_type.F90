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
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
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
!!  FACE_AREA (THIS) returns the areas of all faces in the enclosure as a
!!    rank-1 real array.  This procedure must be called after the ecnlosure
!!    THIS is initialized with READ_ENCL. This procedure performs no
!!    communication.
!!

#include "f90_assert.fpp"

module re_encl_type

  use kinds, only: r8
  use scl
  use rad_encl_file_type
  implicit none
  private

  public :: bcast_encl
  public :: read_encl, write_encl

  integer, parameter :: MAX_NAME_LEN = 32

  type, public :: encl
    character(len=MAX_NAME_LEN) :: name
    !! Bare surface representation
    integer :: nnode=0, nface=0
    integer,  allocatable :: xface(:), fnode(:)
    real(r8), allocatable :: x(:,:)
    !! Face group attributes
    integer,  allocatable :: gnum(:), group_id_list(:)
    !! Connection to the source mesh
    integer, allocatable :: src_elem(:), src_side(:)
    !! Surface symmetries
    logical :: mirror(3)
    integer :: rot_axis, num_rot
  end type

contains

  subroutine bcast_encl (this)

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

  subroutine write_encl(this, path)

    type(encl), intent(in) :: this
    character(*), intent(in) :: path

    type(rad_encl_file) :: file

    if (scl_rank() == 1) then
      call file%create(path)
      call file%put_surface(this%xface, this%fnode, this%x)
      call file%put_group_info(this%gnum, this%group_id_list)
      call file%put_symmetry_info(this%mirror, this%rot_axis, this%num_rot)
      call file%put_source_info(this%src_elem, this%src_side)
      call file%close
    end if

  end subroutine write_encl

  subroutine read_encl(this, path, has_vf)

    type(encl), intent(out) :: this
    character(*), intent(in) :: path
    logical, intent(out), optional :: has_vf

    type(rad_encl_file) :: file
    integer, allocatable :: fsize(:)
    integer :: j

    if (scl_rank() == 1) then
      call file%open_ro(path)
      call file%get_surface(fsize, this%fnode, this%x)
      this%nnode = size(this%x,dim=2)
      this%nface = size(fsize)
      allocate(this%xface(this%nface+1))
      this%xface(1) = 1
      do j = 1, size(fsize)
        this%xface(j+1) = this%xface(j) + fsize(j)
      end do
      deallocate(fsize)
      call file%get_group_info(this%gnum, this%group_id_list)
      call file%get_symmetry_info(this%mirror, this%rot_axis, this%num_rot)
      call file%get_source_info(this%src_elem, this%src_side)
      INSIST(size(this%gnum) == this%nface)
      INSIST(size(this%src_elem) == this%nface)
      INSIST(size(this%src_side) == this%nface)
      if (present(has_vf)) has_vf = file%has_vf_data()
      call file%close
    end if

    call bcast_encl(this)
    if (present(has_vf)) call scl_bcast(has_vf)

  end subroutine read_encl

end module re_encl_type
