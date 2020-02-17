!!
!! RE_ENCL_TYPE
!!
!! This module provides a derived type for describing the enclosing surface of
!! a radiation enclosure and methods that operate on instances of this type.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 3 Apr 2008; updated February 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The ENCL type describes the meshed surface of a radiation enclosure. It has
!! the following type bound subroutines.
!!
!!  BCAST() broadcasts the contents of the object on process rank 1 to all
!!    other process ranks, replicating the enclosure on all processes.
!!
!!  READ(PATH [,HAS_VF]) reads the radiation enclosure dataset PATH and
!!    initializes the object using the data.  The optional logical argument
!!    HAS_VF is assigned the value true if the dataset contains view factor
!!    data.  This is a collective procedure. Input occurs on process rank 1
!!    and the enclosure data is then replicated on all other processes.
!!
!!  WRITE(PATH) writes the enclosure object to the radiation enclosure dataset
!!    PATH. Any previous file is overwritten.  This is a collective procedure,
!!    but only for the purposes of error handling.  Output occurs on process
!!    rank 1 using its data; the object is ignored on all other processes.
!!
!! The ENCL_LIST type extends the ENCL type with additional data that defines a
!! list of related enclosure surfaces.  The surfaces have the same logical mesh
!! topology and symmetries, differing only in a displacement of certain surface
!! components.  Objects of this type are instantiated using the INIT_ENCL_LIST
!! procedure from the RE_EXODUS_ENCL module.  The ENCL components of the object
!! describe the current enclosure of the list, and the object is positioned at
!! the first enclosure in the list immediately after initialization.  The type
!! has the following type bound methods.
!!
!!  NUM_ENCL() returns the number of enclosures in the list.
!!
!!  NEXT_ENCL() advances the object to the next enclosure in the list.
!!
!!  THIS_LABEL() returns a unique character string label for the current
!!    enclosure which can be used in forming a unique output file name. The
!!    labels will all consist of the same number of nonblank characters. In
!!    some cases this is just the numeric position in the list (left padded
!!    with 0s).  In other cases, part of the hash of the displacement.
!!

#include "f90_assert.fpp"

module re_encl_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scl
  use rad_encl_file_type
  implicit none
  private

  type, public :: encl
    character(:), allocatable :: name
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
  contains
    procedure :: bcast => encl_bcast
    procedure :: write => encl_write
    procedure :: read  => encl_read
  end type encl

  type, extends(encl), public :: encl_list
    integer :: n = 0, index = 0
    real(r8), allocatable :: x0(:,:)
    logical,  allocatable :: mask(:)
    real(r8), allocatable :: dx(:,:)
    character(:), allocatable :: label(:)
  contains
    procedure :: next_encl
    procedure :: this_label
    procedure :: num_encl
    procedure :: bcast => encl_list_bcast
  end type encl_list

contains

  !! Broadcast the rank 1 ENCL_LIST components to the other ranks.
  subroutine encl_list_bcast(this)
    class(encl_list), intent(inout) :: this
    call this%encl%bcast
    call scl_bcast(this%n)
    call scl_bcast(this%index)
    call scl_bcast_alloc(this%mask)
    call scl_bcast_alloc(this%dx)
    call scl_bcast_alloc(this%x0)
    call scl_bcast_alloc(this%label)
  end subroutine encl_list_bcast

  !! Advance to the next enclosure.
  subroutine next_encl(this)
    class(encl_list), intent(inout) :: this
    integer :: j
    if (this%index >= this%n) return
    this%index = this%index + 1
    do j = 1, this%nnode
      if (this%mask(j)) this%x(:,j) = this%x0(:,j) + this%dx(:,this%index)
    end do
  end subroutine next_encl

  !! Return the unique string label for the current enclosure.
  function this_label(this) result(label)
    class(encl_list), intent(in) :: this
    character(:), allocatable :: label
    label = this%label(this%index)
  end function this_label

  !! Return the number of enclosures in the list.
  integer function num_encl(this)
    class(encl_list), intent(in) :: this
    num_encl = this%n
  end function num_encl

  !! Broadcast the rank 1 ENCL components to the other ranks.
  subroutine encl_bcast(this)
    class(encl), intent(inout) :: this
    call scl_bcast_alloc(this%name)
    call scl_bcast(this%nnode)
    call scl_bcast(this%nface)
    call scl_bcast_alloc(this%xface)
    call scl_bcast_alloc(this%fnode)
    call scl_bcast_alloc(this%x)
    call scl_bcast_alloc(this%group_id_list)
    call scl_bcast_alloc(this%gnum)
    call scl_bcast_alloc(this%src_elem)
    call scl_bcast_alloc(this%src_side)
    call scl_bcast(this%mirror)
    call scl_bcast(this%rot_axis)
    call scl_bcast(this%num_rot)
  end subroutine encl_bcast

  !! Write a radiation enclosure file PATH with the enclosure data from rank 1.
  subroutine encl_write(this, path)

    class(encl), intent(in) :: this
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

  end subroutine encl_write

  !! Initialize the object with data read from the radiation enclosure data
  !! set PATH.  The file is read on rank 1 and the data broadcast to other
  !! ranks.  HAS_VF returns true if the file contains view factor data.
  subroutine encl_read(this, path, has_vf)

    class(encl), intent(out) :: this
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

    call this%bcast
    if (present(has_vf)) call scl_bcast(has_vf)

  end subroutine encl_read

end module re_encl_type
