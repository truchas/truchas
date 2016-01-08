!!
!! BASE_MESH_CLASS
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module base_mesh_class

  use kinds, only: r8
  use parallel_communication
  use index_partitioning
  use bitfield_type
  implicit none
  private

  type, abstract, public :: base_mesh
    integer :: nnode=0, nface=0, ncell=0

    !! Relationship to external numbering.
    integer, allocatable :: xnode(:)  ! external node number
    integer, allocatable :: xcell(:)  ! external cell number

    !! Cell set data
    integer, allocatable :: cell_set_id(:)
    integer, allocatable :: cell_set_mask(:)

    !! Face set data
    integer, allocatable :: face_set_id(:)
    type(bitfield), allocatable :: face_set_mask(:)

    !! Node set data
    integer, allocatable :: node_set_id(:)
    integer, allocatable :: node_set_mask(:)

    real(r8), allocatable :: x(:,:)
    real(r8), allocatable :: area(:)
    real(r8), allocatable :: volume(:)

    !! Partitioning and inter-process communication data.
    integer :: nnode_onP=0, nface_onP=0, ncell_onP=0
    type(ip_desc) :: node_ip, face_ip, cell_ip
  contains
    procedure :: get_face_set_ids
    procedure :: get_cell_set_bitmask
    procedure :: get_face_set_bitmask
    procedure :: get_global_x_array
    procedure :: get_global_volume_array
    procedure(wp), deferred :: write_profile
  end type base_mesh

  abstract interface
    subroutine wp (this)
      import base_mesh
      class(base_mesh), intent(in) :: this
    end subroutine
  end interface

contains

  !! Given a list of local face indices (on each process), this subroutine
  !! returns the list of IDs of face sets to which one or more of the faces
  !! belong.  This is a global procedure, returning the same result on all
  !! processes.

  subroutine get_face_set_ids (this, faces, setids)

    class(base_mesh), intent(in) :: this
    integer, intent(in) :: faces(:)
    integer, allocatable, intent(out) :: setids(:)

    integer :: j, n
    type(bitfield) :: bitmask

    bitmask = ZERO_BITFIELD
    do j = 1, size(faces)
      bitmask = ior(bitmask, this%face_set_mask(faces(j)))
    end do
    bitmask = ibclr(bitmask, pos=0) ! clear the boundary flag
    bitmask = global_ior(bitmask)

    !! Create the list of involved side set IDS.
    n = 0 ! count first to allocate
    do j = 1, size(this%face_set_id)
      if (btest(bitmask,j)) n = n + 1
    end do
    allocate(setids(n))
    n = 0 ! now store the data
    do j = 1, size(this%face_set_id)
      if (btest(bitmask,j)) then
        n = n + 1
        setids(n) = this%face_set_id(j)
      end if
    end do

  end subroutine get_face_set_ids

  !! Returns a scalar bit mask for use in bit operations with the cell_set_mask
  !! array component.  The corresponding bit is set for each cell set ID given
  !! in the array SETIDS.  STAT returns a non-zero value if an unknown cell set
  !! ID is specified, and the optional allocatable deferred-length character
  !! ERRMSG is assigned an explanatory message if present.

  subroutine get_cell_set_bitmask (this, setids, bitmask, stat, errmsg)
    use string_utilities, only: i_to_c
    class(base_mesh), intent(in) :: this
    integer, intent(in) :: setids(:)
    integer(kind(this%cell_set_mask)), intent(out) :: bitmask
    integer, intent(out) :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    integer :: i, j
    bitmask = 0
    do i = 1, size(setids)
      do j = size(this%cell_set_id), 1, -1
        if (setids(i) == this%cell_set_id(j)) exit
      end do
      if (j == 0) then
        stat = 1
        if (present(errmsg)) errmsg = 'unknown cell set ID: ' // i_to_c(setids(i))
        return
      end if
      bitmask = ibset(bitmask, j)
    end do
    stat = 0
  end subroutine get_cell_set_bitmask

  !! Returns a scalar bit mask for use in bit operations with the face_set_mask
  !! array component.  The corresponding bit is set for each face set ID given
  !! in the array SETIDS.  STAT returns a non-zero value if an unknown face set
  !! ID is specified, and the optional allocatable deferred-length character
  !! ERRMSG is assigned an explanatory message if present.

  subroutine get_face_set_bitmask (this, setids, bitmask, stat, errmsg)
    use string_utilities, only: i_to_c
    class(base_mesh), intent(in) :: this
    integer, intent(in) :: setids(:)
    type(bitfield), intent(out) :: bitmask
    integer, intent(out) :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    integer :: i, j
    bitmask = ZERO_BITFIELD
    do i = 1, size(setids)
      do j = size(this%face_set_id), 1, -1
        if (setids(i) == this%face_set_id(j)) exit
      end do
      if (j == 0) then
        stat = 1
        if (present(errmsg)) errmsg = 'unknown face set ID: ' // i_to_c(setids(i))
        return
      end if
      bitmask = ibset(bitmask, j)
    end do
    stat = 0
  end subroutine get_face_set_bitmask

  !! Returns the global node coordinate array for the mesh. The collated coord
  !! array is returned on the IO processor and a 0-sized array on all others.

  subroutine get_global_x_array (this, x)
    class(base_mesh), intent(in) :: this
    real(r8), allocatable, intent(out) :: x(:,:)
    allocate(x(size(this%x,1),merge(this%node_ip%global_size(),0,is_IOP)))
    call collate (x, this%x(:,:this%nnode_onP))
  end subroutine get_global_x_array

  !! Returns the global cell volume array for the mesh.  The collated volume
  !! array is returned on the IO processor and a 0-sized array on all others.

  subroutine get_global_volume_array (this, volume)
    class(base_mesh), intent(in) :: this
    real(r8), allocatable, intent(out) :: volume(:)
    allocate(volume(merge(this%cell_ip%global_size(),0,is_IOP)))
    call collate (volume, this%volume(:this%ncell_onP))
  end subroutine get_global_volume_array

end module base_mesh_class
