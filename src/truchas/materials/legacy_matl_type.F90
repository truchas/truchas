!!
!! LEGACY_MATL_TYPE
!!
!! This is a modern reimplementation of the legacy MATL data structure. It is
!! faithful to the original structure in order to serve as a baseline reference
!! for new alternative designs of storing and managing the material cell volume
!! fractions. The modernization included using allocatable components instead
!! of pointers, and wrapping the original MATL structure array in a new derived
!! type LEGACY_MATL with type bound procedures. Truchas will use a singleton of
!! this type. In addition, the original collection of MATL procedures was pruned
!! to just those currently being used.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! July 2022
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module legacy_matl_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type :: material
    integer  :: id = 0
    real(r8) :: vof = 0.0_r8
  end type

  type :: matl_slot
    type(material), allocatable :: cell(:)
  end type

  type, public :: legacy_matl
    private
    integer :: nmat ! protected
    integer :: ncell
    type(matl_slot), allocatable :: slot(:)
  contains
    procedure :: init
    procedure :: gather_vof
    procedure :: get_vof
    procedure :: get_cell_vof
    procedure :: set_vof
    procedure, private :: resize
  end type

contains

  subroutine init(this, nmat, ncell)
    class(legacy_matl), intent(out) :: this
    integer, intent(in) :: nmat, ncell
    this%nmat = nmat
    this%ncell = ncell
    allocate(this%slot(1))
    allocate(this%slot(1)%cell(ncell))
  end subroutine

  !! This returns a rank-1 array VOF of volume fractions of the given material
  !! index M for all cells.

  subroutine gather_vof(this, m, vof)
    class(legacy_matl), intent(in) :: this
    integer, intent(in) :: m
    real(r8), intent(inout) :: vof(:)
    integer :: j, s
    ASSERT(size(vof) == this%ncell)
    vof = 0.0_r8
    do s = 1, size(this%slot)
      associate (cell => this%slot(s)%cell)
        do j = 1, size(vof)
          if (cell(j)%id == m) vof(j) = cell(j)%vof
        end do
      end associate
    end do
  end subroutine

  !! This returns a rank-2 array VOF of all material volume fractions
  !! for all cells -- the uncompressed form of the internal data structure.

  subroutine get_vof(this, vof)
    class(legacy_matl), intent(in) :: this
    real(r8), intent(inout) :: vof(:,:)
    integer :: j, s
    ASSERT(size(vof,dim=1) == this%nmat)
    ASSERT(size(vof,dim=2) == this%ncell)
    vof = 0.0_r8
    do s = 1, size(this%slot)
      associate (cell => this%slot(s)%cell)
        do j = 1, size(vof,dim=2)
          if (cell(j)%id > 0) vof(cell(j)%id,j) = cell(j)%vof
        end do
      end associate
    end do
  end subroutine

  !! This returns a rank-1 array VOF of all material volume fractions
  !! for the given cell index J.

  subroutine get_cell_vof(this, j, vof)
    class(legacy_matl), intent(in) :: this
    integer, intent(in) :: j
    real(r8), intent(out) :: vof(:)
    integer :: s
    ASSERT(size(vof) == this%nmat)
    vof = 0.0_r8
    do s = 1, size(this%slot)
      associate (cell_j => this%slot(s)%cell(j))
        if (cell_j%id > 0) vof(cell_j%id) = cell_j%vof
      end associate
    end do
  end subroutine

  !! This defines the internal compressed storage structure given the
  !! uncompressed rank-2 array of volume fractions for all materials.

  subroutine set_vof(this, vof)
    use parallel_communication, only: global_maxval
    class(legacy_matl), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:)
    integer :: j, m , s, nslot
    ASSERT(size(vof,dim=1) <= this%nmat)
    ASSERT(size(vof,dim=2) == this%ncell)
    nslot = global_maxval(count(vof > 0, dim=1))
    call this%resize(nslot)
    do j = 1, this%ncell
      s = 1
      do m = 1, size(vof,dim=1)
        if (vof(m,j) <= 0.0_r8) cycle
        associate (cell_j => this%slot(s)%cell(j))
          cell_j%id  = m
          cell_j%vof = vof(m,j)
        end associate
        s = s + 1
      end do
      do s = s, size(this%slot)  ! fill remainder with null values
        associate (cell_j => this%slot(s)%cell(j))
          cell_j%id  = 0
          cell_j%vof = 0.0_r8
        end associate
      end do
    end do

  end subroutine

  !! This resizes the internal storage to the specified number of slots.
  !! Allocated slots and their data are preserved where possible. Newly
  !! Allocated slots are default initialized with null data.

  subroutine resize(this, nslot)
    class(legacy_matl), intent(inout) :: this
    integer, intent(in) :: nslot
    integer :: s
    type(matl_slot), allocatable :: old_slot(:)
    if (nslot /= size(this%slot)) then
      call move_alloc(this%slot, old_slot)
      allocate(this%slot(nslot))
      do s = 1, min(nslot, size(old_slot))
        call move_alloc(old_slot(s)%cell, this%slot(s)%cell)
      end do
      do s = size(old_slot)+1, nslot
        allocate(this%slot(s)%cell(this%ncell))
      end do
    end if
  end subroutine

end module legacy_matl_type
