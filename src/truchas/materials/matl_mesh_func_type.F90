!! MATL_MESH_FUNC
!!
!! This module defines a derived type to represent the distribution of
!! materials across a mesh.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! September 2018. A F2003 rewrite of an earlier version.
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
!! An instance of the derived type MATL_MESH_FUNC describes how materials are
!! distributed across a mesh. Conceptually the mesh cells are partitioned into
!! regions, and a set of one or more materials associated with each region. A
!! cell may contain any mixture of the materials associated with the region to
!! which it belongs. Cells in a single-material region obviously contain only
!! that material. Multi-material regions include material volume fraction data
!! for each cell that describe the mixture of materials. Materials are
!! associated with a non-negative integer ID. ID 0 has special significance
!! as a "void" material in several procedures, but is otherwise no different.
!!
!! A MATL_MESH_FUNC object is defined incrementally using the following type
!! bound subroutines:
!!
!!  INIT(MESH) initializes the object to begin receiving the specification of
!!    mesh material regions. MESH is a BASE_MESH class object on which the data
!!    is defined. The object maintains an internal read-only reference to MESH.
!!
!!  DEFINE_REGION(SETIDS, MATIDS, STAT, ERRMSG) defines a material region that
!!    consists of the cells belonging to the cell sets with IDs specified by
!!    rank-1 integer array SETIDS. The materials that may be present in the
!!    region are specified by the rank-1 integer array MATIDS of their IDs.
!!    A nonzero value is assigned to the integer STAT if an error occurs, and
!!    an explanatory message assigned to the deferred length allocatable
!!    character argument ERRMSG.  The returned STAT and ERRMSG values are
!!    collective across all processes. Possible errors are referring to an
!!    unknown cell set ID, or attempting to define overlapping regions. This
!!    subroutine can be called multiple times after the initial call to INIT.
!!
!!  DEFINE_COMPLETE(STAT, ERRMSG) performs the final configuration of the
!!    object after all the desired calls to DEFINE_REGION have been made.
!!    A nonzero value is assigned to the integer STAT if an error occurs,
!!    and an explanatory message assigned to the deferred length allocatable
!!    character argument ERRMSG.  The returned STAT and ERRMSG values are
!!    collective across all processes. It is an error if the defined regions
!!    do not completely cover mesh.
!!
!! Once defined, the following type bound procedures may be used.
!!
!!  NUM_REG() returns the number of regions.
!!
!!  NUM_REG_MATL(N) returns the number of materials in region N.
!!
!!  REG_MATL(N) returns a pointer to a rank-1 integer array of the material IDs
!!    in region N. The size of the target array is NUM_REG_MATL(N).  The array
!!    is ordered.
!!
!!  REG_CELLS(N) returns a pointer to a rank-1 integer array containing the
!!    indices of the cells that comprise region N. The index list is ordered.
!!
!!  REG_VOL_FRAC(N) returns a pointer to a rank-2 real array containing the
!!    material volume fractions for region N. For single-material regions a
!!    null pointer is returned; in this case the volume fraction is 1 for all
!!    cells. If
!!
!!      MATID => REG_MATL(N), CELL => REG_CELLS(N), VF => REG_VOL_FRAC(N),
!!
!!    then VF(j,k) is the volume fraction of material MATID(k) in cell CELL(j).
!!    VF has shape [SIZE(CELL), SIZE(MATID)]. (NB: should we swap dimensions?)
!!    Application code is expected to define the values of the target array;
!!    the object merely manages the storage.
!!
!!  MESH_PTR() returns a pointer to the BASE_MESH class mesh on which the
!!    object is defined.
!!
!!  GET_ALL_MATL(MATIDS [,DROP_VOID]) returns the list of all material IDS that
!!    are present in any region in the allocatable integer array MATIDS. The
!!    list is ordered and there are no duplicates. If the optional argument
!!    DROP_VOID is present with value true, then the void material index (0) is
!!    not included in the list.
!!
!!  REG_HAS_MATL(N, MATID) returns true if region N contains material MATID.
!!
!!  HAS_MATL(MATID) returns true if some region contains material MATID.
!!
!!  GET_MATL_VOL_FRAC(MATID, VOL_FRAC) returns the volume fractions for
!!    material MATID for all mesh cells in the rank-1 real array VOL_FRAC.
!!    The volume fraction for cells in regions not containing the material is 0.
!!
!!  GET_VOID_VOL_FRAC(VOL_FRAC) returns the volume fraction of void (ID 0) for
!!    all mesh cells in the rank-1 real array VOL_FRAC.
!!

#include "f90_assert.fpp"

module matl_mesh_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use base_mesh_class
  use cell_group_builder_type
  use integer_set_type
  implicit none
  private

  type :: region
    integer,  allocatable :: mlist(:)   ! list of material IDs for the region
    integer,  allocatable :: index(:)   ! list of cell indices comprising the region
    real(r8), allocatable :: vfrac(:,:) ! volume fraction array
  end type region

  type, public :: matl_mesh_func
    private
    class(base_mesh), pointer :: mesh => null() ! reference only -- do not own
    type(region), allocatable :: reg(:)
    ! temporaries used during construction
    type(cell_group_builder), allocatable :: builder
    type(mset_list), pointer :: list => null()
  contains
    procedure :: init
    procedure :: define_region
    procedure :: define_complete
    procedure :: num_reg
    procedure :: num_reg_mat
    procedure :: reg_matl
    procedure :: reg_cell
    procedure :: reg_vol_frac
    procedure :: mesh_ptr
    procedure :: get_all_matl
    procedure :: reg_has_matl
    procedure :: has_matl
    procedure :: get_matl_vol_frac
    procedure :: get_void_vol_frac
    final :: matl_mesh_func_delete
  end type

  type :: mset_list
    type(integer_set) :: mset
    type(mset_list), pointer :: next => null()
  end type mset_list

contains

  !! Final subroutine for MATL_MESH_FUNC objects
  subroutine matl_mesh_func_delete(this)
    type(matl_mesh_func), intent(inout) :: this
    type(mset_list), pointer :: first
    do while (associated(this%list))
      first => this%list
      this%list => first%next
      deallocate(first)
    end do
  end subroutine matl_mesh_func_delete

  subroutine init(this, mesh)
    class(matl_mesh_func), intent(out) :: this
    class(base_mesh), intent(in), target :: mesh
    this%mesh => mesh
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init

  subroutine define_region(this, setids, matids, stat, errmsg)

    class(matl_mesh_func), intent(inout) :: this
    integer,   intent(in) :: setids(:), matids(:)
    integer,   intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(mset_list), pointer :: first

    call this%builder%add_cell_group(setids, stat, errmsg)
    if (stat /= 0) return

    INSIST(size(matids) > 0)
    INSIST(all(matids >= 0))

    !! Prepend the material ID set to LIST
    allocate(first)
    call first%mset%add(matids)
    first%next => this%list
    this%list => first

  end subroutine define_region

  subroutine define_complete(this, stat, errmsg)

    use string_utilities, only: i_to_c

    class(matl_mesh_func), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n, ncell, nmat, nreg
    integer, allocatable :: index(:), xgroup(:)
    type(mset_list), pointer :: first

    call this%builder%get_cell_groups(nreg, xgroup, index)

    !! Verify every cell has been associated with a region.
    n = this%mesh%ncell - size(index)
    if (n /= 0) then
      stat = 1
      errmsg = i_to_c(n) // ' cells not associated with a region'
      return
    end if

    !! Initialize %REG structure-array component using the data accumulated
    !! by the preceding calls to DEFINE_REGION. Note that %LIST is used as a
    !! stack (LIFO) and so we must fill the %REG array from the end as we pop
    !! data off the top of %LIST.
    allocate(this%reg(nreg))
    do n = nreg, 1, -1
      this%reg(n)%mlist = this%list%mset
      this%reg(n)%index = index(xgroup(n):xgroup(n+1)-1)
      nmat = size(this%reg(n)%mlist)
      ncell = size(this%reg(n)%index)
      if (nmat > 1) allocate(this%reg(n)%vfrac(ncell,nmat))
      !! Unlink the first element in the list of material sets
      first => this%list
      this%list => first%next
      deallocate(first)
    end do

    stat = 0

  end subroutine define_complete

  pure integer function num_reg(this)
    class(matl_mesh_func), intent(in) :: this
    num_reg = size(this%reg)
  end function num_reg

  pure integer function num_reg_mat(this, n)
    class(matl_mesh_func), intent(in) :: this
    integer, intent(in) :: n
    num_reg_mat = size(this%reg(n)%mlist)
  end function num_reg_mat

  function reg_matl(this, n) result(mlist)
    class(matl_mesh_func), intent(in), target :: this
    integer, intent(in) :: n
    integer, pointer :: mlist(:)
    mlist => this%reg(n)%mlist
  end function reg_matl

  pure logical function reg_has_matl(this, n, matid)
    class(matl_mesh_func), intent(in) :: this
    integer, intent(in) :: n, matid
    integer :: j
    reg_has_matl = .true.
    do j = 1, size(this%reg(n)%mlist)
      if (matid == this%reg(n)%mlist(j)) return
      if (matid < this%reg(n)%mlist(j)) exit ! list is sorted
    end do
    reg_has_matl = .false.
  end function reg_has_matl

  pure logical function has_matl(this, matid)
    class(matl_mesh_func), intent(in) :: this
    integer, intent(in) :: matid
    integer :: n
    has_matl = .true.
    do n = 1, size(this%reg)
      if (reg_has_matl(this, n, matid)) return
    end do
    has_matl = .false.
  end function has_matl

  subroutine get_all_matl(this, matids, drop_void)
    class(matl_mesh_func), intent(in) :: this
    integer, allocatable :: matids(:)
    logical, optional :: drop_void
    integer :: n
    type(integer_set) :: set
    do n = 1, size(this%reg)
      call set%add(this%reg(n)%mlist)
    end do
    if (present(drop_void)) then
      if (drop_void) call set%remove(0)
    end if
    matids = set
  end subroutine get_all_matl

  function reg_cell(this, n) result(list)
    class(matl_mesh_func), intent(in), target :: this
    integer, intent(in) :: n
    integer, pointer :: list(:)
    list => this%reg(n)%index
  end function reg_cell

  function reg_vol_frac(this, n) result(vfrac)
    class(matl_mesh_func), intent(in), target :: this
    integer, intent(in) :: n
    real(r8), pointer :: vfrac(:,:)
    if (allocated(this%reg(n)%vfrac)) then
      vfrac => this%reg(n)%vfrac
    else  ! single material region
      vfrac => null()
    end if
  end function reg_vol_frac

  function mesh_ptr(this)
    class(matl_mesh_func), intent(in) :: this
    class(base_mesh), pointer :: mesh_ptr
    mesh_ptr => this%mesh
  end function mesh_ptr

  subroutine get_void_vol_frac(this, vol_frac)
    class(matl_mesh_func), intent(in) :: this
    real(r8), intent(out) :: vol_frac(:)
    call get_matl_vol_frac(this, 0, vol_frac)
  end subroutine get_void_vol_frac

  subroutine get_matl_vol_frac(this, matid, vol_frac)

    class(matl_mesh_func), intent(in) :: this
    integer, intent(in) :: matid
    real(r8), intent(out) :: vol_frac(:)

    integer :: n, i

    ASSERT(size(vol_frac) == this%mesh%ncell)

    vol_frac = 0.0_r8
    do n = 1, size(this%reg)
      associate (mlist => this%reg(n)%mlist)
        do i = 1, size(mlist)
          if (mlist(i) == matid) then
            if (allocated(this%reg(n)%vfrac)) then
              vol_frac(this%reg(n)%index) = this%reg(n)%vfrac(:,i)
            else  ! single material region
              vol_frac(this%reg(n)%index) = 1.0_r8
            end if
            exit
          end if
        end do
      end associate
    end do

  end subroutine get_matl_vol_frac

end module matl_mesh_func_type
