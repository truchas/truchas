!!
!! PEC_BNDRY_FUNC_TYPE
!!
!! This module defines an implementation of the base class BNDRY_FUNC1 that
!! implements the perfect electric conductor (PEC) EM boundary condition for
!! the edge-based electric field unknowns on a SIMPL_MESH type mesh. Since
!! the boundary values are all 0, this is relatively trivial.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!

#include "f90_assert.fpp"

module pec_bndry_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_func1_class
  use simpl_mesh_type
  use bndry_edge_group_builder_type
  implicit none
  private

  type, extends(bndry_func1), public :: pec_bndry_func
    ! temporaries used during construction
    type(bndry_edge_group_builder), allocatable :: builder
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
  end type

contains

  subroutine init(this, mesh)
    class(pec_bndry_func), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine

  subroutine add(this, setids, stat, errmsg)
    class(pec_bndry_func), intent(inout) :: this
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_face_group(setids, stat, errmsg)
  end subroutine

  subroutine add_complete(this, stat)
    class(pec_bndry_func), intent(inout) :: this
    integer, intent(out) :: stat
    integer :: ngroup
    integer, allocatable :: xgroup(:)
    ASSERT(allocated(this%builder))
    call this%builder%get_edge_groups(ngroup, xgroup, this%index, stat)
    deallocate(this%builder)
    allocate(this%value(size(this%index)))
    this%value = 0.0_r8
  end subroutine

  subroutine compute(this, t)
    class(pec_bndry_func), intent(inout) :: this
    real(r8), intent(in) :: t
    ! All values are 0, set by add_complete.
  end subroutine

end module pec_bndry_func_type
