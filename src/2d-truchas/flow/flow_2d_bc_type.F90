!!
!! FLOW_2D_BC_TYPE
!!
!! This module defines FLOW_2D_BC, the boundary-condition data used by
!! two-dimensional flow. It owns old-style sparse boundary functions created
!! by FLOW_2D_BC_FACTORY and selects a single pressure reference face when
!! all pressure boundaries are homogeneous Neumann conditions.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use parameter_list_type
  use bndry_func1_class
  use bndry_vfunc_class
  use flow_2d_bc_factory_type
  use parallel_communication
  implicit none
  private

  type, public :: flow_2d_bc
    class(bndry_func1), allocatable :: pressure_dirichlet
    class(bndry_func1), allocatable :: pressure_correction_dirichlet
    class(bndry_func1), allocatable :: pressure_neumann
    class(bndry_func1), allocatable :: velocity_zero_normal
    class(bndry_vfunc), allocatable :: velocity_dirichlet
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_initial
    procedure :: pressure_pin_face
  end type

contains

  subroutine init(this, mesh, params, stat, errmsg)
    class(flow_2d_bc), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(flow_2d_bc_factory) :: factory

    call factory%init(mesh, params)
    call factory%alloc_dir_vel_bc(this%velocity_dirichlet, stat, errmsg)
    if (stat /= 0) return
    call factory%alloc_zero_vn_bc(this%velocity_zero_normal, stat, errmsg)
    if (stat /= 0) return
    call factory%alloc_dir_prs_bc(this%pressure_dirichlet, stat, errmsg)
    if (stat /= 0) return
    call factory%alloc_dir_prs_bc(this%pressure_correction_dirichlet, stat, errmsg)
    if (stat /= 0) return
    call factory%alloc_neu_prs_bc(this%pressure_neumann, stat, errmsg)
    if (stat /= 0) return
    call apply_default(this, mesh)
  end subroutine


  subroutine compute(this, time, dt)
    class(flow_2d_bc), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), optional, intent(in) :: dt

    real(r8) :: dt_

    dt_ = 0.0_r8
    if (present(dt)) dt_ = dt
    call this%velocity_dirichlet%compute(time)
    call this%velocity_zero_normal%compute(time)
    call this%pressure_dirichlet%compute(time)
    call this%pressure_neumann%compute(time)
    call this%pressure_correction_dirichlet%compute(time + dt_)
    this%pressure_correction_dirichlet%value = this%pressure_correction_dirichlet%value - &
        this%pressure_dirichlet%value
  end subroutine


  subroutine compute_initial(this, time)
    class(flow_2d_bc), intent(inout) :: this
    real(r8), intent(in) :: time

    call this%velocity_dirichlet%compute(time)
    call this%velocity_zero_normal%compute(time)
    call this%pressure_dirichlet%compute(time)
    call this%pressure_neumann%compute(time)
    this%pressure_correction_dirichlet%value = 0.0_r8
  end subroutine


  !! Return the local boundary face at which a zero pressure reference is to
  !! be imposed. All ranks must call this collective function. A value of zero
  !! indicates that a pressure Dirichlet condition already supplies a reference.
  function pressure_pin_face(this) result(face)
    class(flow_2d_bc), intent(in) :: this
    integer :: face

    integer :: pin_pe
    logical :: is_candidate, candidate(nPE)

    face = 0
    if (global_any(size(this%pressure_dirichlet%index) > 0)) return
    is_candidate = size(this%pressure_neumann%index) > 0
    call gather(is_candidate, candidate)
    if (is_IOP) pin_pe = findloc(candidate, .true., dim=1)
    call broadcast(pin_pe)
    INSIST(pin_pe > 0 .and. pin_pe <= nPE)
    if (this_PE == pin_pe) face = this%pressure_neumann%index(1)
  end function


  subroutine apply_default(this, mesh)
    use bndry_face_func_type
    use scalar_func_class
    use scalar_func_factories, only: alloc_const_scalar_func

    class(flow_2d_bc), intent(inout) :: this
    type(unstr_2d_mesh), intent(in) :: mesh

    class(scalar_func), allocatable :: func
    integer, allocatable :: faces(:)
    integer :: f, nface

    allocate(faces(mesh%nface_onP))
    nface = 0
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) /= 0) cycle
      if (any(this%pressure_dirichlet%index == f)) cycle
      if (any(this%pressure_neumann%index == f)) cycle
      nface = nface + 1
      faces(nface) = f
    end do
    if (nface > 0) then
      call alloc_const_scalar_func(func, 0.0_r8)
      select type (bndry => this%pressure_neumann)
      type is (bndry_face_func)
        call bndry%add_face_list(func, faces(:nface))
      class default
        ASSERT(.false.)
      end select
    end if

    nface = 0
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) /= 0) cycle
      if (any(this%velocity_dirichlet%index == f)) cycle
      if (any(this%velocity_zero_normal%index == f)) cycle
      ! A pressure Dirichlet boundary is an open boundary: its normal velocity
      ! is determined by the pressure solve, not prescribed as free slip.
      if (any(this%pressure_dirichlet%index == f)) cycle
      nface = nface + 1
      faces(nface) = f
    end do
    if (nface == 0) return
    call alloc_const_scalar_func(func, 0.0_r8)
    select type (bndry => this%velocity_zero_normal)
    type is (bndry_face_func)
      call bndry%add_face_list(func, faces(:nface))
    class default
      ASSERT(.false.)
    end select
  end subroutine

end module flow_2d_bc_type
