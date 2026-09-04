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
  use flow_domain_types
  use simulation_environment_type
  use parallel_communication
  implicit none
  private

  type, public :: flow_2d_bc
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    class(bndry_func1), allocatable :: pressure_dirichlet
    class(bndry_func1), allocatable :: pressure_correction_dirichlet
    class(bndry_func1), allocatable :: pressure_neumann
    class(bndry_func1), allocatable :: velocity_zero_normal
    class(bndry_vfunc), allocatable :: velocity_dirichlet
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_initial
    procedure :: check_velocity_flux
    procedure :: pressure_pin_face
  end type

contains

  subroutine init(this, env, mesh, params, stat, errmsg)
    class(flow_2d_bc), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(flow_2d_bc_factory) :: factory
    integer :: i
    logical :: overlap

    this%mesh => mesh
    call factory%init(mesh, params)
    call factory%alloc_dir_vel_bc(this%velocity_dirichlet, env, stat, errmsg, report=.true.)
    if (stat /= 0) return
    call factory%alloc_zero_vn_bc(this%velocity_zero_normal, env, stat, errmsg, report=.true.)
    if (stat /= 0) return
    call factory%alloc_dir_prs_bc(this%pressure_dirichlet, env, stat, errmsg, report=.true.)
    if (stat /= 0) return
    call factory%alloc_dir_prs_bc(this%pressure_correction_dirichlet, env, stat, errmsg)
    if (stat /= 0) return
    call factory%alloc_neu_prs_bc(this%pressure_neumann, env, stat, errmsg)
    if (stat /= 0) return
    call apply_default(this, mesh)

    overlap = .false.
    do i = 1, size(this%pressure_dirichlet%index)
      if (any(this%velocity_zero_normal%index == this%pressure_dirichlet%index(i)) .or. &
          any(this%velocity_dirichlet%index == this%pressure_dirichlet%index(i))) then
        overlap = .true.
        exit
      end if
    end do
    if (global_any(overlap)) then
      stat = 1
      errmsg = 'pressure Dirichlet boundary overlaps a velocity boundary condition'
    end if
  end subroutine


  !! Check the compatibility condition for a closed incompressible domain.
  !! If a pressure Dirichlet boundary is present, its normal velocity is not
  !! prescribed and the pressure solve supplies the compensating flux.
  subroutine check_velocity_flux(this, stat, errmsg)
    class(flow_2d_bc), intent(in) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: f, i
    real(r8) :: local_flux, flux, scale, tolerance

    stat = 0
    if (global_any(size(this%pressure_dirichlet%index) > 0)) return
    local_flux = 0.0_r8
    scale = 0.0_r8
    do f = 1, this%mesh%nface_onP
      if (this%mesh%fcell(2,f) /= 0) cycle
      i = findloc(this%velocity_dirichlet%index, f, dim=1)
      if (i == 0) cycle  ! The default or explicit free-slip value is zero.
      flux = this%mesh%area(f)*dot_product(this%mesh%unit_normal(:,f), &
          this%velocity_dirichlet%value(:,i))
      local_flux = local_flux + flux
      scale = scale + abs(flux)
    end do
    flux = global_sum(local_flux)
    scale = global_sum(scale)
    tolerance = 100.0_r8*epsilon(1.0_r8)*max(1.0_r8, scale)
    if (abs(flux) > tolerance) then
      stat = 1
      errmsg = 'incompatible prescribed velocity boundary flux'
    end if
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
  !! be imposed. If FACE_T is present, only currently regular faces are
  !! considered. All ranks must call this collective function. A value of zero
  !! indicates that a pressure Dirichlet condition already supplies a reference.
  function pressure_pin_face(this, face_t) result(face)
    class(flow_2d_bc), intent(in) :: this
    integer, optional, intent(in) :: face_t(:)
    integer :: face

    integer :: i, pin_pe
    logical :: is_candidate, candidate(nPE)

    face = 0
    if (global_any(size(this%pressure_dirichlet%index) > 0)) return
    pin_pe = 0
    if (present(face_t)) then
      is_candidate = .false.
      do i = 1, size(this%pressure_neumann%index)
        if (this%pressure_neumann%index(i) <= size(face_t)) &
            is_candidate = is_candidate .or. face_t(this%pressure_neumann%index(i)) == regular_t
      end do
    else
      is_candidate = size(this%pressure_neumann%index) > 0
    end if
    call gather(is_candidate, candidate)
    if (is_IOP) pin_pe = findloc(candidate, .true., dim=1)
    call broadcast(pin_pe)
    if (pin_pe == 0) return
    if (this_PE == pin_pe) then
      if (present(face_t)) then
        do i = 1, size(this%pressure_neumann%index)
          if (this%pressure_neumann%index(i) <= size(face_t)) then
            if (face_t(this%pressure_neumann%index(i)) == regular_t) then
              face = this%pressure_neumann%index(i)
              exit
            end if
          end if
        end do
      else
        face = this%pressure_neumann%index(1)
      end if
    end if
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
