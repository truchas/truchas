!!
!! T2D_FLOW_BC_TYPE
!!
!! This module defines T2D_FLOW_BC, the boundary-condition data used by
!! two-dimensional flow. It owns old-style sparse boundary functions created
!! by T2D_FLOW_BC_FACTORY and selects a single pressure reference face when
!! all pressure boundaries are homogeneous Neumann conditions.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "t2d_assert.inc"

module t2d_flow_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use t2d_unstr_mesh_type
  use parameter_list_type
  use bndry_func1_class
  use bndry_vfunc_class
  use t2d_flow_bc_factory_type
  use flow_domain_types
  use simulation_environment_type
  use parallel_communication
  implicit none
  private

  type, public :: t2d_flow_inflow_material
    character(:), allocatable :: name
    integer, allocatable :: face(:)
  end type

  type, public :: t2d_flow_bc
    type(t2d_unstr_mesh), pointer :: mesh => null()  ! unowned reference
    class(bndry_func1), allocatable :: pressure_dirichlet
    class(bndry_func1), allocatable :: pressure_correction_dirichlet
    class(bndry_func1), allocatable :: pressure_neumann
    class(bndry_func1), allocatable :: velocity_zero_normal
    class(bndry_vfunc), allocatable :: velocity_dirichlet
    type(t2d_flow_inflow_material), allocatable :: inflow_material(:)
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_initial
    procedure :: check_velocity_flux
    procedure :: pressure_pin_face
  end type

contains

  subroutine init(this, env, mesh, params, stat, errmsg)
    class(t2d_flow_bc), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(t2d_unstr_mesh), target, intent(in) :: mesh
    type(parameter_list), target, intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(t2d_flow_bc_factory) :: factory
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
    call read_inflow_material(this, mesh, params, stat, errmsg)
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


  !! Extract static material-inflow data from velocity and pressure boundary
  !! conditions.  The tracker consumes the resolved local face lists after
  !! the flow material mapping has established its reduced slot ordering.
  subroutine read_inflow_material(this, mesh, params, stat, errmsg)
    use bitfield_type, only: bitfield, ZERO_BITFIELD, iand, operator(/=)
    use string_utilities, only: lower_case

    class(t2d_flow_bc), intent(inout) :: this
    type(t2d_unstr_mesh), intent(in) :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: iter
    type(parameter_list), pointer :: plist
    type(bitfield) :: mask
    integer, allocatable :: setids(:)
    character(:), allocatable :: bc_type
    integer :: f, n

    n = 0
    iter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.iter%at_end())
      plist => iter%sublist()
      if (plist%is_parameter('inflow-material')) n = n + 1
      call iter%next()
    end do
    if (n == 0) then
      stat = 0
      errmsg = ''
      return
    end if
    allocate(this%inflow_material(n))

    n = 0
    iter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.iter%at_end())
      plist => iter%sublist()
      if (.not.plist%is_parameter('inflow-material')) then
        call iter%next()
        cycle
      end if
      call plist%get('type', bc_type, stat, errmsg)
      if (stat /= 0) exit
      if (lower_case(bc_type) /= 'velocity' .and. lower_case(bc_type) /= 'pressure') then
        stat = 1
        errmsg = 'T2D_FLOW_BC[' // iter%name() // ']: "inflow-material" requires a velocity or pressure boundary'
        exit
      end if
      n = n + 1
      call plist%get('inflow-material', this%inflow_material(n)%name, stat, errmsg)
      if (stat /= 0) exit
      call plist%get('face-set-ids', setids, stat, errmsg)
      if (stat /= 0) exit
      call mesh%get_face_set_bitmask(setids, mask, stat, errmsg)
      if (stat /= 0) exit
      this%inflow_material(n)%face = pack([(f, f=1,mesh%nface_onP)], &
          iand(mesh%face_set_mask(:mesh%nface_onP), mask) /= ZERO_BITFIELD)
      call iter%next()
    end do
    if (stat /= 0) errmsg = 'T2D_FLOW_BC[' // iter%name() // ']: ' // errmsg
  end subroutine


  !! Check the compatibility condition for a closed incompressible domain.
  !! If a pressure Dirichlet boundary is present, its normal velocity is not
  !! prescribed and the pressure solve supplies the compensating flux.
  subroutine check_velocity_flux(this, stat, errmsg)
    class(t2d_flow_bc), intent(in) :: this
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
    class(t2d_flow_bc), intent(inout) :: this
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
    class(t2d_flow_bc), intent(inout) :: this
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
  !! indicates that a pressure Dirichlet condition or a fluid/VOID interface
  !! already supplies a reference.
  function pressure_pin_face(this, face_t) result(face)
    class(t2d_flow_bc), intent(in) :: this
    integer, optional, intent(in) :: face_t(:)
    integer :: face

    integer :: i, pin_pe
    logical :: is_candidate, candidate(nPE)

    face = 0
    if (global_any(size(this%pressure_dirichlet%index) > 0)) return
    if (present(face_t)) then
      if (global_any(face_t == void_t)) return
    end if
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

    class(t2d_flow_bc), intent(inout) :: this
    type(t2d_unstr_mesh), intent(in) :: mesh

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

end module t2d_flow_bc_type
