!!
!! VISCOPLASTIC_SOLVER_TYPE
!!
!! This type owns the viscoplastic model and controls the choice of integrator.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! May 2022
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module viscoplastic_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use viscoplastic_model_type
  use bdf2_integrator
  use idaesol_type
  use sm_material_model_type
  use integration_geometry_type
  use truchas_timers
  use truchas_logging_services
  implicit none
  private

  type, public :: viscoplastic_solver
    private
    real(r8) :: strain_limit, atol, rtol, ntol
    logical :: use_bdf2

    type(viscoplastic_model), public, pointer :: model => null()
    type(integration_geometry), pointer :: ig => null() ! unowned reference

    type(idaesol) :: integrator
    class(idaesol_model), pointer :: integ_model => null()

    type(bdf2_control) :: bdf2control
    type(bdf2_state) :: bdf2state
  contains
    procedure :: init
    procedure :: compute_plastic_strain
    procedure :: compute_precon
    procedure, private :: integrate
    final :: viscoplastic_solver_finalize
  end type viscoplastic_solver

contains

  !! WARNING: Not pure, due to deallocation of polymorphic integ_model. Do not
  !! make local variables of type viscoplastic_solver in pure subroutines, it
  !! will leak memory. These variables must be pointers so they are valid
  !! targets in idaesol, according to NAG.
  impure elemental subroutine viscoplastic_solver_finalize(this)
    type(viscoplastic_solver), intent(inout) :: this
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%integ_model)) deallocate(this%integ_model)
  end subroutine


  subroutine init(this, params, ig, matl_model)

    use parameter_list_type
    use viscoplastic_jacob_idaesol_model_type
    use viscoplastic_jfree_idaesol_model_type

    class(viscoplastic_solver), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    type(integration_geometry), intent(in), target :: ig
    type(sm_material_model), intent(in), target :: matl_model

    integer :: stat
    character(:), allocatable :: errmsg, solver_type
    real(r8) :: atol(6)

    this%ig => ig

    allocate(this%model)
    call this%model%init(matl_model)

    if (.not.params%is_parameter("nlk-tol")) call params%set("nlk-tol", 1d-2) ! default, passed to idaesol
    call params%get('strain-limit', this%strain_limit, default=1d-10)
    call params%get('atol', this%atol, default=1d-12)
    call params%get('rtol', this%rtol, default=1d-3)
    call params%get('nlk-tol', this%ntol)

    call params%get('solver', solver_type, default="bdf2")
    select case (solver_type)
    case ("bdf2")
      atol = this%atol
      call bdf2_create_state(this%bdf2state, 6)
      call bdf2_set_param(this%bdf2control, atol=atol, rtol=this%rtol, mvec=0, ntol=this%ntol)
      this%use_bdf2 = .true.
    case ("jacobian")
      this%use_bdf2 = .false.
      call params%set('nlk-max-vec', 0)
      call params%set('nlk-max-itr', 50)
      call params%set('pc-freq', 1)
      allocate(viscoplastic_jacob_idaesol_model :: this%integ_model)
    case ("jfree")
      this%use_bdf2 = .false.
      call params%set('nlk-max-vec', 0)
      call params%set('nlk-max-itr', 10)
      call params%set('pc-freq', 1)
      allocate(viscoplastic_jfree_idaesol_model :: this%integ_model)
    case default
      call TLS_fatal("Invalid selection for viscoplastic_solver")
    end select

    if (.not.this%use_bdf2) then
      call init_integ_model(this%integ_model, this%model, this%atol, this%rtol)
      call this%integrator%init(this%integ_model, params, stat, errmsg)
      if (stat /= 0) call TLS_fatal('Failed to build viscoplastic integrator: '//errmsg)
    end if

  contains

    !! Wrapper around type-specific init methods
    subroutine init_integ_model(integ_model, model, atol, rtol)
      class(idaesol_model) :: integ_model
      type(viscoplastic_model), target :: model
      real(r8), intent(in) :: atol, rtol
      select type (integ_model)
      type is (viscoplastic_jacob_idaesol_model)
        call integ_model%init(model, atol, rtol)
      type is (viscoplastic_jfree_idaesol_model)
        call integ_model%init(model, atol, rtol)
      end select
    end subroutine init_integ_model

  end subroutine init


  subroutine compute_plastic_strain(this, dt, temperature, lame1, lame2, vof, strain_pc, &
      strain_total_old, strain_thermal_old, strain_plastic_old, dstrain_plastic_dt_old, &
      strain_total_new, strain_thermal_new, strain_plastic_new, dstrain_plastic_dt_new)

    use ieee_arithmetic, only: ieee_is_normal

    class(viscoplastic_solver), intent(inout), target :: this
    real(r8), intent(in) :: dt, temperature(:), lame1(:), lame2(:)
    real(r8), intent(in), target :: vof(:,:), strain_pc(:,:), &
        strain_thermal_old(:,:), strain_total_old(:,:), &
        strain_thermal_new(:,:), strain_total_new(:,:)
    real(r8), intent(in) :: strain_plastic_old(:,:), dstrain_plastic_dt_old(:,:)
    real(r8), intent(out) :: strain_plastic_new(:,:), dstrain_plastic_dt_new(:,:)

    integer :: p, j
    real(r8) :: strain_rate_old, strain_rate_new, max_strain!, rate_change

    call start_timer("viscoplasticity")

    do p = 1, this%ig%npt
      j = this%ig%pcell(p)
      associate(u0 => strain_plastic_old(:,p), udot0 => dstrain_plastic_dt_old(:,p), &
          &     u1 => strain_plastic_new(:,p), udot1 => dstrain_plastic_dt_new(:,p))
        call this%model%set_state(dt, temperature(j), lame1(j), lame2(j), vof(:,j), strain_pc(:,j), &
            strain_total_old(:,p), strain_thermal_old(:,j), &
            strain_total_new(:,p), strain_thermal_new(:,j))

        if (this%model%is_elastic()) then
          u1 = 0
          udot1 = 0
          cycle
        end if

        strain_rate_old = this%model%strain_rate(0.0_r8, u0)
        strain_rate_new = this%model%strain_rate(dt, u0)
        max_strain = max(strain_rate_old, strain_rate_new) * dt
        !rate_change = max(strain_rate_new / strain_rate_old, strain_rate_old / strain_rate_new)

        ! The rate_change threshold was present in the legacy code, and it seems to
        ! moderately speed up calculations, but the result is less stable. I've
        ! disabled it for now. -zjibben
        if (dt == 0 .or. max_strain < this%strain_limit) then ! .or. rate_change < 1.1_r8) then
          u1 = u0 + (0.5_r8 * dt) * udot0
          call this%model%compute_udot(0.5_r8 * dt, u1, udot1)
          u1 = u1 + (0.5_r8 * dt) * udot1
        else
          call this%integrate(dt, strain_rate_old, u0, udot0, u1)
        end if

        call this%model%compute_udot(dt, u1, udot1)
      end associate
    end do

    ASSERT(all(ieee_is_normal(strain_plastic_new)))
    ASSERT(all(ieee_is_normal(dstrain_plastic_dt_new)))

    call stop_timer("viscoplasticity")

  end subroutine compute_plastic_strain


  subroutine integrate(this, dt, strain_rate_old, strain_plastic_old, dstrain_plastic_dt_old, &
      strain_plastic_new)

    class(viscoplastic_solver), intent(inout) :: this
    real(r8), intent(in) :: dt, strain_rate_old, strain_plastic_old(:), dstrain_plastic_dt_old(:)
    real(r8), intent(out) :: strain_plastic_new(:)

    integer :: stat
    character(128) :: msg
    real(r8) :: idt

    call start_timer("integrate")

    idt = min(this%strain_limit / strain_rate_old, dt) ! / 10

    if (this%use_bdf2) then
      call bdf2_init_state(this%bdf2state, strain_plastic_old, 0.0_r8, hstart=idt)
      call bdf2_integrate(this%bdf2state, this%bdf2control, tout=dt, stat=stat, rhs=dstrain_dt)
      if (stat /= SOLN_AT_TOUT) then
        write(msg,"(a,i6)") "Viscoplastic strain integrator failed: ", stat
        call TLS_fatal(trim(msg))
      end if
      strain_plastic_new = bdf2_interpolate_solution(this%bdf2state, dt)
    else
      call this%integrator%set_initial_state(0.0_r8, strain_plastic_old, dstrain_plastic_dt_old)
      call this%integrator%integrate(idt, stat, tout=dt)
      if (stat /= 1) then
        write(msg,"(a,i6)") "Viscoplastic strain integrator failed: ", stat
        call TLS_fatal(trim(msg))
      end if

      !call this%integrator%get_last_state_copy(strain_plastic_new)
      call this%integrator%get_interpolated_state(dt, strain_plastic_new)
    end if

    call stop_timer("integrate")

  contains

    subroutine dstrain_dt(t, u, udot)
      real(r8), intent(in) :: t, u(:)
      real(r8), intent(out) :: udot(:)
      call this%model%compute_udot(t, u, udot)
    end subroutine dstrain_dt

  end subroutine integrate


  !! Compute the Jacobian contribution for the stress-strain solver. Derivative of plastic strain
  !! rate wrt plastic strain.
  subroutine compute_precon(this, temperature, lame1, lame2, vof, strain_total, strain_thermal, &
      strain_pc, strain_plastic, precon)

    class(viscoplastic_solver), intent(inout) :: this
    real(r8), intent(in) :: temperature(:), lame1(:), lame2(:)
    real(r8), intent(in), target :: vof(:,:), strain_total(:,:), strain_thermal(:,:), strain_pc(:,:)
    real(r8), intent(in) :: strain_plastic(:,:)
    real(r8), intent(out) :: precon(:,:,:)

    real(r8), parameter :: zeros(6) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
    real(r8), parameter :: dt = huge(1.0_r8)
    integer :: p, j

    call start_timer("viscoplasticity")

    do p = 1, this%ig%npt
      j = this%ig%pcell(p)
      call this%model%set_state(dt, temperature(j), lame1(j), lame2(j), vof(:,j), strain_pc(:,j), &
          strain_total(:,p), strain_thermal(:,j), zeros, zeros)

      if (this%model%is_elastic()) then
        precon(:,:,p) = 0
        cycle
      end if

      call this%model%compute_precon(0.0_r8, strain_plastic(:,p), precon(:,:,p))
    end do

    call stop_timer("viscoplasticity")

  end subroutine compute_precon

end module viscoplastic_solver_type
