!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module alloy_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use alloy_vector_type
  use alloy_model_type
  use alloy_precon_type
  use alloy_norm_type
  use unstr_mesh_type
  !use alloy_idaesol_model_type
  use new_idaesol_type
  use parameter_list_type
  implicit none
  private

  type, extends(idaesol_model), public :: alloy_solver
    type(unstr_mesh), pointer :: mesh => null()
    !type(alloy_idaesol_model) :: integ_model
    type(idaesol) :: integ
    logical :: state_is_pending = .false.
    !! Pending state
    real(r8) :: t, dt
    type(alloy_vector) :: u
    type(alloy_model),  pointer :: model  => null()
    type(alloy_precon), pointer :: precon => null()
    type(alloy_norm),   pointer :: norm   => null()
    type(parameter_list), pointer :: ic_params => null()
    real(r8), allocatable :: Clast(:,:), C(:,:), Cdot(:,:)
    integer :: num_comp
  contains
    procedure :: init
    procedure :: step
    procedure :: commit_pending_state
    procedure :: get_cell_heat_copy, get_cell_heat_view
    procedure :: get_cell_temp_copy, get_cell_temp_view
    procedure :: get_face_temp_copy, get_face_temp_view
    procedure :: get_cell_conc_view
    procedure :: avg_conc
    procedure :: get_C_liq, get_C_sol
    procedure :: get_liq_frac_view
    procedure :: get_cell_temp_grad
    procedure :: get_stepping_stats
    procedure :: last_step_size
    procedure :: last_time
    procedure :: set_initial_state
    procedure :: restart
    procedure :: set_adv_heat
    generic   :: set_cdot => set_cdot_one, set_cdot_all
    procedure, private :: set_cdot_one, set_cdot_all
    !! Deferred procedures from IDAESOL_MODEL
    procedure :: alloc_vector
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
    procedure :: normalize
    final :: alloy_solver_delete
  end type

contains

  !! Final subroutine for alloy_solver objects.
  subroutine alloy_solver_delete (this)
    type(alloy_solver), intent(inout) :: this
    if (associated(this%precon)) deallocate(this%precon)
    if (associated(this%norm)) deallocate(this%norm)
    if (associated(this%ic_params)) deallocate(this%ic_params)
  end subroutine

  subroutine init(this, model, params, stat, errmsg)

    use parallel_communication, only: is_IOP

    class(alloy_solver), intent(out), target :: this
    type(alloy_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, output_unit
    type(parameter_list), pointer :: plist
    logical :: verbose_stepping
    real(r8), allocatable :: C0(:)

    this%model => model
    this%mesh  => model%mesh

    this%num_comp = this%model%num_comp

    allocate(this%C(this%num_comp,this%mesh%ncell))
    allocate(this%Clast, this%Cdot, mold=this%C)

    !! Initial uniform solute concentration
    call params%get('concentration', C0, stat, errmsg)
    if (stat /= 0) return
    if (size(C0) /= this%num_comp) then
      stat = 1
      errmsg = 'invalid concentration size'
      return
    else if (any(C0 < 0)) then
      stat = 1
      errmsg = 'invalid concentration values'
      return
    end if

    do j = 1, this%mesh%ncell
      this%C(:,j) = C0
      this%Clast(:,j) = C0
    end do
    this%Cdot = 0

    allocate(this%norm)
    plist => params%sublist('norm')
    call this%norm%init(model, plist, stat, errmsg)
    if (stat /= 0) return

    allocate(this%precon)
    plist => params%sublist('precon')
    call this%precon%init(model, plist, stat, errmsg)
    if (stat /= 0) return

    call this%model%init_vector(this%u)

    block
      integer :: stat
      character(:), allocatable :: errmsg
      call this%integ%init(this, params, stat, errmsg)
      INSIST(stat == 0)
    end block

    call params%get('verbose-stepping', verbose_stepping)
    if (is_IOP .and. verbose_stepping) then
      call params%get('output-unit', output_unit)
      call this%integ%set_verbose_stepping(output_unit)
    end if

    !! Grab parameters for alloy_ic_solver%init
    allocate(this%ic_params)
    block
    real(r8) :: rval
    plist => params%sublist('norm')
    call plist%get('abs-t-tol', rval)
    call this%ic_params%set('atol-temp', 0.01_r8 * rval)
    call plist%get('rel-t-tol', rval)
    call this%ic_params%set('rtol-temp', 0.01_r8 * rval)
    call this%ic_params%set('max-iter', 50)
    call this%ic_params%set('method', 'SSOR')
    plist => this%ic_params%sublist('params')
    call plist%set('num-cycles', 1)
    end block

  end subroutine init

  !! Set the rate of change of solute concentration due to advection
  subroutine set_cdot_one(this, n, cdot)
    class(alloy_solver), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: cdot(:)
    ASSERT(n > 0 .and. n <= size(this%cdot,1))
    ASSERT(size(cdot) == size(this%cdot,2))
    this%cdot(n,:) = cdot
  end subroutine

  subroutine set_cdot_all(this, cdot)
    class(alloy_solver), intent(inout) :: this
    real(r8), intent(in) :: cdot(:,:)
    ASSERT(all(shape(cdot) == shape(this%cdot)))
    this%cdot(:,:) = cdot
  end subroutine

  subroutine set_adv_heat(this, q)
    class(alloy_solver), intent(inout) :: this
    real(r8), intent(in) :: q(:)
    call this%model%set_heat_source(q)
  end subroutine

  !! Get a reference to the pending/current liquid fractions
  subroutine get_liq_frac_view(this, view)
    class(alloy_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    view => this%u%lf
  end subroutine

  !! Get a reference to the pending/current cell temperatures
  subroutine get_cell_temp_view(this, view)
    class(alloy_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    view => this%u%tc
  end subroutine

  !! Get a copy of the pending/current cell temperatures
  subroutine get_cell_temp_copy(this, copy)
    class(alloy_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    n = min(size(copy), size(this%u%tc))
    copy(1:n) = this%u%tc(1:n)
  end subroutine

  !! Get a reference to the pending/current face temperatures
  subroutine get_face_temp_view(this, view)
    class(alloy_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    view => this%u%tf
  end subroutine

  !! Get a copy of the pending/current face temperatures
  subroutine get_face_temp_copy(this, copy)
    class(alloy_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    n = min(size(copy), size(this%u%tf))
    copy(1:n) = this%u%tf(1:n)
  end subroutine

  !! Get a reference to the pending/current cell enthalpies
  subroutine get_cell_heat_view(this, view)
    class(alloy_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    view => this%u%hc
  end subroutine

  !! Get a copy of the pending/current cell enthalpies
  subroutine get_cell_heat_copy(this, copy)
    class(alloy_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    n = min(size(copy), size(this%u%hc))
    copy(1:n) = this%u%hc(1:n)
  end subroutine

  !! Get a copy of the Nth solute concentration in the liquid phase
  subroutine get_C_liq(this, n, C_liq)
    class(alloy_solver), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(inout) :: C_liq(:)
    call this%model%get_liq_conc(this%C, this%u, n, C_liq)
  end subroutine

  !! Get a copy of the Nth solute concentration in the solid phase
  subroutine get_C_sol(this, n, C_sol)
    class(alloy_solver), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(inout) :: C_sol(:)
    call this%model%get_sol_conc(this%C, this%u, n, C_sol)
  end subroutine

  subroutine get_cell_conc_view(this, view)
    class(alloy_solver), intent(in), target :: this
    real(r8), pointer, intent(out) :: view(:,:)
    view => this%C
  end subroutine

  function avg_conc(this) result(C_avg)
    use parallel_communication, only: global_sum
    class(alloy_solver), intent(in) :: this
    real(r8), allocatable :: C_avg(:)
    integer :: i
    associate (n => this%mesh%ncell_onp)
      C_avg = matmul(this%C(:,:n), this%mesh%volume(:n))
      do i = 1, size(C_avg)
        C_avg(i) = global_sum(C_avg(i))
      end do
      C_avg = C_avg / global_sum(this%mesh%volume(:n))
    end associate
  end function

  !! Compute a cell-based approximation to the gradient of the pending/current
  !! cell temperatures. Note that this derived quantity is not used in the
  !! discretization of the heat equation.

  subroutine get_cell_temp_grad(this, tgrad)
    use mfd_disc_type
    class(alloy_solver), intent(inout) :: this
    real(r8), intent(out) :: tgrad(:,:)
    INSIST(size(tgrad,1) == 3)
    INSIST(size(tgrad,2) == this%model%mesh%ncell_onP)
    call this%model%mesh%face_imap%gather_offp(this%u%tf)
    call this%model%disc%compute_cell_grad(this%u%tf, tgrad)
  end subroutine

  !! Return the time of the current (committed) solution
  function last_time(this) result(t)
    class(alloy_solver), intent(in) :: this
    real(r8) :: t
    t = this%integ%last_time()
  end function

  !! Return the size of the last successful (committed) time step
  function last_step_size(this) result(h)
    class(alloy_solver), intent(in) :: this
    real(r8) :: h
    h = this%integ%last_step_size()
  end function

  subroutine get_stepping_stats(this, counters)
    class(alloy_solver), intent(in) :: this
    integer, intent(out) :: counters(:)
    ASSERT(size(counters) == 6)
    call this%integ%get_stepping_statistics(counters)
  end subroutine

  !! Advance the current solution state by a single time step from the current
  !! time to time T. STAT returns 0 if the step was successful. The advanced
  !! state is only provisional, however, until COMMIT_PENDING_STATE is called
  !! to commit it. If STAT returns a nonzero value the step has failed (STAT=1
  !! if the nonlinear iteration failed using a fresh preconditioner; STAT=2 if
  !! if the step succeeded but the predictor error exceeded its tolerance).
  !! HNEXT returns the suggested next time step size: if STAT == 0 this is the
  !! estimated step size needed to keep the estimated predictor error in its
  !! target range; otherwise this is a reduced step size to use when retrying
  !! the time step.

  subroutine step(this, t, hnext, stat)
    class(alloy_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hnext
    integer, intent(out) :: stat

    ! Update the average solute concentration (mass or volume fraction)
    this%C = this%Clast + (t - this%integ%last_time())*this%Cdot

    call this%integ%step(t, this%u, hnext, stat)
    if (stat == 0) then
      this%t = t
      this%state_is_pending = .true.
      call this%u%gather_offp !TODO: Can the be made unnecessary?
    else
      this%C = this%Clast
      call this%integ%get_last_state_copy(this%u)
      this%state_is_pending = .false.
    end if
  end subroutine

  subroutine commit_pending_state(this)
    class(alloy_solver), intent(inout) :: this
    if (this%state_is_pending) then
      this%Clast = this%C
      call this%integ%commit_state(this%t, this%u)
      this%state_is_pending = .false.
    end if
  end subroutine

  subroutine set_initial_state(this, t, temp, dt)

    use alloy_ic_solver_type

    class(alloy_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: temp(:)
    real(r8), intent(in) :: dt

    type(alloy_vector) :: udot
    type(alloy_ic_solver) :: ic

    INSIST(associated(this%model))

    this%t = t
    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute(t, this%C, this%Cdot, temp, this%u, udot)
    call this%integ%set_initial_state(t, this%u, udot)

  end subroutine set_initial_state

  subroutine restart(this, dt)

    use alloy_ic_solver_type

    class(alloy_solver), intent(inout) :: this
    real(r8), intent(in) :: dt

    type(alloy_vector) :: udot
    type(alloy_ic_solver) :: ic

    INSIST(associated(this%model))

    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute_udot(this%t, this%C, this%Cdot, this%u, udot)
    call this%integ%set_initial_state(this%t, this%u, udot)

  end subroutine

!!!! DEFERRED PROCEDURES FROM IDAESOL_MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_vector(this, vec)
    use vector_class
    class(alloy_solver), intent(in) :: this
    class(vector), allocatable, intent(out) :: vec
    type(alloy_vector), allocatable :: tmp
    allocate(tmp)
    call this%model%init_vector(tmp)
    call move_alloc(tmp, vec)
  end subroutine

  subroutine compute_f(this, t, u, udot, f)
    use vector_class
    class(alloy_solver) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u, udot ! data is intent(in)
    class(vector), intent(inout) :: f       ! data is intent(out)
    select type (u)
    class is (alloy_vector)
      select type (udot)
      class is (alloy_vector)
        select type (f)
        class is (alloy_vector)
          call this%model%compute_f(this%C, this%Cdot, t, u, udot, f)
        end select
      end select
    end select
  end subroutine

  subroutine apply_precon(this, t, u, f)
    use vector_class
    class(alloy_solver) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u  ! data is intent(in)
    class(vector), intent(inout) :: f  ! data is intent(inout)
    select type (u)
    class is (alloy_vector)
      select type (f)
      class is (alloy_vector)
        call this%precon%apply(t, u, f)
      end select
    end select
  end subroutine

  subroutine compute_precon(this, t, u, udot, dt)
    use vector_class
    class(alloy_solver) :: this
    real(r8), intent(in) :: t, dt
    class(vector), intent(inout) :: u, udot
    select type (u)
    class is (alloy_vector)
      select type (udot)
      class is (alloy_vector)
        call this%precon%compute(this%C, this%Cdot, t, u, udot, dt)
      end select
    end select
  end subroutine

  subroutine du_norm(this, t, u, du, error)
    use vector_class
    class(alloy_solver) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: u, du
    real(r8), intent(out) :: error
    select type (u)
    class is (alloy_vector)
      select type (du)
      class is (alloy_vector)
        call this%norm%compute(t, u, du, error)
      end select
    end select
  end subroutine

  subroutine normalize(this, u)
    use vector_class
    class(alloy_solver), intent(in) :: this
    class(vector), intent(inout) :: u
    integer :: j
    select type (u)
    class is (alloy_vector)
      do j = 1, this%mesh%ncell
        if (u%lf(j) > 1) u%lf(j) = 1
        if (u%lf(j) < 0) u%lf(j) = 0
      end do
    end select
  end subroutine
end module alloy_solver_type
