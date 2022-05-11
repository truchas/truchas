!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ht_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_vector_type
  use ht_model_type
  use ht_precon_type
  use ht_norm_type
  use matl_mesh_func_type
  use unstr_mesh_type
  use ht_idaesol_model_type
  use new_idaesol_type
  use parameter_list_type
  implicit none
  private

  type, public :: ht_solver
    type(matl_mesh_func), pointer :: mmf => null()
    type(unstr_mesh), pointer :: mesh => null()
    type(ht_idaesol_model) :: integ_model
    type(idaesol) :: integ
    logical :: state_is_pending = .false.
    !! Pending state
    real(r8) :: t, dt
    type(ht_vector) :: u
    type(ht_model),  pointer :: model  => null()
    type(ht_precon), pointer :: precon => null()
    type(ht_norm),   pointer :: norm   => null()
    type(parameter_list), pointer :: ic_params => null()
  contains
    procedure :: init
    procedure :: step
    procedure :: commit_pending_state
    procedure :: get_cell_heat_copy, get_cell_heat_view
    procedure :: get_cell_temp_copy, get_cell_temp_view
    procedure :: get_face_temp_copy, get_face_temp_view
    procedure :: get_cell_temp_grad
    procedure :: get_stepping_stats
    procedure :: last_step_size
    procedure :: last_time
    procedure :: set_initial_state
    procedure :: restart
    final :: ht_solver_delete
  end type

contains

  !! Final subroutine for HT_SOLVER objects.
  subroutine ht_solver_delete (this)
    type(ht_solver), intent(inout) :: this
    if (associated(this%precon)) deallocate(this%precon)
    if (associated(this%norm)) deallocate(this%norm)
    !! The solver owns the void cell and face mask components of the model.
    if (associated(this%model)) then  !TODO: SOLVER MUST BE FINALIZED BEFORE MODEL!
      if (associated(this%model%void_cell)) deallocate(this%model%void_cell)
      if (associated(this%model%void_face)) deallocate(this%model%void_face)
    end if
    if (associated(this%ic_params)) deallocate(this%ic_params)
  end subroutine

!TODO: model should hold a reference to mmf

  subroutine init(this, mmf, model, params, stat, errmsg)

    use parallel_communication, only: is_IOP

    class(ht_solver), intent(out), target :: this
    type(matl_mesh_func), intent(in), target :: mmf
    type(ht_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: output_unit
    type(parameter_list), pointer :: plist
    logical :: verbose_stepping

    this%mmf   => mmf
    this%model => model
    this%mesh  => model%mesh

    allocate(this%norm)
    plist => params%sublist('norm')
    call this%norm%init(model, plist, stat, errmsg)
    if (stat /= 0) return

    allocate(this%precon)
    plist => params%sublist('precon')
    call this%precon%init(model, plist, stat, errmsg)
    if (stat /= 0) return

    if (associated(this%model%vf_rad_prob)) then
      call this%u%init(this%mesh, this%model%vf_rad_prob%size())
    else
      call this%u%init(this%mesh)
    end if

    call this%integ_model%init(this%model, this%precon, this%norm)

    block
      integer :: stat
      character(:), allocatable :: errmsg
      call this%integ%init(this%integ_model, params, stat, errmsg)
      INSIST(stat == 0)
    end block

    call params%get('verbose-stepping', verbose_stepping)
    if (is_IOP .and. verbose_stepping) then
      call params%get('output-unit', output_unit)
      call this%integ%set_verbose_stepping(output_unit)
    end if

    !! Grab parameters for ht_ic_solver%init
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

  !! Get a reference to the pending/current cell temperatures
  subroutine get_cell_temp_view(this, view)
    class(ht_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    view => this%u%tc
  end subroutine

  !! Get a copy of the pending/current cell temperatures
  subroutine get_cell_temp_copy(this, copy)
    class(ht_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    n = min(size(copy), size(this%u%tc))
    copy(1:n) = this%u%tc(1:n)
  end subroutine

  !! Get a reference to the pending/current face temperatures
  subroutine get_face_temp_view(this, view)
    class(ht_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    view => this%u%tf
  end subroutine

  !! Get a copy of the pending/current face temperatures
  subroutine get_face_temp_copy(this, copy)
    class(ht_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    n = min(size(copy), size(this%u%tf))
    copy(1:n) = this%u%tf(1:n)
  end subroutine

  !! Get a reference to the pending/current cell enthalpies
  subroutine get_cell_heat_view(this, view)
    class(ht_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    view => this%u%hc
  end subroutine

  !! Get a copy of the pending/current cell enthalpies
  subroutine get_cell_heat_copy(this, copy)
    class(ht_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    n = min(size(copy), size(this%u%hc))
    copy(1:n) = this%u%hc(1:n)
  end subroutine

  !! Compute a cell-based approximation to the gradient of the pending/current
  !! cell temperatures. Note that this derived quantity is not used in the
  !! discretization of the heat equation.

  subroutine get_cell_temp_grad(this, tgrad)
    use mfd_disc_type
    class(ht_solver), intent(inout) :: this
    real(r8), intent(out) :: tgrad(:,:)
    INSIST(size(tgrad,1) == 3)
    INSIST(size(tgrad,2) == this%model%mesh%ncell_onP)
    call this%model%mesh%face_imap%gather_offp(this%u%tf)
    call this%model%disc%compute_cell_grad(this%u%tf, tgrad)
  end subroutine

  !! Return the time of the current (committed) solution
  function last_time(this) result(t)
    class(ht_solver), intent(in) :: this
    real(r8) :: t
    t = this%integ%last_time()
  end function

  !! Return the size of the last successful (committed) time step
  function last_step_size(this) result(h)
    class(ht_solver), intent(in) :: this
    real(r8) :: h
    h = this%integ%last_step_size()
  end function

  subroutine get_stepping_stats(this, counters)
    class(ht_solver), intent(in) :: this
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
    class(ht_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hnext
    integer, intent(out) :: stat
    call this%integ%step(t, this%u, hnext, stat)
    if (stat == 0) then
      this%t = t
      this%state_is_pending = .true.
      call this%u%gather_offp !TODO: Can the be made unnecessary?
    else
      call this%integ%get_last_state_copy(this%u)
      this%state_is_pending = .false.
    end if
  end subroutine

  subroutine commit_pending_state(this)
    class(ht_solver), intent(inout) :: this
    if (this%state_is_pending) then
      call this%integ%commit_state(this%t, this%u)
      this%state_is_pending = .false.
    end if
  end subroutine

  subroutine set_initial_state(this, t, temp, dt)

    use ht_ic_solver_type
    use parallel_communication, only: global_count

    class(ht_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: temp(:)
    real(r8), intent(in) :: dt

    integer :: j
    real(r8), allocatable :: void_vol_frac(:)
    type(ht_vector) :: udot
    type(ht_ic_solver) :: ic

    INSIST(associated(this%model))

    !! Establish the void cell and face masks in the model.
    !! N.B. We assume no one else allocates and defines these arrays.
    allocate(void_vol_frac(this%mesh%ncell), this%model%void_cell(this%mesh%ncell))
    call this%mmf%get_void_vol_frac(void_vol_frac)
    this%model%void_cell = (void_vol_frac > 1.0_r8 - 2*epsilon(1.0_r8))
    deallocate(void_vol_frac)
    if (global_count(this%model%void_cell) > 0) then
      allocate(this%model%void_face(this%mesh%nface))
      this%model%void_face = .true.
      do j = 1, this%mesh%ncell
        if (.not.this%model%void_cell(j)) then
          associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
            this%model%void_face(cface) = .false.
          end associate
        end if
      end do
      call this%mesh%face_imap%gather_offp(this%model%void_face)
    else
      deallocate(this%model%void_cell)
      this%model%void_face => null()
    end if

    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute(t, temp, this%u, udot)
    call this%integ%set_initial_state(t, this%u, udot)

  end subroutine set_initial_state

  subroutine restart(this, dt)

    use ht_ic_solver_type

    class(ht_solver), intent(inout) :: this
    real(r8), intent(in) :: dt

    type(ht_vector) :: udot
    type(ht_ic_solver) :: ic

    INSIST(associated(this%model))

    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute_udot(this%t, this%u, udot)
    call this%integ%set_initial_state(this%t, this%u, udot)

  end subroutine

end module ht_solver_type
