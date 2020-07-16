!TODO: finish documentation
!! HT_2D_SOLVER_TYPE
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! July 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HT_2d_solver_type

  use kinds
  use unstr_2d_mesh_type
  use HT_2d_model_type
  use HT_2d_precon_type
  use HT_2d_norm_type
  use HT_2d_idaesol_model_type
  use idaesol_type
  use parallel_communication
  use truchas_logging_services
  implicit none
  private

  type, public:: HT_2d_solver
    private
    type(unstr_2d_mesh), pointer :: mesh => null()  ! reference only -- do not own
    type(HT_2d_model), pointer :: model => null()   ! reference only -- do not own
    !TODO: should precon and norm be allocatable?
    type(HT_2d_precon), pointer :: precon => null()
    type(HT_2d_norm), pointer :: norm => null()
    type(HT_2d_idaesol_model) :: integ_model
    type(idaesol) :: integ
    integer :: lun = 0  ! logical unit for integrator output
    !! Pending/current state
    real(r8) :: t
    real(r8), allocatable :: u(:)
    logical :: state_is_pending = .false.
    real(r8) :: hmin
    integer :: max_step_tries
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: integrate
    procedure :: step
    procedure :: commit_pending_state
    procedure :: time
    procedure :: get_cell_heat_soln
    procedure :: get_cell_temp_soln
    procedure :: write_metrics
    procedure :: compute_initial_state
    final :: HT_2d_solver_delete
  end type HT_2d_solver

contains

  subroutine HT_2d_solver_delete(this)
    type(HT_2d_solver), intent(inout) :: this
    if (associated(this%precon)) deallocate(this%precon)
    if (associated(this%norm)) deallocate(this%norm)
  end subroutine HT_2d_solver_delete


  subroutine init(this, model, params)

    use parameter_list_type

    class(HT_2d_solver), intent(out), target :: this
    type(HT_2d_model), intent(in), target :: model
    type(parameter_list) :: params

    type(parameter_list), pointer :: plist
    character(:), allocatable :: context, errmsg
    integer :: stat
    logical :: verbose

    this%mesh => model%mesh
    this%model => model

    allocate(this%u(this%model%num_dof()))

    !! Create the preconditioner
    context = 'processing ' // params%name() // ': '
    if (params%is_sublist('preconditioner')) then
      plist => params%sublist('preconditioner')
      allocate(this%precon)
      call this%precon%init(this%model, plist)
    else
      call TLS_fatal(context//'missing "preconditioner" sublist parameter')
    end if

    !! Create the error norm
    if (params%is_sublist('error-norm')) then
      allocate(this%norm)
      plist => params%sublist('error-norm')
      call this%norm%init(this%model, plist)
    else
      call TLS_fatal(context//'missing "error-norm" sublist parameter')
    end if

    !! Create the IDAESOL model
    call this%integ_model%init(this%model, this%precon, this%norm)

    !! Create the IDAESOL integrator
    !TODO: where to write verbose output?
    if (params%is_sublist('integrator')) then
      plist => params%sublist('integrator')
      call this%integ%init(this%integ_model, plist, stat, errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('verbose-stepping', verbose, default=.false.)
      if (verbose .and. is_IOP) then
        open(newunit=this%lun,file='bdf2.log')
        call this%integ%set_verbose_stepping(this%lun)
      end if
    else
      call TLS_fatal(context//'missing "integrator" sublist parameter')
    end if

    !! BDF2 control parameters
    !TODO: default values for hmin and max_step_tries?  idaesol%bdf2_step_driver suggests
    !   max_try = 10
    !   hmin = tiny(1.0_r8)
    call params%get('hmin', this%hmin, stat=stat, errmsg=errmsg)
    if (stat /= 0) call TLS_fatal(context//errmsg)
    call params%get('max_step_tries', this%max_step_tries, stat=stat, errmsg=errmsg)
    if (stat /= 0) call TLS_fatal(context//errmsg)

  end subroutine init


  !TODO: is there a better way to set/get rel_tol, max_itr
  subroutine set_initial_state(this, t, dt, temp, rel_tol, max_itr)
    class(HT_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt, temp(:), rel_tol
    integer, intent(in) :: max_itr
    character(:), allocatable :: errmsg
    real(r8), allocatable :: udot(:)
    integer :: stat
    this%t = t
    allocate(udot(size(this%u)))
    call this%compute_initial_state(t, dt, temp, this%u, udot, rel_tol, max_itr, stat, errmsg)
    if (stat /= 0) call TLS_fatal('HT_2D_SOLVER%SET_INITIAL_STATE: '//errmsg)
    call this%integ%set_initial_state(t, this%u, udot)
  end subroutine set_initial_state

  !! Returns the current integration time.

  real(r8) function time(this)
    class(HT_2d_solver), intent(in) :: this
    time = this%integ%last_time()
  end function time

  !! Returns the current cell enthalpy solution.

  subroutine get_cell_heat_soln(this, enth)
    class(HT_2d_solver), intent(in) :: this
    real(r8), intent(inout) :: enth(:)
    ASSERT(size(enth) == this%mesh%ncell_onP)
    call this%model%get_cell_heat_copy(this%u, enth)
  end subroutine get_cell_heat_soln

  !! Returns the current cell temperature solution.

  subroutine get_cell_temp_soln(this, temp)
    class(HT_2d_solver), intent(in) :: this
    real(r8), intent(inout) :: temp(:)
    ASSERT(size(temp) == this%mesh%ncell_onP)
    call this%model%get_cell_temp_copy(this%u, temp)
  end subroutine get_cell_temp_soln

  subroutine write_metrics(this, string)
    class(HT_2d_solver), intent(in) :: this
    character(*), intent(out) :: string(:)
    ASSERT(size(string) == 2)
    call this%integ%write_metrics(string)
  end subroutine write_metrics

  !! This delegates to the IDAESOL integration driver.  A target time (TOUT)
  !! and/or (maximum) number of steps (NSTEP) is specified and the driver
  !! integrates until the target time or number of steps has been reached.
  !! The driver will adjust the time step as needed, and attempt to recover
  !! from failed steps by decreasing the time step if necessary.  The minimum
  !! and maximum step sizes (HMIN/HMAX) can be specified; if not, there is no
  !! limit.  The maximum number of attempts (MTRY) at a time step can also be
  !! specified; it defaults to a reasonable value.  The integration status is
  !! returned in STATUS; the possible values from IDAESOL are exported (see
  !! above).  The input value of HNEXT is the initial time step the driver
  !! will attempt to use.  Its return value is the time step the driver would
  !! use on the next step if it were continuing to integrate.  For the first
  !! call, HNEXT should be set to the (user-specified) initial time step, but
  !! thereafter the return value should normally be used for the next call.
  !! It permissible to change it, but there is little reason to do so in this
  !! multi-step driver scenario.

  subroutine integrate(this, hnext, status, nstep, tout, hmin, hmax, mtry)
    class(HT_2d_solver), intent(inout) :: this
    real(r8), intent(inout) :: hnext
    integer, intent(out) :: status
    integer,  intent(in), optional :: nstep, mtry
    real(r8), intent(in), optional :: tout, hmin, hmax
    call this%integ%integrate(hnext, status, nstep, tout, hmin, hmax, mtry)
    call this%integ%get_last_state_copy(this%u)
  end subroutine integrate

  !! This delegates to the IDAESOL single step subroutine.  It takes a step from
  !! the current time to time T and writes the resulting solution to THIS%U
  !! along with a suggestion for the next step size HNEXT.  STAT returns a
  !! non-zero value if the step was unsuccessful.

  subroutine step(this, t, hnext, stat)

    class(HT_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hnext
    integer, intent(out) :: stat

    real(r8) :: h

    h = t - this%time()
    call this%integ%step(h, this%hmin, this%max_step_tries, this%u, hnext, stat)
    if (stat == 0) then
      this%t = t
      this%state_is_pending = .true.
    else ! failed -- restore last good state before returning
      call this%integ%get_last_state_copy(this%u)
      this%state_is_pending = .false.
    end if

  end subroutine step

  !! This delegates to the IDAESOL subroutine. It commits the current solution
  !! the new current state for the DAE system effectively advancing by a time
  !! step.

  subroutine commit_pending_state(this)
    class(HT_2d_solver), intent(inout) :: this
    INSIST(this%state_is_pending)
    call this%integ%commit_state(this%t, this%u)
    this%state_is_pending = .false.
  end subroutine commit_pending_state

  !! This auxiliary procedure computes the consistent initial state (u, du/dt)
  !! given the initial cell temperatures.  For a typical explicit ODE system
  !! du/dt = F(t,u) this is trivial; u is given and F evaluated to get du/dt.
  !! However for our implicit index-1 DAE system F(t,u,du/dt) = 0 this is much
  !! more involved.  We are only given part of u; the remaining part must be
  !! obtained by solving the algebraic equation portion of the DAE system.
  !! Furthermore F=0 only defines du/dt for the cell enthalpies; the remaining
  !! time derivatives must be solved for (by differentiating F=0 with respect
  !! to time) or approximated (which we do here).

  !TODO: is there a better way to set/get rel_tol, max_itr
  subroutine compute_initial_state(this, t, dt, temp, u, udot, rel_tol, max_itr, stat, errmsg)

    class(HT_2d_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt, temp(:), rel_tol
    real(r8), intent(out),  target :: u(:)
    real(r8), intent(out), target :: udot(:)
    integer, intent(in) :: max_itr
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: state(:,:), Hcell(:)
    real(r8), pointer :: ucell(:), uface(:)

    ASSERT(size(temp) == this%mesh%ncell_onP)
    ASSERT(size(u) == this%model%num_dof())
    ASSERT(size(udot) == size(u))

    call this%model%get_cell_temp_view(u, ucell)  ! cell temp
    call this%model%get_face_temp_view(u, uface)  ! face temp

    ucell = temp ! set the cell temperatures from the input

    !! Compute cell enthalpy
    allocate(Hcell(this%mesh%ncell))
    call this%model%new_state_array(u, state)
    call this%model%H_of_T%compute_value(state, Hcell)
    call this%model%set_cell_heat(Hcell, u)
    deallocate(state, Hcell)

    !! Solve for the face temperatures.
    uface = 0.0_r8 ! initial guess (we could do much better)
    !TODO: better to pass state to compute_face_temp?
    call this%model%compute_face_temp(t, u, rel_tol, max_itr, stat, errmsg)
    if (stat /= 0) return

    !! Compute udot given consistent state u
    call this%model%compute_udot(t, dt, u, udot, rel_tol, max_itr, stat, errmsg)
    if (stat /= 0) return

  end subroutine compute_initial_state

end module HT_2d_solver_type
