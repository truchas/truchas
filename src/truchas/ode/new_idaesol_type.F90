!!
!! IDAESOL_TYPE
!!
!! A solver for index-1 DAE in implicit form using the BDF2 method.
!!
!! This F2008 version is adapted from much earlier F95 and F77 implementations
!! and was obtained from http://sourceforge.net/projects/calliope.mfeproject.p/
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2013  Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module new_idaesol_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use new_state_history_type
  use new_nka_type
  use vector_class
  implicit none
  private

  type, public :: idaesol
    private
    class(idaesol_model), pointer :: model => null()
    integer  :: seq = -1            ! number of steps taken
    real(r8) :: hlast               ! last step size
    real(r8) :: hpc                 ! step size built into the current preconditioner
    logical  :: usable_pc = .false. ! whether the current preconditioner is usable
    integer  :: pc_age, pc_freq
    integer  :: freeze_count = 0    ! don't increase step size for this number of steps
    integer  :: mitr = 5            ! maximum number of nonlinear iterations
    real(r8) :: ntol = 0.1_r8       ! nonlinear solver error tolerance (relative to 1)
    type(nka), allocatable :: nka   ! nonlinear Krylov accelerator
    type(state_history) :: uhist    ! solution history structure

    !! Persistent temporary workspace
    class(vector), allocatable :: u, u0, up
    class(vector), allocatable :: du, udot  ! local to bce_step

    !! Perfomance counters
    integer :: pcfun_calls = 0      ! number of calls to PCFUN
    integer :: updpc_calls = 0      ! number of calls to UPDPC
    integer :: updpc_failed = 0     ! number of UPDPC calls returning an error
    integer :: retried_bce = 0      ! number of retried BCE steps
    integer :: failed_bce = 0       ! number of completely failed BCE steps
    integer :: rejected_steps = 0   ! number of steps rejected on error tolerance
    real(r8) :: hmin = huge(1.0_r8) ! minimum step size used on a successful step
    real(r8) :: hmax = tiny(1.0_r8) ! maximum step size used on a successful step

    !! Diagnostics
    integer :: unit = 0
    logical :: verbose = .false.
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: integrate => bdf2_step_driver
    procedure :: step
    procedure :: commit_state
    procedure :: get_interpolated_state
    procedure :: get_last_state_copy
    procedure :: get_last_state_view
    procedure :: last_time
    procedure :: last_step_size
    procedure :: set_verbose_stepping
    procedure :: set_quiet_stepping
    generic :: write_metrics => write_metrics_unit, write_metrics_string
    procedure, private :: write_metrics_unit
    procedure, private :: write_metrics_string
    procedure :: get_stepping_statistics
  end type idaesol

  type, abstract, public :: idaesol_model
  contains
    procedure(alloc_vector), deferred :: alloc_vector
    procedure(compute_f), deferred :: compute_f
    procedure(apply_precon), deferred :: apply_precon
    procedure(compute_precon), deferred :: compute_precon
    procedure(du_norm), deferred :: du_norm
  end type

  abstract interface
    subroutine alloc_vector(this, vec)
      import idaesol_model, vector
      class(idaesol_model), intent(in) :: this
      class(vector), allocatable, intent(out) :: vec
    end subroutine
    subroutine compute_f(this, t, u, udot, f)
      import idaesol_model, vector, r8
      class(idaesol_model) :: this
      real(r8), intent(in) :: t
      class(vector), intent(inout) :: u, udot
      class(vector), intent(inout) :: f
    end subroutine
    subroutine apply_precon(this, t, u, f)
      import idaesol_model, vector, r8
      class(idaesol_model) :: this
      real(r8), intent(in) :: t
      class(vector), intent(inout) :: u
      class(vector), intent(inout) :: f
    end subroutine
    subroutine compute_precon(this, t, u, dt)
      import idaesol_model, vector, r8
      class(idaesol_model) :: this
      real(r8), intent(in) :: t, dt
      class(vector), intent(inout) :: u
    end subroutine
    subroutine du_norm(this, t, u, du, error)
      import :: idaesol_model, vector, r8
      class(idaesol_model) :: this
      real(r8), intent(in) :: t
      class(vector), intent(in) :: u, du
      real(r8), intent(out) :: error
    end subroutine
  end interface

  real(r8), parameter, private :: RMIN = 0.25_r8
  real(r8), parameter, private :: RMAX = 4.0_r8
  real(r8), parameter, private :: MARGIN = 3.0_r8

  !! Successful STATUS return codes:
  integer, parameter, public :: SOLVED_TO_TOUT = 1
  integer, parameter, public :: SOLVED_TO_NSTEP = 2

  !! Unsuccessful STATUS return codes:
  integer, parameter, public :: BAD_INPUT = -1
  integer, parameter, public :: STEP_FAILED = -2
  integer, parameter, public :: STEP_SIZE_TOO_SMALL = -3

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(idaesol), intent(out) :: this
    class(idaesol_model), intent(in), target :: model
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context
    integer :: maxv
    real(r8) :: vtol

    this%model => model

    context = 'processing ' // params%path() // ': '
    call params%get('nlk-max-iter', this%mitr, stat, errmsg, default=5)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%mitr < 2) then
      stat = 1
      errmsg = context//'"nlk-max-iter" must be > 1'
      return
    end if

    call params%get('nlk-tol', this%ntol, stat, errmsg, default=0.1_r8)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%ntol <= 0.0_r8 .or. this%ntol > 1.0_r8) then
      stat = 1
      errmsg = context//'"nlk-tol" must be > 0.0 and <= 1.0'
      return
    end if

    call params%get('nlk-max-vec', maxv, stat, errmsg, default=this%mitr-1)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (maxv < 0) then
      stat = 1
      errmsg = context//'"nlk-max-vec" must be >= 0'
      return
    end if
    maxv = min(maxv, this%mitr-1)

    call params%get('nlk-vec-tol', vtol, stat, errmsg, default=0.01_r8)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (vtol <= 0.0_r8) then
      stat = 1
      errmsg = context//'"nlk-vec-tol" must be > 0.0'
      return
    end if

    call params%get('pc-freq', this%pc_freq, stat, errmsg, default=huge(this%pc_freq))
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if
    if (this%pc_freq < 1) then
      stat = 1
      errmsg = context//'"pc-freq" must be > 0'
      return
    end if

    call model%alloc_vector(this%u)
    call model%alloc_vector(this%u0)
    call model%alloc_vector(this%up)
    call model%alloc_vector(this%du)
    call model%alloc_vector(this%udot)

    !! Initialize the NKA structure.
    if (maxv > 0) then
      allocate(this%nka)
      call this%nka%init(this%u, maxv)
      call this%nka%set_vec_tol(vtol)
    end if

    !! We need to maintain 3 solution vectors for quadratic extrapolation.
    call this%uhist%init(3, this%u)

  end subroutine init


  subroutine set_initial_state(this, t, u, udot)
    class(idaesol), intent(inout) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: u, udot
    call this%uhist%flush(t, u, udot)
    this%seq = 0
    this%usable_pc = .false.
  end subroutine set_initial_state


  subroutine bdf2_step_driver(this, hnext, status, nstep, tout, hmin, hmax, mtry)

    class(idaesol), intent(inout) :: this
    real(r8), intent(inout) :: hnext
    integer,  intent(out)   :: status
    integer,  intent(in) :: nstep, mtry
    real(r8), intent(in) :: tout, hmin, hmax
    optional :: nstep, tout, hmin, hmax, mtry

    integer  :: max_step, max_try, step, errc
    real(r8) :: tlast, h, t_out, h_min, h_max

    ASSERT(this%seq >= 0)

   !!!
   !!! PROCESS THE INPUT AND CHECK IT FOR CORRECTNESS

    status = 0

    !! Set the maximum number of time steps; default is unlimited.
    max_step = huge(1)
    if (present(nstep)) max_step = nstep
    if (max_step < 1) status = BAD_INPUT

    !! Set the target integration time; default is +infinity.
    t_out = huge(1.0_r8)
    if (present(tout)) t_out = tout
    if (t_out <= this%uhist%last_time()) status = BAD_INPUT

    !! Verify that at least one of NSTEP and TOUT were specified.
    if (.not.present(nstep) .and. .not.present(tout)) status = BAD_INPUT

    !! Set the bounds on the step size; default is none.
    h_min = tiny(1.0_r8)
    h_max = huge(1.0_r8)
    if (present(hmin)) h_min = hmin
    if (present(hmax)) h_max = hmax
    if (h_min < 0.0_r8 .or. hnext < h_min .or. hnext > h_max) status = BAD_INPUT

    !! Set the maximum number of attempts allowed for a step.
    max_try = 10
    if (present(mtry)) max_try = mtry
    if (max_try < 1) status = BAD_INPUT

    if (status == BAD_INPUT) return

   !!!
   !!! BEGIN TIME STEPPING

    step = 0
    STEP_LOOP: do

      h = hnext
      tlast = this%uhist%last_time()

      !! Check for a normal return before proceeding.
      if (t_out <= tlast) then
        status = SOLVED_TO_TOUT
        exit STEP_LOOP
      end if
      if (step >= max_step) then
        status = SOLVED_TO_NSTEP
        exit STEP_LOOP
      end if

      step = step + 1

      call bdf2_step(this, h, h_min, max_try, this%u, hnext, errc)
      if (errc /= 0) then
        status = STEP_FAILED
        exit STEP_LOOP   ! step failed
      end if

      !! BDF2 step was successful; commit the solution.
      call commit_state(this, tlast+h, this%u)

      !! Set the next step size.
      hnext = min(h_max, hnext)
      if (this%verbose) write(this%unit,fmt=1) hnext/h

    end do STEP_LOOP

    1 format(2x,'Changing H by a factor of ',f6.3)

  end subroutine bdf2_step_driver


  subroutine commit_state(this, t, u)

    class(idaesol), intent(inout) :: this
    real(r8), intent(in) :: t
    class(vector), intent(in) :: u

    real(r8) :: h

    h = t - this%uhist%last_time()
    ASSERT(h > 0.0_r8)
    call this%uhist%record_state(t, u)

    this%hlast = h
    this%seq = this%seq + 1
    this%freeze_count = max(0, this%freeze_count - 1)

    this%hmin = min(h, this%hmin)
    this%hmax = max(h, this%hmax)

  end subroutine commit_state


  subroutine bdf2_step(this, h, hmin, mtry, u, hnext, errc)

    type(idaesol), intent(inout) :: this
    real(r8), intent(inout) :: h
    real(r8), intent(in)    :: hmin
    integer,  intent(in)    :: mtry
    class(vector), intent(inout) :: u
    real(r8), intent(out)   :: hnext
    integer,  intent(out)   :: errc

    integer :: try

    errc = 0
    if (hmin < 0.0_r8 .or. h < hmin) errc = BAD_INPUT
    if (mtry < 1) errc = BAD_INPUT
    if (errc == BAD_INPUT) return

    try = 0
    do
      try = try + 1

      !! Check for too many attempts at a single step.
      if (try > mtry) then
        errc  = STEP_FAILED
        exit
      end if

      !! Check for a too-small step size.
      if (h < hmin) then
        errc  = STEP_SIZE_TOO_SMALL
        exit
      end if

      !! Attempt a BDF2 step.
      call step(this, this%uhist%last_time()+h, u, hnext, errc)
      if (errc == 0) exit

      !! Step failed; try again with the suggested step size.
      if (this%verbose) write(this%unit,fmt=1) hnext/h
      h = hnext

    end do

    1 format(2x,'Changing H by a factor of ',f6.3)

  end subroutine bdf2_step

  subroutine step(this, t, u, hnext, errc)

    class(idaesol), intent(inout) :: this
    real(r8), intent(in)  :: t
    class(vector), intent(inout) :: u ! data is intent(out)
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: errc

    real(r8) :: eta, etah, h, t0, tlast, perr, dt(3)
    logical  :: fresh_pc, predictor_error

    ASSERT(this%seq >= 0)

    tlast = this%uhist%last_time()
    h = t - tlast
    INSIST(h > 0)

    !! Predicted solution and base point for BCE step.
    if (this%uhist%depth() == 2) then ! trapezoid method
      if (this%verbose) write(this%unit,fmt=1) this%seq+1, tlast, h
      etah = 0.5_r8 * h
      t0 = tlast + etah
      call this%uhist%interp_state(t,  this%up, order=1)
      call this%uhist%interp_state(t0, this%u0, order=1)
    else  ! BDF2
      if (this%verbose) write(this%unit,fmt=2) this%seq+1, tlast, h
      eta = (this%hlast + h) / (this%hlast + 2.0_r8 * h)
      etah = eta * h
      t0 = tlast + (1.0_r8 - eta)*h
      call this%uhist%interp_state(t,  this%up, order=2)
      call this%uhist%interp_state(t0, this%u0, order=1)
    end if

    fresh_pc = .false.

    if (this%usable_pc) then
      this%pc_age = this%pc_age + 1
      if (this%pc_age >= this%pc_freq) this%usable_pc = .false.
    end if

    !! If the PC step size is too different than the current step size we tag
    !! it as unusable in order to preempt a possible nonlinear solve failure.
    if (this%usable_pc) then
      if (this%hpc/etah > 1.0_r8 + MARGIN) this%usable_pc = .false.
      if (etah/this%hpc > 1.0_r8 + MARGIN) this%usable_pc = .false.
    end if

    BCE: do

      !! Update the preconditioner if necessary.
      if (.not.this%usable_pc) then
        this%updpc_calls = this%updpc_calls + 1
        call this%model%compute_precon(t, this%up, etah)
        if (this%verbose) write(this%unit,fmt=3) t
        this%hpc = etah
        this%usable_pc = .true.
        this%pc_age = 0
        fresh_pc = .true.
      end if

      !! Solve the nonlinear BCE system.
      call u%copy(this%up) ! Initial solution guess is the predictor.
      call bce_step(this, t, etah, this%u0, u, errc)
      if (errc == 0) exit BCE ! the BCE step was successful.

      if (fresh_pc) then ! preconditioner was fresh; cut h and return error condition.
        this%failed_bce = this%failed_bce + 1
        hnext = 0.5_r8 * h
        this%freeze_count = 1
        errc = 1
        return
      else ! update the preconditioner and retry the nonlinear solve.
        this%retried_bce = this%retried_bce + 1
        this%usable_pc = .false.
        cycle BCE
      end if

    end do BCE

    predictor_error = (this%seq >= 3)

    if (predictor_error) then

      !! Predictor error control.
      !du = u - up
      call this%du%copy(u)
      call this%du%update(-1.0_r8, this%up)
      call this%model%du_norm(t, u, this%du, perr)
      if (perr < 4.0_r8) then ! accept the step.
        if (this%verbose) write(this%unit,fmt=4) perr
        errc = 0
      else ! reject the step; cut h and return error condition.
        this%rejected_steps = this%rejected_steps + 1
        if (this%verbose) write(this%unit,fmt=5) perr
        hnext = 0.5_r8 * h
        this%freeze_count = 1
        errc = 2
        return
      end if

      !! Select the next step size based on the predictor error and past step
      !! size history, but don't change the step size by too great a factor.
      dt(1) = h
      dt(2:) = h + this%uhist%time_deltas()
      call select_step_size(dt, perr, hnext)
      hnext = max(RMIN*h, min(RMAX*h, hnext))
      if (this%freeze_count /= 0) hnext = min(h, hnext)

    else

      if (this%verbose) write(this%unit,fmt=6)
      hnext = h
      errc = 0

    end if

    1 format(/,'TRAP step ',i6,': T=',es22.15,', H=',es22.15)
    2 format(/,'BDF2 step ',i6,': T=',es22.15,', H=',es22.15)
    3 format(2x,'Preconditioner updated at T=',es22.15)
    4 format(2x,'Step accepted: perr=',es22.15)
    5 format(2x,'Step REJECTED: perr=',es22.15)
    6 format(2x,'Step accepted: no local error control')

  end subroutine step

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  SELECT_STEP_SIZE -- Choose a new time step.
 !!

  subroutine select_step_size(dt, perr, h)

    real(r8), intent(in)  :: dt(:), perr
    real(r8), intent(out) :: h

    real(r8), parameter :: tol = 0.001_r8
    real(r8) :: a, dh, phi, dphi

    ASSERT( size(dt) == 3 )

    a = 0.5_r8*dt(1)*dt(2)*dt(3)/max(perr,0.001_r8)
    h = dt(1)

    do ! until converged -- DANGEROUS!

      phi  = h*(h + dt(1))*(h + dt(2)) - a
      dphi = (2.0_r8*h + dt(1))*(h + dt(2)) + h*(h + dt(1))

      dh = phi / dphi
      h = h - dh
      if (abs(dh) / h < tol) exit

    end do

  end subroutine select_step_size

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SOLVE_BCE_AIN -- Solve the Backward Cauchy-Euler system using AIN.
 !!
 !! The backward Cauchy-Euler (BCE) method applied to the implicit DAE
 !!
 !!     f(t,u,u') = 0
 !!
 !! yields a nonlinear system of equations for advancing the solution from a
 !! given state u0 at time t - h to the unknown solution u at time t,
 !!
 !!     f(t,u,(u-u0)/h) = 0.
 !!
 !! This subroutine solves this nonlinear system using an accelerated fixed
 !! point iteration [1] for the preconditioned system
 !! g(u) = pc(f(t,u,(u-u0)/h)) = 0:
 !!
 !!    u given
 !!    Do until converged:
 !!      du <-- g(u)
 !!      du <-- NKA(du)
 !!      u  <-- u - du
 !!    End do
 !!
 !! The procedure NKA uses information about g' gleaned from the unaccelerated
 !! correction du=g(u) and previous g values to compute an improved correction.
 !! The preconditioning function pc() is typically an approximate solution
 !! of the Newton correction equation  J*du = f(t,u,(u-u0)/h) where J is an
 !! approximation to the Jacobian of f(t,u,(u-u0)/h) as a function of u.  Thus
 !! this method can be regarded as an accelerated inexact Newton (AIN) method.
 !!
 !! The dummy procedure PCFUN evaluates the preconditioned function g.
 !!
 !! [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
 !!     weighted moving finite element code I: in one dimension", SIAM J.
 !!     Sci. Comput;, 19 (1998), pp. 728-765..
 !!

  subroutine bce_step(this, t, h, u0, u, errc)

    type(idaesol), intent(inout) :: this
    real(r8), intent(in)    :: t, h
    class(vector), intent(in)    :: u0
    class(vector), intent(inout) :: u
    integer,  intent(out)   :: errc

    integer  :: itr
    real(r8) :: error

    if (allocated(this%nka)) call this%nka%restart

    itr = 0
    do

      if (itr >= this%mitr) then  ! too many nonlinear iterations
        if (this%verbose) write(this%unit,fmt=1) itr, error
        errc = 1
        exit
      end if

      itr = itr + 1

      !! Evaluate the nonlinear function and precondition it.
      this%pcfun_calls = this%pcfun_calls + 1

      !call this%model%compute_f(t, u, (u-u0)/h, du)
      call this%udot%copy(u)
      call this%udot%update(-1.0_r8, u0)
      call this%udot%scale(1.0_r8/h)
      call this%model%compute_f(t, u, this%udot, this%du)
      call this%model%apply_precon(t, u, this%du)

      !! NKA accelerated correction.
      if (allocated(this%nka)) call this%nka%accel_update(this%du)

      !! Next solution iterate.
      !u  = u - du
      call u%update(-1.0_r8, this%du)

      !! Error estimate.
      call this%model%du_norm(t, u, this%du, error)
      if (this%verbose) write(this%unit,fmt=3) itr, error

      !! Check for convergence.
      if (((error < this%ntol) .and. (itr > 1)) .or. (error < 0.01_r8 * this%ntol)) then
        if (this%verbose) write(this%unit,fmt=2) itr, error
        errc = 0
        exit
      end if

    end do

    1 format(2x,'NLK BCE solve FAILED: ',i3,' iterations (max), error=',es22.15)
    2 format(2x,'NLK BCE solve succeeded: ',i3,' iterations, error=',es22.15)
    3 format(2x,i3,': error=',es22.15)

  end subroutine bce_step

  subroutine set_verbose_stepping(this, unit)
    class(idaesol), intent(inout) :: this
    integer, intent(in) :: unit
    this%unit = unit
    this%verbose = .true.
  end subroutine

  subroutine set_quiet_stepping(this)
    class(idaesol), intent(inout) :: this
    this%verbose = .false.
  end subroutine

  subroutine get_last_state_view(this, view)
    class(idaesol), intent(in) :: this
    class(vector), pointer :: view
    call this%uhist%get_last_state_view(view)
  end subroutine

  subroutine get_last_state_copy(this, copy)
    class(idaesol), intent(in) :: this
    class(vector), intent(inout) :: copy  ! data intent(out)
    call this%uhist%get_last_state_copy(copy)
  end subroutine

  function last_time(this) result(t)
    class(idaesol), intent(in) :: this
    real(r8) :: t
    t = this%uhist%last_time()
  end function

  function last_step_size(this) result(h)
    class(idaesol), intent(in) :: this
    real(r8) :: h
    h = this%hlast
  end function

  subroutine get_interpolated_state(this, t, u)
    class(idaesol), intent(in) :: this
    real(r8), intent(in) :: t
    class(vector), intent(inout) :: u ! data intent(out)
    call this%uhist%interp_state(t, u)
  end subroutine

  subroutine write_metrics_unit(this, unit)
    class(idaesol), intent(in) :: this
    integer, intent(in) :: unit
    character(80) :: string(2)
    call write_metrics_string(this, string)
    write(unit,'(/,(a))') string
  end subroutine

  subroutine write_metrics_string(this, string)
    class(idaesol), intent(in) :: this
    character(*), intent(out) :: string(:)
    ASSERT(size(string) == 2)
    write(string(1),fmt='(a,i6,a,es11.5,a,es9.3)') &
      'STEP=', this%seq, ', T=', this%uhist%last_time(), ', H=', this%hlast
    write(string(2),fmt='(a,i7.7,":",i5.5,a,5(i4.4,:,":"))') &
      'NFUN:NPC=', this%pcfun_calls, this%updpc_calls, &
      ', NPCF:NNR:NNF:NSR=', this%updpc_failed, &
      this%retried_bce, this%failed_bce, this%rejected_steps
  end subroutine

  subroutine get_stepping_statistics(this, counters)
    class(idaesol), intent(in) :: this
    integer, intent(out) :: counters(:)
    ASSERT(size(counters) == 6)
    counters(1) = this%pcfun_calls
    counters(2) = this%updpc_calls
    counters(3) = this%updpc_failed
    counters(4) = this%retried_bce
    counters(5) = this%failed_bce
    counters(6) = this%rejected_steps
  end subroutine

end module new_idaesol_type
