!!
!! FHT_SOLVER_TYPE
!!
!! This module defines the concrete solver for thermal transport coupled to
!! free-surface flow with dynamic void. It owns the FHT model, preconditioner,
!! residual norm, state history, nonlinear acceleration state, dynamic void
!! masks, and current/pending thermal state used by the common
!! thermal_species_solver interface.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!
!! Notes
!!
!! External state getters currently expose the pending state when one exists.
!! That matches normal use after initialization and accepted steps, but it is
!! not suitable for a future rejected-step/redo path, which would need access
!! to the last committed state.
!!

#include "f90_assert.fpp"

module fht_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use thermal_species_solver_class
  use fht_model_type
  use fht_vector_type
  use fht_norm_type
  use fht_precon_type
  use matl_mesh_func_type
  use new_state_history_type, only: state_history
  use new_nka_type, only: nka
  use vector_class
  use unstr_mesh_type
  use parallel_communication
  use enthalpy_advector_class
  use parameter_list_type
  implicit none
  private

  type, public, extends(thermal_species_solver) :: fht_solver
    type(matl_mesh_func), pointer :: mmf => null() ! unowned reference
    type(fht_model) :: model
    logical :: model_is_initialized = .false.
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(fht_precon) :: precon
    type(fht_norm) :: norm
    type(parameter_list) :: ic_params
    class(enthalpy_advector), allocatable :: hadv
    
    real(r8) :: epsilon
    integer  :: seq = 0
    logical  :: verbose
    integer  :: unit
    integer  :: num_residual = 0
    integer  :: num_precon_compute = 0
    integer  :: num_precon_apply = 0
    
    !! Last state
    real(r8), allocatable :: Hlast(:)
    logical,  allocatable :: last_void_cell(:)
    logical,  allocatable :: last_void_face(:)
    logical,  allocatable :: last_tot_void_cell(:)
    
    !! Pending state
    logical  :: state_is_pending = .false.
    real(r8) :: t
    type(fht_vector) :: u
    real(r8), allocatable :: H(:)
    logical,  allocatable :: void_cell(:)
    logical,  allocatable :: void_face(:)
    logical,  allocatable :: tot_void_cell(:)
    
    !! ODE integrator data
    type(state_history) :: uhist
    type(nka) :: accel
    integer  :: max_itr
  contains
    procedure :: init
    procedure :: step
    procedure :: commit_step
    procedure :: restart
    procedure :: set_initial_state
    procedure :: last_time
    procedure :: log_step_stats
    procedure :: get_cell_temp_copy
    procedure :: get_cell_heat_copy
    procedure :: get_face_temp_copy
    procedure :: get_face_temp_view
    procedure :: get_cell_temp_grad
    procedure :: set_ext_enthalpy_rate
    procedure :: update_moving_vf
    procedure :: add_moving_vf_events
  end type fht_solver
  
contains

  subroutine init(this, mesh, mmf, params, hadv, stat, errmsg)

    use truchas_env, only: output_file_name

    !! TARGET is needed because assemble_solver initializes persistent
    !! pointer references to components of this solver.
    class(fht_solver), intent(inout), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout), target :: params
    class(enthalpy_advector), allocatable, intent(inout) :: hadv
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: solver_params, model_params, norm_params
    integer :: lun
    logical :: verbose_stepping

    solver_params => params%sublist('solver')
    model_params => params%sublist('model')

    call this%model%init(mesh, mmf, model_params, stat, errmsg)
    if (stat /= 0) return
    this%model_is_initialized = .true.

    call solver_params%get('verbose-stepping', verbose_stepping)
    if (verbose_stepping) then
      lun = -1
      if (is_IOP) open(newunit=lun,file=output_file_name('bdf2.out'),position='rewind',action='write')
      call solver_params%set('unit', lun)
      norm_params => solver_params%sublist('norm')
      call norm_params%set('unit', lun)
    end if

    call assemble_solver(this, mmf, solver_params, stat, errmsg)
    if (stat /= 0) return
    call move_alloc(hadv, this%hadv)
    stat = 0

  end subroutine init

  subroutine assemble_solver(this, mmf, params, stat, errmsg)

    !! TARGET allows the preconditioner and norm to retain references to
    !! solver components such as model.
    type(fht_solver), intent(inout), target :: this
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    
    type(parameter_list), pointer :: plist
  
    ASSERT(this%model_is_initialized)

    this%mmf   => mmf
    this%mesh  => this%model%mesh
    
    !! Storage for the various state vectors
    call this%model%init_vector(this%u)
    allocate(this%H(this%mesh%ncell_onP), this%Hlast(this%mesh%ncell_onP))
    allocate(this%void_cell(this%mesh%ncell), this%last_void_cell(this%mesh%ncell))
    allocate(this%void_face(this%mesh%nface), this%last_void_face(this%mesh%nface))
    allocate(this%tot_void_cell(this%mesh%ncell), this%last_tot_void_cell(this%mesh%ncell))
    
    plist => params%sublist('precon')
    call this%precon%init(this%model, plist, stat, errmsg)
    if (stat /= 0) return
    plist => params%sublist('norm')
    call this%norm%init(this%model, plist, stat, errmsg)
    if (stat /= 0) return
    
    !! Parameters for fht_ic_solver. These govern only the auxiliary
    !! consistent face-temperature solve used during initial-state setup.
    block
      type(parameter_list), pointer :: plist
      call this%ic_params%set('atol-temp', 1.0e-6_r8)
      call this%ic_params%set('rtol-temp', 1.0e-6_r8)
      call this%ic_params%set('max-iter', 50)
      call this%ic_params%set('method', 'SSOR')
      plist => this%ic_params%sublist('params')
      call plist%set('num-cycles', 1)
    end block
    
    !! Setup the backward-Euler integrator.
    block
      integer  :: nlk_max_vec
      real(r8) :: nlk_vec_tol
      type(fht_vector) :: mold
      call params%get('nlk-max-vec', nlk_max_vec)
      call this%model%init_vector(mold)
      call this%accel%init(mold, nlk_max_vec)
      call params%get('nlk-vec-tol', nlk_vec_tol)
      call this%accel%set_vec_tol(nlk_vec_tol)
      call this%uhist%init(2, mold)
      call params%get('nlk-max-iter', this%max_itr)

      call params%get('epsilon', this%epsilon)
      call params%get('verbose-stepping', this%verbose)
      if (this%verbose) call params%get('unit', this%unit)
      this%verbose = this%verbose .and. is_IOP
    end block

    stat = 0
    
  end subroutine assemble_solver
  
  !! Tentatively advance the thermal state. The computed state must be
  !! committed by commit_pending_state to be accepted and retained.
  !!
  !! Advancing the state is more involved here because material advection can
  !! change the void classification of cells and faces from one step to the
  !! next. The thermal solve uses dummy equations for void cells and faces,
  !! advances enthalpy density in essentially void cells from advection alone,
  !! and then imposes dummy temperatures on totally void cells and void faces.
  !!
  !! TODO: Consider moving the computation of dQ outside this procedure.
  !! Keeping it here makes this module explicitly dependent on flow code.

  subroutine advance_state(this, t, stat)
  
    type(fht_solver), intent(inout), target :: this
    real(r8), intent(in) :: t
    integer, intent(out) :: stat
  
    integer :: j, k, n
    class(vector), pointer :: ulast
    real(r8) :: dt, void_vol_frac(this%mesh%ncell)
    real(r8), allocatable :: dQ(:), Tmin(:), Tmax(:)
    real(r8), allocatable :: Tcell(:)
    integer, allocatable :: fnbr(:,:)

    !! Update the void and totally-void domain masks.
    call this%mmf%get_void_vol_frac(void_vol_frac)
    this%void_cell = (void_vol_frac > 1.0_r8 - this%epsilon)
    this%tot_void_cell = (void_vol_frac == 1.0_r8)
    
    !! Advance the enthalpy density due to advection. dQ is the advected
    !! enthalpy increment, and Tmin and Tmax give the min and max temperatures
    !! of the enthalpy parcels advected into (or remaining in) the cells. These
    !! are used later when solving for temperature given the enthalpy density.
    call this%uhist%get_last_state_view(ulast)
    select type (ulast)
    type is (fht_vector)
      allocate(dQ(this%mesh%ncell_onP), Tmin(this%mesh%ncell_onP), Tmax(this%mesh%ncell_onP))
      call this%hadv%get_advected_enthalpy(ulast%tc(:this%mesh%ncell_onP), dQ, Tmin, Tmax)
    end select
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        this%H(j) = this%Hlast(j)
      else
        this%H(j) = this%Hlast(j) + dQ(j) / this%mesh%volume(j)
      end if
    end do
    deallocate(dQ)
    
    !! Compute the predicted state for the thermal transport solve.
        
    !! Linear extrapolation from previous states.  This is the desired choice
    !! for components that have not switched type from the previous step.
    call this%uhist%interp_state(t, this%u, order=1)
    
    !! Set the predicted temperature at cells that have switched type.
    do j = 1, this%mesh%ncell_onP
      if (this%void_cell(j)) then
        this%u%tc(j) = this%model%void_temp ! dummy value
      else if (this%last_void_cell(j)) then ! newly non-void cell
        !! Use the temperature corresponding to the advected enthalpy density.
        call this%model%T_of_H%compute(j, this%H(j), Tmin(j), Tmax(j), this%u%tc(j))
      end if
    end do
    
    !! Face-to-non-void-cell data structure; correct for on-process faces only.
    allocate(fnbr(2,this%mesh%nface))
    fnbr = 0
    do j = 1, this%mesh%ncell
      if (this%void_cell(j)) cycle
      associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
        do k = 1, size(cface)
          n = 1
          if (fnbr(1,cface(k)) /= 0) n = 2
          fnbr(n,cface(k)) = j
        end do
      end associate
    end do
    
    !! Void face mask; correct for all faces.
    this%void_face = (fnbr(1,:) == 0)
    call this%mesh%face_imap%gather_offp(this%void_face)
    
    !! Set the predicted temperature at faces that have switched type.
    allocate(Tcell(this%mesh%ncell))
    Tcell = this%u%tc
    call this%mesh%cell_imap%gather_offp(Tcell)
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j)) then
        this%u%tf(j) = this%model%void_temp ! dummy value
      else if (this%last_void_face(j)) then ! newly non-void face
        !! Use the average of neighboring non-void predicted cell temperatures.
        if (fnbr(2,j) == 0) then  ! boundary of non-void region
          this%u%tf(j) = Tcell(fnbr(1,j))
        else  ! interior of non-void region
          this%u%tf(j) = 0.5_r8 * (Tcell(fnbr(1,j)) + Tcell(fnbr(2,j)))
        end if
      end if
    end do
    deallocate(Tcell, fnbr)
    
    !! Correct predicted face temperatures with Dirichlet boundary data.
    if (allocated(this%model%thermal%bc_dir)) then
      call this%model%thermal%bc_dir%compute(t)
      do j = 1, size(this%model%thermal%bc_dir%index)
        n = this%model%thermal%bc_dir%index(j)
        if (n <= this%mesh%nface_onP) this%u%tf(n) = this%model%thermal%bc_dir%value(j)
      end do
    end if
    
    !! Set the current void context for the thermal transport model.
    this%model%void_cell => this%void_cell
    this%model%void_face => this%void_face
    
    !! Solve the backward Euler time step system to advance the temperatures
    !! at non-void cells and faces.  Void cell and face values are dummies.
    dt = t - this%uhist%last_time()
    call backward_euler_solve(this, t, dt, this%H, this%u, stat)
    if (stat /= 0) return
    
    !! Advance the cell temperature for void cells.  For essentially-void cells
    !! this is the temperature corresponding to the advanced enthalpy density.
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        this%u%tc(j) = this%model%void_temp
      else if (this%void_cell(j)) then  ! essentially void cell
        call this%model%T_of_H%compute(j, this%H(j), Tmin(j), Tmax(j), this%u%tc(j))
      end if
    end do
    deallocate(Tmin, Tmax)
    
    !! Evaluate the advanced enthalpy density using the advanced cell
    !! temperatures. We do this even for essentially void cells to ensure the
    !! enthalpy densities and cell temperatures are exactly consistent.
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        this%H(j) = 0.0_r8
      else
        call this%model%thermal%H_of_T%compute_value(j, this%u%tc(j), this%H(j))
      end if
    end do
    
    !! Set the advanced temperature for void faces.
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j)) this%u%tf(j) = this%model%void_temp
    end do
    
    this%t = t
    
    this%state_is_pending = .true.
    
    !! Diagnostic output
    
    ! We've computed new cell and face temperatures, but the residual of the
    ! system isn't exactly zero.  Consider doing the following.
    ! 1. Evaluate heat fluxes using the computed temperatures, using averaged
    !    values at the faces to get consistent fluxes, because the algebraic
    !    face temperature equations were not exactly solved. (Consider actually
    !    changing the system to use averaged values for the conservation eqn;
    !    makes the system non-symmetric, but ...)
    ! 2. Advance the enthalpy densities using this flux.  This guarantees exact
    !    conservation.
    ! 3. Solve for the corresponding cell temperatures.  The resulting cell
    !    and face temperatures don't exactly solve the system still, but is
    !    the residual significantly worse than before?
    ! One issue is that the resulting cell temperatures and enthalpy densities
    ! aren't necessarily exactly consistent, unless given H we solve the H(T)
    ! relation for T to machine precision.  Probably practical to do but not
    ! sure we want to do this.
    
  end subroutine
  
  !! Accept and record the pending state computed by advance_state.
  
  subroutine commit_pending_state(this)
  
    type(fht_solver), intent(inout) :: this
    
    integer :: j
    logical :: reset_udot
    
    INSIST(this%state_is_pending)
    
    !! Commit the pending state vector to the history.
    this%seq = this%seq + 1
    call this%uhist%record_state(this%t, this%u)
    
    !! Certain state components have switched type (totally-void, essentially-
    !! void, non-void) from the last state.  For some of those components we
    !! want to reset the time derivative (first divided difference) to zero
    !! to reflect the absence of a valid history.  This is mainly for computing
    !! reasonable interpolated solutions in the last time interval.  Those
    !! components of the extrapolated predicted state that is computed on the
    !! next step are generally ignored anyway.  The exception is a non-void
    !! cell that was previously essentially void.
    
    do j = 1, this%mesh%ncell_onP
      reset_udot = (this%tot_void_cell(j) .neqv. this%last_tot_void_cell(j))
      !! We may want to also reset the time derivative of non-void cells that
      !! were previously essentially void.  This catches that case.
      !reset_udot = reset_udot .or (.not.this%void_cell(j) .and. this%last_void_cell(j))
      if (reset_udot) then
        call this%uhist%revise(revise_cell_temp, xdot=0.0_r8)
      end if
    end do
    
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j) .neqv. this%last_void_face(j)) then
        call this%uhist%revise(revise_face_temp, xdot=0.0_r8)
      end if
    end do
      
    this%Hlast = this%H
    this%last_void_cell = this%void_cell
    this%last_void_face = this%void_face
    this%last_tot_void_cell = this%tot_void_cell
    
    this%state_is_pending = .false.
    
    !! Write some diagnostic information
    call write_solver_diagnostics

  contains

    subroutine revise_cell_temp(vec, value)
      class(vector), intent(inout) :: vec
      real(r8), intent(in) :: value
      select type (vec)
      type is (fht_vector)
        call vec%set_cell_temp(j, value)
      end select
    end subroutine revise_cell_temp

    subroutine revise_face_temp(vec, value)
      class(vector), intent(inout) :: vec
      real(r8), intent(in) :: value
      select type (vec)
      type is (fht_vector)
        call vec%set_face_temp(j, value)
      end select
    end subroutine revise_face_temp
  
    subroutine write_solver_diagnostics
      use truchas_logging_services
      integer :: n1, n2, n3, max_itr, max_adj
      real avg_itr, rec_rate, avg_adj
      character(80) :: msg
      n1 = global_count(this%tot_void_cell(:this%mesh%ncell_onP))
      n2 = global_count(this%void_cell(:this%mesh%ncell_onP)) - n1
      n3 = global_sum(this%mesh%ncell_onP) - n1 - n2
      write(msg,'(a,3(i0,:,"/"))') 'DS: totally/essentially/non-void cell counts = ', n1, n2, n3
      call TLS_info(msg)
      call this%model%T_of_H%get_metrics(avg_itr, max_itr, rec_rate, avg_adj, max_adj)
      write(msg,'(a,f5.2,a,i0,a)') 'DS: T(H) iterations: ', avg_itr, '(avg), ', max_itr, '(max)'
      call TLS_info(msg)
      write(msg,'(a,f5.3,a,f5.2,a,i0,a)') 'DS: T(H) salvage rate = ', rec_rate, &
          '; interval adjustments = ', avg_adj, '(avg), ', max_adj, '(max)'
      call TLS_info(msg)
    end subroutine
    
  end subroutine commit_pending_state

  !! Construct the initial thermal state and the initial divided difference
  !! used by the BDF-style state history.
  
  subroutine set_solver_initial_state(this, t, temp, dt)

    use fht_ic_solver_type
  
    type(fht_solver), intent(inout), target :: this
    real(r8), intent(in) :: t, temp(:), dt
    
    integer :: j
    real(r8), allocatable :: void_vol_frac(:)
    type(fht_vector) :: udot
    type(fht_ic_solver) :: ic
    
    ASSERT(size(temp) == this%mesh%ncell_onP)
    call this%model%set_initial_time(t)

    allocate(void_vol_frac(this%mesh%ncell))
    call this%mmf%get_void_vol_frac(void_vol_frac)
    this%void_cell = (void_vol_frac > 1.0_r8 - this%epsilon)
    this%tot_void_cell = (void_vol_frac == 1.0_r8)
    deallocate(void_vol_frac)
    
    this%void_face = .true.
    do j = 1, this%mesh%ncell
      if (.not.this%void_cell(j)) then
        associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
          this%void_face(cface) = .false.
        end associate
      end if
    end do
    call this%mesh%face_imap%gather_offp(this%void_face)
    
    !! Set the current void context for the thermal transport model.
    this%model%void_cell => this%void_cell
    this%model%void_face => this%void_face
    
    !! Compute the initial U and UDOT on non-void cells.
    call this%model%init_vector(udot)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute(t, temp, this%u, udot)
    
    !! Set the temperature on all void cells.
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        this%u%tc(j) = this%model%void_temp
      else if (this%void_cell(j)) then ! essentially void cell
        this%u%tc(j) = temp(j)
      end if
    end do
    
    !! Zero the temperature on all void faces. 
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j)) this%u%tf(j) = this%model%void_temp
    end do
    
    !! Zero the temperature time derivative on all void cells. 
    do j = 1, this%mesh%ncell_onP
      if (this%void_cell(j)) udot%tc(j) = 0.0_r8
    end do
    
    !! Zero the temperature time derivative on all void faces. 
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j)) udot%tf(j) = 0.0_r8
    end do
    
    !! Compute the enthalpy density on all cells.
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        this%H(j) = 0.0_r8
      else
        call this%model%thermal%H_of_T%compute_value(j, this%u%tc(j), this%H(j))
      end if
    end do
    
    !! Commit the initial state.
    this%seq = 0
    this%t = t
    call this%uhist%flush(t, this%u, udot)
    this%Hlast = this%H
    this%last_void_cell = this%void_cell
    this%last_void_face = this%void_face
    this%last_tot_void_cell = this%tot_void_cell
    this%state_is_pending = .false.
    
  end subroutine set_solver_initial_state
  
  !! Solve the backward Euler nonlinear system for the thermal transport
  !! problem on the current non-void domain.

  subroutine backward_euler_solve(this, t, dt, Hlast, u, stat)
  
    type(fht_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt, Hlast(:)
    type(fht_vector), intent(inout) :: u
    integer, intent(out) :: stat
  
    integer :: j, itr
    real(r8) :: H, error
    real(r8), allocatable :: Hdot(:)
    type(fht_vector) :: fvec
  
    call this%model%init_vector(fvec)

    !! Update the preconditioner.
    this%num_precon_compute = this%num_precon_compute + 1
    call this%precon%compute(t, u, dt)
    
    allocate(Hdot(this%mesh%ncell_onP))
    
    !! Compute the initial residual and norm.
    this%num_residual = this%num_residual + 1
    do j = 1, this%mesh%ncell_onP
      call this%model%thermal%H_of_T%compute_value(j, u%tc(j), H)
      Hdot(j) = (H - Hlast(j)) / dt
    end do
    call this%model%residual(t, u, Hdot, fvec)
    call this%norm%fnorm(t, u, Hdot, fvec)
    
    itr = 0
    call this%accel%restart
    do ! until converged
    
      itr = itr + 1
      
      !! Compute the next solution iterate.
      this%num_precon_apply = this%num_precon_apply + 1
      call this%precon%apply(t, u, fvec)
      call this%accel%accel_update(fvec)
      call u%update(-1.0_r8, fvec)
      
      !! Compute the residual and norm.
      this%num_residual = this%num_residual + 1
      do j = 1, this%mesh%ncell_onP
        call this%model%thermal%H_of_T%compute_value(j, u%tc(j), H)
        Hdot(j) = (H - Hlast(j)) / dt
      end do
      call this%model%residual(t, u, Hdot, fvec)
      call this%norm%fnorm(t, u, Hdot, fvec, error)
      if (this%verbose) write(this%unit,fmt=3) itr, error
      
      !! Convergence check and iteration control.
      if (itr >= 2 .and. error <= 1.0_r8) then
        if (this%verbose) write(this%unit,fmt=2) itr, error
        stat = 0
        exit
      else if (itr >= this%max_itr) then
        if (this%verbose) write(this%unit,fmt=1) itr, error
        stat = 1
        exit
      end if
      
    end do
    
    1 format(2x,'NLK BCE solve FAILED: ',i3,' iterations (max), error=',es12.5)
    2 format(2x,'NLK BCE solve succeeded: ',i3,' iterations, error=',es12.5)
    3 format(2x,i3,': error=',es12.5)

  end subroutine backward_euler_solve
  
  subroutine get_cell_temp_view(this, view)
    type(fht_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    !! External copies and views expose only on-process entries; the
    !! off-process vector entries are solver scratch space.
    view => this%u%tc(:this%mesh%ncell_onP)
  end subroutine
  
  subroutine get_cell_temp_copy(this, copy)
    class(fht_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    !! External copies and views expose only on-process entries; the
    !! off-process vector entries are solver scratch space.
    n = min(size(copy), this%mesh%ncell_onP)
    copy(1:n) = this%u%tc(1:n)
  end subroutine
  
  subroutine get_face_temp_view(this, view)
    class(fht_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    !! External copies and views expose only on-process entries; the
    !! off-process vector entries are solver scratch space.
    view => this%u%tf(:this%mesh%nface_onP)
  end subroutine
  
  subroutine get_face_temp_copy(this, copy)
    class(fht_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    !! External copies and views expose only on-process entries; the
    !! off-process vector entries are solver scratch space.
    n = min(size(copy), this%mesh%nface_onP)
    copy(1:n) = this%u%tf(1:n)
  end subroutine
  
  subroutine get_cell_heat_view(this, view)
    type(fht_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    view => this%H
  end subroutine
  
  subroutine get_cell_heat_copy(this, copy)
    class(fht_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    ASSERT(size(copy) >= size(this%H))
    copy(:size(this%H)) = this%H
  end subroutine
  
  subroutine get_void_cell_view(this, view)
    type(fht_solver), intent(in), target :: this
    logical, pointer :: view(:)
    view => this%void_cell(:size(this%H))
  end subroutine
  
  subroutine get_void_cell_copy(this, copy)
    type(fht_solver), intent(in) :: this
    logical, intent(inout) :: copy(:)
    ASSERT(size(copy) >= size(this%H))
    copy(:size(this%H)) = this%void_cell(:size(this%H))
  end subroutine
  
  function last_time(this) result(t)
    class(fht_solver), intent(in) :: this
    real(r8) :: t
    t = this%uhist%last_time()
  end function
  
  function last_step_size(this) result(dt)
    type(fht_solver), intent(in) :: this
    real(r8) :: dt
    real(r8), allocatable :: deltas(:)
    deltas = this%uhist%time_deltas()
    dt = deltas(1)
  end function
  
  subroutine get_stepping_stats(this, counters)
    type(fht_solver), intent(in) :: this
    integer, intent(inout) :: counters(:)
    ASSERT(size(counters) >= 4)
    counters(1) = this%seq
    counters(2) = this%num_residual
    counters(3) = this%num_precon_compute
    counters(4) = this%num_precon_apply
  end subroutine
  
  subroutine get_cell_temp_grad(this, tgrad)
    class(fht_solver), intent(inout) :: this
    real(r8), intent(out) :: tgrad(:,:)
    INSIST(size(tgrad,1) == 3)
    INSIST(size(tgrad,2) == this%model%mesh%ncell_onP)
    !! This getter needs off-process face temperatures for the reconstruction.
    !! They are gathered directly into the solver's scratch entries.
    call this%model%mesh%face_imap%gather_offp(this%u%tf)
    call this%model%disc%compute_cell_grad(this%u%tf, tgrad)
  end subroutine

  subroutine step(this, t, hnext, errc)
    class(fht_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hnext
    integer, intent(out) :: errc
    real(r8) :: h
    h = t - this%last_time()
    call advance_state(this, t, errc)
    hnext = merge(h/2, huge(1.0_r8), errc /= 0)
  end subroutine

  !! Keep this type-bound override as the target-free common solver interface
  !! adapter. The FHT-specific pending-state bookkeeping lives in
  !! commit_pending_state.
  subroutine commit_step(this)
    class(fht_solver), intent(inout) :: this
    call commit_pending_state(this)
  end subroutine

  subroutine restart(this, dt)
    class(fht_solver), intent(inout) :: this
    real(r8), intent(in) :: dt
    !! TODO: Implement the FHT phase-event restart/reset path. FHT carries
    !! state history in uhist; restart should recompute a derivative for the
    !! current state using DT and flush uhist, analogous to the IDAESOL-based
    !! solvers.
  end subroutine

  subroutine set_initial_state(this, t, dt, temp, conc)
    class(fht_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional :: temp(:), conc(:,:)
    ASSERT(present(temp))
    call set_solver_initial_state(this, t, temp, dt)
  end subroutine

  subroutine log_step_stats(this)
    use truchas_logging_services, only: TLS_info
    class(fht_solver), intent(in) :: this
    real(r8) :: hlast
    integer :: counters(6)
    character(len=80) :: message(2)
    hlast = last_step_size(this)
    call get_stepping_stats(this, counters)
    write(message,1) hlast, counters(2:4)
    1 format(/,'DS: dt=',es9.3,', NRES:NPC:NPA=',3(i7.7,:,':'))
    call TLS_info(message)
  end subroutine

  subroutine set_ext_enthalpy_rate(this, enthalpy_rate)
    class(fht_solver), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)
    call this%model%set_ext_enthalpy_rate(enthalpy_rate)
  end subroutine

  subroutine update_moving_vf(this)
    class(fht_solver), intent(inout) :: this
    call this%model%update_moving_vf
  end subroutine

  subroutine add_moving_vf_events(this, eventq, rank)
    use sim_event_queue_type, only: sim_event_queue
    class(fht_solver), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    integer, intent(in), optional :: rank
    call this%model%add_moving_vf_events(eventq, rank)
  end subroutine

end module fht_solver_type
