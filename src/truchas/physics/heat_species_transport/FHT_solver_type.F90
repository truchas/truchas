!!
!! FHT_SOLVER_TYPE
!!
!! This is a special-purpose solver for heat transfer that is coupled to free
!! surface flows involving void.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! CRUCIAL BEHAVIOR THAT NEEDS TO BE FIXED: THE STATE 'GET' ROUTINES ALL
!! RETURN THE PENDING STATE.  THIS EQUALS THE LAST STATE AFTER SETTING THE
!! INITIAL STATE OR AFTER COMMITTING THE PENDING STATE, WHICH IS NORMALLY
!! FINE BEHAVIOR.  HOWEVER IF THE INTENT IS TO SCRAP THE PENDING STATE AND
!! DO THE TIME STEP OVER, THESE ROUTINES ARE RETURNING THE WRONG THING --
!! THEY OUGHT TO BE RETURNING THE LAST STATE.  THIS WILL NOT CAUSE ANY
!! PROBLEMS UNTIL WE ACTUALLY IMPLEMENT STEP REDOS IN TRUCHAS.
!!

#include "f90_assert.fpp"

!#define GMV_DIAGNOSTICS

module FHT_solver_type

  use kinds, only: r8
  use FHT_model_type
  use FHT_norm_type
  use FHT_precon_type
  use TofH_type
  use matl_mesh_func_type
  use property_mesh_function
  use solution_history
  use nka_type
  use unstr_mesh_type
  use parallel_communication
  use index_partitioning
  use enthalpy_advector_class
  implicit none
  private

  type, public :: FHT_solver
    type(matl_mesh_func), pointer :: mmf => null()
    type(FHT_model),  pointer :: model => null()
    type(unstr_mesh), pointer :: mesh => null()
    type(FHT_precon) :: precon
    type(FHT_norm) :: norm
    type(TofH) :: T_of_H
    class(enthalpy_advector), allocatable :: hadv
    
    real(r8) :: epsilon
    integer  :: seq = 0
    logical  :: verbose
    integer  :: unit
    integer  :: num_compute_f = 0
    integer  :: num_precon_compute = 0
    integer  :: num_precon_apply = 0
    
    !! Last state
    real(r8), pointer :: Hlast(:) => null()
    logical,  pointer :: last_void_cell(:) => null()
    logical,  pointer :: last_void_face(:) => null()
    logical,  pointer :: last_tot_void_cell(:) => null()
    
    !! Pending state
    logical  :: state_is_pending = .false.
    real(r8) :: t
    real(r8), pointer :: u(:) => null()
    real(r8), pointer :: H(:) => null()
    logical,  pointer :: void_cell(:) => null()
    logical,  pointer :: void_face(:) => null()
    logical,  pointer :: tot_void_cell(:) => null()
    
    !! ODE integrator data
    type(history) :: uhist
    type(nka) :: accel
    integer  :: max_itr
  end type FHT_solver
  
  type, public :: FHT_solver_params
    type(FHT_norm_params) :: norm_params
    type(FHT_precon_params) :: precon_params
    real(r8) :: epsilon  ! material vol frac conduction threshold, > 0
    integer  :: nlk_max_itr
    integer  :: nlk_max_vec
    real(r8) :: nlk_vec_tol
    logical  :: verbose
    integer  :: unit
    !! Parameters for the TofH solver
    real(r8) :: TofH_tol      ! absolute temperature convergence tolerance, >= 0
    real(r8) :: TofH_delta    ! initial endpoint shift when seeking bracket, > 0
    integer  :: TofH_max_try  ! max tries at seeking a bracketing interval, >= 0
  end type FHT_solver_params
  public :: FHT_norm_params, FHT_precon_params
  public :: diff_precon_params, ssor_precon_params, boomer_amg_precon_params
  
  public :: FHT_solver_init
  public :: FHT_solver_delete
  public :: FHT_solver_advance_state
  public :: FHT_solver_commit_pending_state
  public :: FHT_solver_set_initial_state
  public :: FHT_solver_get_cell_temp_view
  public :: FHT_solver_get_cell_temp_copy
  public :: FHT_solver_get_face_temp_view
  public :: FHT_solver_get_face_temp_copy
  public :: FHT_solver_get_cell_heat_view
  public :: FHT_solver_get_cell_heat_copy
  public :: FHT_solver_get_void_cell_view
  public :: FHT_solver_get_void_cell_copy
  public :: FHT_solver_last_time
  public :: FHT_solver_last_step_size
  public :: FHT_solver_get_stepping_stats
  public :: FHT_solver_get_cell_temp_grad

contains

  subroutine FHT_solver_init (this, mmf, model, params)
  
    type(FHT_solver), intent(out) :: this
    type(matl_mesh_func), intent(in), target :: mmf
    type(FHT_model), intent(in), target :: model
    type(FHT_solver_params), intent(inout) :: params
    
    integer :: n
    procedure(pardp), pointer :: dp
  
    this%mmf   => mmf
    this%model => model
    this%mesh  => model%mesh
    
    !! Storage for the various state vectors
    allocate(this%u(FHT_model_size(model)))
    allocate(this%H(this%mesh%ncell_onP), this%Hlast(this%mesh%ncell_onP))
    allocate(this%void_cell(this%mesh%ncell), this%last_void_cell(this%mesh%ncell))
    allocate(this%void_face(this%mesh%nface), this%last_void_face(this%mesh%nface))
    allocate(this%tot_void_cell(this%mesh%ncell), this%last_tot_void_cell(this%mesh%ncell))
    
    call FHT_precon_init (this%precon, this%model, params%precon_params)
    call FHT_norm_init (this%norm, this%model, params%norm_params)
    
    call this%T_of_H%init (model%H_of_T, eps=params%TofH_tol, &
        max_try=params%TofH_max_try, delta=params%TofH_delta)
    
    !! Setup the backward-Euler integrator.
    n = FHT_model_size(model)
    call this%accel%init (n, params%nlk_max_vec)
    call this%accel%set_vec_tol (params%nlk_vec_tol)
    dp => pardp ! NB: in F2008 can make pardp an internal sub and pass directly
    call this%accel%set_dot_prod (dp)
    call create_history (this%uhist, 2, n)
    this%max_itr = params%nlk_max_itr
    
    this%epsilon = params%epsilon
    this%verbose = params%verbose .and. is_IOP
    this%unit = params%unit
    
  end subroutine FHT_solver_init
  
  function pardp (a, b) result (dp)
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    dp = global_dot_product(a, b)
  end function pardp
  
  subroutine FHT_solver_delete (this)
    type(FHT_solver), intent(inout) :: this
    call FHT_precon_delete (this%precon)
    call destroy (this%uhist)
  end subroutine FHT_solver_delete
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! ADVANCE
 !!
 !! This procedure tentatively advances the state of the heat transfer system.
 !! The computed state must be committed by COMMIT_PENDING_STATE in order to
 !! be accepted and retained.
 !!
 !! Advancing the state is a complex process in the context of flow with void.
 !! The mesh cells are partitioned into three domains: totally-void cells are
 !! those with void volume fraction equal 1; essentially-void cells are those
 !! with a void volume fraction that is nearly 1; and non-void cells are the
 !! remaining cells.  A void cell is one that is either totally or essentially
 !! void.  In the essentially-void and non-void domains we model conservation
 !! of heat and must account for increments of heat due to material advection,
 !! but we only model heat conduction in the non-void domain.  For that latter
 !  model we also have temperature variables on mesh faces.  Non-void faces are
 !! those adjacent to a non-void cell, and hence those that enter into the heat
 !  conduction system.  The remaining faces are void faces.  The heat conduction
 !! system is solved mesh-wide using dummy equations for all void cells and
 !! faces.  Afterwards the enthalpy and temperature of essentially-void cells
 !! must be advanced due to the advected heat only.  Finally dummy values are
 !! imposed for all totally-void cells and void faces.
 !!
 !! From one call to the next, the underlying distribution of material due to
 !! advection will change, particularly void.  This means that the partitioning
 !! of cells and faces into totally-, essentially-, and non-void domains is
 !! different with each call.  This is particularly relevant to the heat
 !! conduction system which computes a predicted state based on past history
 !! and so needs to take equations that have switched type between non-void
 !! and void into special account.
 !!
 !! NOTES
 !! * Consider moving the computation of dQ outside of this procedure.
 !!   Keeping it here makes this module explicitly dependent on flow code.
 !!
 
  subroutine FHT_solver_advance_state (this, t, stat)
  
    use boundary_data
  
    type(FHT_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    integer, intent(out) :: stat
  
    integer :: j, k, n
    real(r8), pointer :: ulast(:), Tcell(:), Tface(:)
    real(r8) :: dt, void_vol_frac(this%mesh%ncell)
    real(r8), allocatable :: dQ(:), Tmin(:), Tmax(:)
    integer, allocatable :: fnbr(:,:)

    !! Update the void and totally-void domain masks.
    call this%mmf%get_void_vol_frac(void_vol_frac)
    this%void_cell = (void_vol_frac > 1.0_r8 - this%epsilon)
    this%tot_void_cell = (void_vol_frac == 1.0_r8)
    
    !! Advance the heat density due to advection.  dQ is the advected heat
    !! increment, and Tmin and Tmax are give the min and max temperatures
    !! of the heat parcels advected into (or remaining in) the cells.  These
    !! are used later when solving for temperature given the heat density.
    ulast => most_recent_solution(this%uhist)
    call FHT_model_get_cell_temp_view (this%model, ulast, Tcell)
    allocate(dQ(size(Tcell)), Tmin(size(Tcell)), Tmax(size(Tcell)))
    call this%hadv%get_advected_enthalpy(Tcell, dQ, Tmin, Tmax)
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        this%H(j) = this%Hlast(j)
      else
        this%H(j) = this%Hlast(j) + dQ(j) / this%mesh%volume(j)
      end if
    end do
    deallocate(dQ)
    
    !! COMPUTE THE PREDICTED STATE FOR THE HEAT CONDUCTION SOLVE
        
    !! Linear extrapolation from previous states.  This is the desired choice
    !! for components that have not switched type from the previous step.
    call interpolate_solution (this%uhist, t, this%u, order=1)
    !call interpolate_solution (this%uhist, t, this%u, order=0)
    
    !! Set the predicted temperature at cells that have switched type.
    call FHT_model_get_cell_temp_view (this%model, this%u, Tcell)
    do j = 1, this%mesh%ncell_onP
      if (this%void_cell(j)) then
        Tcell(j) = 0.0_r8 ! dummy value
      else if (this%last_void_cell(j)) then ! newly non-void cell
        !! Use the temperature corresponding to the advected heat density.
        call this%T_of_H%compute (j, this%H(j), Tmin(j), Tmax(j), Tcell(j))
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
    call gather_boundary (this%mesh%face_ip, this%void_face)
    
    !! Set the predicted temperature at faces that have switched type.
    call FHT_model_get_face_temp_view (this%model, this%u, Tface)
    allocate(Tcell(this%mesh%ncell))
    call FHT_model_get_cell_temp_copy (this%model, this%u, Tcell)
    call gather_boundary (this%mesh%cell_ip, Tcell)
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j)) then
        Tface(j) = 0.0_r8 ! dummy value
      else if (this%last_void_face(j)) then ! newly non-void face
        !! Use the average of neighboring non-void predicted cell temperatures.
        if (fnbr(2,j) == 0) then  ! boundary of non-void region
          Tface(j) = Tcell(fnbr(1,j))
        else  ! interior of non-void region
          Tface(j) = 0.5_r8 * (Tcell(fnbr(1,j)) + Tcell(fnbr(2,j)))
        end if
      end if
    end do
    deallocate(Tcell, fnbr)
    
    !! Correct predicted face temperatures with Dirichlet boundary data.
    call bd_data_eval (this%model%bc_dir, t)
    do j = 1, size(this%model%bc_dir%faces)
      n = this%model%bc_dir%faces(j)
      if (n <= this%mesh%nface_onP) Tface(n) = this%model%bc_dir%values(1,j)
    end do
    
    !! Set the current void context for the heat transfer model.
    this%model%void_cell => this%void_cell
    this%model%void_face => this%void_face
    
    !! Solve the backward Euler time step system to advance the temperatures
    !! at non-void cells and faces.  Void cell and face values are dummies.
    dt = t - most_recent_time(this%uhist)
    call backward_euler_solve (this, t, dt, this%H, this%u, stat)
    if (stat /= 0) return
    
    !! Advance the cell temperature for void cells.  For essentially-void cells
    !! this is the temperature corresponding to the advanced heat density.
    call FHT_model_get_cell_temp_view (this%model, this%u, Tcell)
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        Tcell(j) = 0.0_r8
      else if (this%void_cell(j)) then  ! essentially void cell
        call this%T_of_H%compute (j, this%H(j), Tmin(j), Tmax(j), Tcell(j))
      end if
    end do
    deallocate(Tmin, Tmax)
    
    !! Evaluate the advanced heat density using the advanced cell temperatures.
    !! We do this even for essentially-void cells to ensure the heat densities
    !! and cell temperatures are exactly consistent.
    call FHT_model_get_cell_temp_view (this%model, this%u, Tcell)
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        this%H(j) = 0.0_r8
      else
        call pmf_eval (this%model%H_of_T, j, Tcell(j:j), this%H(j))
      end if
    end do
    
    !! Set the advanced temperature for void faces.
    call FHT_model_get_face_temp_view (this%model, this%u, Tface)
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j)) Tface(j) = 0.0_r8
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
    ! 2. Advance the heat densities using this flux.  This guarantees exact
    !    conservation.
    ! 3. Solve for the corresponding cell temperatures.  The resulting cell
    !    and face temperatures don't exactly solve the system still, but is
    !    the residual significantly worse than before?
    ! One issue is that the resulting cell temperatures and heat densities
    ! aren't necessarily exactly consistent, unless given H we solve the H(T)
    ! relation for T to machine precision.  Probably practical to do but not
    ! sure we want to do this.
    
  end subroutine FHT_solver_advance_state
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! COMMIT_PENDING_STATE
 !!
  
  subroutine FHT_solver_commit_pending_state (this)
  
    type(FHT_solver), intent(inout) :: this
    
    integer :: j, n
    logical :: reset_udot
    
    INSIST(this%state_is_pending)
    
    !! Commit the pending state vector to the history.
    this%seq = this%seq + 1
    call record_solution (this%uhist, this%t, this%u)
    
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
        n = FHT_model_cell_temp_index(this%model, j)
        call revise_history (this%uhist, n, xdot=0.0_r8)
      end if
    end do
    
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j) .neqv. this%last_void_face(j)) then
        n = FHT_model_face_temp_index(this%model, j)
        call revise_history (this%uhist, n, xdot=0.0_r8)
      end if
    end do
      
    this%Hlast = this%H
    this%last_void_cell = this%void_cell
    this%last_void_face = this%void_face
    this%last_tot_void_cell = this%tot_void_cell
    
    this%state_is_pending = .false.
    
    !! Write some diagnostic information
    call write_solver_diagnostics

#ifdef GMV_DIAGNOSTICS
    call gmv_write_state (this)
#endif
  contains
  
    subroutine write_solver_diagnostics ()
      use truchas_logging_services
      integer :: n1, n2, n3, max_itr, max_adj
      real avg_itr, rec_rate, avg_adj
      character(80) :: msg
      n1 = global_count(this%tot_void_cell(:this%mesh%ncell_onP))
      n2 = global_count(this%void_cell(:this%mesh%ncell_onP)) - n1
      n3 = global_sum(this%mesh%ncell_onP) - n1 - n2
      write(msg,'(a,3(i0,:,"/"))') 'DS: totally/essentially/non-void cell counts = ', n1, n2, n3
      call TLS_info (msg)
      call this%T_of_H%get_metrics (avg_itr, max_itr, rec_rate, avg_adj, max_adj)
      write(msg,'(a,f5.2,a,i0,a)') 'DS: T(H) iterations: ', avg_itr, '(avg), ', max_itr, '(max)'
      call TLS_info (msg)
      write(msg,'(a,f5.3,a,f5.2,a,i0,a)') 'DS: T(H) salvage rate = ', rec_rate, &
          '; interval adjustments = ', avg_adj, '(avg), ', max_adj, '(max)'
      call TLS_info (msg)
    end subroutine
    
  end subroutine FHT_solver_commit_pending_state
  
#ifdef GMV_DIAGNOSTICS
  subroutine gmv_write_state (this)
    use unstr_mesh_gmv
    type(FHT_solver), intent(in) :: this
    character(63) :: filename
    real(r8), pointer :: Tcell(:)
    real(r8), allocatable :: array(:)
    integer :: j
    write(filename,'(a,i4.4)') 'FHT_solver-gmv-', this%seq
    call gmv_open (trim(filename))
    call gmv_write_unstr_mesh (this%mesh)
    call gmv_begin_variables (this%t, this%seq)
    call FHT_model_get_cell_temp_view (this%model, this%u, Tcell)
    call gmv_write_dist_cell_var (this%mesh, Tcell, "T")
    call gmv_write_dist_cell_var (this%mesh, this%H, "H")
    allocate(array(size(Tcell)))
    do j = 1, size(Tcell)
      if (this%void_cell(j)) then
        if (this%tot_void_cell(j)) then
          array(j) = 0
        else
          array(j) = 1
        end if
      else
        array(j) = 2
      end if
    end do
    call gmv_write_dist_cell_var (this%mesh, array, "void")
    call gmv_end_variables ()
    call gmv_close ()
  end subroutine gmv_write_state
  subroutine gmv_write_field (this, array, name)
    use unstr_mesh_gmv
    type(FHT_solver), intent(in) :: this
    real(r8), intent(in) :: array(:)
    character(*), intent(in) :: name
    character(63) :: filename
    write(filename,'(a,i4.4)') 'FHT_solver-' // trim(name) // '-gmv-', this%seq
    call gmv_open (trim(filename))
    call gmv_write_unstr_mesh (this%mesh)
    call gmv_begin_variables (this%t, this%seq)
    call gmv_write_dist_cell_var (this%mesh, array, trim(name))
    call gmv_end_variables ()
    call gmv_close ()
  end subroutine gmv_write_field
  subroutine gmv_write_face_field (this, array, name)
    use unstr_mesh_gmv
    use gmvwrite_c_binding
    type(FHT_solver), intent(in) :: this
    real(r8), intent(in) :: array(:)
    character(*), intent(in) :: name
    integer, allocatable :: matnum(:)
    integer :: j
    character(63) :: filename
    write(filename,'(a,i4.4)') 'FHT_solver-face-' // trim(name) // '-gmv-', this%seq
    call gmv_open (trim(filename))

    call gmvwrite_node_data (size(this%mesh%x,dim=2), this%mesh%x(1,:), this%mesh%x(2,:), this%mesh%x(3,:))
    call gmvwrite_cell_header (size(this%mesh%fnode,dim=2))
    do j = 1, size(this%mesh%fnode,dim=2)
      call gmvwrite_cell_type ('quad', 4, this%mesh%fnode(:,j))
    end do
    call gmvwrite_material_header (3, CELLDATA)
    call gmvwrite_material_name ('x-faces')
    call gmvwrite_material_name ('y-faces')
    call gmvwrite_material_name ('z-faces')
    
    allocate(matnum(this%mesh%nface))
    do j = 1, this%mesh%nface
      if (abs(this%mesh%normal(1,j)) > 0.8*this%mesh%area(j)) then
        matnum(j) = 1
      else if (abs(this%mesh%normal(2,j)) > 0.8*this%mesh%area(j)) then
        matnum(j) = 2
      else if (abs(this%mesh%normal(3,j)) > 0.8*this%mesh%area(j)) then
        matnum(j) = 3
      else
        matnum(j) = 0
      end if
    end do
    call gmvwrite_material_ids (matnum, CELLDATA)
    
    call gmv_begin_variables (this%t, this%seq)
    call gmvwrite_variable_name_data (CELLDATA, name, array)
    call gmv_end_variables ()
    call gmv_close ()
  end subroutine gmv_write_face_field
#endif
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SET_INITIAL_STATE
 !!
  
  subroutine FHT_solver_set_initial_state (this, t, temp)
  
    type(FHT_solver), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:)
    
    integer :: j
    real(r8), pointer :: Tcell(:), Tface(:)
    real(r8), allocatable :: void_vol_frac(:), udot(:)
    target :: udot
    
    ASSERT(size(temp) == this%mesh%ncell_onP)
    
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
    call gather_boundary (this%mesh%face_ip, this%void_face)
    
    !! Set the current void context for the heat transfer model.
    this%model%void_cell => this%void_cell
    this%model%void_face => this%void_face
    
    !! Compute the initial U and UDOT on non-void cells.
    allocate(udot(size(this%u)))
    call FHT_model_compute_udot (this%model, t, temp, this%u, udot)
    
    !! Set the temperature on all void cells.
    call FHT_model_get_cell_temp_view (this%model, this%u, Tcell)
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        Tcell(j) = 0.0_r8
      else if (this%void_cell(j)) then ! essentially void cell
        Tcell(j) = temp(j)
      end if
    end do
    
    !! Zero the temperature on all void faces. 
    call FHT_model_get_face_temp_view (this%model, this%u, Tface)
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j)) Tface(j) = 0.0_r8
    end do
    
    !! Zero the temperature time derivative on all void cells. 
    call FHT_model_get_cell_temp_view (this%model, udot, Tcell)
    do j = 1, this%mesh%ncell_onP
      if (this%void_cell(j)) Tcell(j) = 0.0_r8
    end do
    
    !! Zero the temperature time derivative on all void faces. 
    call FHT_model_get_face_temp_view (this%model, udot, Tface)
    do j = 1, this%mesh%nface_onP
      if (this%void_face(j)) Tface(j) = 0.0_r8
    end do
    
    !! Compute the enthalpy density on all cells.
    call FHT_model_get_cell_temp_view (this%model, this%u, Tcell)
    do j = 1, this%mesh%ncell_onP
      if (this%tot_void_cell(j)) then
        this%H(j) = 0.0_r8
      else
        call pmf_eval (this%model%H_of_T, j, Tcell(j:j), this%H(j))
      end if
    end do
    
    !! Commit the initial state.
    this%seq = 0
    this%t = t
    call flush_history (this%uhist, t, this%u, udot)
    this%Hlast = this%H
    this%last_void_cell = this%void_cell
    this%last_void_face = this%void_face
    this%last_tot_void_cell = this%tot_void_cell
    this%state_is_pending = .false.
    
#ifdef GMV_DIAGNOSTICS
    call gmv_write_state (this)
#endif
  end subroutine FHT_solver_set_initial_state
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! COMPUTE_ADVECTED_HEAT
 !!
 !! Given cell temperatures, this auxillary procedure returns the heat
 !! increment for each cell due to material phase advection.  It also returns
 !! the minimum and maximum temperature of the heat parcels advected into (or
 !! remaining in) each cell.  The actual computation is done by a procedure
 !! from ADVECTION_MODULE; this procedure primarily handles the mapping of
 !! fields between old and new meshes.
 !!
  
  subroutine compute_advected_heat (this, T, dQ, Tmin, Tmax)
  
    use legacy_mesh_api, only: ncells
    use advection_module, only: compute_advected_enthalpy

    type(FHT_solver), intent(in) :: this
    real(r8), intent(in) :: T(:)
    real(r8), intent(out) :: dQ(:), Tmin(:), Tmax(:)
    
    real(r8), dimension(ncells) :: T_t, dQ_t, Tmin_t, Tmax_t
    
    T_t(:this%mesh%ncell_onP) = T
    T_t(this%mesh%ncell_onP+1:) = 0.0_r8
    call compute_advected_enthalpy (T_t, dQ_t, Tmin_t, Tmax_t)
    dQ = dQ_t(:this%mesh%ncell_onP)
    Tmin = Tmin_t(:this%mesh%ncell_onP)
    Tmax = Tmax_t(:this%mesh%ncell_onP)
    ASSERT(all(Tmin <= T .and. T <= Tmax))
    
  end subroutine compute_advected_heat
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BACKWARD_EULER_SOLVE
 !!
 !! This procedure solves the backward Euler time step system for the heat
 !! conduction problem.
 !!

  subroutine backward_euler_solve (this, t, dt, Hlast, u, stat)
  
#ifdef GMV_DIAGNOSTICS
    use unstr_mesh_gmv
#endif
  
    type(FHT_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt, Hlast(:)
    real(r8), intent(inout), target :: u(:)
    integer, intent(out) :: stat
  
    integer :: j, itr
    real(r8) :: f(size(u)), H, error
    real(r8), pointer :: Tcell(:)
    real(r8), allocatable :: Hdot(:)
  
    !! Update the preconditioner.
    this%num_precon_compute = this%num_precon_compute + 1
    call FHT_precon_compute (this%precon, t, u, dt)
    
    call FHT_model_get_cell_temp_view (this%model, u, Tcell)
    allocate(Hdot(size(Tcell)))
    
    !! Compute the initial residual and norm.
    this%num_compute_f = this%num_compute_f + 1
    do j = 1, size(Tcell)
      call pmf_eval (this%model%H_of_T, j, Tcell(j:j), H)
      Hdot(j) = (H - Hlast(j)) / dt
    end do
    call FHT_model_compute_f (this%model, t, u, Hdot, f)
    call FHT_norm_fnorm (this%norm, t, u, Hdot, f)
    
    itr = 0
    call this%accel%restart
    do ! until converged
    
      itr = itr + 1
      
      !! Compute the next solution iterate.
      this%num_precon_apply = this%num_precon_apply + 1
      call FHT_precon_apply (this%precon, t, u, f)
      call this%accel%accel_update (f)
      u = u - f
      
      !! Compute the residual and norm.
      this%num_compute_f = this%num_compute_f + 1
      do j = 1, size(Tcell)
        call pmf_eval (this%model%H_of_T, j, Tcell(j:j), H)
        Hdot(j) = (H - Hlast(j)) / dt
      end do
      call FHT_model_compute_f (this%model, t, u, Hdot, f)
      call FHT_norm_fnorm (this%norm, t, u, Hdot, f, error)
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
  
  subroutine FHT_solver_get_cell_temp_view (this, view)
    type(FHT_solver), intent(in) :: this
    real(r8), pointer :: view(:)
    call FHT_model_get_cell_temp_view (this%model, this%u, view)
  end subroutine FHT_solver_get_cell_temp_view
  
  subroutine FHT_solver_get_cell_temp_copy (this, copy)
    type(FHT_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    call FHT_model_get_cell_temp_copy (this%model, this%u, copy)
  end subroutine FHT_solver_get_cell_temp_copy
  
  subroutine FHT_solver_get_face_temp_view (this, view)
    type(FHT_solver), intent(in) :: this
    real(r8), pointer :: view(:)
    call FHT_model_get_face_temp_view (this%model, this%u, view)
  end subroutine FHT_solver_get_face_temp_view
  
  subroutine FHT_solver_get_face_temp_copy (this, copy)
    type(FHT_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    call FHT_model_get_face_temp_copy (this%model, this%u, copy)
  end subroutine FHT_solver_get_face_temp_copy
  
  subroutine FHT_solver_get_cell_heat_view (this, view)
    type(FHT_solver), intent(in) :: this
    real(r8), pointer :: view(:)
    view => this%H
  end subroutine FHT_solver_get_cell_heat_view
  
  subroutine FHT_solver_get_cell_heat_copy (this, copy)
    type(FHT_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    ASSERT(size(copy) >= size(this%H))
    copy(:size(this%H)) = this%H
  end subroutine FHT_solver_get_cell_heat_copy
  
  subroutine FHT_solver_get_void_cell_view (this, view)
    type(FHT_solver), intent(in) :: this
    logical, pointer :: view(:)
    view => this%void_cell(:size(this%H))
  end subroutine FHT_solver_get_void_cell_view
  
  subroutine FHT_solver_get_void_cell_copy (this, copy)
    type(FHT_solver), intent(in) :: this
    logical, intent(inout) :: copy(:)
    ASSERT(size(copy) >= size(this%H))
    copy(:size(this%H)) = this%void_cell(:size(this%H))
  end subroutine FHT_solver_get_void_cell_copy
  
  function FHT_solver_last_time (this) result (t)
    type(FHT_solver), intent(in) :: this
    real(r8) :: t
    t = most_recent_time(this%uhist)
  end function FHT_solver_last_time
  
  function FHT_solver_last_step_size (this) result (dt)
    type(FHT_solver), intent(in) :: this
    real(r8) :: dt
    dt = last_time_delta(this%uhist)
  end function FHT_solver_last_step_size
  
  subroutine FHT_solver_get_stepping_stats (this, counters)
    type(FHT_solver), intent(in) :: this
    integer, intent(inout) :: counters(:)
    ASSERT(size(counters) >= 4)
    counters(1) = this%seq
    counters(2) = this%num_compute_f
    counters(3) = this%num_precon_compute
    counters(4) = this%num_precon_apply
  end subroutine FHT_solver_get_stepping_stats
  
  subroutine FHT_solver_get_cell_temp_grad (this, tgrad)
    use mfd_disc_type
    type(FHT_solver), intent(in) :: this
    real(r8), intent(out) :: tgrad(:,:)
    real(r8) :: tface(this%model%mesh%nface)
    INSIST(size(tgrad,1) == 3)
    INSIST(size(tgrad,2) == this%model%mesh%ncell_onP)
    call FHT_model_get_face_temp_copy (this%model, this%u, tface)
    call gather_boundary (this%model%mesh%face_ip, tface)
    call this%model%disc%compute_cell_grad (tface, tgrad)
  end subroutine FHT_solver_get_cell_temp_grad

end module FHT_solver_type
