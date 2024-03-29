!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HTSD_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use HTSD_model_type
  use HTSD_precon_type
  use HTSD_norm_type
  use matl_mesh_func_type
  use unstr_mesh_type
  use htsd_idaesol_model_type
  use idaesol_type
  use parameter_list_type
  implicit none
  private
  
  type, public :: HTSD_solver
    type(matl_mesh_func), pointer :: mmf => null()
    type(unstr_mesh), pointer :: mesh => null()
    type(htsd_idaesol_model) :: integ_model
    type(idaesol) :: integ
    logical :: state_is_pending = .false.
    !! Pending state
    real(r8) :: t
    real(r8), pointer :: u(:) => null()
    type(HTSD_model),  pointer :: model  => null()
    type(HTSD_precon), pointer :: precon => null()
    type(HTSD_norm),   pointer :: norm   => null()
    real(r8) :: hmin
    integer :: max_step_tries
    !integer :: step_method
    type(parameter_list), pointer :: ic_params
  end type HTSD_solver
  
  public :: HTSD_solver_init
  public :: HTSD_solver_delete
  public :: HTSD_solver_set_initial_state
  public :: HTSD_solver_advance_state
  public :: HTSD_solver_commit_pending_state
  public :: HTSD_solver_restart
  public :: HTSD_solver_get_soln_view, HTSD_solver_get_soln_copy
  public :: HTSD_solver_get_cell_temp_view, HTSD_solver_get_cell_temp_copy
  public :: HTSD_solver_get_face_temp_view, HTSD_solver_get_face_temp_copy
  public :: HTSD_solver_get_cell_heat_view, HTSD_solver_get_cell_heat_copy
  public :: HTSD_solver_get_cell_conc_view, HTSD_solver_get_cell_conc_copy
  public :: HTSD_solver_get_cell_temp_grad
  public :: HTSD_solver_get_stepping_stats  
  public :: HTSD_solver_last_step_size
  public :: HTSD_solver_last_time
   
contains

  subroutine HTSD_solver_init (this, mmf, model, params)
    type(HTSD_solver), intent(out), target :: this
    type(matl_mesh_func), intent(in), target :: mmf
    type(HTSD_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer :: n, output_unit
    type(parameter_list), pointer :: plist
    logical :: verbose_stepping
    this%mmf   => mmf
    this%model => model
    this%mesh  => model%mesh
    allocate(this%norm)
    plist => params%sublist('norm')
    call HTSD_norm_init (this%norm, model, plist)
    allocate(this%precon)
    plist => params%sublist('precon')
    call HTSD_precon_init (this%precon, model, plist)
    n = HTSD_model_size(model)
    allocate(this%u(n))

    call this%integ_model%init(this%model, this%precon, this%norm)

    block
      integer :: stat
      character(:), allocatable :: errmsg
      call this%integ%init(this%integ_model, params, stat, errmsg)
      INSIST(stat == 0)
    end block

    call params%get('verbose-stepping', verbose_stepping)
    if (verbose_stepping) then
      call params%get('output-unit', output_unit)
      call this%integ%set_verbose_stepping(output_unit)
    end if
    call params%get('hmin', this%hmin)
    call params%get('max-step-tries', this%max_step_tries)
    !this%step_method = params%step_method
    !! Grab parameters for HTSD_init_cond%init
    allocate(this%ic_params)
    if (associated(model%ht)) then
      block
      real(r8) :: rval
      plist => params%sublist('norm')
      call plist%get('abs-t-tol', rval)
      call this%ic_params%set ('atol-temp', 0.01_r8 * rval)
      call plist%get('rel-t-tol', rval)
      call this%ic_params%set ('rtol-temp', 0.01_r8 * rval)
      call this%ic_params%set ('max-iter', 50)
      call this%ic_params%set ('method', 'SSOR')
      plist => this%ic_params%sublist('params')
      call plist%set ('num-cycles', 1)
      end block
    end if
  end subroutine HTSD_solver_init
  
  subroutine HTSD_solver_delete (this)
    type(HTSD_solver), intent(inout) :: this
    if (associated(this%u)) deallocate(this%u)
    if (associated(this%precon)) then
      call HTSD_precon_delete (this%precon)
      deallocate(this%precon)
    end if
    if (associated(this%norm)) deallocate(this%norm)
    !! The solver owns the void cell and face mask components of the model.
    if (associated(this%model%void_cell)) deallocate(this%model%void_cell)
    if (associated(this%model%void_face)) deallocate(this%model%void_face)
    if (associated(this%ic_params)) deallocate(this%ic_params)
  end subroutine HTSD_solver_delete
  
  subroutine HTSD_solver_get_soln_view (this, view)
    type(HTSD_solver), intent(in) :: this
    real(r8), pointer :: view(:)
    view => this%u
  end subroutine HTSD_solver_get_soln_view
  
  subroutine HTSD_solver_get_soln_copy (this, copy)
    type(HTSD_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    ASSERT(size(copy) >= size(this%u))
    copy(:size(this%u)) = this%u
  end subroutine HTSD_solver_get_soln_copy
  
  subroutine HTSD_solver_get_cell_temp_view (this, view)
    type(HTSD_solver), intent(in) :: this
    real(r8), pointer :: view(:)
    call HTSD_model_get_cell_temp_view (this%model, this%u, view)
  end subroutine HTSD_solver_get_cell_temp_view
  
  subroutine HTSD_solver_get_cell_temp_copy (this, copy)
    type(HTSD_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    call HTSD_model_get_cell_temp_copy (this%model, this%u, copy)
  end subroutine HTSD_solver_get_cell_temp_copy
  
  subroutine HTSD_solver_get_face_temp_view (this, view)
    type(HTSD_solver), intent(in) :: this
    real(r8), pointer :: view(:)
    call HTSD_model_get_face_temp_view (this%model, this%u, view)
  end subroutine HTSD_solver_get_face_temp_view
  
  subroutine HTSD_solver_get_face_temp_copy (this, copy)
    type(HTSD_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    call HTSD_model_get_face_temp_copy (this%model, this%u, copy)
  end subroutine HTSD_solver_get_face_temp_copy
  
  subroutine HTSD_solver_get_cell_heat_view (this, view)
    type(HTSD_solver), intent(in) :: this
    real(r8), pointer :: view(:)
    call HTSD_model_get_cell_heat_view (this%model, this%u, view)
  end subroutine HTSD_solver_get_cell_heat_view
  
  subroutine HTSD_solver_get_cell_heat_copy (this, copy)
    type(HTSD_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    call HTSD_model_get_cell_heat_copy (this%model, this%u, copy)
  end subroutine HTSD_solver_get_cell_heat_copy
  
  subroutine HTSD_solver_get_cell_conc_view (this, n, view)
    type(HTSD_solver), intent(in) :: this
    integer, intent(in) :: n
    real(r8), pointer :: view(:)
    call HTSD_model_get_cell_conc_view (this%model, n, this%u, view)
  end subroutine HTSD_solver_get_cell_conc_view
  
  subroutine HTSD_solver_get_cell_conc_copy (this, n, copy)
    type(HTSD_solver), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(inout) :: copy(:)
    call HTSD_model_get_cell_conc_copy (this%model, n, this%u, copy)
  end subroutine HTSD_solver_get_cell_conc_copy
  
  subroutine HTSD_solver_get_cell_temp_grad (this, tgrad)
    use mfd_disc_type
    type(HTSD_solver), intent(in) :: this
    real(r8), intent(out) :: tgrad(:,:)
    real(r8) :: tface(this%model%mesh%nface)
    INSIST(size(tgrad,1) == 3)
    INSIST(size(tgrad,2) == this%model%mesh%ncell_onP)
    INSIST(associated(this%model%ht))
    call HTSD_model_get_face_temp_copy (this%model, this%u, tface)
    call this%model%mesh%face_imap%gather_offp(tface)
    call this%model%disc%compute_cell_grad (tface, tgrad)
  end subroutine HTSD_solver_get_cell_temp_grad
    
  subroutine HTSD_solver_get_stepping_stats (this, counters)
    type(HTSD_solver), intent(in) :: this
    integer, intent(out) :: counters(:)
    ASSERT(size(counters) == 6)
    call this%integ%get_stepping_statistics(counters)
  end subroutine HTSD_solver_get_stepping_stats
  
  function HTSD_solver_last_time (this) result (t)
    type(HTSD_solver), intent(in) :: this
    real(r8) :: t
    t = this%integ%last_time()
  end function HTSD_solver_last_time
  
  function HTSD_solver_last_step_size (this) result (h)
    type(HTSD_solver), intent(in) :: this
    real(r8) :: h
    h = this%integ%last_step_size()
  end function HTSD_solver_last_step_size
  
  subroutine HTSD_solver_advance_state(this, t, hnext, stat)
    type(HTSD_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hnext
    integer, intent(out) :: stat
    call this%integ%step(t, this%u, hnext, stat)
    if (stat == 0) then
      this%t = t
      this%state_is_pending = .true.
    else ! failed -- restore last good state before returning
      call this%integ%get_last_state_copy(this%u)
      this%state_is_pending = .false.
    end if
  end subroutine HTSD_solver_advance_state

  subroutine HTSD_solver_commit_pending_state (this)
    type(HTSD_solver), intent(inout) :: this
    INSIST(this%state_is_pending)
    call this%integ%commit_state(this%t, this%u)
    this%state_is_pending = .false.
  end subroutine HTSD_solver_commit_pending_state
  
  subroutine HTSD_solver_set_initial_state (this, t, dt, temp, conc)
  
    use parallel_communication, only: global_count
    use HTSD_init_cond_type
  
    type(HTSD_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional :: temp(:), conc(:,:)
    
    integer :: j
    real(r8), allocatable :: void_vol_frac(:), udot(:)
    type(HTSD_init_cond) :: ic
    
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
    
    allocate(udot(size(this%u)))
    call this%ic_params%set ('dt', dt)
    call ic%init (this%model, this%ic_params)
    call ic%compute (t, temp, conc, u=this%u, udot=udot)
    call this%integ%set_initial_state(t, this%u, udot)
    deallocate(udot)
    
  end subroutine HTSD_solver_set_initial_state

  subroutine HTSD_solver_restart (this, dt)

    use HTSD_init_cond_type

    type(HTSD_solver), intent(inout) :: this
    real(r8), intent(in) :: dt

    real(r8) :: udot(size(this%u))
    type(HTSD_init_cond) :: ic

    INSIST(associated(this%model))

    call this%ic_params%set ('dt', dt)
    call ic%init (this%model, this%ic_params)
    call ic%compute_udot (this%t, this%u, udot)
    call this%integ%set_initial_state(this%t, this%u, udot)
    !call set_initial_state (this%bdf2_state, this%t, this%u, udot)

  end subroutine HTSD_solver_restart

end module HTSD_solver_type
