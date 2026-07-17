!!
!! HTSD_SOLVER_TYPE
!!
!! This module defines the concrete solver for coupled thermal and species
!! transport. It owns the coupled transport model, preconditioner, norm,
!! integrator adapter, integrator, and current solution vector, and coordinates
!! initialization, time stepping, restart, and access to thermal and species
!! state data through the common thermal_species_solver interface.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module htsd_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use thermal_species_solver_class
  use htsd_vector_type
  use htsd_model_type
  use htsd_precon_type
  use htsd_norm_type
  use matl_mesh_func_type
  use unstr_mesh_type
  use htsd_idaesol_model_type
  use new_idaesol_type
  use parameter_list_type
  implicit none
  private
  
  type, public, extends(thermal_species_solver) :: htsd_solver
    type(matl_mesh_func), pointer :: mmf => null() ! unowned reference
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(htsd_idaesol_model) :: integ_model
    type(idaesol) :: integ
    logical :: state_is_pending = .false.
    !! Pending state
    real(r8) :: t
    type(htsd_vector) :: u
    type(htsd_model) :: model
    logical :: model_is_initialized = .false.
    type(htsd_precon) :: precon
    type(htsd_norm)   :: norm
    type(parameter_list) :: ic_params
  contains
    procedure :: init
    procedure :: step
    procedure :: commit_step
    procedure :: restart
    procedure :: set_initial_state
    procedure :: set_solver_initial_state
    procedure :: last_time
    procedure :: last_step_size
    procedure :: get_stepping_stats
    procedure :: log_step_stats
    procedure :: get_cell_temp_copy
    procedure :: get_cell_heat_copy
    procedure :: get_cell_conc_copy
    procedure :: get_face_temp_copy
    procedure :: get_face_temp_view
    procedure :: get_cell_temp_grad
    procedure :: set_ext_enthalpy_rate
    procedure :: set_ext_species_rate
    procedure :: update_moving_vf
    procedure :: add_moving_vf_events
  end type htsd_solver
  
contains

  subroutine init(this, mesh, mmf, params, stat, errmsg)

    use parallel_communication, only: is_IOP
    use truchas_env, only: output_file_name

    !! TARGET is needed because assemble_solver initializes persistent
    !! pointer references to components of this solver.
    class(htsd_solver), intent(inout), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: solver_params, model_params
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
      call solver_params%set('output-unit', lun)
    end if

    call assemble_solver(this, mmf, solver_params, stat, errmsg)

  end subroutine init

  subroutine assemble_solver(this, mmf, params, stat, errmsg)
    use parallel_communication, only: is_IOP

    !! TARGET allows IDAESOL and its model adapter to retain references to
    !! solver components such as model, precon, norm, and integ_model.
    class(htsd_solver), intent(inout), target :: this
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer :: output_unit
    type(parameter_list), pointer :: plist
    logical :: verbose_stepping

    ASSERT(this%model_is_initialized)

    this%mmf   => mmf
    this%mesh  => this%model%mesh
    plist => params%sublist('norm')
    call this%norm%init(this%model, plist, stat, errmsg)
    if (stat /= 0) return
    plist => params%sublist('precon')
    call this%precon%init(this%model, plist, stat, errmsg)
    if (stat /= 0) return
    call this%model%init_vector(this%u)

    call this%integ_model%init(this%model, this%precon, this%norm)

    call this%integ%init(this%integ_model, params, stat, errmsg)
    if (stat /= 0) return

    call params%get('verbose-stepping', verbose_stepping)
    if (is_IOP .and. verbose_stepping) then
      call params%get('output-unit', output_unit)
      call this%integ%set_verbose_stepping(output_unit)
    end if
    !! Grab parameters for htsd_ic_solver%init
    block
      real(r8) :: rval
      real(r8), allocatable :: cval(:)
      plist => params%sublist('norm')
      call plist%get('abs-t-tol', rval)
      call this%ic_params%set('atol-temp', 0.01_r8 * rval)
      call plist%get('rel-t-tol', rval)
      call this%ic_params%set('rtol-temp', 0.01_r8 * rval)
      call plist%get('abs-c-tol', cval)
      call this%ic_params%set('atol-conc', 0.01_r8 * cval)
      call plist%get('rel-c-tol', cval)
      call this%ic_params%set('rtol-conc', 0.01_r8 * cval)
      call this%ic_params%set('max-iter', 50)
      call this%ic_params%set('method', 'SSOR')
      plist => this%ic_params%sublist('params')
      call plist%set('num-cycles', 1)
    end block

  end subroutine assemble_solver
  
  subroutine get_cell_temp_copy(this, copy)
    class(htsd_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    !! External copies and views expose only on-process entries; the
    !! off-process vector entries are solver scratch space.
    n = min(size(copy), this%mesh%ncell_onP)
    copy(:n) = this%u%tc(:n)
  end subroutine get_cell_temp_copy
  
  subroutine get_face_temp_view(this, view)
    class(htsd_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    !! External copies and views expose only on-process entries; the
    !! off-process vector entries are solver scratch space.
    view => this%u%tf(:this%mesh%nface_onP)
  end subroutine get_face_temp_view
  
  subroutine get_face_temp_copy(this, copy)
    class(htsd_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    !! External copies and views expose only on-process entries; the
    !! off-process vector entries are solver scratch space.
    n = min(size(copy), this%mesh%nface_onP)
    copy(:n) = this%u%tf(:n)
  end subroutine get_face_temp_copy
  
  subroutine get_cell_heat_copy(this, copy)
    class(htsd_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    integer :: n
    !! External copies expose only on-process entries; the off-process vector
    !! entries are solver scratch space.
    n = min(size(copy), this%mesh%ncell_onP)
    copy(:n) = this%u%hc(:n)
  end subroutine get_cell_heat_copy
  
  subroutine get_cell_conc_copy(this, n, copy)
    class(htsd_solver), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(inout) :: copy(:)
    integer :: ncopy
    ASSERT(n > 0 .and. n <= this%model%num_comp)
    !! External copies expose only on-process entries; the off-process vector
    !! entries are solver scratch space.
    ncopy = min(size(copy), this%mesh%ncell_onP)
    copy(:ncopy) = this%u%cc(:ncopy,n)
  end subroutine get_cell_conc_copy
  
  subroutine get_cell_temp_grad(this, tgrad)
    class(htsd_solver), intent(inout) :: this
    real(r8), intent(out) :: tgrad(:,:)
    INSIST(size(tgrad,1) == 3)
    INSIST(size(tgrad,2) == this%model%mesh%ncell_onP)
    !! This getter needs off-process face temperatures for the reconstruction.
    !! They are gathered directly into the solver's scratch entries.
    call this%model%mesh%face_imap%gather_offp(this%u%tf)
    call this%model%disc%compute_cell_grad(this%u%tf, tgrad)
  end subroutine get_cell_temp_grad
    
  subroutine get_stepping_stats(this, counters)
    class(htsd_solver), intent(in) :: this
    integer, intent(out) :: counters(:)
    ASSERT(size(counters) == 6)
    call this%integ%get_stepping_statistics(counters)
  end subroutine get_stepping_stats
  
  function last_time(this) result(t)
    class(htsd_solver), intent(in) :: this
    real(r8) :: t
    t = this%integ%last_time()
  end function last_time
  
  function last_step_size(this) result(h)
    class(htsd_solver), intent(in) :: this
    real(r8) :: h
    h = this%integ%last_step_size()
  end function last_step_size
  
  subroutine step(this, t, hnext, errc)
    class(htsd_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hnext
    integer, intent(out) :: errc
    call this%integ%step(t, this%u, hnext, errc)
    if (errc == 0) then
      this%t = t
      this%state_is_pending = .true.
    else
      call this%integ%get_last_state_copy(this%u)
      this%state_is_pending = .false.
    end if
  end subroutine step

  subroutine commit_step(this)
    class(htsd_solver), intent(inout) :: this
    if (this%state_is_pending) then
      call this%integ%commit_state(this%t, this%u)
      this%state_is_pending = .false.
    end if
  end subroutine commit_step

  subroutine set_initial_state(this, t, dt, temp, conc)
    class(htsd_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional :: temp(:), conc(:,:)
    ASSERT(present(temp))
    ASSERT(present(conc))
    call this%set_solver_initial_state(t, temp, conc, dt)
  end subroutine set_initial_state

  subroutine set_solver_initial_state(this, t, temp, conc, dt)

    use htsd_ic_solver_type

    class(htsd_solver), intent(inout), target :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: temp(:), conc(:,:)
    real(r8), intent(in) :: dt
    
    type(htsd_vector) :: udot
    type(htsd_ic_solver) :: ic

    INSIST(this%model_is_initialized)
    call this%model%set_initial_time(t)

    !! Static void masks are model-owned and defined once from the initial MMF.
    call this%model%define_void_masks(this%mmf)
    
    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute(t, temp, conc, u=this%u, udot=udot)
    call this%integ%set_initial_state(t, this%u, udot)
    
  end subroutine set_solver_initial_state

  subroutine restart(this, dt)
    class(htsd_solver), intent(inout) :: this
    real(r8), intent(in) :: dt
    call restart_target(this, dt)
  end subroutine restart

  !! Keep this helper separate from the type-bound restart override. The
  !! abstract restart interface does not have TARGET on THIS, but the local
  !! IC solver initialized here stores a reference to THIS%model. NAG dangling
  !! pointer checking requires the dummy that supplies that target to have the
  !! TARGET attribute.
  subroutine restart_target(this, dt)

    use htsd_ic_solver_type

    class(htsd_solver), intent(inout), target :: this
    real(r8), intent(in) :: dt

    type(htsd_vector) :: udot
    type(htsd_ic_solver) :: ic

    INSIST(this%model_is_initialized)

    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call udot%init(this%u)
    call ic%compute_udot(this%t, this%u, udot)
    call this%integ%set_initial_state(this%t, this%u, udot)

  end subroutine restart_target

  subroutine log_step_stats(this)
    use truchas_logging_services, only: TLS_info
    class(htsd_solver), intent(in) :: this
    real(r8) :: hlast
    integer :: counters(6)
    character(len=80) :: message(2)
    hlast = this%last_step_size()
    call this%get_stepping_stats(counters)
    write(message,1) hlast, counters(1:2), counters(4:6)
    1 format(/,'DS: dt=',es9.3,', NRES:NPC=',i7.7,':',i5.5,', NNR:NNF:NSR=',3(i4.4,:,':'))
    call TLS_info(message)
  end subroutine log_step_stats

  subroutine set_ext_enthalpy_rate(this, enthalpy_rate)
    class(htsd_solver), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)
    call this%model%set_ext_enthalpy_rate(enthalpy_rate)
  end subroutine set_ext_enthalpy_rate

  subroutine set_ext_species_rate(this, n, species_rate)
    class(htsd_solver), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: species_rate(:)
    call this%model%set_ext_species_rate(n, species_rate)
  end subroutine set_ext_species_rate

  subroutine update_moving_vf(this)
    class(htsd_solver), intent(inout) :: this
    call this%model%update_moving_vf
  end subroutine update_moving_vf

  subroutine add_moving_vf_events(this, eventq, rank)
    use sim_event_queue_type, only: sim_event_queue
    class(htsd_solver), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    integer, intent(in), optional :: rank
    call this%model%add_moving_vf_events(eventq, rank)
  end subroutine add_moving_vf_events

end module htsd_solver_type
