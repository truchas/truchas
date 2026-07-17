!!
!! SD_SOLVER_TYPE
!!
!! This module defines the concrete solver for standalone species transport.
!! It owns the species transport model, preconditioner, norm, integrator
!! adapter, integrator, and current solution vector, and coordinates
!! initialization, time stepping, restart, and access to species state data
!! through the common thermal_species_solver interface.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module sd_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use thermal_species_solver_class
  use sd_vector_type
  use sd_model_type
  use sd_precon_type
  use sd_norm_type
  use matl_mesh_func_type
  use unstr_mesh_type
  use sd_idaesol_model_type
  use new_idaesol_type
  use parameter_list_type
  implicit none
  private

  type, public, extends(thermal_species_solver) :: sd_solver
    type(matl_mesh_func), pointer :: mmf => null() ! unowned reference
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(sd_idaesol_model) :: integ_model
    type(idaesol) :: integ
    logical :: state_is_pending = .false.
    real(r8) :: t
    type(sd_vector) :: u
    type(sd_model) :: model
    logical :: model_is_initialized = .false.
    type(sd_precon) :: precon
    type(sd_norm) :: norm
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
    procedure :: get_cell_conc_copy
    procedure :: set_ext_species_rate
  end type

contains

  subroutine init(this, mesh, mmf, params, stat, errmsg)

    use parallel_communication, only: is_IOP
    use truchas_env, only: output_file_name

    !! TARGET is needed because assemble_solver initializes persistent
    !! pointer references to components of this solver.
    class(sd_solver), intent(inout), target :: this
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
    class(sd_solver), intent(inout), target :: this
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

    call this%u%init(this%mesh, this%model%num_comp)

    call this%integ_model%init(this%model, this%precon, this%norm)

    call this%integ%init(this%integ_model, params, stat, errmsg)
    if (stat /= 0) return

    call params%get('verbose-stepping', verbose_stepping)
    if (is_IOP .and. verbose_stepping) then
      call params%get('output-unit', output_unit)
      call this%integ%set_verbose_stepping(output_unit)
    end if
    !! Grab parameters for sd_ic_solver%init
    block
      real(r8), allocatable :: rval(:)
      plist => params%sublist('norm')
      call plist%get('abs-c-tol', rval)
      call this%ic_params%set('atol-conc', 0.01_r8 * rval)
      call plist%get('rel-c-tol', rval)
      call this%ic_params%set('rtol-conc', 0.01_r8 * rval)
      call this%ic_params%set('max-iter', 50)
      call this%ic_params%set('method', 'SSOR')
      plist => this%ic_params%sublist('params')
      call plist%set('num-cycles', 1)
    end block

  end subroutine assemble_solver

  subroutine step(this, t, hnext, errc)
    class(sd_solver), intent(inout) :: this
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
  end subroutine

  subroutine commit_step(this)
    class(sd_solver), intent(inout) :: this
    if (this%state_is_pending) then
      call this%integ%commit_state(this%t, this%u)
      this%state_is_pending = .false.
    end if
  end subroutine

  subroutine restart(this, dt)
    class(sd_solver), intent(inout) :: this
    real(r8), intent(in) :: dt
    call restart_target(this, dt)
  end subroutine restart

  !! Keep this helper separate from the type-bound restart override. The
  !! abstract restart interface does not have TARGET on THIS, but the local
  !! IC solver initialized here stores a reference to THIS%model. NAG dangling
  !! pointer checking requires the dummy that supplies that target to have the
  !! TARGET attribute.
  subroutine restart_target(this, dt)
    use sd_ic_solver_type
    class(sd_solver), intent(inout), target :: this
    real(r8), intent(in) :: dt
    type(sd_vector) :: udot
    type(sd_ic_solver) :: ic
    INSIST(this%model_is_initialized)
    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute_udot(this%t, this%u, udot)
    call this%integ%set_initial_state(this%t, this%u, udot)
  end subroutine restart_target

  subroutine set_initial_state(this, t, dt, temp, conc)
    class(sd_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional :: temp(:), conc(:,:)
    ASSERT(present(conc))
    call this%set_solver_initial_state(t, conc, dt)
  end subroutine set_initial_state

  subroutine set_solver_initial_state(this, t, conc, dt)

    use sd_ic_solver_type

    class(sd_solver), intent(inout), target :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: conc(:,:)
    real(r8), intent(in) :: dt

    type(sd_vector) :: udot
    type(sd_ic_solver) :: ic

    INSIST(this%model_is_initialized)

    !! Static void masks are model-owned and defined once from the initial MMF.
    call this%model%define_void_masks(this%mmf)

    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute(t, conc, this%u, udot)
    call this%integ%set_initial_state(t, this%u, udot)

  end subroutine set_solver_initial_state


  function last_time(this) result(t)
    class(sd_solver), intent(in) :: this
    real(r8) :: t
    t = this%integ%last_time()
  end function

  function last_step_size(this) result(h)
    class(sd_solver), intent(in) :: this
    real(r8) :: h
    h = this%integ%last_step_size()
  end function

  subroutine get_stepping_stats(this, counters)
    class(sd_solver), intent(in) :: this
    integer, intent(out) :: counters(:)
    ASSERT(size(counters) == 6)
    call this%integ%get_stepping_statistics(counters)
  end subroutine

  subroutine log_step_stats(this)
    use truchas_logging_services, only: TLS_info
    class(sd_solver), intent(in) :: this
    real(r8) :: hlast
    integer :: counters(6)
    character(len=80) :: message(2)
    hlast = this%last_step_size()
    call this%get_stepping_stats(counters)
    write(message,1) hlast, counters(1:2), counters(4:6)
    1 format(/,'SD: dt=',es9.3,', NRES:NPC=',i7.7,':',i5.5,', NNR:NNF:NSR=',3(i4.4,:,':'))
    call TLS_info(message)
  end subroutine

  subroutine get_cell_conc_copy(this, n, copy)
    class(sd_solver), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(inout) :: copy(:)
    integer :: ncopy
    ASSERT(n > 0 .and. n <= this%model%num_comp)
    !! External copies expose only on-process entries; the off-process vector
    !! entries are solver scratch space.
    ncopy = min(size(copy), this%mesh%ncell_onP)
    copy(:ncopy) = this%u%cc(:ncopy,n)
  end subroutine

  subroutine set_ext_species_rate(this, n, species_rate)
    class(sd_solver), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: species_rate(:)
    call this%model%set_ext_species_rate(n, species_rate)
  end subroutine

end module sd_solver_type
