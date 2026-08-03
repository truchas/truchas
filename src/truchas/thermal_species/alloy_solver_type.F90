!!
!! ALLOY_SOLVER_TYPE
!!
!! Concrete thermal/phase-change solver for multicomponent alloy
!! solidification.  Solute concentrations are explicit transport state while
!! liquid fraction, enthalpy, temperature, and face temperature are solved
!! implicitly by the common IDAESOL integrator.
!!

#include "f90_assert.fpp"

module alloy_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use thermal_species_solver_class
  use alloy_vector_type
  use alloy_model_type
  use alloy_precon_type
  use alloy_norm_type
  use alloy_idaesol_model_type
  use matl_mesh_func_type
  use unstr_mesh_type
  use new_idaesol_type
  use parameter_list_type
  implicit none
  private

  type, extends(thermal_species_solver), public :: alloy_solver
    type(matl_mesh_func), pointer :: mmf => null() ! unowned reference
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(alloy_idaesol_model) :: integ_model
    type(idaesol) :: integ
    logical :: state_is_pending = .false.
    real(r8) :: t
    type(alloy_vector) :: u
    type(alloy_model) :: model
    type(alloy_precon) :: precon
    type(alloy_norm) :: norm
    type(parameter_list), pointer :: ic_params => null()
    real(r8), pointer :: clast(:,:) => null(), c(:,:) => null(), cdot(:,:) => null()
    integer :: num_comp = 0
  contains
    procedure :: init
    procedure :: step
    procedure :: commit_step
    procedure :: get_cell_heat_copy
    procedure :: get_cell_temp_copy
    procedure :: get_face_temp_copy
    procedure :: get_face_temp_view
    procedure :: get_liq_frac_view
    procedure :: get_cell_conc_copy
    procedure :: get_cell_temp_grad
    procedure :: get_stepping_stats
    procedure :: last_step_size
    procedure :: last_time
    procedure :: set_initial_state
    procedure :: set_solver_initial_state
    procedure :: restart
    procedure :: log_step_stats
    procedure :: set_ext_enthalpy_rate
    procedure :: set_ext_species_rate
  end type

contains

  subroutine get_liq_frac_view(this, view)
    class(alloy_solver), intent(in), target :: this
    real(r8), pointer, intent(out) :: view(:)
    view => this%u%lf
  end subroutine get_liq_frac_view

  subroutine init(this, mesh, mmf, params, stat, errmsg)

    use material_model_driver, only: matl_model
    use material_class, only: material
    use parallel_communication, only: is_IOP
    use truchas_env, only: output_file_name

    class(alloy_solver), intent(inout), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    class(material), pointer :: matl
    type(parameter_list), pointer :: plist
    character(:), allocatable :: name
    real(r8), allocatable :: c0(:)
    integer :: j, mid, lun
    logical :: verbose_stepping

    this%mmf => mmf
    this%mesh => mesh

    call params%get('material', name, stat, errmsg)
    if (stat /= 0) return
    mid = matl_model%matl_index(name)
    if (mid == 0) then
      stat = 1
      errmsg = 'unknown alloy material: ' // name
      return
    end if
    call matl_model%get_matl_ref(mid, matl)
    call this%model%init(mesh, matl, params, stat, errmsg)
    if (stat /= 0) return

    this%num_comp = this%model%num_comp
    allocate(this%c(this%num_comp, mesh%ncell), this%clast(this%num_comp, mesh%ncell), &
        this%cdot(this%num_comp, mesh%ncell))
    call params%get('concentration', c0, stat, errmsg)
    if (stat /= 0) return
    if (size(c0) /= this%num_comp .or. any(c0 < 0.0_r8)) then
      stat = 1
      errmsg = 'invalid initial alloy concentration'
      return
    end if
    do j = 1, mesh%ncell
      this%c(:,j) = c0
    end do
    this%clast = this%c
    this%cdot = 0.0_r8

    plist => params%sublist('norm')
    call this%norm%init(this%model, plist, stat, errmsg)
    if (stat /= 0) return
    plist => params%sublist('precon')
    call this%precon%init(this%model, plist, stat, errmsg)
    if (stat /= 0) return

    call this%model%init_vector(this%u)
    call this%integ_model%init(this%model, this%precon, this%norm, this%c, this%cdot)
    call this%integ%init(this%integ_model, params, stat, errmsg)
    if (stat /= 0) return

    call params%get('verbose-stepping', verbose_stepping)
    if (verbose_stepping) then
      lun = -1
      if (is_IOP) open(newunit=lun, file=output_file_name('bdf2.out'), position='rewind', action='write')
      call params%set('output-unit', lun)
      if (is_IOP) call this%integ%set_verbose_stepping(lun)
    end if

    block
      real(r8) :: rval
      allocate(this%ic_params)
      plist => params%sublist('norm')
      call plist%get('abs-t-tol', rval)
      call this%ic_params%set('atol-temp', 0.01_r8*rval)
      call plist%get('rel-t-tol', rval)
      call this%ic_params%set('rtol-temp', 0.01_r8*rval)
      call this%ic_params%set('max-iter', 50)
      call this%ic_params%set('method', 'SSOR')
      plist => this%ic_params%sublist('params')
      call plist%set('num-cycles', 1)
    end block

  end subroutine

  subroutine step(this, t, hnext, errc)
    class(alloy_solver), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hnext
    integer, intent(out) :: errc
    this%c = this%clast + (t - this%integ%last_time())*this%cdot
    call this%integ%step(t, this%u, hnext, errc)
    if (errc == 0) then
      this%t = t
      this%state_is_pending = .true.
    else
      this%c = this%clast
      call this%integ%get_last_state_copy(this%u)
      this%state_is_pending = .false.
    end if
  end subroutine

  subroutine commit_step(this)
    class(alloy_solver), intent(inout) :: this
    if (this%state_is_pending) then
      this%clast = this%c
      call this%integ%commit_state(this%t, this%u)
      this%state_is_pending = .false.
    end if
  end subroutine

  subroutine set_initial_state(this, t, dt, temp, conc)
    class(alloy_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional :: temp(:), conc(:,:)
    call set_solver_initial_state(this, t, dt, temp, conc)
  end subroutine

  subroutine set_solver_initial_state(this, t, dt, temp, conc)
    use alloy_ic_solver_type
    class(alloy_solver), intent(inout), target :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), optional :: temp(:), conc(:,:)
    type(alloy_vector) :: udot
    type(alloy_ic_solver) :: ic
    ASSERT(present(temp))
    if (present(conc)) then
      ASSERT(all(shape(conc) == [this%mesh%ncell_onP, this%num_comp]))
      this%c(:, :this%mesh%ncell_onP) = transpose(conc)
      this%clast = this%c
    end if
    this%t = t
    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute(t, this%c, this%cdot, temp, this%u, udot)
    call this%integ%set_initial_state(t, this%u, udot)
  end subroutine

  subroutine restart(this, dt)
    use alloy_ic_solver_type
    class(alloy_solver), intent(inout) :: this
    real(r8), intent(in) :: dt
    type(alloy_vector) :: udot
    type(alloy_ic_solver) :: ic
    call udot%init(this%u)
    call this%ic_params%set('dt', dt)
    call ic%init(this%model, this%ic_params)
    call ic%compute_udot(this%t, this%c, this%cdot, this%u, udot)
    call this%integ%set_initial_state(this%t, this%u, udot)
  end subroutine

  subroutine get_cell_temp_copy(this, copy)
    class(alloy_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    copy(:min(size(copy), this%mesh%ncell_onP)) = this%u%tc(:min(size(copy), this%mesh%ncell_onP))
  end subroutine

  subroutine get_cell_heat_copy(this, copy)
    class(alloy_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    copy(:min(size(copy), this%mesh%ncell_onP)) = this%u%hc(:min(size(copy), this%mesh%ncell_onP))
  end subroutine

  subroutine get_face_temp_copy(this, copy)
    class(alloy_solver), intent(in) :: this
    real(r8), intent(inout) :: copy(:)
    copy(:min(size(copy), this%mesh%nface_onP)) = this%u%tf(:min(size(copy), this%mesh%nface_onP))
  end subroutine

  subroutine get_face_temp_view(this, view)
    class(alloy_solver), intent(in), target :: this
    real(r8), pointer :: view(:)
    view => this%u%tf(:this%mesh%nface_onP)
  end subroutine

  subroutine get_cell_conc_copy(this, n, copy)
    class(alloy_solver), intent(in) :: this
    integer, intent(in) :: n
    real(r8), intent(inout) :: copy(:)
    ASSERT(n >= 1 .and. n <= this%num_comp)
    copy(:min(size(copy), this%mesh%ncell_onP)) = this%c(n, :min(size(copy), this%mesh%ncell_onP))
  end subroutine

  subroutine get_cell_temp_grad(this, tgrad)
    class(alloy_solver), intent(inout) :: this
    real(r8), intent(out) :: tgrad(:,:)
    call this%mesh%face_imap%gather_offp(this%u%tf)
    call this%model%disc%compute_cell_grad(this%u%tf, tgrad)
  end subroutine

  subroutine set_ext_enthalpy_rate(this, enthalpy_rate)
    class(alloy_solver), intent(inout) :: this
    real(r8), intent(in) :: enthalpy_rate(:)
    call this%model%set_heat_source(enthalpy_rate)
  end subroutine

  subroutine set_ext_species_rate(this, n, species_rate)
    class(alloy_solver), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: species_rate(:)
    ASSERT(n >= 1 .and. n <= this%num_comp)
    ASSERT(size(species_rate) == this%mesh%ncell_onP)
    this%cdot(n, :this%mesh%ncell_onP) = species_rate
  end subroutine

  function last_time(this) result(t)
    class(alloy_solver), intent(in) :: this
    real(r8) :: t
    t = this%integ%last_time()
  end function

  function last_step_size(this) result(h)
    class(alloy_solver), intent(in) :: this
    real(r8) :: h
    h = this%integ%last_step_size()
  end function

  subroutine get_stepping_stats(this, counters)
    class(alloy_solver), intent(in) :: this
    integer, intent(out) :: counters(:)
    call this%integ%get_stepping_statistics(counters)
  end subroutine

  subroutine log_step_stats(this)
    use truchas_logging_services, only: TLS_info
    class(alloy_solver), intent(in) :: this
    integer :: counters(6)
    character(len=96) :: message
    call this%get_stepping_stats(counters)
    write(message, 1) this%last_step_size(), counters(1:2), counters(4:6)
    1 format('ALLOY: dt=',es10.3,', NRES:NPC=',i7.7,':',i5.5,', NNR:NNF:NSR=',3(i4.4,:,':'))
    call TLS_info(message)
  end subroutine

end module alloy_solver_type
