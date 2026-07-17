!!
!! THERMAL_VIEW_FACTOR_COUPLING_TYPE
!!
!! This module defines the adapter that couples the enthalpy transport models
!! to view-factor radiation problems. It owns the heat-transport-specific
!! boundary-condition bookkeeping and exposes only the radiation operations
!! needed by the HT, HTSD, and FHT residual and preconditioner assembly code.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, July 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module thermal_view_factor_coupling_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use mfd_diff_matrix_type
  use parameter_list_type
  use rad_problem_type
  use thermal_component_type
  use unstr_mesh_type
  implicit none
  private

  public :: VFR_JAC, VFR_FGS, VFR_BGS, VFR_FAC
  public :: vfr_precon_coupling_from_string

  !! Methods of coupling heat conduction and radiosity preconditioning.
  enum, bind(c)
    enumerator :: VFR_JAC = 1 ! Jacobi (radiosity and conduction decoupled)
    enumerator :: VFR_FGS = 2 ! Forward Gauss-Seidel (radiosity, then conduction)
    enumerator :: VFR_BGS = 3 ! Backward Gauss-Seidel (conduction, then radiosity)
    enumerator :: VFR_FAC = 4 ! Factorization (radiosity, conduction, radiosity)
  end enum

  type :: view_factor_enclosure_coupling
    type(rad_problem) :: problem
    logical, allocatable :: heat_flux_mask(:)
    integer, allocatable :: faces(:)
  contains
    procedure, private :: init
    procedure, private :: suppress_heat_flux_on_faces
    procedure, public :: size => coupling_size
    procedure, public :: solve_radiosity
    procedure, public :: compute_radiosity_residual
    procedure, public :: add_heat_flux_residual
    procedure, public :: add_heat_precon_deriv
    procedure, public :: apply_rad_precon
    procedure, public :: apply_rad_precon_matvec1
    procedure, public :: compute_rhs
    procedure, public :: add_rhs_deriv_times_face_vector
    procedure, public :: add_heat_flux_to_residual
    procedure, public :: set_initial_time
    procedure, public :: update_moving_vf
    procedure, public :: add_moving_vf_events
  end type view_factor_enclosure_coupling

  type, public :: thermal_view_factor_coupling
    private
    type(view_factor_enclosure_coupling), allocatable :: encl(:)
    integer, allocatable :: offset(:)
  contains
    procedure, public :: init => coupling_init
    procedure, public :: validate_bc
    procedure, public :: is_active
    procedure, public :: size => total_size
    procedure, private :: num_enclosures
    procedure, private :: enclosure_range
    procedure, public :: solve_radiosity => coupling_solve_radiosity
    procedure, public :: compute_radiosity_residual => coupling_compute_radiosity_residual
    procedure, public :: relative_residual_norm => coupling_relative_residual_norm
    procedure, public :: add_heat_flux_residual => coupling_add_heat_flux_residual
    procedure, public :: add_heat_precon_deriv => coupling_add_heat_precon_deriv
    procedure, public :: apply_rad_precon => coupling_apply_rad_precon
    procedure, public :: apply_rad_precon_matvec1 => coupling_apply_rad_precon_matvec1
    procedure, public :: compute_rhs => coupling_compute_rhs
    procedure, public :: add_rhs_deriv_times_face_vector => coupling_add_rhs_deriv_times_face_vector
    procedure, public :: add_heat_flux_to_residual => coupling_add_heat_flux_to_residual
    procedure, public :: set_initial_time => coupling_set_initial_time
    procedure, public :: update_moving_vf => coupling_update_moving_vf
    procedure, public :: add_moving_vf_events => coupling_add_moving_vf_events
  end type thermal_view_factor_coupling

contains

  integer function vfr_precon_coupling_from_string(method) result(coupling)
    character(*), intent(in) :: method

    select case (method)
    case ('jacobi')
      coupling = VFR_JAC
    case ('forward-gs')
      coupling = VFR_FGS
    case ('backward-gs')
      coupling = VFR_BGS
    case ('factorization')
      coupling = VFR_FAC
    case default
      INSIST(.false.)
      coupling = VFR_BGS
    end select
  end function

  subroutine init(this, mesh, name, params)
    class(view_factor_enclosure_coupling), intent(out) :: this
    type(unstr_mesh), intent(in) :: mesh
    character(len=*), intent(in) :: name
    type(parameter_list), intent(inout) :: params
    call this%problem%init(mesh, name, params)
    this%faces = this%problem%faces
    allocate(this%heat_flux_mask(size(this%faces)), source=.true.)
  end subroutine

  integer function coupling_size(this)
    class(view_factor_enclosure_coupling), intent(in) :: this
    coupling_size = this%problem%size()
  end function

  subroutine solve_radiosity(this, time, Tface, qrad, stat, numitr, error)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:)
    real(r8), intent(inout) :: qrad(:)
    integer, intent(out) :: stat
    integer, intent(out) :: numitr
    real(r8), intent(out) :: error
    call this%problem%solve_radiosity(time, Tface(this%faces), qrad, stat, numitr, error)
  end subroutine

  subroutine compute_radiosity_residual(this, time, Tface, qrad, r)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:)
    real(r8), intent(in) :: qrad(:)
    real(r8), intent(out) :: r(:)
    call this%problem%residual(time, qrad, Tface(this%faces), r)
  end subroutine

  subroutine add_heat_flux_residual(this, time, Tface, qrad, area, Fface)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:), qrad(:), area(:)
    real(r8), intent(inout) :: Fface(:)
    real(r8), allocatable :: flux(:)
    allocate(flux(size(this%faces)))
    call this%problem%heat_flux(time, qrad, Tface(this%faces), flux)
    call this%add_heat_flux_to_residual(area, flux, Fface)
  end subroutine

  subroutine add_heat_precon_deriv(this, time, Tface, area, matrix)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:), area(:)
    type(mfd_diff_matrix), intent(inout) :: matrix
    real(r8), allocatable :: values(:)
    allocate(values(size(this%faces)))
    call this%problem%rhs_deriv(time, Tface(this%faces), values)
    where (.not.this%heat_flux_mask) values = 0
    call matrix%incr_face_diag(this%faces, area(this%faces) * values)
  end subroutine

  subroutine apply_rad_precon(this, time, z)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(inout) :: z(:)
    call this%problem%precon(time, z)
  end subroutine

  subroutine apply_rad_precon_matvec1(this, time, z)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(inout) :: z(:)
    call this%problem%precon_matvec1(time, z)
  end subroutine

  subroutine compute_rhs(this, time, Tface, rhs)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:)
    real(r8), intent(out) :: rhs(:)
    call this%problem%rhs(time, Tface(this%faces), rhs)
  end subroutine

  subroutine add_rhs_deriv_times_face_vector(this, time, Tface, face_vector, qrad_rhs)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:), face_vector(:)
    real(r8), intent(inout) :: qrad_rhs(:)
    real(r8), allocatable :: values(:)
    allocate(values(size(this%faces)))
    call this%problem%rhs_deriv(time, Tface(this%faces), values)
    qrad_rhs = qrad_rhs + values * face_vector(this%faces)
  end subroutine

  subroutine add_heat_flux_to_residual(this, area, flux, Fface)
    class(view_factor_enclosure_coupling), intent(in) :: this
    real(r8), intent(in) :: area(:), flux(:)
    real(r8), intent(inout) :: Fface(:)
    integer :: j
    do j = 1, size(this%faces)
      if (this%heat_flux_mask(j)) &
          Fface(this%faces(j)) = Fface(this%faces(j)) + area(this%faces(j)) * flux(j)
    end do
  end subroutine

  subroutine suppress_heat_flux_on_faces(this, faces)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    integer, intent(in) :: faces(:)
    integer :: j
    do j = 1, size(this%faces)
      if (findloc(faces, this%faces(j), dim=1) > 0) this%heat_flux_mask(j) = .false.
    end do
  end subroutine

  subroutine set_initial_time(this, time)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    call this%problem%set_initial_time(time)
  end subroutine

  subroutine update_moving_vf(this)
    class(view_factor_enclosure_coupling), intent(inout) :: this
    call this%problem%update_moving_vf
  end subroutine

  subroutine add_moving_vf_events(this, eventq, rank)
    use sim_event_queue_type, only: sim_event_queue
    class(view_factor_enclosure_coupling), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    integer, intent(in), optional :: rank
    call this%problem%add_moving_vf_events(eventq, rank)
  end subroutine

  logical function is_active(this)
    class(thermal_view_factor_coupling), intent(in) :: this
    if (allocated(this%encl)) then
      is_active = (size(this%encl) > 0)
    else
      is_active = .false.
    end if
  end function

  integer function total_size(this)
    class(thermal_view_factor_coupling), intent(in) :: this
    ASSERT(allocated(this%encl))
    ASSERT(allocated(this%offset))
    total_size = this%offset(size(this%offset)) - 1
  end function

  integer function num_enclosures(this)
    class(thermal_view_factor_coupling), intent(in) :: this
    ASSERT(allocated(this%encl))
    num_enclosures = size(this%encl)
  end function

  subroutine enclosure_range(this, n, first, last)
    class(thermal_view_factor_coupling), intent(in) :: this
    integer, intent(in) :: n
    integer, intent(out) :: first, last
    ASSERT(allocated(this%encl))
    ASSERT(allocated(this%offset))
    first = this%offset(n)
    last = this%offset(n+1) - 1
  end subroutine


  subroutine coupling_solve_radiosity(this, time, Tface, qrad, stat, numitr, error)

    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:)
    real(r8), intent(inout) :: qrad(:)
    integer, intent(out) :: stat
    integer, intent(out) :: numitr
    real(r8), intent(out) :: error

    integer :: n, first, last, stat1, numitr1
    real(r8) :: error1

    ASSERT(allocated(this%encl))
    stat = 0
    numitr = 0
    error = 0.0_r8
    do n = 1, this%num_enclosures()
      call this%enclosure_range(n, first, last)
      call this%encl(n)%solve_radiosity(time, Tface, qrad(first:last), stat1, numitr1, error1)
      if (stat1 /= 0) stat = stat1
      numitr = numitr + numitr1
      error = max(error, error1)
    end do

  end subroutine coupling_solve_radiosity


  subroutine coupling_compute_radiosity_residual(this, time, Tface, qrad, r)

    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:)
    real(r8), intent(in) :: qrad(:)
    real(r8), intent(out) :: r(:)

    integer :: n, first, last

    ASSERT(allocated(this%encl))
    do n = 1, this%num_enclosures()
      call this%enclosure_range(n, first, last)
      call this%encl(n)%compute_radiosity_residual(time, Tface, qrad(first:last), r(first:last))
    end do

  end subroutine coupling_compute_radiosity_residual


  subroutine coupling_relative_residual_norm(this, time, Tface, qrad, error)

    use parallel_communication, only: global_sum

    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:)
    real(r8), intent(in) :: qrad(:)
    real(r8), intent(out) :: error

    real(r8) :: res2, rhs2
    real(r8), allocatable :: res(:), rhs(:)

    ASSERT(allocated(this%encl))
    allocate(res(this%size()), rhs(this%size()))
    call this%compute_radiosity_residual(time, Tface, qrad, res)
    call this%compute_rhs(time, Tface, rhs)

    if (size(res) > 0) then
      res2 = norm2(res)**2
      rhs2 = norm2(rhs)**2
    else
      res2 = 0.0_r8
      rhs2 = 0.0_r8
    end if
    res2 = global_sum(res2)
    rhs2 = global_sum(rhs2)

    if (rhs2 > 0.0_r8) then
      error = sqrt(res2/rhs2)
    else
      error = sqrt(res2)
    end if

  end subroutine coupling_relative_residual_norm


  subroutine coupling_add_heat_flux_residual(this, time, Tface, qrad, area, Fface)
    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:), qrad(:), area(:)
    real(r8), intent(inout) :: Fface(:)
    integer :: n, first, last
    ASSERT(allocated(this%encl))
    do n = 1, this%num_enclosures()
      call this%enclosure_range(n, first, last)
      call this%encl(n)%add_heat_flux_residual(time, Tface, qrad(first:last), area, Fface)
    end do
  end subroutine

  subroutine coupling_add_heat_precon_deriv(this, time, Tface, area, matrix)
    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:), area(:)
    type(mfd_diff_matrix), intent(inout) :: matrix
    integer :: n
    ASSERT(allocated(this%encl))
    do n = 1, this%num_enclosures()
      call this%encl(n)%add_heat_precon_deriv(time, Tface, area, matrix)
    end do
  end subroutine

  subroutine coupling_apply_rad_precon(this, time, z)
    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(inout) :: z(:)
    integer :: n, first, last
    ASSERT(allocated(this%encl))
    do n = 1, this%num_enclosures()
      call this%enclosure_range(n, first, last)
      call this%encl(n)%apply_rad_precon(time, z(first:last))
    end do
  end subroutine

  subroutine coupling_apply_rad_precon_matvec1(this, time, z)
    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(inout) :: z(:)
    integer :: n, first, last
    ASSERT(allocated(this%encl))
    do n = 1, this%num_enclosures()
      call this%enclosure_range(n, first, last)
      call this%encl(n)%apply_rad_precon_matvec1(time, z(first:last))
    end do
  end subroutine

  subroutine coupling_compute_rhs(this, time, Tface, rhs)
    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:)
    real(r8), intent(out) :: rhs(:)
    integer :: n, first, last
    ASSERT(allocated(this%encl))
    do n = 1, this%num_enclosures()
      call this%enclosure_range(n, first, last)
      call this%encl(n)%compute_rhs(time, Tface, rhs(first:last))
    end do
  end subroutine

  subroutine coupling_add_rhs_deriv_times_face_vector(this, time, Tface, face_vector, qrad_rhs)
    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: Tface(:), face_vector(:)
    real(r8), intent(inout) :: qrad_rhs(:)
    integer :: n, first, last
    ASSERT(allocated(this%encl))
    do n = 1, this%num_enclosures()
      call this%enclosure_range(n, first, last)
      call this%encl(n)%add_rhs_deriv_times_face_vector(time, Tface, face_vector, qrad_rhs(first:last))
    end do
  end subroutine

  subroutine coupling_add_heat_flux_to_residual(this, area, flux, Fface)
    class(thermal_view_factor_coupling), intent(in) :: this
    real(r8), intent(in) :: area(:), flux(:)
    real(r8), intent(inout) :: Fface(:)
    integer :: n, first, last
    ASSERT(allocated(this%encl))
    do n = 1, this%num_enclosures()
      call this%enclosure_range(n, first, last)
      call this%encl(n)%add_heat_flux_to_residual(area, flux(first:last), Fface)
    end do
  end subroutine

  subroutine coupling_set_initial_time(this, time)
    class(thermal_view_factor_coupling), intent(inout) :: this
    real(r8), intent(in) :: time
    integer :: n
    if (.not.this%is_active()) return
    do n = 1, this%num_enclosures()
      call this%encl(n)%set_initial_time(time)
    end do
  end subroutine

  subroutine coupling_update_moving_vf(this)
    class(thermal_view_factor_coupling), intent(inout) :: this
    integer :: n
    if (.not.this%is_active()) return
    do n = 1, this%num_enclosures()
      call this%encl(n)%update_moving_vf
    end do
  end subroutine

  subroutine coupling_add_moving_vf_events(this, eventq, rank)
    use sim_event_queue_type, only: sim_event_queue
    class(thermal_view_factor_coupling), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    integer, intent(in), optional :: rank
    integer :: n
    if (.not.this%is_active()) return
    do n = 1, this%num_enclosures()
      call this%encl(n)%add_moving_vf_events(eventq, rank)
    end do
  end subroutine

  subroutine coupling_init(this, mesh, params, stat, errmsg)

    use bitfield_type, only: btest
    use parallel_communication, only: global_any, global_all

    class(thermal_view_factor_coupling), intent(out) :: this
    type(unstr_mesh), intent(in) :: mesh
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n, j, precon_iter
    real(r8) :: solve_tol
    character(:), allocatable :: precon_method
    logical, allocatable :: mask(:)
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist

    stat = 0

    !! Initialize the enclosure radiation (view factor) problems, if any.
    piter = parameter_list_iterator(params, sublists_only=.true.)
    n = piter%count()
    if (n == 0) return

    call params%get('solve-tol', solve_tol, default=1.0e-3_r8)
    call params%get('precon-method', precon_method, default='JACOBI')
    call params%get('precon-iter', precon_iter, default=1)

    allocate(this%encl(n))
    allocate(this%offset(n+1))
    this%offset(1) = 1
    do j = 1, n
      plist => piter%sublist()
      call plist%set('error-tol', solve_tol)
      call plist%set('precon-method', precon_method)
      call plist%set('precon-iter', precon_iter)
      call this%encl(j)%init(mesh, piter%name(), plist)
      this%offset(j+1) = this%offset(j) + this%encl(j)%size()
      !! Verify that these enclosure faces are boundary faces.
      if (.not.global_all(btest(mesh%face_set_mask(this%encl(j)%faces),0))) then
        stat = -1
        errmsg = 'some enclosure faces are not boundary faces'
        exit
      else if (n > 1) then  ! multiple enclosures
        if (j == 1) then
          allocate(mask(mesh%nface))
          mask = .false.
        !! Verify that these faces don't overlap with other enclosure faces.
        else if (global_any(mask(this%encl(j)%faces))) then
          stat = -1
          errmsg = 'some enclosure faces belong to other enclosures'
          exit
        !! Tag the faces belonging to this enclosure.
        else
          mask(this%encl(j)%faces) = .true.
        end if
      end if
      call piter%next
    end do
    if (allocated(mask)) deallocate(mask)

  end subroutine coupling_init


  subroutine validate_bc(this, mesh, thermal, stat, errmsg)

    use bitfield_type
    use parallel_communication, only: global_any, global_count
    use string_utilities, only: i_to_c

    class(thermal_view_factor_coupling), intent(inout) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(thermal_component), intent(inout) :: thermal
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n
    logical, allocatable :: mask(:), rmask(:)
    integer, allocatable :: setids(:)
    character(160) :: string1, string2
    type(bitfield) :: bitmask

    allocate(mask(mesh%nface), rmask(mesh%nface))

    mask = .false.
    rmask = .false.

    !! Mark all faces covered by standard thermal boundary/interface conditions.
    if (allocated(thermal%ic_htc)) then
      mask(thermal%ic_htc%index(1,:)) = .true.
      mask(thermal%ic_htc%index(2,:)) = .true.
    end if
    if (allocated(thermal%ic_rad)) then
      mask(thermal%ic_rad%index(1,:)) = .true.
      mask(thermal%ic_rad%index(2,:)) = .true.
    end if
    if (allocated(thermal%ic_htc) .or. allocated(thermal%ic_rad)) then
      call mesh%face_imap%gather_offp(mask)
    end if
    if (allocated(thermal%bc_flux)) mask(thermal%bc_flux%index) = .true.
    if (allocated(thermal%bc_vflux)) then
      do j = 1, size(thermal%bc_vflux%index)
        mask(thermal%bc_vflux%index(j)) = .true.
      end do
    end if
    if (allocated(thermal%bc_htc)) mask(thermal%bc_htc%index) = .true.
    if (allocated(thermal%bc_rad)) mask(thermal%bc_rad%index) = .true.
    if (allocated(thermal%bc_dir)) mask(thermal%bc_dir%index) = .true.

    !! Mark all faces covered by enclosure radiation.
    if (this%is_active()) then
      do n = 1, this%num_enclosures()
        rmask(this%encl(n)%faces) = .true.
      end do
      call mesh%face_imap%gather_offp(rmask)
    end if

    !! Simple radiation BCs and enclosure radiation are distinct models and
    !! cannot both contribute on the same face.
    if (allocated(thermal%bc_rad)) then
      if (global_any(rmask(thermal%bc_rad%index))) then
        stat = -1
        errmsg = 'temperature radiation boundary condition overlaps with enclosure radiation'
        return
      end if
    end if

    !! Temperature Dirichlet conditions suppress the heat-equation flux
    !! contribution from overlapping enclosure radiation faces.
    if (this%is_active() .and. allocated(thermal%bc_dir)) then
      do n = 1, this%num_enclosures()
        call this%encl(n)%suppress_heat_flux_on_faces(thermal%bc_dir%index)
      end do
    end if

    !! Finally verify that a condition has been applied to every boundary face.
    mask = (mask .or. rmask) .neqv. btest(mesh%face_set_mask,0)
    if (global_any(mask)) then
      call mesh%get_face_set_ids(pack([(j,j=1,mesh%nface)], mask), setids)
      if (size(setids) == 0) then
        string1 = '(none)'
      else
        write(string1,'(i0,*(:,", ",i0))') setids
      end if
      call mesh%get_link_set_ids(mask, setids)
      if (size(setids) == 0) then
        string2 = '(none)'
      else
        write(string2,'(i0,*(:,", ",i0))') setids
      end if
      errmsg = 'incomplete temperature boundary/interface specification;' // &
          ' remaining boundary faces belong to face sets ' // trim(string1) // &
          '; and interface link sets ' // trim(string2)
      bitmask = ibset(ZERO_BITFIELD, 0)
      mask = mask .and. (mesh%face_set_mask == bitmask)
      mask(mesh%lface(1,:)) = .false.
      mask(mesh%lface(2,:)) = .false.
      n = global_count(mask(:mesh%nface_onP))
      if (n > 0) errmsg = errmsg // '; ' // i_to_c(n) // ' faces belong to neither'
      stat = -1
      return
    end if

    stat = 0

  end subroutine validate_bc

end module thermal_view_factor_coupling_type
