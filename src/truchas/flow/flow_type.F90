#include "f90_assert.fpp"

module flow_type

  use kinds, only: r8
  use flow_domain_types
  use constants_module
  use truchas_logging_services
  use truchas_timers
  use flow_projection_type
  use flow_prediction_type
  use flow_props_type
  use flow_bc_type
  use unstr_mesh_type
  use index_partitioning
  use parallel_communication
  implicit none
  private

  public :: flow

  type :: flow
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(flow_props), pointer :: props => null() ! unowned reference
    real(r8), allocatable :: vel_cc(:,:) ! cell-centered velocity (dims, ncell)
    real(r8), allocatable :: vel_cc_n(:,:) ! cell-centered velocity (dims, ncell)
    real(r8), allocatable :: vel_fn(:) ! outward oriented face-normal velocity
    real(r8), allocatable :: P_cc(:) ! cell-centered pressure
    real(r8), allocatable :: grad_p_rho_cc_n(:,:) ! dynamic pressure gradient over rho
    real(r8) :: body_force(3)
    type(flow_projection) :: proj
    type(flow_prediction) :: pred
    type(flow_bc), pointer :: bc
    logical :: inviscid, stokes, prescribed
    real(r8) :: viscous_number
    real(r8) :: courant_number
  contains
    procedure :: init
    procedure :: step
    procedure :: accept
    procedure :: correct_non_regular_cells
    procedure :: timestep
    ! views into data
    procedure :: vel_cc_view
    procedure :: vel_fn_view
    procedure :: P_cc_view
  end type flow

contains

  function vel_cc_view(this) result(p)
    class(flow), target, intent(in) :: this
    real(r8), pointer :: p(:,:)
    p => this%vel_cc
  end function vel_cc_view

  function P_cc_view(this) result(p)
    class(flow), target, intent(in) :: this
    real(r8), pointer :: p(:)
    p => this%P_cc
  end function P_cc_view

  function vel_fn_view(this) result(p)
    class(flow), target, intent(in) :: this
    real(r8), pointer :: p(:)
    p => this%vel_fn
  end function vel_fn_view

  subroutine timestep(this, dtc, dtv)
    use parallel_communication

    class(flow), intent(in) :: this
    real(r8), intent(out) :: dtc, dtv
    real(r8) :: v
    integer :: i, j

    dtc = huge(1.0_r8)
    dtv = huge(1.0_r8)

    ! try something different (simpler) compared to truchas
    associate (vof => this%props%vof, rho => this%props%rho_cc, mu => this%props%mu_cc, m => this%mesh)
      do j = 1, m%nface_onP
        dtc = min(dtc, m%face_normal_dist(j)/(epsilon(1.0_r8)+abs(this%vel_fn(j))))
      end do

      if (.not.this%inviscid) then
        do j = 1, m%ncell_onP
          if (this%props%cell_t(j) == regular_t .and. mu(j) > 0.0_r8) then
            v = (m%volume(j)*vof(j))**(1.0_r8/3.0_r8)
            dtv = min(dtv, v**2*rho(j)/mu(j))
          end if
        end do
      end if
    end associate

    dtc = this%courant_number * global_minval(dtc)
    dtv = this%viscous_number * global_minval(dtv)

  end subroutine timestep

  subroutine init(this, mesh, props, prescribed, vel_cc, P_cc, vel_fn, params)

    use parameter_list_type

    class(flow), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(flow_props), intent(in), target :: props
    logical, optional, intent(in) :: prescribed
    real(r8), optional, intent(in) :: vel_cc(:,:), P_cc(:), vel_fn(:)
    type(parameter_list), intent(inout) :: params

    integer :: i
    real(r8) :: xc, yc, sigma, x, y, r2, theta, w, vel(3)
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: array(:)

    this%mesh => mesh
    this%props => props

    plist => params%sublist('options')
    call plist%get('inviscid', this%inviscid, default=.false.)
    call plist%get('stokes', this%stokes, default=.false.)
    call plist%get('viscous number', this%viscous_number, default=0.1_r8)
    call plist%get('courant number', this%courant_number, default=0.5_r8)
    call plist%get('body force', array, default=[0.0_r8, 0.0_r8, 0.0_r8])
    this%body_force = array

    associate (nc => mesh%ncell, nf => mesh%nface)
      allocate(this%vel_cc(3, nc))
      allocate(this%vel_cc_n(3, nc))
      allocate(this%P_cc(nc))
      allocate(this%vel_fn(nf))
      allocate(this%grad_p_rho_cc_n(3,nc))
    end associate

    if (present(vel_cc)) then
      do i = 1, mesh%ncell_onP
        this%vel_cc(:,i) = vel_cc(:,i)
      end do
      call gather_boundary(mesh%cell_ip, this%vel_cc)
      ! set the face velocities
      if (present(vel_fn)) then
        ASSERT(size(vel_fn) == this%mesh%nface_onP)
        this%vel_fn(:this%mesh%nface_onP) = vel_fn
      else
        ! WARN: This will only give reasonable results when the initial velocity
        !       is uniform. A real solution should get velocities from an input
        !       function and ensure the face velocities are discretely solenoidal.
        do i = 1, mesh%nface_onP
          vel = this%vel_cc(:,mesh%fcell(1,i))
          if (mesh%fcell(2,i) /= 0) then
            vel = vel + this%vel_cc(:,mesh%fcell(2,i))
            vel = vel / 2
          end if
          this%vel_fn(i) = dot_product(mesh%normal(:,i), vel)
        end do
      end if
      call gather_boundary(mesh%face_ip, this%vel_fn)
    else
      this%vel_cc = 0.0_r8
      this%vel_fn = 0.0_r8
    end if
    this%vel_cc_n = this%vel_cc
    if (present(P_cc)) then
      INSIST(size(p_cc) == mesh%ncell_onP)
      this%P_cc(:mesh%ncell_onP) = P_cc
      call gather_boundary(mesh%cell_ip, this%p_cc)
    else
      this%P_cc = 0.0_r8
    end if

    this%grad_p_rho_cc_n = 0.0_r8

    ! Disable BC initialization and the prediction and
    ! projection solvers if we have prescribed flow.
    if (present(prescribed)) then
      this%prescribed = prescribed
    else
      this%prescribed = .false.
    end if
    if (this%prescribed) return

    allocate(this%bc)
    call this%bc%init(mesh, params)
    call this%pred%init(mesh, this%bc, this%inviscid, this%stokes, params)
    call this%proj%init(mesh, this%bc, params)

    !call this%proj%grad_p_rho

  end subroutine init

  subroutine correct_non_regular_cells(this)
    class(flow), intent(inout) :: this
    !-
    integer :: i

    associate (cell_t => this%props%cell_t, face_t => this%props%face_t)
      do i = 1, this%mesh%ncell
        if (cell_t(i) /= regular_t) then
          this%vel_cc(:,i) = 0.0_r8
          this%P_cc(i) = 0.0_r8
        end if
      end do

      do i = 1, this%mesh%nface_onP
        if (face_t(i) > regular_t) then
          this%vel_fn(i) = 0.0_r8
        end if
      end do
    end associate
  end subroutine correct_non_regular_cells

  ! note 1: On the initial pass, calculate initial pressure
  !         and pressure gradients, but revert cell-centered
  !         and face-centered velocities to initial values.
  !         Some issues can arise here with non-uniform, but
  !         analytically solenoidal velocity fields, which
  !         are difficult to set discretely on faces of an
  !         unstructured mesh. One pass of the flow solver
  !         calculate discretely divergence-free face
  !         velocities, but one timestep offset from the
  !         cell-centered velocities, which introduces an
  !         inconsistency in unsteady flow because the
  !         predictor uses the face velocity as the donor-cell
  !         velocity along boundaries. Since this term uses
  !         a mix of cell- and face-centered velocities,
  !         the initial mismatch introduces a momentum sync
  !         in the first timestep. By reverting the
  !         face-velocities to initial values, we avoid this
  !         issue and produce the correct result for initially
  !         uniform velocity fields. Further consideration is
  !         necessary for non-uniform initial velocity fields.
  subroutine step(this, t, dt, flux_volumes, initial)
    class(flow), intent(inout) :: this
    real(r8), intent(in) :: t, dt, flux_volumes(:,:)
    logical, optional, intent(in) ::  initial
    !-
    real(r8) :: p_max, p_min
    real(r8), allocatable :: vel_fn_n(:)
    integer :: j

    ! a useless copy when the previous solution hase been accepted but
    ! probably not a performance bottleneck
    this%vel_cc = this%vel_cc_n

    if (present(initial)) then
      if (initial) then
        vel_fn_n = this%vel_fn
      end if
    end if

    call this%bc%compute(t, dt, initial=initial)

    call start_timer("prediction")
    call this%pred%setup(dt, this%props, this%vel_cc, initial=initial)
    call this%pred%solve(dt, this%props, this%grad_p_rho_cc_n, flux_volumes, this%vel_fn, this%vel_cc, initial=initial)
    call stop_timer("prediction")

    call start_timer("projection")
    call this%proj%setup(dt, this%props, this%grad_p_rho_cc_n, this%body_force, this%vel_cc, this%P_cc, this%vel_fn, initial=initial)
    call this%proj%solve(dt, this%props, this%grad_p_rho_cc_n, this%vel_cc, this%P_cc, this%vel_fn, initial=initial)
    call stop_timer("projection")

    call this%correct_non_regular_cells()
    !p_min = global_minval(this%p_cc)
    !p_max = global_maxval(this%p_cc)

    if (present(initial)) then
      if (initial) then
        this%vel_cc = this%vel_cc_n
        this%vel_fn = vel_fn_n
      end if
    end if
  end subroutine step

  subroutine accept(this)
    class(flow), intent(inout) :: this

    this%vel_cc_n = this%vel_cc

    ! disable the prediction and projection solvers if we have prescribed flow
    if (this%prescribed) return

    call this%pred%accept()
    call this%proj%accept(this%grad_p_rho_cc_n)

  end subroutine accept

end module flow_type
