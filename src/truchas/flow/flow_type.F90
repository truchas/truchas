#define ASDF 0
#define Q 717
module flow_type

  use kinds, only: r8
  use flow_domain_types
  use constants_module
  use truchas_logging_services
  use truchas_timers
  use flow_mesh_type
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
    type(flow_mesh), pointer :: mesh => null() ! unowned reference
    real(r8), allocatable :: vel_cc(:,:) ! cell-centered velocity (dims, ncell)
    real(r8), allocatable :: vel_cc_n(:,:) ! cell-centered velocity (dims, ncell)
    real(r8), allocatable :: vel_fn(:) ! outward oriented face-normal velocity
    real(r8), allocatable :: P_cc(:) ! cell-centered pressure
    real(r8), allocatable :: grad_p_rho_cc_n(:,:) ! dynamic pressure gradient over rho
    real(r8) :: body_force(3)
    type(flow_projection) :: proj
    type(flow_prediction) :: pred
    type(flow_bc), pointer :: bc
    logical :: inviscid, stokes
    real(r8) :: viscous_number
    real(r8) :: courant_number
  contains
    procedure :: read_params
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

  subroutine timestep(this, props, dtc, dtv)
    use parallel_communication

    class(flow), intent(in) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(out) :: dtc, dtv
    real(r8) :: v
    integer :: i, j

    dtc = huge(1.0_r8)
    dtv = huge(1.0_r8)

    ! try something different (simpler) compared to truchas
    associate( m => this%mesh%mesh, vof => props%vof, &
        rho => props%rho_cc, mu => props%mu_cc)

      do j = 1, m%ncell_onP
        if (vof(j) > 0.0_r8) then
          v = (m%volume(j)*vof(j))**(1.0_r8/3.0_r8)
          dtc = min(dtc, v/(epsilon(1.0_r8)+norm2(this%vel_cc(:,j))))
          if (.not.this%inviscid .and. mu(j) > 0.0_r8) &
              dtv = min(dtv, v**2*rho(j)/mu(j))
        end if
      end do
    end associate

    dtc = this%courant_number * global_minval(dtc)
    dtv = this%viscous_number * global_minval(dtv)

  end subroutine timestep

  subroutine read_params(this, p)
    use parameter_list_type
    class(flow), intent(inout) :: this
    type(parameter_list), pointer, intent(in) :: p
    type(parameter_list), pointer :: pp
    real(r8), allocatable :: body_force(:)
    integer :: stat
    !

    pp => p%sublist("options")

    call pp%get('inviscid', this%inviscid, .false.)
    call pp%get('stokes', this%stokes, .false.)
    call pp%get('viscous number', this%viscous_number, 0.1_r8)
    call pp%get('courant number', this%courant_number, 0.5_r8)
    call pp%get('body force', body_force, stat=stat)
    if (stat == 0) then
      this%body_force = body_force
    else
      this%body_force = 0.0_r8
    end if

    allocate(this%bc)
    call this%bc%read_params(p)
    call this%proj%read_params(p)
    call this%pred%read_params(p)

  end subroutine read_params


  subroutine init(this, m, vel_cc, P_cc)
    class(flow), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: m
    real(r8), optional, intent(in) :: vel_cc(:,:), P_cc
    !-
    integer :: i
    real(r8) :: xc, yc, sigma, x, y, r2, theta, w, vel(3)

    this%mesh => m

    call this%bc%init(m)

    associate (nc => m%mesh%ncell, nf => m%mesh%nface)
      allocate(this%vel_cc(3, nc))
      allocate(this%vel_cc_n(3, nc))
      allocate(this%P_cc(nc))
      allocate(this%vel_fn(nf))
      allocate(this%grad_p_rho_cc_n(3,nc))
    end associate


    if (present(vel_cc)) then
      ! vortex centered at x=0.5, y=0.5
      xc = 0.5_r8
      yc = 0.5_r8
      sigma = 0.1_r8
#if ASDF
      print* , ">>>> HIGHJACKING FLOW INTIAL CONDITIONS <<<<"

      do i = 1, m%mesh%ncell
        x = m%cell_centroid(1,i) - xc
        y = m%cell_centroid(2,i) - yc
        !print*, y, x
        theta = atan2(y, x)
        r2 = (sqrt(x**2+y**2) - 0.25_r8)**2
        ! gaussian
        w = exp(-r2/(2.0_r8*sigma**2))/sqrt(2.0_r8*pi*sigma)
        this%vel_cc(:,i) = [-w*sin(theta), w*cos(theta), 0.0_r8]
      end do
#else
      do i = 1, m%mesh%ncell_onP
        this%vel_cc(:,i) = vel_cc(:,i)
      end do
      call gather_boundary(m%mesh%cell_ip, this%vel_cc)
#endif
      ! set the face velocities
      ! WARN: This will only give reasonable results when the initial velocity
      !       is uniform. A real solution should get velocities from an input
      !       function and ensure the face velocities are discretely solenoidal.
      do i = 1, m%mesh%nface_onP
        vel = this%vel_cc(:,m%fcell(2,i))
        if (m%fcell(1,i) /= 0) then
          vel = vel + this%vel_cc(:,m%fcell(1,i))
          vel = vel / 2
        end if
        this%vel_fn(i) = dot_product(m%mesh%normal(:,i), vel)
      end do
      call gather_boundary(m%mesh%face_ip, this%vel_fn)
    else
      this%vel_cc = 0.0_r8
      this%vel_fn = 0.0_r8
    end if
    this%vel_cc_n = this%vel_cc
    if (present(P_cc)) then
      this%P_cc = P_cc
    else
      this%P_cc = 0.0_r8
    end if

    this%grad_p_rho_cc_n = 0.0_r8

    call this%pred%init(m, this%bc, this%inviscid, this%stokes)
    call this%proj%init(m, this%bc)

  end subroutine init

  subroutine correct_non_regular_cells(this, props)
    class(flow), intent(inout) :: this
    type(flow_props), intent(in) :: props
    !-
    integer :: i

    associate (m => this%mesh%mesh, cell_t => props%cell_t, face_t => props%face_t)
      do i = 1, m%ncell
        if (cell_t(i) /= regular_t) then
          this%vel_cc(:,i) = 0.0_r8
          this%P_cc(i) = 0.0_r8
        end if
      end do

      do i = 1, m%nface_onP
        if (face_t(i) /= regular_t) then
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
  subroutine step(this, t, dt, props, flux_volumes, initial)
    class(flow), intent(inout) :: this
    real(r8), intent(in) :: t, dt, flux_volumes(:,:)
    type(flow_props), intent(inout) :: props
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

#if ASDF
    write(*, '("Pre-Predictor u[",i4,"]: ", 3es20.12)') Q, this%vel_cc(:,Q)
    write(*, '("Pre-Predictor P[",i4,"]: ", es20.12)') Q, this%P_cc(Q)
#endif
    call start_timer("prediction")
    call this%pred%setup(dt, props, this%vel_cc, initial=initial)
    call this%pred%solve(dt, props, this%grad_p_rho_cc_n, flux_volumes, this%vel_fn, this%vel_cc, initial=initial)
    call stop_timer("prediction")
#if ASDF
    write(*, '("Post-Predictor u[",i4,"]: ", 3es20.12)') Q, this%vel_cc(:,Q)
    write(*, '("Post-Predictor P[",i4,"]: ", es20.12)') Q, this%P_cc(Q)
#endif

    call start_timer("projection")
    call this%proj%setup(dt, props, this%grad_p_rho_cc_n, this%body_force, this%vel_cc, this%P_cc, this%vel_fn, initial=initial)
    call this%proj%solve(dt, props, this%grad_p_rho_cc_n, this%vel_cc, this%P_cc, this%vel_fn, initial=initial)
    call stop_timer("projection")
#if ASDF
    write(*, '("Post-Projection u[",i4,"]: ", 3es20.12)') Q, this%vel_cc(:,Q)
    write(*, '("Post-Projection P[",i4,"]: ", es20.12)') Q, this%P_cc(Q)

    associate (m=>this%mesh%mesh)
      write(*, '("Post-Projection vel_fn[",i4,"]: ", 6es20.12)') Q, this%vel_fn(m%cface(m%xcface(Q):m%xcface(Q+1)-1))
    end associate
#endif

    call this%correct_non_regular_cells(props)
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

    call this%pred%accept()
    call this%proj%accept(this%grad_p_rho_cc_n)

  end subroutine accept

end module flow_type
