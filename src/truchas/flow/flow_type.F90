module flow_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use flow_mesh_type
  use flow_projection_type
  use flow_prediction_type
  use flow_props_type
  use flow_bc_type

  use unstr_mesh_type
  use index_partitioning
  implicit none
  private

  public :: flow

  type :: flow
    type(flow_mesh), pointer :: mesh ! unowned reference
    real(r8), allocatable :: vel_cc(:,:) ! cell-centered velocity (dims, ncell)
    real(r8), allocatable :: vel_fc(:,:) ! face-centered velocity (dims, nface)
    real(r8), allocatable :: vel_fn(:) ! outward oriented face-normal velocity
    real(r8), allocatable :: P_cc(:) ! cell-centered pressure
    real(r8), allocatable :: grad_p_rho_cc(:,:) ! dynamic pressure gradient over rho
    real(r8) :: body_force(3)
    type(flow_projection) :: proj
    type(flow_prediction) :: pred
    type(flow_bc), pointer :: bc
    logical :: active
  contains
    procedure :: read_params
    procedure :: init
    procedure :: step
    procedure :: accept
    procedure :: zero_out_solid_velocities
  end type flow

contains

  subroutine read_params(this, p)
    use parameter_list_type
    class(flow), intent(inout) :: this
    type(parameter_list), pointer, intent(in) :: p
    !
    real(r8), allocatable :: body_force(:)
    integer :: stat
    type(parameter_list), pointer :: sub

    call p%get('active', this%active, .true.)

    call p%get('body-force', body_force, stat=stat)
    if (stat == 0) then
      this%body_force = body_force
    else
      this%body_force = 0.0_r8
    end if

    sub => p%sublist("projection")
    call this%proj%read_params(sub)

    sub => p%sublist("prediction")
    call this%pred%read_params(sub)

  end subroutine read_params


  subroutine init(this, m, bc, vel_cc, P_cc)
    class(flow), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: m
    type(flow_bc), pointer, intent(inout) :: bc
    real(r8), optional, intent(in) :: vel_cc, P_cc
    !-
    this%mesh => m
    this%bc => bc

    associate (nc => m%mesh%ncell, nf => m%mesh%nface)
      allocate(this%grad_p_rho_cc(3,nc))
      allocate(this%vel_cc(3, nc))
      allocate(this%vel_fc(3, nf))
      allocate(this%P_cc(nc))
      allocate(this%vel_fn(nf))
    end associate

    if (present(vel_cc)) then
      this%vel_cc = vel_cc
    else
      this%vel_cc = 0.0_r8
    end if

    if (present(P_cc)) then
      this%P_cc = P_cc
    else
      this%P_cc = 0.0_r8
    end if

    call this%pred%init(m, this%bc)
    call this%proj%init(m, this%bc)
  end subroutine init

  subroutine zero_out_solid_velocities(this, props)
    class(flow), intent(inout) :: this
    type(flow_props), intent(in) :: props

    print *, 'zero out all velocities on solid cells and faces'
  end subroutine zero_out_solid_velocities

  subroutine step(this, t, dt, props)
    class(flow), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(flow_props), intent(inout) :: props

    call this%bc%compute(t)

    call this%pred%setup(dt, props, this%vel_cc)
    call this%pred%solve(...)

    call this%proj%setup(dt, props, this%body_force, this%vel_cc, this%P_cc, this%vel_fn)
    call this%proj%solve(dt, props, this%vel_cc, this%P_cc, this%vel_fn)

  end subroutine step

  subroutine accept(this)
    class(flow), intent(inout) :: this

    call props%accept()
    call this%pred%accept()
    call this%proj%accept()

  end subroutine accept

end module flow_type
