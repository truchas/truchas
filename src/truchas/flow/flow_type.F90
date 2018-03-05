module flow_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use flow_mesh_type
  use flow_projection_type
  use unstr_mesh_type
  use index_partitioning
  implicit none
  private

  public :: flow

  type :: flow
    type(flow_mesh), pointer :: mesh ! unowned reference
    real(r8), allocatable :: vel_cc(:,:) ! cell-centered velocity (dims, ncells)
    real(r8), allocatable :: vel_fn(:) ! outward oriented face-normal velocity
    real(r8), allocatable :: P_cc(:) ! cell-centered pressure
    type(flow_projection) :: fp
  contains
    procedure :: read_params
    procedure :: init
  end type flow

contains

  subroutine read_params(this, p)
    use parameter_list_type
    class(flow), intent(inout) :: this
    type(parameter_list), intent(inout) :: p

    print *, 'figure out flow parameters'
    call fp%read_params(p)
  end subroutine read_params


  subroutine init(this, m)
    class(flow), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: m
    !-
    type(unstr_mesh), pointer :: mesh

    this%mesh => m
    mesh => m%mesh

    allocate(this%vel_cc(3, mesh%ncell))
    allocate(this%P_cc(mesh%ncell))
    allocate(this%vel_fn(size(m%cface)))

    print *, 'allocation required in flow%init'

    print *, 'initialization required in flow%init'
    call fp%init(m)
  end subroutine init

  subroutine zero_out_solid_velocities(this, props)
    class(flow), intent(inout) :: this
    type(flow_props), intent(in) :: props

    print *, 'zero out all velocities on solid cells and faces'
  end subroutine zero_out_solid_velocities

  subroutine step(this, t, dt, props)
    class(flow), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(flow_props), intent(in) :: props

    ! predictor

    ! corrector
    call this%fp%setup(props, this%vel, dt)

  end subroutine step

end module flow_type
