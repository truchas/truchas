module flow_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use unstr_mesh_type
  use index_partitioning
  implicit none
  private

  public :: flow

  type :: flow
    type(unstr_mesh), pointer :: mesh ! unowned reference
    real(r8), allocatable :: vel_cc(:,:) ! cell-centered velocity (dims, ncells)

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
  end subroutine read_params


  subroutine init(this, m)
    class(flow), intent(inout) :: this
    type(unstr_mesh), pointer, intent(in) :: m

    this%mesh => m

    print *, 'allocation required in flow%init'

    print *, 'initialization required in flow%init'
  end subroutine init


  subroutine step(this, t, dt)
    class(flow), intent(inout) :: this
    real(r8), intent(in) :: t, dt

    !=================================================================
    ! outline of routine from truchas below
    !=================================================================
    fluidRho = 0.0_r8

    ! =================================================================
    ! ?? move fluid property calculations to driver...
    ! Do we need to track things like isPureImmobile, Solid_face...
    ! Or can they be hidden from the flow solver and handled in material property
    ! calculations only...
    ! =================================================================
    ! Evaluate cell properties excluding immobile materials, and
    ! check that there are at least some flow equations to solve
    call FLUID_PROPERTIES (abort, t)

    if (.not. abort) then

       fluid_to_move = .true.
       ! Predictor Step
        call PREDICTOR()

        ! Projection Step
        call PROJECTION()

        if(cycle_number == 0) then
           ! Special operations required during the prepass
           prelim_projection_iterations = mac_projection_iterations
           prelim_viscous_iterations = viscous_iterations
           do n = 1,ndim
              Zone%Vc(n) = Zone%Vc_Old(n)
           end do
        endif

    else

       ! Everything solid; set velocities equal to zero, and check again in the
       ! next timestep.
       fluid_to_move = .false.
       Fluxing_Velocity = 0
       do n = 1,ndim
          Zone%Vc(n) = 0
          Zone%Vc_Old(n) = 0
       end do

       if(cycle_number == 0) then
          prelim_projection_iterations = 0
          prelim_viscous_iterations = 0
       endif

    endif
  end subroutine step

end module flow_type
