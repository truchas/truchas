MODULE FLUID_FLOW_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures needed to solve the Navier-Stokes equations.
  !
  !   Public Interface:
  !
  !     call NAVIER_STOKES
  !
  !       Main driver for incrementing the Navier-Stokes equations
  !       by one time step.
  !
  ! Contains: NAVIER_STOKES
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  implicit none

  ! Private Module
  private

  ! Public Subroutines
  public :: NAVIER_STOKES, FLUID_FLOW_DRIVER


  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE FLUID_FLOW_DRIVER ()

    use fluid_data_module, only: fluid_flow, applyflow, fluid_to_move

    if (fluid_flow) then
      if (applyflow) then
        fluid_to_move = .true.
        call appliedflowfield()
      else
        call NAVIER_STOKES()
      end if
    end if

  END SUBROUTINE FLUID_FLOW_DRIVER

  SUBROUTINE NAVIER_STOKES ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Navier-Stokes (NS) driver: increment NS equations by one time step.
    !
    !======================================================================
    use constants_module,       only: zero
    use fluid_data_module,      only: fluid_flow, fluidRho, Solid_Face, &
                                      isPureImmobile, fluidDeltaRho,      &
                                      fluid_to_move, Fluxing_Velocity
    use kind_module,            only: int_kind, real_kind, log_kind
    use parameter_module,       only: ncells, nfc, ndim
    use predictor_module,       only: PREDICTOR
    use projection_data_module, only: Face_Density, mac_projection_iterations, &
                                      prelim_projection_iterations
    use projection_module,      only: PROJECTION
    use property_module,        only: FLUID_PROPERTIES
    use time_step_module,       only: cycle_number
    use viscous_data_module,    only: prelim_viscous_iterations, &
                                      viscous_iterations
    use zone_module,            only: Zone
    use truchas_logging_services

    implicit none

    ! Local Variables
    integer                                      :: status
    integer(int_kind)                            :: n
    logical(log_kind)                            :: abort

    ! Argument List

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (.not. fluid_flow) return

    ALLOCATE (Solid_Face(nfc,ncells), STAT = status)
    if (status /= 0) &
         call TLS_panic ('FLUID_FLOW: Solid_Face(nfc,ncells) allocation failed')
    ALLOCATE (isPureImmobile(ncells), STAT = status)
    if (status /= 0) &
         call TLS_panic ('FLUID_FLOW: isPureImmobile(ncells) allocation failed')
    ALLOCATE (fluidDeltaRho(ncells), STAT = status)
    if (status /= 0) &
         call TLS_panic ('FLUID_FLOW: fluidDeltaRho(ncells) allocation failed')
    ALLOCATE (Face_Density(nfc,ncells), STAT = status)
    if (status /= 0) &
       call TLS_panic ('FLUID_FLOW: Face_Density(nfc,ncells) allocation failed')

    fluidRho = 0.0_real_kind

    ! Evaluate cell properties excluding immobile materials, and
    ! check that there are at least some flow equations to solve
    call FLUID_PROPERTIES (abort)

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
       Fluxing_Velocity = zero
       do n = 1,ndim
          Zone%Vc(n) = zero
          Zone%Vc_Old(n) = zero
       end do

       if(cycle_number == 0) then
          prelim_projection_iterations = 0
          prelim_viscous_iterations = 0
       endif
       
    endif

    DEALLOCATE (Solid_Face)
    DEALLOCATE (isPureImmobile)
    DEALLOCATE (fluidDeltaRho)
    DEALLOCATE (Face_Density)

    return

  END SUBROUTINE NAVIER_STOKES

  SUBROUTINE appliedflowfield ()

    !=======================================================================
    ! Purpose(s):
    !
    !   Navier-Stokes (NS) driver: increment NS equations by one time step.
    !
    !======================================================================
    use fluid_data_module,      only: fluid_flow, Solid_Face, isPureImmobile, fluidDeltaRho
    use kind_module,            only: log_kind
    use parameter_module,       only: ncells, nfc
    use projection_data_module, only: Face_Density
    use property_module,        only: FLUID_PROPERTIES
    use overwrite_module,       only: PRESCRIBE_VELOCITY
    use truchas_logging_services

    implicit none

    ! Local Variables
    integer                                      :: status
    logical(log_kind)                            :: abort

    ! Argument List

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (.not. fluid_flow) return

    ALLOCATE (Solid_Face(nfc,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('FLUID_FLOW: Solid_Face(nfc,ncells) allocation failed')
    ALLOCATE (isPureImmobile(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('FLUID_FLOW: isPureImmobile(ncells) allocation failed')
    ALLOCATE (fluidDeltaRho(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('FLUID_FLOW: fluidDeltaRho(ncells) allocation failed')
    ALLOCATE (Face_Density(nfc,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('FLUID_FLOW: Face_Density(nfc,ncells) allocation failed')

    call fluid_properties(abort)


    call PRESCRIBE_VELOCITY('diagxz')

    DEALLOCATE (Solid_Face)
    DEALLOCATE (isPureImmobile)
    DEALLOCATE (fluidDeltaRho)
    DEALLOCATE (Face_Density)

    return

  END SUBROUTINE appliedflowfield

END MODULE FLUID_FLOW_MODULE
