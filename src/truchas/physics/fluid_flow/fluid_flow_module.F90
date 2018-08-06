!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: NAVIER_STOKES, FLUID_FLOW_DRIVER

CONTAINS

  SUBROUTINE FLUID_FLOW_DRIVER (t)

    use fluid_data_module, only: fluid_flow, applyflow, fluid_to_move

    real(r8), intent(in) :: t

    if (fluid_flow) then
      if (applyflow) then
        fluid_to_move = .true.
        call appliedflowfield(t)
      else
        call NAVIER_STOKES(t)
      end if
    end if

  END SUBROUTINE FLUID_FLOW_DRIVER

  SUBROUTINE NAVIER_STOKES (t)
    !=======================================================================
    ! Purpose(s):
    !
    !   Navier-Stokes (NS) driver: increment NS equations by one time step.
    !
    !======================================================================
    use fluid_data_module,      only: fluid_flow, fluidRho, Solid_Face, &
        isPureImmobile, fluidDeltaRho,      &
        fluid_to_move, Fluxing_Velocity
    use legacy_mesh_api,        only: ncells, nfc, ndim
    use predictor_module,       only: PREDICTOR
    use projection_data_module, only: Face_Density, mac_projection_iterations, &
        prelim_projection_iterations
    use projection_module,      only: PROJECTION
    use property_module,        only: FLUID_PROPERTIES
    use time_step_module,       only: cycle_number
    use viscous_data_module,    only: prelim_viscous_iterations, &
        viscous_iterations
    use zone_module,            only: Zone

    ! Local Variables
    integer :: status
    integer :: n
    logical :: abort

    ! Argument List
    real(r8), intent(in) :: t

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

    fluidRho = 0.0_r8

    ! Evaluate cell properties excluding immobile materials, and
    ! check that there are at least some flow equations to solve
    call FLUID_PROPERTIES (abort, t)

    if (.not. abort) then

      fluid_to_move = .true.
      print *, "<< Pre Predictor"
      do n = 1, ncells
        write(*,'("cell ",i4, " vel_cc: ", 2es15.5, " P_cc: ", es15.5)') &
            n, Zone(n)%Vc(1:2), Zone(n)%P
      end do
      ! Predictor Step
      call PREDICTOR()
      print *, "<< Post Predictor"
      do n = 1, ncells
        write(*,'("cell ",i4, " vel_cc: ", 2es15.5, " P_cc: ", es15.5)') &
            n, Zone(n)%Vc(1:2), Zone(n)%P
      end do
      ! Projection Step
      call PROJECTION()
      print *, "<< Post Projection"
      do n = 1, ncells
        write(*,'("cell ",i4, " vel_cc: ", 2es15.5, " P_cc: ", es15.5)') &
            n, Zone(n)%Vc(1:2), Zone(n)%P
      end do


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

    DEALLOCATE (Solid_Face)
    DEALLOCATE (isPureImmobile)
    DEALLOCATE (fluidDeltaRho)
    DEALLOCATE (Face_Density)

  END SUBROUTINE NAVIER_STOKES

  SUBROUTINE appliedflowfield (t)

    !=======================================================================
    ! Purpose(s):
    !
    !   Navier-Stokes (NS) driver: increment NS equations by one time step.
    !
    !======================================================================
    use fluid_data_module,      only: fluid_flow, Solid_Face, isPureImmobile, fluidDeltaRho
    use fluid_data_module,      only: fluxing_velocity
    use legacy_mesh_api,        only: ncells, nfc, cell
    use projection_data_module, only: Face_Density
    use property_module,        only: FLUID_PROPERTIES
    use overwrite_module,       only: PRESCRIBE_VELOCITY
    use zone_module, only: zone
    use advection_velocity_namelist, only: adv_vel

    real(r8), intent(in) :: t

    ! Local Variables
    integer :: j, k, status
    logical :: abort
    real(r8) :: args(0:3), vel(3)

    if (.not. fluid_flow) return

    ALLOCATE (Solid_Face(nfc,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('FLUID_FLOW: Solid_Face(nfc,ncells) allocation failed')
    ALLOCATE (isPureImmobile(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('FLUID_FLOW: isPureImmobile(ncells) allocation failed')
    ALLOCATE (fluidDeltaRho(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('FLUID_FLOW: fluidDeltaRho(ncells) allocation failed')
    ALLOCATE (Face_Density(nfc,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('FLUID_FLOW: Face_Density(nfc,ncells) allocation failed')

    call fluid_properties(abort, t)

    call PRESCRIBE_VELOCITY('deforming_sphere', t)

!!$    args(0) = t
!!$    do j = 1, ncells
!!$      args(1:3) = cell(j)%centroid
!!$      zone(j)%vc = adv_vel%eval(args)
!!$      do k = 1, nfc
!!$        args(1:3) = cell(j)%face_centroid(:,k)
!!$        vel = adv_vel%eval(args)
!!$        fluxing_velocity(k,j) = dot_product(cell(j)%face_normal(:,k), vel)
!!$      end do
!!$    end do
!!$
    DEALLOCATE (Solid_Face)
    DEALLOCATE (isPureImmobile)
    DEALLOCATE (fluidDeltaRho)
    DEALLOCATE (Face_Density)

  END SUBROUTINE appliedflowfield

END MODULE FLUID_FLOW_MODULE
