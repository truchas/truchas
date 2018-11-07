!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! NNC, Nov 2018. The legacy flow parts of time_step_module were moved here
!! with only minimal necessary changes.

module flow_time_step_module

  use kinds, only: r8
  implicit none
  private

  public :: legacy_flow_time_step

  real(r8), save, public :: courant_number ! Fluid Flow
  real(r8), save, public :: viscous_number ! Viscous Stress
  real(r8), save, public :: surften_number ! surface tension

contains

  subroutine legacy_flow_time_step(t, dt, dt_constraint)

    use legacy_mesh_api, only: ncells, nfc
    use fluid_data_module, only: FluidDeltaRho, Solid_face, isPureImmobile
    use property_module, only: fluid_properties
    use surface_tension_module, only: surface_tension
    use string_utilities, only: i_to_c

    real(r8), intent(in)  :: t
    real(r8), intent(out) :: dt
    character(:), allocatable, intent(inout) :: dt_constraint

    logical  :: abort

    real(r8) :: dt_courant, dt_viscous, dt_surften
    integer  :: min_dt_courant_cell(1), min_dt_viscous_cell(1), min_dt_surften_cell(1)

    dt = huge(dt)

    allocate(FluidDeltaRho(ncells),Solid_Face(nfc,ncells),isPureImmobile(ncells))
    call FLUID_PROPERTIES(abort, t)
    if(.not.abort) then
      call TIME_STEP_COURANT(dt_courant, min_dt_courant_cell)
      call TIME_STEP_VISCOUS(dt_viscous, min_dt_viscous_cell)
      dt = MIN(dt_courant, dt_viscous)
    else
      dt_courant = 1.0d10
      dt_viscous = 1.0d10
    endif
    deallocate(FluidDeltaRho,Solid_Face,isPureImmobile)

    if (surface_tension) then
      call TIME_STEP_SURFACE_TENSION(dt_surften, min_dt_surften_cell)
      dt = MIN(dt, dt_surften)
    else
      dt_surften = 1.0d10
    endif

    if (dt == dt_courant) then
      dt_constraint = 'courant ['//i_to_c(min_dt_courant_cell(1))//']'
    else if (dt == dt_viscous) then
      dt_constraint = 'viscous ['//i_to_c(min_dt_viscous_cell(1))//']'
    else if (dt == dt_surften) then
      dt_constraint = 'surface tension ['//i_to_c(min_dt_surften_cell(1))//']'
    end if

  end subroutine legacy_flow_time_step

  SUBROUTINE TIME_STEP_COURANT (dt_courant, min_dt_courant_cell)
    !=======================================================================
    ! Purpose(s):
    !   Compute the time step limit due to explicit advection. This
    !   time step restriction is limited by a CFL <= 1.0, based on
    !   Fluxing_Velocity, which is the face-normal component of the
    !   face flux velocities.
    !=======================================================================
    use legacy_mesh_api,   only: ncells, ndim, nfc, Cell
    use cutoffs_module,    only: alittle
    use fluid_data_module, only: fluxing_velocity, fluidvof, fluidRho, courant
    use PGSLib_module,     only: PGSLib_GLOBAL_MINLOC, PGSLib_GLOBAL_MINVAL, &
                                 PGSLib_GLOBAL_MAXLOC, PGSLib_GLOBAL_MAXVAL

    ! Argument List
    real(r8), intent(out) :: dt_courant
    integer,  intent(out) :: min_dt_courant_cell(1)

    ! Local Variables
    integer :: l, n, r, f
    real(r8) :: dt_next, max_Flux_Sum
    real(r8), dimension(ndim,ncells) :: Dl
    real(r8), dimension(ncells)      :: Local_Dt, Flux_Sum
    real(r8), dimension(ncells)      :: V_r, V_l

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Variables
    dt_courant          = 1.0d10
    min_dt_courant_cell = 0

    ! Compute Dl
    call TIME_STEP_DISTANCE (Dl)

    ! Courant Time Step: dt = C*dl/V
    ! (Courant Number: C = V*dt/dl)
    do n = 1, ndim
       ! Array Indices
       r = 2*n ; l = r-1

       ! Outward Normal Face Velocities
       !V_r = MERGE(0.0_r8, Fluxing_Velocity(r,:), Fluxing_Velocity(r,:) < 0.0_r8)
       !V_l = MERGE(0.0_r8, Fluxing_Velocity(l,:), Fluxing_Velocity(l,:) < 0.0_r8)
       V_r = max(0.0_r8, Fluxing_Velocity(r,:))
       V_l = max(0.0_r8, Fluxing_Velocity(l,:))

       ! Local Courant Time Step
       Local_Dt = courant_number * Dl(n,:) / (V_r + V_l + alittle)

!!$       ! The timestep should not be based on cells that are filled with void
!!$       where (FluidRho == 0.0_r8) Local_Dt = 1.0d10

       ! Minimum Courant Time Step
       dt_next = PGSLib_GLOBAL_MINVAL(Local_Dt)

       ! Global Courant Time Step
       dt_courant = MIN(dt_courant,dt_next)
       if (dt_courant == dt_next) then
          min_dt_courant_cell = PGSLib_GLOBAL_MINLOC(Local_Dt)
       end if
    end do

    ! The above routine does allow for the possibility that the sum of
    ! positive (outward) Fluxing_Velocity * Face_Area * dt_courant exceeds
    ! the volume of the cell; need to check, and reduce dt_courant
    ! accordingly.

    Flux_Sum = 0.0_r8
    do f = 1, nfc
       Flux_Sum(:) = Flux_Sum(:) + MAX(Fluxing_Velocity(f,:),0.0_r8)*Cell(:)%Face_Area(f)*dt_courant
    end do
    where (fluidVof > 0.0_r8) Flux_Sum = Flux_Sum / (Cell%Volume * fluidVof)
    where (FluidVof == 0.0_r8) Flux_Sum = 0.0_r8
    max_Flux_Sum = PGSLib_GLOBAL_MAXVAL(Flux_Sum)
    if (max_Flux_Sum > 1.0_r8) then
       dt_courant = 0.99*dt_courant/max_Flux_Sum
       min_dt_courant_cell = PGSLib_GLOBAL_MAXLOC(Flux_Sum)
    end if

    do n = 1, ndim
       ! Array Indices
       r = 2*n ; l = r-1
       ! Outward Normal Face Velocities
       !V_r = MERGE(0.0_r8, Fluxing_Velocity(r,:), Fluxing_Velocity(r,:) < 0.0_r8)
       !V_l = MERGE(0.0_r8, Fluxing_Velocity(l,:), Fluxing_Velocity(l,:) < 0.0_r8)
       V_r = max(0.0_r8, Fluxing_Velocity(r,:))
       V_l = max(0.0_r8, Fluxing_Velocity(l,:))

       where (fluidRho > 0.0_r8)
          courant = max(courant,(V_r+V_l+alittle)*dt_courant/Dl(n,:))
       elsewhere
          courant = 0.0_r8
       end where

    end do

  END SUBROUTINE TIME_STEP_COURANT

  SUBROUTINE TIME_STEP_VISCOUS (dt_viscous, min_dt_viscous_cell)
    !=======================================================================
    ! Purpose(s):
    !   Compute the time step limit due to explicit viscous stress terms
    !   in the Navier-Stokes Equations. The time step is computed using
    !   the Viscous Number Vn = nu*dt/(dl)**2 where the kinematic viscosity
    !   is given by nu = mu/rho.
    !=======================================================================
    use fluid_data_module, only: FluidRho
    use legacy_mesh_api,   only: ncells, ndim
    use PGSLib_module,     only: PGSLib_GLOBAL_MINLOC, PGSLib_GLOBAL_MINVAL
    use property_module,   only: get_viscosity
    use viscous_data_module,    only: inviscid

    real(r8), intent(out) :: dt_viscous
    integer,  intent(out) :: min_dt_viscous_cell(1)

    ! Local Variables
    integer :: n
    real(r8) :: dt_next

    real(r8), dimension(ndim,ncells) :: Dl_Squared
    real(r8), dimension(ncells)      :: Local_Dt
    real(r8), dimension(ncells)      :: Mu_Cell

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Variables
    dt_viscous          = 1.0d10
    min_dt_viscous_cell = 0

    if(.not.inviscid .and. viscous_number > 0.0_r8) then
      ! Update Properties
      call get_viscosity (Mu_Cell)

      ! Compute Dl_Squared
      call TIME_STEP_DISTANCE_SQUARED (Dl_Squared)


      ! Viscous Time Step: dt = Nv*(dl)**2/(Nu)
      ! (Viscous Number: Vn = Nu*dt/(dl)**2)
      do n = 1, ndim
       ! Local Viscous Time Step
       where(Mu_Cell > 0.0_r8)
           Local_Dt = viscous_number * Dl_Squared(n,:) * FluidRho / Mu_Cell
       elsewhere
           Local_Dt = 1.0d10
       endwhere

       ! Minimum Viscous Time Step
       dt_next = PGSLib_GLOBAL_MINVAL(Local_Dt)

       ! Global Viscous Time Step
       dt_viscous = MIN(dt_viscous, dt_next)
       if (dt_viscous == dt_next) then
          min_dt_viscous_cell = PGSLib_GLOBAL_MINLOC(Local_Dt)
       end if
      end do
    endif

  END SUBROUTINE TIME_STEP_VISCOUS

  SUBROUTINE TIME_STEP_DISTANCE (Dl)
    !=======================================================================
    ! Purpose(s):
    !   Compute the square of a characteristic length for a cell used in
    !   calculating a time step constraint.
    !=======================================================================
    use legacy_mesh_api, only: ncells, ndim, Cell, orthogonal_mesh

    ! Argument List
    real(r8), dimension(ndim,ncells), intent(OUT) :: Dl

    ! Local Variables
    integer :: l, n, r

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (orthogonal_mesh) then
       dl = 0.0_r8
       do n = 1, ndim
          r = 2*n ; l = r-1             ! Array Indices
          Dl(n,:) = Cell%HalfWidth(r) + Cell%HalfWidth(l) ! Face-to-Face Distance
       end do
    else
       call time_step_distance_squared(dl)
       dl(:,:) = sqrt(dl(:,:))
    end if

  END SUBROUTINE TIME_STEP_DISTANCE

  SUBROUTINE TIME_STEP_DISTANCE_SQUARED (Dl_Squared)
    !=======================================================================
    ! Purpose(s):
    !   Compute the square of a characteristic length for a cell
    !   used in calculating a time step constraint.
    !=======================================================================
    use legacy_mesh_api, only: ncells, ndim, Cell, orthogonal_mesh, Mesh, GAP_ELEMENT_1

    ! Argument List
    real(r8), dimension(ndim,ncells), intent(OUT) :: Dl_Squared

    ! Local Variables
    integer :: l, n, r, s

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Variables
    Dl_Squared = 0.0_r8

    ! Compute Dl_Squared
    do n = 1, ndim
       ! Array Indices
       r = 2*n ; l = r-1

       ! Dl_Squared
       if (orthogonal_mesh) then
          ! Orthogonal Mesh
          Dl_Squared(n,:) = (Cell%HalfWidth(r) + Cell%HalfWidth(l))**2
       else
          ! Non-Orthogonal Mesh
          do s = 1, ndim
             Dl_Squared(n,:) = Dl_Squared(n,:) + (Cell%Face_Centroid(s,r) - Cell%Face_Centroid(s,l))**2
          end do
       end if
       ! Do not use gap element dimensions (some of which are zero) for time step control
       where (Mesh(:)%Cell_Shape >= GAP_ELEMENT_1)
          Dl_Squared(n,:) = 1.0e12
       end where
    end do

  END SUBROUTINE TIME_STEP_DISTANCE_SQUARED

  SUBROUTINE TIME_STEP_SURFACE_TENSION (dt_surften, min_dt_surften_cell)
    !=======================================================================
    ! Purpose(s):
    !   Compute the time step limit due to surface tension forces
    !   Here note it is only due to the normal force component
    !   but it is currently used for both normal and tangential forces cases
    !   dt = coeff*sqrt(rho*dx^3/(2*pi*sigma)) with coeff between 0 and 1
    !
    !  Author(s)name : Marianne M. Francois (mmfran@lanl.gov)
    !           date : June 2005
    !  Reviewer name : Jim Sicilian (sicilian@lanl.gov)
    !           date : ??/??/2005
    !
    !=======================================================================
    use legacy_mesh_api,             only: ncells, ndim
    use constants_module,            only: pi, big
    use fluid_data_module,           only: fluidRho
    use PGSLib_module,               only: PGSLib_GLOBAL_MINLOC, PGSLib_GLOBAL_MINVAL
    use zone_module,                 only: Zone
    use surface_tension_module,      only: sigma_func, csf_boundary

    real(r8), intent(out) :: dt_surften
    integer,  intent(out) :: min_dt_surften_cell(1)

    ! Local Variables
    real(r8) :: state(1)
    real(r8), dimension(ncells) :: Sigma
    real(r8), dimension(ndim,ncells) :: Dl
    real(r8), dimension(ncells) :: Dl_min
    real(r8), dimension(ncells) :: Dt_local
    integer :: n

    !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get the cell-centered surface tension coefficient.
    do n = 1, ncells
      state(1) = Zone(n)%Temp_Old
      Sigma(n) = sigma_func%eval(state)
    end do

    ! Compute Dl
    call TIME_STEP_DISTANCE (Dl)
    ! find minimum
    do n=1,ncells
      Dl_min(n) = MIN (Dl(1,n),Dl(2,n),Dl(3,n))
    enddo

    ! Compute Surface Tension time step
    do n=1,ncells
      if (FluidRho(n) == 0.0_r8) then
        Dt_local(n)=big
      else
        Dt_local(n)=sqrt((FluidRho(n)*Dl_min(n)**3)/(2.0_r8*pi*Sigma(n)))
      endif
    enddo

    dt_surften = surften_number * PGSLib_GLOBAL_MINVAL(Dt_local)
    min_dt_surften_cell = PGSLib_GLOBAL_MINLOC(Dt_local)

    ! WARNING!!! Experimental: Deactivate/override time step limit due to
    ! tangential surface tension only for the special experimental case of
    ! boundary-applied surface tension force (csf_boundary)
    if (csf_boundary) dt_surften = 1.0d10

  END SUBROUTINE TIME_STEP_SURFACE_TENSION

end module flow_time_step_module
