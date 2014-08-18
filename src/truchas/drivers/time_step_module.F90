MODULE TIME_STEP_MODULE 
  !======================================================================= 
  ! Purpose(s): 
  ! 
  !   Time step and cycle information. 
  ! 
  !   Public Interface(s): 
  ! 
  !     * call TIME_STEP () 
  ! 
  !         Compute the new time step. 
  ! 
  ! 
  ! Contains: TIME_STEP 
  !           TIME_STEP_COURANT 
  !           TIME_STEP_VISCOUS 
  ! 
  !           TIME_STEP_DISTANCE 
  !           TIME_STEP_DISTANCE_SQUARED 
  ! 
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov) 
  ! 
  !======================================================================= 
  use kinds, only: r8
  implicit none 
  private 

  ! Public Subroutines 
  public :: TIME_STEP, TIME_STEP_COURANT 

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 

  ! Namelist Input Variables 
  ! Current Simulation Time 
  real(r8), save, public :: t = 0.0_r8 

  ! Cycle Numbers 
  integer, save, public :: cycle_max 
  integer, save, public :: cycle_number 

  ! Time Step Numbers 
  real(r8), save, public :: dt_constant 
  real(r8), save, public :: dt_init, dt_grow
  real(r8), save, public :: dt_max, dt_min 

  ! Time Step Limit Numbers 
  real(r8), save, public :: courant_number ! Fluid Flow 
  real(r8), save, public :: viscous_number ! Viscous Stress 
!  real(r8), save, public :: strain_limit   ! plastic strain 

  real(r8), save, public :: surften_number ! surface tension 

  ! Derived Quantities 
  character(LEN = 80),      save, public :: dt_constraint 
  logical, save, public :: constant_dt ! Constant dt flag 

  ! Heat transfer derived quantities. 
  real(r8), save, public :: dx_avg    ! avg cell width 

  integer,               save, public :: cycle_number_restart 
  integer, dimension(1), save, public :: min_dt_courant_cell 
  integer, dimension(1), save, public :: min_dt_viscous_cell 
  integer, dimension(1), save, public :: min_dt_surften_cell
  integer, dimension(1), save, public :: min_dt_overall_cell 

  real(r8), save, public :: dt            ! Current time step 
  real(r8), save, public :: dt_old        ! Previous time step
  real(r8), save, public :: dt_courant    ! Courant time step limit 
  real(r8), save, public :: dt_viscous    ! Viscous time step limit 
  real(r8), save, public :: dt_strain     ! plastic strain time step limit 
  real(r8), save, public :: dt_ds = huge(1.0d0)        ! diffusion solver time step limit 
  real(r8), save, public :: t1, t2        ! Pre and Post-Cycle times 
  real(r8), save, public :: dt_surften    ! surface tension time step limit

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 

CONTAINS 

  SUBROUTINE TIME_STEP () 
    !======================================================================= 
    ! Purpose(s): 
    !   Compute the time step. The time step is limited by constraints for 
    !   fluid-flow (advection and viscosity) and heat-transfer (conduction). 
    !======================================================================= 
    use fluid_data_module,        only: fluid_flow, FluidDeltaRho, & 
                                        Solid_face, isPureImmobile, Fluxing_Velocity
    use matl_module,              only: Matl 
    use parameter_module,         only: mat_slot, ncells, ndim, nfc 
    use property_module,          only: fluid_properties 
    use restart_variables,        only: restart
    use timing_tree
    use zone_module,              only: Zone
    use surface_tension_module,   only: surface_tension
    use diffusion_solver_data,    only: ds_enabled
    use truchas_logging_services

    ! Local Variables 
    integer :: n, s, status
    logical :: abort
    real(r8) :: dt_next, dt_growth
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 

    ! Start Time-Step Timer 
    call start_timer("Time Step")

    ! Initialize all old time data. 
    do n = 1, ndim 
       Zone%Vc_Old(n) = Zone%Vc(n) 
    end do 
    Zone%Rho_Old      = Zone%Rho 
    Zone%Temp_Old     = Zone%Temp 
    Zone%Enthalpy_Old = Zone%Enthalpy 
    do s = 1,mat_slot 
       Matl(s)%Cell%Vof_Old = Matl(s)%Cell%Vof 
    end do 

    ! Initialize Variables 
    min_dt_overall_cell = 0 
    dt_old = dt

    ! Time Step Limits 
    ! For the next time step, take the minimum of all constraints. 

    ! For the next time step, take the minimum of all constraints. 
    dt_growth = dt_grow*dt
    dt_next = MIN(dt_growth, dt_max)

    ! Diffusion solver time step limit.  The DT_DS value was returned by
    ! the solver on last step with a huge default initialized value to start.
    if (ds_enabled) dt_next = min(dt_next, dt_ds)

    ! Conduction and solidification front time steps.
    ! NNC: removed -- old HT

    ! Fluid-Flow Time Step: This must be done here because the ENTHALPY module 
    !                       may have changed the fluid properties after NAVIER_STOKES 
    !                       completed 
    if(fluid_flow) then 
        ! Evaluate cell properties excluding immobile materials, and 
        ! check that there are at least some flow equations to solve 
        ALLOCATE (Solid_Face(nfc,ncells), STAT = status) 
          if (status /= 0) call TLS_panic ('TIME_STEP: Solid_Face(nfc,ncells) allocation failed') 
        ALLOCATE (isPureImmobile(ncells), STAT = status) 
          if (status /= 0) call TLS_panic ('TIME_STEP: isPureImmobile(ncells) allocation failed') 
        ALLOCATE (fluidDeltaRho(ncells), STAT = status) 
          if (status /= 0) call TLS_panic ('TIME_STEP: fluidDeltaRho(ncells) allocation failed') 

        call FLUID_PROPERTIES (abort) 
        if(.not.abort) then 
            call TIME_STEP_COURANT(Fluxing_Velocity) 
            call TIME_STEP_VISCOUS 
            dt_next = MIN(dt_next,       & 
                          dt_courant,    & ! Courant 
                          dt_viscous)      ! Viscous 
          else 
            dt_courant = 1.0d10 
            dt_viscous = 1.0d10 
        endif 
      
        if (surface_tension) then
          call TIME_STEP_SURFACE_TENSION
          dt_next = MIN(dt_next, dt_surften)
        else
          dt_surften = 1.0d10 
        endif    
   

        DEALLOCATE (Solid_Face) 
        DEALLOCATE (isPureImmobile) 
        DEALLOCATE (fluidDeltaRho) 
    endif 

    ! Time step limit for plasticity 
!    if(solid_mechanics) then 
!       call TIME_STEP_PLASTICITY 
!       dt_next = MIN(dt_next,       & 
!            dt_strain) 
!    else 
       dt_strain = 1.0d10 
!    end if 

    ! Set a character string according to the constraint that's been activated. 
    if (dt_next == dt_growth) then 
       dt_constraint = 'growth' 
    else if (dt_next == dt_max) then 
       dt_constraint = 'maximum' 
    else if (dt_next == dt_ds) then
       dt_constraint = 'diffusion solver'
    else if (dt_next == dt_courant) then 
       min_dt_overall_cell = min_dt_courant_cell 
       dt_constraint = 'courant' 
    else if (dt_next == dt_viscous) then 
       min_dt_overall_cell = min_dt_viscous_cell 
       dt_constraint = 'viscous' 
    else if (dt_next == dt_strain) then 
       min_dt_overall_cell = min_dt_viscous_cell 
       dt_constraint = 'plastic strain' 
    else if (dt_next == dt_surften) then
       min_dt_overall_cell = min_dt_surften_cell
       dt_constraint = 'surface tension'
    end if 

    ! (Non) Constant Time Step 
    ! Constant Time Step 
    if (constant_dt) then 

       if (dt_next < dt_constant) then 
          write(message,10) dt_constant, TRIM(dt_constraint), dt_next 
10        format ('Constant time step of ',1pe13.5,' > ', & 
                  a,' time step constraint of ',1pe13.5)
          call TLS_warn (message) 
       end if 
       dt = dt_constant             ! Constant Time Step 
       dt_constraint = 'constant'   ! Time Step Constraint 

    ! Non-Constant Time Step 
    else 
       if (cycle_number - 1 == cycle_number_restart .and. .not.restart) then 
          ! First cycle; use initial time step 
          dt = dt_init 
          dt_constraint = 'initial' 
!          if (solid_mechanics) then 
!             if (dt_init > dt_strain) then 
!                Output_String =  blank_line 
!                write (Output_String, 15) dt, dt_strain 
!15              format (/,9x,'FATAL: Initial time step too large: dt = ',1pe13.5,' > dt_strain',/ & 
!                         ' restart with initial time step smaller than ',1pe13.5) 
!                call PUNT (Output_String, 'TIME_STEP') 
!             end if 
!          end if 
           if (surface_tension) then
            if (dt_init > dt_surften) then
              write (message, 17) dt, dt_surften
17            format ('Initial time step too large: dt = ',1pe13.5,' > dt_surften (',es13.5,')')
              call TLS_fatal (message)
            end if
          end if
       else 
          ! Non-First Cycle; use computed time step 
          dt = dt_next 
       end if 

    end if 

    ! Minimum Time Step 
    if (dt < dt_min) then 

       write (message, 30) dt_courant, min_dt_courant_cell
       call TLS_info (message)
30     format (15x,' dt_courant = ',1p,e12.3,5x,'at cell ',i8)

       write (message, 40) dt_viscous, min_dt_viscous_cell
       call TLS_info (message)
40     format (15x,' dt_viscous = ',1p,e12.3,5x,'at cell ',i8)

       write (message, 70) dt_surften, min_dt_surften_cell
       call TLS_info (message)
70     format (15x,' dt_surften = ',1p,e12.3,5x,'at cell ',i8)

       write (message, 80) dt_ds
       call TLS_info (message)
80     format (15x,' dt_ds = ',1p,e12.3)

       write (message, 20) dt 
20     format ('Time step too small: dt = ',1pe13.5,' < dt_min') 
       call TLS_fatal (message) 
    end if 

    ! Stop Time Step Timer 
    call stop_timer("Time Step") 

  END SUBROUTINE TIME_STEP 

  SUBROUTINE TIME_STEP_COURANT (Fluxing_Velocity) 
    !======================================================================= 
    ! Purpose(s): 
    !   Compute the time step limit due to explicit advection. This 
    !   time step restriction is limited by a CFL <= 1.0, based on 
    !   Fluxing_Velocity, which is the face-normal component of the 
    !   face flux velocities. 
    !======================================================================= 
    use mesh_module,       only: Cell 
    use cutoffs_module,    only: alittle 
    use fluid_data_module, only: fluidvof, fluidRho, courant
    use parameter_module,  only: ncells, ndim, nfc 
    use PGSLib_module,     only: PGSLib_GLOBAL_MINLOC, PGSLib_GLOBAL_MINVAL, & 
                                 PGSLib_GLOBAL_MAXLOC, PGSLib_GLOBAL_MAXVAL                    

    ! Argument List 
    real(r8), dimension(nfc,ncells), intent(IN) :: Fluxing_Velocity 

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

!!$       ! The timestep shouldn't be based on cells that are filled with void. 
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
    ! positive (outward) Fluxing_Velocity's * Face_Area * dt_courant exceeds 
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

  SUBROUTINE TIME_STEP_VISCOUS () 
    !======================================================================= 
    ! Purpose(s): 
    !   Compute the time step limit due to explicit viscous stress terms 
    !   in the Navier-Stokes Equations. The time step is computed using 
    !   the Viscous Number Vn = nu*dt/(dl)**2 where the kinematic viscosity 
    !   is given by nu = mu/rho. 
    !======================================================================= 
    use fluid_data_module, only: FluidRho 
    use parameter_module,  only: ncells, ndim
    use PGSLib_module,     only: PGSLib_GLOBAL_MINLOC, PGSLib_GLOBAL_MINVAL 
    use property_module,   only: get_viscosity
    use viscous_data_module,    only: inviscid 

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

!  SUBROUTINE TIME_STEP_PLASTICITY 
!    !======================================================================= 
!    ! Purpose(s): 
!    !   Compute the time step limit based on the plastic strain increment. 
!    !   The user specifies the (approximate) maximum change in 
!    !   elastic strain for a time step, a.  The plastic strain rate at  
!    !   the beginning of the time step, edot, is used to compute the maximum  
!    !   delta t by delta t = a/edot 
!    !    
!    !======================================================================= 
!    use constants_module, only: big 
!    Use node_operator_module, Only: nipc 
!    use solid_mechanics_data, only: SMech_IP 
!    use PGSLib_module,        only: PGSLib_GLOBAL_MAXLOC, PGSLib_GLOBAL_MAXVAL                    
! 
!    implicit none 
! 
!    ! Local variables 
!    real(r8) :: maxrate, testrate 
!    integer :: ip 
!     
!    maxrate = 0.0 
!    ! If strain rate is effectively zero 
!    do ip = 1,nipc 
!       testrate = PGSLib_GLOBAL_MAXVAL(SMech_IP(ip)%Plastic_Strain_Rate) 
!       if (testrate > maxrate) maxrate = testrate 
!    end do 
!    if (maxrate > 1e-12) then 
!          dt_strain = strain_limit / maxrate 
!    else 
!       dt_strain = big 
!    end if 
! 
!  END SUBROUTINE TIME_STEP_PLASTICITY 


  SUBROUTINE TIME_STEP_DISTANCE (Dl) 
    !======================================================================= 
    ! Purpose(s): 
    !   Compute the square of a characteristic length for a cell used in 
    !   calculating a time step constraint. 
    !======================================================================= 
    use mesh_module,      only: Cell, orthogonal_mesh 
    use parameter_module, only: ncells, ndim 

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
    use mesh_module,          only: Cell, orthogonal_mesh, Mesh, GAP_ELEMENT_1 
    use parameter_module,     only: ncells, ndim

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
  
  SUBROUTINE TIME_STEP_SURFACE_TENSION
    !=======================================================================
    ! Purpose(s):
    !   Compute the time step limit due to surface tension forces
    !   Here note it's only due to the normal force component 
    !   but it's currently used for both normal and tangential forces cases
    !   dt = coeff*sqrt(rho*dx^3/(2*pi*sigma)) with coeff between 0 and 1
    !
    !  Author(s)name : Marianne M. Francois (mmfran@lanl.gov)
    !           date : June 2005
    !  Reviewer name : Jim Sicilian (sicilian@lanl.gov)
    !           date : ??/??/2005
    ! 
    !=======================================================================
    use parameter_module,            only: ncells,ndim
    use constants_module,            only: pi, big
    use fluid_data_module,           only: fluidRho
    use PGSLib_module,               only: PGSLib_GLOBAL_MINLOC, PGSLib_GLOBAL_MINVAL
    use zone_module,                 only: Zone
    use surface_tension_module,      only: sigma_func

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

  END SUBROUTINE TIME_STEP_SURFACE_TENSION

END MODULE TIME_STEP_MODULE 












