#include "f90_assert.fpp"

MODULE ADVECTION_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define the procedures for computing the advection term in
  !   mass, momentum, enthalpy, and species conservation equations.
  !
  ! Public Interface(s):
  !
  !   * call ADVECT_MASS
  !   * call ADVECT_MOMENTUM
  !   * call ADVECT_ENTHALPY
  !
  !     Main advection driver.
  !
  ! Contains: ADVECT_MASS
  !           UPDATE_MASS
  !           ADVECT_MOMENTUM
  !           ADVECT_MOMENTUM_ACCUMULATION
  !           ADVECT_MOMENTUM_DC
  !           ADVECT_MOMENTUM_HO
  !           ADVECT_ENTHALPY_ONCE
  !           ADVECT_ENTHALPY_DC
  !           ADVECT_ENTHALPY_HO
  !           ADVECT_ENTHALPY
  !           ADVECT_SPECIES
  !           ADVECT_SPECIES_DC
  !           ADVECT_SPECIES_HO
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !            Markus Bussmann (bussmann@lanl.gov)
  !            Edward D. Dendy (dendy@lanl.gov)
  !            Mark A. Christon (christon@lanl.gov)
  !
  !=======================================================================
  use kind_module, only: real_kind

  implicit none

  private
              
  ! Public Procedures
  public :: ADVECT_MASS, ADVECT_MOMENTUM
  public :: compute_advected_enthalpy, advected_phi

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE ADVECT_MASS
    !=======================================================================
    ! Purpose(s):
    !
    !   Call the volume tracker from here, and then update mass and
    !   concentration distributions.
    !
    !=======================================================================
    use advect_volume_module, only: ADVECT_VOLUME
    use advection_data,       only: Volume_Flux
    use bc_module,            only: IN_FLOW, OUT_FLOW
    use constants_module,     only: zero
    use fluid_data_module,    only: qin, qout, fluid_to_move, Fluxing_Velocity
    use kind_module,          only: int_kind, log_kind
    use matl_utilities,       only: MATL_GET_VOF, MATL_SET_VOF
    use mesh_module,          only: Cell
    use parameter_module,     only: ncells, nfc, nmat
    use pgslib_module,        only: PGSLib_Global_ANY, PGSLib_Global_SUM
    use time_step_module,     only: dt    
    use timing_tree
    use vof_data_module,      only: volume_track_interfaces, VT_Interface_Mask
    use truchas_logging_services

    ! Local Variables
    integer                                          :: status
    integer(int_kind)                                :: f
    logical(log_kind), dimension(:),     allocatable :: Mask
    real(real_kind),   dimension(:),     allocatable :: Tmp
    real(real_kind),   dimension(:,:),   allocatable :: Vof, Vof_n
   
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (.not.fluid_to_move) return

    ! Start Volume Advection Timer
    call start_timer("Mass Advection")
    
    ! Allocate working arrays
    ALLOCATE (Mask(ncells),                     &
              Tmp(ncells),                      &
              Vof(nmat,ncells),                 &
              Vof_n(nmat,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('ADVECT_MASS: allocation failed')

    if (volume_track_interfaces .and. .not. ASSOCIATED(VT_Interface_Mask)) then
       ALLOCATE (VT_Interface_Mask(nfc,ncells), STAT = status)
       if (status /= 0) call TLS_panic ('ADVECT_MASS: VT_Interface_Mask(nfc,ncells) allocation failed')
    end if

    if (.not. ALLOCATED(Volume_Flux)) then
       ALLOCATE (Volume_Flux(nmat,nfc,ncells), STAT=status)
       if (status /= 0) call TLS_panic ('ADVECT_MASS: Volume_Flux(nmat,nfc,ncells) allocation failed')
    end if

    ! Grab a copy of the volume fractions from Matl; will return time n+1 values to Matl at the end
    ! of this subroutine, via MATL_SET_VOF
    call MATL_GET_VOF (Vof)

    ! Will also need a copy of the time n Vof values.
    Vof_n = Vof

    ! Advect material volumes.
    call ADVECT_VOLUME (Vof, Vof_n, Fluxing_Velocity, Volume_Flux)

    ! Update the mass & concentration distributions.
    if (volume_track_interfaces) then 
       call UPDATE_MASS (Fluxing_Velocity, Vof, Vof_n, VT_Interface_Mask)
    else
       call UPDATE_MASS (Fluxing_Velocity, Vof, Vof_n)
    end if

    ! Return Vof values at this point back into the Matl structure.
    call MATL_SET_VOF (Vof)

    ! Accumulate Momentum Advection Array
    call ADVECT_MOMENTUM_ACCUMULATION ()

    ! Find and store the accumulated inflow and outflow mass.
    do f = 1, nfc
       Mask = IN_FLOW (f, Fluxing_Velocity) .OR. OUT_FLOW(f, Fluxing_Velocity)

       if (PGSLib_Global_ANY(Mask)) then
          ! Accumulate inflow volume - Fluxing_Velocity < zero
          Tmp = MIN(dt*Cell%Face_Area(f)*Fluxing_Velocity(f,:),zero)
          qin = qin + ABS(PGSLib_Global_SUM(Tmp, MASK = Mask))

          ! Accumulate outflow volume - Fluxing_Velocity > zero
          Tmp = MAX(dt*Cell%Face_Area(f)*Fluxing_Velocity(f,:),zero)
          qout = qout + PGSLib_Global_SUM(Tmp, MASK = Mask)

       end if
    end do

    DEALLOCATE (Mask)
    DEALLOCATE (Tmp)
    DEALLOCATE (Vof)
    DEALLOCATE (Vof_n)

    ! Stop the Volume Advection Timer
    call stop_timer("Mass Advection")

    return

  END SUBROUTINE ADVECT_MASS

  SUBROUTINE UPDATE_MASS (Fluxing_Velocity, Vof, Vof_n, VT_Interface_Mask)
    !=======================================================================
    ! Purpose(s):
    !
    !   Update the volume fractions.
    !
    !=======================================================================
    use constants_module,       only: zero
    use fluid_data_module,      only: fluidVof, fluidVof_n, fluidRho, fluidRho_n, isImmobile, &
                                      Cell_isnt_Void, Ngbr_isnt_Void
    use gs_module,              only: EE_GATHER
    use kind_module,            only: int_kind, log_kind, real_kind
    use parameter_module,       only: ncells, nfc, nmat
    use property_module,        only: DENSITY_MATERIAL
    use zone_module,            only: Zone
    use truchas_logging_services

    implicit none

    ! Arguments
    real(real_kind),   dimension(nfc,ncells),                intent(IN)    :: Fluxing_Velocity
    real(real_kind),   dimension(nmat,ncells),               intent(INOUT) :: Vof
    real(real_kind),   dimension(nmat,ncells),               intent(INOUT) :: Vof_n
    logical(log_kind), dimension(nfc,ncells),      optional, intent(INOUT) :: VT_Interface_Mask

    ! Local Variables
    integer                                        :: status
    integer(int_kind)                              :: m, n, f
    logical(log_kind), dimension(:),   allocatable :: Mask
    logical(log_kind), dimension(:,:), allocatable :: VT_Interface_Mask_Ngbr
    real(real_kind)                                :: m_density

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (Mask(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('UPDATE_MASS: Mask(ncells) allocation failed')

    ! Update Zone%Rho, and fluidRho and fluidVof at old and new times.
    Zone%Rho = zero
    fluidVof_n = zero
    fluidVof = zero
    fluidRho_n = zero
    fluidRho = zero

    do m = 1, nmat
       m_density = DENSITY_MATERIAL(m)
       where (Vof(m,:) > zero) Zone%Rho = Zone%Rho + Vof(m,:)*m_density
       if (isImmobile(m)) cycle
       do n = 1, ncells
          if (Vof_n(m,n) > zero) then
             fluidRho_n(n) = fluidRho_n(n) + Vof_n(m,n)*m_density
             fluidVof_n(n) = fluidVof_n(n) + Vof_n(m,n)
          end if
          if (Vof(m,n) > zero) then
             fluidRho(n) = fluidRho(n) + Vof(m,n)*m_density
             fluidVof(n) = fluidVof(n) + Vof(m,n)
          end if
       end do
    end do

    where (fluidVof_n > zero) fluidRho_n = fluidRho_n / fluidVof_n
    where (fluidVof   > zero) fluidRho   = fluidRho   / fluidVof

    ! Set void cell indicator arrays.
    Cell_isnt_Void = Zone%Rho > zero
    call EE_GATHER(Ngbr_isnt_Void, Cell_isnt_Void)

    ! Calculate VT_Interface_Mask, a flag for faces near a volume
    ! tracked interface, that have a positive Fluxing_Velocity.
    ! (This isn't the only definition we can imagine, but we define
    ! the mask as all faces of those cells that no longer contain
    ! the same single material.)

    if (PRESENT(VT_Interface_Mask)) then

       ALLOCATE (VT_Interface_Mask_Ngbr(nfc,ncells), STAT = status)
       if (status /= 0) call TLS_panic ('UPDATE_MASS: VT_Interface_Mask_Ngbr(nfc,ncells) allocation failed')
       Mask = .false.

       ! In cells where a material volume fraction has changed ...
       do m = 1,nmat
          where (Vof(m,:)-Vof_n(m,:) /= zero) Mask = .true.
       end do

       ! ... all faces should be masked the same way.
       do f = 1,nfc
          VT_Interface_Mask(f,:) = Mask
       end do

       ! The mask should also be the same from either side of a face.
       call EE_GATHER (VT_Interface_Mask_Ngbr, VT_Interface_Mask)
       do f = 1,nfc
          where (VT_Interface_Mask_Ngbr(f,:)) VT_Interface_Mask(f,:) = .true.
       end do

       ! Now limit the .true.'s to faces with a positive Fluxing_Velocity.
       do f = 1,nfc
          where (Fluxing_Velocity(f,:) <= zero) VT_Interface_Mask(f,:) = .false.
       end do

       DEALLOCATE (VT_Interface_Mask_Ngbr)

    end if

    DEALLOCATE (Mask)

    return

  END SUBROUTINE UPDATE_MASS
 
  SUBROUTINE ADVECT_MOMENTUM_ACCUMULATION ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the liquid momentum advection term.
    !
    !=======================================================================
    use advection_data,       only: advection_order_momentum, Momentum_Delta
    use time_step_module,     only: dt
    use timing_tree

    implicit none

    ! Local Varialbes

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start Timer
    call  start_timer("Momentum Advection")

    ! Advect Momentum.
    if (advection_order_momentum > 1) then
       ! call high order routine...
       call ADVECT_MOMENTUM_HO (Momentum_Delta)
    else
       ! call Markus Bussmann's first order routine...
       call ADVECT_MOMENTUM_DC (dt, Momentum_Delta)
    end if

    ! Stop Timer
    call stop_timer ("Momentum Advection")

    return

  END SUBROUTINE ADVECT_MOMENTUM_ACCUMULATION

  SUBROUTINE ADVECT_MOMENTUM_DC (dt, Momentum_Delta)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the liquid momentum advection term using donor-cell
    !
    !=======================================================================
    use advection_data,       only: Volume_Flux
    use bc_module,            only: BC_Vel, BC_Mat, BC_Temp, IN_FLOW
    use constants_module,     only: zero, preset
    use fluid_data_module,    only: Fluxing_Velocity
    use gs_module,            only: EE_GATHER
    use kind_module,          only: int_kind, log_kind
    use mesh_module,          only: Cell
    use parameter_module,     only: ncells, ndim, nfc, nmat
    use pgslib_module,        only: PGSLib_Global_ANY
    use time_step_module,     only: cycle_number
    use zone_module,          only: Zone
    use property_module,      only: DENSITY_MATERIAL
    use truchas_logging_services

    implicit none

    ! Argument List
    real(real_kind),  intent(IN)   :: dt
    real(real_kind),  intent(OUT)  :: Momentum_Delta(:,:)

    ! Local Variables
    integer                                        :: status
    integer(int_kind)                              :: f, m, n, nc
    logical(log_kind), dimension(:),   allocatable :: Mask
    real(real_kind),   dimension(:),   allocatable :: Inflow_Density, Inflow_Temp, Inflow_Vel, Tmp
    real(real_kind),   dimension(:),   allocatable :: Momentum_Delta_Component
    real(real_kind),   dimension(:,:), allocatable :: Momentum, Vc_Ngbr
    real(real_kind), dimension(:,:,:), allocatable :: Mass_Flux
    real(real_kind) :: rhom

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (Mask(ncells),                 &
              Momentum_Delta_Component(ncells),       &
              Inflow_Density(ncells),       &
              Inflow_Temp(ncells),          &
              Inflow_Vel(ncells),           &
              Momentum(nfc,ncells),         &
              Tmp(ncells),                  &
              Vc_Ngbr(nfc,ncells),          &
              Mass_Flux(nmat,nfc,ncells),   &
              STAT = status)
    call TLS_fatal_if_any (status /= 0, 'ADVECT_MOMENTUM_DC: allocation failed')

    Momentum_Delta = zero

    ! Calculate Mass_Flux using the previously computed Volume_Flux.
    do m = 1,nmat
       Mass_Flux(m,:,:) = Volume_Flux(m,:,:) * DENSITY_MATERIAL(m)
    end do
 
    ! Loop over each velocity component.
    DIMENSIONS: do n = 1, ndim

       Momentum = zero
       Tmp = Zone%Vc(n)
       call EE_GATHER(Vc_Ngbr, Tmp)

          ! Calculate Momentum change from Mass_Flux.

          do nc = 1,ncells
             do f = 1,nfc
                do m = 1,nmat
                   if(Mass_Flux(m,f,nc)>zero) then
                      Momentum(f,nc) = Momentum(f,nc) + Mass_Flux(m,f,nc)*Zone(nc)%Vc(n)
                   else
                      Momentum(f,nc) = Momentum(f,nc) + Mass_Flux(m,f,nc)*Vc_Ngbr(f,nc)
                   endif
                end do
             end do
          end do

       INFLOW_LOOP: do f = 1,nfc

          ! Find inflow faces at mesh boundaries.
          Mask = IN_FLOW (f, Fluxing_Velocity)

          if (PGSLIB_GLOBAL_ANY(Mask)) then

             Inflow_Temp = Zone%Temp
             if (ASSOCIATED(BC_Temp)) then
                where (BC_Temp(f,:) /= preset) Inflow_Temp = BC_Temp(f,:)
             end if

             Inflow_Density = Zone%Rho_Old
             do m = 1,nmat
                rhom = DENSITY_MATERIAL(m)
                where (BC_Mat(f,:) == m) Inflow_Density = rhom
             end do

             Inflow_Vel = Zone%Vc(n)
             where (BC_Vel(n,f,:) /= preset) Inflow_Vel = BC_Vel(n,f,:)


             ! fix up for 1D_plug_flow problem
             if (cycle_number == 1) then
                
                where (Mask .and. Zone(:)%Rho_old>zero) 
                    Momentum(f,:) = dt*Fluxing_Velocity(f,:)*Cell%Face_Area(f)*Zone%Vc_Old(n)*Zone%Rho_Old
                  elsewhere (Mask)
                    Momentum(f,:) = dt*Fluxing_Velocity(f,:)*Cell%Face_Area(f)*Inflow_Vel*Inflow_Density
                endwhere

             else
                where (Mask) Momentum(f,:) = dt*Fluxing_Velocity(f,:)*Cell%Face_Area(f)*Inflow_Vel*Inflow_Density
             end if

          end if

       end do INFLOW_LOOP

       ! Loop over faces to calculate Momentum_Delta_Component
       Momentum_Delta_Component = zero
       do f = 1,nfc
          Momentum_Delta_Component = Momentum_Delta_Component - Momentum(f,:)
       end do

       ! Increment the momentum delta by advection.
       Momentum_Delta(n,:) = Momentum_Delta(n,:) + Momentum_Delta_Component/Cell%Volume

    end do DIMENSIONS

    DEALLOCATE (Mask,             &
                Momentum_Delta_Component,   &
                Inflow_Density,   &
                Inflow_Temp,      &
                Inflow_Vel,       &
                Momentum,         &
                Tmp,              &
                Vc_Ngbr,          &
                Mass_Flux)

  END SUBROUTINE ADVECT_MOMENTUM_DC

  SUBROUTINE ADVECT_MOMENTUM_HO (Momentum_Delta)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the liquid momentum advection term using second-order scheme 
    !
    !=======================================================================
    use advection_data,       only: Volume_Flux
    use bc_module,            only: BC_Vel, IN_FLOW
    use constants_module,     only: zero, preset
    use fluid_data_module,    only: Fluxing_Velocity
    use gs_module,            only: EE_GATHER
    use kind_module,          only: int_kind, log_kind
    use parameter_module,     only: ncells, ndim, nfc, nmat
    use pgslib_module,        only: PGSLib_Global_ANY
    use time_step_module,     only: cycle_number
    use zone_module,          only: Zone
    use property_module,      only: DENSITY_MATERIAL
    use hoadvection,          only: ADVECT_SCALAR
    use truchas_logging_services

    implicit none

    ! Argument List
    real(real_kind),  intent(OUT)  :: Momentum_Delta(:,:)

    ! Local Variables
    integer                                        :: status
    integer(int_kind)                              :: f, m, n, nc
    real(real_kind)                                :: sum
    real(real_kind),   dimension(:),   allocatable :: Momentum_Delta_Component, Inflow_Density, Inflow_Temp, Inflow_Vel, Velocity
    real(real_kind),   dimension(:,:), allocatable :: Momentum, Vc_Ngbr
    real(real_kind), dimension(:,:,:), allocatable :: Mass_Flux
    real(real_kind),   dimension(:,:), allocatable :: ModFluxVolume

    logical(log_kind),   dimension(:,:), allocatable :: InflowMask
    real(real_kind),   dimension(:,:), allocatable :: InflowPhi

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (Momentum_Delta_Component(ncells),       &
              Inflow_Density(ncells),       &
              Inflow_Temp(ncells),          &
              Inflow_Vel(ncells),           &
              Momentum(nfc,ncells),         &
              Velocity(ncells),             &
              Vc_Ngbr(nfc,ncells),          &
              Mass_Flux(nmat,nfc,ncells),   &
              ModFluxVolume(nfc,ncells),    &
              InflowMask(nfc,ncells),       &
              InflowPhi(nfc,ncells),        &
              STAT = status)
    call TLS_fatal_if_any (status /= 0, 'ADVECT_MOMENTUM_HO: allocation failed')

    ! Calculate Mass_Flux using the previously computed Volume_Flux.
    do m = 1,nmat
       Mass_Flux(m,:,:) = Volume_Flux(m,:,:) * DENSITY_MATERIAL(m)
    end do

    Momentum_Delta = zero

    ! Loop over each velocity component.
    DIMENSIONS: do n = 1, ndim

       Momentum = zero
       Velocity = Zone%Vc(n)
       call EE_GATHER(Vc_Ngbr, Velocity)

       ! to set InflowMask and InflowPhi...
       INFLOW_LOOP: do f = 1,nfc

          InflowMask(f,:) = IN_FLOW(f, Fluxing_Velocity)

          if (PGSLIB_GLOBAL_ANY(InflowMask)) then
             InflowPhi(f,:) = Zone%Vc(n)
             where (BC_Vel(n,f,:) /= preset) InflowPhi(f,:) =  BC_Vel(n,f,:)
             if (cycle_number == 1) then
               where (InflowMask(f,:) .and. Zone(:)%Rho_old > zero) 
                 InflowPhi(f,:) = Zone%Vc(n)
               end where
             end if
          end if

       end do INFLOW_LOOP

       do nc = 1,ncells
          do f = 1,nfc
             sum = zero
             do m = 1,nmat             
                sum = sum + Mass_Flux(m,f,nc)
             end do
             ModFluxVolume(f,nc) = sum
          end do
       end do

       call ADVECT_SCALAR (Velocity, Zone%Rho_Old, Zone%Rho, ModFluxVolume, Fluxing_Velocity, InflowPhi, InflowMask)

       ! now to fit in with the rest of code calculate momentum_delta (MOM_NP1 - MOM_N)

       Momentum_Delta_Component = zero
       Momentum_Delta_Component = Velocity*Zone%Rho - Zone%Vc(n)*Zone%Rho_Old

       ! Increment the momentum delta by advection.
       Momentum_Delta(n,:) = Momentum_Delta(n,:) + Momentum_Delta_Component

    end do DIMENSIONS

    DEALLOCATE (Momentum_Delta_Component,       &
                Inflow_Density,       &
                Inflow_Temp,          &
                Inflow_Vel,           &
                Momentum,             &
                Velocity,             &
                Vc_Ngbr,              &
                Mass_Flux,            &
                ModFluxVolume,        &
                InflowMask,           &
                InflowPhi)

  END SUBROUTINE ADVECT_MOMENTUM_HO

  SUBROUTINE ADVECT_MOMENTUM (Mom_Delta)
    !=======================================================================
    ! Purpose(s):
    !
    !   Return the change in momentum due to advection.  This ncells
    !   quantity was computed during the advection part of the calculation, and
    !   is stored in Momentum_Delta, defined at the top of this module.
    !
    !=======================================================================
    use advection_data,       only: Momentum_Delta

    implicit none

    ! Argument List
    real(real_kind), dimension(:,:), intent(INOUT) :: Mom_Delta

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    Mom_Delta = Mom_Delta + Momentum_Delta

    return

  END SUBROUTINE ADVECT_MOMENTUM

!---------------mf----------------
  subroutine advected_phi(phi,delta_phi)
    
! advect scalar phi using DC
! dt 
    use constants_module,  only : zero
    use fluid_data_module, only : Fluxing_Velocity
    use gs_module,         only : EE_GATHER
    use kind_module,       only : int_kind,log_kind
    use parameter_module,  only : ncells,nfc
    use pgslib_module,     only : PGSLib_Global_Any
    use mesh_module,       only : Cell
    use time_step_module,  only : dt
    use truchas_logging_services

    implicit none

    real(real_kind), intent(in)  :: phi(:)
    real(real_kind), intent(out) :: delta_phi(:)


    integer(int_kind) :: status
    integer(int_kind) :: j,f 
    logical(log_kind), dimension(:),     allocatable :: Mask
    real(real_kind),   dimension(:),     allocatable :: Inflow_Conc
    real(real_kind),   dimension(:,:),   allocatable :: Conc_Flux
    real(real_kind),   dimension(:,:),   allocatable :: phi_ngbr 

    Allocate (Mask(ncells),            &
              Inflow_Conc(ncells),     &
              Conc_Flux(nfc,ncells),   & 
              phi_ngbr(nfc,ncells),    &
              STAT=status)
    call TLS_fatal_if_any (status /= 0, 'ADVECTED_PHI: allocation failed')

    ! get the neighbor values of phi 
    call EE_GATHER(phi_ngbr,phi)

    Conc_Flux = zero

    do j=1,ncells
      do f=1,nfc
        if (Fluxing_Velocity(f,j)>zero) then
          Conc_Flux(f,j) = Conc_Flux(f,j) + dt*Fluxing_Velocity(f,j)*phi(j)*Cell(j)%Face_Area(f)
        else
          Conc_Flux(f,j) = Conc_Flux(f,j) + dt*Fluxing_Velocity(f,j)*phi_ngbr(f,j)*Cell(j)%Face_Area(f)
        endif
      enddo !nfc
    enddo !ncells

        
!    TI_INFLOW_LOOP: do f = 1,nfc
!       ! Find inflow faces.
!       Mask = IN_FLOW (f, Fluxing_Velocity)
!       if (PGSLIB_GLOBAL_ANY(Mask)) then
!         Inflow_Conc = phi
!         if (ASSOCIATED(BC_Conc)) then
!           where (BC_Conc(f,:) /= preset) 
!              Inflow_Conc = BC_Conc(f,:)
!           end where
!         end if
!         where (Mask) 
!           Conc_Flux(f,:) = &
!                dt*Fluxing_Velocity(f,:)*Cell%Face_Area(f)*Inflow_Conc
!         end where
!       endif
!    end do TI_INFLOW_LOOP


    delta_phi = zero  
    do f=1,nfc
      delta_phi(:) = delta_phi(:) - Conc_Flux(f,:) !/Cell%Volume 
    enddo

    Deallocate(Mask, Inflow_Conc, Conc_Flux, phi_ngbr)

  end subroutine advected_phi

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! COMPUTE_ADVECTED_ENTHALPY
 !!
 !! A replacement for ADVECT_ENTHALPY for use by the diffusion solver only.
 !! This computes the net change in enthalpy for each cell due to the
 !! advection of material between cells as described by the VOLUME_FLUX
 !! module array.  In addition it returns the min and max temperatures of the
 !! the material parcels advected into (or remaining in) each cell.  These
 !! should give bounds on the new cell temperature.
 !!
 !! The procedure used is intended to replicate that of ADVECT_ENTHALPY_DC.
 !!
 !! NOTES:
 !! (1) BC are handled differently here.  The previously computed VOLUME_FLUX
 !!     array accounts for this (FLUX_BC subroutine), but the original
 !!     ADVECT_ENTHALPY_DC code replicated its BC handling for some reason.
 !!     I believe the result here should be identical.
 !! (2) I don't know how void is treated in VOLUME_FLUX, but I've been
 !!     careful to skip doing anything for void, which has no enthalpy anyway.
 !! (3) TODO: I'm continuing to treat the boundary inflow temperature as
 !!     before, but I think this needs to be changed to use the diffusion
 !!     solvers heat equation dirichlet boundary condition data.
 !! 

  subroutine compute_advected_enthalpy (Tcell, Hdelta, Tmin, Tmax)
  
    use kinds, only: r8
    use parameter_module, only: ncells, nfc, nmat
    use constants_module, only: preset
    use bc_module, only: IN_FLOW, BC_Temp
    use fluid_data_module, only: Fluxing_Velocity
    use advection_data, only: Volume_Flux
    use gs_module, only: EE_GATHER
    use material_interop, only: ds_enthalpy_density, void_material_index
    
    real(r8), intent(in) :: Tcell(:)
    real(r8), intent(out) :: Hdelta(:)
    real(r8), intent(out), optional :: Tmin(:), Tmax(:)
    
    integer :: k, j, m
    logical, allocatable :: mask(:)
    real(r8), allocatable :: Tnbr(:,:)
    real(r8) :: state(1), sum
    
    ASSERT(size(Tcell) == ncells)
    ASSERT(size(Hdelta) == ncells)
    ASSERT(present(Tmin) .eqv. present(Tmax))
    if (present(Tmin)) then
      ASSERT(size(Tmin) == ncells)
      ASSERT(size(Tmax) == ncells)
    end if
    
    !! Caution!  Volume_Flux may not yet be allocated because ADVECT_MASS
    !! has not yet been called with fluid_to_move==.true.
    if (.not.allocated(Volume_Flux)) then
      Hdelta = 0.0_r8
      if (present(Tmin)) then
        Tmin = Tcell
        Tmax = Tcell
      end if
      return
    end if
    
    !! Gather the neighbor cell temperatures.  I'm not sure what EE_GATHER does
    !! when there is no neighbor -- I'd guess it uses the cell temperature --
    !! but this data is only used when VOLUME_FLUX < 0, and for those faces we
    !! explicitly overwrite the value below (assuming VOLUME_FLUX < 0 when and
    !! only when FLUXING_VELOCITY < 0).
    allocate(Tnbr(nfc,ncells))
    call EE_GATHER(Tnbr, Tcell)
    
    !! At boundary inflow faces, overwrite the neighbor temperature with
    !! the boundary temperature, when specified, else the cell temperature.
    allocate(mask(ncells))
    do k = 1, nfc
      mask = IN_FLOW(k, Fluxing_Velocity)
      do j = 1, ncells
        if (mask(j)) then
          Tnbr(k,j) = Tcell(j)
          if (associated(BC_Temp)) then
            if (BC_Temp(k,j) /= preset) Tnbr(k,j) = BC_Temp(k,j)
          end if
        end if
      end do
    end do
    deallocate(mask)
    
    !! Accumulate the net change in enthalpy due to fluxed material.
    if (present(Tmin)) then
      do j = 1, ncells
        Tmin(j) = Tcell(j)
        Tmax(j) = Tcell(j)
        sum = 0.0_r8  ! net outflux of heat
        do k = 1, nfc
          do m = 1, nmat
            if (m == void_material_index) cycle
            if (Volume_Flux(m,k,j) == 0.0_r8) cycle
            if (Volume_Flux(m,k,j) < 0.0_r8) then ! influx
              state(1) = Tnbr(k,j)
              Tmin(j) = min(Tmin(j),Tnbr(k,j))
              Tmax(j) = max(Tmax(j),Tnbr(k,j))
            else  ! outflux
              state(1) = Tcell(j)
            end if
            sum = sum + Volume_Flux(m,k,j)*ds_enthalpy_density(m,state)
          end do
        end do
        Hdelta(j) = -sum
      end do
    else
      do j = 1, ncells
        sum = 0.0_r8  ! net outflux of heat
        do k = 1, nfc
          do m = 1, nmat
            if (m == void_material_index) cycle
            if (Volume_Flux(m,k,j) == 0.0_r8) cycle
            if (Volume_Flux(m,k,j) < 0.0_r8) then ! influx
              state(1) = Tnbr(k,j)
            else  ! outflux
              state(1) = Tcell(j)
            end if
            sum = sum + Volume_Flux(m,k,j)*ds_enthalpy_density(m,state)
          end do
        end do
        Hdelta(j) = -sum
      end do
    end if
    deallocate(Tnbr)
    
  end subroutine compute_advected_enthalpy

END MODULE ADVECTION_MODULE
