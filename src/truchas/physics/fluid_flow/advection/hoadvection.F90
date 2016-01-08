!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE HOADVECTION
  !=======================================================================
  ! Purpose(s):
  !
  !   Define the procedures for computing the advectivej update for a
  !   general scalar such as velocity, energy, or species in a manner that
  !   is compatible with the mass.
  !
  ! Public Interface(s):
  !
  ! Contains:
  !           ADVECT_GEOM
  !           ADVECT_OBJECT
  !           ADVECT_SCALAR
  !           UPDATE_OBJECT
  !
  !           Author(s): Edward D. Dendy (dendy@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: ADVECT_SCALAR

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! PHYSICS namelist variable: Flag for running in scalar advection mode.
  logical, public, save :: scalar_advection
  ! PHYSICS namelist variable: Flag for running in donor cell or first
  ! order advection mode
  logical, public, save :: donor_cell_advection

  character(80), public, save :: limiter_type
  character(80), public, save :: limiter_name

CONTAINS

  SUBROUTINE fluxedQuantity (Phi, Zone, Fluxing_Velocity, &
                             Flux_Phi, ModFluxVolume, ModVolume, phimin, phimax, InflowMask, &
                             limitingtype,limitingname)
    !=======================================================================
    ! Purpose(s):
    !
    !   Calculate a fluix quantity which is compatible with the mass
    !
    !======================================================================= 
    use gs_module,                 only: EE_GATHER
    use parameter_module,          only: ncells, ndim, nfc
    use time_step_module,          only: dt
    use zone_module,               only: CELL_AVG
    use limiter,                   only: fluxLimiterThuburn, limitGradient
    use discrete_op_module,        only: GRADIENT_CELL

    ! Argument List
    real(r8), dimension(:), intent(INOUT) :: Phi
    type(CELL_AVG), dimension(:), intent(INOUT) :: Zone
    real(r8), dimension(:,:), intent(IN)    :: Fluxing_Velocity
    real(r8), dimension(:,:), intent(INOUT) :: Flux_Phi
    real(r8), dimension(:,:), intent(IN)    :: ModFluxVolume
    real(r8), dimension(:),   intent(IN)    :: ModVolume
    real(r8), dimension(:),   intent(INOUT) :: phimin, phimax
    character(*), intent(IN) :: limitingtype
    character(*), intent(IN) :: limitingname
    logical, dimension(:,:), intent(IN) :: InflowMask

    ! Local Variables
    integer :: status
    real(r8), dimension(:,:),   allocatable :: Ngbr_Flux_vel
    real(r8), dimension(:,:),   allocatable :: Face_flx_vols
    real(r8), dimension(:,:,:), allocatable :: Face_flx_cents
    real(r8), dimension(:),     allocatable :: dPhiV
    real(r8), dimension(:),     allocatable :: othervof
    real(r8), dimension(:),     allocatable :: vfdotn
    real(r8), dimension(:,:),   allocatable :: GradPhi

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (Ngbr_Flux_vel(nfc,ncells),            &
              Face_flx_vols(nfc,ncells),            &
              Face_flx_cents(ndim,nfc,ncells),      &
              dPhiV(ncells),                        &
              othervof(ncells),                     &
              vfdotn(ncells),                       &
              GradPhi(ndim,ncells),                 &
              STAT = status)
    call TLS_fatal_if_any (status /= 0, 'FLUXEDQUANTITY: allocation failed')

    ! *************************************************
    ! the text string names used in limiter_type are
    ! fluxlimit
    ! slopelimit
    ! the text string names used in limiter_name are
    ! Thuburn
    ! Barth
    ! Venkat
    ! *************************************************

    limiter_type = limitingtype
    limiter_name = limitingname

    if (limiter_type .eq. 'donorcell') then
      donor_cell_advection = .true.
    end if

    Ngbr_Flux_vel    = 0.0_r8

    ! Gather face neighbor fluxing velocities 
    call EE_GATHER(Ngbr_Flux_vel(:,:), Fluxing_Velocity(:,:))
       
    ! get geometric quantities for advection
    call ADVECT_GEOM(Zone, Fluxing_Velocity, Face_flx_vols, &
                     Face_flx_cents,  dt)
    ! calculate a raw un-limited gradient
    call GRADIENT_CELL(GradPhi, Phi)

    ! if we are slope limiting do so now
    if (limiter_type .eq. 'slopelimit') then
      call limitGradient(Phi, gradPhi, limiter_name)
    end if

    ! use slopes to calculate value of phi at the centroid of the flux volume
    call ADVECT_OBJECT(Zone, Fluxing_Velocity, Face_flx_vols, Face_flx_cents,  &
                       Phi, Flux_Phi, GradPhi, InflowMask)

    ! if we are flux limiting do so now
    if (limiter_type .eq. 'fluxlimit') then

       select case(limiter_name)

       case ('nolimiting')
          ! no limiting...

       case ('Thuburn')
          ! return in place bounded/montone Flux_Phi
          call fluxLimiterThuburn(Phi, Flux_Phi, Face_flx_vols, Fluxing_Velocity, &
               ModFluxVolume, ModVolume, InflowMask, phimin, phimax)

        case default
          status = 1
          call TLS_fatal_if_any (status /= 0, 'FLUXEDQUANTITY: no flux limiter chosen')

       end select

    end if

    DEALLOCATE (Ngbr_Flux_vel, Face_flx_vols, Face_flx_cents, dPhiV, othervof, vfdotn, GradPhi)

  END SUBROUTINE fluxedQuantity

  SUBROUTINE ADVECT_OBJECT (Zone, Fluxing_Velocity, Face_flx_vols, Face_flx_cents,  &
                            Phi, Flux_Phi, GradPhi, InflowMask)                           
    !=======================================================================
    ! Purpose(s):
    !
    !   Calculates the cell gradient and uses it to calculate phi at the
    !   center of the characteristic based flux volume.  Gradient can
    !   be limited or not.
    !
    !======================================================================= 
    use gs_module,                 only: EE_GATHER
    use parameter_module,          only: ncells, ndim, nfc
    use zone_module,               only: CELL_AVG

    ! Argument List
    type(CELL_AVG), dimension(:), intent(INOUT) :: Zone
    real(r8), dimension(:,:),   intent(IN)    :: Fluxing_Velocity
    real(r8), dimension(:,:),   intent(IN)    :: Face_flx_vols
    real(r8), dimension(:,:,:), intent(IN)    :: Face_flx_cents
    real(r8), dimension(:),     intent(IN)    :: Phi
    real(r8), dimension(:,:),   intent(INOUT) :: Flux_Phi
    real(r8), dimension(:,:),   intent(IN)    :: GradPhi
    logical,  dimension(:,:),   intent(IN)    :: InflowMask

    ! Local Variables
    integer :: f,n,i,status
    real(r8), dimension(:,:), allocatable :: Ngbr_flux_Phi
    real(r8), dimension(:),   allocatable :: Psi, Tmp1, Tmp2
    real(r8), dimension(:,:), allocatable :: Ngbr_Phi

    ALLOCATE(Ngbr_flux_Phi(nfc,ncells),     &
             Psi(ncells),                   &
             Tmp1(ncells),                  &
             Tmp2(ncells),                  &
             Ngbr_Phi(nfc,ncells),          &
             STAT = status)
    call TLS_fatal_if_any (status /= 0, 'ADVECT_OBJECT: allocation failed')

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! initialize phi to its cell value...
    do f =1,nfc
       ! make sure we do not overwrite Flux_Phi at inflows
       where (.not.InflowMask(f,:))
          Flux_Phi(f,:) = Phi(:)
       end where
    enddo

    if(.not. donor_cell_advection)then

       call EE_GATHER(Ngbr_Phi, Phi)

       do i = 1,ndim
          do f = 1,nfc
             do n=1,ncells
                ! if fluxing velocity is less than zero do not calc since
                ! it would be unecessary and would overwrite boundary 
                ! initialization the int facePhi.
                if (Fluxing_Velocity(f,n) < 0.0_r8) cycle
                Flux_Phi(f,n) = Flux_Phi(f,n)+GradPhi(i,n)*Face_flx_cents(i,f,n) 
             end do
          end do
       end do

    end if

    DEALLOCATE(Ngbr_flux_Phi, Psi, Tmp1, Tmp2, Ngbr_Phi)

  END SUBROUTINE ADVECT_OBJECT                              

  SUBROUTINE ADVECT_GEOM (Zone, Fluxing_Velocity, Face_flx_vols, &
                     Face_flx_cents, dtcyc)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute evaluation point for calculating fluxed quantity.
    !
    !======================================================================= 
    use mesh_module,               only: Cell
    use parameter_module,          only: ncells, ndim, nfc
    use zone_module,               only: CELL_AVG

    ! Argument List
    type(CELL_AVG),  dimension(:), intent(INOUT) :: Zone
    real(r8), intent(IN) :: dtcyc
    real(r8), dimension(:,:),   intent(IN)  :: Fluxing_Velocity
    real(r8), dimension(:,:),   intent(OUT) :: Face_flx_vols
    real(r8), dimension(:,:,:), intent(OUT) :: Face_flx_cents
    
    ! Local Variables
    integer :: f,n,i,d,d1,status
    real(r8), dimension(:,:,:), allocatable :: vf

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE(vf(ndim,nfc,ncells), &
             STAT = status)
    call TLS_fatal_if_any (status /= 0, 'ADVECT_GEOM: allocation failed')
 
    Face_flx_vols = 0.0_r8
    Do f = 1,nfc
       Face_flx_vols(f,:) = Fluxing_Velocity(f,:) * Cell%Face_Area(f) * dtcyc
    End do


    !************************************************************************
    ! individual components of the face velocity from the Fluxing_Velocity
    ! and the upwind cell centered velocity.
    ! uf = (uf.n)*n + uc(upwind) - (uc(upwind).n)*n or in truchas notation
    ! uf = (Fluxing_Velocity)*n + uc(upwind) - (uc(upwind).n)*n

    do i=1,ncells
       do f=1,nfc
          do d=1,ndim
             vf(d,f,i) = Fluxing_Velocity(f,i) * Cell(i)%Face_Normal(d,f)
          end do

          do d=1,ndim
             vf(d,f,i) = vf(d,f,i) + Zone(i)%Vc(d)
             do d1=1,ndim
                vf(d,f,i) = vf(d,f,i) -                 &
                     Zone(i)%Vc(d1)*             &
                     Cell(i)%Face_Normal(d1,f)*  &
                     Cell(i)%Face_Normal(d,f)
             end do
          end do
       end do
    end do
    !************************************************************************

    If (.not. donor_cell_advection) then

       do i = 1,ndim
          do f = 1,nfc
             do n=1,ncells
                ! if fluxing velocity is less than zero do not calc since
                ! it would be unnecessary work.
                if (Fluxing_Velocity(f,n) < 0.0_r8) cycle
                Face_flx_cents(i,f,n) =   Cell(n)%Face_Centroid(i,f)     & 
                     - dtcyc * 0.5_r8                 & 
                     * vf(i,f,n)                      &
                     - Cell(n)%Centroid(i)
             end do
          end do
       end do

    End if

    DEALLOCATE(vf)

  END SUBROUTINE ADVECT_GEOM

  SUBROUTINE ADVECT_SCALAR (Phi, RhoN, RhoNP1, ModFluxVolume, Fluxing_Velocity, InflowPhi, InflowMask)
    !=======================================================================
    ! Purpose(s):
    !
    !           takes a scalar quantity from time n to time n+1.
    !
    !=======================================================================
    use gs_module,            only: EE_GATHER
    use mesh_module,          only: Cell
    use parameter_module,     only: ncells, nfc
    use zone_module,          only: Zone

    ! Argument List
    real(r8), dimension(:),   intent(INOUT) :: Phi
    real(r8), dimension(:,:), intent(IN)    :: Fluxing_Velocity
    real(r8), dimension(:),   intent(IN)    :: RhoN
    real(r8), dimension(:),   intent(IN)    :: RhoNP1
    real(r8), dimension(:,:), intent(IN)    :: ModFluxVolume
    real(r8), dimension(:,:), intent(IN)    :: InflowPhi
    logical,  dimension(:,:), intent(IN)    :: InflowMask

    ! Local Variables
    integer :: status
    integer :: f, nc
    real(r8) :: t1, t2, t3, vol
   
    real(r8), dimension(:,:), allocatable :: Flux_Phi, Flux_Phi_Ngbr
    real(r8), dimension(:),   allocatable :: ModVolume, phimin, phimax

    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
      ALLOCATE (Flux_Phi(nfc,ncells),      &
                Flux_Phi_Ngbr(nfc,ncells), &
                ModVolume(ncells),         &
                phimin(ncells),            &
                phimax(ncells),            &
                STAT = status)
      call TLS_fatal_if_any (status /= 0, 'ADVECT_SCALAR: allocation failed')

      ! modifiy the flux volume with the specified time n+1 density.
      ModVolume = RhoNP1*Cell%Volume

      ! Initialize the fluxed quantity to account for boundaries
      FLux_Phi(:,:) = 0.0_r8
      do f = 1,nfc
         where (InflowMask(f,:))
            Flux_Phi(f,:)  = InflowPhi(f,:)
         end where
      end do

      call fluxedQuantity (Phi, Zone, Fluxing_Velocity, &
           Flux_Phi(:,:), ModFluxVolume, ModVolume, phimin, phimax, InflowMask, 'fluxlimit', 'Thuburn')

      call EE_GATHER (Flux_Phi_Ngbr(:,:), Flux_Phi(:,:))

      ! need to set Flux_Phi_Ngbr = InflowPhi for boundary cell facess.
      do f = 1,nfc
         where (InflowMask(f,:))
            Flux_Phi_Ngbr(f,:) = InflowPhi(f,:)
         end where
      end do

      ! Now for the update...
      do nc = 1,ncells
         vol  = Cell(nc)%Volume
         t3 = 0.0_r8
         do f = 1,nfc
            if(Fluxing_Velocity(f,nc) > 0.0_r8 ) then
               t3 = t3 + ModFluxVolume(f,nc)*Flux_Phi(f,nc)
            elseif (Fluxing_Velocity(f,nc) < 0.0_r8) then
               t3 = t3 + ModFluxVolume(f,nc)*Flux_Phi_Ngbr(f,nc)
            endif
         end do
         ! this is formulated so that when mass is goes to zero
         ! the time advanced value is just the time n value.
         t1 = (Phi(nc)*vol*(RhoN(nc) - RhoNP1(nc))-t3)*RhoNP1(nc)*vol
         t2 = (RhoNP1(nc)*vol)*(RhoNP1(nc)*vol) + 1.0e-32 
         Phi(nc) =  Phi(nc) + t1 / t2

         ! this would be the standard advection update...
         !  Phi(nc) = Phi(nc)*RhoN(nc)/(RhoNP1(nc)) - t3 / (RhoNP1(nc) * vol)

      end do

      ! now need to make sure the future value is bounded
      ! by limiting future values from the Thuburn limiter.

      do nc=1,ncells
         Phi(nc) = min(phimax(nc), Phi(nc))
         Phi(nc) = max(phimin(nc), Phi(nc))
      end do

      ! cleanup

      DEALLOCATE (Flux_Phi, Flux_Phi_Ngbr, ModVolume, phimin, phimax)

  END SUBROUTINE ADVECT_SCALAR

END MODULE HOADVECTION
