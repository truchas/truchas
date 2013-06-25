MODULE LIMITER
  !=======================================================================
  ! Purpose(s):
  !
  !   Encapsulate all scalars, arrays, and procedures
  !   related to the slope-limiting algorithms
  !
  ! Public Interface(s):
  !
  !   * call LIMITER (Extrapolated_Value, Cell_Value, Neighbor_Values, &
  !                   Slope_Limiter, element, location, method)
  !
  ! Contains: LIMITER
  !
  !
  !=======================================================================
  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: slopeLimiter, fluxLimiterThuburn, limitGradient

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Limiter namelist variables
  ! character(LEN = 80), public, save :: limiter_type

CONTAINS

  SUBROUTINE limitGradient(Phi, gradPhi, method)

    !=======================================================================
    ! PURPOSE -
    !           Given the cell centered values of phi and its gradient
    !           evaluate the taylor series on faces then call slopelimiter
    !           to determine the multiplier which will ensure the slope
    !           is monotone.
    !      
    !=======================================================================

    use gs_module,                 only: EE_GATHER
    use kind_module,               only: real_kind, int_kind
    use mesh_module,               only: Cell, Mesh
    use parameter_module,          only: ncells, ndim, nfc
    use truchas_logging_services

    implicit none

    ! Argument List
    real(real_kind), dimension(:),    intent(IN)    :: Phi
    real(real_kind), dimension(:,:),  intent(INOUT) :: gradPhi
    character(LEN = *), optional,     intent(IN)    :: method

    ! Local Variables
    integer(int_kind)                             :: f,i,status
    real(real_kind), dimension(:,:), allocatable        :: Flux_Phi
    real(real_kind), dimension(:,:), allocatable        :: Ngbr_flux_Phi
    real(real_kind), dimension(:),   allocatable        :: Psi, Tmp1, Tmp2
    real(real_kind), dimension(:,:), allocatable        :: Ngbr_Phi

    ALLOCATE (Flux_Phi(nfc,ncells),      &
              Ngbr_flux_Phi(nfc,ncells), &
              Psi(ncells),               &
              Tmp1(ncells),              &
              Tmp2(ncells),              &
              Ngbr_Phi(nfc,ncells),      &
              STAT = status)
    call TLS_fatal_if_any (status /= 0, 'LIMITGRADIENT: allocation failed')


      call EE_GATHER(Ngbr_Phi, Phi)

      ! Compute the limiter PSI
      Psi = HUGE(Psi)
      do f = 1, nfc
        !Tmp1 = non limited face-extrapolated scalar value
        Tmp1 = Phi
        do i = 1,ndim

         Tmp1 = Tmp1 + (Cell%Face_Centroid(i,f) - Cell%Centroid(i))*gradPhi(i,:)

        end do
        where (Mesh%Ngbr_Cell(f) == 0)Ngbr_Phi(f,:) = Tmp1
        !The limiter is returned in Tmp2
        call slopeLimiter (Tmp1, Phi, Ngbr_Phi, Tmp2, f, 'face', method)
        !The final limiter is taken as the minimum for all faces
        Psi = MIN(Tmp2, Psi)
      end do
      ! Apply the limiter uniformly to the gradient
      do i = 1,ndim
        gradPhi(i,:)= Psi*gradPhi(i,:)
      end do

    DEALLOCATE (Flux_Phi, Ngbr_flux_Phi, Psi, Tmp1, Tmp2, Ngbr_Phi)

  END SUBROUTINE limitGradient

  SUBROUTINE slopeLimiter (Extrapolated_Value, Cell_Value, Neighbor_Values, &
                      Slope_Limiter, element, location, method)
    !=======================================================================
    ! PURPOSE -
    !           Given an value of a cell cnetered variable calculated using
    !           a Taylor series expansion calculate a multiplier to scale
    !           back the slope used so that it is bounded by the surrounding
    !           cell centered values.  Here we only look at face neighbors.
    !           Barth looked at vertices which for our purposes my be overly
    !           limiting.
    !      
    !=======================================================================
    use constants_module, only: one, one_third, zero, two
    use cutoffs_module,   only: alittle
    use kind_module,      only: int_kind, log_kind, real_kind
    use parameter_module, only: ncells
    use mesh_module,      only: Cell
    use truchas_logging_services
  
    implicit none

    ! Arguments
    real(real_kind), dimension(:),     intent(IN)  :: Extrapolated_Value
    real(real_kind), dimension(:),     intent(IN)  :: Cell_Value
    real(real_kind), dimension(:,:),   intent(IN)  :: Neighbor_Values
    real(real_kind), dimension(:),     intent(OUT) :: Slope_Limiter
    integer(int_kind),                 intent(IN)  :: element
    character(LEN = *), optional,      intent(IN)  :: location
    character(LEN = *), optional,      intent(IN)  :: method

    ! Local Variables
    logical(log_kind) :: specified_method, specified_location
    character(LEN = 80) :: limiter_method = 'Venkat', &
                           limiter_location = 'face'
    integer(int_kind)                  :: status
    real(real_kind)                    :: Venkat_constant = one_third
    real(real_kind), dimension(:), allocatable :: Tmp1, Tmp2, Tmp3, Tmp4

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

    ALLOCATE (Tmp1(ncells),      &
              Tmp2(ncells),      &
              Tmp3(ncells),      &
              Tmp4(ncells),      &
              STAT = status)
    call TLS_fatal_if_any (status /= 0, 'SLOPELIMITER: allocation failed')

    ! Process the incoming arguments
    specified_method   = PRESENT(method)
    specified_location = PRESENT(location)
    if (specified_method) limiter_method = method
    if (specified_location) limiter_location = location

    ! Compute the limiter according to the specified location and method
    select case(limiter_method)

       ! No limiting; set slope limiter to unity
       case default

          Slope_Limiter = one

       ! Venkatakrishnan's Limiter (AIAA Paper #AIAA-93-0880)
       case ('Venkat')

          Tmp1 = MAX(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp2 = MIN(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp3 = Extrapolated_Value - Cell_Value
          Tmp3 = SIGN(one,Tmp3)*(ABS(Tmp3) + alittle)
          Tmp4 = (Cell%Volume)**one_third  ! Approximation of cell length
          Tmp4 = (Venkat_constant*Tmp4)**3
          where (Tmp3 > alittle)
             Slope_Limiter = (Tmp3*(Tmp1**2 + Tmp4) + two*Tmp1*Tmp3**2) / &
                         (Tmp1**2 + two*Tmp3**2 + Tmp1*Tmp3 + Tmp4 + alittle)
          elsewhere
             Slope_Limiter = (Tmp3*(Tmp2**2 + Tmp4) + two*Tmp2*Tmp3**2) / &
                         (Tmp2**2 + two*Tmp3**2 + Tmp2*Tmp3 + Tmp4 + alittle)
          end where
          Slope_Limiter = MERGE(one, Slope_Limiter/Tmp3, ABS(Tmp3) <= alittle)

       ! Barth's Limiter (AIAA Paper #AIAA-89-0366)
       case ('Barth')

          Tmp1 = MAX(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp2 = MIN(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp3 = Extrapolated_Value - Cell_Value
          Tmp3 = MERGE(alittle, Tmp3, ABS(Tmp3) <= alittle)
          where (Tmp3 >= alittle)
             Slope_Limiter = MIN(one, Tmp1/Tmp3)
          elsewhere
             Slope_Limiter = MIN(one, Tmp2/Tmp3)
          end where
          Slope_Limiter = MERGE(one, Slope_Limiter, ABS(Tmp3) <= alittle)

    end select

    ! Make sure the limiter is bounded by zero and one
    Slope_Limiter = MIN(one, MAX(zero, Slope_Limiter))

    DEALLOCATE (Tmp1, Tmp2, Tmp3, Tmp4)

  END SUBROUTINE slopeLimiter

  SUBROUTINE fluxLimiterThuburn (phi, facePhi, Face_flx_vols, Fluxing_Velocity, &
                               ModfluxVolume, ModVolume, InflowMask, phimin, phimax)

    !=======================================================================
    ! PURPOSE -
    !           Given a flux quantity make sure that it is physically
    !           bounded and does not introduce new extrema.  The limiting
    !           is done such that the fluxed quantity is compatible with
    !           the mass.
    !     
    !=======================================================================
    use constants_module, only: zero
    use kind_module,      only: int_kind, log_kind, real_kind
    use parameter_module, only: ncells, nfc
    use gs_module,        only: EE_GATHER
    use truchas_logging_services

    implicit none

    ! Arguments
    real(real_kind), dimension(:),     intent(IN)     :: phi
    real(real_kind), dimension(:,:),   intent(INOUT)  :: facePhi
    real(real_kind), dimension(:,:),   intent(IN)     :: Face_flx_vols
    real(real_kind), dimension(:,:),   intent(IN)     :: Fluxing_Velocity
    real(real_kind), dimension(:,:),   intent(IN)     :: ModFluxVolume
    real(real_kind), dimension(:),     intent(IN)     :: ModVolume
    real(real_kind), dimension(:),   intent(INOUT)    :: phimin
    real(real_kind), dimension(:),   intent(INOUT)    :: phimax

    logical(log_kind), dimension(:,:), intent(IN)   :: InflowMask

    ! Local Variables
    real(real_kind)          :: phiInMin, phiInMax, fluxVolume
    integer(int_kind)        :: f,i,status

    ! memory for neighbor values of cell quantities...    

    real(real_kind), dimension(:,:), allocatable :: phiUpMin_Ngbr
    real(real_kind), dimension(:,:), allocatable :: phiUpMax_Ngbr
    real(real_kind), dimension(:,:), allocatable :: facePhi_Ngbr

    real(real_kind), dimension(:),   allocatable :: phiUpMin   
    real(real_kind), dimension(:),   allocatable :: phiUpMax   
    real(real_kind), dimension(:),   allocatable :: phiMp1Min  
    real(real_kind), dimension(:),   allocatable :: phiMp1Max  
    real(real_kind), dimension(:),   allocatable :: sumVolIn   
    real(real_kind), dimension(:),   allocatable :: sumVolOut  
    real(real_kind), dimension(:),   allocatable :: sumPhiVolInMin
    real(real_kind), dimension(:),   allocatable :: sumPhiVolInMax
    real(real_kind), dimension(:),   allocatable :: PhiOutMin
    real(real_kind), dimension(:),   allocatable :: PhiOutMax
    real(real_kind), dimension(:,:),   allocatable :: Phi_Ngbr


    ALLOCATE (phiUpMin_Ngbr(nfc,ncells),  &
              phiUpMax_Ngbr(nfc,ncells),  &
              facePhi_Ngbr(nfc,ncells),   &
              phiUpMin(ncells),        &
              phiUpMax(ncells),        &
              phiMp1Min(ncells),       &
              phiMp1Max(ncells),       &
              sumVolIn(ncells),        &
              sumVolOut(ncells),       &
              sumPhiVolInMin(ncells),  &
              sumPhiVolInMax(ncells),  &
              PhiOutMin(ncells),       &
              PhiOutMax(ncells),       &
              Phi_Ngbr(nfc,ncells),    &
              STAT = status)
    call TLS_fatal_if_any (status /= 0, 'FLUXLIMITER: allocation failed')

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

    ! initialize values to cell value...
    phiUpMin    = Phi
    phiUpMax    = Phi
    phiMp1Min   = Phi
    phiMp1Max   = Phi
    ! initialize values to zero...
    sumVolIn       = zero
    sumVolOut      = zero
    sumPhiVolInMin = zero
    sumPhiVolInMax = zero

    call EE_GATHER (Phi_Ngbr, Phi)

!!$    Flux_Phi coming into the hoadvection module should have been set the the incoming
!!$    fluxed quantity for inflow bc faces so that in the next step we are using the inflow
!!$    value to set Phi_ngbr.

    do f = 1,nfc
       where (InflowMask(f,:))
          Phi_Ngbr(f,:) = facePhi(f,:)
       end where
    end do
    
    do f = 1,nfc
       where (Fluxing_Velocity(f,:) < zero)
          phiUpMin(:) = min(Phi_Ngbr(f,:), phiUpMin(:))
          phiUpMax(:) = max(Phi_Ngbr(f,:), phiUpMax(:))
       end where
    end do ! nfc loop...

    ! do our EE_Gathers up front...
    phiUpMin_Ngbr = zero
    phiUpMax_Ngbr = zero
    facePhi_Ngbr  = zero
    call EE_GATHER (phiUpMin_Ngbr,phiUpMin)
    call EE_GATHER (phiUpMax_Ngbr,phiUpMax)
    call EE_GATHER (facePhi_Ngbr(:,:), facePhi(:,:))

!!$ Make sure that phiup min and max are set to inflow
!!$ value on inflow bc faces.
    do f = 1,nfc
       where (InflowMask(f,:))
          phiUpMin_Ngbr(f,:) = facePhi(f,:)
          phiUpMax_Ngbr(f,:) = facePhi(f,:)
          facePhi_Ngbr(f,:)  = facePhi(f,:)
       end where
    end do

    ! first look all faces as inflow faces...
    do i=1,ncells

       do f=1,nfc

          ! is face an inflow face...
          if (Fluxing_Velocity(f,i) < zero) then

             ! phi in (to the accecptor cell) min and max on this face...
             phiInMin = min(phiUpMin_Ngbr(f,i), phi(i))
             phiInMax = max(phiUpMax_Ngbr(f,i), phi(i))

             ! set down stream side of face value to upstream value...
             facePhi(f,i) = facePhi_Ngbr(f,i)

             ! now limit this face value based on phiInMax and phiInMin...
             facePhi(f,i) = min(facePhi(f,i),phiInMax)
             facePhi(f,i) = max(facephi(f,i),phiInMin)

             ! min and max future values...
             phiMp1Min(i) = min(phiMp1Min(i), phiInMin)
             phiMp1Max(i) = max(phiMp1Max(i), phiInMax)

             ! second refinement for courner coupling...
             phiInMin = min(phiUpMin_Ngbr(f,i),facePhi(f,i))
             phiInMax = max(phiUpMax_Ngbr(f,i),facePhi(f,i))

             ! fluxVolume = abs(Face_flx_vols(f,i))

             fluxVolume = abs(ModFluxVolume(f,i))

             sumVolIn(i)       = sumVolIn(i)       + fluxVolume
             sumPhiVolInMin(i) = sumPhiVolInMin(i) + fluxVolume * phiInMin;
             sumPhiVolInMax(i) = sumPhiVolInMax(i) + fluxVolume * phiInMax;

          end if

       end do ! nfc loop...
    end do ! ncells loop...



    ! now make sure facePhi is the same on both sides of a cell face...
    call EE_GATHER (facePhi_Ngbr(:,:), facePhi(:,:))

    ! now for sumVolOut...
    do i=1,ncells
       do f=1,nfc
          ! is face an inflow face...
          if (Fluxing_Velocity(f,i) > zero) then
             sumVolOut(i)          = sumVolOut(i)    +  abs(Face_flx_vols(f,i))
             ! make sure both sides of face have the same inflow limited value...
             facePhi(f,i)          = facePhi_Ngbr(f, i)
          end if
       end do ! nfc loop...
    end do ! ncells loop...

    ! calculate max and min outflow values...
    phiOutMin(:)=(sumPhiVolInMax(:)-phi(:)*(sumVolIn(:)-sumVolOut(:))- &
                 (phiMp1Max(:)-phi(:))*ModVolume(:))/                  &
                 (sumVolOut(:)+1.0e-32)**2*sumVolOut(:)
    phiOutMax(:)=(sumPhiVolInMin(:)-phi(:)*(sumVolIn(:)-sumVolOut(:))- &
                 (phiMp1Min(:)-phi(:))*ModVolume(:))/                  &
                 (sumVolOut(:)+1.0e-32)**2*sumVolOut(:)


    ! now look all faces as outflow faces...
    do i=1,ncells
       do f=1,nfc
          ! is face an outflow face...
          if (Fluxing_Velocity(f,i) > zero) then
             facePhi(f,i) = min(facePhi(f,i),phiOutMax(i))
             facePhi(f,i) = max(facePhi(f,i),phiOutMin(i))
          end if
       end do ! nfc loop...
    end do ! ncells loop...

    phimin = phiMp1Min
    phimax = phiMp1Max

    DEALLOCATE (phiUpMin_Ngbr,  &
         phiUpMax_Ngbr,  &
         facePhi_Ngbr,   &
         phiUpMin,        &
         phiUpMax,        &
         phiMp1Min,       &
         phiMp1Max,       &
         sumVolIn,        &
         sumVolOut,       &
         sumPhiVolInMin,  &
         sumPhiVolInMax,  &
         PhiOutMin,       &
         PhiOutMax,       &
         Phi_Ngbr)

  END SUBROUTINE fluxLimiterThuburn

END MODULE LIMITER
