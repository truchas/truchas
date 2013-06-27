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
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: slopeLimiter, fluxLimiterThuburn, limitGradient

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

    use gs_module, only: EE_GATHER
    use mesh_module, only: Cell, Mesh
    use parameter_module, only: ncells, ndim, nfc

    ! Argument List
    real(r8), dimension(:), intent(IN) :: Phi
    real(r8), dimension(:,:), intent(INOUT) :: gradPhi
    character(*), optional, intent(IN) :: method

    ! Local Variables
    integer :: f,i,status
    real(r8), dimension(:,:), allocatable :: Flux_Phi
    real(r8), dimension(:,:), allocatable :: Ngbr_flux_Phi
    real(r8), dimension(:),   allocatable :: Psi, Tmp1, Tmp2
    real(r8), dimension(:,:), allocatable :: Ngbr_Phi

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
    use cutoffs_module,   only: alittle
    use parameter_module, only: ncells
    use mesh_module,      only: Cell

    ! Arguments
    real(r8), dimension(:),   intent(IN)  :: Extrapolated_Value
    real(r8), dimension(:),   intent(IN)  :: Cell_Value
    real(r8), dimension(:,:), intent(IN)  :: Neighbor_Values
    real(r8), dimension(:),   intent(OUT) :: Slope_Limiter
    integer, intent(IN)  :: element
    character(*), optional, intent(IN) :: location
    character(*), optional, intent(IN) :: method

    ! Local Variables
    logical :: specified_method, specified_location
    character(80) :: limiter_method = 'Venkat', limiter_location = 'face'
    integer :: status
    real(r8) :: Venkat_constant = 1.0_r8 / 3.0_r8
    real(r8), dimension(:), allocatable :: Tmp1, Tmp2, Tmp3, Tmp4

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

          Slope_Limiter = 1.0_r8

       ! Venkatakrishnan's Limiter (AIAA Paper #AIAA-93-0880)
       case ('Venkat')

          Tmp1 = MAX(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp2 = MIN(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp3 = Extrapolated_Value - Cell_Value
          Tmp3 = SIGN(1.0_r8,Tmp3)*(ABS(Tmp3) + alittle)
          Tmp4 = (Cell%Volume)**(1.0_r8/3.0_r8)  ! Approximation of cell length
          Tmp4 = (Venkat_constant*Tmp4)**3
          where (Tmp3 > alittle)
             Slope_Limiter = (Tmp3*(Tmp1**2 + Tmp4) + 2.0_r8*Tmp1*Tmp3**2) / &
                         (Tmp1**2 + 2.0_r8*Tmp3**2 + Tmp1*Tmp3 + Tmp4 + alittle)
          elsewhere
             Slope_Limiter = (Tmp3*(Tmp2**2 + Tmp4) + 2.0_r8*Tmp2*Tmp3**2) / &
                         (Tmp2**2 + 2.0_r8*Tmp3**2 + Tmp2*Tmp3 + Tmp4 + alittle)
          end where
          Slope_Limiter = MERGE(1.0_r8, Slope_Limiter/Tmp3, ABS(Tmp3) <= alittle)

       ! Barth's Limiter (AIAA Paper #AIAA-89-0366)
       case ('Barth')

          Tmp1 = MAX(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp2 = MIN(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp3 = Extrapolated_Value - Cell_Value
          Tmp3 = MERGE(alittle, Tmp3, ABS(Tmp3) <= alittle)
          where (Tmp3 >= alittle)
             Slope_Limiter = MIN(1.0_r8, Tmp1/Tmp3)
          elsewhere
             Slope_Limiter = MIN(1.0_r8, Tmp2/Tmp3)
          end where
          Slope_Limiter = MERGE(1.0_r8, Slope_Limiter, ABS(Tmp3) <= alittle)

    end select

    ! Make sure the limiter is bounded by zero and one
    Slope_Limiter = MIN(1.0_r8, MAX(0.0_r8, Slope_Limiter))

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
    use parameter_module, only: ncells, nfc
    use gs_module,        only: EE_GATHER

    ! Arguments
    real(r8), dimension(:),     intent(IN)     :: phi
    real(r8), dimension(:,:),   intent(INOUT)  :: facePhi
    real(r8), dimension(:,:),   intent(IN)     :: Face_flx_vols
    real(r8), dimension(:,:),   intent(IN)     :: Fluxing_Velocity
    real(r8), dimension(:,:),   intent(IN)     :: ModFluxVolume
    real(r8), dimension(:),     intent(IN)     :: ModVolume
    real(r8), dimension(:),   intent(INOUT)    :: phimin
    real(r8), dimension(:),   intent(INOUT)    :: phimax

    logical, dimension(:,:), intent(IN) :: InflowMask

    ! Local Variables
    real(r8) :: phiInMin, phiInMax, fluxVolume
    integer :: f,i,status

    ! memory for neighbor values of cell quantities...    

    real(r8), dimension(:,:), allocatable :: phiUpMin_Ngbr
    real(r8), dimension(:,:), allocatable :: phiUpMax_Ngbr
    real(r8), dimension(:,:), allocatable :: facePhi_Ngbr

    real(r8), dimension(:),   allocatable :: phiUpMin   
    real(r8), dimension(:),   allocatable :: phiUpMax   
    real(r8), dimension(:),   allocatable :: phiMp1Min  
    real(r8), dimension(:),   allocatable :: phiMp1Max  
    real(r8), dimension(:),   allocatable :: sumVolIn   
    real(r8), dimension(:),   allocatable :: sumVolOut  
    real(r8), dimension(:),   allocatable :: sumPhiVolInMin
    real(r8), dimension(:),   allocatable :: sumPhiVolInMax
    real(r8), dimension(:),   allocatable :: PhiOutMin
    real(r8), dimension(:),   allocatable :: PhiOutMax
    real(r8), dimension(:,:),   allocatable :: Phi_Ngbr


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
    sumVolIn       = 0.0_r8
    sumVolOut      = 0.0_r8
    sumPhiVolInMin = 0.0_r8
    sumPhiVolInMax = 0.0_r8

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
       where (Fluxing_Velocity(f,:) < 0.0_r8)
          phiUpMin(:) = min(Phi_Ngbr(f,:), phiUpMin(:))
          phiUpMax(:) = max(Phi_Ngbr(f,:), phiUpMax(:))
       end where
    end do ! nfc loop...

    ! do our EE_Gathers up front...
    phiUpMin_Ngbr = 0.0_r8
    phiUpMax_Ngbr = 0.0_r8
    facePhi_Ngbr  = 0.0_r8
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
          if (Fluxing_Velocity(f,i) < 0.0_r8) then

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
          if (Fluxing_Velocity(f,i) > 0.0_r8) then
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
          if (Fluxing_Velocity(f,i) > 0.0_r8) then
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
