MODULE POROUS_DRAG_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables and routines necessary for modeling
  !   momentum loss due to porous media.
  !
  ! Contains: POROUS_DRAG 
  !
  ! Author(s): Doug Kothe (dbk@lanl.gov)
  !
  !=======================================================================
! use kind_module, only: log_kind

  implicit none

  ! Private Module
  private

  ! Public Subroutines
  public :: POROUS_DRAG

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! PHYSICS namelist variables -
  ! Flag for enabling/disabling the porous drag model.
! logical(KIND = log_kind), public, save :: porous_flow

CONTAINS

  SUBROUTINE POROUS_DRAG (dt, Mom_Delta)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute momentum change due to flow through a mushy zone, following
    !   the simple model presented by V.R. Voller and C. Prakash, "A Fixed
    !   Grid Numerical Modelling Methodology for Convection-Diffusion Mushy
    !   Region Phase-Change Problems,", Int. J. Heat Mass Transfer, 30(8),
    !   pp. 1709-1719.
    !
    !   The current model imparts a momentum change given by:
    !
    !         Mom_Delta = -dt * A * velocity
    !
    !   where A = - Permeability_Constant * (1 - porosity)**2 / porosity**3
    !
    !  and where we allow the Permeability_Constant to be a function of 
    !  orientation.
    !
    !======================================================================= 
    use constants_module,     only: zero, one
    use fluid_data_module,    only: fluidVof, isImmobile, Drag_Coefficient
    use kind_module,          only: real_kind, int_kind
    use matl_module,          only: GATHER_VOF
    use parameter_module,     only: nmat, ncells, ndim
    use porous_drag_data,     only: porous_implicitness
    use property_data_module, only: Permeability_Constant
    use timing_tree
    use zone_module,          only: Zone
    use truchas_logging_services

    implicit none

    ! Argument List
    real(KIND=real_kind),                         intent(IN)    :: dt
    real(KIND=real_kind), dimension(ndim,ncells), intent(INOUT) :: Mom_Delta

    ! Local Variables
    integer                                           :: status
    integer(KIND=int_kind)                            :: m, n
    real(KIND=real_kind)                              :: tweight
    real(KIND=real_kind), dimension(:),   allocatable :: Solid_Vof
    real(KIND=real_kind), dimension(:,:), allocatable :: C

    ! Cutoff value for the drag coefficient
    real(KIND=real_kind) :: drag_cutoff = 1.0e+6

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start Timer
    call start_timer("Porous Drag")

    ALLOCATE (Solid_Vof(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('POROUS_DRAG: Solid_Vof(ncells) allocation failed')
    ALLOCATE (C(ndim,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('POROUS_DRAG: C(ndim,ncells) allocation failed')

    C = zero

    tweight = one - porous_implicitness

    ! This is a hack (for the moment), but because the code allows
    ! Permeability_Constant's to be a function of material, I'll average
    ! over different materials in each cell to get an average C(ndim,ncells)
    do m = 1,nmat
       if (.not.isImmobile(m)) CYCLE
       call GATHER_VOF (m, Solid_Vof)
       do n = 1,ndim
          C(n,:) = C(n,:) + Solid_Vof(:) * Permeability_Constant(n,m)
       end do
    end do
    do n = 1,ndim
       where (fluidVof < one) C(n,:) = C(n,:)/(one-fluidVof(:))
    end do

    ! Calculate A, according to Voller & Prakash
    ! Increment Velocity Delta for porous drag effects.
    do n = 1,ndim
       where (fluidVof > zero) 
          Drag_Coefficient(n,:) = &
               MIN(C(n,:) * (one - fluidVof(:))**2 / fluidVof(:)**3, drag_cutoff)
          Mom_Delta(n,:) =  Mom_Delta(n,:) &
               - tweight* dt*fluidVof(:)*Drag_Coefficient(n,:)*Zone%Vc_old(n)
       endwhere
    end do

    DEALLOCATE (Solid_Vof, C)

    ! Stop Timer
    call stop_timer("Porous Drag")

    return

  END SUBROUTINE POROUS_DRAG

END MODULE POROUS_DRAG_MODULE
