!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  use truchas_logging_services
  use kinds, only: r8
  implicit none
  private

  public :: POROUS_DRAG

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
    use fluid_data_module,    only: fluidVof, isImmobile, Drag_Coefficient
    use matl_module,          only: GATHER_VOF
    use parameter_module,     only: nmat, ncells, ndim
    use porous_drag_data,     only: porous_implicitness
    use property_data_module, only: Permeability_Constant
    use timing_tree
    use zone_module,          only: Zone

    ! Argument List
    real(r8), intent(IN) :: dt
    real(r8), dimension(ndim,ncells), intent(INOUT) :: Mom_Delta

    ! Local Variables
    integer :: status
    integer :: m, n
    real(r8) :: tweight
    real(r8), dimension(:),   allocatable :: Solid_Vof
    real(r8), dimension(:,:), allocatable :: C

    ! Cutoff value for the drag coefficient
    real(r8) :: drag_cutoff = 1.0e+6

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start Timer
    call start_timer("Porous Drag")

    ALLOCATE (Solid_Vof(ncells), STAT = status)
    if (status /= 0) call TLS_panic ('POROUS_DRAG: Solid_Vof(ncells) allocation failed')
    ALLOCATE (C(ndim,ncells), STAT = status)
    if (status /= 0) call TLS_panic ('POROUS_DRAG: C(ndim,ncells) allocation failed')

    C = 0.0_r8

    tweight = 1.0_r8 - porous_implicitness

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
       where (fluidVof < 1.0_r8) C(n,:) = C(n,:)/(1.0_r8-fluidVof(:))
    end do

    ! Calculate A, according to Voller & Prakash
    ! Increment Velocity Delta for porous drag effects.
    do n = 1,ndim
       where (fluidVof > 0.0_r8) 
          Drag_Coefficient(n,:) = &
               MIN(C(n,:) * (1.0_r8 - fluidVof(:))**2 / fluidVof(:)**3, drag_cutoff)
          Mom_Delta(n,:) =  Mom_Delta(n,:) &
               - tweight* dt*fluidVof(:)*Drag_Coefficient(n,:)*Zone%Vc_old(n)
       endwhere
    end do

    DEALLOCATE (Solid_Vof, C)

    ! Stop Timer
    call stop_timer("Porous Drag")

  END SUBROUTINE POROUS_DRAG

END MODULE POROUS_DRAG_MODULE
