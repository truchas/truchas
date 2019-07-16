!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! NNC, October 2015.  The original base_types_module.F90 has been split
!! into two parts, one dealing with data belonging to the original mesh
!! data structure, and the other dealing with everything else (this file).
!! Changes have been kept to a minimum otherwise.

MODULE BASE_TYPES_A_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures to allocate, deallocate, and default base types
  !
  ! Public Interface(s):
  !
  !   * call BASE_TYPES_A_ALLOCATE ()
  !
  !     Allocate the base types.
  !
  !   * call BASE_TYPES_A_DEALLOCATE ()
  !
  !     Deallocate the base types.
  !
  ! Contains: BASE_TYPES_A_ALLOCATE
  !           BASE_TYPES_A_DEALLOCATE
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !            Robert Ferrell (CPCA, Ltd., ferrell@cpca.com)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  ! Public Procedures
  public :: BASE_TYPES_A_ALLOCATE, BASE_TYPES_A_DEALLOCATE

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE BASE_TYPES_A_ALLOCATE ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Allocate the Cell, Matl, and Zone base types and Fluxing_Velocity
    !
    !=======================================================================
    use bc_module,              only: BC
    use matl_module,            only: SLOT_INCREASE, Matl
    use legacy_mesh_api,        only: nfc, ncells
    use parameter_module,       only: mat_slot, mat_slot_new, nmat
    use zone_module,            only: Zone
    use fluid_data_module,      only: Fluxing_Velocity
    use solid_mechanics_module, only: SOLID_MECHANICS_ALLOCATE
    use solid_mechanics_input,  only: solid_mechanics
    use turbulence_module,      only: TURBULENCE_ALLOCATE
    use EM_data_proxy,          only: EM_is_on

    ! Arguments

    ! Local Variables
    integer :: memstat
    integer :: i
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Inform the user of allocation.
    call TLS_info ('')
    call TLS_info ('Allocating base derived types A ...', advance=.false.)

    ! Allocate the BC derived type.
    ALLOCATE (BC(ncells), STAT = memstat)
    if (memstat /= 0) call TLS_panic ('BASE_TYPES_A_ALLOCATE: BC derived type memory allocation error')

    ! BC base type.
    BC%Flag = 0
    BC%Internal = 0

    ! Allocate the Zone derived type.
    ALLOCATE (Zone(ncells), STAT = memstat)
    if (memstat /= 0) call TLS_panic ('BASE_TYPES_A_ALLOCATE: Zone derived type memory allocation error')

    !Allocate Fluxing Velocity variable
    ALLOCATE (Fluxing_Velocity(nfc,ncells), STAT = memstat)
    if (memstat /= 0) call TLS_panic ('BASE_TYPES_A_ALLOCATE: Fluxing_Velocity derived type memory allocation error')

    ! Allocate the Matl derived type.
    call SLOT_INCREASE (Matl, mat_slot, mat_slot_new)

    ! Allocate the material property, displacement, strain, and stress arrays
    call SOLID_MECHANICS_ALLOCATE ()

    ! allocate arrays for the turbulence model
    call TURBULENCE_ALLOCATE ()

    ! Set the new arrays to their defaults
    call BASE_TYPES_DEFAULT ()

    ! Inform the user of allocation.
    call TLS_info ('done.')

  END SUBROUTINE BASE_TYPES_A_ALLOCATE

  SUBROUTINE BASE_TYPES_A_DEALLOCATE ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Deallocate the base types.
    !
    !=======================================================================
    use bc_module,         only: BC
    use matl_module,       only: SLOT_DECREASE, Matl
    use parameter_module,  only: mat_slot, mat_slot_new
    use zone_module,       only: ZONE
    use fluid_data_module, only: Fluxing_Velocity

    ! Local Variables
    integer :: i
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! deallocate the BC derived type.
    if (ASSOCIATED(BC)) DEALLOCATE (BC)

    ! Deallocate the Zone derived type.
    if (ASSOCIATED(Zone)) DEALLOCATE (Zone)

    ! Deallocate Fluxing_Velocity.
    if (ASSOCIATED(Fluxing_Velocity)) DEALLOCATE (Fluxing_Velocity)

    ! Deallocate the Matl derived type by decreasing to zero slots.
    mat_slot_new = 0
    call SLOT_DECREASE (Matl, mat_slot, mat_slot_new)

  END SUBROUTINE BASE_TYPES_A_DEALLOCATE

  SUBROUTINE BASE_TYPES_DEFAULT ()
    !=======================================================================
    ! Purpose:
    !
    !   Default the base types.
    !
    !=======================================================================
    use parameter_module,  only: ndim
    use legacy_mesh_api,   only: nfc
    use zone_module,       only: Zone

    ! Local variables
    integer :: n
    
    ! Zone base type.
    Zone%Rho          = 0.0_r8
    Zone%Rho_old      = 0.0_r8
    Zone%Temp         = 0.0_r8
    Zone%Temp_old     = 0.0_r8
    Zone%Enthalpy     = 0.0_r8
    Zone%Enthalpy_old = 0.0_r8
    Zone%P            = 0.0_r8
    do n = 1,ndim
       Zone%Vc(n)     = 0.0_r8
       Zone%Vc_old(n) = 0.0_r8
    end do

  END SUBROUTINE BASE_TYPES_DEFAULT

END MODULE BASE_TYPES_A_MODULE
