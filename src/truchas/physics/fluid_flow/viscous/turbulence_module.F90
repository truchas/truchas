!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE TURBULENCE_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define turbulence related variables and routines.
  !
  ! Public Interface(s):
  !
  !   * call TURBULENCE (turbulence_model)
  !
  ! Contains: TURBULENCE
  !
  !           ALG_TURB_MODEL
  !
  ! Author(s): Kin L. Lam, LANL ESA-EA (klam@lanl.gov)
  !           
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: string_len
  use legacy_mesh_api, only: ncells
  implicit none
  private

  public :: TURBULENCE, TURBULENCE_ALLOCATE, read_turbulence_namelist_for_legacy

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Turbulent diffusivity at cell centers
  real(r8), allocatable, save, public :: Nu_Turb(:)

  ! Physics Namelist Variables
  character(string_len), save, public :: turbulence_model = 'none'
  real(r8), save, public :: turbulence_length
  real(r8), save, public :: turbulence_cmu
  real(r8), save, public :: turbulence_ke_fraction

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE TURBULENCE (specified_model)
    !=======================================================================
    ! Purpose(s):
    !
    !   Driver routine for all turbulence models.
    !
    !======================================================================= 

    ! Argument List
    character(*), optional, intent(IN) :: specified_model

    ! Local Variables
    character(string_len) :: model
    logical :: model_is_specified

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set the model to be applied.
    model_is_specified = PRESENT(specified_model)

    if (model_is_specified) then
       model = ADJUSTL(specified_model)
    else   
       model = 'alg' ! Algebraic model if optional argument not supplied
    end if

    ! Call the appropriate turbulence model
    select case (TRIM(model))

       case ('none')

          ! No turbulence; do nothing.

       case default

          ! Default: do nothing if specified model is neither 'none' nor 'alg'

       case ('alg')

          ! Algebraic model
          ! Initialize the turbulent diffusivity array
          Nu_Turb = 0.0_r8
          call ALG_TURB_MODEL ()

    end select

  END SUBROUTINE TURBULENCE

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE ALG_TURB_MODEL ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Calculate cell-centered turbulent diffusivity based on the
    !   algebraic turbulence model:
    !
    !   Nu_Turb = cmu * turbulence_length * turbulence_velocity
    !
    !   where turbulence_length is specified via input, and
    !   turbulence_velocity is determined from algebraic relations with the
    !   local velocities approximating the turbulence intensities.  No
    !   solution of transport equations is required.
    !            
    !   The calculated Nu_Turb (public variable accessible through use of
    !   turbulence_module) is to be used in appropriate places elsewhere
    !   in the code as follows:
    !
    !   Effective (Apparent) Diffusivity = Nu_Turb + Molecular Diffusivity
    !
    !   In the case of momentum transport, the diffusivity is also called
    !   kinematic viscosity (nu), which is given by mu / rho.
    !======================================================================= 
    use legacy_mesh_api,  only: ncells, ndim
    use zone_module,      only: Zone

    ! Local Variables
    integer :: n
    real(r8), dimension(ncells) :: Vel_Sq

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Calculate sum of cell-centered velocity squared
    Vel_Sq = 0.0_r8
    do n = 1, ndim
       Vel_Sq = Vel_Sq + Zone%Vc_old(n)**2
    end do

    ! Calculate turbulent diffusivity
    Nu_Turb =  turbulence_cmu * turbulence_length &
               * SQRT(0.5_r8 * turbulence_ke_fraction * Vel_Sq)

  END SUBROUTINE ALG_TURB_MODEL

  !-----------------------------------------------------------------------------

  subroutine turbulence_allocate ()

     ! use turbulence_model as a flag to determine whether Nu_Turb should be allocated

     if (turbulence_model == 'alg') then
        allocate(Nu_Turb(ncells))
        Nu_Turb = 0.0_r8
     end if

  end subroutine turbulence_allocate

  !! NNC, Oct 2018. Extracted out the namelist read to an independent module
  !! for use by the new flow implementation, and we piggyback on it here.

  subroutine read_turbulence_namelist_for_legacy(lun)
    use turbulence_namelist
    integer, intent(in) :: lun
    call read_turbulence_namelist(lun)
    if (associated(params)) then
      turbulence_model = 'alg'
      call params%get('length', turbulence_length)
      call params%get('cmu', turbulence_cmu, default=0.05_r8)
      call params%get('ke fraction', turbulence_ke_fraction, default=0.1_r8)
    end if
  end subroutine read_turbulence_namelist_for_legacy

END MODULE TURBULENCE_MODULE
