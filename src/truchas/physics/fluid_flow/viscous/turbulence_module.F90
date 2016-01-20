!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
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

  public :: TURBULENCE, TURBULENCE_ALLOCATE, read_turbulence_namelist

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Turbulent diffusivity at cell centers
  real(r8), dimension(:), pointer, save, public :: Nu_Turb

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
     use ArrayAllocate_Module, only: ARRAYCREATE

     ! use turbulence_model as a flag to determine whether Nu_Turb should be allocated

     nullify (Nu_Turb)

     if (turbulence_model == 'alg') then
        call ARRAYCREATE (Nu_Turb, 1, ncells, 'Allocation of Nu_Turb(ncells) failed')
        Nu_Turb = 0.0_r8
     end if

  end subroutine turbulence_allocate

  !!
  !! NNC, Jan 2012.  Moved turbulence model parameter input out of the PHYSICS
  !! namelist into its own namelist.  I've tried to maintained the original
  !! behavior, though I think some of it needs to be reconsidered. This routine
  !! should be called to read the TURBULENCE namelist only when viscous flow is
  !! enabled.  If the namelist is found, the algebraic turbulence model is used;
  !! otherwise it is not.
  !! 

  subroutine read_turbulence_namelist (lun)
  
    use kinds, only: r8
    use input_utilities
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use truchas_logging_services
    
    integer, intent(in) :: lun
    
    integer :: ios
    logical :: found
    character(128) :: string
    
    namelist /turbulence/ turbulence_length, turbulence_cmu, turbulence_ke_fraction
    
    !! Locate the TURBULENCE namelist (first occurence).
    if (is_IOP) then
      rewind lun
      call seek_to_namelist (lun, 'TURBULENCE', found, iostat=ios)
    end if

    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast (found)

    if (.not.found) return  ! the namelist is optional

    call TLS_info ('')
    call TLS_info ('Reading TURBULENCE namelist ...')

    !! Read the namelist.
    if (is_IOP) then
      turbulence_length = NULL_R
      turbulence_cmu = NULL_R
      turbulence_ke_fraction = NULL_R
      read(lun,nml=turbulence,iostat=ios)
    end if

    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('Error reading TURBULENCE namelist.')

    !! Broadcast the namelist variables.
    call broadcast (turbulence_length)
    call broadcast (turbulence_cmu)
    call broadcast (turbulence_ke_fraction)

    turbulence_model = 'alg'

    if (turbulence_length == NULL_R) then
      call TLS_fatal ('TURBULENCE_LENGTH must be assigned a value.')
    else if (turbulence_length <= 0.0_r8) then
      call TLS_fatal ('TURBULENCE_LENGTH must be > 0')
    end if

    if (turbulence_cmu == NULL_R) then
      turbulence_cmu = 0.05_r8
      write(string, '(a,es12.5)') '  using default TURBULENCE_CMU value: ', turbulence_cmu
      call TLS_info (string)
    else if (turbulence_cmu <= 0.0_r8) then
      call TLS_fatal ('TURBULENCE_CMU must be > 0')
    end if

    if (turbulence_ke_fraction == NULL_R) then
      turbulence_ke_fraction = 0.1_r8
      write(string, '(a,es12.5)') '  using default TURBULENCE_KE_FRACTION value: ', turbulence_ke_fraction
    else if (turbulence_ke_fraction <= 0.0_r8 .or. turbulence_ke_fraction >= 1.0_r8) then
      call TLS_fatal ('TURBULENCE_KE_FRACTION must be in (0,1)')
    end if

  end subroutine read_turbulence_namelist

END MODULE TURBULENCE_MODULE
