!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

MODULE INTERFACES_INPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define procedures for the input of interfaces namelist parameters.
  !
  ! Public Interface(s):
  !
  !   * call INTERFACES_INPUT()
  !
  !     Read and process the INTERFACES namelist.
  !
  ! Contains: INTERFACES_DEFAULT
  !           INTERFACES_INPUT
  !           INTERFACES_INPUT_PARALLEL
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use truchas_logging_services
  implicit none
  private

  public :: INTERFACES_INPUT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE INTERFACES_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default INTERFACES namelist.
    !
    !=======================================================================
    use interfaces_module, only: background_body, int_particles, &
                                 isrftype, Matnum, nsurf, vof_particles, &
                                 vof_method, &
                                 vof_tolerance, &
                                 vof_max_recursion
    use input_utilities,   only: NULL_R, NULL_I

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Constants
    isrftype = 0
    Matnum   = 0
    nsurf    = 0

    ! vof initialization method
    vof_method = 'default'

    ! division tolerance
    vof_tolerance = NULL_R

    ! maximum recursion depth
    vof_max_recursion = NULL_I

    ! Interface locator particles
    int_particles = 5

    ! Volume fraction particles
    vof_particles = 5

    ! Background body
    background_body = 0

  END SUBROUTINE INTERFACES_DEFAULT

  !-----------------------------------------------------------------------------

  SUBROUTINE INTERFACES_INPUT (lun)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   Read INTERFACES namelist. This routine sets body and interfaces 
    !   defaults. It reads and sets the variables listed in the /INTERFACES/
    !   namelist. It also calls body_inp, which reads and sets the body
    !   namelist.
    !
    !---------------------------------------------------------------------------
    use body_input_module,      only: BODY_CHECK, BODY_DEFAULT, BODY_INPUT
    use interfaces_module,      only: background_body,   &
                                      Matnum,            &
                                      nbody,             &
                                      vof_method,        &
                                      vof_tolerance,     &
                                      vof_max_recursion, &
                                      vof_particles,     &
                                      int_particles
    use input_utilities,        only: seek_to_namelist, NULL_R, NULL_I
    use parallel_info_module,   only: p_info
    use parameter_module,       only: mbody
    use property_data_module,   only: background_material
    use pgslib_module,          only: PGSLib_GLOBAL_ANY, pgslib_bcast
    use property_module,        only: Get_User_Material_ID
    use string_utilities,       only: lower_case
    use string_utilities, only: i_to_c

    integer, intent(in) :: lun

    ! local variables
    logical :: body_namelist, fatal, found
    integer :: ios, ib
    character(128) :: message

    ! interfaces namelist
    namelist /INTERFACES/ vof_method, vof_tolerance, vof_max_recursion, int_particles, vof_particles

    !---------------------------------------------------------------------------

    ! assume success
    fatal = .false.

    call TLS_info ('')
    call TLS_info ('Reading INTERFACES namelist ...')

    ! prepare to read
    call INTERFACES_DEFAULT  ()

    !! Locate the INTERFACES namelist (first occurence).
    if (p_info%IOP) then
      rewind lun
      call seek_to_namelist (lun, 'INTERFACES', found, iostat=ios)
    end if

    call pgslib_bcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading input file; iostat=' // i_to_c(ios))

    call pgslib_bcast (found)
    if (found) then
      if (p_info%IOP) read(lun, nml=interfaces, iostat=ios)
      call pgslib_bcast (ios)
      if (ios /= 0) call TLS_fatal ('error reading INTERFACES namelist; iostat=' // i_to_c(ios))
    else
      call TLS_info ('  INTERFACES namelist not found; using defaults')
    end if

    ! broadcast data
    call INTERFACES_INPUT_PARALLEL ()

    ! INTERFACES namelist error checking

    ! check for method
    select case (trim(lower_case(vof_method)))

    case ('default')
       ! this is where the default method is set
       vof_method = 'points'            ! default method

    case ('points')
       vof_method = 'points'
       ! should check for valid numbers in int_particles and vof_particles here
       ! but that's a code change that's not part of this exercise

    case ('divide')
       vof_method = 'divide'

    case default
       call TLS_fatal ('unknown value for VOF_METHOD: "' // trim(vof_method) // '"')

    end select

    ! method is known, do method-specific error checking
    select case (trim(vof_method))

    case ('points')
       ! make sure int_particles <= vof_particles, if not, correct
       if (vof_particles < int_particles) then
          write (message, 115)
115       format ('VOF_PARTICLES < INT_PARTICLES, reducing INT_PARTICLES (= VOF_PARTICLES)')
          call TLS_warn (message)
          int_particles = vof_particles
       end if

    case ('divide')

       ! if vof_max_recursion wasn't specified, set to default of 100
       if (vof_max_recursion == NULL_I) then
          vof_max_recursion = 100
       else if (vof_max_recursion < 0 .or. vof_max_recursion > 1000) then
          call TLS_fatal ('VOF_MAX_RECURSION must be between 0 and 1000')
       end if

       ! if it's still =NULL_R, we'll default it to something reasonable later
       if (vof_tolerance /= NULL_R .and. vof_tolerance < 0.0d0) then
          call TLS_fatal ('VOF_TOLERANCE must be >= 0.0')
       end if

    case default
       ! this should never happen
       INSIST(.false.)
    end select

    ! *** BODY namelists ***

    ! Initialize the number of bodies
    nbody = 0

    call BODY_DEFAULT ()

    if (p_info%IOP) rewind lun

    ! Read BODY geometry
    do ib = 1, mbody
       call BODY_INPUT (lun, found)
       if (.not. found) exit

       call BODY_CHECK (fatal)
       call TLS_fatal_if_any (fatal, 'error in BODY namelist #' // i_to_c(ib))

       if (Matnum(ib) == background_material) background_body = ib
    end do

    ! Make sure one body namelist has been read
    if (nbody == 0) call TLS_fatal ('no BODY namelists found')

    ! BODY namelist error checking

    ! Make sure the background body has been assigned
    if (background_body == 0) then
      write(message,'(a,i0,a)') 'background material ', &
          Get_User_Material_ID(background_material), ' does not have a BODY namelist'
      call TLS_fatal (message)
    end if

    ! If no fatal errors, then we got this far
    write (message, 210) background_body, Get_User_Material_ID(background_material)
210 format (9x,'BODY Namelist number ',i2,' will be used for background material (',i2,')')
    call TLS_info (message)

  END SUBROUTINE INTERFACES_INPUT

  !-----------------------------------------------------------------------------

  SUBROUTINE INTERFACES_INPUT_PARALLEL ()
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast components of the interfaces namelist to all PE's. If the
    !   interfaces namelist changes in interfaces_input, this routine and
    !   calls to it must be updated appropriately.
    !
    !======================================================================
    use interfaces_module,    only: vof_method,        &
                                    vof_tolerance,     &
                                    vof_max_recursion, &
                                    vof_particles,     &
                                    int_particles
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLib_BCAST

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast Data
    if (.not. p_info%UseGlobalServices) then
       call PGSLib_BCAST (vof_method)
       call PGSLib_BCAST (vof_tolerance)
       call PGSLib_BCAST (vof_max_recursion)
       call PGSLib_BCAST (int_particles)
       call PGSLib_BCAST (vof_particles)
    endif

  END SUBROUTINE INTERFACES_INPUT_PARALLEL

END MODULE INTERFACES_INPUT_MODULE
